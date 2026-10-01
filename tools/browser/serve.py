"""
Local browser front end for generate_pMHC_SCT.py.

Serves one page on 127.0.0.1 where you upload a peptide CSV, pick the HLA /
leader / vector, enter a seed, and download the generated SCT constructs.

The page is only a front end. Every sequence is produced by calling
generate_sct_constructs() from tools/generate_pMHC_SCT.py, so the CSV
downloaded here is what the CLI writes for the same inputs. Nothing about the
construct is reimplemented in the browser.

Run:
    python tools/browser/serve.py

then open http://127.0.0.1:8765 (opened for you unless --no-browser is passed).

Stdlib only: no dependencies beyond what generate_pMHC_SCT.py already needs.
"""

from __future__ import annotations

import argparse
import io
import json
import re
import sys
import threading
import traceback
import uuid
import webbrowser
from contextlib import redirect_stdout
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

BROWSER_DIR = Path(__file__).resolve().parent
TOOLS_DIR = BROWSER_DIR.parent
ROOT_DIR = TOOLS_DIR.parent

# generate_pMHC_SCT imports config from the repo root and SCT_constants from
# tools/, so both have to be importable before it is imported.
for path in (ROOT_DIR, TOOLS_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import generate_pMHC_SCT as sct
from SCT_constants import HLA_nuc, destination_vector_nuc, leader_peptide_aa

PAGE = BROWSER_DIR / "index.html"
PREVIEW_ROWS = 50
MAX_UPLOAD_BYTES = 8 * 1024 * 1024

# numpy's global RNG is what the codon optimisation draws from, so two runs
# must never overlap in this process: concurrent reseeding would give each run
# sequences that neither seed reproduces.
_generation_lock = threading.Lock()

_jobs: dict[str, dict] = {}
_jobs_lock = threading.Lock()

AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


class InputError(Exception):
    """A problem with the uploaded CSV, worth reporting without a traceback."""


def _check_peptides(peptides) -> None:
    """Reject peptides the optimiser would only fail on obscurely.

    dnachisel reverse translates residue by residue, so a stray character
    surfaces as KeyError('3') and a blank cell as a TypeError about floats,
    neither of which points at the row to fix. Checked here rather than in
    generate_pMHC_SCT.py so the generator stays as it is on the command line.
    """
    blank, bad = [], []
    for position, (name, peptide) in enumerate(
            zip(peptides["name"], peptides["peptide_aa"]), start=1):
        # an unnamed row reads back from pandas as the string "nan"
        clean = str(name).strip()
        label = (f"row {position} ({clean!r})"
                 if clean and clean.lower() != "nan" else f"row {position}")
        if not isinstance(peptide, str) or not peptide.strip() or peptide == "NAN":
            blank.append(label)
            continue
        unknown = sorted(set(peptide) - AMINO_ACIDS)
        if unknown:
            bad.append(f"{label}: {', '.join(unknown)}")

    if blank:
        raise InputError(
            f"{len(blank)} row(s) have no peptide sequence: "
            f"{'; '.join(blank[:5])}"
            + (" ..." if len(blank) > 5 else "")
            + ". Remove them or fill the second column."
        )
    if bad:
        raise InputError(
            f"{len(bad)} row(s) contain characters that are not amino acids - "
            f"{'; '.join(bad[:5])}"
            + (" ..." if len(bad) > 5 else "")
            + ". Only the 20 standard residues (ACDEFGHIKLMNPQRSTVWY) can be "
            "reverse translated."
        )


class _ProgressLog(io.TextIOBase):
    """Collect what the generator prints, counting progress as it goes.

    convert_nucleotide() prints one "seed used:" line per sequence it
    optimises, which is a progress signal that costs nothing and leaves
    generate_pMHC_SCT.py untouched. Every other printed line (renamed
    columns, removed duplicates) is kept to show in the page.
    """

    def __init__(self, job: dict):
        self._job = job
        self._lines: list[str] = []
        self._partial = ""

    def write(self, text: str) -> int:
        self._partial += text
        while "\n" in self._partial:
            line, self._partial = self._partial.split("\n", 1)
            self._record(line)
        return len(text)

    def _record(self, line: str) -> None:
        if line.startswith("seed used:"):
            with _jobs_lock:
                self._job["done"] += 1
            return
        if line.strip():
            self._lines.append(line.rstrip())
            with _jobs_lock:
                self._job["log"] = "\n".join(self._lines)

    def flush(self) -> None:
        if self._partial:
            self._record(self._partial)
            self._partial = ""


def _preview(frame) -> dict:
    """First few rows of the result, for the table in the page."""
    head = frame.head(PREVIEW_ROWS)
    return {
        "columns": [str(column) for column in frame.columns],
        "rows": json.loads(head.to_json(orient="records")),
        "total_rows": int(len(frame)),
        "shown_rows": int(len(head)),
    }


def _run_job(job: dict, csv_text: str, seed: int, hla: str, leader: str,
             vector: str | None) -> None:
    log = _ProgressLog(job)
    try:
        with _generation_lock:
            with _jobs_lock:
                job["state"] = "running"
            with redirect_stdout(log):
                # Parsed up front so the peptide count is known before the
                # slow part starts; generate_sct_constructs takes the frame.
                peptides = sct.load_peptides(io.StringIO(csv_text))
                _check_peptides(peptides)
                with _jobs_lock:
                    # +1: the leader peptide is codon optimised once too
                    job["total"] = len(peptides) + 1
                frame = sct.generate_sct_constructs(
                    input_csv=peptides,
                    seed=seed,
                    hla=hla,
                    leader=leader,
                    vector=vector,
                    output_csv=None,
                )
            log.flush()
        with _jobs_lock:
            job["csv"] = frame.to_csv(index=False)
            job["preview"] = _preview(frame)
            job["done"] = job["total"]
            job["state"] = "done"
    except Exception as exc:  # surfaced in the page rather than the terminal
        log.flush()
        with _jobs_lock:
            job["state"] = "error"
            job["error"] = (str(exc) if isinstance(exc, InputError)
                            else f"{type(exc).__name__}: {exc}")
            job["traceback"] = (None if isinstance(exc, InputError)
                                else traceback.format_exc())


def _validate(payload: dict) -> tuple[dict | None, str | None]:
    """Check a generate request, returning (kwargs, error message)."""
    csv_text = payload.get("csv_text") or ""
    if not csv_text.strip():
        return None, "No CSV content received. Choose a peptide CSV first."

    seed_raw = payload.get("seed")
    if seed_raw is None or str(seed_raw).strip() == "":
        return None, ("A seed is required: the sequences cannot be reproduced "
                      "without one.")
    try:
        seed = int(str(seed_raw).strip())
    except ValueError:
        return None, f"Seed must be a whole number, got {str(seed_raw)!r}."
    # the range numpy's legacy global seeding accepts
    if not 0 <= seed <= 2**32 - 1:
        return None, "Seed must be between 0 and 4294967295."

    hla = payload.get("hla")
    if hla not in HLA_nuc:
        return None, f"Unknown HLA {hla!r}. Available: {sorted(HLA_nuc)}."

    leader = payload.get("leader")
    if leader not in leader_peptide_aa:
        return None, (f"Unknown leader {leader!r}. "
                      f"Available: {sorted(leader_peptide_aa)}.")

    vector = payload.get("vector") or None
    if vector is not None and vector not in destination_vector_nuc:
        return None, (f"Unknown vector {vector!r}. "
                      f"Available: {sorted(destination_vector_nuc)}.")

    return {
        "csv_text": csv_text,
        "seed": seed,
        "hla": hla,
        "leader": leader,
        "vector": vector,
    }, None


def _download_name(job: dict) -> str:
    parts = ["SCT_to_order", job["hla"], job["leader"]]
    if job["vector"]:
        parts.append(job["vector"])
    parts.append(f"seed{job['seed']}")
    stem = re.sub(r"[^A-Za-z0-9_.-]+", "-", "_".join(parts))
    return f"{stem}.csv"


class Handler(BaseHTTPRequestHandler):
    server_version = "pMHC-SCT-viewer"
    protocol_version = "HTTP/1.1"

    # -- plumbing ---------------------------------------------------------
    def _send(self, code: int, body: bytes, ctype: str,
              extra: dict[str, str] | None = None) -> None:
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        for key, value in (extra or {}).items():
            self.send_header(key, value)
        self.end_headers()
        self.wfile.write(body)

    def _json(self, code: int, payload: dict) -> None:
        self._send(code, json.dumps(payload).encode("utf-8"),
                   "application/json; charset=utf-8")

    def _job_from_query(self, query: dict) -> tuple[dict | None, str | None]:
        job_id = (query.get("job") or [""])[0]
        with _jobs_lock:
            job = _jobs.get(job_id)
            snapshot = dict(job) if job else None
        if snapshot is None:
            return None, f"Unknown job {job_id!r}."
        return snapshot, None

    def log_message(self, fmt: str, *args) -> None:
        if self.path.startswith("/api/status"):
            return  # polled twice a second; would bury everything else
        sys.stderr.write(f"  {self.command} {self.path}\n")

    # -- routes -----------------------------------------------------------
    def do_GET(self) -> None:
        parsed = urlparse(self.path)
        route = parsed.path
        query = parse_qs(parsed.query)

        if route in ("/", "/index.html"):
            try:
                body = PAGE.read_bytes()
            except OSError:
                self._json(500, {"error": f"Cannot read {PAGE}."})
                return
            self._send(200, body, "text/html; charset=utf-8")
            return

        if route == "/api/options":
            self._json(200, {
                "hla": sorted(HLA_nuc),
                "leaders": sorted(leader_peptide_aa),
                "vectors": sorted(destination_vector_nuc),
                "leader_sequences": dict(leader_peptide_aa),
                "defaults": {"hla": "HLA-A2", "leader": "B2M"},
                "preview_rows": PREVIEW_ROWS,
            })
            return

        if route == "/api/status":
            job, error = self._job_from_query(query)
            if error:
                self._json(404, {"error": error})
                return
            self._json(200, {
                "state": job["state"],
                "done": job["done"],
                "total": job["total"],
                "log": job["log"],
                "error": job.get("error"),
                "traceback": job.get("traceback"),
                "preview": job.get("preview"),
                "filename": _download_name(job),
            })
            return

        if route == "/api/result":
            job, error = self._job_from_query(query)
            if error:
                self._json(404, {"error": error})
                return
            if job["state"] != "done":
                self._json(409, {"error": "That run has not finished."})
                return
            disposition = f'attachment; filename="{_download_name(job)}"'
            self._send(200, job["csv"].encode("utf-8"),
                       "text/csv; charset=utf-8",
                       {"Content-Disposition": disposition})
            return

        self._json(404, {"error": f"No route {route!r}."})

    def do_POST(self) -> None:
        if urlparse(self.path).path != "/api/generate":
            self._json(404, {"error": f"No route {self.path!r}."})
            return

        length = int(self.headers.get("Content-Length") or 0)
        if length > MAX_UPLOAD_BYTES:
            self._json(413, {"error": "That CSV is larger than 8 MB."})
            return
        try:
            payload = json.loads(self.rfile.read(length) or b"{}")
        except (ValueError, UnicodeDecodeError) as exc:
            self._json(400, {"error": f"Could not read the request: {exc}"})
            return

        kwargs, error = _validate(payload)
        if error:
            self._json(400, {"error": error})
            return

        job_id = uuid.uuid4().hex
        job = {
            "state": "queued",
            "done": 0,
            "total": 0,
            "log": "",
            "seed": kwargs["seed"],
            "hla": kwargs["hla"],
            "leader": kwargs["leader"],
            "vector": kwargs["vector"],
        }
        with _jobs_lock:
            _jobs[job_id] = job
        threading.Thread(target=_run_job, args=(job,), kwargs=kwargs,
                         daemon=True).start()
        self._json(202, {"job": job_id, "busy": _generation_lock.locked()})


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--port", type=int, default=8765,
                        help="Port to serve on.")
    parser.add_argument("--no-browser", action="store_true",
                        help="Do not open a browser window on start.")
    args = parser.parse_args(argv)

    # 127.0.0.1 only: this is a local tool and the page runs the optimiser
    address = ("127.0.0.1", args.port)
    server = ThreadingHTTPServer(address, Handler)
    url = f"http://127.0.0.1:{args.port}"
    print(f"pMHC SCT viewer on {url}")
    print(f"  HLA:     {', '.join(sorted(HLA_nuc))}")
    print(f"  leaders: {', '.join(sorted(leader_peptide_aa))}")
    print(f"  vectors: {', '.join(sorted(destination_vector_nuc))}")
    print("  ctrl-c to stop")
    if not args.no_browser:
        threading.Timer(0.4, webbrowser.open, args=(url,)).start()
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")
    finally:
        server.server_close()


if __name__ == "__main__":
    main()
