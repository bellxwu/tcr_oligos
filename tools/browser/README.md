# pMHC SCT browser viewer

A local page for running `tools/generate_pMHC_SCT.py` without the command line:
upload a peptide CSV, pick the HLA / leader / vector, enter a seed, download the
constructs.

## Run it

```bash
conda activate tcr-oligos-env
python tools/browser/serve.py
```

That serves <http://127.0.0.1:8765> and opens it. `--port N` to move it,
`--no-browser` to not open a window, ctrl-c to stop.

No extra dependencies: the server is stdlib `http.server`, and the page is one
HTML file with no build step and no CDN.

## What it does and does not do

The page is a front end only. Every sequence comes from calling
`generate_sct_constructs()` in `tools/generate_pMHC_SCT.py`, so a CSV downloaded
here is byte for byte what the CLI writes for the same inputs:

```bash
# equivalent to choosing HLA-A2 / B2M / pHIV-EGFP / seed 1231 in the page
python tools/generate_pMHC_SCT.py --input-csv peptides.csv --seed 1231 \
    --hla HLA-A2 --leader B2M --vector pHIV-EGFP --output-csv out.csv
```

Nothing about the construct is reimplemented in JavaScript, and that is
deliberate: the codon optimisation is dnachisel driving numpy's global RNG, and
that seeded RNG is the whole reason a run can be reproduced. A JavaScript
reimplementation would quietly return different sequences for the same seed.

Two consequences:

- **It has to be served, not opened as a file.** Double-clicking `index.html`
  gives a page that cannot reach the optimiser.
- **Only one run executes at a time.** numpy's RNG is process-wide, so
  overlapping runs would reseed each other and produce sequences that neither
  seed reproduces. A second request waits, and the page says it is queued.

## Choices and seeds

The dropdowns are read from `tools/SCT_constants.py` at page load, so adding an
allele to `HLA_nuc`, a leader to `leader_peptide_aa` or a vector to
`destination_vector_nuc` makes it appear here on reload, with no change to this
directory.

The seed is typed in by hand, as on the CLI, and is written into the output's
`seed` column. To rebuild a past construct, re-enter its seed with the same
peptide and the same choices.

## Files

- `serve.py` - localhost server: options, job running, progress, CSV download.
  Also pre-checks peptides so a stray residue names the offending row instead
  of failing as `KeyError: '3'` from inside dnachisel.
- `index.html` - the page.
