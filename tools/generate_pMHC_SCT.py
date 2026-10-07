"""
Build pMHC single-chain trimer (SCT) nucleotide sequences from a peptide CSV.

Takes a CSV of peptide amino acid sequences, codon optimises each one, and
concatenates it into the full SCT construct ready for ordering:

    5'addon - signal peptide - peptide - L1 - B2M - L2 - HLA heavy chain - 3'addon

Pulls variables from SCT constants. Changing this file will change the sequences used.


Input CSV must have a "name" column and a "peptide_aa" column:
    name,peptide_aa
    MART-1,AAGIGILTV
    NY-ESO-1,SLLMWITQC

Example for CLI use:
    python analyses/06_HLA_test.py \
        --input-csv all_data/peptides.csv \
        --seed 1231 \
        --output-csv all_data/06_SCT_to_order.csv
"""

import numpy as np
from dnachisel import (DnaOptimizationProblem, CodonOptimize,
                       EnforceTranslation, AvoidPattern, EnforceGCContent)
from dnachisel.biotools import reverse_translate, translate
import argparse
import sys
from pathlib import Path
import pandas as pd
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import ALL_DATA_DIR
from SCT_constants import (B2M_nuc, HLA_aa, destination_vector_nuc,
                           leader_peptide_aa, linker_nuc)
import contextlib
import io
from Bio.SeqUtils import gc_fraction


def require_vector(vector):
    """
    Return the Gibson arms for a vector, refusing a missing or unknown one.

    params:
        vector: key into destination_vector_nuc. None, an empty string and
                unrecognised keys all raise instead of quietly giving back an
                unflanked construct.

    return:
        That vector's {"5_addon": ..., "3_addon": ...} mapping.
    """
    if vector is None or str(vector).strip() == "":
        raise ValueError(
            "a vector is required: without one the construct has no Gibson "
            "homology arms and cannot be assembled into a destination vector, "
            "so pass an explicit vector. available: "
            f"{sorted(destination_vector_nuc)}"
        )
    if vector not in destination_vector_nuc:
        raise KeyError(
            f"unknown vector {vector!r}; available: {sorted(destination_vector_nuc)}"
        )
    return destination_vector_nuc[vector]


def convert_nucleotide(
        aa_seq, seed, species="h_sapiens",
        avoid_patterns=None,
        gc_mini=0.35, gc_maxi=0.65, gc_window=50):
    """
    Codon-optimise an amino acid sequence into a nucleotide sequence.
    dnachisel draws all of its randomness from numpy's global RNG (both the
    randomised back-translation and the optimisation itself), so seeding that
    RNG is what makes a run reproducible. The seed is therefore required: no
    sequence is generated without one, so every sequence on record can be
    regenerated later. The seed used is printed at the end.

    Parameters
    ----------
    aa_seq
        Amino acid sequence to reverse translate.
    seed
        Required integer seed. Passing None raises a ValueError.
    species
        Species for the CodonOptimize objective.
    avoid_patterns
        Optional patterns to forbid. Accepts a single pattern ("BsaI_site")
        or a list of them. If None, no patterns are avoided.
    gc_mini, gc_maxi, gc_window
        Local GC content constraint, applied within each sliding window.

    Returns
    -------
    The optimised nucleotide sequence as a string.
    """
    if seed is None:
        raise ValueError(
            "a seed is required: the sequence cannot be reproduced without "
            "one, so pass an explicit seed (e.g. seed=1234) before running. "
            "Draw one with int(np.random.randint(0, 1_000_000)) and record it."
        )

    aa_seq = "".join(aa_seq.split()).upper()

    if avoid_patterns is None:
        avoid_patterns = []
    elif isinstance(avoid_patterns, str):
        avoid_patterns = [avoid_patterns]

    np.random.seed(seed)

    problem = DnaOptimizationProblem(
        sequence=reverse_translate(aa_seq, randomize_codons=True),
        constraints=[EnforceTranslation(),
                     *[AvoidPattern(pattern) for pattern in avoid_patterns],
                     EnforceGCContent(mini=gc_mini, maxi=gc_maxi,
                                      window=gc_window)],
        objectives=[CodonOptimize(species=species)],
        logger=None,
    )
    problem.resolve_constraints()
    problem.optimize()

    dna = problem.sequence
    assert translate(dna) == aa_seq, "translation does not match input protein"

    print(f"seed used: {seed}")
    return dna


def load_peptides(input_csv):
    """
    Read a peptide CSV into a dataframe of name / peptide_aa pairs.

    params:
        input_csv: path to a CSV with "name" and "peptide_aa" columns.

    return:
        Dataframe with name and peptide_aa columns, whitespace stripped and
        amino acid sequences upper-cased.
    """
    peptides = input_csv.copy() if isinstance(input_csv, pd.DataFrame) else pd.read_csv(input_csv)

    # strip whitespace/BOM from names
    peptides.columns = peptides.columns.str.replace("﻿", "", regex=False).str.strip()

    # first two columns are taken as name / peptide_aa, whatever they are called
    renames = {
        old: new
        for old, new in zip(peptides.columns[:2], ["name", "peptide_aa"])
        if old != new
    }
    if renames:
        peptides = peptides.rename(columns=renames)
        for old, new in renames.items():
            print(f"renamed column {old!r} -> {new!r}")

    # strips white space from values inside columns
    peptides = peptides[["name", "peptide_aa"]].copy()
    peptides["name"] = peptides["name"].astype(str).str.strip()
    peptides["peptide_aa"] = (
        peptides["peptide_aa"].astype(str).str.replace(r"\s+", "", regex=True).str.upper()
    )
    # removes empty values
    peptides = peptides[peptides["peptide_aa"] != ""]
    
    # duplicates are dropped rather than raised, keeping the first occurrence

    # dup_names = peptides.loc[peptides["name"].duplicated(), "name"]
    # for name in dup_names:
    #     print(f"removed duplicate name: {name!r}")
    # peptides = peptides[~peptides["name"].duplicated()]

    # checked after names, so only rows that survived above are considered
    dup_seqs = peptides[peptides["peptide_aa"].duplicated()]
    for name, seq in zip(dup_seqs["name"], dup_seqs["peptide_aa"]):
        print(f"removed duplicate sequence: {name!r} ({seq})")
    peptides = peptides[~peptides["peptide_aa"].duplicated()]

    return peptides.reset_index(drop=True)


def build_sct(peptide_nuc, signal_peptide_nuc, vector, hla="HLA-A2"):
    """
    Join the SCT components into a single nucleotide sequence.

    Chain order is 5'addon - signal peptide - peptide - L1 - B2M - L2 - HLA -
    3'addon, matching SCT_CHAIN_ORDER.

    params:
        peptide_nuc: codon optimised nucleotide sequence for the peptide.
        signal_peptide_nuc: codon optimised nucleotide sequence for the leader.
        vector: required key into destination_vector_nuc; its homology arms
                flank the construct.
        hla: key into HLA_aa naming the heavy chain to use.

    return:
        The assembled nucleotide sequence as a string.
    """
    addons = require_vector(vector)

    if hla not in HLA_nuc:
        raise KeyError(f"unknown HLA {hla!r}; available: {sorted(HLA_nuc)}")
    
    # join the sequences together
    return "".join([
        addons["5_addon"],
        signal_peptide_nuc,
        peptide_nuc,
        linker_nuc["L1"],
        B2M_nuc,
        linker_nuc["L2"],
        HLA_nuc[hla],
        addons["3_addon"],
    ])


def generate_sct_constructs(input_csv, seed, vector, hla="HLA-A2",
                            leader="B2M", output_csv=None):
    """
    Convert a peptide CSV into full SCT constructs ready for ordering.

    Each peptide is codon optimised with the same seed, so any single
    construct can be regenerated by rerunning this with the same seed, the
    same peptide and the same vector.

    params:
        input_csv: path to the peptide CSV, or a dataframe (see load_peptides).
        seed: required integer seed for the codon optimisation.
        vector: required destination vector whose homology arms are added.
        hla: key into HLA_nuc naming the heavy chain to use.
        leader: key into leader_peptide_aa naming the signal peptide.
        output_csv: optionally, a path to write the result to.

    return:
        Dataframe with name, peptide_aa, peptide_nuc and to_order columns.
    """
    # all three are dictionary lookups, so they go before the slow part: a bad
    # vector or leader should not cost a full codon optimisation to discover
    require_vector(vector)
    if leader not in leader_peptide_aa:
        raise KeyError(
            f"unknown leader {leader!r}; available: {sorted(leader_peptide_aa)}"
        )
    if hla not in HLA_nuc:
        raise KeyError(f"unknown HLA {hla!r}; available: {sorted(HLA_nuc)}")
    # load peptide csv
    peptides = load_peptides(input_csv)
    if peptides.empty:
        raise ValueError(f"no peptide sequences found in {input_csv}")
    # convert the amino acid of the leader peptide into nucleotide
    signal_peptide_nuc = convert_nucleotide(leader_peptide_aa[leader], seed=seed)
    
    # conver peptides nucleotides into nucleotides
    peptide_nucs = []
    for peptide_aa in peptides["peptide_aa"]:
        # use entire peptide sequence as window
        gc_window = len(peptide_aa) * 3
        peptide_nucs.append(
            convert_nucleotide(peptide_aa, seed=seed, gc_window=gc_window)
        )
    # build the sct by concatating the sequence length together
    peptides["peptide_nuc"] = peptide_nucs
    peptides["to_order"] = [
        build_sct(peptide_nuc, signal_peptide_nuc, vector, hla=hla)
        for peptide_nuc in peptides["peptide_nuc"]
    ]

    # add metadata columns
    peptides["hla"] = hla
    peptides["leader"] = leader
    peptides["vector"] = vector
    peptides["seed"] = seed

    # add GC fraction
    peptides["GC_fraction"] = [
        round(gc_fraction(seq), 3) for seq in peptides["to_order"]
    ]

    if output_csv is not None:
        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        peptides.to_csv(output_csv, index=False)
        print(f"wrote {len(peptides)} constructs to {output_csv}")

    return peptides

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Build pMHC single-chain trimer sequences from a peptide CSV.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input-csv", required=True,
                        help="CSV with 'name' and 'peptide_aa' columns.")
    parser.add_argument("--seed", type=int, required=True,
                        help="Seed for the codon optimisation; record it to reproduce the run.")
    parser.add_argument("--vector", required=True, choices=sorted(destination_vector_nuc),
                        help="Destination vector whose Gibson arms flank the construct.")
    parser.add_argument("--hla", default="HLA-A2", choices=sorted(HLA_nuc),
                        help="HLA heavy chain to build into the trimer.")
    parser.add_argument("--leader", default="B2M", choices=sorted(leader_peptide_aa),
                        help="Signal peptide to put at the N terminus.")
    parser.add_argument("--output-csv", default=str(ALL_DATA_DIR / "06_SCT_to_order.csv"),
                        help="Where to write the constructs.")
    return parser.parse_args(argv)


# every one of these is refused before a single codon is optimised
for bad in (None, "", "   ", "pHIV_EGFP", "pUC19"):
    try:
        require_vector(bad)
    except (ValueError, KeyError) as exc:
        print(f"{bad!r:>12} -> {type(exc).__name__}")
    else:
        print(f"{bad!r:>12} -> ACCEPTED, which it should not be")

print(f"{'pHIV-EGFP'!r:>12} -> accepted, arms {sorted(require_vector('pHIV-EGFP'))}")


# and the CLI will not parse without one; argparse prints its own usage to
# stderr on the way out, which is swallowed here to keep the output in order
with contextlib.redirect_stderr(io.StringIO()):
    try:
        parse_args(["--input-csv", "x.csv", "--seed", "1"])
        verdict = "argparse accepted it, which it should not"
    except SystemExit:
        verdict = "argparse exits when it is missing"
print(f"{'--vector':>12} -> {verdict}")

if __name__ == "__main__":
    args = parse_args()
    generate_sct_constructs(
        input_csv=args.input_csv,
        seed=args.seed,
        hla=args.hla,
        leader=args.leader,
        vector=args.vector,
        output_csv=args.output_csv,
    )
