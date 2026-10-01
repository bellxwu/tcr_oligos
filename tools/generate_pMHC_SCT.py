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
from SCT_constants import (B2M_nuc, HLA_nuc, destination_vector_nuc,
                           leader_peptide_aa, linker_nuc)
from Bio.SeqUtils import gc_fraction



def convert_nucleotide(
        aa_seq, seed, species="h_sapiens",
        avoid_patterns=None,
        gc_mini=0.35, gc_maxi=0.65, gc_window=50):
    """Codon-optimise an amino acid sequence into a nucleotide sequence.

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


def build_sct(peptide_nuc, signal_peptide_nuc, hla="HLA-A2", vector=None):
    """
    Join the SCT components into a single nucleotide sequence.

    Chain order is signal peptide - peptide - L1 - B2M - L2 - HLA.

    params:
        peptide_nuc: codon optimised nucleotide sequence for the peptide.
        signal_peptide_nuc: codon optimised nucleotide sequence for the leader.
        hla: key into HLA_nuc naming the heavy chain to use.
        vector: optionally, a key into destination_vector_nuc. If given, the
                construct is flanked by that vector's Gibson homology arms.

    return:
        The assembled nucleotide sequence as a string.
    """
    if hla not in HLA_nuc:
        raise KeyError(f"unknown HLA {hla!r}; available: {sorted(HLA_nuc)}")

    sequence = "".join([
        signal_peptide_nuc,
        peptide_nuc,
        linker_nuc["L1"],
        B2M_nuc,
        linker_nuc["L2"],
        HLA_nuc[hla],
    ])

    if vector is not None:
        if vector not in destination_vector_nuc:
            raise KeyError(
                f"unknown vector {vector!r}; available: {sorted(destination_vector_nuc)}"
            )
        addons = destination_vector_nuc[vector]
        sequence = addons["5_addon"] + sequence + addons["3_addon"]

    return sequence


def generate_sct_constructs(input_csv, seed, hla="HLA-A2", leader="B2M",
                            vector=None, output_csv=None):
    """
    Convert a peptide CSV into full SCT constructs ready for ordering.

    Each peptide is codon optimised with the same seed, so any single
    construct can be regenerated by rerunning this with the same seed and the
    same peptide.

    params:
        input_csv: path to the peptide CSV (see load_peptides).
        seed: required integer seed for the codon optimisation.
        hla: key into HLA_nuc naming the heavy chain to use.
        leader: key into leader_peptide_aa naming the signal peptide.
        vector: optionally, a destination vector whose homology arms are added.
        output_csv: optionally, a path to write the result to.

    return:
        Dataframe with name, peptide_aa, peptide_nuc and to_order columns.
    """
    # load leader peptide
    if leader not in leader_peptide_aa:
        raise KeyError(
            f"unknown leader {leader!r}; available: {sorted(leader_peptide_aa)}"
        )
    # load peptide
    peptides = load_peptides(input_csv)
    if peptides.empty:
        raise ValueError(f"no peptide sequences found in {input_csv}")
    # convert leader peptide aa into nucleotide sequence
    signal_peptide_nuc = convert_nucleotide(leader_peptide_aa[leader], seed=seed)

    # convert peptide aa list into nucleotide sequence append to empty list
    peptide_nucs = []
    for peptide_aa in peptides["peptide_aa"]:
        # use entire peptide sequence as window
        gc_window = len(peptide_aa)*3
        peptide_nucs.append(
            convert_nucleotide(peptide_aa, seed=seed, gc_window=gc_window)
        )
    
    # add nucleotide to original loaded dataframe
    peptides["peptide_nuc"] = peptide_nucs
    # concat all nucleotide sequnces together
    peptides["to_order"] = [
        build_sct(peptide_nuc, signal_peptide_nuc, hla=hla, vector=vector)
        for peptide_nuc in peptides["peptide_nuc"]
    ]

    # add metadata
    peptides["hla"] = hla
    peptides["leader"] = leader
    peptides["vector"] = vector
    peptides["seed"] = seed

    # add GC fraction
    peptides["GC_fraction"] = [round(gc_fraction(seq), 3) for seq in peptides["to_order"]]

    if output_csv is not None:
        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        peptides.to_csv(output_csv, index=False)
        print(f"wrote {len(peptides)} constructs to {output_csv}")

    return peptides


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-csv", required=True,
                        help="CSV with 'name' and 'peptide_aa' columns.")
    parser.add_argument("--seed", type=int, required=True,
                        help="Seed for the codon optimisation; record it to reproduce the run.")
    parser.add_argument("--hla", default="HLA-A2", choices=sorted(HLA_nuc),
                        help="HLA heavy chain to build into the trimer.")
    parser.add_argument("--leader", default="B2M", choices=sorted(leader_peptide_aa),
                        help="Signal peptide to put at the N terminus.")
    parser.add_argument("--vector", default=None, choices=sorted(destination_vector_nuc),
                        help="Optionally flank the construct with this vector's Gibson arms.")
    parser.add_argument("--output-csv", default=str(ALL_DATA_DIR / "06_SCT_to_order.csv"),
                        help="Where to write the constructs.")
    return parser.parse_args(argv)


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
