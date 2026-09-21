#!/usr/bin/env python
"""
Reusable script version of analyses/02_Standardize_TCR.ipynb.

Takes a TCR list (in AS_top_96 or positive_ctrls format), standardizes it to
the TCRAFT input format, generates + validates CDR3 oligos, and writes out
oligo-pool CSVs ready for ordering (e.g. from IDT).

Example:
    python analyses/generate_oligopool.py \
        --input-csv all_data/AS_Top_96_TCRs.csv \
        --prefix 01 \
        --name AS_TCRAFT_96 \
        --select 93
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import ALL_DATA_DIR

from TCRAFT.Generate import generate_cdr3_oligos
from TCRAFT.Validate import validate_cdr3_oligos


def standardize_to_sample_format(df, select=None, output_dir=None, filename="standardized_TCRs.csv"):
    """Convert a df in AS_top_96 or positive_ctrls format into sample_df format.

    select: optionally, an integer number of TCRs to randomly select from the result.
    output_dir: optionally, a directory to save the standardized df to as a CSV.
                The directory is created if it doesn't already exist.
    filename: name of the CSV file to write into output_dir (ignored if output_dir is None).
    """
    df = df.copy()
    df.columns = df.columns.str.replace("﻿", "", regex=False).str.strip()

    # AS_top_96 uses short lowercase column names; rename to sample_df's convention.
    # positive_ctrls already uses sample_df's column names (plus extra columns).
    column_map = {
        "CLONE_ID": "Name",
        "name": "Name",
        "va": "V_alpha",
        "vb": "V_beta",
        "cdr3a": "CDR3_alpha",
        "cdr3b": "CDR3_beta",
        "ja": "J_alpha",
        "jb": "J_beta",
    }
    df = df.rename(columns=column_map)

    # A concatenated input can carry both conventions at once (e.g. va and V_alpha,
    # or CLONE_ID and Name), so renaming collapses them onto duplicate labels. Merge
    # each set of duplicates into a single column, keeping the first non-null value
    # per row. Left alone, the duplicates would also make df[sample_cols] below
    # return every matching column rather than one per name.
    # if df.columns.duplicated().any():
    #     df = pd.concat(
    #         [df.loc[:, df.columns == col].bfill(axis=1).iloc[:, 0].rename(col)
    #          for col in df.columns.unique()],
    #         axis=1,
    #     )

    sample_cols = ["Name", "V_alpha", "V_beta", "CDR3_alpha", "CDR3_beta", "J_alpha", "J_beta"]
    df = df[sample_cols].dropna(how="all").reset_index(drop=True)

    if select is not None:
        df = df.sample(n=select).reset_index(drop=True)

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_dir / filename, index=False)

    return df


def generate_oligos(tcraft_csv, prefix, name, data_dir=ALL_DATA_DIR):
    """Generate CDR3 oligos from a TCRAFT-formatted TCR list CSV.

    params:
        tcraft_csv: path to a CSV with columns V_alpha, V_beta, CDR3_alpha,
                    CDR3_beta, J_alpha, J_beta (e.g. from standardize_to_sample_format).
        prefix: run identifier used to namespace the output directory (e.g. "01").
        name: descriptive name for this TCR set (e.g. "AS_TCRAFT_96").
        data_dir: directory to write the output directory into.

    return:
        Path to the generated CDR3 oligos CSV (All_<N>_CDR3_oligos.csv).
    """
    data_dir = Path(data_dir)

    oligos_dir = data_dir / f"{prefix}_oligos_{name}"
    generate_cdr3_oligos(
        input_csv_path=str(tcraft_csv), 
        output_dir=str(oligos_dir)
    )

    oligos_csv_matches = sorted(oligos_dir.glob("All_*_CDR3_oligos.csv"))
    if not oligos_csv_matches:
        raise FileNotFoundError(f"No generated oligos CSV found in {oligos_dir}")

    return oligos_csv_matches[-1]


def validate_oligos(oligos_csv, prefix, name, data_dir=ALL_DATA_DIR):
    """Validate the assembly of generated CDR3 oligos by simulating Golden Gate.

    params:
        oligos_csv: path to the CDR3 oligos CSV produced by generate_oligos.
        prefix: run identifier used to namespace the output directory (e.g. "01").
        name: descriptive name for this TCR set (e.g. "AS_TCRAFT_96").
        data_dir: directory to write the output directory into.

    return:
        Path to the validation output directory.
    """
    data_dir = Path(data_dir)

    validation_dir = data_dir / f"{prefix}_validation_{name}"
    validate_cdr3_oligos(input_csv_path=str(oligos_csv), output_dir=str(validation_dir))

    return validation_dir


def select_sequences_by_pool(oligos_csv, prefix, data_dir=ALL_DATA_DIR):
    """Split a generated CDR3 oligos CSV into per-pool sequence CSVs ready for ordering.

    params:
        oligos_csv: path to the CDR3 oligos CSV produced by generate_oligos.
        prefix: run identifier used to namespace output files (e.g. "01").
        data_dir: directory to write output files into.

    return:
        dict mapping pool label ("all", "A", "B", ...) to the CSV path written for it.
    """
    data_dir = Path(data_dir)
    cdr3_oligos = pd.read_csv(oligos_csv)

    oligopool_paths = {}

    combined_path = data_dir / f"{prefix}_IDT_oligos.csv"
    cdr3_oligos[["Name", "Sequence"]].to_csv(combined_path, index=False)
    oligopool_paths["all"] = combined_path

    for pool_label, pool_df in cdr3_oligos.groupby("Pool"):
        pool_path = data_dir / f"{prefix}_Pool{pool_label}_IDT_oligos.csv"
        pool_df[["Name", "Sequence"]].to_csv(pool_path, index=False)
        oligopool_paths[pool_label] = pool_path

    return oligopool_paths