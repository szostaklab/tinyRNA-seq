#!/usr/bin/env python3
"""
dno_filter.py

Minimal filter for DNO (non-VCG) experiments: drop reads with a zero-length
insert (VCG_Len == 0). No sequence matching against a reference list is
performed — all non-empty inserts pass through.

Outputs:
    vcg_filtered.tsv  — reads with a non-empty insert
    vcg_excluded.tsv  — dropped rows with an 'exclude_reason' column

Called by Snakemake (rule dno_filter) or directly:
    python scripts/dno_filter.py \\
        --parsed parsed_reads.tsv \\
        --filtered vcg_filtered.tsv \\
        --excluded vcg_excluded.tsv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def dno_filter(
    parsed_path: str,
    filtered_path: str,
    excluded_path: str,
    sample: str = "",
) -> None:
    Path(filtered_path).parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(parsed_path, sep="\t", dtype=str)
    df["VCG_Len"] = pd.to_numeric(df["VCG_Len"], errors="coerce").fillna(0).astype(int)

    n_total = len(df)

    zero_mask = df["VCG_Len"] == 0
    excluded = df[zero_mask].copy()
    excluded["exclude_reason"] = "insert_len_zero"
    df = df[~zero_mask]

    n_after = len(df)
    print(f"[dno_filter:{sample}] Dropped {n_total - n_after:,} rows with zero-length insert "
          f"({(n_total - n_after) / max(n_total, 1):.2%})")
    print(f"[dno_filter:{sample}] Remaining: {n_after:,} reads")

    df.to_csv(filtered_path, sep="\t", index=False)
    excluded.to_csv(excluded_path, sep="\t", index=False)

    print(f"[dno_filter:{sample}] Filtered  -> {filtered_path}")
    print(f"[dno_filter:{sample}] Excluded  -> {excluded_path}")


if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    dno_filter(
        parsed_path   = snakemake.input.parsed,
        filtered_path = snakemake.output.filtered,
        excluded_path = snakemake.output.excluded,
        sample        = snakemake.wildcards.sample,
    )

elif __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Filter DNO reads: drop zero-length inserts.")
    ap.add_argument("--parsed",   required=True)
    ap.add_argument("--filtered", required=True)
    ap.add_argument("--excluded", required=True)
    args = ap.parse_args()
    dno_filter(args.parsed, args.filtered, args.excluded)
