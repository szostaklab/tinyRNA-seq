#!/usr/bin/env python3
"""
normalise_dno.py

CPM normalisation for DNO (non-VCG) experiments.

Unlike normalise.py, there is no reference sequence list. Norm_Reads is
counts per million of total unique molecules observed in this sample.
No zero-filling: only inserts that were actually observed appear in the output.

Output columns:
    VCG_Len, VCG (RNA), VCG_DNA, Count, Norm_Reads

Called by Snakemake (rule normalise_dno) or directly:
    python scripts/normalise_dno.py \\
        --collapsed umi_collapsed.tsv \\
        --normalised normalised.tsv \\
        --csv sample_Norm.csv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def normalise_dno(
    collapsed_path: str,
    normalised_path: str,
    sample: str = "",
    csv_path: str | None = None,
) -> None:
    Path(normalised_path).parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(collapsed_path, sep="\t", dtype=str)

    if "VCG" not in df.columns or "VCG_Len" not in df.columns:
        raise ValueError(f"Expected 'VCG' and 'VCG_Len' columns in {collapsed_path}.")

    df["VCG_Len"] = pd.to_numeric(df["VCG_Len"], errors="coerce").astype(int)

    grouped = (
        df.groupby(["VCG_Len", "VCG"], as_index=False)
        .size()
        .rename(columns={"size": "Count"})
    )

    total_counts = grouped["Count"].sum()
    grouped["Norm_Reads"] = (grouped["Count"] / total_counts * 1_000_000).round(4)

    print(f"[normalise_dno:{sample}] Total collapsed reads:    {total_counts:,}")
    print(f"[normalise_dno:{sample}] Unique insert sequences:  {len(grouped):,}")

    grouped["VCG_DNA"] = grouped["VCG"]
    grouped["VCG"]     = grouped["VCG"].str.replace("T", "U", regex=False)

    out = grouped[["VCG_Len", "VCG", "VCG_DNA", "Count", "Norm_Reads"]]
    out = out.sort_values(["VCG_Len", "VCG"], ascending=[True, True])
    out.to_csv(normalised_path, sep="\t", index=False)
    print(f"[normalise_dno:{sample}] Normalised -> {normalised_path}")

    if csv_path:
        out.to_csv(csv_path, index=False)
        print(f"[normalise_dno:{sample}] Normalised -> {csv_path}")


if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    normalise_dno(
        collapsed_path  = snakemake.input.collapsed,
        normalised_path = snakemake.output.normalised,
        sample          = snakemake.wildcards.sample,
        csv_path        = snakemake.output.csv,
    )

elif __name__ == "__main__":
    ap = argparse.ArgumentParser(description="CPM normalise DNO reads.")
    ap.add_argument("--collapsed",  required=True)
    ap.add_argument("--normalised", required=True)
    ap.add_argument("--csv",        default=None)
    args = ap.parse_args()
    normalise_dno(args.collapsed, args.normalised, csv_path=args.csv)
