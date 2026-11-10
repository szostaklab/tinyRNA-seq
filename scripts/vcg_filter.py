#!/usr/bin/env python3
"""
vcg_filter.py

Filter parsed reads:
  1. Drop rows where VCG_Len == 0 (empty VCG after trimming)
  2. Drop rows where VCG is not in the reference list (VCG_Oligos_Labelled.csv)

The FASTQ read_id column is preserved throughout.

Outputs:
    vcg_filtered.tsv  — reads with a valid, non-empty VCG in the reference set
    vcg_excluded.tsv  — dropped rows with an 'exclude_reason' column

Called by Snakemake (rule vcg_filter) or directly:
    python scripts/vcg_filter.py \\
        --parsed parsed_reads.tsv \\
        --vcg-list references/synthetic/VCG_Oligos_Labelled.csv \\
        --filtered vcg_filtered.tsv \\
        --excluded vcg_excluded.tsv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


# ── VCG list loading ──────────────────────────────────────────────────────────

_NUCS = set("ACGTU")


def load_vcg_list(vcg_list_path: str) -> set[str]:
    """
    Load VCG reference sequences in DNA form (to match the DNA-mode parsed reads).

    Accepts either the labelled CSV (VCG_Oligos_Labelled.csv: VCG_Code,sequence,...
    where the sequence is column 2) or a plain one-per-line list (optional leading
    'p'). Handles CRLF, converts RNA→DNA (U→T), and skips header/non-sequence rows.
    """
    seqs = set()
    with open(vcg_list_path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            token = line.split(",")[1].strip() if "," in line else line.lstrip("pP")
            token = token.upper().replace("U", "T")     # RNA→DNA, normalise case
            if token and not (set(token) - _NUCS):      # skip headers / junk rows
                seqs.add(token)
    return seqs


# ── Core ──────────────────────────────────────────────────────────────────────

def vcg_filter(
    parsed_path: str,
    vcg_list_path: str,
    filtered_path: str,
    excluded_path: str,
    sample: str = "",
) -> None:
    Path(filtered_path).parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(parsed_path, sep="\t", dtype=str)
    df["VCG_Len"] = pd.to_numeric(df["VCG_Len"], errors="coerce").fillna(0).astype(int)

    n_total = len(df)

    # ── Step 1: drop zero-length VCG ─────────────────────────────────────────
    zero_mask = df["VCG_Len"] == 0
    excluded  = df[zero_mask].copy()
    excluded["exclude_reason"] = "vcg_len_zero"
    df = df[~zero_mask]

    n_after_len = len(df)
    print(f"[vcg_filter:{sample}] Dropped {n_total - n_after_len:,} rows with VCG_Len == 0 "
          f"({(n_total - n_after_len) / max(n_total, 1):.2%})")

    # ── Step 2: drop VCGs not in reference list ───────────────────────────────
    vcg_set  = load_vcg_list(vcg_list_path)
    not_in   = ~df["VCG"].isin(vcg_set)
    excluded2 = df[not_in].copy()
    excluded2["exclude_reason"] = "vcg_not_in_reference"
    df = df[~not_in]

    n_after_ref = len(df)
    print(f"[vcg_filter:{sample}] Dropped {n_after_len - n_after_ref:,} rows with VCG not in reference "
          f"({(n_after_len - n_after_ref) / max(n_after_len, 1):.2%})")
    print(f"[vcg_filter:{sample}] Remaining: {n_after_ref:,} reads")
    print(f"[vcg_filter:{sample}] Reference VCG list: {len(vcg_set):,} sequences")

    # ── Write outputs ─────────────────────────────────────────────────────────
    all_excluded = pd.concat([excluded, excluded2], ignore_index=True)

    df.to_csv(filtered_path, sep="\t", index=False)
    all_excluded.to_csv(excluded_path, sep="\t", index=False)

    print(f"[vcg_filter:{sample}] Filtered → {filtered_path}")
    print(f"[vcg_filter:{sample}] Excluded → {excluded_path}")


# ── Entry points ──────────────────────────────────────────────────────────────

if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    vcg_filter(
        parsed_path   = snakemake.input.parsed,
        vcg_list_path = snakemake.input.vcg_list,
        filtered_path = snakemake.output.filtered,
        excluded_path = snakemake.output.excluded,
        sample        = snakemake.wildcards.sample,
    )

elif __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Filter parsed reads by VCG length and reference list.")
    ap.add_argument("--parsed",   required=True)
    ap.add_argument("--vcg-list", required=True, dest="vcg_list")
    ap.add_argument("--filtered", required=True)
    ap.add_argument("--excluded", required=True)
    args = ap.parse_args()
    vcg_filter(args.parsed, args.vcg_list, args.filtered, args.excluded)
