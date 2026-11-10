#!/usr/bin/env python3
"""
normalise.py

Group UMI-collapsed reads by (VCG_Len, VCG), count unique reads per VCG,
and compute normalised read counts.

Normalisation:
    Norm_Reads = count / (total_counts / n_reference_seqs)

Where n_reference_seqs is the number of UNIQUE sequences in the VCG reference
list, derived at runtime. This means results are comparable across experiments
that use different sized reference pools.

The VCG column is converted from DNA to RNA (T→U) in the output.
All upstream processing uses DNA notation.

Output columns:
    VCG_Len, VCG (RNA), VCG_DNA, Count, Norm_Reads

Called by Snakemake (rule normalise) or directly:
    python scripts/normalise.py \\
        --collapsed umi_collapsed.tsv \\
        --vcg-list references/synthetic/VCG_Oligos_Labelled.csv \\
        --normalised normalised.tsv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


# ── VCG reference ──────────────────────────────────────────────────────────────

_NUCS = set("ACGTU")


def _load_vcg_reference(path: str) -> set[str]:
    """Return the set of VCG reference sequences in DNA form.

    Accepts the labelled CSV (VCG_Oligos_Labelled.csv: code,sequence,... → column 2)
    or a plain one-per-line list (optional leading 'p'). Handles CRLF line endings,
    converts U→T, and skips header/non-sequence rows.
    """
    seqs: set[str] = set()
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            token = line.split(",")[1].strip() if "," in line else line.lstrip("pP")
            token = token.upper().replace("U", "T")
            if token and not (set(token) - _NUCS):   # skip headers / junk rows
                seqs.add(token)
    return seqs


def load_unique_vcg_set(vcg_list_path: str) -> set[str]:
    """
    Load the unique set of VCG reference sequences in DNA form.

    Accepts either the labelled CSV (VCG_Oligos_Labelled.csv: VCG_Code,sequence,...
    where the sequence is column 2) or a plain one-per-line list (optional leading
    'p'). Sequences are back-transcribed RNA→DNA (U→T) to match upstream notation.
    """
    return _load_vcg_reference(vcg_list_path)


# ── Core ──────────────────────────────────────────────────────────────────────

def normalise(
    collapsed_path: str,
    vcg_list_path: str,
    normalised_path: str,
    sample: str = "",
    csv_path: str | None = None,
) -> None:
    Path(normalised_path).parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(collapsed_path, sep="\t", dtype=str)

    if "VCG" not in df.columns or "VCG_Len" not in df.columns:
        raise ValueError(f"Expected 'VCG' and 'VCG_Len' columns in {collapsed_path}.")

    df["VCG_Len"] = pd.to_numeric(df["VCG_Len"], errors="coerce").astype(int)

    # ── Group and count ───────────────────────────────────────────────────────
    # Each row in umi_collapsed.tsv is one unique UMI cluster (unique molecule).
    # Count rows per VCG to get unique-molecule counts.
    # Works for both UMI-deduped layouts (Count column present) and pass-through
    # layouts (no Count column) — in both cases we want number of unique reads.
    grouped = (
        df.groupby(["VCG_Len", "VCG"], as_index=False)
        .size()
        .rename(columns={"size": "Count"})
    )

    total_counts = grouped["Count"].sum()
    vcg_set      = load_unique_vcg_set(vcg_list_path)
    n_ref_seqs   = len(vcg_set)
    norm_factor  = total_counts / n_ref_seqs

    grouped["Norm_Reads"] = (grouped["Count"] / norm_factor).round(4)

    print(f"[normalise:{sample}] Total collapsed reads:     {total_counts:,}")
    print(f"[normalise:{sample}] Unique reference VCG seqs: {n_ref_seqs:,}")
    print(f"[normalise:{sample}] Normalisation factor:      {norm_factor:.4f}")
    print(f"[normalise:{sample}] Unique VCGs observed:      {len(grouped):,}")

    # ── Zero-fill: include all reference VCGs even if not observed ───────────
    ref_df = pd.DataFrame({"VCG": sorted(vcg_set)})
    ref_df["VCG_Len"] = ref_df["VCG"].str.len()
    grouped = ref_df.merge(grouped, on=["VCG", "VCG_Len"], how="left")
    grouped["Count"]      = grouped["Count"].fillna(0).astype(int)
    grouped["Norm_Reads"] = grouped["Norm_Reads"].fillna(0.0)

    # ── Convert VCG to RNA notation for output ────────────────────────────────
    grouped["VCG_DNA"] = grouped["VCG"]
    grouped["VCG"]     = grouped["VCG"].str.replace("T", "U", regex=False)

    # ── Reorder and write ─────────────────────────────────────────────────────
    out = grouped[["VCG_Len", "VCG", "VCG_DNA", "Count", "Norm_Reads"]]
    out = out.sort_values(["VCG_Len", "VCG"], ascending=[True, True])
    out.to_csv(normalised_path, sep="\t", index=False)
    print(f"[normalise:{sample}] Normalised → {normalised_path}")

    if csv_path:
        out.to_csv(csv_path, index=False)
        print(f"[normalise:{sample}] Normalised → {csv_path}")


# ── Entry points ──────────────────────────────────────────────────────────────

if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    normalise(
        collapsed_path  = snakemake.input.collapsed,
        vcg_list_path   = snakemake.input.vcg_list,
        normalised_path = snakemake.output.normalised,
        sample          = snakemake.wildcards.sample,
        csv_path        = snakemake.output.csv,
    )

elif __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Group, count, and normalise UMI-collapsed reads.")
    ap.add_argument("--collapsed",  required=True)
    ap.add_argument("--vcg-list",   required=True, dest="vcg_list")
    ap.add_argument("--normalised", required=True)
    args = ap.parse_args()
    normalise(args.collapsed, args.vcg_list, args.normalised)
