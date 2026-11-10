#!/usr/bin/env python3
"""
bio_collapse_to_fasta.py

Take parsed_reads.tsv (output of parse_reads.py) for a biological RNA sample,
apply directional Hamming-1 UMI clustering, and write the collapsed inserts
as FASTA for Bowtie1 alignment (-f flag).

One FASTA entry = one unique molecule (post UMI-collapse).
Sequence = VCG column (the biological RNA insert, DNA form).

For layouts with no UMI dedup columns (e.g. FLUSH_NORMAL), all reads are
written through unchanged.

Called by Snakemake (rule bio_collapse_to_fasta) or directly:
    python scripts/bio_collapse_to_fasta.py \\
        --parsed results/sample/parsed_reads.tsv \\
        --fasta  results/sample/collapsed.fasta \\
        --stats  results/sample/collapse_stats.tsv \\
        --read-layout FOUR_N_MADAP \\
        --schemas adaptor_schemas.yaml
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd
import yaml


BASES = "ACGT"


# ── Schema loading ────────────────────────────────────────────────────────────

def load_dedup_columns(schemas_file: str, read_layout: str) -> list:
    with open(schemas_file) as fh:
        s = yaml.safe_load(fh)
    layouts = s.get("read_layouts", {})
    if read_layout not in layouts:
        raise ValueError(
            f"Unknown read_layout '{read_layout}'. Available: {list(layouts)}"
        )
    return layouts[read_layout].get("umi_dedup_columns", [])


# ── Directional UMI clustering ────────────────────────────────────────────────

def hamming1_neighbors(umi: str) -> list:
    neighbors = []
    for i, base in enumerate(umi):
        for b in BASES:
            if b != base:
                neighbors.append(umi[:i] + b + umi[i + 1:])
    return neighbors


def directional_cluster(umi_counts: dict) -> dict:
    """
    Returns {umi: cluster_root} for every UMI.
    UMIs processed descending by count; UMI B absorbed into A if
    count(A) >= 2*count(B) - 1 and Hamming(A,B) == 1.
    """
    umis_sorted = sorted(umi_counts, key=lambda u: -umi_counts[u])
    cluster_of: dict = {}
    for umi in umis_sorted:
        if umi in cluster_of:
            continue
        cluster_of[umi] = umi
        count_root = umi_counts[umi]
        for neighbor in hamming1_neighbors(umi):
            if neighbor in umi_counts and neighbor not in cluster_of:
                if count_root >= 2 * umi_counts[neighbor] - 1:
                    cluster_of[neighbor] = umi
    return cluster_of


def collapse_group(group_df: pd.DataFrame, umi_random_cols: list) -> tuple:
    """One representative row per cluster root. Returns (row_list, stats_dict)."""
    n_before = len(group_df)
    umi_strings = group_df[umi_random_cols].fillna("").agg("".join, axis=1).tolist()

    umi_counts: dict = defaultdict(int)
    for u in umi_strings:
        umi_counts[u] += 1

    cluster_of = directional_cluster(dict(umi_counts))

    seen_roots: set = set()
    result_rows = []
    for row_idx, umi in enumerate(umi_strings):
        root = cluster_of.get(umi, umi)
        if root not in seen_roots:
            seen_roots.add(root)
            result_rows.append(group_df.iloc[row_idx])

    n_after   = len(result_rows)
    n_removed = n_before - n_after
    stats = {
        "reads_before":  n_before,
        "reads_after":   n_after,
        "reads_removed": n_removed,
        "pct_collapsed": round(100 * n_removed / max(n_before, 1), 2),
    }
    return result_rows, stats


# ── Core ──────────────────────────────────────────────────────────────────────

def collapse_to_fasta(
    parsed_path: str,
    fasta_path: str,
    stats_path: str,
    read_layout: str,
    schemas_file: str,
    sample: str = "",
) -> None:
    Path(fasta_path).parent.mkdir(parents=True, exist_ok=True)
    tag = f"[bio_collapse:{sample}]" if sample else "[bio_collapse]"

    dedup_cols = load_dedup_columns(schemas_file, read_layout)
    df = pd.read_csv(parsed_path, sep="\t", dtype=str)
    n_before = len(df)

    if "VCG" not in df.columns:
        raise ValueError(
            f"'VCG' column not found in {parsed_path}. Available: {list(df.columns)}"
        )

    if not dedup_cols:
        print(f"{tag} Layout '{read_layout}' has no UMI dedup columns — passing through.")
        result_rows = [df.iloc[i] for i in range(len(df))]
        pd.DataFrame([{
            "VCG": "*", "reads_before": n_before, "reads_after": n_before,
            "reads_removed": 0, "pct_collapsed": 0.0,
        }]).to_csv(stats_path, sep="\t", index=False)
    else:
        umi_random_cols = [c for c in dedup_cols if c != "VCG"]
        missing = [c for c in umi_random_cols if c not in df.columns]
        if missing:
            raise ValueError(f"UMI columns {missing} not found in {parsed_path}.")

        result_rows = []
        stats_rows = []
        for vcg, group in df.groupby("VCG", sort=True):
            rows, stats = collapse_group(group.reset_index(drop=True), umi_random_cols)
            result_rows.extend(rows)
            vcg_len = group["VCG_Len"].iloc[0] if "VCG_Len" in group.columns else ""
            stats_rows.append({"VCG": vcg, "VCG_Len": vcg_len, **stats})
        pd.DataFrame(stats_rows).to_csv(stats_path, sep="\t", index=False)

    n_after   = len(result_rows)
    n_removed = n_before - n_after

    # Write FASTA
    written = 0
    with open(fasta_path, "w") as fh:
        for i, row in enumerate(result_rows):
            seq = str(row["VCG"]) if hasattr(row, "__getitem__") else ""
            if seq and seq != "nan":
                fh.write(f">read_{i + 1}\n{seq}\n")
                written += 1

    print(f"{tag} Read layout:         {read_layout}")
    print(f"{tag} Reads before:        {n_before:,}")
    print(f"{tag} Reads after:         {n_after:,}  ({100 * n_after / max(n_before, 1):.1f}%)")
    print(f"{tag} Reads removed:       {n_removed:,}  ({100 * n_removed / max(n_before, 1):.1f}%)")
    print(f"{tag} FASTA entries:       {written:,}")
    print(f"{tag} FASTA  → {fasta_path}")
    print(f"{tag} Stats  → {stats_path}")


# ── Entry points ──────────────────────────────────────────────────────────────

if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    collapse_to_fasta(
        parsed_path  = snakemake.input.parsed,
        fasta_path   = snakemake.output.fasta,
        stats_path   = snakemake.output.stats,
        read_layout  = snakemake.params.read_layout,
        schemas_file = snakemake.input.schemas,
        sample       = snakemake.wildcards.sample,
    )

elif __name__ == "__main__":
    ap = argparse.ArgumentParser(description="UMI-collapse parsed reads to FASTA for Bowtie.")
    ap.add_argument("--parsed",      required=True)
    ap.add_argument("--fasta",       required=True)
    ap.add_argument("--stats",       required=True)
    ap.add_argument("--read-layout", required=True)
    ap.add_argument("--schemas",     default="adaptor_schemas.yaml")
    args = ap.parse_args()
    collapse_to_fasta(args.parsed, args.fasta, args.stats, args.read_layout, args.schemas)
