#!/usr/bin/env python3
"""
umi_collapse_dno.py

PCR duplicate removal for DNO (non-VCG) experiments using directional UMI
clustering (UMI-tools method). DNO variant of umi_collapse.py.

Workflow:
  1. Group reads by VCG (species) — deduplication is always within-species.
  2. Within each group, build the full UMI string by concatenating all random
     base columns defined in umi_dedup_columns (excluding VCG).
  3. Apply directional clustering on the UMI strings:
       - Two UMIs are neighbours if Hamming distance <= 1.
       - UMI B is absorbed into cluster A if count(A) >= 2*count(B) - 1.
       - UMIs processed in descending count order so the most abundant
         always becomes the cluster root.
  4. Emit one representative row per cluster (the root UMI's read) with a
     Count column recording how many reads collapsed into that cluster.
  5. Write a per-VCG collapse stats TSV (collapse_stats.tsv) recording
     reads_before, reads_after, reads_removed, and pct_collapsed per species.

If umi_dedup_columns is empty for the layout, all reads are passed through
unchanged (no deduplication).

Called by Snakemake (rule umi_collapse) or directly:
    python scripts/umi_collapse.py \\
        --filtered vcg_filtered.tsv \\
        --collapsed umi_collapsed.tsv \\
        --collapse-stats collapse_stats.tsv \\
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

def load_dedup_columns(schemas_file: str, read_layout: str) -> list[str]:
    with open(schemas_file) as fh:
        s = yaml.safe_load(fh)
    layouts = s.get("read_layouts", {})
    if read_layout not in layouts:
        raise ValueError(
            f"Unknown read_layout '{read_layout}'. Available: {list(layouts)}"
        )
    return layouts[read_layout].get("umi_dedup_columns", [])


# ── Directional clustering ────────────────────────────────────────────────────

def hamming1_neighbors(umi: str) -> list[str]:
    """All strings at Hamming distance exactly 1 from umi."""
    neighbors = []
    for i, base in enumerate(umi):
        for b in BASES:
            if b != base:
                neighbors.append(umi[:i] + b + umi[i + 1:])
    return neighbors


def directional_cluster(umi_counts: dict[str, int]) -> dict[str, str]:
    """
    Directional UMI clustering.
    Returns {umi: cluster_root} for every UMI in umi_counts.

    UMIs are processed in descending count order. Each unassigned UMI
    becomes a cluster root and absorbs any 1-error neighbour B where
    count(root) >= 2*count(B) - 1.
    """
    umis_sorted = sorted(umi_counts, key=lambda u: -umi_counts[u])
    cluster_of: dict[str, str] = {}

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


# ── Per-VCG collapse ──────────────────────────────────────────────────────────

def collapse_group(
    group_df: pd.DataFrame,
    umi_random_cols: list[str],
) -> tuple[pd.DataFrame, dict]:
    """
    Apply directional UMI clustering within one VCG group.
    Returns (collapsed_df, stats_dict).
    collapsed_df has one representative row per cluster with a Count column.
    stats_dict has reads_before, reads_after, reads_removed, pct_collapsed.
    """
    n_before = len(group_df)
    umi_strings = group_df[umi_random_cols].fillna("").agg("".join, axis=1).tolist()

    umi_counts: dict[str, int] = defaultdict(int)
    for u in umi_strings:
        umi_counts[u] += 1

    cluster_of = directional_cluster(dict(umi_counts))

    # Sum raw read counts per cluster root
    cluster_totals: dict[str, int] = defaultdict(int)
    for u, root in cluster_of.items():
        cluster_totals[root] += umi_counts[u]

    # One representative row per cluster root (first occurrence of root UMI)
    seen_roots: set[str] = set()
    result_rows = []
    for row_idx, umi in enumerate(umi_strings):
        root = cluster_of.get(umi, umi)
        if root not in seen_roots:
            seen_roots.add(root)
            row = group_df.iloc[row_idx].to_dict()
            row["Count"] = cluster_totals[root]
            row["UMI"] = root
            result_rows.append(row)

    n_after   = len(result_rows)
    n_removed = n_before - n_after
    stats = {
        "reads_before":   n_before,
        "reads_after":    n_after,
        "reads_removed":  n_removed,
        "pct_collapsed":  round(100 * n_removed / max(n_before, 1), 2),
    }
    return pd.DataFrame(result_rows), stats


# ── Core ──────────────────────────────────────────────────────────────────────

def umi_collapse(
    filtered_path: str,
    collapsed_path: str,
    stats_path: str,
    read_layout: str,
    schemas_file: str,
    sample: str = "",
) -> None:
    Path(collapsed_path).parent.mkdir(parents=True, exist_ok=True)

    dedup_cols = load_dedup_columns(schemas_file, read_layout)
    df = pd.read_csv(filtered_path, sep="\t", dtype=str)
    n_before = len(df)

    if not dedup_cols:
        print(f"[umi_collapse:{sample}] Layout '{read_layout}' has no UMI dedup columns — passing through.")
        df = df.drop(columns=["read_id"], errors="ignore")
        df.to_csv(collapsed_path, sep="\t", index=False)
        # Write a minimal stats file
        pd.DataFrame([{
            "VCG": "*", "VCG_Len": "*",
            "reads_before": n_before, "reads_after": n_before,
            "reads_removed": 0, "pct_collapsed": 0.0,
        }]).to_csv(stats_path, sep="\t", index=False)
        print(f"[umi_collapse:{sample}] Reads retained: {n_before:,}")
        return

    missing = [c for c in dedup_cols if c not in df.columns]
    if missing:
        raise ValueError(
            f"Dedup columns {missing} not found in {filtered_path}. "
            f"Available: {list(df.columns)}"
        )

    # VCG = grouping key; remaining dedup_cols = UMI random bases for clustering
    umi_random_cols = [c for c in dedup_cols if c != "VCG"]

    parts       = []
    stats_rows  = []

    for vcg, group in df.groupby("VCG", sort=True):
        collapsed_group, stats = collapse_group(group.reset_index(drop=True), umi_random_cols)
        parts.append(collapsed_group)
        vcg_len = group["VCG_Len"].iloc[0] if "VCG_Len" in group.columns else ""
        stats_rows.append({
            "VCG":           vcg,
            "VCG_Len":       vcg_len,
            "reads_before":  stats["reads_before"],
            "reads_after":   stats["reads_after"],
            "reads_removed": stats["reads_removed"],
            "pct_collapsed": stats["pct_collapsed"],
        })

    df_collapsed = pd.concat(parts, ignore_index=True) if parts else df.iloc[:0].copy()
    n_after   = len(df_collapsed)
    n_removed = n_before - n_after

    # Console summary
    print(f"[umi_collapse:{sample}] Read layout:  {read_layout}")
    print(f"[umi_collapse:{sample}] UMI columns:  {umi_random_cols}")
    print(f"[umi_collapse:{sample}] Before: {n_before:,}  |  Removed: {n_removed:,}  ({n_removed / max(n_before, 1):.2%})")
    print(f"[umi_collapse:{sample}] After:  {n_after:,}")

    # Per-VCG stats TSV — sorted by reads_removed descending so heaviest
    # collapsers appear first
    stats_df = (
        pd.DataFrame(stats_rows)
        .sort_values("reads_removed", ascending=False)
        .reset_index(drop=True)
    )
    stats_df.to_csv(stats_path, sep="\t", index=False)
    print(f"[umi_collapse:{sample}] Stats →     {stats_path}")

    df_collapsed = df_collapsed.drop(columns=["read_id"], errors="ignore")
    df_collapsed.to_csv(collapsed_path, sep="\t", index=False)
    print(f"[umi_collapse:{sample}] Collapsed → {collapsed_path}")


# ── Entry points ──────────────────────────────────────────────────────────────

if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    umi_collapse(
        filtered_path  = snakemake.input.filtered,
        collapsed_path = snakemake.output.collapsed,
        stats_path     = snakemake.output.collapse_stats,
        read_layout    = snakemake.params.read_layout,
        schemas_file   = snakemake.input.schemas,
        sample         = snakemake.wildcards.sample,
    )

elif __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Directional UMI clustering within VCG species groups.")
    ap.add_argument("--filtered",       required=True)
    ap.add_argument("--collapsed",      required=True)
    ap.add_argument("--collapse-stats", required=True, dest="collapse_stats")
    ap.add_argument("--read-layout",    required=True)
    ap.add_argument("--schemas",        default="adaptor_schemas.yaml")
    args = ap.parse_args()
    umi_collapse(args.filtered, args.collapsed, args.collapse_stats, args.read_layout, args.schemas)
