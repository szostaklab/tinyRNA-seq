#!/usr/bin/env python3
"""
ligation_junction.py  (data only)

Per-VCG nucleotide-frequency analysis of the random bases flanking the ligation
junction, to detect sequence bias introduced by the T4 ligase. This is the
DATA half: it computes the per-VCG base-composition stats and writes

    ligation_junction_summary.tsv     — per-VCG stats for the junction + control columns

The sequence-logo figures built from the same base-frequency matrices live in
figures/ligation_junction_plots.py. The stats computation (compute_per_vcg_stats)
was factored out of the old plotting routine so the summary table can be produced
without importing matplotlib / logomaker.

Junction column (immediately 3' of the VCG insert):
    three_prime MADAP  → "TRR"   ·  FOUR_N → "N4"   ·  NORMAL → none
Control column (5' random tag, ligated differently — expected less biased):
    five_prime FOUR_N / EIGHT_N → "5_Adap"   ·  DEFINED / FLUSH → none
"Abnormal" = any position where one base exceeds max_prob_thresh (0.51) or falls
below min_prob_thresh (0.07).

Snakemake script, or standalone:
    python scripts/ligation_junction.py \\
        --filtered vcg_filtered.tsv --schemas adaptor_schemas.yaml \\
        --read-layout FOUR_N_MADAP --summary-out ligation_junction_summary.tsv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml

BASES = ["A", "C", "G", "T"]

MAX_PROB_THRESH = 0.51
MIN_PROB_THRESH = 0.07

# Which column sits immediately 3' of VCG at the ligation junction
JUNCTION_COL = {"NORMAL": None, "FOUR_N": "N4", "MADAP": "TRR"}
# Which column is the 5' random tag (control)
CONTROL_COL = {"DEFINED": None, "FLUSH": None, "FOUR_N": "5_Adap", "EIGHT_N": "5_Adap"}


# ── Schema loading ─────────────────────────────────────────────────────────────

def load_columns(schemas_file: str, read_layout: str) -> tuple[str | None, str | None]:
    """Return (junction_col, control_col) for this layout, either may be None."""
    with open(schemas_file) as fh:
        s = yaml.safe_load(fh)
    layout     = s["read_layouts"][read_layout]
    five_type  = layout["five_prime"]
    three_type = layout["three_prime"]
    return JUNCTION_COL[three_type], CONTROL_COL[five_type]


# ── Base frequency ─────────────────────────────────────────────────────────────

def base_freq_df(series: pd.Series) -> pd.DataFrame:
    """Per-position base probability matrix. Always columns [A, C, G, T]; missing → 0."""
    series = series.dropna().astype(str)
    if len(series) == 0:
        return pd.DataFrame(columns=BASES)
    seq_len = len(series.iloc[0])
    probs = {b: [] for b in BASES}
    for pos in range(seq_len):
        counts = series.str[pos].value_counts() / len(series)
        for b in BASES:
            probs[b].append(float(counts.get(b, 0.0)))
    return pd.DataFrame(probs)


def is_abnormal(freq: pd.DataFrame,
                max_thresh: float = MAX_PROB_THRESH,
                min_thresh: float = MIN_PROB_THRESH) -> bool:
    if freq.empty:
        return False
    arr = freq.to_numpy()
    return bool((arr.max() > max_thresh) or (arr.min() < min_thresh))


def compute_per_vcg_stats(df: pd.DataFrame, col: str) -> dict:
    """Per-VCG base-composition stats (no plotting).

    Returns {vcg: {VCG_Len, n_reads, max_prob, min_prob, abnormal}} for all VCGs.
    Factored out of the original plot_per_vcg so the summary can be built without
    importing any plotting libraries.
    """
    stats = {}
    for (vcg_len, vcg), group in df.groupby(["VCG_Len", "VCG"]):
        freq = base_freq_df(group[col])
        arr  = freq.to_numpy()
        stats[vcg] = {
            "VCG_Len":  int(vcg_len),
            "n_reads":  len(group),
            "max_prob": round(float(arr.max()), 4) if arr.size else float("nan"),
            "min_prob": round(float(arr.min()), 4) if arr.size else float("nan"),
            "abnormal": is_abnormal(freq),
        }
    return stats


# ── Summary ────────────────────────────────────────────────────────────────────

def build_summary(junc_stats: dict | None, ctrl_stats: dict | None,
                  junc_col: str | None, ctrl_col: str | None) -> pd.DataFrame:
    """Merge junction and control stats into a single per-VCG summary DataFrame."""
    all_vcgs = sorted(
        set(junc_stats or {}) | set(ctrl_stats or {}),
        key=lambda v: ((junc_stats or ctrl_stats or {}).get(v, {}).get("VCG_Len", 0), v),
    )
    records = []
    for vcg in all_vcgs:
        row: dict = {"VCG": vcg}
        base = (junc_stats or {}).get(vcg) or (ctrl_stats or {}).get(vcg) or {}
        row["VCG_Len"] = base.get("VCG_Len", "")
        row["n_reads"] = base.get("n_reads", "")
        if junc_stats and vcg in junc_stats:
            j = junc_stats[vcg]
            row[f"{junc_col}_max_prob"] = j["max_prob"]
            row[f"{junc_col}_min_prob"] = j["min_prob"]
            row[f"{junc_col}_abnormal"] = j["abnormal"]
        if ctrl_stats and vcg in ctrl_stats:
            c = ctrl_stats[vcg]
            row[f"{ctrl_col}_max_prob"] = c["max_prob"]
            row[f"{ctrl_col}_min_prob"] = c["min_prob"]
            row[f"{ctrl_col}_abnormal"] = c["abnormal"]
        records.append(row)
    return pd.DataFrame(records).sort_values(["VCG_Len", "VCG"]).reset_index(drop=True)


# ── Core ───────────────────────────────────────────────────────────────────────

def ligation_junction(filtered_path: str, schemas_file: str, read_layout: str,
                      summary_out: str, sample: str = "") -> None:
    junc_col, ctrl_col = load_columns(schemas_file, read_layout)
    Path(summary_out).parent.mkdir(parents=True, exist_ok=True)

    if junc_col is None and ctrl_col is None:
        print(f"[ligation_junction:{sample}] Layout {read_layout} has no random columns — nothing to do.")
        pd.DataFrame().to_csv(summary_out, sep="\t", index=False)
        return

    df = pd.read_csv(filtered_path, sep="\t", dtype=str)
    n  = len(df)
    print(f"[ligation_junction:{sample}] Layout: {read_layout}  |  Reads: {n:,}")
    print(f"[ligation_junction:{sample}] Junction column: {junc_col}  |  Control column: {ctrl_col}")

    junc_stats = ctrl_stats = None
    if junc_col and junc_col in df.columns:
        junc_stats = compute_per_vcg_stats(df, junc_col)
        n_biased = sum(v["abnormal"] for v in junc_stats.values())
        print(f"[ligation_junction:{sample}]   {n_biased}/{len(junc_stats)} VCGs show junction bias")
    elif junc_col:
        print(f"[ligation_junction:{sample}] Warning — junction column '{junc_col}' not found in TSV")

    if ctrl_col and ctrl_col in df.columns:
        ctrl_stats = compute_per_vcg_stats(df, ctrl_col)
        n_biased = sum(v["abnormal"] for v in ctrl_stats.values())
        print(f"[ligation_junction:{sample}]   {n_biased}/{len(ctrl_stats)} VCGs show 5' bias")
    elif ctrl_col:
        print(f"[ligation_junction:{sample}] Warning — control column '{ctrl_col}' not found in TSV")

    summary = build_summary(junc_stats, ctrl_stats, junc_col, ctrl_col)
    summary.to_csv(summary_out, sep="\t", index=False)
    print(f"[ligation_junction:{sample}] Summary → {summary_out}")


# ── Entry points ───────────────────────────────────────────────────────────────

if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    ligation_junction(
        filtered_path = snakemake.input.filtered,
        schemas_file  = snakemake.input.schemas,
        read_layout   = snakemake.params.read_layout,
        summary_out   = snakemake.output.summary,
        sample        = snakemake.wildcards.sample,
    )

elif __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Ligation junction base-bias analysis (data only).")
    ap.add_argument("--filtered",    required=True)
    ap.add_argument("--schemas",     default="adaptor_schemas.yaml")
    ap.add_argument("--read-layout", required=True)
    ap.add_argument("--summary-out", required=True)
    args = ap.parse_args()
    ligation_junction(args.filtered, args.schemas, args.read_layout, args.summary_out)
