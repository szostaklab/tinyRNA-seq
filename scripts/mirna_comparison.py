#!/usr/bin/env python3
"""
mirna_comparison.py  (data only)

Pool per-sample miRNA CPM and correlate against an external reference. Writes
the tables listed below.

Outputs:
    miRNA_all_experiments.csv / .tsv                 pooled miRNA rows, all samples
    spearman/spearman_summary.csv                    Spearman rho per sample vs reference
    spearman/rank_delta_full.csv                     per-miRNA rank delta vs reference

Snakemake script: entry point.
"""

import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    counts_files  = list(snakemake.input.counts)
    out_csv       = Path(snakemake.output.table)
    ref_csv_path  = Path(snakemake.input.ref_csv)
    spearman_out  = Path(snakemake.output.spearman)
    cpm_threshold = float(getattr(snakemake.params, "cpm_threshold", 1.0))
    display_names = dict(snakemake.config.get("display_names", {}))
else:
    raise RuntimeError("mirna_comparison.py must be run via Snakemake script: directive")

RESULT_SUFFIX = "_results"

OUTDIR = out_csv.parent
OUTDIR.mkdir(parents=True, exist_ok=True)
SPEARMAN_DIR = spearman_out.parent
SPEARMAN_DIR.mkdir(parents=True, exist_ok=True)
OUT_COMBINED_TABLE = OUTDIR / "miRNA_all_experiments.tsv"


# ── 1. Pooled miRNA table ──────────────────────────────────────────────────────

def load_mirna_tables(counts_files: Iterable) -> pd.DataFrame:
    files = sorted(set(Path(p) for p in counts_files))
    if not files:
        raise RuntimeError("No input files provided.")
    print(f"Found {len(files)} experiment result files:")
    for f in files:
        print(" -", f)

    dfs: List[pd.DataFrame] = []
    for file_path in files:
        if not file_path.is_file():
            print(f"[WARN] Missing file: {file_path}")
            continue
        df = pd.read_csv(file_path, sep="\t")
        if "gene_type" not in df.columns or "CPM" not in df.columns:
            df = pd.read_csv(file_path, sep="\t", header=None)
            ncol = df.shape[1]
            df.columns = (
                ["Geneid", "counts", "CPM", "gene_name", "gene_type"]
                + [f"extra_{i}" for i in range(ncol - 5)]
            )
        df = df[df["gene_type"] == "miRNA"].copy()
        exp_name = file_path.parent.name
        if exp_name.endswith(RESULT_SUFFIX):
            exp_name = exp_name[: -len(RESULT_SUFFIX)]
        exp_name = display_names.get(exp_name, exp_name)
        df["experiment"] = exp_name
        dfs.append(df)

    if not dfs:
        raise RuntimeError("No miRNA tables found.")

    df_all = pd.concat(dfs, ignore_index=True)
    df_all.to_csv(OUT_COMBINED_TABLE, sep="\t", index=False)
    df_all.to_csv(out_csv, index=False)
    print(f"\nTotal pooled miRNA rows: {df_all.shape[0]}")
    print(f"Combined table → {OUT_COMBINED_TABLE}")
    return df_all


df_all = load_mirna_tables(counts_files)


# ── 2. Spearman rank correlation vs external reference ─────────────────────────

SAMPLES = {Path(f).parent.name: Path(f) for f in counts_files}

# Group samples by protocol prefix (first '-'-delimited token).
_grp_map: Dict[str, List[str]] = defaultdict(list)
for s in SAMPLES:
    _grp_map[s.split("-")[0]].append(s)
SAMPLE_GROUPS = dict(_grp_map)
print(f"\nSample groups for rank analysis: {SAMPLE_GROUPS}")


def norm_external(name: str):
    """miRBase mature name → Ensembl gene-symbol style (mmu-miR-451a → Mir451a)."""
    if not isinstance(name, str) or name == "-":
        return None
    m = re.match(r"^mmu-let-(.+)", name)
    if m:
        base = re.sub(r"-(3p|5p)$", "", m.group(1))
        base = re.sub(r"-\d+$", "", base)
        return "Mirlet" + base
    m = re.match(r"^mmu-miR-(.+)", name)
    if m:
        base = re.sub(r"-(3p|5p)$", "", m.group(1))
        base = re.sub(r"-\d+$", "", base)
        return "Mir" + base
    return None


def norm_bcrun(name: str) -> str:
    """Strip paralogue copy-number suffix: Mir133a-1 → Mir133a."""
    return re.sub(r"-\d+$", "", name)


def load_sample_mirna(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    mirna = df[(df["gene_type"] == "miRNA") & (df["CPM"].astype(float) >= cpm_threshold)].copy()
    mirna["norm_name"] = mirna["gene_name"].map(norm_bcrun)
    return mirna.groupby("norm_name", as_index=False)["CPM"].sum()


def fuzzy_merge(ref_gene: pd.DataFrame, sample: pd.DataFrame) -> pd.DataFrame:
    """Exact match + 'a'-suffix fallback for miRNA name mismatches."""
    sample_set = set(sample["norm_name"])
    exact = ref_gene[ref_gene["norm_name"].isin(sample_set)].copy()
    exact["match_type"] = "exact"
    unmatched = ref_gene[~ref_gene["norm_name"].isin(sample_set)]
    fallback_rows = []
    for _, row in unmatched.iterrows():
        n = row["norm_name"]
        if n + "a" in sample_set:
            fallback_rows.append({**row, "norm_name": n + "a", "match_type": "append_a"})
        elif n.endswith("a") and n[:-1] in sample_set:
            fallback_rows.append({**row, "norm_name": n[:-1], "match_type": "strip_a"})
    combined = pd.concat([exact, pd.DataFrame(fallback_rows)], ignore_index=True) if fallback_rows else exact
    return combined.merge(sample, on="norm_name", how="inner")


def get_fuzzy_key(ref_df: pd.DataFrame, sample: pd.DataFrame) -> dict:
    sample_set = set(sample["norm_name"])
    mapping = {}
    for _, row in ref_df.iterrows():
        n = row["norm_name"]
        if n in sample_set:
            mapping[n] = n
        elif n + "a" in sample_set:
            mapping[n] = n + "a"
        elif n.endswith("a") and n[:-1] in sample_set:
            mapping[n] = n[:-1]
    return mapping


ref = pd.read_csv(ref_csv_path)
ref["RNA_seq"] = ref["RNA_seq"].str.replace("T", "U", regex=False)
ref_mirna = ref[ref["type"] == "miRNA"].copy()
ref_mirna["norm_name"] = ref_mirna["RNA_name"].map(norm_external)
ref_mirna = ref_mirna.dropna(subset=["norm_name"])
ref_gene = (
    ref_mirna.groupby("norm_name", as_index=False)["expression"]
    .sum()
    .rename(columns={"expression": "ref_expression"})
)
print(f"\nExternal reference: {len(ref_mirna)} rows → {len(ref_gene)} unique genes")

sp_results = []
for sample, path in SAMPLES.items():
    if not path.exists():
        print(f"  SKIP {sample}: file not found")
        continue
    s_df = load_sample_mirna(path)
    merged = fuzzy_merge(ref_gene, s_df)
    nonzero = merged[(merged["ref_expression"] > 0) & (merged["CPM"] > 0)]
    rho, pval = spearmanr(merged["ref_expression"], merged["CPM"])
    rho_nz, pval_nz = spearmanr(nonzero["ref_expression"], nonzero["CPM"])
    sp_results.append({
        "sample": sample, "n_matched": len(merged), "n_nonzero": len(nonzero),
        "rho_all": rho, "pval_all": pval,
        "rho_nonzero": rho_nz, "pval_nonzero": pval_nz,
    })
    print(f"{sample:40s}  n={len(merged):4d}  rho={rho:.3f}  (p={pval:.2e})")

results_df = pd.DataFrame(sp_results)
results_df.to_csv(spearman_out, index=False)
print(f"\nSpearman summary → {spearman_out}")


# ── 3. Per-miRNA rank-delta table ──────────────────────────────────────────────

canonical_seq = ref_mirna.sort_values("expression", ascending=False).groupby("norm_name")["RNA_seq"].first()
ref_gene["seq"]      = ref_gene["norm_name"].map(canonical_seq)
ref_gene["first_nt"] = ref_gene["seq"].str[0]
ref_gene["last_nt"]  = ref_gene["seq"].str[-1]
ref_gene_ranked = ref_gene.sort_values("ref_expression", ascending=False).reset_index(drop=True)
ref_gene_ranked["ref_rank"] = ref_gene_ranked.index + 1

rank_abs  = ref_gene_ranked[["ref_rank", "norm_name", "ref_expression", "seq", "first_nt", "last_nt"]].set_index("norm_name")
delta_tbl = rank_abs[["ref_rank", "ref_expression", "first_nt", "last_nt"]].copy()

for sample, path in SAMPLES.items():
    if not path.exists():
        continue
    s_df   = load_sample_mirna(path)
    merged = fuzzy_merge(ref_gene_ranked, s_df)
    merged["bcrun_rank"] = merged["CPM"].rank(ascending=False, method="min").astype(int)
    key_map     = get_fuzzy_key(ref_gene_ranked, s_df)
    rank_lookup = merged.set_index("norm_name")["bcrun_rank"]
    rank_abs[sample]  = rank_abs.index.map(
        lambda n: rank_lookup.get(key_map.get(n, n), np.nan) if n in key_map else np.nan)
    delta_tbl[sample] = rank_abs.index.map(
        lambda n: rank_lookup.get(key_map.get(n, n), np.nan) - rank_abs.loc[n, "ref_rank"]
        if n in key_map else np.nan)

for grp, members in SAMPLE_GROUPS.items():
    cols = [s for s in members if s in delta_tbl.columns]
    if cols:
        delta_tbl[f"delta_{grp}"] = delta_tbl[cols].mean(axis=1)

delta_tbl.reset_index().sort_values("ref_rank").to_csv(SPEARMAN_DIR / "rank_delta_full.csv", index=False)
print(f"rank_delta_full.csv → {SPEARMAN_DIR / 'rank_delta_full.csv'}")
print("\nAll miRNA comparison tables saved. DONE!")
