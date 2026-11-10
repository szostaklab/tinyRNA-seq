#!/usr/bin/env python3
"""
smallrna_comparison.py  (data only)

Pool counts.smallRNA.withCPM.tsv from all experiments into one combined table:
    smallRNA_all_experiments.csv   (comma)
    smallRNA_all_experiments.tsv   (tab)

Snakemake script: entry point.
"""

import sys
from pathlib import Path
from typing import List

import pandas as pd

if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    counts_files  = list(snakemake.input)
    out_csv       = Path(snakemake.output.csv)
    display_names = dict(snakemake.config.get("display_names", {}))
else:
    raise RuntimeError("smallrna_comparison.py must be run via Snakemake script: directive")

RESULT_SUFFIX      = "_results"
OUTDIR             = out_csv.parent
OUTDIR.mkdir(parents=True, exist_ok=True)
OUT_COMBINED_TABLE = OUTDIR / "smallRNA_all_experiments.tsv"

files = sorted(set(Path(p) for p in counts_files))
print(f"Found {len(files)} experiment result files:")
for f in files:
    print(" -", f)

dfs: List[pd.DataFrame] = []
for fp in files:
    if not fp.is_file():
        print(f"[WARN] Missing: {fp}")
        continue
    df = pd.read_csv(fp, sep="\t")
    if "gene_type" not in df.columns or "CPM" not in df.columns:
        df = pd.read_csv(fp, sep="\t", header=None)
        ncol = df.shape[1]
        df.columns = (
            ["Geneid", "counts", "CPM", "gene_name", "gene_type"]
            + [f"extra_{i}" for i in range(ncol - 5)]
        )
    exp_name = fp.parent.name
    if exp_name.endswith(RESULT_SUFFIX):
        exp_name = exp_name[: -len(RESULT_SUFFIX)]
    exp_name = display_names.get(exp_name, exp_name)
    df["experiment"] = exp_name
    dfs.append(df)

if not dfs:
    raise RuntimeError("No smallRNA tables found.")

df_all = pd.concat(dfs, ignore_index=True)
df_all.to_csv(OUT_COMBINED_TABLE, sep="\t", index=False)
df_all.to_csv(out_csv, index=False)
print(f"Total pooled smallRNA rows: {df_all.shape[0]}")
print(f"Combined table → {OUT_COMBINED_TABLE}")
print(f"Combined table → {out_csv}")
