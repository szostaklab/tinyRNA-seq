#!/usr/bin/env python3
"""
mirna_dispersion.py

Compute miRNA dispersion metrics (Gini, HHI_norm, CV, top-k shares) and write
them as a table. Data only — Lorenz-curve and rank-CPM plots live in the
figures layer (figures/mirna_dispersion_plots.py).

Snakemake script: entry point.
"""

import sys

import numpy as np
import pandas as pd
from pathlib import Path

if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    infile        = snakemake.input.counts
    out_txt       = snakemake.output.disp
    sample        = snakemake.wildcards.sample
    cpm_threshold  = float(getattr(snakemake.params, "cpm_threshold", 1.0))
else:
    raise RuntimeError("mirna_dispersion.py must be run via Snakemake script: directive")

Path(out_txt).parent.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(infile, sep="\t")
mir = df[df["gene_type"] == "miRNA"].copy()
mir["CPM"] = pd.to_numeric(mir["CPM"], errors="coerce")
mir = mir.dropna(subset=["CPM"])
mir = mir[mir["CPM"] >= cpm_threshold].sort_values("CPM", ascending=False).reset_index(drop=True)
vals = mir["CPM"].values


def gini(x):
    x = np.asarray(x, float)
    if x.size == 0:
        return 0.0
    x = np.where(x < 0, 0, x)
    mean = x.mean()
    if mean == 0:
        return 0.0
    mad = np.abs(np.subtract.outer(x, x)).mean()
    return 0.5 * mad / mean


def hhi_norm(x):
    x = np.asarray(x, float)
    if x.size == 0:
        return 0.0
    s = x.sum()
    if s == 0:
        return 0.0
    shares = x / s
    h = np.sum(shares ** 2)
    N = len(x)
    return 0.0 if N <= 1 else (h - 1.0 / N) / (1.0 - 1.0 / N)


def topk(x, k):
    x = np.asarray(x, float)
    s = x.sum()
    return 0.0 if s == 0 else np.sort(x)[::-1][:k].sum() / s


g   = gini(vals)
h   = hhi_norm(vals)
cv  = float(np.std(vals) / (np.mean(vals) + 1e-12))
t10 = topk(vals, 10) * 100
t50 = topk(vals, 50) * 100

with open(out_txt, "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"Gini\t{g:.6f}\n")
    f.write(f"HHI_norm\t{h:.6f}\n")
    f.write(f"CV\t{cv:.6f}\n")
    f.write(f"Top10_share_percent\t{t10:.2f}\n")
    f.write(f"Top50_share_percent\t{t50:.2f}\n")

print(
    f"[mirna_dispersion:{sample}] "
    f"Gini={g:.4f}  HHI={h:.4f}  CV={cv:.4f}  "
    f"Top10={t10:.1f}%  Top50={t50:.1f}%  n_miRNA={len(vals)}"
)
