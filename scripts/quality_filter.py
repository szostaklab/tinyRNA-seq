#!/usr/bin/env python3
"""
quality_filter.py

Filter reads by average Phred quality score. Writes high- and low-quality
reads to separate gzipped FASTQ files.

Called by Snakemake (rule quality_filter) or directly:
    python scripts/quality_filter.py \\
        --input raw.fastq[.gz] \\
        --high high_quality.fastq.gz \\
        --low  low_quality.fastq.gz \\
        --phred 30
"""

import argparse
import gzip
import math
import os
import sys
from pathlib import Path

from Bio import SeqIO


# ── Helpers ───────────────────────────────────────────────────────────────────

def average_phred(scores: list[int]) -> float:
    """Convert per-base Phred scores to a single average Phred score."""
    n = len(scores)
    mean_p = sum(10 ** (-q / 10) for q in scores) / n
    return -10 * math.log10(mean_p)


def approx_num_reads(path: str, sample_lines: int = 40) -> int:
    """Estimate read count by sampling average line length vs file size."""
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt") as fh:
            lengths = [len(fh.readline()) for _ in range(sample_lines)]
        avg = sum(lengths) / len(lengths)
        size = os.path.getsize(path)
        return max(1, round(size / avg / 4))
    except Exception:
        return 0


def open_fastq(path: str):
    """Open a plain or gzipped FASTQ for reading."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


# ── Core ──────────────────────────────────────────────────────────────────────

def filter_reads(
    input_path: str,
    high_path: str,
    low_path: str,
    phred_cutoff: int,
    sample: str = "",
) -> dict:
    """
    Split reads in input_path into high/low quality FASTQ.gz files.
    Returns a summary dict with read counts.
    """
    n_total = n_high = n_low = 0
    with (
        open_fastq(input_path) as fh_in,
        gzip.open(high_path, "wt") as fh_high,
        gzip.open(low_path, "wt") as fh_low,
    ):
        for record in SeqIO.parse(fh_in, "fastq"):
            n_total += 1
            scores = record.letter_annotations["phred_quality"]
            if average_phred(scores) >= phred_cutoff:
                n_high += 1
                SeqIO.write(record, fh_high, "fastq")
            else:
                n_low += 1
                SeqIO.write(record, fh_low, "fastq")

    return {
        "total":    n_total,
        "high":     n_high,
        "low":      n_low,
        "pct_high": f"{n_high / n_total:.2%}" if n_total else "N/A",
        "pct_low":  f"{n_low  / n_total:.2%}" if n_total else "N/A",
    }


# ── Entry points ──────────────────────────────────────────────────────────────

def run(input_path: str, high_path: str, low_path: str, phred: int, sample: str = "") -> None:
    Path(high_path).parent.mkdir(parents=True, exist_ok=True)
    Path(low_path).parent.mkdir(parents=True, exist_ok=True)

    tag = f"[quality_filter:{sample}]" if sample else "[quality_filter]"
    print(f"{tag} Input:  {input_path}")
    print(f"{tag} Phred cutoff: {phred}")

    summary = filter_reads(input_path, high_path, low_path, phred, sample=sample)

    print(f"{tag} Total reads:      {summary['total']:,}")
    print(f"{tag} High quality:     {summary['high']:,}  ({summary['pct_high']})")
    print(f"{tag} Low quality:      {summary['low']:,}  ({summary['pct_low']})")
    print(f"{tag} High → {high_path}")
    print(f"{tag} Low  → {low_path}")


# Snakemake entry point
if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    run(
        input_path=snakemake.input.fastq,
        high_path=snakemake.output.high,
        low_path=snakemake.output.low,
        phred=int(snakemake.params.phred),
        sample=snakemake.wildcards.sample,
    )

# CLI entry point
elif __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Filter FASTQ reads by average Phred score.")
    ap.add_argument("--input",  required=True, help="Input FASTQ (.fastq or .fastq.gz)")
    ap.add_argument("--high",   required=True, help="Output high-quality FASTQ.gz path")
    ap.add_argument("--low",    required=True, help="Output low-quality FASTQ.gz path")
    ap.add_argument("--phred",  type=int, default=30, help="Phred cutoff (default: 30)")
    args = ap.parse_args()
    run(args.input, args.high, args.low, args.phred)
