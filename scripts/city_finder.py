#!/usr/bin/env python3
"""
city_finder.py

Scan a collapsed FASTA for reads 12–15 nt long and check whether each
read aligns anywhere along a known miRNA sequence (sliding window).

For each unique short read, every miRNA is scanned with a sliding window
of the read length.  The best match (fewest mismatches, ties broken by
earliest position) is recorded per (read, miRNA) pair.

The canonical miRNA sequence per name is the isoform with the highest
expression value in the reference CSV.

Optionally joins collapse_stats.tsv to add per-VCG read counts.

Outputs a TSV with one row per matching (unique_seq, miRNA) pair.
Output columns:
    sequence, length, read_count, mirna_name, mirna_seq,
    mirna_match_pos, match_window, mismatches

Called by Snakemake or directly:
    python scripts/city_finder.py \\
        --fasta   results/SAMPLE/SAMPLE.collapsed.fasta \\
        --mirna   mm39.brain.mirna.filtered.csv \\
        --out     results/SAMPLE/SAMPLE.city_finder.tsv \\
        [--stats  results/SAMPLE/SAMPLE.collapse_stats.tsv] \\
        [--min-len 12] [--max-len 15] [--mismatches 0]
"""

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd


# ── Mismatch counting ─────────────────────────────────────────────────────────

def count_mismatches(a: str, b: str) -> int:
    """Hamming distance between two equal-length strings."""
    return sum(x != y for x, y in zip(a, b))


# ── miRNA reference ────────────────────────────────────────────────────────────

def load_mirna_seqs(csv_path: str) -> dict:
    """
    Returns {mirna_name: full_seq} using the highest-expression isoform
    per miRNA name.  Sequences are returned in DNA form (U→T).
    """
    df = pd.read_csv(csv_path)
    required = {"RNA_seq", "RNA_name", "expression"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"miRNA CSV missing columns: {missing}")

    canonical = (
        df.sort_values("expression", ascending=False)
        .drop_duplicates(subset="RNA_name")
        .set_index("RNA_name")["RNA_seq"]
        .str.upper()
        .str.replace("U", "T", regex=False)
    )
    return canonical.to_dict()


# ── Input parsing ─────────────────────────────────────────────────────────────

def iter_fasta(path: str):
    """Yield (header, sequence) pairs from a FASTA file."""
    header = None
    seq_parts = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
    if header is not None:
        yield header, "".join(seq_parts)


def iter_fastq(path: str):
    """Yield (header, sequence) pairs from a FASTQ or gzipped FASTQ file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        while True:
            header = fh.readline().rstrip()
            if not header:
                break
            seq  = fh.readline().rstrip().upper()
            fh.readline()   # +
            fh.readline()   # quality
            yield header[1:], seq


def iter_reads(path: str):
    """Dispatch to FASTA or FASTQ parser based on file extension."""
    p = str(path).lower()
    if ".fastq" in p or ".fq" in p:
        yield from iter_fastq(path)
    else:
        yield from iter_fasta(path)


def count_unique_seqs(path: str, min_len: int, max_len: int) -> dict:
    """
    For FASTQ input (no collapse_stats), count occurrences of each unique
    sequence directly from the reads.  Returns {seq: count}.
    """
    counts: dict = {}
    for _, seq in iter_reads(path):
        if min_len <= len(seq) <= max_len:
            counts[seq] = counts.get(seq, 0) + 1
    return counts


# ── Core ───────────────────────────────────────────────────────────────────────

def best_window_match(seq: str, mirna_seq: str, max_mismatches: int) -> tuple | None:
    """
    Slide a window of len(seq) across mirna_seq.
    Returns (mismatches, pos, window) for the best match within budget,
    or None if no window qualifies.
    Ties on mismatches broken by earliest position.
    """
    L = len(seq)
    best = None
    for pos in range(len(mirna_seq) - L + 1):
        window = mirna_seq[pos : pos + L]
        mm = count_mismatches(seq, window)
        if mm <= max_mismatches:
            if best is None or mm < best[0]:
                best = (mm, pos, window)
    return best


def run_city_finder(
    fasta_path: str,
    mirna_csv: str,
    out_path: str,
    stats_path: str | None = None,
    min_len: int = 12,
    max_len: int = 15,
    max_mismatches: int = 0,
    sample: str = "",
) -> None:
    tag = f"[city_finder:{sample}]" if sample else "[city_finder]"

    print(f"{tag} Loading miRNA sequences from {mirna_csv}", flush=True)
    mirna_seqs = load_mirna_seqs(mirna_csv)
    print(f"{tag} {len(mirna_seqs)} canonical miRNAs loaded", flush=True)

    # Optional: VCG → read count from collapse_stats
    vcg_counts: dict = {}
    if stats_path and Path(stats_path).exists():
        stats = pd.read_csv(stats_path, sep="\t", dtype=str)
        if {"VCG", "reads_before"}.issubset(stats.columns):
            vcg_counts = dict(
                zip(stats["VCG"], stats["reads_before"].astype(int))
            )
            print(f"{tag} Loaded read counts for {len(vcg_counts)} VCGs", flush=True)

    is_fastq = any(x in str(fasta_path).lower() for x in (".fastq", ".fq"))
    print(f"{tag} Input type:           {'FASTQ' if is_fastq else 'FASTA'}", flush=True)
    print(f"{tag} Scanning {fasta_path} for {min_len}–{max_len} nt reads …", flush=True)

    # Collect unique sequences in the length range
    if is_fastq:
        # Count directly from FASTQ — no collapse_stats needed
        seen_seqs = count_unique_seqs(fasta_path, min_len, max_len)
        n_scanned = sum(1 for _ in iter_reads(fasta_path))
    else:
        seen_seqs = {}
        n_scanned = 0
        for header, seq in iter_fasta(fasta_path):
            n_scanned += 1
            if not (min_len <= len(seq) <= max_len):
                continue
            if seq not in seen_seqs:
                seen_seqs[seq] = vcg_counts.get(seq, None)

    n_unique = len(seen_seqs)
    print(f"{tag} Reads scanned:       {n_scanned:,}", flush=True)
    print(f"{tag} Unique seqs in range: {n_unique:,}", flush=True)

    rows = []
    for seq, count in seen_seqs.items():
        L = len(seq)
        for mirna_name, mirna_seq in mirna_seqs.items():
            if len(mirna_seq) < L:
                continue
            hit = best_window_match(seq, mirna_seq, max_mismatches)
            if hit is not None:
                mm, pos, window = hit
                rows.append({
                    "sequence":       seq,
                    "length":         L,
                    "read_count":     count,
                    "mirna_name":     mirna_name,
                    "mirna_seq":      mirna_seq,
                    "mirna_match_pos": pos,
                    "match_window":   window,
                    "mismatches":     mm,
                })

    print(f"{tag} Matches found:        {len(rows):,}", flush=True)

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    if rows:
        out_df = pd.DataFrame(rows).sort_values(
            ["mismatches", "read_count", "mirna_name"],
            ascending=[True, False, True],
        )
    else:
        out_df = pd.DataFrame(columns=[
            "sequence", "length", "read_count", "mirna_name", "mirna_seq",
            "mirna_match_pos", "match_window", "mismatches",
        ])
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"{tag} Output → {out_path}", flush=True)


# ── Entry points ───────────────────────────────────────────────────────────────

if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    _stats = snakemake.input.get("stats", None)
    if isinstance(_stats, list):
        _stats = _stats[0] if _stats else None
    run_city_finder(
        fasta_path     = snakemake.input.fasta,
        mirna_csv      = snakemake.input.mirna_csv,
        out_path       = snakemake.output.tsv,
        stats_path     = _stats,
        min_len        = snakemake.params.get("min_len", 12),
        max_len        = snakemake.params.get("max_len", 15),
        max_mismatches = snakemake.params.get("mismatches", 0),
        sample         = snakemake.wildcards.sample,
    )

elif __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Find short reads (12–15 nt) matching anywhere in miRNA sequences."
    )
    ap.add_argument("--fasta",      required=True, help="collapsed FASTA input")
    ap.add_argument("--mirna",      required=True, help="miRNA CSV (RNA_seq, RNA_name, expression)")
    ap.add_argument("--out",        required=True, help="output TSV path")
    ap.add_argument("--stats",      default=None,  help="collapse_stats.tsv (adds read counts)")
    ap.add_argument("--min-len",    type=int, default=12)
    ap.add_argument("--max-len",    type=int, default=15)
    ap.add_argument("--mismatches", type=int, default=0, help="max allowed mismatches")
    args = ap.parse_args()
    run_city_finder(
        fasta_path     = args.fasta,
        mirna_csv      = args.mirna,
        out_path       = args.out,
        stats_path     = args.stats,
        min_len        = args.min_len,
        max_len        = args.max_len,
        max_mismatches = args.mismatches,
    )
