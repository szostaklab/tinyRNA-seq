#!/usr/bin/env python3
"""
assign_piRNA.py

Post-featureCounts rescue step for piRNA.

Reads that aligned to the genome but received no feature assignment
(Unassigned_NoFeatures in the featureCounts per-read file) are extracted
from the BAM and their sequences checked against one or more piRNA reference
FASTAs (e.g. mmu.gold.fa.gz, piRNAdb.mmu.v1_7_6.fa.gz).

Matched reads are counted per piRNA name and written in with_meta.tsv format:
    gene_id  count  name  biotype
so they slot directly into the existing filter_smallrna step without any
changes to the downstream pipeline.

No double-counting risk: by construction, Unassigned_NoFeatures reads have
no overlap with anything already in with_meta.tsv.

Called by Snakemake or directly:
    python scripts/assign_piRNA.py \\
        --bam   results/SAMPLE/SAMPLE.NH.sorted.bam \\
        --fc    results/SAMPLE/SAMPLE.NH.sorted.bam.featureCounts \\
        --refs  mmu.gold.fa.gz piRNAdb.mmu.v1_7_6.fa.gz \\
        --out   results/SAMPLE/SAMPLE.pirna_matched.tsv \\
        [--sample SAMPLE]
"""

import argparse
import gzip
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path


# ── Reference loading ─────────────────────────────────────────────────────────

def load_ref(path: str) -> dict:
    """Return {sequence_upper: piRNA_name} from a (gzipped) FASTA."""
    ref = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        name = None
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:]
            elif name:
                seq = line.upper().replace("U", "T")
                if seq not in ref:          # first entry wins on collision
                    ref[seq] = name
    return ref


def merge_refs(paths: list) -> dict:
    """Merge multiple reference FASTAs; first file takes priority on collision."""
    merged = {}
    for p in paths:
        for seq, name in load_ref(p).items():
            if seq not in merged:
                merged[seq] = name
    return merged


# ── featureCounts per-read file ───────────────────────────────────────────────

def candidate_read_names(fc_path: str) -> set:
    """Return read names that are candidates for piRNA rescue:
    - Unassigned_NoFeatures, OR
    - Assigned exclusively to RM_* (RepeatMasker) features with no overlap
      with any annotated gene — these are the reads swallowed by repeat
      annotations that would otherwise be piRNAs.
    """
    names = set()
    with open(fc_path) as fh:
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) < 2:
                continue
            status = parts[1]
            if status == "Unassigned_NoFeatures":
                names.add(parts[0])
            elif status == "Assigned" and len(parts) >= 4:
                features = parts[3].split(",")
                if all(f.startswith("RM_") for f in features if f):
                    names.add(parts[0])
    return names


# ── BAM sequence extraction ───────────────────────────────────────────────────

_RC = str.maketrans("ACGT", "TGCA")

def revcomp(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def extract_sequences(bam_path: str, read_names: set) -> dict:
    """
    Use samtools view -N to pull sequences for a set of read names.
    Reverse-complements reads with FLAG & 16 so sequences match the
    original strand (consistent with the collapsed FASTA and reference sets).
    Returns {read_name: sequence}.
    """
    seqs = {}
    with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as tmp:
        tmp.write("\n".join(read_names))
        tmp_path = tmp.name

    cmd = ["samtools", "view", "-N", tmp_path, bam_path]
    result = subprocess.run(cmd, capture_output=True, text=True)
    Path(tmp_path).unlink(missing_ok=True)

    if result.returncode != 0:
        raise RuntimeError(f"samtools view failed:\n{result.stderr}")

    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 10:
            rname = fields[0]
            flag  = int(fields[1])
            seq   = fields[9].upper()
            if flag & 16:
                seq = revcomp(seq)
            if rname not in seqs:
                seqs[rname] = seq
    return seqs


# ── Core ──────────────────────────────────────────────────────────────────────

def run(
    bam_path: str,
    fc_path: str,
    ref_paths: list,
    out_path: str,
    sample: str = "",
    pirna_min_len: int = 24,
    pirna_max_len: int = 32,
) -> None:
    tag = f"[assign_piRNA:{sample}]" if sample else "[assign_piRNA]"

    print(f"{tag} Loading piRNA references …", flush=True)
    ref = merge_refs(ref_paths)
    print(f"{tag} {len(ref):,} unique piRNA sequences loaded from {len(ref_paths)} file(s)",
          flush=True)
    print(f"{tag} piRNA length gate: {pirna_min_len}–{pirna_max_len} nt", flush=True)

    print(f"{tag} Parsing candidate reads from {fc_path} …", flush=True)
    candidates = candidate_read_names(fc_path)
    print(f"{tag} {len(candidates):,} candidate reads "
          f"(Unassigned_NoFeatures + RM_-only assigned)", flush=True)

    if not candidates:
        print(f"{tag} Nothing to rescue — writing empty output.", flush=True)
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)
        Path(out_path).write_text("")
        return

    print(f"{tag} Extracting sequences from BAM …", flush=True)
    seqs = extract_sequences(bam_path, candidates)
    print(f"{tag} {len(seqs):,} sequences retrieved", flush=True)

    n_length_fail = 0
    counts: dict = defaultdict(int)
    for seq in seqs.values():
        if not (pirna_min_len <= len(seq) <= pirna_max_len):
            n_length_fail += 1
            continue
        if seq in ref:
            counts[ref[seq]] += 1
    print(f"{tag} {n_length_fail:,} reads excluded by length gate", flush=True)

    n_matched = sum(counts.values())
    print(f"{tag} {n_matched:,} reads matched across {len(counts):,} unique piRNAs",
          flush=True)

    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        for pirna_name, count in sorted(counts.items(), key=lambda x: -x[1]):
            # Format: gene_id  count  name  biotype  (matches with_meta.tsv)
            fh.write(f"{pirna_name}\t{count}\t{pirna_name}\tpiRNA\n")

    print(f"{tag} Output → {out_path}", flush=True)


# ── Entry points ──────────────────────────────────────────────────────────────

if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    run(
        bam_path      = snakemake.input.bam,
        fc_path       = snakemake.input.fc,
        ref_paths     = snakemake.input.refs,
        out_path      = snakemake.output[0],
        sample        = snakemake.wildcards.sample,
        pirna_min_len = snakemake.params.get("pirna_min_len", 24),
        pirna_max_len = snakemake.params.get("pirna_max_len", 32),
    )

elif __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Rescue unassigned piRNA reads post-featureCounts."
    )
    ap.add_argument("--bam",    required=True,        help="NH.sorted.bam")
    ap.add_argument("--fc",     required=True,        help="*.bam.featureCounts per-read file")
    ap.add_argument("--refs",   required=True, nargs="+", help="piRNA reference FASTA(s)")
    ap.add_argument("--out",    required=True,        help="output TSV path")
    ap.add_argument("--sample", default="",           help="sample name for log tags")
    args = ap.parse_args()
    run(args.bam, args.fc, args.refs, args.out, args.sample)
