#!/usr/bin/env python3
"""
add_sequences.py

Appends five columns to counts.smallRNA.withCPM.tsv and adds a header row:
  ref_seq     — canonical reference RNA sequence (T→U) from genome FASTA + GTF
  top_bam_seq — most abundant read sequence in BAM for each gene (DNA)
  top_bam_pct — % of BAM reads for that gene carrying top_bam_seq
  top_vcg     — most abundant VCG from parsed_reads.tsv (FASTQ-start only; NA otherwise)
  top_vcg_pct — % of pre-dedup reads for that gene carrying top_vcg (FASTQ-start only; NA otherwise)

For SAM/BAM-start samples (no parsed_reads.tsv), top_vcg and top_vcg_pct are NA.
For FASTQ-start samples, top_vcg counts are from the pre-UMI-collapse parsed_reads.tsv;
  top_bam_seq/pct are from the post-collapse BAM.

Requires:
  - samtools on PATH (for view + faidx)
  - genome FASTA indexed with samtools faidx (genome.fa.fai must exist)
  - featureCounts run with -R flag to produce the read-assignment file

Called by Snakemake rule add_sequences.
"""

import gzip
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


# ── GTF parsing ───────────────────────────────────────────────────────────────

def parse_gtf_gene_coords(gtf_path: str) -> dict:
    """
    Return {gene_id: (chrom, start_0based, end_1based, strand)}.
    Tries 'gene' features first; falls back to first 'exon' hit per gene_id
    for GTFs that only have exon-level entries (tRNA, piRNA, RepeatMasker).
    """
    primary: dict = {}
    fallback: dict = {}
    opener = gzip.open if str(gtf_path).endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feat = fields[2]
            if feat not in ("gene", "exon"):
                continue
            m = re.search(r'gene_id "([^"]+)"', fields[8])
            if not m:
                continue
            gid = m.group(1)
            entry = (fields[0], int(fields[3]) - 1, int(fields[4]), fields[6])
            if feat == "gene":
                primary[gid] = entry
            elif gid not in fallback:
                fallback[gid] = entry
    # merge: gene features take priority
    merged = {**fallback, **primary}
    return merged


# ── Genome sequence extraction ────────────────────────────────────────────────

_RC = str.maketrans("ACGTacgt", "TGCAtgca")


def _revcomp(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def fetch_ref_seq(genome_fasta: str, chrom: str, start_0: int, end_1: int, strand: str) -> str:
    """Extract sequence with samtools faidx. Returns RNA string (T→U)."""
    region = f"{chrom}:{start_0 + 1}-{end_1}"   # samtools faidx is 1-based inclusive
    result = subprocess.run(
        ["samtools", "faidx", genome_fasta, region],
        capture_output=True, text=True,
    )
    if result.returncode != 0 or not result.stdout.strip():
        return ""
    lines = result.stdout.strip().split("\n")
    seq = "".join(lines[1:]).upper()
    if strand == "-":
        seq = _revcomp(seq)
    return seq.replace("T", "U")


# ── featureCounts -R assignment file ─────────────────────────────────────────

def read_fc_assignments(assign_path: str) -> dict:
    """
    Return {read_name: [gene_id, ...]} for Assigned reads.
    featureCounts -R output format (tab-separated, lines starting with '#' are comments):
      ReadName  Status  NumFeatures  Gene(s)
    Gene(s) may be comma-separated for multi-overlap reads.
    """
    assignments: dict = {}
    with open(assign_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            # Skip header row if present
            if parts[0] == "ReadName" or parts[1] == "Status":
                continue
            if parts[1] == "Assigned" and parts[0] not in assignments:
                assignments[parts[0]] = parts[3].split(",")
    return assignments


# ── BAM read sequences ────────────────────────────────────────────────────────

def read_bam_seqs(bam_path: str) -> dict:
    """
    Return {read_name: seq} from BAM (primary alignments only).
    Uses samtools view; first occurrence of each read name is kept.
    """
    seqs: dict = {}
    result = subprocess.run(
        ["samtools", "view", "-F", "260", bam_path],  # exclude unmapped + secondary
        capture_output=True, text=True,
    )
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) >= 10 and fields[0] not in seqs:
            seqs[fields[0]] = fields[9]
    return seqs


# ── Gene-level sequence aggregation ──────────────────────────────────────────

def build_gene_seq_counters(assignments: dict, bam_seqs: dict) -> dict:
    """Return {gene_id: Counter({seq: count})} from post-dedup BAM reads."""
    counters: dict = defaultdict(Counter)
    for rname, gids in assignments.items():
        seq = bam_seqs.get(rname, "")
        if seq and seq != "*":
            for gid in gids:
                counters[gid][seq] += 1
    return dict(counters)


def top_seq_pct(counter: Counter) -> tuple:
    """Return (top_seq, pct_rounded) for the most abundant sequence in counter."""
    if not counter:
        return ("", 0.0)
    total = sum(counter.values())
    seq, cnt = counter.most_common(1)[0]
    return seq, round(100.0 * cnt / total, 2)


# ── Pre-dedup VCG counts (FASTQ-start only) ───────────────────────────────────

def build_vcg_gene_data(parsed_path: str, gene_seq_counters: dict) -> dict:
    """
    Return {gene_id: (top_vcg, top_vcg_pct)} using pre-collapse VCG counts.

    Strategy:
      1. Count all VCG occurrences in parsed_reads.tsv (includes PCR duplicates).
      2. For each gene, the set of VCG sequences is taken from the post-collapse BAM
         (gene_seq_counters), since BAM_seq == VCG for FASTQ-start samples.
      3. For each VCG in a gene's set, look up its raw (pre-dedup) count.
      4. Report the VCG with the highest raw count and its % of the gene's total raw reads.
    """
    vcg_raw: Counter = Counter()
    df = pd.read_csv(parsed_path, sep="\t", dtype=str, usecols=["VCG"])
    for vcg in df["VCG"].dropna():
        vcg_raw[vcg] += 1

    result: dict = {}
    for gid, seq_counter in gene_seq_counters.items():
        gene_vcgs = list(seq_counter.keys())
        gene_raw = {v: vcg_raw.get(v, 0) for v in gene_vcgs}
        total_raw = sum(gene_raw.values())
        if total_raw == 0:
            result[gid] = ("", 0.0)
            continue
        top_vcg = max(gene_raw, key=lambda v: gene_raw[v])
        pct = round(100.0 * gene_raw[top_vcg] / total_raw, 2)
        result[gid] = (top_vcg, pct)
    return result


# ── Main ──────────────────────────────────────────────────────────────────────

def main(
    counts_path: str,
    bam_path: str,
    assign_path: str,
    gtf_path: str,
    genome_fasta: str,
    parsed_path,          # str | None
    out_path: str,
    strip_versions: bool = False,
) -> None:
    tag = "[add_sequences]"

    # Load the core counts/CPM table (headered: gene_id count CPM gene_name gene_type)
    df = pd.read_csv(counts_path, sep="\t", dtype={"gene_name": str, "gene_type": str})
    print(f"{tag} Loaded {len(df)} genes from {counts_path}")

    # ── BAM-derived sequence counters ────────────────────────────────────────
    print(f"{tag} Reading featureCounts assignments …")
    assignments = read_fc_assignments(assign_path)
    print(f"{tag}   Assigned reads: {len(assignments):,}")

    print(f"{tag} Reading BAM sequences …")
    bam_seqs = read_bam_seqs(bam_path)
    print(f"{tag}   Unique read names in BAM: {len(bam_seqs):,}")

    gene_counters = build_gene_seq_counters(assignments, bam_seqs)
    print(f"{tag} Sequence counters built for {len(gene_counters):,} genes")

    top_bam = {gid: top_seq_pct(gene_counters.get(gid, Counter())) for gid in df["gene_id"]}
    df["top_bam_seq"] = df["gene_id"].map(lambda g: top_bam[g][0])
    df["top_bam_pct"] = df["gene_id"].map(lambda g: top_bam[g][1])

    # ── VCG columns (FASTQ-start only) ───────────────────────────────────────
    if parsed_path:
        print(f"{tag} Reading parsed_reads.tsv for pre-dedup VCG counts …")
        vcg_data = build_vcg_gene_data(parsed_path, gene_counters)
        df["top_vcg"]     = df["gene_id"].map(lambda g: vcg_data.get(g, ("", 0.0))[0])
        df["top_vcg_pct"] = df["gene_id"].map(lambda g: vcg_data.get(g, ("", 0.0))[1])
    else:
        df["top_vcg"]     = pd.NA
        df["top_vcg_pct"] = pd.NA

    # ── Reference sequences from genome + GTF ────────────────────────────────
    print(f"{tag} Parsing GTF gene coordinates …")
    gene_coords = parse_gtf_gene_coords(gtf_path)

    if strip_versions:
        gene_coords = {re.sub(r"\.[0-9]+$", "", k): v for k, v in gene_coords.items()}

    print(f"{tag} Extracting reference sequences (this may take a moment) …")
    ref_seqs: dict = {}
    n_found = 0
    for gid in df["gene_id"]:
        coords = gene_coords.get(gid)
        if coords:
            ref_seqs[gid] = fetch_ref_seq(genome_fasta, *coords)
            n_found += 1
        else:
            ref_seqs[gid] = ""
    print(f"{tag}   Reference seqs found: {n_found}/{len(df)}")

    df["ref_seq"] = df["gene_id"].map(ref_seqs)

    # ── Write with header ────────────────────────────────────────────────────
    cols = [
        "gene_id", "count", "CPM", "gene_name", "gene_type",
        "ref_seq", "top_bam_seq", "top_bam_pct", "top_vcg", "top_vcg_pct",
    ]
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    df[cols].to_csv(out_path, sep="\t", index=False)
    print(f"{tag} Written → {out_path}")


# ── Entry points ─────────────────────────────────────────────────────────────

if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    parsed_raw = snakemake.input.get("parsed", [])
    if not parsed_raw:
        parsed_path = None
    elif isinstance(parsed_raw, str):
        parsed_path = parsed_raw
    else:
        parsed_path = str(parsed_raw[0])
    main(
        counts_path    = snakemake.input.counts,
        bam_path       = snakemake.input.bam,
        assign_path    = snakemake.input.read_assign,
        gtf_path       = snakemake.input.gtf,
        genome_fasta   = snakemake.input.genome,
        parsed_path    = parsed_path,
        out_path       = snakemake.output[0],
        strip_versions = str(snakemake.params.get("strip_versions", False)).lower()
                         in ("true", "1", "yes"),
    )

elif __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Add ref/BAM/VCG sequence columns to withCPM.tsv.")
    ap.add_argument("--counts",        required=True)
    ap.add_argument("--bam",           required=True)
    ap.add_argument("--read-assign",   required=True, dest="read_assign")
    ap.add_argument("--gtf",           required=True)
    ap.add_argument("--genome",        required=True)
    ap.add_argument("--parsed",        default=None)
    ap.add_argument("--output",        required=True)
    ap.add_argument("--strip-versions", action="store_true", dest="strip_versions")
    args = ap.parse_args()
    main(
        args.counts, args.bam, args.read_assign, args.gtf, args.genome,
        args.parsed, args.output, args.strip_versions,
    )
