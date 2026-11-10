"""
apply_hierarchy.py

Hierarchical conflict resolution for featureCounts output.

Reads the featureCounts -R CORE per-read assignment file and gene_info TSV.
For each assigned read, looks up the gene_type and applies a priority tier.
For reads assigned to multiple features (multi-overlap via -O), only features
at the highest-priority tier are kept; count is split equally among tied winners.

Output: two-column TSV (gene_id, hier_count) — same format as
allRNA.gene.slim.tsv but with hierarchically resolved counts.

Tier priority:
  0: miRNA
  1: rRNA, Mt_rRNA
  2: tRNA, Mt_tRNA
  3: piRNA
  4: snoRNA, scaRNA, snRNA
  5: misc_RNA, Y_RNA, vaultRNA, scRNA, srpRNA
  6: RepeatMasker-derived (SINE, LINE, LTR, DNA, Satellite, Simple_repeat, Retroposon,
     RC, Low_complexity, Unknown, Unspecified)
  7: protein_coding and everything else

Notes:
  - With featureCounts CORE format, multi-overlap reads may appear as multiple
    rows with the same ReadName (one row per feature). This script collects all
    assignments per read before applying hierarchy.
  - If a read has a single assignment (the common case with -k 1 bowtie), it
    passes through unchanged regardless of tier.
  - Run on results/{sample}/{sample}.NH.sorted.bam.featureCounts
"""

import sys
from collections import defaultdict

TIERS = [
    {"miRNA"},
    {"rRNA", "Mt_rRNA"},
    {"tRNA", "Mt_tRNA"},
    {"piRNA"},
    {"snoRNA", "scaRNA", "snRNA"},
    {"misc_RNA", "Y_RNA", "vaultRNA", "scRNA", "srpRNA"},   # real ncRNA (incl. new scRNA/srpRNA)
    # RepeatMasker-derived classes — reported, but lowest priority so a read overlapping a
    # real ncRNA above and a repeat here goes to the ncRNA (fixes misc_RNA/scRNA dilution).
    {"SINE", "LINE", "LTR", "DNA", "Satellite", "Simple_repeat", "Retroposon", "RC",
     "Low_complexity", "Unknown", "Unspecified"},
]


def gene_tier(gene_type: str) -> int:
    for i, tier_set in enumerate(TIERS):
        if gene_type in tier_set:
            return i
    return len(TIERS)


def load_gene_info(path: str) -> dict:
    """Return {gene_id: gene_type} from gene_info.mm39.tsv (gene_id, gene_name, gene_type)."""
    d = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                d[parts[0]] = parts[2]
    return d


def load_assignments(path: str) -> dict:
    """
    Return {read_name: [gene_id, ...]} for Assigned reads.
    Collects ALL gene_ids per read to handle multi-overlap rows.
    """
    assignments = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            if parts[0] in ("ReadName",) or (len(parts) > 1 and parts[1] == "Status"):
                continue
            if parts[1] == "Assigned":
                # featureCounts -R CORE puts multi-overlap targets comma-separated in
                # the last column (e.g. "GENE_A,GENE_B"); split them into separate genes.
                for gene in parts[-1].split(","):
                    gene = gene.strip()
                    if gene and gene not in assignments[parts[0]]:
                        assignments[parts[0]].append(gene)
    return dict(assignments)


def apply_hierarchy(assignments: dict, gene_info: dict) -> dict:
    """
    For each read, keep only the features at the best (lowest) tier.
    Split count equally among tied winners.
    Returns {gene_id: float_count}.
    """
    counts = defaultdict(float)
    n_single = 0
    n_multi = 0
    n_demoted = 0

    for read_name, genes in assignments.items():
        if len(genes) == 1:
            counts[genes[0]] += 1.0
            n_single += 1
            continue

        n_multi += 1
        tiers = [gene_tier(gene_info.get(g, "")) for g in genes]
        best = min(tiers)
        kept = [g for g, t in zip(genes, tiers) if t == best]
        if len(kept) < len(genes):
            n_demoted += 1
        weight = 1.0 / len(kept)
        for g in kept:
            counts[g] += weight

    return dict(counts), n_single, n_multi, n_demoted


if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")

    gene_info = load_gene_info(snakemake.input.gene_info)
    print(f"[apply_hierarchy] Loaded {len(gene_info):,} gene type entries", flush=True)

    assignments = load_assignments(snakemake.input.read_assign)
    print(f"[apply_hierarchy] Loaded {len(assignments):,} assigned reads from CORE file", flush=True)

    counts, n_single, n_multi, n_demoted = apply_hierarchy(assignments, gene_info)

    print(f"[apply_hierarchy] Single-assignment reads: {n_single:,}", flush=True)
    print(f"[apply_hierarchy] Multi-assignment reads:  {n_multi:,}", flush=True)
    print(f"[apply_hierarchy] Reads with tier demotion: {n_demoted:,}", flush=True)
    print(f"[apply_hierarchy] Genes with counts: {len(counts):,}", flush=True)

    with open(snakemake.output[0], "w") as fh:
        for gene_id, count in sorted(counts.items()):
            fh.write(f"{gene_id}\t{count:.6f}\n")

    print(f"[apply_hierarchy] Written → {snakemake.output[0]}", flush=True)
