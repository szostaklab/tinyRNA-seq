#!/usr/bin/env python3
"""
make_input_tsv.py — Stage 1 of the pipeline: auto-discover input FASTQs.

Scans a raw-reads directory (one subfolder per experiment) and writes a starting
experiments.tsv in the unified schema. You then REVIEW/EDIT the TSV (set the correct
experiment_type, organism, adaptors, etc.) before running snakemake.

Workflow:
    1. Drop your reads into raw_reads/<experiment_name>/*.fastq[.gz]
    2. python scripts/make_input_tsv.py        # writes experiments.tsv
    3. Review experiments.tsv (fix experiment_type / organism / adaptors)
    4. snakemake --cores 8

Track detection:
    - Immediate subfolders of raw_reads named vcg / dno / bio are treated as track
      groupings: raw_reads/bio/<experiment>/ → experiment_type=bio, etc.
    - Any other immediate subfolder is one experiment with experiment_type=--type.
    So you can keep everything flat (all one --type) or nest by track, or mix both.

FASTQ detection per experiment folder:
    - A file whose name contains "trimmed" → trimmed_fastq  (pre-trimmed / supplier-kit
      path; see docs/pretrimming.md). R1 preferred.
    - Otherwise the raw FASTQ → fastq_file. R1 preferred.

Paths are written RELATIVE to the current directory by default (portable across
machines); pass --absolute for absolute paths.

Usage:
    python scripts/make_input_tsv.py
    python scripts/make_input_tsv.py --type vcg --vcg references/synthetic/VCG_Oligos_Labelled.csv
    python scripts/make_input_tsv.py --type bio --organism mouse
    python scripts/make_input_tsv.py --raw-dir /path/to/reads --absolute
"""

import argparse
import os
import re
from pathlib import Path

# ── Defaults (all overridable per row in the TSV afterwards) ───────────────────

DEF_RAW_DIR     = "raw_reads"
DEF_TYPE        = "bio"                         # experiment_type for un-nested folders
DEF_ORGANISM    = "mouse"                       # applied to bio rows
DEF_VCG         = "references/synthetic/VCG_Oligos_Labelled.csv"     # VCG reference (code,sequence,strand); vcg/dno rows
DEF_ADAPTOR_3P  = "AGATCGGAAGAGCACACGTCT"
DEF_READ_LAYOUT = "FOUR_N_MADAP"
DEF_PHRED       = 30

TRACK_DIRS = {"vcg", "dno", "bio"}

COLUMNS = [
    "experiment_name", "experiment_type", "organism",
    "fastq_file", "trimmed_fastq", "sam_path", "sorted_bam",
    "read_layout", "adaptor_3p", "vcg_list_path", "phred_cutoff",
    "cofold_analysis", "cofold_temperature", "adaptor_5p_ext", "adaptor_3p_tail",
]

R1_RE      = re.compile(r"R1", re.IGNORECASE)
TRIMMED_RE = re.compile(r"trim", re.IGNORECASE)

FASTQ_GLOBS = ("*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz")


# ── Helpers ────────────────────────────────────────────────────────────────────

def _fmt_path(p: Path, absolute: bool) -> str:
    ap = Path(p).resolve()
    return str(ap) if absolute else os.path.relpath(ap, start=os.getcwd())


def _pick(cands, r1_only=False):
    """Prefer an R1 file, else the first (alphabetical)."""
    if not cands:
        return None
    r1 = [p for p in cands if R1_RE.search(p.name)]
    if r1:
        return r1[0]
    return None if r1_only else cands[0]


def find_fastqs(folder: Path):
    """Return (fastq_file, trimmed_fastq) Path|None for an experiment folder.

    A trimmed FASTQ takes precedence: its presence signals the pre-trimmed / supplier-kit
    workflow, so we fill trimmed_fastq and leave fastq_file blank (the Snakefile routes on
    fastq_file first, so filling both would wrongly force the raw path).
    """
    all_fq = sorted(
        (p for g in FASTQ_GLOBS for p in folder.glob(g)),
        key=lambda p: p.name.lower(),
    )
    trimmed = [p for p in all_fq if TRIMMED_RE.search(p.name)]
    if trimmed:
        return None, _pick(trimmed)
    return _pick(all_fq), None


def _iter_experiments(raw_dir: Path, default_type: str):
    """Yield (experiment_folder, experiment_type) honouring vcg/dno/bio track dirs."""
    for child in sorted(p for p in raw_dir.iterdir() if p.is_dir()):
        if child.name.lower() in TRACK_DIRS:
            for exp in sorted(p for p in child.iterdir() if p.is_dir()):
                yield exp, child.name.lower()
        else:
            yield child, default_type


def build_rows(raw_dir, default_type, organism, vcg_path, adaptor_3p,
               read_layout, phred, absolute) -> list[dict]:
    rows = []
    for folder, etype in _iter_experiments(raw_dir, default_type):
        raw_fq, trimmed_fq = find_fastqs(folder)
        is_bio = etype == "bio"
        rows.append({
            "experiment_name": folder.name,
            "experiment_type": etype,
            "organism":        organism if is_bio else "",
            "fastq_file":      _fmt_path(raw_fq, absolute) if raw_fq else "",
            "trimmed_fastq":   _fmt_path(trimmed_fq, absolute) if trimmed_fq else "",
            "sam_path":        "",
            "sorted_bam":      "",
            "read_layout":     read_layout,
            "adaptor_3p":      adaptor_3p,
            "vcg_list_path":   ("" if is_bio else (_fmt_path(Path(vcg_path), absolute)
                                                   if vcg_path and Path(vcg_path).exists() else vcg_path)),
            "cofold_analysis": "no",
            "cofold_temperature": "",
            "phred_cutoff":    phred,
            "adaptor_5p_ext":  "",
            "adaptor_3p_tail": "",
        })
    return rows


def write_tsv(rows, out: Path) -> None:
    with open(out, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")


# ── Main ────────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description="Stage 1: discover FASTQs → experiments.tsv (unified schema).")
    ap.add_argument("--raw-dir", default=DEF_RAW_DIR,
                    help=f"Directory of per-experiment FASTQ folders (default: {DEF_RAW_DIR})")
    ap.add_argument("-o", "--out", default="experiments.tsv", help="Output TSV (default: experiments.tsv)")
    ap.add_argument("--type", default=DEF_TYPE, choices=sorted(TRACK_DIRS),
                    help=f"experiment_type for folders not under a vcg/dno/bio track dir (default: {DEF_TYPE})")
    ap.add_argument("--organism", default=DEF_ORGANISM, help=f"organism for bio rows (default: {DEF_ORGANISM})")
    ap.add_argument("--vcg", default=DEF_VCG, help=f"vcg_list_path for vcg/dno rows (default: {DEF_VCG})")
    ap.add_argument("--adaptor-3p", default=DEF_ADAPTOR_3P, help=f"default 3' adaptor (default: {DEF_ADAPTOR_3P})")
    ap.add_argument("--read-layout", default=DEF_READ_LAYOUT, help=f"default read_layout (default: {DEF_READ_LAYOUT})")
    ap.add_argument("--phred", type=int, default=DEF_PHRED, help=f"default phred_cutoff (default: {DEF_PHRED})")
    ap.add_argument("--absolute", action="store_true", help="write absolute paths instead of relative")
    ap.add_argument("-f", "--force", action="store_true", help="overwrite an existing output TSV")
    args = ap.parse_args()

    raw_dir = Path(args.raw_dir)
    if not raw_dir.is_dir():
        raise SystemExit(f"Raw reads directory not found: {raw_dir}")

    out = Path(args.out)
    if out.exists() and not args.force:
        raise SystemExit(f"{out} already exists — use --force to overwrite (you may have edited it).")

    rows = build_rows(
        raw_dir=raw_dir, default_type=args.type, organism=args.organism,
        vcg_path=args.vcg, adaptor_3p=args.adaptor_3p, read_layout=args.read_layout,
        phred=args.phred, absolute=args.absolute,
    )
    if not rows:
        raise SystemExit(f"No experiment subfolders found in {raw_dir}")

    missing = [r["experiment_name"] for r in rows if not r["fastq_file"] and not r["trimmed_fastq"]]
    write_tsv(rows, out)

    print(f"Wrote {len(rows)} experiment(s) to {out}")
    if missing:
        print(f"  WARNING: no FASTQ found for: {', '.join(missing)}")
    print("  NEXT: review experiments.tsv — check experiment_type, organism, read_layout,")
    print(f"        and adaptor_3p per row (defaults: type={args.type}, adaptor_3p={args.adaptor_3p}).")
    print("  Then run:  snakemake --cores 8 -n   (dry run)   and   snakemake --cores 8")


if __name__ == "__main__":
    main()
