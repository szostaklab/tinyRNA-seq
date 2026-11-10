"""Reference-file resolution and startup validation for the bio (genome-alignment) track.

Reference files are not shipped with the repo; the user downloads them and points the
`organisms:` block in config.yaml at their location (see references/README.md). This module
resolves an organism name to its reference paths and fails fast, before any job runs, if a
required file is missing or a genome FASTA is not indexed.

Used two ways:
  1. From the Snakefile at parse time:
         from references import validate_references, organism_refs
         validate_references(exp_df, config)          # raises SystemExit if anything missing
  2. Standalone sanity check:
         python scripts/references.py                 # uses config.yaml + experiments.tsv
         python scripts/references.py --organism mouse
"""

from __future__ import annotations

import glob
from pathlib import Path

# Files a bio sample always needs, and the optional ones (feature disabled if key absent).
BOWTIE_SUFFIXES = [".1", ".2", ".3", ".4", ".rev.1", ".rev.2"]   # <prefix><suffix>.ebwt
OPTIONAL_KEYS   = ("pirna_refs", "mirna_csv")


def resolve_gtf_path(configured):
    """Locate the GTF by presence, not by exact filename.

    The `gtf:` value in config.yaml points at a *location*, not a required name — users
    may call their annotation anything. Resolution order:
      - if it is an existing file → use it;
      - else look in its directory (or the directory it names) for a single .gtf / .gtf.gz
        and use that.

    Returns (path, None) on success, or (None, message) describing exactly what is wrong
    (not set / directory missing / none found / several found → disambiguate).
    """
    if not configured:
        return None, "gtf: not set in config.yaml"
    p = Path(configured)
    if p.is_file():
        return str(p), None
    search = p if p.is_dir() else p.parent
    if not search.is_dir():
        return None, f"gtf directory not found: {search}/ (config gtf: {configured})"
    hits = sorted(glob.glob(str(search / "*.gtf")) + glob.glob(str(search / "*.gtf.gz")))
    if len(hits) == 1:
        return hits[0], None
    if not hits:
        return None, f"no .gtf / .gtf.gz file found in {search}/ (config gtf: {configured})"
    names = ", ".join(Path(h).name for h in hits)
    return None, (f"multiple GTF files in {search}/ — set 'gtf:' to the exact file. Found: {names}")


def organism_refs(organism: str, config: dict) -> dict:
    """Return the reference dict for `organism`, or raise KeyError with valid names listed."""
    organisms = config.get("organisms", {}) or {}
    if organism not in organisms:
        valid = ", ".join(sorted(organisms)) or "(none defined)"
        raise KeyError(
            f"organism '{organism}' is not defined in config.yaml under 'organisms:'. "
            f"Defined organisms: {valid}"
        )
    return organisms[organism] or {}


def _bowtie_index_missing(prefix: str) -> list[str]:
    """Return the list of missing bowtie1 index files for `prefix` (empty if complete).

    Accepts both small (.ebwt) and large (.ebwtl) index formats: an index is complete if
    every one of the 6 index parts exists in at least one of the two extensions.
    """
    missing = []
    for suf in BOWTIE_SUFFIXES:
        small = Path(f"{prefix}{suf}.ebwt")
        large = Path(f"{prefix}{suf}.ebwtl")
        if not small.exists() and not large.exists():
            missing.append(f"{prefix}{suf}.ebwt")
    return missing


def _check_organism(organism: str, refs: dict,
                    need_pirna: bool, need_mirna: bool, need_genome: bool) -> list[str]:
    """Return a list of human-readable problems for one organism (empty == all good).

    Only the reference sets the enabled modules actually use are required:
      bowtie_index + gtf  → always (core alignment + counting)
      genome_fasta        → only if sequence_annotation module is on (need_genome)
      pirna_refs          → only if pirna_matching module is on (need_pirna)
      mirna_csv           → only if city_finder or comparisons module is on (need_mirna)
    """
    problems: list[str] = []

    # ── Always required: bowtie index (prefix -> 6 files) ──
    idx = refs.get("bowtie_index")
    if not idx:
        problems.append("bowtie_index: not set in config.yaml")
    else:
        for f in _bowtie_index_missing(idx):
            problems.append(f"bowtie index file not found: {f}")

    # ── Always required: GTF (located by presence, any filename) ──
    _gtf_path, _gtf_err = resolve_gtf_path(refs.get("gtf"))
    if _gtf_err:
        problems.append(_gtf_err)

    # ── Genome FASTA + samtools faidx index (only for sequence_annotation) ──
    if need_genome:
        fa = refs.get("genome_fasta")
        if not fa:
            problems.append("genome_fasta: sequence_annotation is on but genome_fasta not set in config.yaml")
        elif not Path(fa).exists():
            problems.append(f"genome FASTA not found: {fa}")
        elif not Path(f"{fa}.fai").exists():
            problems.append(f"genome FASTA not indexed: {fa}  ->  run: samtools faidx {fa}")

    # ── Optional: piRNA FASTAs (only if any sample needs piRNA rescue) ──
    if need_pirna:
        pirna = refs.get("pirna_refs") or []
        if isinstance(pirna, str):
            pirna = [pirna]
        if not pirna:
            problems.append("pirna_refs: piRNA matching requested but no pirna_refs set in config.yaml")
        for p in pirna:
            if not Path(p).exists():
                problems.append(f"piRNA reference not found: {p}")

    # ── Optional: miRNA CSV (city_finder + miRNA comparison) ──
    if need_mirna:
        csv = refs.get("mirna_csv")
        if not csv:
            problems.append("mirna_csv: miRNA/city_finder analysis requested but mirna_csv not set in config.yaml")
        elif not Path(csv).exists():
            problems.append(f"miRNA reference CSV not found: {csv}")

    return problems


def modules_to_needs(config) -> dict:
    """Translate the config `modules:` block into the reference sets it requires."""
    m = config.get("modules", {}) or {}
    def on(k, default=False):
        return bool(m.get(k, default))
    return {
        "need_pirna":  on("pirna_matching", True),   # on by default (mirrors Snakefile)
        "need_genome": on("sequence_annotation"),
        "need_mirna":  on("city_finder") or on("comparisons"),
    }


def validate_references(exp_df, config, *, need_pirna=False, need_genome=False, need_mirna=False) -> None:
    """Fail fast (SystemExit) if any bio sample's organism references are missing.

    Only the reference sets required by the enabled analysis modules are checked (the
    need_* flags; typically derived via modules_to_needs(config)). bowtie_index + gtf are
    always required for bio. VCG/DNO samples need no references and are ignored, so a
    synthetic-only or lean-core run passes with a near-empty references/ directory.
    """
    if "experiment_type" not in exp_df.columns:
        return  # nothing to route; leave manifest validation to the Snakefile
    bio_rows = exp_df[exp_df["experiment_type"].astype(str).str.lower() == "bio"]
    if bio_rows.empty:
        return

    if "organism" not in bio_rows.columns or (bio_rows["organism"].astype(str).str.strip() == "").any():
        raise SystemExit(
            "\nReference check failed: one or more bio samples have no 'organism' set in "
            "experiments.tsv.\nEvery bio row needs an organism matching a block under "
            "'organisms:' in config.yaml.\n"
        )

    organisms = sorted({str(o).strip() for o in bio_rows["organism"]})
    report: list[str] = []
    for organism in organisms:
        try:
            refs = organism_refs(organism, config)
        except KeyError as e:
            report.append(f"organism '{organism}':\n  - {e}")
            continue
        problems = _check_organism(
            organism, refs,
            need_pirna=need_pirna, need_mirna=need_mirna, need_genome=need_genome,
        )
        if problems:
            lines = "\n".join(f"  - {p}" for p in problems)
            report.append(f"organism '{organism}':\n{lines}")

    if report:
        raise SystemExit(
            "\n"
            + "=" * 70
            + "\nReference check failed. The bio (genome-alignment) track needs reference\n"
            "files that are not shipped with this repo. Download them and place them where\n"
            "config.yaml points (see references/README.md for the expected layout).\n"
            + "=" * 70
            + "\n\n"
            + "\n\n".join(report)
            + "\n"
        )


# ── Standalone check ──────────────────────────────────────────────────────────

def _main() -> None:
    import argparse
    import sys

    import pandas as pd
    import yaml

    ap = argparse.ArgumentParser(description="Validate bio reference files against config.yaml.")
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--tsv", default="experiments.tsv")
    ap.add_argument("--organism", help="Check just this organism (ignores the TSV).")
    args = ap.parse_args()

    with open(args.config) as fh:
        config = yaml.safe_load(fh) or {}

    if args.organism:
        # Manual check verifies every reference set regardless of module config.
        refs = organism_refs(args.organism, config)
        problems = _check_organism(args.organism, refs,
                                   need_pirna=True, need_mirna=True, need_genome=True)
        if problems:
            print(f"organism '{args.organism}':")
            for p in problems:
                print(f"  - {p}")
            sys.exit(1)
        print(f"organism '{args.organism}': all reference files present.")
        return

    exp_df = pd.read_csv(args.tsv, sep="\t", dtype=str, keep_default_na=False)
    validate_references(exp_df, config, **modules_to_needs(config))   # raises SystemExit on failure
    print("Reference check passed for all bio samples (given enabled modules).")


if __name__ == "__main__":
    _main()
