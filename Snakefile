# Snakefile for tinyRNA-seq pipeline
#
# Use the manifest (experiments.tsv) to define inputs for the workflow. Each sample is routed by its
# `experiment_type` column:
#
#   vcg  → synthetic VCG track          (no genome alignment)
#   dno  → synthetic random-RNA (DNO)   (no genome alignment)
#   bio  → biological small-RNA track   (genome alignment via bowtie1)
#
# Shared front-end for vcg / dno / bio-from-raw-fastq:  quality_filter → parse_reads.
# Bio reference files are resolved per sample from the `organism` column against the
# organisms: block in config.yaml, and validated at startup (scripts/references.py).
#
# Usage:
#   Stage 1 (discover inputs):  python scripts/make_input_tsv.py   # → experiments.tsv
#             then review experiments.tsv (experiment_type / organism / adaptors)
#   Stage 2 (run):
#     snakemake --cores 8 -n                       # dry run
#     snakemake --cores 8
#     snakemake --cores 8 --config tsv=other.tsv

import sys
import re
import itertools
from pathlib import Path

import pandas as pd

configfile: "config.yaml"

# reference resolution + startup validation (scripts/references.py)
sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))
from references import validate_references, organism_refs, resolve_gtf_path

# ── Load manifest (Stage 1: auto-discover from raw_reads if empty/missing) ─────

import make_input_tsv as _mk

TSV     = config.get("tsv", "experiments.tsv")
RAW_DIR = config.get("raw_dir", "raw_reads")


def _tsv_has_samples(path) -> bool:
    """True if the TSV exists and has at least one non-blank data row (beyond the header)."""
    p = Path(path)
    if not p.exists():
        return False
    with open(p) as fh:
        return any(ln.strip() for i, ln in enumerate(fh) if i > 0)


# An empty (header-only) or missing manifest auto-populates from raw_reads/ and stops
# so you can review it, then re-run. Set --config tsv=... to point elsewhere.
if not _tsv_has_samples(TSV):
    _raw = Path(RAW_DIR)
    _rows = _mk.build_rows(
        raw_dir=_raw, default_type=_mk.DEF_TYPE, organism=_mk.DEF_ORGANISM,
        vcg_path=_mk.DEF_VCG, adaptor_3p=_mk.DEF_ADAPTOR_3P,
        read_layout=_mk.DEF_READ_LAYOUT, phred=_mk.DEF_PHRED, absolute=False,
    ) if _raw.is_dir() else []
    if not _rows:
        raise SystemExit(
            f"\n'{TSV}' has no samples and nothing was found in '{RAW_DIR}/'.\n"
            f"Add your reads as {RAW_DIR}/<experiment_name>/*.fastq[.gz] and re-run.\n"
        )
    _mk.write_tsv(_rows, Path(TSV))
    _missing = [r["experiment_name"] for r in _rows if not r["fastq_file"] and not r["trimmed_fastq"]]
    raise SystemExit(
        f"\nAuto-discovered {len(_rows)} experiment(s) from '{RAW_DIR}/' → wrote '{TSV}'.\n"
        + (f"  WARNING: no FASTQ found for: {', '.join(_missing)}\n" if _missing else "")
        + "REVIEW experiments.tsv (experiment_type / organism / adaptors), then re-run snakemake.\n"
    )

exp_df = pd.read_csv(TSV, sep="\t", dtype=str, keep_default_na=False)
if "experiment_name" not in exp_df.columns:
    raise SystemExit("experiments.tsv must have an 'experiment_name' column.")
exp_df = exp_df.set_index("experiment_name", drop=False)
SAMPLES = exp_df["experiment_name"].tolist()

DEFAULTS     = config.get("defaults", {})
SCHEMAS_FILE = "adaptor_schemas.yaml"
RESULTS      = "results"

# ── Optional analysis modules (lean core: default off) ─────────────────────────
MODULES = config.get("modules", {})


def _mod(name):
    return bool(MODULES.get(name, False))


DO_PIRNA        = bool(MODULES.get("pirna_matching", True))   # on by default (see config.yaml)
DO_ANNOTATION   = _mod("sequence_annotation")
DO_METRICS      = _mod("per_sample_metrics")   # by_type composition + miRNA_dispersion (data only)
DO_CITY         = bool(MODULES.get("city_finder", True))  # on by default; city_finder.tsv (needs mirna_csv)
DO_COMPARISON   = _mod("comparisons")          # cross-sample comparison TABLES (data only; plots in figures/)
DO_DS           = _mod("comparison_downsample") # rarefied equal-depth comparison supplement (opt-in)

# NOTE: this release produces tabular data only — no plotting. Figure/analysis code
# (per-sample plots, cross-sample comparisons, VCG graphs, cofold/bias analysis, Excel
# export) lives in figures/ and is run manually against an existing results/ tree; it is
# intentionally not wired into this workflow. See figures/ and the README "Figures" note.


# ── Manifest helpers ───────────────────────────────────────────────────────────

def _row(sample):
    return exp_df.loc[sample]


def cfg(sample, col):
    """Required value — raises if the column is missing or the cell is blank."""
    val = str(_row(sample).get(col, "")).strip()
    if not val:
        raise SystemExit(f"Sample '{sample}' missing required column '{col}'")
    return val


def cfg_opt(sample, col, default=None):
    """Optional value with layered fallback.

    blank cell   → config defaults[col], else `default`
    cell 'None'  → "" (explicit empty override; e.g. force no adaptor tail)
    otherwise    → the cell value
    """
    val = str(_row(sample).get(col, "")).strip()
    if val == "None":
        return ""
    if val == "":
        return DEFAULTS.get(col, default)
    return val


def _truthy(v) -> bool:
    return str(v).strip().lower() in ("yes", "true", "1", "y")


# ── Intermediate-file retention ────────────────────────────────────────────────
# Lean by default: large regenerable intermediates (uncompressed SAM, pre-sort BAMs,
# collapsed FASTA, quality-passed FASTQ) are marked temp() and deleted once their
# downstream rule has consumed them. Set keep_intermediates: true in config.yaml
# (or --config keep_intermediates=true) to retain everything. Requesting such a file
# directly as a target also keeps it, regardless of this flag.
KEEP_INTERMEDIATES = _truthy(config.get("keep_intermediates", False))


def _tmp(path):
    """Wrap an output path in temp() unless the user opted to keep intermediates."""
    return path if KEEP_INTERMEDIATES else temp(path)


def _has(sample, col) -> bool:
    return bool(str(_row(sample).get(col, "")).strip())


def has_random_regions(sample: str) -> bool:
    """True if the layout has any randomised bases (drives ligation_junction eligibility)."""
    import yaml
    with open(SCHEMAS_FILE) as fh:
        s = yaml.safe_load(fh)
    layout = s["read_layouts"][cfg(sample, "read_layout")]
    no_random_5p = {"DEFINED", "FLUSH"}
    return layout["five_prime"] not in no_random_5p or layout["three_prime"] != "NORMAL"


def bio_ref(sample, key, default=None):
    """Resolve a reference-file path for a bio sample from its organism block in config."""
    return organism_refs(cfg(sample, "organism"), config).get(key, default)


def bio_gtf(sample):
    """Resolve the GTF for a bio sample by presence, not filename (see resolve_gtf_path)."""
    path, err = resolve_gtf_path(bio_ref(sample, "gtf"))
    if err:
        raise SystemExit(f"[gtf] organism '{cfg(sample, 'organism')}': {err}")
    return path


def _etype(sample) -> str:
    return str(_row(sample).get("experiment_type", "")).strip().lower()


# ── Manifest validation (checks values against allowed sets; reports all at once) ─
# Valid read_layouts come from adaptor_schemas.yaml and valid organisms from config.yaml,
# so error messages enumerate the live options (never drift from the README).

def validate_manifest():
    import yaml
    with open(SCHEMAS_FILE) as fh:
        layouts = sorted(yaml.safe_load(fh).get("read_layouts", {}))
    organisms = sorted(config.get("organisms", {}) or {})
    entry_cols = ("fastq_file", "trimmed_fastq", "sam_path", "sorted_bam")
    _layout_list = ", ".join(layouts)
    _org_list = ", ".join(organisms) or "(none defined in config.yaml organisms:)"

    problems = []
    for s in SAMPLES:
        et = _etype(s)
        if et not in ("vcg", "dno", "bio", ""):
            problems.append(f"{s}: experiment_type '{et}' invalid — use vcg / dno / bio")
            continue
        is_synth = et in ("vcg", "dno", "")
        is_bio   = et == "bio"
        bio_raw  = is_bio and _has(s, "fastq_file")

        # front-end samples (synthetic, or bio starting from raw fastq) need layout + adaptor
        if is_synth or bio_raw:
            rl = str(_row(s).get("read_layout", "")).strip()
            if not rl:
                problems.append(f"{s}: read_layout is required — one of: {_layout_list}")
            elif rl not in layouts:
                problems.append(f"{s}: read_layout '{rl}' unknown — one of: {_layout_list}")
            if not _has(s, "adaptor_3p"):
                problems.append(f"{s}: adaptor_3p is required (the 3' sequencing adaptor)")

        if is_synth:
            if not _has(s, "fastq_file"):
                problems.append(f"{s}: fastq_file is required for {et or 'vcg'} samples")
            if et in ("vcg", "") and not _has(s, "vcg_list_path"):
                problems.append(f"{s}: vcg_list_path is required for vcg samples")

        if is_bio:
            org = str(_row(s).get("organism", "")).strip()
            if not org:
                problems.append(f"{s}: organism is required for bio samples — one of: {_org_list}")
            elif org not in organisms:
                problems.append(f"{s}: organism '{org}' not defined in config.yaml — one of: {_org_list}")
            if not any(_has(s, c) for c in entry_cols):
                problems.append(f"{s}: bio sample needs one of fastq_file / trimmed_fastq / sam_path / sorted_bam")

        ph = str(_row(s).get("phred_cutoff", "")).strip()
        if ph and not ph.isdigit():
            problems.append(f"{s}: phred_cutoff '{ph}' must be a whole number")
        cf = str(_row(s).get("cofold_analysis", "")).strip().lower()
        if cf and cf not in ("yes", "no", "true", "false", "y", "n", "1", "0"):
            problems.append(f"{s}: cofold_analysis '{cf}' should be yes or no")

    if problems:
        raise SystemExit(
            "\nexperiments.tsv has problems (fix these rows, then re-run):\n"
            + "\n".join(f"  - {p}" for p in problems) + "\n"
        )


validate_manifest()


# ── Reference preflight (bio samples only; no-op for a vcg/dno-only run) ────────

validate_references(
    exp_df, config,
    need_pirna=DO_PIRNA,
    need_genome=DO_ANNOTATION,
    need_mirna=DO_CITY or DO_COMPARISON,
)


# ── Sample routing ─────────────────────────────────────────────────────────────

_VALID_TYPES = {"vcg", "dno", "bio", ""}
_bad = {s: _etype(s) for s in SAMPLES if _etype(s) not in _VALID_TYPES}
if _bad:
    raise SystemExit(
        "Unknown experiment_type values in experiments.tsv (use vcg / dno / bio):\n"
        + "\n".join(f"  {s}: '{t}'" for s, t in _bad.items())
    )

VCG_SAMPLES = [s for s in SAMPLES if _etype(s) in ("vcg", "")]   # blank defaults to vcg
DNO_SAMPLES = [s for s in SAMPLES if _etype(s) == "dno"]
BIO_SAMPLES = [s for s in SAMPLES if _etype(s) == "bio"]

# Bio entry points (mutually exclusive per sample, checked in priority order)
BIO_FASTQ   = [s for s in BIO_SAMPLES if _has(s, "fastq_file")]
BIO_SAM     = [s for s in BIO_SAMPLES if not _has(s, "fastq_file") and _has(s, "sam_path")]
BIO_TRIMMED = [s for s in BIO_SAMPLES
               if not _has(s, "fastq_file") and not _has(s, "sam_path")
               and not _has(s, "sorted_bam") and _has(s, "trimmed_fastq")]
BIO_BAM     = [s for s in BIO_SAMPLES
               if not _has(s, "fastq_file") and not _has(s, "sam_path")
               and not _has(s, "trimmed_fastq") and _has(s, "sorted_bam")]
BIO_SAMCONV = BIO_FASTQ + BIO_SAM + BIO_TRIMMED         # need sam → sorted bam
BIO_CITY    = BIO_FASTQ + BIO_TRIMMED                   # city_finder eligible

_bio_unrouted = [s for s in BIO_SAMPLES
                 if s not in BIO_FASTQ + BIO_SAM + BIO_TRIMMED + BIO_BAM]
if _bio_unrouted:
    raise SystemExit(
        "These bio samples have no input set (need one of fastq_file / trimmed_fastq / "
        "sam_path / sorted_bam):\n  " + "\n  ".join(_bio_unrouted)
    )

# Shared front-end (quality_filter + parse_reads)
FRONTEND_SAMPLES = VCG_SAMPLES + DNO_SAMPLES + BIO_FASTQ

# Cofold / bias analysis: opt-in per sample via cofold_analysis=yes (default no), VCG-only.
# Produces summary tables.
COFOLD_SAMPLES  = [s for s in VCG_SAMPLES if _truthy(cfg_opt(s, "cofold_analysis"))]
LIGJUNC_SAMPLES = [s for s in COFOLD_SAMPLES if has_random_regions(s)]


def _pat(samples):
    return "|".join(re.escape(s) for s in samples) if samples else "NOMATCH_EMPTY"


# ── bio helpers (biotypes, small-RNA regex) ────────────────────────────────────

def _biotypes(sample):
    val = cfg_opt(sample, "plot_biotypes", "")
    if val:
        if isinstance(val, list):
            return [str(b) for b in val]
        return [b.strip() for b in str(val).split(",") if b.strip()]
    return ["miRNA", "piRNA", "snoRNA", "tRNA"]


def _small_rna_regex(sample):
    val = str(_row(sample).get("small_rna_types", "")).strip()
    if val and val != "None":
        types = [t.strip() for t in val.split(",") if t.strip()]
    else:
        default_types = DEFAULTS.get("small_rna_types", [
            "miRNA", "piRNA", "snoRNA", "snRNA", "tRNA", "Mt_tRNA",
            "rRNA", "Mt_rRNA", "Y_RNA", "vaultRNA", "scaRNA", "misc_RNA",
        ])
        types = [str(t) for t in default_types]
    return "|".join(types)


ALL_BIOTYPES = sorted({bt for s in BIO_SAMPLES for bt in _biotypes(s)})

# miRNA reference for the cross-sample comparison spearman (organism of first bio sample)
BIO_MIRNA_CSV = bio_ref(BIO_SAMPLES[0], "mirna_csv") if (BIO_SAMPLES and DO_COMPARISON) else None


# ── Target assembly ────────────────────────────────────────────────────────────

# VCG + DNO share the same terminal outputs
SYNTH_SAMPLES = VCG_SAMPLES + DNO_SAMPLES
SYNTH_TARGETS = (
    # core: normalised counts (CPM) + concatenated log
    [f"{RESULTS}/{s}/{s}.normalised.tsv" for s in SYNTH_SAMPLES]
    + [f"{RESULTS}/{s}/{s}_Norm.csv" for s in SYNTH_SAMPLES]
    + [f"{RESULTS}/{s}/{s}.pipeline.log" for s in SYNTH_SAMPLES]
    # cofold / bias summary TABLES — opt-in per sample (cofold_analysis=yes)
    + [f"{RESULTS}/{s}/{s}.first_ligation_reads.tsv" for s in COFOLD_SAMPLES]
    + [f"{RESULTS}/{s}/{s}.secondary_struct_reads.tsv" for s in COFOLD_SAMPLES]
    + [f"{RESULTS}/{s}/{s}.ligation_junction_summary.tsv" for s in LIGJUNC_SAMPLES]
)

# ── bio targets ──
_TRIALS = config.get("trials", "")
TRIAL_TARGETS = []
if _TRIALS in ("cpm_mapped", "both"):
    TRIAL_TARGETS += [f"{RESULTS}/{s}/{s}.counts.smallRNA.withCPM.mapped.tsv" for s in BIO_SAMPLES]
if _TRIALS in ("hier", "both"):
    TRIAL_TARGETS += [f"{RESULTS}/{s}/{s}.counts.smallRNA.withCPM.hier.tsv" for s in BIO_SAMPLES]

COUNTS_CPM = [f"{RESULTS}/{s}/{s}.counts.smallRNA.withCPM.tsv" for s in BIO_SAMPLES]

BIO_SAMPLE_TARGETS = []
for _s in BIO_SAMPLES:
    # core: indexed BAM, featureCounts summary, counts/CPM table, log
    BIO_SAMPLE_TARGETS += [
        f"{RESULTS}/{_s}/{_s}.NH.sorted.bam.bai",
        f"{RESULTS}/{_s}/{_s}.allRNA.gene.txt.summary",
        f"{RESULTS}/{_s}/{_s}.counts.smallRNA.withCPM.tsv",
        f"{RESULTS}/{_s}/{_s}.pipeline_bio.log",
    ]
    if DO_ANNOTATION:
        BIO_SAMPLE_TARGETS.append(f"{RESULTS}/{_s}/{_s}.counts.smallRNA.withCPM.annotated.tsv")
    if DO_METRICS:
        BIO_SAMPLE_TARGETS += [
            f"{RESULTS}/{_s}/{_s}.counts.smallRNA.by_type.tsv",
            f"{RESULTS}/{_s}/{_s}.miRNA_dispersion.txt",
        ]
    if DO_CITY and _s in BIO_CITY:
        BIO_SAMPLE_TARGETS.append(f"{RESULTS}/{_s}/{_s}.city_finder.tsv")

# Cross-sample comparison TABLES (data only; the plots are in figures/). Opt-in via
# the `comparisons` module. Needs >= 1 bio sample and the organism's mirna_csv.
MIRNA_CMP_OUT    = f"{RESULTS}/comparison/miRNA/miRNA_all_experiments.csv"
SPEARMAN_OUT     = f"{RESULTS}/comparison/miRNA/spearman/spearman_summary.csv"
SMALLRNA_CMP_OUT = f"{RESULTS}/comparison/smallRNA/smallRNA_all_experiments.csv"
RUN_SUMMARY_OUT  = f"{RESULTS}/comparison/run_summary.tsv"   # tidy long-format (data)
RUN_SUMMARY_TXT  = f"{RESULTS}/comparison/run_summary.txt"   # rendered human-readable

BIO_COMPARISON_TARGETS = []
if BIO_SAMPLES and DO_COMPARISON:
    BIO_COMPARISON_TARGETS = [
        MIRNA_CMP_OUT, SPEARMAN_OUT, SMALLRNA_CMP_OUT, RUN_SUMMARY_OUT, RUN_SUMMARY_TXT,
    ]

# Rarefied (equal-depth) comparison supplement — opt-in (module comparison_downsample);
# needs comparisons on and >= 2 bio samples. Full-depth comparison above is untouched.
RAREFIED_OUT = f"{RESULTS}/comparison/rarefied/unique_species_by_biotype.rarefied.tsv"
if BIO_COMPARISON_TARGETS and DO_DS and len(BIO_SAMPLES) >= 2:
    BIO_COMPARISON_TARGETS.append(RAREFIED_OUT)


# ── Target ─────────────────────────────────────────────────────────────────────

rule all:
    input:
        SYNTH_TARGETS + BIO_SAMPLE_TARGETS + BIO_COMPARISON_TARGETS + TRIAL_TARGETS


wildcard_constraints:
    sample=_pat(SAMPLES),
    bt=_pat(ALL_BIOTYPES),


# ══════════════════════════════════════════════════════════════════════════════
#  SHARED FRONT-END  (vcg + dno + bio-from-raw-fastq)
# ══════════════════════════════════════════════════════════════════════════════

rule quality_filter:
    input:
        fastq=lambda wc: cfg(wc.sample, "fastq_file"),
    output:
        high=_tmp(f"{RESULTS}/{{sample}}/{{sample}}.high_quality.fastq.gz"),
        low=f"{RESULTS}/{{sample}}/{{sample}}.low_quality.fastq.gz",
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.quality_filter.log",
    params:
        phred=lambda wc: int(cfg_opt(wc.sample, "phred_cutoff", 30)),
    wildcard_constraints:
        sample=_pat(FRONTEND_SAMPLES),
    threads: 1
    script:
        "scripts/quality_filter.py"


rule parse_reads:
    input:
        fastq=f"{RESULTS}/{{sample}}/{{sample}}.high_quality.fastq.gz",
        schemas=SCHEMAS_FILE,
    output:
        parsed=f"{RESULTS}/{{sample}}/{{sample}}.parsed_reads.tsv",
        dropped=f"{RESULTS}/{{sample}}/{{sample}}.bad_reads.tsv",
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.parse_reads.log",
    params:
        read_layout=lambda wc: cfg(wc.sample, "read_layout"),
        adaptor_3p=lambda wc: cfg(wc.sample, "adaptor_3p"),
    wildcard_constraints:
        sample=_pat(FRONTEND_SAMPLES),
    threads: 1
    script:
        "scripts/parse_reads.py"


# ══════════════════════════════════════════════════════════════════════════════
#  VCG TRACK  (synthetic; no genome alignment)
# ══════════════════════════════════════════════════════════════════════════════

rule vcg_filter:
    input:
        parsed=f"{RESULTS}/{{sample}}/{{sample}}.parsed_reads.tsv",
        vcg_list=lambda wc: cfg(wc.sample, "vcg_list_path"),
    output:
        filtered=f"{RESULTS}/{{sample}}/{{sample}}.vcg_filtered.tsv",
        excluded=f"{RESULTS}/{{sample}}/{{sample}}.vcg_excluded.tsv",
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.vcg_filter.log",
    wildcard_constraints:
        sample=_pat(VCG_SAMPLES),
    threads: 1
    script:
        "scripts/vcg_filter.py"


# ── Cofold analysis (opt-in: cofold_analysis=yes; summary tables) ──

rule first_ligation_cofold:
    input:
        filtered=f"{RESULTS}/{{sample}}/{{sample}}.vcg_filtered.tsv",
        schemas=SCHEMAS_FILE,
    output:
        reads=f"{RESULTS}/{{sample}}/{{sample}}.first_ligation_reads.tsv",
        summary=f"{RESULTS}/{{sample}}/{{sample}}.first_ligation_summary.tsv",
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.first_ligation_cofold.log",
    params:
        read_layout=lambda wc: cfg(wc.sample, "read_layout"),
        max_reads_per_vcg=config.get("max_reads_per_vcg", 2000),
        adaptor_3p_tail=lambda wc: cfg_opt(wc.sample, "adaptor_3p_tail"),
        cofold_temperature=lambda wc: float(cfg_opt(wc.sample, "cofold_temperature", 16.0)),
    wildcard_constraints:
        sample=_pat(COFOLD_SAMPLES),
    threads: 1
    script:
        "scripts/first_ligation_cofold.py"


rule secondary_structures:
    input:
        filtered=f"{RESULTS}/{{sample}}/{{sample}}.vcg_filtered.tsv",
        schemas=SCHEMAS_FILE,
    output:
        reads=f"{RESULTS}/{{sample}}/{{sample}}.secondary_struct_reads.tsv",
        summary=f"{RESULTS}/{{sample}}/{{sample}}.secondary_struct_summary.tsv",
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.secondary_structures.log",
    params:
        read_layout=lambda wc: cfg(wc.sample, "read_layout"),
        max_reads_per_vcg=config.get("max_reads_per_vcg", 2000),
        adaptor_5p_ext=lambda wc: cfg_opt(wc.sample, "adaptor_5p_ext"),
        adaptor_3p_tail=lambda wc: cfg_opt(wc.sample, "adaptor_3p_tail"),
        cofold_temperature=lambda wc: float(cfg_opt(wc.sample, "cofold_temperature", 16.0)),
    wildcard_constraints:
        sample=_pat(COFOLD_SAMPLES),
    threads: 1
    script:
        "scripts/secondary_structures.py"


rule ligation_junction:
    input:
        filtered=f"{RESULTS}/{{sample}}/{{sample}}.vcg_filtered.tsv",
        schemas=SCHEMAS_FILE,
    output:
        summary=f"{RESULTS}/{{sample}}/{{sample}}.ligation_junction_summary.tsv",
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.ligation_junction.log",
    params:
        read_layout=lambda wc: cfg(wc.sample, "read_layout"),
    wildcard_constraints:
        sample=_pat(LIGJUNC_SAMPLES),
    threads: 1
    script:
        "scripts/ligation_junction.py"


rule umi_collapse:
    input:
        filtered=f"{RESULTS}/{{sample}}/{{sample}}.vcg_filtered.tsv",
        schemas=SCHEMAS_FILE,
    output:
        collapsed=f"{RESULTS}/{{sample}}/{{sample}}.umi_collapsed.tsv",
        collapse_stats=f"{RESULTS}/{{sample}}/{{sample}}.collapse_stats.tsv",
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.umi_collapse.log",
    params:
        read_layout=lambda wc: cfg(wc.sample, "read_layout"),
    wildcard_constraints:
        sample=_pat(VCG_SAMPLES),
    threads: 1
    script:
        "scripts/umi_collapse.py"


rule normalise:
    input:
        collapsed=f"{RESULTS}/{{sample}}/{{sample}}.umi_collapsed.tsv",
        vcg_list=lambda wc: cfg(wc.sample, "vcg_list_path"),
    output:
        normalised=f"{RESULTS}/{{sample}}/{{sample}}.normalised.tsv",
        csv=f"{RESULTS}/{{sample}}/{{sample}}_Norm.csv",
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.normalise.log",
    wildcard_constraints:
        sample=_pat(VCG_SAMPLES),
    threads: 1
    script:
        "scripts/normalise.py"


VCG_STEP_LOG_NAMES = [
    "quality_filter.log", "parse_reads.log", "vcg_filter.log",
    "first_ligation_cofold.log", "secondary_structures.log", "ligation_junction.log",
    "umi_collapse.log", "normalise.log",
]


rule concat_logs:
    input:
        normalised=f"{RESULTS}/{{sample}}/{{sample}}.normalised.tsv",
    output:
        f"{RESULTS}/{{sample}}/{{sample}}.pipeline.log",
    wildcard_constraints:
        sample=_pat(VCG_SAMPLES),
    run:
        sample = wildcards.sample
        with Path(output[0]).open("w") as out:
            out.write(f"Pipeline log for sample: {sample}\n" + "=" * 60 + "\n\n")
            for name in VCG_STEP_LOG_NAMES:
                logfile = Path(f"{RESULTS}/{sample}/logs/{sample}.{name}")
                out.write("─" * 60 + f"\n  {sample}.{name}\n" + "─" * 60 + "\n")
                out.write(logfile.read_text() if logfile.exists() and logfile.stat().st_size
                          else "(not run for this sample)\n")
                out.write("\n")


# ══════════════════════════════════════════════════════════════════════════════
#  DNO TRACK  (synthetic random RNA; shares front-end)
# ══════════════════════════════════════════════════════════════════════════════

DNO_STEP_LOG_NAMES = [
    "quality_filter.log", "parse_reads.log", "dno_filter.log",
    "umi_collapse.log", "normalise_dno.log",
]


rule dno_filter:
    input:
        parsed=f"{RESULTS}/{{sample}}/{{sample}}.parsed_reads.tsv",
    output:
        filtered=f"{RESULTS}/{{sample}}/{{sample}}.vcg_filtered.tsv",
        excluded=f"{RESULTS}/{{sample}}/{{sample}}.vcg_excluded.tsv",
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.dno_filter.log",
    wildcard_constraints:
        sample=_pat(DNO_SAMPLES),
    threads: 1
    script:
        "scripts/dno_filter.py"


rule umi_collapse_dno:
    input:
        filtered=f"{RESULTS}/{{sample}}/{{sample}}.vcg_filtered.tsv",
        schemas=SCHEMAS_FILE,
    output:
        collapsed=f"{RESULTS}/{{sample}}/{{sample}}.umi_collapsed.tsv",
        collapse_stats=f"{RESULTS}/{{sample}}/{{sample}}.collapse_stats.tsv",
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.umi_collapse.log",
    params:
        read_layout=lambda wc: cfg(wc.sample, "read_layout"),
    wildcard_constraints:
        sample=_pat(DNO_SAMPLES),
    threads: 1
    script:
        "scripts/umi_collapse_dno.py"


rule normalise_dno:
    input:
        collapsed=f"{RESULTS}/{{sample}}/{{sample}}.umi_collapsed.tsv",
    output:
        normalised=f"{RESULTS}/{{sample}}/{{sample}}.normalised.tsv",
        csv=f"{RESULTS}/{{sample}}/{{sample}}_Norm.csv",
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.normalise_dno.log",
    wildcard_constraints:
        sample=_pat(DNO_SAMPLES),
    threads: 1
    script:
        "scripts/normalise_dno.py"


rule concat_logs_dno:
    input:
        normalised=f"{RESULTS}/{{sample}}/{{sample}}.normalised.tsv",
    output:
        f"{RESULTS}/{{sample}}/{{sample}}.pipeline.log",
    wildcard_constraints:
        sample=_pat(DNO_SAMPLES),
    run:
        sample = wildcards.sample
        with Path(output[0]).open("w") as out:
            out.write(f"Pipeline log for sample: {sample}\n" + "=" * 60 + "\n\n")
            for name in DNO_STEP_LOG_NAMES:
                logfile = Path(f"{RESULTS}/{sample}/logs/{sample}.{name}")
                out.write("─" * 60 + f"\n  {sample}.{name}\n" + "─" * 60 + "\n")
                out.write(logfile.read_text() if logfile.exists() and logfile.stat().st_size
                          else "(not run for this sample)\n")
                out.write("\n")


# ══════════════════════════════════════════════════════════════════════════════
#  BIO TRACK  (biological small RNA; genome alignment)
# ══════════════════════════════════════════════════════════════════════════════

if BIO_FASTQ:
    rule bio_collapse_to_fasta:
        input:
            parsed=f"{RESULTS}/{{sample}}/{{sample}}.parsed_reads.tsv",
            schemas=SCHEMAS_FILE,
        output:
            fasta=_tmp(f"{RESULTS}/{{sample}}/{{sample}}.collapsed.fasta"),
            stats=f"{RESULTS}/{{sample}}/{{sample}}.collapse_stats.tsv",
        params:
            read_layout=lambda wc: cfg(wc.sample, "read_layout"),
        log:
            f"{RESULTS}/{{sample}}/logs/{{sample}}.bio_collapse.log",
        wildcard_constraints:
            sample=_pat(BIO_FASTQ),
        script:
            "scripts/bio_collapse_to_fasta.py"

    rule bowtie_align:
        input:
            fasta=f"{RESULTS}/{{sample}}/{{sample}}.collapsed.fasta",
        output:
            sam=_tmp(f"{RESULTS}/{{sample}}/{{sample}}.mm39.sam"),
            short=f"{RESULTS}/{{sample}}/{{sample}}.short_reads.fasta",
        params:
            index=lambda wc: bio_ref(wc.sample, "bowtie_index"),
            extra=lambda wc: cfg_opt(wc.sample, "bowtie_extra",
                                     "-n 0 -l 5 -v 1 --best --strata -k 1 -m 1000"),
            min_len=lambda wc: int(cfg_opt(wc.sample, "min_read_length", 18)),
        threads:
            lambda wc: int(cfg_opt(wc.sample, "threads", 8))
        log:
            f"{RESULTS}/{{sample}}/logs/{{sample}}.bowtie_align.log",
        wildcard_constraints:
            sample=_pat(BIO_FASTQ),
        shell: r"""
            > {output.short}
            awk -v min_len={params.min_len} -v short="{output.short}" \
                '/^>/{{h=$0; next}}
                 length($0) >= min_len {{print h; print $0}}
                 length($0) <  min_len {{print h >> short; print $0 >> short}}' \
                {input.fasta} \
            | bowtie -f -p {threads} -x {params.index} \
              {params.extra} - -S {output.sam} > {log} 2>&1
        """


if BIO_SAM:
    rule stage_sam:
        input:
            lambda wc: cfg(wc.sample, "sam_path"),
        output:
            _tmp(f"{RESULTS}/{{sample}}/{{sample}}.mm39.sam"),
        wildcard_constraints:
            sample=_pat(BIO_SAM),
        shell:
            "cp {input} {output}"


# Pre-trimmed bio path: for supplier kits whose adapters are
# removed BEFORE the pipeline using scripts/trimming/*.sh. Set the trimmed R1 in the
# `trimmed_fastq` column (leave fastq_file blank). See docs/pretrimming.md.
# Skips quality_filter / parse_reads / bio_collapse; aligns the trimmed FASTQ directly.
if BIO_TRIMMED:
    rule bowtie_align_trimmed:
        input:
            fastq=lambda wc: cfg(wc.sample, "trimmed_fastq"),
        output:
            sam=_tmp(f"{RESULTS}/{{sample}}/{{sample}}.mm39.sam"),
            short=f"{RESULTS}/{{sample}}/{{sample}}.short_reads.fastq",
        params:
            index=lambda wc: bio_ref(wc.sample, "bowtie_index"),
            extra=lambda wc: cfg_opt(wc.sample, "bowtie_extra",
                                     "-n 0 -l 5 -v 1 --best --strata -k 1 -m 1000"),
            min_len=lambda wc: int(cfg_opt(wc.sample, "min_read_length", 18)),
        threads:
            lambda wc: int(cfg_opt(wc.sample, "threads", 8))
        log:
            f"{RESULTS}/{{sample}}/logs/{{sample}}.bowtie_align.log",
        wildcard_constraints:
            sample=_pat(BIO_TRIMMED),
        shell: r"""
            case "{input.fastq}" in
              *.gz) _reader="gunzip -c" ;;
              *)    _reader="cat" ;;
            esac
            > {output.short}
            $_reader {input.fastq} \
            | awk -v min_len={params.min_len} -v short="{output.short}" \
              'NR%4==1{{h=$0}} NR%4==2{{s=$0}} NR%4==3{{p=$0}}
               NR%4==0 && length(s) >= min_len {{print h; print s; print p; print $0}}
               NR%4==0 && length(s) <  min_len {{print h >> short; print s >> short; print p >> short; print $0 >> short}}' \
            | bowtie -q -p {threads} -x {params.index} \
              {params.extra} - -S {output.sam} > {log} 2>&1
        """


if BIO_BAM:
    rule stage_sorted_bam:
        input:
            lambda wc: cfg(wc.sample, "sorted_bam"),
        output:
            f"{RESULTS}/{{sample}}/{{sample}}.mm39.sorted.bam",
        wildcard_constraints:
            sample=_pat(BIO_BAM),
        shell:
            "cp {input} {output}"


if BIO_SAMCONV:
    rule sam_to_bam:
        input:
            f"{RESULTS}/{{sample}}/{{sample}}.mm39.sam",
        output:
            _tmp(f"{RESULTS}/{{sample}}/{{sample}}.mm39.bam"),
        threads:
            lambda wc: int(cfg_opt(wc.sample, "threads", 8))
        wildcard_constraints:
            sample=_pat(BIO_SAMCONV),
        shell:
            "samtools view -@ {threads} -b -o {output} {input}"

    rule coord_sort_bam:
        input:
            f"{RESULTS}/{{sample}}/{{sample}}.mm39.bam",
        output:
            f"{RESULTS}/{{sample}}/{{sample}}.mm39.sorted.bam",
        threads:
            lambda wc: int(cfg_opt(wc.sample, "threads", 8))
        wildcard_constraints:
            sample=_pat(BIO_SAMCONV),
        shell:
            "samtools sort -@ {threads} -o {output} {input}"


rule name_sort:
    input:
        f"{RESULTS}/{{sample}}/{{sample}}.mm39.sorted.bam",
    output:
        _tmp(f"{RESULTS}/{{sample}}/{{sample}}.name.bam"),
    threads:
        lambda wc: int(cfg_opt(wc.sample, "threads", 8))
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    shell:
        "samtools sort -@ {threads} -n -o {output} {input}"


rule add_NH_and_resort:
    input:
        f"{RESULTS}/{{sample}}/{{sample}}.name.bam",
    output:
        f"{RESULTS}/{{sample}}/{{sample}}.NH.sorted.bam",
    threads:
        lambda wc: int(cfg_opt(wc.sample, "threads", 8))
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.add_NH.log",
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    shell: r"""
        set -euo pipefail
        samtools view -h {input} \
        | awk 'BEGIN{{FS=OFS="\t"}}
            /^@/ {{print; next}}
            {{
              if($1!=qn && qn!=""){{
                for(i=1;i<=n;i++){{
                  line=buf[i]; sub(/\tNH:i:[0-9]+/,"",line);
                  print line "\tNH:i:" n
                }}; n=0
              }}
              qn=$1; buf[++n]=$0
            }}
            END{{
              if(n>0) for(i=1;i<=n;i++){{
                line=buf[i]; sub(/\tNH:i:[0-9]+/,"",line);
                print line "\tNH:i:" n
              }}
            }}' \
        | samtools view -@ {threads} -bS - \
        | samtools sort -@ {threads} -o {output} - 2> {log}
    """


rule index_nh:
    input:
        f"{RESULTS}/{{sample}}/{{sample}}.NH.sorted.bam",
    output:
        f"{RESULTS}/{{sample}}/{{sample}}.NH.sorted.bam.bai",
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    threads: 1
    shell:
        "samtools index {input}"


rule featureCounts_gene:
    input:
        bam=f"{RESULTS}/{{sample}}/{{sample}}.NH.sorted.bam",
        gtf=lambda wc: bio_gtf(wc.sample),
    output:
        txt=f"{RESULTS}/{{sample}}/{{sample}}.allRNA.gene.txt",
        summary=f"{RESULTS}/{{sample}}/{{sample}}.allRNA.gene.txt.summary",
        read_assign=f"{RESULTS}/{{sample}}/{{sample}}.NH.sorted.bam.featureCounts",
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.featureCounts.log",
    threads:
        lambda wc: int(cfg_opt(wc.sample, "threads", 8))
    params:
        sflag=lambda wc: int(cfg_opt(wc.sample, "strandedness", 0)),
        fc_flags=lambda wc: cfg_opt(wc.sample, "fc_flags", "-O --fraction"),
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    shell: r"""
        featureCounts -T {threads} \
          -a {input.gtf} -t exon -g gene_id \
          -s {params.sflag} \
          {params.fc_flags} \
          -R CORE \
          -o {output.txt} {input.bam} > {log} 2>&1
    """


rule slim_counts:
    input:
        f"{RESULTS}/{{sample}}/{{sample}}.allRNA.gene.txt",
    output:
        f"{RESULTS}/{{sample}}/{{sample}}.allRNA.gene.slim.tsv",
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    shell:
        r"""awk 'BEGIN{{FS=OFS="\t"}} NR>2 {{print $1, $NF}}' {input} > {output}"""


rule gene_info_from_gtf:
    input:
        gtf=lambda wc: bio_gtf(wc.sample),
    output:
        f"{RESULTS}/{{sample}}/{{sample}}.gene_info.mm39.tsv",
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    shell: r"""
        case "{input.gtf}" in
          *.gz) _cat="gzip -dc" ;;
          *)    _cat="cat" ;;
        esac
        $_cat {input.gtf} \
        | perl -F'\t' -lane 'next if /^#/; next unless $F[2] eq "gene" || $F[2] eq "tRNA";
            ($id)   = $F[8] =~ /gene_id "([^"]+)"/;
            ($name) = $F[8] =~ /gene_name "([^"]+)"/;
            ($type) = $F[8] =~ /gene_type "([^"]+)"/;
            print join("\t",$id||"",$name||"",$type||"") if $id;' \
        | sort -u > {output}
    """


rule normalize_ids_if_needed:
    input:
        counts=f"{RESULTS}/{{sample}}/{{sample}}.allRNA.gene.slim.tsv",
        info=f"{RESULTS}/{{sample}}/{{sample}}.gene_info.mm39.tsv",
    output:
        counts_sorted=f"{RESULTS}/{{sample}}/{{sample}}.allRNA.gene.slim.sorted.tsv",
        info_sorted=f"{RESULTS}/{{sample}}/{{sample}}.gene_info.mm39.sorted.tsv",
    params:
        strip=lambda wc: cfg_opt(wc.sample, "strip_versions", False),
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    run:
        if str(params.strip).lower() in ("true", "1", "yes"):
            shell(r"awk 'BEGIN{{FS=OFS=\"\t\"}} {{sub(/\.[0-9]+$/, \"\", $1)}}1' "
                  r"{input.counts} | sort -k1,1 > {output.counts_sorted}")
            shell(r"awk 'BEGIN{{FS=OFS=\"\t\"}} {{sub(/\.[0-9]+$/, \"\", $1)}}1' "
                  r"{input.info} | sort -k1,1 > {output.info_sorted}")
        else:
            shell(r"sort -k1,1 {input.counts} > {output.counts_sorted}")
            shell(r"sort -k1,1 {input.info} > {output.info_sorted}")


rule join_counts_meta:
    input:
        counts=f"{RESULTS}/{{sample}}/{{sample}}.allRNA.gene.slim.sorted.tsv",
        info=f"{RESULTS}/{{sample}}/{{sample}}.gene_info.mm39.sorted.tsv",
    output:
        f"{RESULTS}/{{sample}}/{{sample}}.with_meta.tsv",
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    shell:
        r"""join -t $'\t' -1 1 -2 1 {input.counts} {input.info} > {output}"""


rule assign_piRNA:
    input:
        bam=f"{RESULTS}/{{sample}}/{{sample}}.NH.sorted.bam",
        fc=f"{RESULTS}/{{sample}}/{{sample}}.NH.sorted.bam.featureCounts",
        refs=lambda wc: bio_ref(wc.sample, "pirna_refs"),
    output:
        f"{RESULTS}/{{sample}}/{{sample}}.pirna_matched.tsv",
    params:
        pirna_min_len=lambda wc: int(cfg_opt(wc.sample, "pirna_min_len", 24)),
        pirna_max_len=lambda wc: int(cfg_opt(wc.sample, "pirna_max_len", 32)),
        sample=lambda wc: wc.sample,
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.assign_piRNA.log",
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    script:
        "scripts/assign_piRNA.py"


rule merge_pirna:
    input:
        with_meta=f"{RESULTS}/{{sample}}/{{sample}}.with_meta.tsv",
        pirna=f"{RESULTS}/{{sample}}/{{sample}}.pirna_matched.tsv",
    output:
        f"{RESULTS}/{{sample}}/{{sample}}.with_meta.plus_pirna.tsv",
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    shell:
        "cat {input.with_meta} {input.pirna} > {output}"


rule filter_smallrna:
    input:
        # piRNA rescue (optional) concatenates rescued piRNAs before biotype filtering
        f"{RESULTS}/{{sample}}/{{sample}}."
        + ("with_meta.plus_pirna.tsv" if DO_PIRNA else "with_meta.tsv"),
    output:
        f"{RESULTS}/{{sample}}/{{sample}}.counts.smallRNA.tsv",
    params:
        types=lambda wc: _small_rna_regex(wc.sample),
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    shell:
        r"""awk 'BEGIN{{FS=OFS="\t"}} $4 ~ /^({params.types})$/' {input} > {output}"""


rule library_size:
    input:
        f"{RESULTS}/{{sample}}/{{sample}}.NH.sorted.bam",
    output:
        f"{RESULTS}/{{sample}}/{{sample}}.libsize.txt",
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    shell:
        "samtools view -c -F 4 {input} > {output}"


# Core CPM table (headered: gene_id count CPM gene_name gene_type). This is the
# bio output that all downstream analysis reads.
rule add_cpm:
    input:
        counts=f"{RESULTS}/{{sample}}/{{sample}}.counts.smallRNA.tsv",
        lib=f"{RESULTS}/{{sample}}/{{sample}}.libsize.txt",
    output:
        f"{RESULTS}/{{sample}}/{{sample}}.counts.smallRNA.withCPM.tsv",
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    shell: r"""
        LIB=$(cat {input.lib})
        {{ printf 'gene_id\tcount\tCPM\tgene_name\tgene_type\n'; \
           awk -v lib="$LIB" 'BEGIN{{FS=OFS="\t"}} \
               {{cpm=($2/lib)*1e6; print $1,$2,cpm,$3,$4}}' {input.counts}; \
        }} > {output}
    """


# Optional: append reference + observed read sequences to the counts table.
# Needs the genome FASTA in references/ . Enable via modules.sequence_annotation.
if DO_ANNOTATION:
    rule add_sequences:
        input:
            counts=f"{RESULTS}/{{sample}}/{{sample}}.counts.smallRNA.withCPM.tsv",
            bam=f"{RESULTS}/{{sample}}/{{sample}}.NH.sorted.bam",
            read_assign=f"{RESULTS}/{{sample}}/{{sample}}.NH.sorted.bam.featureCounts",
            gtf=lambda wc: bio_gtf(wc.sample),
            genome=lambda wc: bio_ref(wc.sample, "genome_fasta"),
            parsed=lambda wc: (f"{RESULTS}/{wc.sample}/{wc.sample}.parsed_reads.tsv"
                               if wc.sample in BIO_FASTQ else []),
        output:
            f"{RESULTS}/{{sample}}/{{sample}}.counts.smallRNA.withCPM.annotated.tsv",
        params:
            strip_versions=lambda wc: cfg_opt(wc.sample, "strip_versions", False),
        log:
            f"{RESULTS}/{{sample}}/logs/{{sample}}.add_sequences.log",
        wildcard_constraints:
            sample=_pat(BIO_SAMPLES),
        script:
            "scripts/add_sequences.py"


rule by_type:
    input:
        counts=f"{RESULTS}/{{sample}}/{{sample}}.counts.smallRNA.tsv",
        lib=f"{RESULTS}/{{sample}}/{{sample}}.libsize.txt",
    output:
        f"{RESULTS}/{{sample}}/{{sample}}.counts.smallRNA.by_type.tsv",
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    shell: r"""
        set -euo pipefail
        LIB=$(cat {input.lib})
        awk 'BEGIN{{FS=OFS="\t"}} {{sum[$4]+=$2}} END{{for (k in sum) print k, sum[k]}}' \
            {input.counts} \
        | sort -k2,2nr \
        | awk -v lib="$LIB" 'BEGIN{{OFS="\t"; print "biotype","count","% of mapped","% of smallRNA"}} \
            {{printf "%s\t%.6f\t%.3f\t0\n", $1, $2, ($2/lib)*100}}' \
        > {output}
        TOT=$(awk 'BEGIN{{FS=OFS="\t"}} NR>1 {{s+=$2}} END{{print s+0}}' {output})
        awk -v tot="$TOT" 'BEGIN{{FS=OFS="\t"}} NR==1{{print;next}} \
            {{print $1,$2,$3,($2/tot)*100}}' {output} > {output}.tmp
        mv {output}.tmp {output}
    """


rule city_finder:
    input:
        fasta=lambda wc: (f"{RESULTS}/{wc.sample}/{wc.sample}.short_reads.fasta"
                          if wc.sample in BIO_FASTQ
                          else f"{RESULTS}/{wc.sample}/{wc.sample}.short_reads.fastq"),
        mirna_csv=lambda wc: bio_ref(wc.sample, "mirna_csv"),
        stats=lambda wc: (f"{RESULTS}/{wc.sample}/{wc.sample}.collapse_stats.tsv"
                         if wc.sample in BIO_FASTQ else []),
    output:
        tsv=f"{RESULTS}/{{sample}}/{{sample}}.city_finder.tsv",
    params:
        min_len=lambda wc: int(cfg_opt(wc.sample, "city_finder_min_len", 12)),
        max_len=lambda wc: int(cfg_opt(wc.sample, "city_finder_max_len", 15)),
        mismatches=lambda wc: int(cfg_opt(wc.sample, "city_finder_mismatches", 0)),
        sample=lambda wc: wc.sample,
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.city_finder.log",
    wildcard_constraints:
        sample=_pat(BIO_CITY),
    script:
        "scripts/city_finder.py"


rule mirna_dispersion:
    input:
        counts=f"{RESULTS}/{{sample}}/{{sample}}.counts.smallRNA.withCPM.tsv",
    output:
        disp=f"{RESULTS}/{{sample}}/{{sample}}.miRNA_dispersion.txt",
    params:
        cpm_threshold=DEFAULTS.get("cpm_threshold", 1.0),
    log:
        f"{RESULTS}/{{sample}}/logs/{{sample}}.mirna_dispersion.log",
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    script:
        "scripts/mirna_dispersion.py"


# ── bio cross-sample comparison TABLES (module: comparisons; data only) ────────

if BIO_SAMPLES and DO_COMPARISON:
    rule mirna_comparison:
        input:
            counts=COUNTS_CPM,
            ref_csv=BIO_MIRNA_CSV,
        output:
            table=MIRNA_CMP_OUT,
            spearman=SPEARMAN_OUT,
        params:
            cpm_threshold=DEFAULTS.get("cpm_threshold", 1.0),
        log:
            f"{RESULTS}/comparison/miRNA/mirna_comparison.log",
        script:
            "scripts/mirna_comparison.py"

    rule smallrna_comparison:
        input:
            COUNTS_CPM,
        output:
            csv=SMALLRNA_CMP_OUT,
        params:
            cpm_threshold=DEFAULTS.get("cpm_threshold", 1.0),
        log:
            f"{RESULTS}/comparison/smallRNA/smallrna_comparison.log",
        script:
            "scripts/smallrna_comparison.py"

    rule run_summary:
        input:
            counts_cpm=COUNTS_CPM,
            spearman=SPEARMAN_OUT,
        output:
            tsv=RUN_SUMMARY_OUT,
            txt=RUN_SUMMARY_TXT,
        params:
            all_samples=BIO_SAMPLES,
            fastq_samples=BIO_FASTQ,
            trimmed_samples=BIO_TRIMMED,
            cpm_threshold=DEFAULTS.get("cpm_threshold", 1.0),
        log:
            f"{RESULTS}/comparison/run_summary.log",
        script:
            "scripts/run_summary.py"


# ── rarefied (equal-depth) bio comparison — supplement (module: comparison_downsample) ──
# Downsamples each bio sample's aligned reads to a common depth and re-counts, for a
# depth-fair copy of the bio comparison. The full-depth comparison above is untouched.
# One self-contained rule (loops per sample internally); off by default.
if DO_COMPARISON and DO_DS and len(BIO_SAMPLES) >= 2:
    # cityRNA rarefaction reuses the native city_finder.tsv (native pools) and the
    # short_reads pools (pools above target); present only when city_finder is on.
    _CITY_TSVS  = ([f"{RESULTS}/{s}/{s}.city_finder.tsv" for s in BIO_CITY]
                   if DO_CITY else [])
    _CITY_SHORT = (([f"{RESULTS}/{s}/{s}.short_reads.fasta" for s in BIO_FASTQ]
                    + [f"{RESULTS}/{s}/{s}.short_reads.fastq" for s in BIO_TRIMMED])
                   if DO_CITY else [])
    rule rarefied_comparison:
        input:
            bams  = expand(f"{RESULTS}/{{s}}/{{s}}.NH.sorted.bam",       s=BIO_SAMPLES),
            libs  = expand(f"{RESULTS}/{{s}}/{{s}}.libsize.txt",         s=BIO_SAMPLES),
            ginfo = expand(f"{RESULTS}/{{s}}/{{s}}.gene_info.mm39.tsv",  s=BIO_SAMPLES),
            cpm   = COUNTS_CPM,
            city  = _CITY_TSVS,
            short = _CITY_SHORT,
        output:
            RAREFIED_OUT,
        params:
            depth  = DEFAULTS.get("comparison_depth", "auto"),
            city   = DEFAULTS.get("city_finder_depth", "auto"),
            seed   = int(DEFAULTS.get("downsample_seed", 13)),
            keep   = "--keep-below" if DEFAULTS.get("comparison_keep_below", True) else "",
            outdir = f"{RESULTS}/comparison/rarefied",
        threads: 8
        log:
            f"{RESULTS}/comparison/rarefied/rarefied_comparison.log",
        shell:
            "python scripts/rarefied_comparison.py --depth {params.depth} "
            "--city-depth {params.city} --seed {params.seed} {params.keep} "
            "--threads {threads} --out {params.outdir} > {log} 2>&1"


# ── bio per-sample concatenated log ────────────────────────────────────────────

def _bio_log_inputs(wildcards):
    s = wildcards.sample
    logs = []
    if s in BIO_FASTQ:
        logs += [
            f"{RESULTS}/{s}/logs/{s}.quality_filter.log",
            f"{RESULTS}/{s}/logs/{s}.parse_reads.log",
            f"{RESULTS}/{s}/logs/{s}.bio_collapse.log",
            f"{RESULTS}/{s}/logs/{s}.bowtie_align.log",
        ]
    elif s in BIO_TRIMMED:
        logs += [f"{RESULTS}/{s}/logs/{s}.bowtie_align.log"]
    logs += [
        f"{RESULTS}/{s}/logs/{s}.add_NH.log",
        f"{RESULTS}/{s}/logs/{s}.featureCounts.log",
    ]
    if DO_PIRNA:
        logs.append(f"{RESULTS}/{s}/logs/{s}.assign_piRNA.log")
    if DO_METRICS:
        logs.append(f"{RESULTS}/{s}/logs/{s}.mirna_dispersion.log")
    if DO_CITY and s in BIO_CITY:
        logs.append(f"{RESULTS}/{s}/logs/{s}.city_finder.log")
    return logs


rule concat_logs_bio:
    input:
        _bio_log_inputs,
    output:
        f"{RESULTS}/{{sample}}/{{sample}}.pipeline_bio.log",
    wildcard_constraints:
        sample=_pat(BIO_SAMPLES),
    shell:
        "cat {input} > {output}"


# ── bio featureCounts trials (--config trials=cpm_mapped|hier|both) ────────────

if _TRIALS in ("cpm_mapped", "both"):
    rule mapped_libsize:
        input:
            f"{RESULTS}/{{sample}}/{{sample}}.NH.sorted.bam",
        output:
            f"{RESULTS}/{{sample}}/{{sample}}.libsize.mapped.txt",
        wildcard_constraints:
            sample=_pat(BIO_SAMPLES),
        shell:
            "samtools view -c -F 4 {input} > {output}"

    rule add_cpm_mapped:
        input:
            counts=f"{RESULTS}/{{sample}}/{{sample}}.counts.smallRNA.tsv",
            lib=f"{RESULTS}/{{sample}}/{{sample}}.libsize.mapped.txt",
        output:
            f"{RESULTS}/{{sample}}/{{sample}}.counts.smallRNA.withCPM.mapped.tsv",
        wildcard_constraints:
            sample=_pat(BIO_SAMPLES),
        shell: r"""
            LIB=$(cat {input.lib})
            awk -v lib="$LIB" 'BEGIN{{FS=OFS="\t"}} \
                {{cpm=($2/lib)*1e6; print $1,$2,cpm,$3,$4}}' \
                {input.counts} > {output}
        """


if _TRIALS in ("hier", "both"):
    rule apply_hierarchy:
        input:
            read_assign=f"{RESULTS}/{{sample}}/{{sample}}.NH.sorted.bam.featureCounts",
            gene_info=f"{RESULTS}/{{sample}}/{{sample}}.gene_info.mm39.tsv",
        output:
            f"{RESULTS}/{{sample}}/{{sample}}.allRNA.gene.hier.slim.tsv",
        log:
            f"{RESULTS}/{{sample}}/logs/{{sample}}.apply_hierarchy.log",
        wildcard_constraints:
            sample=_pat(BIO_SAMPLES),
        script:
            "scripts/apply_hierarchy.py"

    rule sort_hier:
        input:
            f"{RESULTS}/{{sample}}/{{sample}}.allRNA.gene.hier.slim.tsv",
        output:
            f"{RESULTS}/{{sample}}/{{sample}}.allRNA.gene.hier.slim.sorted.tsv",
        params:
            strip=lambda wc: cfg_opt(wc.sample, "strip_versions", False),
        wildcard_constraints:
            sample=_pat(BIO_SAMPLES),
        run:
            if str(params.strip).lower() in ("true", "1", "yes"):
                shell(r"awk 'BEGIN{{FS=OFS=\"\t\"}} {{sub(/\.[0-9]+$/, \"\", $1)}}1' "
                      r"{input} | sort -k1,1 > {output}")
            else:
                shell("sort -k1,1 {input} > {output}")

    rule join_counts_meta_hier:
        input:
            counts=f"{RESULTS}/{{sample}}/{{sample}}.allRNA.gene.hier.slim.sorted.tsv",
            info=f"{RESULTS}/{{sample}}/{{sample}}.gene_info.mm39.sorted.tsv",
        output:
            f"{RESULTS}/{{sample}}/{{sample}}.with_meta.hier.tsv",
        wildcard_constraints:
            sample=_pat(BIO_SAMPLES),
        shell:
            r"""join -t $'\t' -1 1 -2 1 {input.counts} {input.info} > {output}"""

    rule filter_smallrna_hier:
        input:
            f"{RESULTS}/{{sample}}/{{sample}}.with_meta.hier.tsv",
        output:
            f"{RESULTS}/{{sample}}/{{sample}}.counts.smallRNA.hier.tsv",
        params:
            types=lambda wc: _small_rna_regex(wc.sample),
        wildcard_constraints:
            sample=_pat(BIO_SAMPLES),
        shell:
            r"""awk 'BEGIN{{FS=OFS="\t"}} $4 ~ /^({params.types})$/' {input} > {output}"""

    rule add_cpm_hier:
        input:
            counts=f"{RESULTS}/{{sample}}/{{sample}}.counts.smallRNA.hier.tsv",
            lib=f"{RESULTS}/{{sample}}/{{sample}}.libsize.txt",
        output:
            f"{RESULTS}/{{sample}}/{{sample}}.counts.smallRNA.withCPM.hier.tsv",
        wildcard_constraints:
            sample=_pat(BIO_SAMPLES),
        shell: r"""
            LIB=$(cat {input.lib})
            awk -v lib="$LIB" 'BEGIN{{FS=OFS="\t"}} \
                {{cpm=($2/lib)*1e6; print $1,$2,cpm,$3,$4}}' \
                {input.counts} > {output}
        """
