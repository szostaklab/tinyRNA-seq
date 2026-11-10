# tinyRNA-seq 

A workflow for smallRNA sequencing data analysis using Snakemake.

This project is designed to handle three types of input fastqs from illumina sequencing simultaneously:

- **`bio`** — biological small RNA-seq where pre-trimming (scripts provided), a reference genome and annotations are required. 
- **`vcg`** — synthetic *Virtual Circular Genome* (or similar) experiments: match reads to a defined
  oligomer pool, UMI-collapse, normalise.
- **`dno`** — de novo oligo generation experiments synthetic random-RNA experiments: like `vcg` but with no reference pool.

Every experiment is one row in `experiments.tsv`; its `experiment_type` decides how it is processed.

<br>

## Contents

- [Requirements and Installation](#requirements-and-installation)
  - [Other software required](#other-software-required)
  - [Reference files (bio only)](#reference-files-bio-only)
  - [A note on Snakemake](#a-note-on-snakemake)
  - [Quick start](#quick-start)
- [Detailed Explanations](#detailed-explanations)
  - [Repository layout](#repository-layout)
  - [The manifest (`experiments.tsv`)](#the-manifest-experimentstsv)
  - [Read layouts (`adaptor_schemas.yaml`)](#read-layouts-adaptor_schemasyaml)
  - [Optional analysis modules](#optional-analysis-modules-configyaml--modules)
  - [piRNA & cityRNA detection](#pirna--cityrna-detection-sequence-based-recovery)
  - [Pre-trimming bio libraries from commercial kits](#pre-trimming-bio-libraries-from-commercial-kits)
  - [Script reference](#script-reference)
- [Outputs](#outputs-per-experiment-in-resultsexperiment)
- [Output breakdowns](#output-breakdowns)
- [License](#license)

<br>

# Requirements and Installation

**Platform:** Unix-like OS, tested on **macOS**, **Ubuntu 22.04**, and **Windows via WSL2**. Tested with Python 3.11.10 and Snakemake 9.13.7.

**Recommended: one conda command installs everything** (Python stack + the
command-line tools the bio/vcg tracks need):

```bash
conda env create -f environment.yml
conda activate tinyrna-seq
```

**Or pip** for the Python packages only:

```bash
pip install -r requirements.txt
```
or
```bash
pip install -r requirements.lock.txt
``` 

for the exact validated versions. 

## Other software required

| Tool | Needed for | Version used |
|---|---|---|
| `bowtie` (v1) | bio alignment | 1.3.1 |
| `samtools` | bio BAM handling | 1.21 |
| `featureCounts` (Subread) | bio quantification | 2.1.1 |
| `perl` | bio GTF parsing | 5.32.1 |
| `cutadapt` | adaptor trimming | 4.9 |


**Verify your install.** Each of these should print a version without error:

```bash
snakemake --version
bowtie --version
samtools --version
featureCounts -v
cutadapt --version
```


## Reference files (bio only)

Reference files for biological RNAs are not included in this repository please download from OSF.io if not using your own; see Data availability in the manuscript. Place them where `config.yaml → organisms:` points. Each reference is required only by the module that uses it (checked at startup only when that module is on):

| Reference | Required when | Default? |
|---|---|---|
| `bowtie_index` (6 `.ebwt`) | always (alignment) | — |
| `gtf` | always (counting) | — |
| `pirna_refs` (FASTAs) | `pirna_matching` | **on** |
| `mirna_csv` | `city_finder` or `comparisons` | city_finder: **on**, comparisons: **off** |
| `genome_fasta` (+ `.fai`) | `sequence_annotation` | off |

So a **default bio run needs the Bowtie1 index + GTF + piRNA FASTAs + miRNA CSV.** The GTF is located by presence, not filename — any single `*.gtf` in the configured directory is used.

The mouse merged annotation GTF (GENCODE + tRNA, ± piRNA/repeats) is provided on OSF.io as a
prebuilt file; place it in the configured `references/mouse/` directory. Instructions for
assembling a GTF for a new organism will be added in a later release.

Add another organism by copying the `mouse:` block in `config.yaml` and pointing it at
your files; set that name in the `organism` column.

<br>

## A note on Snakemake

[Snakemake](https://snakemake.readthedocs.io) is a workflow manager:
you declare *rules* (each with input files, output files and a command), and it works out
the dependency graph and the order to run them in. You ask for the final files you want and
it runs only the steps needed to produce them. Key things to know:

- The workflow lives in the **`Snakefile`**; per-run settings live in **`config.yaml`** and
  **`experiments.tsv`**. You edit those, not the `Snakefile`.
- Running `snakemake` builds every target for the samples in your manifest. Add `-n`
  (`--dry-run`) to see what *would* run without doing anything.
- It is **resumable and incremental**: if a run is interrupted, or you add a sample and
  re-run, Snakemake regenerates only what is missing or out of date — finished samples are
  left untouched. (`--rerun-incomplete` cleans up after a hard interruption.)
- You can build a single file by naming it as a target, e.g.
  `snakemake --cores 8 results/MySample/MySample.counts.smallRNA.withCPM.tsv`.

<br>

## Quick start

1. **Install tools** (see [Requirements and Installation](#requirements-and-installation)) and, for `bio` samples, place your reference files (see [Reference files](#reference-files-bio-only)).

2. **Add your reads** to `raw_reads/<experiment_name>/*.fastq[.gz]`.

3. **Build the manifest** by running the pipeline once — it auto-discovers your reads, writes `experiments.tsv`, then stops so you can review it:

   ```bash
   snakemake --cores 8
   ```

4. **Review `experiments.tsv`** — fix `experiment_type` / `organism` / `read_layout` / `adaptor_3p` as needed.

5. **Run for real:**

   ```bash
   snakemake --cores 8 -n      # dry run (recommended - check number of outputs matches number of input files)
   snakemake --cores 8         # run analysis
   ```

`--cores N` sets how many CPU cores Snakemake may use: it runs independent jobs in
parallel and hands multi-threaded steps (Bowtie1, samtools, featureCounts) up to `N`
threads each. `--cores 8` is a sensible default; raise it to your core count for more
samples in parallel. It is required — Snakemake
will not run without a core count (`--cores all` uses every available core).

The manifest is checked before every run: unknown values (e.g. a mistyped `read_layout`
or an `organism` not in `config.yaml`) are reported up front, listing the valid options.


<br>
<br>

# Detailed Explanations

## Repository layout

```
tinyRNA-seq_v1.0/
├── Snakefile                  the workflow control file
├── config.yaml                defaults, optional modules, organism reference paths
├── adaptor_schemas.yaml       read-layout / adaptor definitions
├── experiments.tsv            per-experiment manifest - edit this
├── scripts/                   pipeline steps (+ make_input_tsv.py, references.py)
│   └── trimming/              cutadapt wrappers for supplier kits
├── references/                bio reference files
│   ├── mouse/                 per-organism bio references (example)
│   │   ├── mm39.*.ebwt        Bowtie1 index (6 files)
│   │   ├── mm39.gtf           merged annotation GTF (GENCODE + tRNA + repeats)
│   │   ├── *.fa.gz            piRNA reference sets
│   │   └── *.mirna.csv        miRNA reference (city_finder / comparisons)
│   └── synthetic/             VCG_Oligos_Labelled.csv (vcg/dno reference)
├── raw_reads/<experiment>/    your input FASTQs
└── results/<experiment>/      outputs (git-ignored)
```

<br>
<br>
<br>

## The Manifest (`experiments.tsv`)

**This is the key file to edit to add your experiments to.** It's built in two stages: **auto-discovery** fills what it can from `raw_reads/`, then the user 
should review the files. Discovery only determines the experiment name, which FASTQ
(raw vs pre-trimmed) and its path, plus `experiment_type` *if* you nest reads as
`raw_reads/{vcg,dno,bio}/<experiment>/`. Everything else is a default you must check. See `experiments_example.tsv` and below for defaults and change them in `config.yaml`. 

### Columns

| Column | Applies to | Required? | Allowed values / format | Default |
|---|---|---|---|---|
| `experiment_name` | all | yes | short, filesystem-safe id (used throughout) | folder name (auto) |
| `experiment_type` | all | yes | `vcg` \| `dno` \| `bio` (blank → `vcg`) | `bio` (auto) |
| `organism` | bio | bio only | a key under `organisms:` in `config.yaml` (e.g. `mouse`) | `mouse` (auto) |
| `fastq_file` | vcg, dno, bio | one entry point | path to raw R1 FASTQ (`.gz` ok) | discovered path (auto) |
| `trimmed_fastq` | bio | one entry point | path to pre-trimmed R1 (supplier kits) | blank |
| `sam_path` | bio | one entry point | path to an existing SAM | blank |
| `sorted_bam` | bio | one entry point | path to an existing coord-sorted BAM | blank |
| `read_layout` | vcg, dno, bio-from-raw | yes (front-end) | a layout name from `adaptor_schemas.yaml` (see below) | `FOUR_N_MADAP` (auto) |
| `adaptor_3p` | vcg, dno, bio-from-raw | yes (front-end) | 3' sequencing adaptor sequence | `AGATCGGAAGAGCACACGTCT` (auto) |
| `vcg_list_path` | vcg | vcg only | path to the VCG reference | `references/synthetic/VCG_Oligos_Labelled.csv` (auto) |
| `phred_cutoff` | front-end | no | integer | `30` |
| `cofold_analysis` | vcg | no | `yes` \| `no`; opt-in cofold/bias **summary tables** | `no` |
| `cofold_temperature` | vcg (cofold) | no | °C for ViennaRNA folding; per-sample override | `16.0` |
| `adaptor_5p_ext` | vcg (cofold) | no | Leave blank (cofold only) | blank |
| `adaptor_3p_tail` | vcg (cofold) | no | Leave blank (cofold only) | blank |

*"one entry point" = each bio sample sets exactly one of `fastq_file` / `trimmed_fastq` /
`sam_path` / `sorted_bam`.*

*"(auto)" = pre-filled when the manifest is auto-generated from `raw_reads/`. Blank optional
columns fall back to `config.yaml → defaults:` at run time.*

### Minimum columns per track

| Track | Needs |
|---|---|
| `vcg` | `experiment_type=vcg`, `fastq_file`, `read_layout`, `adaptor_3p`, `vcg_list_path` |
| `dno` | `experiment_type=dno`, `fastq_file`, `read_layout`, `adaptor_3p` |
| `bio` (raw FASTQ) | `experiment_type=bio`, `organism`, `fastq_file`, `read_layout`, `adaptor_3p` |
| `bio` (pre-trimmed kit) | `experiment_type=bio`, `organism`, `trimmed_fastq` |
| `bio` (from SAM/BAM) | `experiment_type=bio`, `organism`, `sam_path` **or** `sorted_bam` |

## Read layouts (`adaptor_schemas.yaml`)

`read_layout` must be one of these (5' structure + 3' structure) or define your own, we have provided many 
commonly used combinations in curent smallRNA-seq protocols:

| Layout | Description | Read structure (5' → 3') |
|---|---|---|
| `FLUSH_NORMAL` | Flush 5' adaptor + defined 3' adaptor | `[VCG] \| 3' adaptor` |
| `DEFINED_NORMAL` | TCTA fixed 5' adaptor + defined 3' adaptor | `[TCTA][VCG] \| 3' adaptor` |
| `DEFINED_FOUR_N` | TCTA fixed 5' adaptor + 4N region before 3' adaptor | `[TCTA][VCG][N4] \| 3' adaptor` |
| `DEFINED_MADAP` | TCTA fixed 5' adaptor + MADAP 3' structure | `[TCTA][VCG][TRR][BR][MAC][IRR][CT] \| 3' adaptor` |
| `FOUR_N_NORMAL` | 4N 5' adaptor + defined 3' adaptor | `[5_Adap: 4bp][VCG] \| 3' adaptor` |
| `FOUR_N_FOUR_N` | 4N 5' adaptor + 4N region before 3' adaptor | `[5_Adap: 4bp][VCG][N4] \| 3' adaptor` |
| `FOUR_N_MADAP` | 4N 5' adaptor + full MADAP 3' structure (classic Ext layout) | `[5_Adap: 4bp][VCG][TRR][BR][MAC][IRR][CT] \| 3' adaptor` |
| `SUPPLIER_T_POLYA` | Supplier T polyA kit; reads pre-trimmed externally with `cutadapt -u 3 -a AAAAAAAAAA -m 2`, then point `fastq_file` at the trimmed reads; no 3' adaptor search, no UMI | `[VCG insert] (polyA and 3N removed upstream)` |

| 5' prefix | Meaning | 3' suffix | Meaning |
|---|---|---|---|
| `FLUSH` | read starts at the insert (no 5' tag) | `NORMAL` | plain 3' adaptor |
| `DEFINED` | fixed 5' sequence (e.g. TCTA) | `FOUR_N` | 4N random tag before adaptor |
| `FOUR_N` | 4N random UMI tag | `MADAP` | full MADAP structure |

Layouts with random regions provide UMIs for deduplication; `FLUSH`/`DEFINED` layouts have none and pass reads through the collapse step unchanged.

## Optional analysis modules (`config.yaml → modules:`)

A default run produces the counts/CPM table and little else. Turn modules on to add
more table outputs. Some need extra reference files (checked at startup only when on).

| Module | Default | Adds (tables) | Needs |
|---|---|---|---|
| `pirna_matching` | **on** | piRNA reads assigned by sequence (GENCODE has no piRNA loci) | `pirna_refs` |
| `sequence_annotation` | off | reference/read sequences appended to the counts table | `genome_fasta` (+ `.fai`) |
| `per_sample_metrics` | off | `by_type` biotype-composition table + `miRNA_dispersion` metrics table | — |
| `city_finder` | **on** | ultra-short reads matched within known miRNAs (`city_finder.tsv`) | `mirna_csv` |
| `comparisons` | **on** | cross-sample **tables**: miRNA/smallRNA CPM matrices, spearman-vs-reference, run summary | `mirna_csv` |

### Intermediate files & disk usage

By default the pipeline removes a handful of large intermediate files that are marked `temp()` and deleted the moment the next rule has consumed them. Everything you
would actually inspect is kept. This keeps a run's on-disk footprint small without losing
any output file.

**Keeping everything.** To retain every intermediate on disk:

```yaml
# config.yaml
keep_intermediates: true
```
```bash
# or per-run, without editing config
snakemake --cores 8 --config keep_intermediates=true
```
**Requesting a temp file directly keeps it** — e.g. `snakemake … results/S/S.mm39.sam`
builds and retains that SAM regardless of the flag.

**Auto-deleted (temp) by default**:

| File | Produced by | Why it's safe to drop |
|---|---|---|
| `<s>.high_quality.fastq.gz` | quality_filter | re-derived from the input FASTQ |
| `<s>.collapsed.fasta` | bio_collapse_to_fasta | re-derived from `parsed_reads.tsv` |
| `<s>.mm39.sam` | bowtie_align / stage_sam | uncompressed alignment (largest file; the sorted BAM supersedes it) |
| `<s>.mm39.bam` | sam_to_bam | pre-sort BAM, redundant once sorted |
| `<s>.name.bam` | name_sort | transient name-sorted copy for NH tagging |

**Always kept** (outputs and diagnostics):

| File | What it is |
|---|---|
| `<s>.NH.sorted.bam` (`.bai`) | coordinate-sorted, indexed alignment (featureCounts input; browser-ready) |
| counts tables (`<s>.counts.smallRNA.withCPM.tsv`, module tables) | core outputs |
| `<s>.bad_reads.tsv` | reads rejected at parse, with reason |
| `<s>.low_quality.fastq.gz` | reads failing the quality filter |
| `<s>.collapse_stats.tsv` | UMI dedup stats (in / unique / collapse ratio) |
| `<s>.short_reads.*` | reads too short to align, set aside (cityRNA input) |
| `<s>.allRNA.gene.txt.summary` | featureCounts assignment summary |


## piRNA & cityRNA detection (sequence-based recovery)

Two default-on bio modules recover small-RNA signal that the genome-coordinate counting
step *cannot* see: **piRNA** (GENCODE annotates no piRNA loci) and **cityRNA** (reads too
short to align at all). Both work by **exact sequence matching against a reference set**,
not by genomic coordinates.

### piRNA — `assign_piRNA.py`  (module `pirna_matching`, on by default)

**How it operates**

1. featureCounts is run with `-R CORE`, which writes a per-read assignment file
   (`<sample>.NH.sorted.bam.featureCounts`: `read_id · status · n_targets · targets`).
2. `assign_piRNA` takes only the reads whose status is **`Unassigned_NoFeatures`** — reads
   that aligned to the genome but overlapped no annotated feature — and pulls their
   sequences from the BAM.
3. Each such sequence, if within the piRNA length window
   (`pirna_min_len`–`pirna_max_len`), is compared for **exact, full-length equality** to
   every sequence in the piRNA reference FASTAs (`pirna_refs`, e.g. the piRNAdb + "gold"
   sets). A hit increments that piRNA's integer read count.
4. Output `<sample>.pirna_matched.tsv` has one row per matched piRNA
   (`gene_id · count · name · biotype`); it is appended to the counts as
   `<sample>.with_meta.plus_pirna.tsv`. The run summary reports **unique piRNAs** (distinct
   species = row count) and **reads matched** (sum of counts).


### cityRNA — `city_finder.py`  (module `city_finder`, on by default)

**How it operates**

1. Reads shorter than `min_read_length` (18 nt) never reach bowtie — the aligner rule sets
   them aside as `<sample>.short_reads.{fasta,fastq}`. Many are genuine miRNA degradation
   fragments, so this fraction is not discarded.
2. `city_finder` scans the short reads in the window
   `city_finder_min_len`–`city_finder_max_len` (**12–15 nt**).
3. The reference is `mirna_csv`; the **canonical isoform per miRNA name is the highest
   `expression` row** (sequences taken in DNA form, U→T).
4. For each *unique* short read, a window of the read's length is slid across each mature
   miRNA; the best placement within `city_finder_mismatches` (default **0**) is kept (ties
   broken by earliest position). i.e. the read must be an **exact substring of a mature
   miRNA**.
5. Output `<sample>.city_finder.tsv` has **one row per matching (read, miRNA) pair**
   (`sequence · length · read_count · mirna_name · mirna_seq · match_pos · window · mismatches`).


<br>


## Pre-trimming bio libraries from commercial kits

Libraries prepared from biological RNA using commercial kits should be trimmed before the pipeline with the wrappers in
`scripts/trimming/`; the trimmed R1 goes in the `trimmed_fastq` column and skips the
front-end. 

<br>

## Script reference

Every script in `scripts/`, listed **in execution order** and
numbered by pipeline stage. A run flows top-to-bottom; at stage 4 each sample follows
exactly one of the three tracks (chosen by `experiment_type`), then bio samples continue
to stages 5–6. Tags: *(module: X)* runs only when that module is enabled. Stage 0 is **not**
part of the Snakemake DAG. Plotting is a separate downstream layer, not part of this release.

```
0. pre-trim (kits, optional)
1. discover + validate manifest
        │
        ▼
2. quality filter ─▶ 3. parse reads
        │
        ▼
4. route by experiment_type
   ├─▶ vcg ─▶ filter ─▶ UMI collapse ─▶ normalise
   ├─▶ dno ─▶ filter ─▶ UMI collapse ─▶ normalise
   └─▶ bio ─▶ align ─▶ featureCounts ─▶ 5. per-sample tables ─▶ 6. compare
```

### 0. Supplier-kit pre-processing, before the pipeline *(kit samples only)*

| Script | What it does |
|---|---|
| `cutadapt_supplierT.sh` · `cutadapt_supplierR.sh` · `cutadapt_supplierN.sh` | cutadapt wrappers that trim supplier kit adaptors into a `trimmed_fastq`. |

### 1. Manifest discovery & startup validation

| Script | What it does |
|---|---|
| `make_input_tsv.py` | Auto-discovers input FASTQs under `raw_reads/` and writes/updates `experiments.tsv` (the first thing a run does). |
| `references.py` | Resolves bio reference files and validates them at startup (`validate_references`, name-agnostic GTF resolver). |

### 2–3. Shared front-end *(vcg / dno / bio-from-raw)*

| Script | What it does |
|---|---|
| `quality_filter.py` | **(2)** Filters reads by mean Phred score; writes high-/low-quality gzipped FASTQs. |
| `parse_reads.py` | **(3)** Trims adaptors and slices each read into structural regions per the `read_layout` in `adaptor_schemas.yaml`, applying QC filters. Writes `parsed_reads.tsv` + `bad_reads.tsv`. |

### 4. Track processing — one track per sample, routed by `experiment_type`

**4A · VCG**

| Script | What it does |
|---|---|
| `vcg_filter.py` | Keeps parsed reads whose insert matches a sequence in the VCG reference (`references/synthetic/VCG_Oligos_Labelled.csv`); drops empties and non-matches. |
| `first_ligation_cofold.py` | *(cofold)* Per-read RNA cofold of the first ligation junction → `first_ligation_{reads,summary}.tsv`. Data only. |
| `secondary_structures.py` | *(cofold)* Per-read cofold of the 5′-adaptor × insert+3′ junction → `secondary_struct_{reads,summary}.tsv`. Data only. |
| `ligation_junction.py` | *(cofold, random-region layouts)* Per-VCG junction base-bias → `ligation_junction_summary.tsv`. Data only. |
| `umi_collapse.py` | Directional UMI clustering (UMI-tools method) to remove PCR duplicates within each VCG species. |
| `normalise.py` | Counts unique reads per VCG and normalises to the reference pool (mean = 1). |

**4B · DNO**

| Script | What it does |
|---|---|
| `dno_filter.py` | Drops zero-length inserts; no reference matching (all non-empty inserts pass). |
| `umi_collapse_dno.py` | DNO variant of `umi_collapse.py`. |
| `normalise_dno.py` | CPM normalisation for DNO experiments. |

**4C · Bio** *(genome alignment → featureCounts → CPM; Bowtie1/samtools/featureCounts run as shell rules between these scripts)*

| Script | What it does |
|---|---|
| `bio_collapse_to_fasta.py` | UMI-collapses parsed bio reads and writes collapsed inserts as FASTA for Bowtie1. |
| `assign_piRNA.py` | *(module: pirna_matching)* Assigns piRNA reads by sequence match to `pirna_refs` (GENCODE has no piRNA loci). |
| `apply_hierarchy.py` | *(`--config trials=both`)* Hierarchical resolution of multi-overlap featureCounts assignments by biotype tier. |
| `add_sequences.py` | *(module: sequence_annotation)* Appends reference & observed read sequences to the counts table. |

### 5. Bio per-sample tables *(optional modules)*

| Script | What it does |
|---|---|
| `city_finder.py` | *(module: city_finder)* Finds ultra-short (12–15 nt) reads that align within known miRNAs; writes `city_finder.tsv`. |
| `mirna_dispersion.py` | *(module: per_sample_metrics)* miRNA dispersion metrics (Gini / HHI / CV / top-k) → `miRNA_dispersion.txt`. |

*(The `per_sample_metrics` module also produces `by_type.tsv` via a shell rule, not a script.)*

### 6. Cross-sample comparison tables *(module: comparisons; runs once after all samples)*

| Script | What it does |
|---|---|
| `mirna_comparison.py` | Pools per-sample miRNA CPM → `miRNA_all_experiments.csv/.tsv`; Spearman vs external reference → `spearman_summary.csv`; per-miRNA rank deltas → `rank_delta_full.csv`. Data only (no plots). |
| `smallrna_comparison.py` | Pools all-biotype CPM → `smallRNA_all_experiments.csv/.tsv`. Data only. |
| `run_summary.py` | Per-sample & cross-sample QC report (`run_summary.txt`): read funnel, biotype breakdown, quality metrics. |

<br>

<br>
<br>
<br>

# Outputs (per experiment, in `results/<experiment>/`)

| Track | Core output |
|---|---|
| `vcg` / `dno` | `<s>.normalised.tsv`, `<s>_Norm.csv` (per-VCG counts; `vcg` normalises to the reference pool, `dno` is CPM) |
| `bio` | `<s>.counts.smallRNA.withCPM.tsv` (`gene_id, count, CPM, gene_name, gene_type`), `<s>.NH.sorted.bam(.bai)`, featureCounts summary |
| all | `<s>.pipeline*.log` (concatenated step logs) |

Enabled modules add files alongside these (annotated counts, comparison tables in
`results/comparison/`).


<br>
<br>

# Output breakdowns

### Full `vcg` output breakdown

Every file a `vcg` sample produces, in pipeline order. All tables; the cofold **summary
tables** are opt-in (`cofold_analysis=yes`).

```
results/
  <sample>/
    <sample>.high_quality.fastq.gz          # quality filter — passing reads
    <sample>.low_quality.fastq.gz           # quality filter — failing reads
    <sample>.parsed_reads.tsv               # adaptor parse — one row per accepted read
    <sample>.bad_reads.tsv                   # reads rejected at parse, with reason
    <sample>.vcg_filtered.tsv               # reads matched to the oligomer pool
    <sample>.vcg_excluded.tsv               # reads that did not match the pool
    <sample>.first_ligation_reads.tsv       # (cofold) per-read first-ligation cofold
    <sample>.first_ligation_summary.tsv     # (cofold) per-VCG first-ligation summary
    <sample>.secondary_struct_reads.tsv     # (cofold) per-read secondary-structure cofold
    <sample>.secondary_struct_summary.tsv   # (cofold) per-VCG secondary-structure summary
    <sample>.ligation_junction_summary.tsv  # (cofold, random-region layouts) junction base-bias
    <sample>.umi_collapsed.tsv              # UMI-collapsed unique molecules
    <sample>.collapse_stats.tsv             # dedup stats (in / unique / collapse ratio)
    <sample>.normalised.tsv                # ★ CORE: per-oligo counts, Norm_Reads (pool mean = 1)
    <sample>_Norm.csv                      # ★ CORE: same, CSV
    <sample>.pipeline.log                  # all step logs concatenated in order
    logs/
      <sample>.quality_filter.log
      <sample>.parse_reads.log
      <sample>.vcg_filter.log
      <sample>.first_ligation_cofold.log    # (cofold)
      <sample>.secondary_structures.log     # (cofold)
      <sample>.ligation_junction.log        # (cofold)
      <sample>.umi_collapse.log
      <sample>.normalise.log
```

Cofold outputs only appear when `cofold_analysis=yes` (`ligation_junction` also needs a
random-region layout). `FLUSH` / `DEFINED` layouts carry no UMIs, so `umi_collapsed.tsv`
passes reads through unchanged and `collapse_stats.tsv` reports a 1:1 collapse.

### Full `dno` output breakdown

`dno` shares the front-end with `vcg` but has no reference pool. Every file it produces
(all core):

```
results/
  <sample>/
    <sample>.high_quality.fastq.gz          # quality filter — passing reads
    <sample>.low_quality.fastq.gz           # quality filter — failing reads
    <sample>.parsed_reads.tsv               # adaptor parse — one row per accepted read
    <sample>.bad_reads.tsv                   # reads rejected at parse, with reason
    <sample>.vcg_filtered.tsv               # reads kept (length/quality only — dno_filter reuses this name)
    <sample>.vcg_excluded.tsv               # reads excluded
    <sample>.umi_collapsed.tsv              # UMI-collapsed unique molecules
    <sample>.collapse_stats.tsv             # dedup stats (in / unique / collapse ratio)
    <sample>.normalised.tsv                # ★ CORE: per-sequence counts, CPM
    <sample>_Norm.csv                      # ★ CORE: same, CSV
    <sample>.pipeline.log                  # all step logs concatenated in order
    logs/
      <sample>.quality_filter.log
      <sample>.parse_reads.log
      <sample>.dno_filter.log
      <sample>.umi_collapse.log
      <sample>.normalise_dno.log
```

Note: `dno` reuses the `vcg_filtered.tsv` / `vcg_excluded.tsv` filenames even though there
is no pool matching (the `dno_filter` step writes them).

### Full `bio` output breakdown

Every file a `bio` sample can produce, in pipeline order. Files marked *(module)* only
appear when that optional module is enabled; everything else is core.

```
results/
  <sample>/
    <sample>.high_quality.fastq.gz          # quality filter — passing reads (raw-FASTQ entry only)
    <sample>.low_quality.fastq.gz           # quality filter — failing reads
    <sample>.parsed_reads.tsv               # adaptor parse — one row per accepted read
    <sample>.bad_reads.tsv                   # reads rejected at parse, with reason
    <sample>.collapsed.fasta                # UMI-collapsed unique molecules → Bowtie input
    <sample>.collapse_stats.tsv             # dedup stats (in / unique / collapse ratio)
    <sample>.short_reads.fasta              # reads too short to align, set aside
    <sample>.mm39.sam                        # Bowtie1 alignment (unsorted SAM)
    <sample>.mm39.bam                        # SAM → BAM
    <sample>.mm39.sorted.bam                # coordinate-sorted BAM
    <sample>.name.bam                        # name-sorted (for NH tagging)
    <sample>.NH.sorted.bam                  # NH tag (multimapper count) added → featureCounts input
    <sample>.NH.sorted.bam.bai             # BAM index
    <sample>.allRNA.gene.txt                # featureCounts raw output (all biotypes)
    <sample>.allRNA.gene.txt.summary       # featureCounts assignment summary
    <sample>.NH.sorted.bam.featureCounts   # per-read gene assignments (-R output)
    <sample>.allRNA.gene.slim.tsv          # slimmed to gene_id + count
    <sample>.allRNA.gene.slim.sorted.tsv   # sorted, for join
    <sample>.gene_info.mm39.tsv            # gene_id, gene_name, gene_type (from GTF)
    <sample>.gene_info.mm39.sorted.tsv     # sorted, for join
    <sample>.pirna_matched.tsv             # (pirna_matching, on by default) reads assigned to piRNA by sequence
    <sample>.with_meta.tsv                 # counts joined to gene metadata
    <sample>.with_meta.plus_pirna.tsv      # (pirna_matching) with_meta + sequence-matched piRNA rows
    <sample>.counts.smallRNA.tsv           # filtered to small-RNA biotypes
    <sample>.libsize.txt                    # library size (denominator for CPM)
    <sample>.counts.smallRNA.withCPM.tsv   # ★ CORE: gene_id, count, CPM, gene_name, gene_type (has header)
    <sample>.counts.smallRNA.withCPM.annotated.tsv  # (sequence_annotation) + ref_seq / top BAM+VCG sequence cols
    <sample>.counts.smallRNA.by_type.tsv   # (per_sample_metrics) CPM summed per biotype
    <sample>.miRNA_dispersion.txt          # (per_sample_metrics) Gini / HHI / CV / top-k shares
    <sample>.city_finder.tsv               # (city_finder) short-read "city" miRNA table
    <sample>.pipeline_bio.log              # all step logs concatenated in order
    logs/
      <sample>.<step>.log                  # one per rule (quality_filter, parse_reads, bowtie_align, …)

  comparison/                               # (comparisons module) cross-sample TABLES, once per run
    miRNA/miRNA_all_experiments.csv/.tsv    # pooled miRNA CPM matrix
    miRNA/spearman/spearman_summary.csv     # Spearman rho vs external reference, per sample
    miRNA/spearman/rank_delta_full.csv      # per-miRNA rank delta vs reference
    smallRNA/smallRNA_all_experiments.csv/.tsv  # pooled all-biotype CPM matrix
    run_summary.txt                          # read funnel + biotype breakdown + QC metrics
```

The `comparison/` tables are produced by the `comparisons` module.

With `--config trials=both`, featureCounts is run in extra modes, adding
`<sample>.counts.smallRNA.withCPM.mapped.tsv` (CPM over mapped reads) and
`<sample>.counts.smallRNA.withCPM.hier.tsv` (hierarchical assignment), plus their
`.libsize.mapped.txt` / `.allRNA.gene.hier.*` intermediates.

Handy invocations:
```bash
snakemake --cores 8 -n                                  # dry run
snakemake --cores 8 results/MySample/MySample.counts.smallRNA.withCPM.tsv   # one target
snakemake --cores 8 --config trials=both                # bio featureCounts trials
snakemake --cores 8 --rerun-incomplete                  # resume after interruption
```

<br>


## License

Released under the **MIT License**; see the [`LICENSE`](LICENSE) file for the full text.


<br>
