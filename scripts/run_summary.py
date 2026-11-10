#!/usr/bin/env python3
"""
run_summary.py

Per-sample and cross-sample QC report for the tinyRNA-seq pipeline.

Sections:
  1. Read processing funnel (input → filter → parse → collapse → align → assign)
  2. RNA subtype breakdown (unique species + total CPM per type, per sample)
  3. Quality metrics table (Spearman ρ, Gini, HHI, CV, Top10%, Top50%)
  4. piRNA and cityRNA summary

Snakemake script: entry point.
"""

import re
import sys
from datetime import date
from pathlib import Path

import pandas as pd

if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    all_samples     = list(snakemake.params.all_samples)
    fastq_samples   = set(snakemake.params.fastq_samples)
    trimmed_samples = set(snakemake.params.trimmed_samples)
    cpm_threshold   = float(getattr(snakemake.params, "cpm_threshold", 1.0))
    spearman_csv    = Path(snakemake.input.spearman)
    tsv_path        = Path(snakemake.output.tsv)
    txt_path        = Path(snakemake.output.txt)
else:
    raise RuntimeError("run_summary.py must be run via Snakemake script: directive")

RESULTS = Path("results")


# ── Log parsers ───────────────────────────────────────────────────────────────

def _int(s):
    """Parse comma-formatted integer string."""
    return int(s.replace(",", "").replace(" ", ""))


def parse_quality_filter(sample):
    p = RESULTS / sample / "logs" / f"{sample}.quality_filter.log"
    out = {}
    if not p.exists():
        return out
    text = p.read_text()
    m = re.search(r"Total reads:\s+([\d,]+)", text)
    if m: out["input_reads"] = _int(m.group(1))
    m = re.search(r"High quality:\s+([\d,]+)", text)
    if m: out["hq_reads"] = _int(m.group(1))
    return out


def parse_parse_reads(sample):
    p = RESULTS / sample / "logs" / f"{sample}.parse_reads.log"
    out = {}
    if not p.exists():
        return out
    text = p.read_text()
    m = re.search(r"Parsed \(pass\):\s+([\d,]+)", text)
    if m: out["parsed_reads"] = _int(m.group(1))
    return out


def parse_bio_collapse(sample):
    p = RESULTS / sample / "logs" / f"{sample}.bio_collapse.log"
    out = {}
    if not p.exists():
        return out
    text = p.read_text()
    m = re.search(r"Reads after:\s+([\d,]+)", text)
    if m: out["collapsed_reads"] = _int(m.group(1))
    return out


def parse_bowtie(sample):
    p = RESULTS / sample / "logs" / f"{sample}.bowtie_align.log"
    out = {}
    if not p.exists():
        return out
    text = p.read_text()
    m = re.search(r"reads processed:\s+([\d,]+)", text)
    if m: out["bowtie_input"] = _int(m.group(1))
    m = re.search(r"reads with at least one alignment:\s+([\d,]+)", text)
    if m: out["aligned"] = _int(m.group(1))
    m = re.search(r"reads with alignments suppressed due to -m:\s+([\d,]+)", text)
    if m: out["suppressed"] = _int(m.group(1))
    return out


def parse_libsize(sample):
    p = RESULTS / sample / f"{sample}.libsize.txt"
    if p.exists():
        return {"libsize": int(p.read_text().strip())}
    return {}


def parse_featurecounts(sample):
    p = RESULTS / sample / f"{sample}.allRNA.gene.txt.summary"
    out = {}
    if not p.exists():
        return out
    text = p.read_text()
    m = re.search(r"^Assigned\s+(\d+)", text, re.MULTILINE)
    if m: out["fc_assigned"] = int(m.group(1))
    m = re.search(r"^Unassigned_NoFeatures\s+(\d+)", text, re.MULTILINE)
    if m: out["fc_unassigned_nofeat"] = int(m.group(1))
    return out


def parse_mirna_dispersion(sample):
    p = RESULTS / sample / "logs" / f"{sample}.mirna_dispersion.log"
    out = {}
    if not p.exists():
        return out
    text = p.read_text()
    m = re.search(
        r"Gini=([\d.]+)\s+HHI=([\d.]+)\s+CV=([\d.]+)\s+Top10=([\d.]+)%\s+Top50=([\d.]+)%\s+n_miRNA=(\d+)",
        text,
    )
    if m:
        out["gini"]     = float(m.group(1))
        out["hhi"]      = float(m.group(2))
        out["cv"]       = float(m.group(3))
        out["top10pct"] = float(m.group(4))
        out["top50pct"] = float(m.group(5))
        out["n_mirna"]  = int(m.group(6))
    return out


def parse_assign_pirna(sample):
    p = RESULTS / sample / "logs" / f"{sample}.assign_piRNA.log"
    out = {}
    if not p.exists():
        return out
    text = p.read_text()
    m = re.search(r"(\d+) reads matched across (\d+) unique piRNAs", text)
    if m:
        out["pirna_reads"]  = int(m.group(1))
        out["pirna_unique"] = int(m.group(2))
    return out


def parse_city_finder(sample):
    p = RESULTS / sample / "logs" / f"{sample}.city_finder.log"
    out = {}
    if not p.exists():
        return out
    text = p.read_text()
    m = re.search(r"Reads scanned:\s+([\d,]+)", text)
    if m: out["city_scanned"] = _int(m.group(1))
    m = re.search(r"Unique seqs in range:\s+([\d,]+)", text)
    if m: out["city_unique_in_range"] = _int(m.group(1))
    m = re.search(r"Matches found:\s+([\d,]+)", text)
    if m: out["city_matches"] = _int(m.group(1))
    # 5'-end seed-restricted counts, read from the per-match TSV:
    # keep matches where the read spans the mature-miRNA seed (nt 2-8), i.e. it
    # aligns at the 5' end (mirna_match_pos <= 1) and covers position 7 (0-based).
    tsv = RESULTS / sample / f"{sample}.city_finder.tsv"
    if tsv.exists():
        try:
            ct = pd.read_csv(tsv, sep="\t")
            if {"mirna_match_pos", "length", "mirna_name"}.issubset(ct.columns):
                pos = pd.to_numeric(ct["mirna_match_pos"], errors="coerce")
                ln  = pd.to_numeric(ct["length"], errors="coerce")
                seed = (pos <= 1) & (pos + ln - 1 >= 7)
                out["city_seed_matches"] = int(seed.sum())
                out["city_seed_mirnas"]  = int(ct.loc[seed, "mirna_name"].nunique())
        except Exception:
            pass
    return out


def load_cpm_table(sample, threshold):
    p = RESULTS / sample / f"{sample}.counts.smallRNA.withCPM.tsv"
    if not p.exists():
        return pd.DataFrame()
    df = pd.read_csv(p, sep="\t")
    if "CPM" not in df.columns or "gene_type" not in df.columns:
        return pd.DataFrame()
    df["CPM"] = pd.to_numeric(df["CPM"], errors="coerce")
    return df[df["CPM"] >= threshold]


# ── Collect per-sample data ───────────────────────────────────────────────────

data = {}
for s in all_samples:
    d = {}
    if s in fastq_samples:
        d.update(parse_quality_filter(s))
        d.update(parse_parse_reads(s))
        d.update(parse_bio_collapse(s))
    d.update(parse_bowtie(s))
    d.update(parse_libsize(s))
    d.update(parse_featurecounts(s))
    d.update(parse_mirna_dispersion(s))
    d.update(parse_assign_pirna(s))
    d.update(parse_city_finder(s))
    d["cpm_df"] = load_cpm_table(s, cpm_threshold)
    data[s] = d

spearman_df = pd.read_csv(spearman_csv).set_index("sample") if spearman_csv.exists() else pd.DataFrame()


# ── Build tidy long-format table ──────────────────────────────────────────────
# One row per (section, metric, rna_type, sample, value). `rna_type` is blank for
# everything except the RNA-subtype breakdown. `value` is numeric (blank = N/A) so
# the whole file loads with pd.read_csv(sep="\t") and pivots cleanly, e.g.:
#     df = pd.read_csv("run_summary.tsv", sep="\t")
#     df[df.section=="funnel"].pivot(index="metric", columns="sample", values="value")

records = []   # dicts: section, metric, rna_type, sample, value

def add(section, metric, sample, value, rna_type=""):
    records.append({
        "section":  section,
        "metric":   metric,
        "rna_type": rna_type,
        "sample":   sample,
        "value":    value,
    })


# ── 1. Read processing funnel ─────────────────────────────────────────────────

FUNNEL_KEYS = [
    ("input_reads",          "input_reads"),
    ("high_quality",         "hq_reads"),
    ("parsed_reads",         "parsed_reads"),
    ("post_collapse",        "collapsed_reads"),
    ("bowtie_input",         "bowtie_input"),
    ("aligned",              "aligned"),
    ("suppressed_m",         "suppressed"),
    ("mapped_libsize",       "libsize"),
    ("fc_assigned",          "fc_assigned"),
    ("fc_unassigned",        "fc_unassigned_nofeat"),
]
for s in all_samples:
    for metric, key in FUNNEL_KEYS:
        add("funnel", metric, s, data[s].get(key))
    inp, aln = data[s].get("bowtie_input"), data[s].get("aligned")
    add("funnel", "alignment_rate_pct", s,
        round(aln / inp * 100, 2) if inp and aln else None)
    lib, asg = data[s].get("libsize"), data[s].get("fc_assigned")
    add("funnel", "assignment_rate_pct", s,
        round(asg / lib * 100, 2) if lib and asg else None)


# ── 2. RNA subtype breakdown (unique species + total CPM per type) ────────────

all_types = sorted(set(
    t for s in all_samples if not data[s]["cpm_df"].empty
    for t in data[s]["cpm_df"]["gene_type"].unique()
))
for s in all_samples:
    df = data[s]["cpm_df"]
    for rtype in all_types:
        if df.empty:
            n, cpm = None, None
        else:
            sub = df[df["gene_type"] == rtype]
            n = int(sub["gene_id"].nunique()) if "gene_id" in sub.columns else int(sub.shape[0])
            cpm = round(float(sub["CPM"].sum()), 2)
        add("subtype", "unique_species", s, n, rna_type=rtype)
        add("subtype", "total_cpm",      s, cpm, rna_type=rtype)
    # TOTAL across all types
    if df.empty:
        n_tot, cpm_tot = None, None
    else:
        n_tot = int(df["gene_id"].nunique()) if "gene_id" in df.columns else int(df.shape[0])
        cpm_tot = round(float(df["CPM"].sum()), 2)
    add("subtype", "unique_species", s, n_tot, rna_type="TOTAL")
    add("subtype", "total_cpm",      s, cpm_tot, rna_type="TOTAL")


# ── 3. Quality metrics ────────────────────────────────────────────────────────

for s in all_samples:
    if not spearman_df.empty and s in spearman_df.index:
        add("quality", "spearman_rho",     s, round(float(spearman_df.loc[s, "rho_all"]), 4))
        add("quality", "n_mirnas_matched", s, int(spearman_df.loc[s, "n_matched"]))
    else:
        add("quality", "spearman_rho",     s, None)
        add("quality", "n_mirnas_matched", s, None)
    add("quality", "gini",           s, data[s].get("gini"))
    add("quality", "hhi",            s, data[s].get("hhi"))
    add("quality", "cv",             s, data[s].get("cv"))
    add("quality", "top10_pct_cpm",  s, data[s].get("top10pct"))
    add("quality", "top50_pct_cpm",  s, data[s].get("top50pct"))
    add("quality", "detected_mirnas", s, data[s].get("n_mirna"))


# ── 4. piRNA and cityRNA ──────────────────────────────────────────────────────

for s in all_samples:
    add("pirna", "unique_pirnas",  s, data[s].get("pirna_unique"))
    add("pirna", "reads_matched",  s, data[s].get("pirna_reads"))
    add("cityRNA", "reads_scanned",   s, data[s].get("city_scanned"))
    add("cityRNA", "in_range_unique", s, data[s].get("city_unique_in_range"))
    add("cityRNA", "matches",         s, data[s].get("city_matches"))
    add("cityRNA", "seed_matches",    s, data[s].get("city_seed_matches"))
    add("cityRNA", "seed_miRNAs",     s, data[s].get("city_seed_mirnas"))


# ── Write tidy TSV ────────────────────────────────────────────────────────────

out_df = pd.DataFrame(records, columns=["section", "metric", "rna_type", "sample", "value"])
tsv_path.parent.mkdir(parents=True, exist_ok=True)
out_df.to_csv(tsv_path, sep="\t", index=False)
print(f"Run summary (tidy) ({len(all_samples)} samples, CPM≥{cpm_threshold}) → {tsv_path}  [{len(out_df)} rows]")


# ── Render human-readable report (.txt), driven by the same tidy table ─────────

W  = 14   # per-sample column width
LW = 28   # left label width
TW = 16   # RNA-type label width

# pivots keyed by (metric or (metric, rna_type)) → {sample: value}
_by = {}
for r in records:
    key = (r["metric"], r["rna_type"]) if r["section"] == "subtype" else r["metric"]
    _by.setdefault(key, {})[r["sample"]] = r["value"]

def _cell(v, kind):
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return "N/A".rjust(W)
    if kind == "pct":   return f"{v:.1f}%".rjust(W)
    if kind == "int":   return f"{int(v):,}".rjust(W)
    if kind == "cpm":   return f"{v:,.1f}".rjust(W)
    if kind == "f3":    return f"{v:.3f}".rjust(W)
    if kind == "f4":    return f"{v:.4f}".rjust(W)
    return str(v).rjust(W)

def hrow(label, key, kind="int"):
    vals = _by.get(key, {})
    return label.ljust(LW) + "".join(_cell(vals.get(s), kind) for s in all_samples)

def header():
    return "".ljust(LW) + "".join(s.rjust(W) for s in all_samples)

def bar(w, ch="─"):
    return ch * w

def section(title):
    return f"\n{'═'*80}\n  {title}\n{'═'*80}\n"

L = []
L.append("=" * 80)
L.append("  tinyRNA-seq Pipeline Run Summary")
L.append(f"  Generated : {date.today()}")
L.append(f"  Samples   : {len(all_samples)}  ({', '.join(all_samples)})")
L.append(f"  CPM threshold (filter): {cpm_threshold}")
L.append("=" * 80)

# 1. Funnel
L.append(section("1. READ PROCESSING FUNNEL"))
L.append(header()); L.append(bar(LW + W * len(all_samples)))
L.append(hrow("Input reads",       "input_reads"))
L.append(hrow("  High quality",    "high_quality"))
L.append(hrow("  Parsed reads",    "parsed_reads"))
L.append(hrow("  Post-collapse",   "post_collapse"))
L.append(hrow("Bowtie input",      "bowtie_input"))
L.append(hrow("  Aligned",         "aligned"))
L.append(hrow("  Suppressed (-m)", "suppressed_m"))
L.append(hrow("Mapped (libsize)",  "mapped_libsize"))
L.append(hrow("  FC assigned",     "fc_assigned"))
L.append(hrow("  FC unassigned",   "fc_unassigned"))
L.append("")
L.append(hrow("Alignment rate",    "alignment_rate_pct",  "pct"))
L.append(hrow("Assignment rate",   "assignment_rate_pct", "pct"))
L.append("")
L.append("  Note: FASTQ_SAMPLES (4N-protocol) show full funnel; TRIMMED_SAMPLES")
L.append("  start at Bowtie input (adapter trimming performed externally).")

# 2. RNA subtype breakdown
L.append(section(f"2. RNA SUBTYPE BREAKDOWN  (CPM ≥ {cpm_threshold})"))
tbar = bar(TW + W * len(all_samples))
def type_table(metric, kind, title):
    L.append(f"  {title}"); L.append("")
    L.append("  " + "RNA type".ljust(TW) + "".join(s.rjust(W) for s in all_samples))
    L.append("  " + tbar)
    for rtype in all_types:
        vals = _by.get((metric, rtype), {})
        L.append("  " + rtype.ljust(TW) + "".join(_cell(vals.get(s), kind) for s in all_samples))
    L.append("  " + tbar)
    tot = _by.get((metric, "TOTAL"), {})
    L.append("  " + "TOTAL".ljust(TW) + "".join(_cell(tot.get(s), kind) for s in all_samples))
    L.append("")
type_table("unique_species", "int", "A. Unique species (n genes detected)")
type_table("total_cpm",      "cpm", "B. Total CPM per RNA type")

# 3. Quality metrics
L.append(section("3. QUALITY METRICS"))
L.append(header()); L.append(bar(LW + W * len(all_samples)))
L.append(hrow("Spearman ρ (vs brain ref)", "spearman_rho",     "f3"))
L.append(hrow("  n miRNAs matched",        "n_mirnas_matched", "int"))
L.append("")
L.append(hrow("Gini coefficient",          "gini",            "f3"))
L.append(hrow("HHI",                       "hhi",             "f4"))
L.append(hrow("CV",                        "cv",              "f3"))
L.append(hrow("Top-10 miRNAs (% CPM)",     "top10_pct_cpm",   "pct"))
L.append(hrow("Top-50 miRNAs (% CPM)",     "top50_pct_cpm",   "pct"))
L.append(hrow("Detected miRNAs (n)",       "detected_mirnas", "int"))

# 4. piRNA and cityRNA
L.append(section("4. piRNA  AND  cityRNA"))
L.append("  piRNA (rescued from RM-assigned + unassigned reads):")
L.append(f"  {'Sample':<16}{'Unique piRNAs':>14}{'Reads matched':>14}")
L.append("  " + bar(44))
for s in all_samples:
    u = _by.get("unique_pirnas", {}).get(s)
    r = _by.get("reads_matched", {}).get(s)
    u_s = f"{int(u):,}" if u is not None and not pd.isna(u) else "N/A"
    r_s = f"{int(r):,}" if r is not None and not pd.isna(r) else "N/A"
    L.append(f"  {s:<16}{u_s:>14}{r_s:>14}")
L.append("")
L.append("  cityRNA (12–15 nt reads matching miRNA sequences; seed = read spans miRNA nt 2–8):")
L.append(f"  {'Sample':<16}{'Reads scanned':>14}{'In-range unique':>16}{'Matches':>10}"
         f"{'Seed match':>12}{'Seed miRNAs':>13}")
L.append("  " + bar(81))
for s in all_samples:
    sc = _by.get("reads_scanned", {}).get(s)
    ur = _by.get("in_range_unique", {}).get(s)
    mt = _by.get("matches", {}).get(s)
    sm = _by.get("seed_matches", {}).get(s)
    sn = _by.get("seed_miRNAs", {}).get(s)
    f = lambda v: f"{int(v):,}" if v is not None and not pd.isna(v) else "N/A"
    L.append(f"  {s:<16}{f(sc):>14}{f(ur):>16}{f(mt):>10}{f(sm):>12}{f(sn):>13}")
L.append("")

txt_path.write_text("\n".join(L) + "\n")
print(f"Run summary (rendered) → {txt_path}")
