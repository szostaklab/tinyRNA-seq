#!/usr/bin/env python3
"""
rarefied_comparison.py

Opt-in rarefied (equal-depth) cross-sample comparison for the BIO track. A depth-fair
supplement to the standard full-depth comparison: downsamples each bio sample's aligned
reads to a common depth (fixed seed) and re-counts, so detection/richness is compared at
equal sequencing depth. Runs against an existing results/ tree; the full-depth comparison
is left untouched. Enabled via config modules.comparison_downsample (off by default).

For every bio sample it:
  1. downsamples the aligned reads (NH.sorted.bam) to a common depth (default: the
     smallest bio library, so nothing is upsampled), with a fixed seed,
  2. re-runs featureCounts (same flags as the pipeline) on the downsampled BAM,
  3. re-runs assign_piRNA on the downsampled BAM (piRNA is rescued from the same
     reads, so it is depth-matched for free),
  4. rebuilds the smallRNA + CPM table at the matched depth.

Then it writes the rarefied per-sample tables and a rarefied comparison
(unique species and total CPM per biotype), plus a full-vs-rarefied delta so the
depth effect is visible.

cityRNA is rarefied SEPARATELY, on its own denominator: short_reads never enter the
BAM, so cityRNA richness scales with the number of sub-18 nt reads scanned, not with
mapped depth. It is downsampled to a common SCANNED depth (city_finder_depth: auto =
smallest short-read pool | integer | a bio sample name), the same seed, then the real
city_finder is re-run and matches / seed-matches / seed-miRNAs are re-counted at equal
scanned depth. Reusing the aligned fraction here would be wrong (pool sizes are not
proportional to mapped depth).

Usage (from the pipeline root):
    python scripts/rarefied_comparison.py                 # auto depth = min library
    python scripts/rarefied_comparison.py --depth 250000  # fixed floor; drop samples below it
    python scripts/rarefied_comparison.py --seed 7 --keep-bam
"""
import argparse
import glob
import os
import random
import subprocess
import sys
from pathlib import Path

import pandas as pd

try:
    import yaml
except ImportError:
    sys.exit("needs pyyaml (conda env with snakemake has it)")


# ── config / manifest ──────────────────────────────────────────────────────────

def load_config(path):
    with open(path) as fh:
        return yaml.safe_load(fh)

def bio_samples(experiments_tsv):
    df = pd.read_csv(experiments_tsv, sep="\t", dtype=str).fillna("")
    bio = df[df["experiment_type"] == "bio"]
    return list(zip(bio["experiment_name"], bio["organism"]))

def resolve_gtf(gtf_path):
    if os.path.exists(gtf_path):
        return gtf_path
    d = os.path.dirname(gtf_path) or "."
    hits = sorted(glob.glob(f"{d}/*.gtf") + glob.glob(f"{d}/*.gtf.gz"))
    if len(hits) != 1:
        sys.exit(f"gtf: expected exactly one *.gtf in {d}, found {len(hits)}")
    return hits[0]


# ── per-sample downsample + recount ────────────────────────────────────────────

def sh(cmd, log=None):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        if log:
            Path(log).write_text(r.stdout + r.stderr)
        sys.exit(f"FAILED: {cmd}\n{r.stderr[-2000:]}")
    return r.stdout

def downsample_bam(bam, frac, seed, out_bam, threads):
    sf = f"{seed + frac:.6f}"   # samtools -s SEED.FRACTION (e.g. 13.838360)
    sh(f'samtools view -@ {threads} -b -s {sf} "{bam}" > "{out_bam}"')
    sh(f'samtools index "{out_bam}"')
    return int(sh(f'samtools view -c -F 4 "{out_bam}"').strip())

def featurecounts(bam, gtf, out_txt, threads, log):
    sh(f'featureCounts -T {threads} -a "{gtf}" -t exon -g gene_id '
       f'-s 0 -O --fraction -R CORE -o "{out_txt}" "{bam}"', log=log)

def run_pirna(script, bam, fc_core, refs, out_tsv, sample, log):
    refs_s = " ".join(f'"{r}"' for r in refs)
    sh(f'python "{script}" --bam "{bam}" --fc "{fc_core}" --refs {refs_s} '
       f'--out "{out_tsv}" --sample "{sample}"', log=log)


# ── table building ──────────────────────────────────────────────────────────────

def gene_biotype_map(gene_info_tsv):
    m = {}
    for l in open(gene_info_tsv):
        p = l.rstrip("\n").split("\t")
        if len(p) >= 3:
            m[p[0]] = (p[1], p[2])   # gene_id -> (name, biotype)
    return m

def build_smallrna_cpm(gene_txt, g2nb, pirna_tsv, libsize, small_types):
    rows = []
    started = False
    for l in open(gene_txt):
        if l.startswith("#"):
            continue
        p = l.rstrip("\n").split("\t")
        if not started:            # column-header line
            started = True; continue
        gid = p[0]; cnt = float(p[-1])
        name, bt = g2nb.get(gid, ("", ""))
        if bt in small_types and cnt > 0:
            rows.append((gid, cnt, cnt / libsize * 1e6, name, bt))
    if pirna_tsv and os.path.exists(pirna_tsv):
        for l in open(pirna_tsv):
            p = l.rstrip("\n").split("\t")
            if len(p) >= 4:
                try:
                    c = float(p[1])
                    rows.append((p[0], c, c / libsize * 1e6, p[2], "piRNA"))
                except ValueError:
                    pass
    return pd.DataFrame(rows, columns=["gene_id", "count", "CPM", "gene_name", "gene_type"])

def unique_by_biotype(df, thr=1.0):
    d = df[df["CPM"] >= thr]
    return d.groupby("gene_type")["gene_id"].nunique()


# ── cityRNA rarefaction (separate short-read pool, own common scanned depth) ─────

def find_short_reads(res_dir, sample):
    for ext in ("fasta", "fastq", "fa", "fq"):
        p = res_dir / sample / f"{sample}.short_reads.{ext}"
        if p.exists():
            return str(p)
    return None

def _is_fastq(path):
    return str(path).lower().rstrip(".gz").endswith(("fastq", "fq"))

def count_records(path):
    n = sum(1 for _ in open(path))
    return n // (4 if _is_fastq(path) else 2)

def subsample_reads(src, dst, n, seed):
    """Record-level subsample of a 2-line FASTA / 4-line FASTQ short_reads file."""
    per = 4 if _is_fastq(src) else 2
    lines = open(src).readlines()
    recs = ["".join(lines[i:i + per]) for i in range(0, len(lines), per)]
    if len(recs) <= n:
        keep = recs
    else:
        idx = sorted(random.Random(seed).sample(range(len(recs)), n))
        keep = [recs[i] for i in idx]
    with open(dst, "w") as out:
        out.writelines(keep)
    return len(keep)

def run_city(script, inp, mirna_csv, out_tsv, min_len, max_len, mm, log):
    sh(f'python "{script}" --fasta "{inp}" --mirna "{mirna_csv}" --out "{out_tsv}" '
       f'--min-len {min_len} --max-len {max_len} --mismatches {mm}', log=log)

def city_metrics(tsv):
    """(matches, seed_matches, seed_miRNAs) from a city_finder TSV, run_summary seed rule."""
    d = pd.read_csv(tsv, sep="\t")
    if d.empty:
        return 0, 0, 0
    pos = pd.to_numeric(d["mirna_match_pos"], errors="coerce")
    ln  = pd.to_numeric(d["length"], errors="coerce")
    seed = (pos <= 1) & (pos + ln - 1 >= 7)
    return len(d), int(seed.sum()), int(d.loc[seed, "mirna_name"].nunique())


# ── main ────────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--experiments", default="experiments.tsv")
    ap.add_argument("--results", default="results")
    ap.add_argument("--out", default="results/comparison/rarefied")
    ap.add_argument("--depth", default=None,
                    help="'auto' (min library), an integer, or a bio sample name "
                         "(cap target at that sample's mapped depth)")
    ap.add_argument("--city-depth", dest="city_depth", default=None,
                    help="cityRNA common SCANNED depth: 'auto' (smallest short-read "
                         "pool), an integer, or a bio sample name")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--keep-bam", action="store_true")
    ap.add_argument("--keep-below", action="store_true",
                    help="keep samples below target at their native depth "
                         "(reuse their full-depth table) instead of excluding them")
    args = ap.parse_args()

    cfg = load_config(args.config)
    defaults = cfg.get("defaults", {})
    small_types = set(defaults.get("small_rna_types", []))
    seed = args.seed if args.seed is not None else int(defaults.get("downsample_seed", 13))
    depth_cfg = args.depth if args.depth is not None else str(defaults.get("comparison_depth", "auto"))

    samples = bio_samples(args.experiments)
    if len(samples) < 2:
        sys.exit("need >=2 bio samples to compare")

    RES = Path(args.results)
    OUT = Path(args.out); OUT.mkdir(parents=True, exist_ok=True)
    pirna_script = "scripts/assign_piRNA.py"

    # collect full libsizes + references
    info = {}
    for s, org in samples:
        lib = int((RES / s / f"{s}.libsize.txt").read_text().strip())
        oblk = cfg["organisms"][org]
        info[s] = dict(org=org, full_lib=lib,
                       bam=str(RES / s / f"{s}.NH.sorted.bam"),
                       gene_info=str(RES / s / f"{s}.gene_info.mm39.tsv"),
                       gtf=resolve_gtf(oblk["gtf"]),
                       pirna_refs=oblk.get("pirna_refs", []))

    # target depth: 'auto' (min library), an integer, or a bio sample name
    min_lib = min(v["full_lib"] for v in info.values())
    if depth_cfg == "auto":
        target = min_lib
    elif depth_cfg in info:
        target = info[depth_cfg]["full_lib"]   # cap at a named sample's depth
    else:
        target = int(depth_cfg)
    if args.keep_below:
        kept, dropped = dict(info), []
    else:
        kept = {s: v for s, v in info.items() if v["full_lib"] >= target}
        dropped = [s for s in info if s not in kept]
    for s in dropped:
        print(f"[warn] {s} full depth {info[s]['full_lib']:,} < target {target:,} — EXCLUDED")
    if len(kept) < 2:
        sys.exit("fewer than 2 samples at/above target depth")

    mode = ("auto=min library" if depth_cfg == "auto"
            else f"cap at {depth_cfg}" if depth_cfg in info else "fixed")
    print(f"[downsample] seed={seed}  target={target:,} mapped reads  ({mode})"
          f"{'  [keep-below: native depth]' if args.keep_below else ''}")

    # per-sample downsample + recount
    ds_lib, ds_tables, native = {}, {}, []
    for s, v in kept.items():
        frac = target / v["full_lib"]
        if frac >= 1.0:   # at/below target → keep native depth (reuse full-depth table)
            native.append(s)
            ds_lib[s] = v["full_lib"]
            full_tbl = pd.read_csv(RES / s / f"{s}.counts.smallRNA.withCPM.tsv", sep="\t")
            full_tbl.to_csv(OUT / f"{s}.counts.smallRNA.withCPM.ds.tsv", sep="\t", index=False)
            ds_tables[s] = full_tbl
            print(f"[{s}] native depth {v['full_lib']:,} (<= target) — kept full", flush=True)
            continue
        print(f"[{s}] full={v['full_lib']:,}  frac={frac:.4f}", flush=True)
        dsbam = str(OUT / f"{s}.ds.bam")
        ds_lib[s] = downsample_bam(v["bam"], frac, seed, dsbam, args.threads)
        gtxt = str(OUT / f"{s}.ds.gene.txt")
        featurecounts(dsbam, v["gtf"], gtxt, args.threads, str(OUT / f"{s}.fc.log"))
        pir = ""
        if v["pirna_refs"]:
            pir = str(OUT / f"{s}.ds.pirna.tsv")
            run_pirna(pirna_script, dsbam, f"{dsbam}.featureCounts",
                      v["pirna_refs"], pir, s, str(OUT / f"{s}.pirna.log"))
        tbl = build_smallrna_cpm(gtxt, gene_biotype_map(v["gene_info"]),
                                 pir, ds_lib[s], small_types)
        tbl.to_csv(OUT / f"{s}.counts.smallRNA.withCPM.ds.tsv", sep="\t", index=False)
        ds_tables[s] = tbl
        if not args.keep_bam:
            for ext in ("", ".bai", ".featureCounts"):
                try: os.remove(dsbam + ext)
                except OSError: pass
        print(f"[{s}] done  ds_lib={ds_lib[s]:,}", flush=True)

    # comparison: unique species per biotype, full vs rarefied
    ks = list(kept)
    full_u = {s: unique_by_biotype(
        pd.read_csv(RES / s / f"{s}.counts.smallRNA.withCPM.tsv", sep="\t")) for s in ks}
    down_u = {s: unique_by_biotype(ds_tables[s]) for s in ks}
    down_cpm = {s: ds_tables[s][ds_tables[s]["CPM"] >= 1].groupby("gene_type")["CPM"].sum()
                for s in ks}

    biotypes = sorted(set().union(*[set(full_u[s].index) | set(down_u[s].index) for s in ks]))
    uniq = pd.DataFrame({s: down_u[s] for s in ks}).reindex(biotypes).fillna(0).astype(int)
    uniq.index.name = "gene_type"
    uniq.to_csv(OUT / "unique_species_by_biotype.rarefied.tsv", sep="\t")
    pd.DataFrame({s: down_cpm[s] for s in ks}).reindex(biotypes).round(1).to_csv(
        OUT / "total_cpm_by_biotype.rarefied.tsv", sep="\t")

    # full-vs-rarefied delta (tidy)
    rec = []
    for s in ks:
        for bt in biotypes:
            f = int(full_u[s].get(bt, 0)); d = int(down_u[s].get(bt, 0))
            rec.append(dict(sample=s, gene_type=bt, full=f, rarefied=d, delta=d - f))
    pd.DataFrame(rec).to_csv(OUT / "full_vs_rarefied.unique_species.tsv", sep="\t", index=False)

    # summary
    lines = ["# Rarefied bio comparison (equal-depth supplement)",
             f"seed={seed}  target={target:,}  depth_mode={depth_cfg}",
             "", "sample        full_lib     ds_lib", "-" * 36]
    for s in ks:
        flag = "  (native, kept full)" if s in native else ""
        lines.append(f"{s:12s}{info[s]['full_lib']:>11,}{ds_lib[s]:>11,}{flag}")
    if native:
        lines += ["", "native depth (kept below target): " + ", ".join(native)]
    if dropped:
        lines += ["", "excluded (below target): " + ", ".join(dropped)]
    lines += ["", "unique smallRNA species (CPM>=1), full -> rarefied:"]
    for s in ks:
        tf = int(sum(full_u[s])); td = int(sum(down_u[s]))
        lines.append(f"  {s:12s} {tf:>6} -> {td:>6}  ({(td-tf)/tf*100:+.1f}%)")
    (OUT / "rarefaction_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n[downsample] outputs -> {OUT}/")

    # ── cityRNA rarefaction (own denominator: sub-18 nt reads scanned) ────────────
    city_min = int(defaults.get("city_finder_min_len", 12))
    city_max = int(defaults.get("city_finder_max_len", 15))
    city_mm  = int(defaults.get("city_finder_mismatches", 0))
    city_depth_cfg = (args.city_depth if args.city_depth is not None
                      else str(defaults.get("city_finder_depth", "auto")))
    city_script = "scripts/city_finder.py"

    cinfo = {}
    for s, org in samples:
        sr = find_short_reads(RES, s)
        mcsv = cfg["organisms"][org].get("mirna_csv")
        ntsv = RES / s / f"{s}.city_finder.tsv"
        if sr and mcsv and ntsv.exists():
            cinfo[s] = dict(short_reads=sr, mirna_csv=mcsv,
                            scanned=count_records(sr), native_tsv=str(ntsv))

    if len(cinfo) < 2:
        print("[cityRNA] skipped (need >=2 bio samples with short_reads + mirna_csv)")
        return

    cmin = min(v["scanned"] for v in cinfo.values())
    if city_depth_cfg == "auto":
        ctarget = cmin
    elif city_depth_cfg in cinfo:
        ctarget = cinfo[city_depth_cfg]["scanned"]
    else:
        ctarget = int(city_depth_cfg)
    keep_below = args.keep_below or bool(defaults.get("comparison_keep_below", True))
    if keep_below:
        ckept, cdropped = dict(cinfo), []
    else:
        ckept = {s: v for s, v in cinfo.items() if v["scanned"] >= ctarget}
        cdropped = [s for s in cinfo if s not in ckept]
    for s in cdropped:
        print(f"[cityRNA][warn] {s} scanned {cinfo[s]['scanned']:,} < target {ctarget:,} — EXCLUDED")
    if len(ckept) < 2:
        print("[cityRNA] skipped (fewer than 2 short-read pools at/above target)")
        return

    cmode = ("auto=smallest pool" if city_depth_cfg == "auto"
             else f"cap at {city_depth_cfg}" if city_depth_cfg in cinfo else "fixed")
    print(f"[cityRNA] seed={seed}  target={ctarget:,} scanned reads  ({cmode})"
          f"{'  [keep-below: native depth]' if keep_below else ''}")

    crows = []
    for s, v in ckept.items():
        nat_m, nat_sm, nat_snm = city_metrics(v["native_tsv"])
        if v["scanned"] <= ctarget:            # native — reuse existing city_finder.tsv
            sc, m, sm, snm, is_nat = v["scanned"], nat_m, nat_sm, nat_snm, True
        else:
            sub = str(OUT / f"{s}.short_reads.rare{Path(v['short_reads']).suffix}")
            sc = subsample_reads(v["short_reads"], sub, ctarget, seed)
            dtsv = str(OUT / f"{s}.city_finder.ds.tsv")
            run_city(city_script, sub, v["mirna_csv"], dtsv,
                     city_min, city_max, city_mm, str(OUT / f"{s}.city.log"))
            m, sm, snm = city_metrics(dtsv)
            is_nat = False
            if not args.keep_bam:
                try: os.remove(sub)
                except OSError: pass
        crows.append(dict(sample=s, native=is_nat,
                          scanned_full=v["scanned"], scanned_rarefied=sc,
                          matches_full=nat_m, matches_rarefied=m,
                          seed_matches_full=nat_sm, seed_matches_rarefied=sm,
                          seed_miRNAs_full=nat_snm, seed_miRNAs_rarefied=snm))
        print(f"[cityRNA:{s}] scanned {v['scanned']:,}->{sc:,}  matches {nat_m}->{m}  "
              f"seed {nat_sm}->{sm}  seed_miRNAs {nat_snm}->{snm}", flush=True)

    pd.DataFrame(crows).to_csv(OUT / "cityRNA_by_sample.rarefied.tsv", sep="\t", index=False)
    with open(OUT / "rarefaction_summary.txt", "a") as fh:
        fh.write(f"\ncityRNA (common scanned depth = {ctarget:,}, seed={seed}, "
                 f"mode={city_depth_cfg}):\n")
        fh.write(f"  {'sample':12s} {'scanned':>18}  {'matches':>14}  {'seed_miRNAs':>14}\n")
        for r in crows:
            fh.write(f"  {r['sample']:12s} "
                     f"{r['scanned_full']:>8,}->{r['scanned_rarefied']:>8,}  "
                     f"{r['matches_full']:>6}->{r['matches_rarefied']:>6}  "
                     f"{r['seed_miRNAs_full']:>6}->{r['seed_miRNAs_rarefied']:>6}\n")
    print(f"[cityRNA] outputs -> {OUT}/cityRNA_by_sample.rarefied.tsv")


if __name__ == "__main__":
    main()
