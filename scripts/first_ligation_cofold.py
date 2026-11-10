#!/usr/bin/env python3
"""
first_ligation_cofold.py

Per-read RNA cofold analysis of the FIRST ligation junction — between the VCG
insert (left strand) and the 3' adaptor structure (right strand).

This simulates the ligation competency of each VCG before the library 5'
adaptor is considered.  It answers: does the VCG 3' end present a good
substrate for 3' adaptor ligation?

Left strand  (ligation junction → 5' end of VCG):
    VCG sequence (RNA, T→U)

Right strand (ligation junction → 3' end):
    [MADAP/4N structure in DNA] + 3_Adap + adaptor_3p_tail

    middle_cols follow the same RIGHT_COLS scheme as secondary_structures.py.
    3_Adap is taken from the vcg_filtered.tsv column (actually-sequenced
    adaptor, so error-tolerant vs. the reference sequence).
    adaptor_3p_tail is resolved from the TSV override → schema field →
    DEFAULT_ADAPTOR_3P_TAIL (same resolution order as secondary_structures).

Note: adaptor_5p_ext is NOT used here — the first ligation junction is
defined solely by the VCG and what is downstream (3' adaptor side).

Cofold classes (Zhuang 2012 extended):
    Classes 1–16 — same integer scheme as secondary_structures.py (second
                   ligation).  Same class number = same junction topology.
    Class 17     — right-strand homodimer only; VCG free (adapt_stem)
    Class 18     — left-strand homodimer only; adaptor free (vcg_stem)
    Class 19     — no base-pairing anywhere (open)
    Class 0      — unclassified

    Ligation-competent reads sit in classes where left_junction = PHe
    (VCG 3' end engaged in cross-strand pairing at the junction):
    classes 1, 3, 9, 12 — and, with adaptor internal structure, class 4.

Outputs:
    first_ligation_reads.tsv    — per-read: all input cols + seq_l, seq_r,
                                   dotbracket, mfe, junction_loop_size,
                                   cofold_class + 7 topology columns
    first_ligation_summary.tsv  — per-VCG: n_reads, mean_mfe, mean_loop_size,
                                   cofold_class frequency counts
    first_ligation/             — figures directory

Called by Snakemake (rule first_ligation_cofold) or directly:
    python scripts/first_ligation_cofold.py \\
        --filtered results/EXP1/EXP1.vcg_filtered.tsv \\
        --schemas adaptor_schemas.yaml \\
        --read-layout FOUR_N_MADAP \\
        --reads-out results/EXP1/EXP1.first_ligation_reads.tsv \\
        --summary-out results/EXP1/EXP1.first_ligation_summary.tsv \\
        --figures-dir results/EXP1/first_ligation/
"""

import argparse
import fcntl
import json
import sys
import tempfile
from pathlib import Path

import pandas as pd
import RNA
import yaml

# ── Import shared classifier from secondary_structures.py ─────────────────────
sys.path.insert(0, str(Path(__file__).parent))
from secondary_structures import (
    _classify_single,
    CLASS_DESCRIPTIONS as _ZHUANG_DESCRIPTIONS,
)

# ── Parallel progress table (same mechanism as secondary_structures.py) ───────

_PROGRESS_FILE = Path(tempfile.gettempdir()) / "snakemake_firstlig_progress.json"
_LOCK_FILE     = Path(tempfile.gettempdir()) / "snakemake_firstlig_progress.lock"
_BAR_WIDTH     = 20


def _redraw_progress_table(progress_out, sample, current, total, unit):
    """Update this sample's entry and atomically redraw the full table."""
    if progress_out is None:
        return
    pct = current / max(total, 1)
    with open(_LOCK_FILE, "w") as lf:
        fcntl.flock(lf, fcntl.LOCK_EX)
        try:
            try:
                state = json.loads(_PROGRESS_FILE.read_text()) if _PROGRESS_FILE.exists() else {}
            except Exception:
                state = {}
            prev_lines = state.pop("__lines__", 0)
            state[sample] = {"n": current, "total": total, "pct": pct, "unit": unit}
            rows = []
            for s, info in sorted(state.items()):
                filled = int(info["pct"] * _BAR_WIDTH)
                bar = "█" * filled + "░" * (_BAR_WIDTH - filled)
                rows.append(
                    f"  {s:<15} [{bar}] {info['pct']:>4.0%}"
                    f"  {info['n']:,}/{info['total']:,} {info['unit']}"
                )
            lines = ["1st-ligation cofold:"] + rows
            if prev_lines:
                progress_out.write(f"\033[{prev_lines}A\033[J")
            progress_out.write("\n".join(lines) + "\n")
            progress_out.flush()
            state["__lines__"] = len(lines)
            _PROGRESS_FILE.write_text(json.dumps(state))
        finally:
            fcntl.flock(lf, fcntl.LOCK_UN)


def _finish_progress_entry(progress_out, sample):
    """Remove this sample from the table when done; clear table if last."""
    if progress_out is None:
        return
    with open(_LOCK_FILE, "w") as lf:
        fcntl.flock(lf, fcntl.LOCK_EX)
        try:
            try:
                state = json.loads(_PROGRESS_FILE.read_text()) if _PROGRESS_FILE.exists() else {}
            except Exception:
                state = {}
            prev_lines = state.pop("__lines__", 0)
            state.pop(sample, None)
            remaining = {k: v for k, v in state.items() if not k.startswith("__")}
            if remaining:
                rows = []
                for s, info in sorted(remaining.items()):
                    filled = int(info["pct"] * _BAR_WIDTH)
                    bar = "█" * filled + "░" * (_BAR_WIDTH - filled)
                    rows.append(
                        f"  {s:<15} [{bar}] {info['pct']:>4.0%}"
                        f"  {info['n']:,}/{info['total']:,} {info['unit']}"
                    )
                lines = ["1st-ligation cofold:"] + rows
                if prev_lines:
                    progress_out.write(f"\033[{prev_lines}A\033[J")
                progress_out.write("\n".join(lines) + "\n")
                progress_out.flush()
                state["__lines__"] = len(lines)
                _PROGRESS_FILE.write_text(json.dumps(state))
            else:
                if prev_lines:
                    progress_out.write(f"\033[{prev_lines}A\033[J")
                progress_out.flush()
                _PROGRESS_FILE.unlink(missing_ok=True)
        finally:
            fcntl.flock(lf, fcntl.LOCK_UN)


PALETTE = ["#000000", "#009E73", "#0072B2", "#56B4E9",
           "#F0E442", "#E69F00", "#D55E00", "#CC79A7"]


# ── Cofold class definitions ───────────────────────────────────────────────────
#
# Classes 1–16: Zhuang 2012 (shared with secondary_structures.py).
#   Left strand  = VCG (positions 0 .. amp-1)
#   Right strand = adaptor structure (positions amp+1 .. end)
#
# Extended classes for topologies not covered by the Zhuang lookup:
#   17 — right homodimer only, no het  (adaptor stem; VCG free)
#   18 — left homodimer only, no het   (VCG stem; adaptor free)
#   19 — no pairing anywhere           (open)
#    0 — unclassified
#
# Ligation-competent state = VCG 3' end engaged cross-strand at the junction
# (left_junction = PHe).  Primary classes: 1, 3, 9.  Class 4 also competent
# (right bifurcated: VCG paired at junction, adaptor internally folded).

COFOLD_CLASSES = list(range(20))

COFOLD_CLASS_DESCRIPTIONS = {
    **_ZHUANG_DESCRIPTIONS,
    17: "Adaptor stem only (right homodimer; VCG free)",
    18: "VCG stem only (left homodimer; adaptor free)",
    19: "Open — no base-pairing",
}


def classify_first_ligation(structure: str) -> int:
    """Classify a first-ligation cofold dot-bracket into an integer class.

    Applies _classify_single (Zhuang 2012) and falls back to classes 17/18/19
    for single-strand homodimers and fully open structures that are not
    covered by the Zhuang CLASS_LOOKUP.

    Returns an integer 0–19.
    """
    result = _classify_single(structure)
    cls = result["cofold_class"]
    if cls != 0:
        return cls

    # Fallback for structures that return 0 from the Zhuang lookup.
    left_homo = result["left_homodimer"]
    right_homo = result["right_homodimer"]
    left_het   = result["left_heterodimer"]

    if left_het:
        return 0  # has het but no Zhuang match — genuinely unclassified
    if right_homo and not left_homo:
        return 17  # adaptor stem only
    if left_homo and not right_homo:
        return 18  # VCG stem only
    if not left_homo and not right_homo:
        return 19  # open
    return 0  # both homodimers, no het, not in Zhuang lookup


# ── Schema and layout loading ──────────────────────────────────────────────────

DEFAULT_ADAPTOR_3P_TAIL = "AAUGAATC"   # mirrors secondary_structures.py

DEFAULT_COFOLD_TEMPERATURE = 16.0   # °C — ligation reactions run at 16 °C (override per run/sample)

# Right-strand intermediate column order per 3' adaptor type (all DNA)
RIGHT_COLS = {
    "NORMAL": [],
    "FOUR_N": ["N4"],
    "MADAP":  ["TRR", "BR", "MAC", "IRR", "SP"],
}


def load_layout(schemas_file: str, read_layout: str) -> tuple[dict, str, str, str]:
    """Returns (layout_cfg, five_type, three_type, adaptor_3p_tail)."""
    with open(schemas_file) as fh:
        s = yaml.safe_load(fh)
    layout    = s["read_layouts"][read_layout]
    three_type = s["three_prime_adaptors"][layout["three_prime"]].get("type", layout["three_prime"])
    # Resolve adaptor_3p_tail: schema layout field → DEFAULT
    adaptor_3p_tail = layout.get("adaptor_3p_tail", DEFAULT_ADAPTOR_3P_TAIL)
    return layout, layout["five_prime"], layout["three_prime"], adaptor_3p_tail


# ── Strand reconstruction ──────────────────────────────────────────────────────

def build_strands(
    row: pd.Series,
    three_type: str,
    adaptor_3p_tail: str,
) -> tuple[str, str]:
    """
    Build (seq_l, seq_r) for first-ligation RNA.cofold.

    seq_l = VCG (RNA)
    seq_r = [middle cols as DNA] + 3_Adap + adaptor_3p_tail
    """
    seq_l = str(row["VCG"]).replace("T", "U").replace("t", "u")

    middle_cols = RIGHT_COLS.get(three_type, [])
    middle = "".join(str(row[c]) for c in middle_cols if c in row.index)

    seq_r = middle + str(row.get("3_Adap", "")) + adaptor_3p_tail
    return seq_l, seq_r


# ── Junction analysis ──────────────────────────────────────────────────────────

def junction_loop_size(db: str, junction_index: int) -> int:
    """Count unpaired bases immediately flanking the '&' in the dot-bracket."""
    def dots_left():
        count, i = 0, junction_index - 1
        while i >= 0:
            if db[i] == ".":
                count += 1; i -= 1
            else:
                break
        return count

    def dots_right():
        count, i = 0, junction_index + 1
        while i < len(db):
            if db[i] == ".":
                count += 1; i += 1
            else:
                break
        return count

    if db[junction_index] != ".":
        return dots_left()
    return dots_left() + dots_right() + 1


def _cofold_strands(seq_l: str, seq_r: str) -> dict:
    """Run RNA.cofold on a pre-built strand pair and return result dict."""
    combined = seq_l + "&" + seq_r
    dot_bracket, mfe = RNA.cofold(combined)
    if "&" not in dot_bracket:
        dot_bracket = dot_bracket[:len(seq_l)] + "&" + dot_bracket[len(seq_l):]
    amp_pos   = dot_bracket.index("&")
    loop_size = junction_loop_size(dot_bracket, amp_pos)
    cls_result = _classify_single(dot_bracket)
    cofold_cls = classify_first_ligation(dot_bracket)
    return {
        "seq_l":              seq_l,
        "seq_r":              seq_r,
        "dotbracket":         dot_bracket,
        "mfe":                round(mfe, 4),
        "junction_loop_size": loop_size,
        "left_homodimer":     cls_result["left_homodimer"],
        "left_heterodimer":   cls_result["left_heterodimer"],
        "right_homodimer":    cls_result["right_homodimer"],
        "right_heterodimer":  cls_result["right_heterodimer"],
        "left_junction":      cls_result["left_junction"],
        "right_junction":     cls_result["right_junction"],
        "cofold_class":       cofold_cls,
    }




# ── Summary ────────────────────────────────────────────────────────────────────

def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    base = df.groupby(["VCG_Len", "VCG"], as_index=False).agg(
        n_reads=("mfe", "count"),
        mean_mfe=("mfe", "mean"),
        mean_loop_size=("junction_loop_size", "mean"),
    )
    base["mean_mfe"]       = base["mean_mfe"].round(4)
    base["mean_loop_size"] = base["mean_loop_size"].round(2)

    class_counts = (
        df.groupby(["VCG", "cofold_class"])
        .size()
        .unstack(fill_value=0)
        .add_prefix("class_")
        .reset_index()
    )
    summary = base.merge(class_counts, on="VCG", how="left")
    return summary.sort_values(["VCG_Len", "VCG"]).reset_index(drop=True)


# ── Core ───────────────────────────────────────────────────────────────────────

def first_ligation_cofold(
    filtered_path: str,
    schemas_file: str,
    read_layout: str,
    reads_out: str,
    summary_out: str,
    sample: str = "",
    max_reads_per_vcg: int = 2000,
    progress_out=None,
    adaptor_3p_tail_override: str | None = None,
    cofold_temperature: float = DEFAULT_COFOLD_TEMPERATURE,
) -> None:
    RNA.cvar.temperature = cofold_temperature
    layout_cfg, five_type, three_type, adaptor_3p_tail = load_layout(schemas_file, read_layout)
    if adaptor_3p_tail_override is not None:
        adaptor_3p_tail = adaptor_3p_tail_override

    Path(reads_out).parent.mkdir(parents=True, exist_ok=True)

    tag = f"[first_ligation_cofold:{sample}]"
    print(f"{tag} cofold temperature: {cofold_temperature} °C", flush=True)
    print(f"{tag} loading {filtered_path} ...", flush=True)
    df = pd.read_csv(filtered_path, sep="\t", dtype=str)
    n_total = len(df)

    # Per-VCG subsampling (same threshold as secondary_structures)
    SUBSAMPLE_THRESHOLD = 500_000
    if n_total > SUBSAMPLE_THRESHOLD:
        df = (
            df.groupby("VCG", group_keys=False)
            .apply(lambda g: g.sample(min(len(g), max_reads_per_vcg), random_state=42))
            .reset_index(drop=True)
        )
        n = len(df)
        print(
            f"{tag} Layout: {read_layout}  |  "
            f"Reads: {n:,} (subsampled to {max_reads_per_vcg}/VCG; input was {n_total:,})"
        )
    else:
        n = n_total
        print(f"{tag} Layout: {read_layout}  |  Reads: {n:,}")

    seq_ls, seq_rs = zip(*[
        build_strands(row, three_type, adaptor_3p_tail)
        for _, row in df.iterrows()
    ]) if n > 0 else ([], [])

    if n > 0:
        df["_seq_l"] = seq_ls
        df["_seq_r"] = seq_rs
        unique_pairs = df[["_seq_l", "_seq_r"]].drop_duplicates().reset_index(drop=True)
        n_unique = len(unique_pairs)
        msg = f"{tag} cofolding {n_unique:,} unique pairs (from {n:,} reads)...\n"
        print(msg.strip(), flush=True)
        if progress_out is not None:
            progress_out.write(msg); progress_out.flush()

        cache = {}
        report_every = max(1, n_unique // 20)
        for i, (_, urow) in enumerate(unique_pairs.iterrows()):
            key = (urow["_seq_l"], urow["_seq_r"])
            cache[key] = _cofold_strands(urow["_seq_l"], urow["_seq_r"])
            if (i + 1) % report_every == 0:
                _redraw_progress_table(progress_out, sample, i + 1, n_unique, "pairs")
        _redraw_progress_table(progress_out, sample, n_unique, n_unique, "pairs")

        result_rows = [cache[(row["_seq_l"], row["_seq_r"])] for _, row in df.iterrows()]
        df = df.drop(columns=["_seq_l", "_seq_r"])
    else:
        result_rows = []

    _finish_progress_entry(progress_out, sample)
    if progress_out is not None:
        progress_out.write(f"{tag} cofold complete\n")
        progress_out.flush()

    result_df = pd.concat([df.reset_index(drop=True), pd.DataFrame(result_rows)], axis=1) if n > 0 else df
    if n > 0:
        result_df["mfe"]                = pd.to_numeric(result_df["mfe"])
        result_df["junction_loop_size"] = pd.to_numeric(result_df["junction_loop_size"])
        result_df["VCG_Len"]            = pd.to_numeric(result_df["VCG_Len"])

    if n > 0:
        class_counts = result_df["cofold_class"].value_counts().sort_index()
        print(f"{tag} Cofold class distribution ({n:,} reads):")
        for cls, count in class_counts.items():
            desc = COFOLD_CLASS_DESCRIPTIONS.get(cls, "")
            pct  = f"{count / n:.1%}"
            print(f"{tag}   class {cls:<3d}  {count:>8,}  ({pct})  {desc}")

    result_df.to_csv(reads_out, sep="\t", index=False)
    print(f"{tag} Reads → {reads_out}")

    if n > 0:
        summary = build_summary(result_df)
        summary.to_csv(summary_out, sep="\t", index=False)
        print(f"{tag} Summary → {summary_out}")
    else:
        pd.DataFrame(columns=["VCG_Len", "VCG", "n_reads", "mean_mfe", "mean_loop_size"]).to_csv(
            summary_out, sep="\t", index=False
        )
        print(f"{tag} No reads — empty outputs written.")


# ── Entry points ───────────────────────────────────────────────────────────────

if "snakemake" in dir():
    _term = sys.stderr
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    first_ligation_cofold(
        filtered_path          = snakemake.input.filtered,
        schemas_file           = snakemake.input.schemas,
        read_layout            = snakemake.params.read_layout,
        reads_out              = snakemake.output.reads,
        summary_out            = snakemake.output.summary,
        sample                 = snakemake.wildcards.sample,
        max_reads_per_vcg      = snakemake.params.max_reads_per_vcg,
        progress_out           = _term,
        adaptor_3p_tail_override = snakemake.params.adaptor_3p_tail,
        cofold_temperature     = float(snakemake.params.cofold_temperature),
    )

elif __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Per-read RNA cofold — first ligation junction (VCG vs 3' adaptor).")
    ap.add_argument("--filtered",          required=True)
    ap.add_argument("--schemas",           default="adaptor_schemas.yaml")
    ap.add_argument("--read-layout",       required=True)
    ap.add_argument("--reads-out",         required=True)
    ap.add_argument("--summary-out",       required=True)
    ap.add_argument("--max-reads-per-vcg", type=int, default=2000)
    ap.add_argument("--adaptor-3p-tail",   default=None,
                    help="Override adaptor_3p_tail (use 'None' for explicit empty)")
    ap.add_argument("--cofold-temperature", type=float, default=DEFAULT_COFOLD_TEMPERATURE)
    args = ap.parse_args()
    tail_override = "" if args.adaptor_3p_tail == "None" else args.adaptor_3p_tail
    first_ligation_cofold(
        args.filtered, args.schemas, args.read_layout,
        args.reads_out, args.summary_out,
        max_reads_per_vcg=args.max_reads_per_vcg,
        adaptor_3p_tail_override=tail_override,
        cofold_temperature=args.cofold_temperature,
    )
