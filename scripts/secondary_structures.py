#!/usr/bin/env python3
"""
secondary_structures.py

Per-read RNA cofold analysis of the ligation junction between the 5' library
adaptor and the VCG insert + 3' adaptor structure.

For each read in vcg_filtered.tsv:
  1. Reconstruct the two strands flanking the ligation junction.
  2. Run RNA.cofold to get the dot-bracket secondary structure and MFE.
  3. Measure the junction loop size (unpaired bases at the ligation site).
  4. Classify the dot-bracket into one of 16 structural classes (Zhuang 2012).

Left strand  (5' library adaptor → ligation junction):
    ADAPTOR_5P_EXT + 5_Adap   (random layouts)
    ADAPTOR_5P_EXT + TCTA     (DEFINED layout)

Right strand (ligation junction → 3' end):
    VCG(RNA) + [MADAP/4N structure in DNA] + 3_Adap + ADAPTOR_3P_TAIL

VCG is output as RNA (T→U). All adaptor structure columns remain as DNA.
ViennaRNA handles mixed T/U notation natively.

Outputs:
    secondary_struct_reads.tsv   — per-read: all input cols + dotbracket, mfe,
                                   junction_loop_size, + 7 cofold_class columns
    secondary_struct_summary.tsv — per-VCG: n_reads, mean_mfe, mean_loop_size,
                                   class frequency counts
    secondary_structures/        — figures directory

Called by Snakemake (rule secondary_structures) or directly:
    python scripts/secondary_structures.py \\
        --filtered vcg_filtered.tsv \\
        --schemas adaptor_schemas.yaml \\
        --read-layout FOUR_N_MADAP \\
        --reads-out secondary_struct_reads.tsv \\
        --summary-out secondary_struct_summary.tsv \\
        --figures-dir secondary_structures/
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

# ── Parallel progress table ───────────────────────────────────────────────────
# Each Snakemake job writes its progress to a shared JSON file (locked via
# fcntl) and then redraws a multi-row table on the terminal using ANSI
# escape codes.  The file lock ensures only one process writes to the
# terminal at a time, preventing interleaved output.

_PROGRESS_FILE = Path(tempfile.gettempdir()) / "snakemake_cofold_progress.json"
_LOCK_FILE     = Path(tempfile.gettempdir()) / "snakemake_cofold_progress.lock"
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
            lines = ["Cofolding:"] + rows
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
                lines = ["Cofolding:"] + rows
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


# ── Cofold classifier (inlined from Zhuang 2012 NAR 40(7):e54) ───────────────
#
# Classifies RNAcofold dot-bracket strings into 16 structural classes based on
# the local topology at the RNA–adaptor ligation junction.
# Key: (left_homo, left_het, right_homo, right_het, left_junc, right_junc)

CLASS_LOOKUP = {
    # Dual heterodimer
    (False, True,  False, True,  'PHe', 'PHe'): 1,
    (False, True,  False, True,  'UB',  'UB' ): 5,
    (False, True,  False, True,  'PHe', 'UB' ): 9,
    (False, True,  False, True,  'UB',  'PHe'): 13,
    # Dual homodimer
    (True,  False, True,  False, 'PH',  'PH' ): 2,
    (True,  False, True,  False, 'U',   'U'  ): 6,
    (True,  False, True,  False, 'PH',  'U'  ): 10,
    (True,  False, True,  False, 'U',   'PH' ): 14,
    # Mixed, both ends paired
    (True,  True,  False, True,  'PH',  'PHe'): 3,
    (False, True,  True,  True,  'PHe', 'PH' ): 4,
    # Mixed, both ends free (bifurcated)
    (True,  True,  False, True,  'UB',  'UB' ): 7,
    (False, True,  True,  True,  'UB',  'UB' ): 8,
    # Left bifurcated
    (True,  True,  False, True,  'PH',  'UB' ): 11,
    (True,  True,  False, True,  'UB',  'PHe'): 15,
    # Right bifurcated
    (False, True,  True,  True,  'PHe', 'UB' ): 12,
    (False, True,  True,  True,  'UB',  'PH' ): 16,
}

CLASS_DESCRIPTIONS = {
    1:  "Full heterodimer, blunt junction",
    2:  "Dual homodimer, both ends paired",
    3:  "Left bifurcated, both ends paired",
    4:  "Right bifurcated, both ends paired",
    5:  "Full heterodimer, both ends looped",
    6:  "Dual homodimer, both ends free",
    7:  "Left bifurcated, both ends free",
    8:  "Right bifurcated, both ends free",
    9:  "Heterodimer, right end looped",
    10: "Dual homodimer, left end paired / right free",
    11: "Left bifurcated, left junction in homodimer",
    12: "Right bifurcated, right junction in heterodimer",
    13: "Heterodimer, left end looped",
    14: "Dual homodimer, left free / right end paired",
    15: "Left bifurcated, left junction looped",
    16: "Right bifurcated, right junction in homodimer",
    0:  "Unclassified",
}


def _get_bracket_pairs(structure: str) -> dict:
    stack, pairs = [], {}
    for i, c in enumerate(structure):
        if c in '.&':
            continue
        if c == '(':
            stack.append(i)
        elif c == ')' and stack:
            j = stack.pop()
            pairs[i] = j
            pairs[j] = i
    return pairs


def _is_junction_bifurcated(homo_pairs, het_pairs, strand: str) -> bool:
    if not homo_pairs or not het_pairs:
        return False
    if strand == 'left':
        nearest_het = max(p for p, q in het_pairs)
        return any(p > nearest_het or q > nearest_het for p, q in homo_pairs)
    else:
        nearest_het = min(p for p, q in het_pairs)
        return any(p < nearest_het or q < nearest_het for p, q in homo_pairs)


def _junction_state(pos, own_pos, other_pos, this_het, other_het, pairs, structure) -> str:
    if pos < 0 or pos >= len(structure) or structure[pos] == '&':
        return 'U'
    c = structure[pos]
    if c in '()':
        if pos in pairs:
            partner = pairs[pos]
            if partner in own_pos:
                return 'PH'
            if partner in other_pos:
                return 'PHe'
        return 'PH'
    # c == '.'
    if this_het:
        if other_het:
            return 'UB'
        het_left  = any(p < pos and p in pairs and pairs[p] in other_pos for p in own_pos)
        het_right = any(p > pos and p in pairs and pairs[p] in other_pos for p in own_pos)
        if het_left and het_right:
            return 'UB'
    return 'U'


def _classify_single(structure: str) -> dict:
    structure = structure.strip()
    if '&' not in structure:
        raise ValueError(f"Structure missing '&' junction marker: {structure}")

    amp   = structure.index('&')
    pairs = _get_bracket_pairs(structure)

    left_pos  = set(range(0, amp))
    right_pos = set(range(amp + 1, len(structure)))

    lhet_pairs  = [(p, pairs[p]) for p in left_pos  if p in pairs and pairs[p] in right_pos]
    rhet_pairs  = [(p, pairs[p]) for p in right_pos if p in pairs and pairs[p] in left_pos]
    lhomo_pairs = [(p, pairs[p]) for p in left_pos  if p in pairs and pairs[p] in left_pos  and p < pairs[p]]
    rhomo_pairs = [(p, pairs[p]) for p in right_pos if p in pairs and pairs[p] in right_pos and p < pairs[p]]

    left_het   = bool(lhet_pairs)
    right_het  = bool(rhet_pairs)
    left_homo  = _is_junction_bifurcated(lhomo_pairs, lhet_pairs, 'left')  if left_het  else bool(lhomo_pairs)
    right_homo = _is_junction_bifurcated(rhomo_pairs, rhet_pairs, 'right') if right_het else bool(rhomo_pairs)

    lj = _junction_state(amp - 1, left_pos,  right_pos, left_het,  right_het, pairs, structure)
    rj = _junction_state(amp + 1, right_pos, left_pos,  right_het, left_het,  pairs, structure)

    key          = (left_homo, left_het, right_homo, right_het, lj, rj)
    cofold_class = CLASS_LOOKUP.get(key, 0)

    state_map = {'PH': 'paired_homodimer', 'PHe': 'paired_heterodimer',
                 'UB': 'unpaired_bulge',   'U':   'unpaired'}

    return {
        'left_homodimer':    left_homo,
        'left_heterodimer':  left_het,
        'right_homodimer':   right_homo,
        'right_heterodimer': right_het,
        'left_junction':     state_map[lj],
        'right_junction':    state_map[rj],
        'cofold_class':      cofold_class,
    }

# Default flanking sequences for the original FOUR_N / DEFINED protocol.
# FLUSH_NORMAL layouts (and any other layout with different library chemistry)
# can override these per-layout via adaptor_5p_ext / adaptor_3p_tail fields
# in the read_layouts section of adaptor_schemas.yaml.
DEFAULT_ADAPTOR_5P_EXT  = "GUUCAGAGUUCUACAGUCCGACGAUC"   # 5' library adaptor (RNA)
DEFAULT_ADAPTOR_3P_TAIL = "AAUGAATC"                      # 3' end tail (RNA)

DEFAULT_COFOLD_TEMPERATURE = 16.0   # °C — ligation reactions run at 16 °C (override per run/sample)

# Right-strand column order for each 3' adaptor type (all DNA except VCG)
RIGHT_COLS = {
    "NORMAL": [],
    "FOUR_N": ["N4"],
    "MADAP":  ["TRR", "BR", "MAC", "IRR", "SP"],
}


# ── Schema loading ────────────────────────────────────────────────────────────

def load_layout(schemas_file: str, read_layout: str) -> tuple[dict, dict, dict, str, str]:
    """Returns (layout_cfg, five_cfg, three_cfg, adaptor_5p_ext, adaptor_3p_tail)."""
    with open(schemas_file) as fh:
        s = yaml.safe_load(fh)
    layouts = s["read_layouts"]
    layout    = layouts[read_layout]
    five_cfg  = s["five_prime_adaptors"][layout["five_prime"]]
    three_cfg = s["three_prime_adaptors"][layout["three_prime"]]
    adaptor_5p_ext  = layout.get("adaptor_5p_ext",  DEFAULT_ADAPTOR_5P_EXT)
    adaptor_3p_tail = layout.get("adaptor_3p_tail", DEFAULT_ADAPTOR_3P_TAIL)
    return layout, five_cfg, three_cfg, adaptor_5p_ext, adaptor_3p_tail


# ── Strand reconstruction ─────────────────────────────────────────────────────

def build_strands(
    row: pd.Series,
    five_type: str,
    three_type: str,
    adaptor_5p_ext: str,
    adaptor_3p_tail: str,
) -> tuple[str, str]:
    """
    Build (seq_l, seq_r) for RNA.cofold.

    seq_l = adaptor_5p_ext + 5' tag/fixed bases
        DEFINED  → ext + "TCTA"
        FLUSH    → ext only (no per-read tag; ligation junction is immediately after ext)
        FOUR_N / EIGHT_N → ext + 5_Adap bases

    seq_r = VCG(RNA) + [intermediate structure cols as DNA] + 3_Adap + adaptor_3p_tail
    """
    if five_type == "DEFINED":
        seq_l = adaptor_5p_ext + "TCTA"
    elif five_type == "FLUSH":
        seq_l = adaptor_5p_ext
    else:  # FOUR_N, EIGHT_N
        seq_l = adaptor_5p_ext + str(row.get("5_Adap", ""))

    vcg_rna = str(row["VCG"]).replace("T", "U").replace("t", "u")

    middle_cols = RIGHT_COLS.get(three_type, [])
    middle = "".join(str(row[c]) for c in middle_cols if c in row.index)

    seq_r = vcg_rna + middle + str(row.get("3_Adap", "")) + adaptor_3p_tail

    return seq_l, seq_r


# ── Junction loop size ────────────────────────────────────────────────────────

def junction_loop_size(db: str, junction_index: int) -> int:
    """
    Count unpaired bases at the ligation junction in the dot-bracket string.
    junction_index = len(seq_l), i.e. the position of the '&' separator.
    """
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


# ── Per-read cofold ───────────────────────────────────────────────────────────

def _cofold_strands(seq_l: str, seq_r: str) -> dict:
    """Run RNA.cofold on a pre-built strand pair and return result dict."""
    combined = seq_l + "&" + seq_r
    dot_bracket, mfe = RNA.cofold(combined)
    if "&" not in dot_bracket:
        dot_bracket = dot_bracket[:len(seq_l)] + "&" + dot_bracket[len(seq_l):]
    amp_pos   = dot_bracket.index("&")
    loop_size = junction_loop_size(dot_bracket, amp_pos)
    cls = _classify_single(dot_bracket)
    return {
        "seq_l":              seq_l,
        "seq_r":              seq_r,
        "dotbracket":         dot_bracket,
        "mfe":                round(mfe, 4),
        "junction_loop_size": loop_size,
        **cls,
    }


def cofold_row(
    row: pd.Series,
    five_type: str,
    three_type: str,
    adaptor_5p_ext: str,
    adaptor_3p_tail: str,
) -> dict:
    seq_l, seq_r = build_strands(row, five_type, three_type, adaptor_5p_ext, adaptor_3p_tail)
    return _cofold_strands(seq_l, seq_r)


def _is_no_random(five_type: str, three_type: str) -> bool:
    """True when seq_l is constant and seq_r depends only on VCG + 3_Adap.
    Allows deduplication: fold once per unique strand pair rather than per read."""
    return five_type in {"FLUSH", "DEFINED"} and three_type == "NORMAL"




# ── Summary ───────────────────────────────────────────────────────────────────

def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    base = df.groupby(["VCG_Len", "VCG"], as_index=False).agg(
        n_reads=("mfe", "count"),
        mean_mfe=("mfe", "mean"),
        mean_loop_size=("junction_loop_size", "mean"),
    )
    base["mean_mfe"]       = base["mean_mfe"].round(4)
    base["mean_loop_size"] = base["mean_loop_size"].round(2)

    # Class frequency counts per VCG
    class_counts = (
        df.groupby(["VCG", "cofold_class"])
        .size()
        .unstack(fill_value=0)
        .add_prefix("class_")
        .reset_index()
    )
    summary = base.merge(class_counts, on="VCG", how="left")
    return summary.sort_values(["VCG_Len", "VCG"]).reset_index(drop=True)


# ── Core ──────────────────────────────────────────────────────────────────────

def secondary_structures(
    filtered_path: str,
    schemas_file: str,
    read_layout: str,
    reads_out: str,
    summary_out: str,
    sample: str = "",
    max_reads_per_vcg: int = 2000,
    progress_out=None,
    adaptor_5p_ext_override: str | None = None,
    adaptor_3p_tail_override: str | None = None,
    cofold_temperature: float = DEFAULT_COFOLD_TEMPERATURE,
) -> None:
    RNA.cvar.temperature = cofold_temperature
    layout_cfg, five_cfg, three_cfg, adaptor_5p_ext, adaptor_3p_tail = load_layout(schemas_file, read_layout)
    if adaptor_5p_ext_override is not None:
        adaptor_5p_ext = adaptor_5p_ext_override
    if adaptor_3p_tail_override is not None:
        adaptor_3p_tail = adaptor_3p_tail_override
    five_type  = layout_cfg["five_prime"]
    three_type = layout_cfg["three_prime"]

    Path(reads_out).parent.mkdir(parents=True, exist_ok=True)

    print(f"[secondary_structures:{sample}] cofold temperature: {cofold_temperature} °C", flush=True)
    print(f"[secondary_structures:{sample}] loading {filtered_path} ...", flush=True)
    df = pd.read_csv(filtered_path, sep="\t", dtype=str)
    n_total = len(df)

    # Cap reads per VCG to keep cofold tractable (only when input exceeds threshold)
    SUBSAMPLE_THRESHOLD = 500_000
    if n_total > SUBSAMPLE_THRESHOLD:
        df = (
            df.groupby("VCG", group_keys=False)
            .apply(lambda g: g.sample(min(len(g), max_reads_per_vcg), random_state=42))
            .reset_index(drop=True)
        )
        n = len(df)
        print(
            f"[secondary_structures:{sample}] Layout: {read_layout}  |  "
            f"Reads: {n:,} (subsampled to {max_reads_per_vcg}/VCG; input was {n_total:,})"
        )
    else:
        n = n_total
        print(f"[secondary_structures:{sample}] Layout: {read_layout}  |  Reads: {n:,}")

    if _is_no_random(five_type, three_type):
        # seq_l is constant; seq_r depends only on VCG + 3_Adap.
        # Fold each unique strand pair once and broadcast to all reads.
        seq_ls, seq_rs = zip(*[
            build_strands(row, five_type, three_type, adaptor_5p_ext, adaptor_3p_tail)
            for _, row in df.iterrows()
        ])
        df["_seq_l"] = seq_ls
        df["_seq_r"] = seq_rs
        unique_pairs = df[["_seq_l", "_seq_r"]].drop_duplicates().reset_index(drop=True)
        n_unique = len(unique_pairs)
        msg = f"[secondary_structures:{sample}] cofolding {n_unique:,} unique pairs (from {n:,} reads)...\n"
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
        if progress_out is not None:
            progress_out.write(f"[secondary_structures:{sample}] cofolding {n:,} reads...\n")
            progress_out.flush()

        result_rows = []
        report_every = max(1, n // 20)
        for i, (_, row) in enumerate(df.iterrows()):
            result_rows.append(cofold_row(row, five_type, three_type, adaptor_5p_ext, adaptor_3p_tail))
            if (i + 1) % report_every == 0:
                _redraw_progress_table(progress_out, sample, i + 1, n, "reads")
        _redraw_progress_table(progress_out, sample, n, n, "reads")

    _finish_progress_entry(progress_out, sample)
    if progress_out is not None:
        progress_out.write(f"[secondary_structures:{sample}] cofold complete\n")
        progress_out.flush()

    result_df = pd.concat([df.reset_index(drop=True), pd.DataFrame(result_rows)], axis=1)
    result_df["mfe"]                = pd.to_numeric(result_df["mfe"])
    result_df["junction_loop_size"] = pd.to_numeric(result_df["junction_loop_size"])
    result_df["VCG_Len"]            = pd.to_numeric(result_df["VCG_Len"])

    # Class count summary
    tag = f"[secondary_structures:{sample}]"
    class_counts = result_df["cofold_class"].value_counts().sort_index()
    print(f"{tag} Cofold class distribution ({n:,} reads):")
    for cls, count in class_counts.items():
        desc = CLASS_DESCRIPTIONS.get(cls, "unknown")
        pct  = f"{count / n:.1%}"
        print(f"{tag}   class {cls:>2d}  {count:>8,}  ({pct})  {desc}")

    result_df.to_csv(reads_out, sep="\t", index=False)
    print(f"{tag} Reads → {reads_out}")

    summary = build_summary(result_df)
    summary.to_csv(summary_out, sep="\t", index=False)
    print(f"{tag} Summary → {summary_out}")


# ── Entry points ──────────────────────────────────────────────────────────────

if "snakemake" in dir():
    _term = sys.stderr  # save before redirect — writes go to Snakemake terminal
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    secondary_structures(
        filtered_path         = snakemake.input.filtered,
        schemas_file          = snakemake.input.schemas,
        read_layout           = snakemake.params.read_layout,
        reads_out             = snakemake.output.reads,
        summary_out           = snakemake.output.summary,
        sample                = snakemake.wildcards.sample,
        max_reads_per_vcg     = snakemake.params.max_reads_per_vcg,
        progress_out          = _term,
        adaptor_5p_ext_override  = snakemake.params.adaptor_5p_ext,
        adaptor_3p_tail_override = snakemake.params.adaptor_3p_tail,
        cofold_temperature       = float(snakemake.params.cofold_temperature),
    )

elif __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Per-read RNA cofold ligation junction analysis.")
    ap.add_argument("--filtered",           required=True)
    ap.add_argument("--schemas",            default="adaptor_schemas.yaml")
    ap.add_argument("--read-layout",        required=True)
    ap.add_argument("--reads-out",          required=True)
    ap.add_argument("--summary-out",        required=True)
    ap.add_argument("--max-reads-per-vcg",  type=int, default=2000)  # matches snakefile default
    ap.add_argument("--cofold-temperature", type=float, default=DEFAULT_COFOLD_TEMPERATURE)
    args = ap.parse_args()
    secondary_structures(
        args.filtered, args.schemas, args.read_layout,
        args.reads_out, args.summary_out,
        max_reads_per_vcg=args.max_reads_per_vcg,
        cofold_temperature=args.cofold_temperature,
    )
