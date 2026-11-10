#!/usr/bin/env python3
"""
parse_reads.py

Trim adaptors, slice reads into structural regions, and apply QC filters.
Behaviour is fully controlled by a named read_layout in adaptor_schemas.yaml.

Processing is two-stage:
  1. Strip the 5' adaptor (fixed sequence or random N-mer tag).
  2. Find and parse the 3' adaptor structure on the remainder.

Reads are processed as DNA throughout. VCG→RNA conversion happens only at
the normalisation step. The FASTQ read ID is preserved as 'read_id' in all
outputs until UMI collapse.

The full raw read is preserved as 'Seq' and the matched 3' adaptor sequence
as '3_Adap' so the entire read can be reconstructed from the output columns.

5' random tag bases (FOUR_N / EIGHT_N) are always recorded as the '5_Adap'
column for tracking, but are not used for UMI deduplication.

Outputs:
    parsed_reads.tsv  — successfully parsed reads with all region columns
    bad_reads.tsv     — reads dropped at any QC step, with reason column

Called by Snakemake (rule parse_reads) or directly:
    python scripts/parse_reads.py \\
        --input high_quality.fastq.gz \\
        --parsed parsed_reads.tsv \\
        --dropped bad_reads.tsv \\
        --read-layout FOUR_N_MADAP \\
        --adaptor-3p AGATCGGAAGAGCACACGTCT \\
        --schemas adaptor_schemas.yaml
"""

import argparse
import csv
import gzip
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Optional

import pandas as pd
import yaml
from Bio import SeqIO
from cutadapt import adapters


# ── Schema loading ────────────────────────────────────────────────────────────

def load_layout(schemas_file: str, read_layout: str) -> tuple[dict, dict, dict]:
    """Returns (layout_cfg, five_prime_cfg, three_prime_cfg) for the named layout."""
    with open(schemas_file) as fh:
        s = yaml.safe_load(fh)

    layouts = s.get("read_layouts", {})
    if read_layout not in layouts:
        raise ValueError(
            f"Unknown read_layout '{read_layout}'. "
            f"Available: {list(layouts)}"
        )

    layout     = layouts[read_layout]
    five_cfg   = s["five_prime_adaptors"][layout["five_prime"]]
    three_cfg  = s["three_prime_adaptors"][layout["three_prime"]]
    return layout, five_cfg, three_cfg


# ── Cutadapt adaptor objects ──────────────────────────────────────────────────

def make_3p_adaptor(seq: str) -> adapters.FrontAdapter:
    return adapters.FrontAdapter(
        sequence=seq, max_errors=0, adapter_wildcards=True, indels=False
    )

def make_shift_templates(landmarks: list[str]) -> list[adapters.AnywhereAdapter]:
    return [
        adapters.AnywhereAdapter(
            sequence=lm,
            max_errors=1,
            adapter_wildcards=False,
            indels=True,
            min_overlap=len(lm),
        )
        for lm in landmarks
    ]


# ── Progress bar helper ───────────────────────────────────────────────────────

def approx_num_reads(path: str, sample_lines: int = 40) -> int:
    opener = gzip.open if str(path).endswith(".gz") else open
    try:
        with opener(path, "rt") as fh:
            lengths = [len(fh.readline()) for _ in range(sample_lines)]
        avg = sum(lengths) / len(lengths)
        return max(1, round(os.path.getsize(path) / avg / 4))
    except Exception:
        return 0


# ── Stage 1: 5' adaptor stripping ────────────────────────────────────────────

def strip_five_prime(
    seq: str, read_id: str, five_cfg: dict
) -> tuple[Optional[str], Optional[str], Optional[dict]]:
    """
    Strip the 5' adaptor from seq.
    Returns (five_bases, remainder, bad_row). bad_row is non-None if dropped.
    """
    length = five_cfg["length"]

    if len(seq) <= length:
        return None, None, {"read_id": read_id, "reason": "read_too_short"}

    if five_cfg["type"] == "fixed":
        expected = five_cfg["sequence"]
        if not seq.startswith(expected):
            return None, None, {
                "read_id": read_id,
                "reason": f"bad_5p_adaptor:{seq[:length]}",
            }

    five_bases = seq[:length]
    remainder  = seq[length:]
    return five_bases, remainder, None


# ── Stage 2: 3' adaptor parsing ───────────────────────────────────────────────

def parse_3p_normal(
    seq: str, read_id: str, adaptor_obj: adapters.FrontAdapter
) -> tuple[Optional[dict], Optional[dict]]:
    """VCG is everything left of the 3' adaptor."""
    match = adaptor_obj.match_to(seq)
    if not match:
        return None, {"read_id": read_id, "reason": "no_adaptor"}

    vcg     = seq[:match.rstart]
    adap_3p = seq[match.rstart:match.rstart + 21]

    if "-" in vcg or "_" in vcg:
        return None, {"read_id": read_id, "reason": "delimiter_in_vcg"}

    return {
        "VCG": vcg, "VCG_Len": len(vcg),
        "3_Adap": adap_3p, "Errors": match.errors,
    }, None


def parse_3p_pretrimmed(
    seq: str, read_id: str
) -> tuple[Optional[dict], Optional[dict]]:
    """Reads are already trimmed to the insert (no 3' adaptor to find);
    the whole read IS the VCG. Used by pre-trimmed layouts (e.g. SUPPLIER_T_POLYA)."""
    if not seq:
        return None, {"read_id": read_id, "reason": "empty_insert"}
    if "-" in seq or "_" in seq:
        return None, {"read_id": read_id, "reason": "delimiter_in_vcg"}
    return {
        "VCG": seq, "VCG_Len": len(seq),
        "3_Adap": "", "Errors": 0,
    }, None


def parse_3p_four_n(
    seq: str, read_id: str, adaptor_obj: adapters.FrontAdapter
) -> tuple[Optional[dict], Optional[dict]]:
    """Last 4 bp before adaptor = N4 tag; everything before = VCG."""
    match = adaptor_obj.match_to(seq)
    if not match:
        return None, {"read_id": read_id, "reason": "no_adaptor"}

    pre     = seq[:match.rstart]
    adap_3p = seq[match.rstart:match.rstart + 21]

    if len(pre) < 4:
        return None, {"read_id": read_id, "reason": "read_too_short_for_n4"}

    n4  = pre[-4:]
    vcg = pre[:-4]

    if "-" in vcg or "_" in vcg or "-" in n4 or "_" in n4:
        return None, {"read_id": read_id, "reason": "delimiter_in_region"}

    return {
        "VCG": vcg, "VCG_Len": len(vcg),
        "N4": n4, "3_Adap": adap_3p, "Errors": match.errors,
    }, None


def _find_madap_shift(
    pre_adaptor: str,
    templates: list,
    expected_start: int = 4,
) -> int | str:
    """
    Search for the BR+MADAP landmark in the last 19 bases of pre_adaptor.
    Returns shift offset or 'not found'.

    expected_start=4: in a correctly formed read the landmark begins 4 bases
    from the start of the 19-base window (after the TRR region).
    """
    window = pre_adaptor[-19:]
    for t in templates:
        m = t.match_to(window)
        if m:
            return m.rstart - expected_start
    return "not found"


def parse_3p_madap(
    seq: str,
    read_id: str,
    adaptor_obj: adapters.FrontAdapter,
    shift_templates: list,
    fixed_regions: dict,
) -> tuple[Optional[dict], Optional[dict]]:
    """
    Parse the MADAP 3' structure from seq (5' adaptor already stripped).

    Structure of seq:
        [VCG: variable][TRR: 4bp][BR: AGATCG][MAC: CAA|GTT]
        [IRR: 4bp][SP: CT] | 3' adaptor
    """
    match = adaptor_obj.match_to(seq)
    if not match:
        return None, {"read_id": read_id, "reason": "no_adaptor"}

    errors      = match.errors
    pre_adaptor = seq[:match.rstart]
    adap_3p     = seq[match.rstart:match.rstart + 21]

    shift = _find_madap_shift(pre_adaptor, shift_templates, expected_start=4)
    if shift == "not found":
        return None, {"read_id": read_id, "reason": "no_madap_landmark"}

    spacer     = pre_adaptor[-2:]
    rand_4     = pre_adaptor[-6:-2]
    madap_code = pre_adaptor[-9:-6]
    bind_6     = pre_adaptor[-15:-9]
    n4         = pre_adaptor[-19:-15]
    vcg        = pre_adaptor[:(-19 + shift)]

    for name, val in [("TRR", n4), ("BR", bind_6), ("MAC", madap_code),
                      ("IRR", rand_4), ("VCG", vcg)]:
        if "-" in val or "_" in val:
            return None, {"read_id": read_id, "reason": f"delimiter_in_{name}"}

    exp_bind6  = fixed_regions.get("BR",     {}).get("dna", "AGATCG")
    exp_madaps = fixed_regions.get("MAC", {}).get("dna", ["CAA", "GTT"])
    exp_spacer = fixed_regions.get("SP",     {}).get("dna", "CT")

    if bind_6 != exp_bind6:
        return None, {"read_id": read_id, "reason": f"bad_bind6:{bind_6}"}
    if madap_code not in exp_madaps:
        return None, {"read_id": read_id, "reason": f"bad_madap:{madap_code}"}
    if spacer != exp_spacer:
        return None, {"read_id": read_id, "reason": f"bad_spacer:{spacer}"}

    return {
        "VCG": vcg, "VCG_Len": len(vcg),
        "TRR": n4, "BR": bind_6, "MAC": madap_code,
        "IRR": rand_4, "SP": spacer,
        "3_Adap": adap_3p,
        "Errors": errors, "Shift": shift,
    }, None


# ── Per-record dispatcher ─────────────────────────────────────────────────────

def parse_record(
    record,
    five_cfg: dict,
    three_cfg: dict,
    three_prime_type: str,
    adaptor_obj: adapters.FrontAdapter,
    shift_templates: list,
) -> tuple[Optional[dict], Optional[dict]]:
    """Full parse of one FASTQ record. Returns (row_dict, bad_dict)."""
    read_id = record.id
    seq     = str(record.seq)

    # Stage 1: strip 5' adaptor
    five_bases, remainder, bad = strip_five_prime(seq, read_id, five_cfg)
    if bad:
        return None, bad

    # Stage 2: parse 3' adaptor structure
    if three_prime_type == "NORMAL":
        row, bad = parse_3p_normal(remainder, read_id, adaptor_obj)
    elif three_prime_type == "FOUR_N":
        row, bad = parse_3p_four_n(remainder, read_id, adaptor_obj)
    elif three_prime_type == "MADAP":
        row, bad = parse_3p_madap(
            remainder, read_id, adaptor_obj,
            shift_templates, three_cfg.get("fixed_regions", {}),
        )
    elif three_prime_type == "PRETRIMMED":
        row, bad = parse_3p_pretrimmed(remainder, read_id)
    else:
        raise ValueError(f"Unknown three_prime adaptor type: {three_prime_type}")

    if bad:
        return None, bad

    # Record 5' tag bases if the adaptor has a named column
    if five_cfg.get("column_name"):
        row[five_cfg["column_name"]] = five_bases

    # Prepend read_id and full raw read so they appear as first columns
    row = {"read_id": read_id, "Seq": seq, **row}
    return row, None


# ── Main parsing loop ─────────────────────────────────────────────────────────

def parse_reads(
    input_path: str,
    parsed_path: str,
    dropped_path: str,
    read_layout: str,
    adaptor_3p: str,
    schemas_file: str,
    sample: str = "",
) -> None:
    layout_cfg, five_cfg, three_cfg = load_layout(schemas_file, read_layout)
    three_prime_type = layout_cfg["three_prime"]

    adaptor_obj     = None if three_prime_type == "PRETRIMMED" else make_3p_adaptor(adaptor_3p)
    shift_templates = []
    if three_prime_type == "MADAP":
        landmarks       = three_cfg.get("landmarks_dna", ["AGATCGGTT", "AGATCGCAA"])
        shift_templates = make_shift_templates(landmarks)

    opener = gzip.open if str(input_path).endswith(".gz") else open
    Path(parsed_path).parent.mkdir(parents=True, exist_ok=True)

    all_parsed_rows = []
    dropped_rows    = []
    chunk_size      = 1000
    dropped_writer  = None
    dropped_header_written = False
    n_total         = 0
    reason_counts: Counter = Counter()

    with (
        open(dropped_path, "w", newline="") as fh_d,
        opener(input_path, "rt") as fh_in,
    ):
        for record in SeqIO.parse(fh_in, "fastq"):
            n_total += 1
            row, bad = parse_record(
                record, five_cfg, three_cfg, three_prime_type,
                adaptor_obj, shift_templates,
            )

            if row:
                all_parsed_rows.append(row)
            if bad:
                # Normalise reasons that embed a sequence (e.g. bad_bind6:AGATCC → bad_bind6)
                reason = bad["reason"].split(":")[0]
                reason_counts[reason] += 1
                dropped_rows.append(bad)

            if len(dropped_rows) >= chunk_size:
                if not dropped_header_written:
                    dropped_writer = csv.DictWriter(
                        fh_d, fieldnames=list(dropped_rows[0].keys()), delimiter="\t"
                    )
                    dropped_writer.writeheader()
                    dropped_header_written = True
                dropped_writer.writerows(dropped_rows)
                dropped_rows = []

        if dropped_rows:
            if not dropped_header_written:
                dropped_writer = csv.DictWriter(
                    fh_d, fieldnames=list(dropped_rows[0].keys()), delimiter="\t"
                )
                dropped_writer.writeheader()
            dropped_writer.writerows(dropped_rows)

        if not dropped_header_written:
            fh_d.write("read_id\treason\n")

    # Build df_4N and write TSV
    if all_parsed_rows:
        df_4N = pd.DataFrame(all_parsed_rows).convert_dtypes()
    else:
        df_4N = pd.DataFrame(columns=["read_id", "Seq", "VCG", "VCG_Len", "Errors"])

    df_4N.to_csv(parsed_path, sep="\t", index=False)

    n_parsed  = len(df_4N)
    n_dropped = n_total - n_parsed
    pct_pass  = f"{n_parsed  / n_total:.2%}" if n_total else "N/A"
    pct_drop  = f"{n_dropped / n_total:.2%}" if n_total else "N/A"

    tag = f"[parse_reads:{sample}]"
    print(f"{tag} Layout:           {read_layout}")
    print(f"{tag} Total reads:      {n_total:,}")
    print(f"{tag} Parsed (pass):    {n_parsed:,}  ({pct_pass})")
    print(f"{tag} Dropped (fail):   {n_dropped:,}  ({pct_drop})")

    if reason_counts:
        print(f"{tag} Drop reason breakdown:")
        for reason, count in reason_counts.most_common():
            pct = f"{count / n_total:.2%}"
            print(f"{tag}   {reason:<35s} {count:>8,}  ({pct})")

    # For MADAP layouts: show how many reads have each MADAP code (CAA/GTT)
    if three_prime_type == "MADAP" and "MAC" in df_4N.columns:
        code_counts = df_4N["MAC"].value_counts()
        print(f"{tag} MAC distribution (parsed reads):")
        for code, count in code_counts.items():
            pct = f"{count / n_parsed:.2%}"
            print(f"{tag}   {code}  {count:>8,}  ({pct})")

    print(f"{tag} Parsed  → {parsed_path}")
    print(f"{tag} Dropped → {dropped_path}")


# ── Entry points ──────────────────────────────────────────────────────────────

if "snakemake" in dir():
    sys.stdout = sys.stderr = open(snakemake.log[0], "w", buffering=1, encoding="utf-8")
    parse_reads(
        input_path      = snakemake.input.fastq,
        parsed_path     = snakemake.output.parsed,
        dropped_path    = snakemake.output.dropped,
        read_layout     = snakemake.params.read_layout,
        adaptor_3p      = snakemake.params.adaptor_3p,
        schemas_file    = snakemake.input.schemas,
        sample          = snakemake.wildcards.sample,
    )

elif __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Two-stage adaptor trim and region slice.")
    ap.add_argument("--input",        required=True)
    ap.add_argument("--parsed",       required=True)
    ap.add_argument("--dropped",      required=True)
    ap.add_argument("--read-layout",  required=True)
    ap.add_argument("--adaptor-3p",   required=True)
    ap.add_argument("--schemas",      default="adaptor_schemas.yaml")
    args = ap.parse_args()
    parse_reads(
        args.input, args.parsed, args.dropped,
        args.read_layout, args.adaptor_3p, args.schemas,
    )
