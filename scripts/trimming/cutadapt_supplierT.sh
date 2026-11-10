#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Supplier T small RNA cutadapt wrapper.

Usage:
  cutadapt_supplierT.sh -i input.fastq -o output.fastq [options]

Options (defaults match the old one-liner):
  -i, --input FILE         Input FASTQ (required)
  -o, --output FILE        Output FASTQ (required)
  -a, --adapter SEQ        3' adapter sequence (default: AAAAAAAAAA)
  -m, --min-length N       Minimum length to keep (default: 15)
  -u, --trim-front N       Bases to trim from 5' end (default: 3)
  --cutadapt-args "ARGS"   Extra args to pass to cutadapt verbatim
  -h, --help               Show this help and exit

Example (defaults):
  cutadapt_supplierT.sh -i input.fastq -o output.fastq
USAGE
}

INPUT=""
OUTPUT=""
ADAPTER="AAAAAAAAAA"
MIN_LEN=2
TRIM_FRONT=3
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
  -i | --input)
    INPUT="${2-}"
    shift 2
    ;;
  -o | --output)
    OUTPUT="${2-}"
    shift 2
    ;;
  -a | --adapter)
    ADAPTER="${2-}"
    shift 2
    ;;
  -m | --min-length)
    MIN_LEN="${2-}"
    shift 2
    ;;
  -u | --trim-front)
    TRIM_FRONT="${2-}"
    shift 2
    ;;
  --cutadapt-args)
    EXTRA_ARGS="${2-}"
    shift 2
    ;;
  -h | --help)
    usage
    exit 0
    ;;
  *)
    echo "Unknown option: $1" >&2
    usage
    exit 1
    ;;
  esac
done

if [[ -z "$INPUT" || -z "$OUTPUT" ]]; then
  echo "Error: --input and --output are required." >&2
  usage
  exit 1
fi

if ! command -v cutadapt >/dev/null 2>&1; then
  echo "Error: cutadapt is not in PATH." >&2
  exit 1
fi

cutadapt \
  -m "$MIN_LEN" \
  -u "$TRIM_FRONT" \
  -a "$ADAPTER" \
  -o "$OUTPUT" \
  $EXTRA_ARGS \
  "$INPUT"
