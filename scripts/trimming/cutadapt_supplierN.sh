#!/bin/bash
# Usage: ./trim_reads.sh input.fastq.gz output.fastq
# Example: ./trim_reads.sh sample.fastq.gz sample.trimmed.fastq

# Exit on error or undefined variable
set -euo pipefail

# Check argument count
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <input_R1.fastq.gz> <input_R2.fastq.gz> <output_R1.fastq> <output_R2.fastq>" >&2
  exit 1
fi

# Assign arguments
INPUT_FILE_R1="$1"
INPUT_FILE_R2="$2"
OUTPUT_FILE_R1="$3"
OUTPUT_FILE_R2="$4"

# Run cutadapt
cutadapt \
  --cores 0 \
  --overlap 8 \
  --error-rate 0 \
  --quality-cutoff 30 \
  --discard-untrimmed --pair-filter=any \
  -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -o "$OUTPUT_FILE_R1" \
  -p "$OUTPUT_FILE_R2" \
  "$INPUT_FILE_R1" \
  "$INPUT_FILE_R2"

echo "Trimming complete. Output saved to:"
echo "  R1: $OUTPUT_FILE_R1"
echo "  R2: $OUTPUT_FILE_R2"
