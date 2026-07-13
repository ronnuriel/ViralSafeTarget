#!/usr/bin/env bash
set -euo pipefail

INPUT=${1:-data/processed/virus_sample.fasta}
OUTPUT=${2:-data/processed/virus_aligned.fasta}

mkdir -p "$(dirname "$OUTPUT")"
mafft --auto "$INPUT" > "$OUTPUT"
echo "Alignment written to $OUTPUT"
