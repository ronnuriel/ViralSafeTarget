#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PYTHON_BIN="${PYTHON_BIN:-python}"
DATASETS_BIN="${DATASETS_BIN:-datasets}"
SUMMARY="data/raw/hsv2_population/summary.jsonl"
DOWNLOAD_DIR="reports/hsv2_population_download"
HELDOUT_DIR="reports/hsv2_population_heldout"
ZIP="data/raw/hsv2_population/ncbi_dataset.zip"
DOWNLOADED_FASTA="data/raw/hsv2_population/population_download.fasta"
CANDIDATES="reports/hsv2_genome_wide/genome_wide_candidates_post_human.csv"
REFERENCE="data/raw/hsv2_reference/ncbi_dataset/data/genomic.fna"
DISCOVERY_FASTA="data/processed/hsv2_aligned_25.fasta"
REPORT_DIR="reports/hsv2_population_report_balanced"

mkdir -p "$(dirname "$SUMMARY")" "$DOWNLOAD_DIR"

if [[ ! -s "$SUMMARY" ]]; then
  "$DATASETS_BIN" summary virus genome taxon 10310 --as-json-lines > "$SUMMARY"
fi

"$PYTHON_BIN" scripts/prepare_hsv2_population_panel.py \
  --summary-jsonl "$SUMMARY" \
  --out-dir "$DOWNLOAD_DIR"

if [[ ! -s "$ZIP" ]]; then
  "$DATASETS_BIN" download virus genome accession \
    --inputfile "$DOWNLOAD_DIR/accessions.txt" \
    --include genome \
    --filename "$ZIP" \
    --no-progressbar
fi

if [[ ! -s "$DOWNLOADED_FASTA" ]]; then
  unzip -p "$ZIP" 'ncbi_dataset/data/genomic.fna' > "$DOWNLOADED_FASTA"
fi

"$PYTHON_BIN" scripts/prepare_hsv2_population_panel.py \
  --summary-jsonl "$SUMMARY" \
  --source-fasta "$DOWNLOADED_FASTA" \
  --exclude-fasta "$DISCOVERY_FASTA" \
  --out-dir "$HELDOUT_DIR"

"$PYTHON_BIN" -m viral_safe_target analyze population \
  --population-fasta "$HELDOUT_DIR/population_unique.fasta" \
  --reference-fasta "$REFERENCE" \
  --candidates "$CANDIDATES" \
  --out-dir "$REPORT_DIR"

echo "Population report: $ROOT/$REPORT_DIR/population_validation_report.html"
