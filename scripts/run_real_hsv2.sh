#!/usr/bin/env bash
set -euo pipefail

# Reproducible HSV-2 real-data pilot.
# Downloads public NCBI data, selects a deterministic sample, aligns it, scans
# conserved SpCas9-compatible sites, and optionally prepares a GRCh38
# Cas-OFFinder query. It does not validate editing or viral inactivation.
#
# Usage:
#   bash scripts/run_real_hsv2.sh
#   bash scripts/run_real_hsv2.sh --sample-size 50
#   bash scripts/run_real_hsv2.sh --with-human --sample-size 25

WITH_HUMAN=0
SAMPLE_SIZE=25
MIN_COVERAGE=0.95

while [[ $# -gt 0 ]]; do
  case "$1" in
    --with-human)
      WITH_HUMAN=1
      shift
      ;;
    --sample-size)
      SAMPLE_SIZE="$2"
      shift 2
      ;;
    --min-coverage)
      MIN_COVERAGE="$2"
      shift 2
      ;;
    -h|--help)
      sed -n '1,24p' "$0"
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

for cmd in datasets python mafft unzip curl; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Missing required command: $cmd" >&2
    echo "Create the Conda environment first: conda env create -f environment.yml" >&2
    exit 1
  fi
done

mkdir -p data/raw data/processed reports/real_hsv2

if [[ ! -f data/raw/hsv2_reference.zip ]]; then
  datasets download virus genome accession NC_001798.2 \
    --include genome,cds,protein \
    --filename data/raw/hsv2_reference.zip
fi
if [[ ! -d data/raw/hsv2_reference ]]; then
  unzip -q data/raw/hsv2_reference.zip -d data/raw/hsv2_reference
fi

if [[ ! -f data/raw/hsv2_complete.zip ]]; then
  datasets download virus genome taxon 10310 \
    --complete-only --include genome \
    --filename data/raw/hsv2_complete.zip
fi
if [[ ! -d data/raw/hsv2_complete ]]; then
  unzip -q data/raw/hsv2_complete.zip -d data/raw/hsv2_complete
fi

if [[ ! -s data/raw/hsv2_reference.gb ]]; then
  curl -fL --retry 3 \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001798.2&rettype=gbwithparts&retmode=text" \
    -o data/raw/hsv2_reference.gb
fi

REFERENCE_FASTA=$(find data/raw/hsv2_reference -type f -name 'genomic.fna' | head -n 1)
ALL_FASTA=$(find data/raw/hsv2_complete -type f -name 'genomic.fna' | head -n 1)
if [[ -z "$REFERENCE_FASTA" || -z "$ALL_FASTA" ]]; then
  echo "Could not locate genomic.fna inside the extracted NCBI packages." >&2
  exit 1
fi

python scripts/prepare_real_hsv2.py \
  --reference-fasta "$REFERENCE_FASTA" \
  --reference-genbank data/raw/hsv2_reference.gb \
  --all-genomes-fasta "$ALL_FASTA" \
  --output-fasta data/processed/hsv2_sample_${SAMPLE_SIZE}.fasta \
  --output-gff data/processed/hsv2_reference.gff3 \
  --sample-size "$SAMPLE_SIZE"

mafft --auto --thread -1 data/processed/hsv2_sample_${SAMPLE_SIZE}.fasta \
  > data/processed/hsv2_aligned_${SAMPLE_SIZE}.fasta

HUMAN_ARGS=()
if [[ "$WITH_HUMAN" -eq 1 ]]; then
  if [[ ! -f data/raw/human_GRCh38.zip ]]; then
    datasets download genome accession GCF_000001405.40 \
      --include genome,gff3 \
      --filename data/raw/human_GRCh38.zip
  fi
  if [[ ! -d data/raw/human_GRCh38 ]]; then
    unzip -q data/raw/human_GRCh38.zip -d data/raw/human_GRCh38
  fi
  HUMAN_FASTA=$(find data/raw/human_GRCh38 -type f -name '*_genomic.fna' | head -n 1)
  if [[ -z "$HUMAN_FASTA" ]]; then
    echo "Could not locate the GRCh38 genomic FASTA." >&2
    exit 1
  fi
  HUMAN_DIR=$(dirname "$HUMAN_FASTA")
  HUMAN_ARGS=(--human-fasta-directory "$HUMAN_DIR")
fi

python scripts/generate_real_candidates.py \
  --alignment data/processed/hsv2_aligned_${SAMPLE_SIZE}.fasta \
  --gff data/processed/hsv2_reference.gff3 \
  --reference-id NC_001798.2 \
  --min-coverage "$MIN_COVERAGE" \
  --out-dir reports/real_hsv2 \
  "${HUMAN_ARGS[@]}"

viral-safe-target simulate-pairs \
  --candidates reports/real_hsv2/candidates_pre_human.csv \
  --gff data/processed/hsv2_reference.gff3 \
  --virus-alignment data/processed/hsv2_aligned_${SAMPLE_SIZE}.fasta \
  --reference-id NC_001798.2 \
  --output reports/real_hsv2/simulated_pairs_pre_human.csv \
  --max-candidates 150

python - <<PY
from pathlib import Path
from viral_safe_target import write_run_manifest
write_run_manifest(
    "reports/real_hsv2/run_manifest.json",
    [
        "data/processed/hsv2_sample_${SAMPLE_SIZE}.fasta",
        "data/processed/hsv2_aligned_${SAMPLE_SIZE}.fasta",
        "data/processed/hsv2_reference.gff3",
    ],
    {
        "virus_taxon_id": 10310,
        "reference_accession": "NC_001798.2",
        "sample_size": int("${SAMPLE_SIZE}"),
        "minimum_exact_site_coverage": float("${MIN_COVERAGE}"),
        "human_screen_prepared": bool(int("${WITH_HUMAN}")),
    },
)
PY

echo
echo "Done. Main outputs:"
echo "  reports/real_hsv2/report_pre_human.html"
echo "  reports/real_hsv2/candidates_pre_human.csv"
echo "  reports/real_hsv2/simulated_pairs_pre_human.csv"
echo "  reports/real_hsv2/run_manifest.json"
if [[ "$WITH_HUMAN" -eq 1 ]]; then
  echo
echo "Next, after installing Cas-OFFinder:"
  echo "  cas-offinder reports/real_hsv2/cas_offinder_input.txt C reports/real_hsv2/cas_offinder_output.tsv"
fi
