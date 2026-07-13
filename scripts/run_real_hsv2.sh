#!/usr/bin/env bash
set -euo pipefail

# Resumable public-data HSV-2 workflow. Outputs are computational hypotheses,
# not claims of editing, viral inactivation, safety, or clinical efficacy.

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

WITH_HUMAN=0
SAMPLE_SIZE=25
MIN_COVERAGE=0.95
CONFIG=configs/research_v0.3.yaml

while [[ $# -gt 0 ]]; do
  case "$1" in
    --with-human) WITH_HUMAN=1; shift ;;
    --sample-size) SAMPLE_SIZE="$2"; shift 2 ;;
    --min-coverage) MIN_COVERAGE="$2"; shift 2 ;;
    --config) CONFIG="$2"; shift 2 ;;
    -h|--help)
      sed -n '1,28p' "$0"
      exit 0
      ;;
    *) echo "Unknown argument: $1" >&2; exit 2 ;;
  esac
done

if command -v vst >/dev/null 2>&1; then
  VST=(vst)
elif command -v viral-safe-target >/dev/null 2>&1; then
  VST=(viral-safe-target)
else
  VST=(python -m viral_safe_target)
fi

export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
mkdir -p data/raw data/processed reports/real_hsv2 reports/real_hsv2/.cache

declare -a STAGE_SUMMARY=()
stage_done() { STAGE_SUMMARY+=("$1: $2"); }

require_command() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Missing required command: $1" >&2
    exit 1
  fi
}

require_command python
require_command unzip
require_command curl

valid_zip() { [[ -s "$1" ]] && unzip -tq "$1" >/dev/null 2>&1; }

download_dataset() {
  local output=$1
  shift
  require_command datasets
  if valid_zip "$output"; then
    stage_done "download $(basename "$output")" "cached and validated"
    return
  fi
  rm -f "$output"
  local attempt
  for attempt in 1 2 3; do
    echo "Dataset download attempt ${attempt}/3: $output"
    if datasets download "$@" --filename "$output" && valid_zip "$output"; then
      stage_done "download $(basename "$output")" "completed and validated"
      return
    fi
    rm -f "$output"
  done
  echo "Failed to download a valid archive: $output" >&2
  exit 1
}

download_dataset data/raw/hsv2_reference.zip virus genome accession NC_001798.2 \
  --include genome,cds,protein
if [[ ! -d data/raw/hsv2_reference ]] || \
   [[ -z "$(find data/raw/hsv2_reference -type f -name 'genomic.fna' -print -quit 2>/dev/null)" ]]; then
  rm -rf data/raw/hsv2_reference
  unzip -q data/raw/hsv2_reference.zip -d data/raw/hsv2_reference
  stage_done "extract HSV-2 reference" "completed"
else
  stage_done "extract HSV-2 reference" "cached"
fi

download_dataset data/raw/hsv2_complete.zip virus genome taxon 10310 \
  --complete-only --include genome
if [[ ! -d data/raw/hsv2_complete ]] || \
   [[ -z "$(find data/raw/hsv2_complete -type f -name 'genomic.fna' -print -quit 2>/dev/null)" ]]; then
  rm -rf data/raw/hsv2_complete
  unzip -q data/raw/hsv2_complete.zip -d data/raw/hsv2_complete
  stage_done "extract HSV-2 collection" "completed"
else
  stage_done "extract HSV-2 collection" "cached"
fi

if [[ ! -s data/raw/hsv2_reference.gb ]]; then
  curl -fL --retry 5 --retry-all-errors \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001798.2&rettype=gbwithparts&retmode=text" \
    -o data/raw/hsv2_reference.gb
  stage_done "reference GenBank" "downloaded"
else
  stage_done "reference GenBank" "cached"
fi

REFERENCE_FASTA=$(find data/raw/hsv2_reference -type f -name 'genomic.fna' -print -quit)
ALL_FASTA=$(find data/raw/hsv2_complete -type f -name 'genomic.fna' -print -quit)
SAMPLE=data/processed/hsv2_sample_${SAMPLE_SIZE}.fasta
ALIGNMENT=data/processed/hsv2_aligned_${SAMPLE_SIZE}.fasta
GFF=data/processed/hsv2_reference.gff3
QC_REPORT=reports/real_hsv2/accession_qc.csv
ACCESSIONS=reports/real_hsv2/accessions_used.txt

PREPARE_STAMP=reports/real_hsv2/.cache/prepare_${SAMPLE_SIZE}.json
if python scripts/cache_stage.py check --stamp "$PREPARE_STAMP" \
  --input "$REFERENCE_FASTA" --input data/raw/hsv2_reference.gb --input "$ALL_FASTA" \
  --output "$SAMPLE" --output "$GFF" --output "$QC_REPORT" --output "$ACCESSIONS" \
  --parameter "sample_size=${SAMPLE_SIZE}"; then
  stage_done "prepare and QC" "cached"
else
  python scripts/prepare_real_hsv2.py \
    --reference-fasta "$REFERENCE_FASTA" \
    --reference-genbank data/raw/hsv2_reference.gb \
    --all-genomes-fasta "$ALL_FASTA" \
    --output-fasta "$SAMPLE" \
    --output-gff "$GFF" \
    --sample-size "$SAMPLE_SIZE" \
    --qc-report "$QC_REPORT" \
    --accessions-used "$ACCESSIONS"
  python scripts/cache_stage.py stamp --stamp "$PREPARE_STAMP" \
    --input "$REFERENCE_FASTA" --input data/raw/hsv2_reference.gb --input "$ALL_FASTA" \
    --parameter "sample_size=${SAMPLE_SIZE}"
  stage_done "prepare and QC" "completed"
fi

require_command mafft
MAFFT_STAMP=reports/real_hsv2/.cache/mafft_${SAMPLE_SIZE}.json
if python scripts/cache_stage.py check --stamp "$MAFFT_STAMP" --input "$SAMPLE" \
  --output "$ALIGNMENT" --parameter "mode=auto" --parameter "threads=-1"; then
  stage_done "MAFFT" "cached"
else
  TEMP_ALIGNMENT="${ALIGNMENT}.tmp"
  mafft --auto --thread -1 "$SAMPLE" > "$TEMP_ALIGNMENT"
  [[ -s "$TEMP_ALIGNMENT" ]] || { echo "MAFFT produced an empty alignment" >&2; exit 1; }
  mv "$TEMP_ALIGNMENT" "$ALIGNMENT"
  python scripts/cache_stage.py stamp --stamp "$MAFFT_STAMP" --input "$SAMPLE" \
    --parameter "mode=auto" --parameter "threads=-1"
  stage_done "MAFFT" "completed"
fi

CANDIDATE_STAMP=reports/real_hsv2/.cache/candidates_${SAMPLE_SIZE}.json
RANKED=reports/real_hsv2/candidates_ranked_pre_human.csv
REJECTED=reports/real_hsv2/candidates_rejected_pre_human.csv
if python scripts/cache_stage.py check --stamp "$CANDIDATE_STAMP" \
  --input "$ALIGNMENT" --input "$GFF" --input "$CONFIG" \
  --output "$RANKED" --output "$REJECTED" \
  --parameter "minimum_coverage=${MIN_COVERAGE}"; then
  stage_done "candidate scan and rank" "cached"
else
  "${VST[@]}" scan \
    --virus-alignment "$ALIGNMENT" \
    --reference-id NC_001798.2 \
    --gff "$GFF" \
    --gene-evidence data/curated/hsv2_gene_evidence.csv \
    --out-dir reports/real_hsv2 \
    --min-coverage "$MIN_COVERAGE" \
    --config "$CONFIG"
  python scripts/cache_stage.py stamp --stamp "$CANDIDATE_STAMP" \
    --input "$ALIGNMENT" --input "$GFF" --input "$CONFIG" \
    --parameter "minimum_coverage=${MIN_COVERAGE}"
  stage_done "candidate scan and rank" "completed"
fi

"${VST[@]}" simulate-pairs \
  --candidates "$RANKED" \
  --gff "$GFF" \
  --virus-alignment "$ALIGNMENT" \
  --reference-id NC_001798.2 \
  --out-dir reports/real_hsv2 \
  --config "$CONFIG"
stage_done "pair hypotheses" "completed"

"${VST[@]}" report \
  --candidates "$RANKED" \
  --rejected "$REJECTED" \
  --pairs reports/real_hsv2/pair_hypotheses_same_gene.csv \
  --multi-pairs reports/real_hsv2/pair_hypotheses_multi_target.csv \
  --out-dir reports/real_hsv2 \
  --title "HSV-2 real-data computational target-discovery report"
stage_done "research report" "completed"

if [[ "$WITH_HUMAN" -eq 1 ]]; then
  download_dataset data/raw/human_GRCh38.zip genome accession GCF_000001405.40 \
    --include genome,gff3
  if [[ ! -d data/raw/human_GRCh38 ]] || \
     [[ -z "$(find data/raw/human_GRCh38 -type f -name '*_genomic.fna' -print -quit 2>/dev/null)" ]]; then
    rm -rf data/raw/human_GRCh38
    unzip -q data/raw/human_GRCh38.zip -d data/raw/human_GRCh38
  fi
  HUMAN_FASTA=$(find data/raw/human_GRCh38 -type f -name '*_genomic.fna' -print -quit)
  "${VST[@]}" build-offtarget-input \
    --candidates "$RANKED" \
    --human-fasta-directory "$(dirname "$HUMAN_FASTA")" \
    --output reports/real_hsv2/cas_offinder_input.txt \
    --manifest reports/real_hsv2/offtarget_selected_candidates.csv \
    --config "$CONFIG"
  stage_done "human off-target input" "completed"
fi

python - "$CONFIG" "$SAMPLE_SIZE" "$MIN_COVERAGE" "$WITH_HUMAN" <<'PY'
import json
import sys
from pathlib import Path

import pandas as pd

from viral_safe_target import load_config, write_run_manifest

config_path, sample_size, minimum_coverage, with_human = sys.argv[1:]
config = load_config(config_path)
qc = pd.read_csv("reports/real_hsv2/accession_qc.csv")
accepted = qc.loc[qc["decision"] == "accepted", "accession"].astype(str).tolist()
rejected = qc.loc[qc["decision"] == "rejected", ["accession", "rejection_reason"]].to_dict("records")
outputs = [
    path for path in Path("reports/real_hsv2").glob("*")
    if path.is_file() and path.name != "run_manifest.json"
]
write_run_manifest(
    "reports/real_hsv2/run_manifest.json",
    [
        f"data/processed/hsv2_sample_{sample_size}.fasta",
        f"data/processed/hsv2_aligned_{sample_size}.fasta",
        "data/processed/hsv2_reference.gff3",
    ],
    {
        "virus_taxon_id": 10310,
        "reference_accession": "NC_001798.2",
        "sample_size": int(sample_size),
        "minimum_exact_site_coverage": float(minimum_coverage),
        "human_screen_prepared": bool(int(with_human)),
    },
    config_path=config_path,
    editor_profile=config["editor"],
    accepted_accessions=accepted,
    rejected_accessions=rejected,
    human_assembly_identifier=(
        config["off_target"]["human_assembly_accession"] if int(with_human) else None
    ),
    command_line=[
        "bash",
        "scripts/run_real_hsv2.sh",
        "--sample-size",
        sample_size,
        "--min-coverage",
        minimum_coverage,
        "--config",
        config_path,
        *(["--with-human"] if int(with_human) else []),
    ],
    random_seed=int(config["random_seed"]),
    output_paths=outputs,
)
PY
stage_done "run manifest" "completed"

echo
echo "Stage summary"
for stage in "${STAGE_SUMMARY[@]}"; do echo "- $stage"; done
echo
echo "Main outputs: reports/real_hsv2/"
echo "- candidates_ranked_pre_human.csv"
echo "- candidates_rejected_pre_human.csv"
echo "- pair_hypotheses_same_gene.csv"
echo "- pair_hypotheses_multi_target.csv"
echo "- run_manifest.json"
echo "- report.html"
