#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"

CONFIG=${VST_CONFIG:-configs/hsv2_pilot.yaml}
REPORTS=reports/hsv2_pilot
REAL_REPORTS=reports/real_hsv2
mkdir -p "$REPORTS"

if command -v vst >/dev/null 2>&1; then
  VST=(vst)
elif command -v viral-safe-target >/dev/null 2>&1; then
  VST=(viral-safe-target)
else
  VST=(python -m viral_safe_target)
fi

RANKED=${VST_RANKED_CANDIDATES:-$REAL_REPORTS/candidates_ranked_pre_human.csv}
if [[ ! -s "$RANKED" ]]; then
  echo "No v0.3 ranked candidates found; resuming the real-data workflow."
  bash scripts/run_real_hsv2.sh --sample-size "${VST_SAMPLE_SIZE:-25}" --config "$CONFIG"
fi

HUMAN_DIR=${HUMAN_FASTA_DIRECTORY:-}
if [[ -z "$HUMAN_DIR" ]]; then
  HUMAN_FASTA=$(find data/raw/human_GRCh38 -type f -name '*_genomic.fna' -print -quit 2>/dev/null || true)
  [[ -n "$HUMAN_FASTA" ]] && HUMAN_DIR=$(dirname "$HUMAN_FASTA")
fi
if [[ -z "$HUMAN_DIR" ]]; then
  echo "Human FASTA directory not found. Set HUMAN_FASTA_DIRECTORY or run:" >&2
  echo "  bash scripts/run_real_hsv2.sh --with-human --sample-size ${VST_SAMPLE_SIZE:-25}" >&2
  exit 1
fi

"${VST[@]}" build-offtarget-input \
  --candidates "$RANKED" \
  --human-fasta-directory "$HUMAN_DIR" \
  --output "$REPORTS/cas_offinder_input.txt" \
  --manifest "$REPORTS/offtarget_selected_candidates.csv" \
  --genes UL19 UL30 \
  --max-candidates 200 \
  --config "$CONFIG"

ALIGNMENT=$(find data/processed -type f -name 'hsv2_aligned_*.fasta' -print | sort | tail -n 1)
"${VST[@]}" simulate-pairs \
  --candidates "$REPORTS/offtarget_selected_candidates.csv" \
  --gff data/processed/hsv2_reference.gff3 \
  --virus-alignment "$ALIGNMENT" \
  --reference-id NC_001798.2 \
  --out-dir "$REPORTS" \
  --genes UL19 UL30 \
  --maximum-candidates-per-gene 100 \
  --max-candidates 200 \
  --config "$CONFIG"

CAS_OUTPUT=${CAS_OFFINDER_OUTPUT:-$REPORTS/cas_offinder_output.tsv}
FINAL_CANDIDATES="$REPORTS/offtarget_selected_candidates.csv"
if [[ -s "$CAS_OUTPUT" ]]; then
  HUMAN_GFF=$(find data/raw/human_GRCh38 -type f -name '*genomic.gff' -print -quit 2>/dev/null || true)
  GFF_ARGS=()
  [[ -n "$HUMAN_GFF" ]] && GFF_ARGS=(--human-gff "$HUMAN_GFF")
  "${VST[@]}" summarize-offtargets \
    --candidates "$REPORTS/offtarget_selected_candidates.csv" \
    --cas-output "$CAS_OUTPUT" \
    --manifest "$REPORTS/offtarget_selected_candidates.csv" \
    --out-dir "$REPORTS" \
    --config "$CONFIG" \
    "${GFF_ARGS[@]}"
  FINAL_CANDIDATES="$REPORTS/candidates_ranked_post_human.csv"
else
  echo "Cas-OFFinder output not present; post-human ranking remains pending."
  echo "Run: cas-offinder $REPORTS/cas_offinder_input.txt C $CAS_OUTPUT"
fi

REPORT_ARGS=(
  --candidates "$FINAL_CANDIDATES"
  --rejected "$REAL_REPORTS/candidates_rejected_pre_human.csv"
  --pairs "$REPORTS/pair_hypotheses_same_gene.csv"
  --multi-pairs "$REPORTS/pair_hypotheses_multi_target.csv"
  --out-dir "$REPORTS"
  --title "HSV-2 UL19/UL30 computational pilot"
)
if [[ -s "$REPORTS/predicted_human_hits.csv" ]]; then
  "${VST[@]}" report "${REPORT_ARGS[@]}" \
    --predicted-hits "$REPORTS/predicted_human_hits.csv"
else
  "${VST[@]}" report "${REPORT_ARGS[@]}"
fi

echo "HSV-2 pilot outputs: $REPORTS/"
