#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python}"

if command -v vst >/dev/null 2>&1; then
  VST=(vst)
elif command -v viral-safe-target >/dev/null 2>&1; then
  VST=(viral-safe-target)
else
  VST=("$PYTHON_BIN" -m viral_safe_target)
fi

"${VST[@]}" scan \
  --virus-alignment data/demo/virus_aligned.fasta \
  --reference-id HSV2_demo_ref \
  --gff data/demo/reference.gff3 \
  --small-host-fasta data/demo/human_mini.fasta \
  --out-dir reports/demo \
  --min-coverage 0.80

"${VST[@]}" simulate-pairs \
  --candidates reports/demo/candidates.csv \
  --gff data/demo/reference.gff3 \
  --virus-alignment data/demo/virus_aligned.fasta \
  --reference-id HSV2_demo_ref \
  --output reports/demo/simulated_pairs.csv \
  --max-candidates 100

echo "Demo completed: reports/demo/"
