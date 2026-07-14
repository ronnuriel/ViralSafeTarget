#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

if command -v vst >/dev/null 2>&1; then
  VST=(vst)
else
  VST=(python -m viral_safe_target)
fi

"${VST[@]}" tools doctor

MODE=(--run-cas-offinder)
for argument in "$@"; do
  if [[ "$argument" == "--analysis-only" ]]; then
    MODE=()
    break
  fi
done

"${VST[@]}" discover genome-wide \
  --virus hsv2 \
  --config configs/hsv2_genome_wide.yaml \
  --out-dir reports/hsv2_genome_wide \
  "${MODE[@]}" \
  "$@"
