#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

if [[ -z "${CAS_OFFINDER_PATH:-}" ]]; then
  echo "Set CAS_OFFINDER_PATH to an executable Cas-OFFinder binary." >&2
  exit 2
fi

vst discover genome-wide \
  --virus hsv2 \
  --config configs/hsv2_genome_wide.yaml \
  --out-dir reports/hsv2_genome_wide_exhaustive \
  --exhaustive \
  --confirm-exhaustive \
  --run-cas-offinder
