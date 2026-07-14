#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python}"
"${PYTHON_BIN}" -m viral_safe_target.cli tools benchmark \
  --config configs/benchmarks/hsv2_multitool.yaml
