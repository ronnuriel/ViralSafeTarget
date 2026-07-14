#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"

if command -v vst >/dev/null 2>&1; then
  vst doctor "$@"
elif command -v viral-safe-target >/dev/null 2>&1; then
  viral-safe-target doctor "$@"
else
  PYTHONPATH=src python -m viral_safe_target doctor "$@"
fi
