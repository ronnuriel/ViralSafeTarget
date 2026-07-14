#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

PYTHON_BIN="${PYTHON_BIN:-/Users/ronnuriel/mambaforge/envs/vst-runtime/bin/python}"
EXHAUSTIVE_PID="${EXHAUSTIVE_PID:-}"
CAS_OFFINDER_PATH="${CAS_OFFINDER_PATH:-/tmp/cas-offinder-vst-build/build/cas-offinder}"
LOG="reports/hsv2_genome_wide_exhaustive/finalization.log"

if [[ -n "$EXHAUSTIVE_PID" ]]; then
  while kill -0 "$EXHAUSTIVE_PID" 2>/dev/null; do
    sleep 30
  done
fi

exec > >(tee -a "$LOG") 2>&1
echo "Finalization started: $(date -u +%Y-%m-%dT%H:%M:%SZ)"

if [[ ! -x "$CAS_OFFINDER_PATH" ]]; then
  echo "Cas-OFFinder executable is missing or not executable: $CAS_OFFINDER_PATH" >&2
  exit 1
fi
export CAS_OFFINDER_PATH

all_batches_completed() {
  "$PYTHON_BIN" - <<'PY'
import json
from pathlib import Path

batch_root = Path("reports/hsv2_genome_wide_exhaustive/batches")
manifests = sorted(batch_root.glob("batch_*/manifest.json"))
incomplete = []
for path in manifests:
    record = json.loads(path.read_text(encoding="utf-8"))
    if record.get("status") != "completed":
        incomplete.append(path.parent.name)
print(f"Batch manifests: {len(manifests)}; incomplete: {len(incomplete)}")
if incomplete:
    print("Incomplete batches:", ", ".join(incomplete[:20]))
raise SystemExit(0 if manifests and not incomplete else 1)
PY
}

for attempt in 1 2 3; do
  echo "Exhaustive completion/resume attempt $attempt/3"
  "$PYTHON_BIN" -m viral_safe_target discover genome-wide \
    --virus hsv2 \
    --config configs/hsv2_genome_wide.yaml \
    --out-dir reports/hsv2_genome_wide_exhaustive \
    --exhaustive \
    --confirm-exhaustive \
    --run-cas-offinder
  if all_batches_completed; then
    break
  fi
done

if ! all_batches_completed; then
  echo "Exhaustive screening remains incomplete after three resume attempts." >&2
  exit 1
fi

"$PYTHON_BIN" -m viral_safe_target analyze population \
  --population-fasta reports/hsv2_population_heldout/population_unique.fasta \
  --reference-fasta data/raw/hsv2_reference/ncbi_dataset/data/genomic.fna \
  --candidates reports/hsv2_genome_wide_exhaustive/genome_wide_candidates_post_human.csv \
  --out-dir reports/hsv2_population_report_exhaustive

"$PYTHON_BIN" -m viral_safe_target analyze gene-function \
  --genome-wide-dir reports/hsv2_genome_wide_exhaustive \
  --out-dir reports/hsv2_gene_function_exhaustive

"$PYTHON_BIN" scripts/compare_discovery_modes.py \
  --balanced-dir reports/hsv2_genome_wide \
  --exhaustive-dir reports/hsv2_genome_wide_exhaustive \
  --out-dir reports/hsv2_discovery_mode_comparison

"$PYTHON_BIN" -m pytest -q
"$PYTHON_BIN" -m ruff check .
git diff --check

echo "Finalization completed: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
