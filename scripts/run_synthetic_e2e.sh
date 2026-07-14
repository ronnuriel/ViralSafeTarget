#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$ROOT"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"

if command -v vst >/dev/null 2>&1; then
  VST=(vst)
elif command -v viral-safe-target >/dev/null 2>&1; then
  VST=(viral-safe-target)
else
  VST=(python -m viral_safe_target)
fi

OUT=reports/synthetic_e2e
mkdir -p "$OUT"

"${VST[@]}" scan \
  --virus-alignment data/demo/virus_aligned.fasta \
  --reference-id HSV2_demo_ref \
  --gff data/demo/reference.gff3 \
  --out-dir "$OUT" \
  --min-coverage 0.0

"${VST[@]}" build-offtarget-input \
  --candidates "$OUT/candidates_ranked_pre_human.csv" \
  --human-fasta-directory data/demo \
  --output "$OUT/cas_offinder_input.txt" \
  --manifest "$OUT/offtarget_selected_candidates.csv" \
  --max-candidates 8

"${VST[@]}" summarize-offtargets \
  --candidates "$OUT/candidates_ranked_pre_human.csv" \
  --cas-output tests/fixtures/cas_offinder_output.tsv \
  --manifest "$OUT/offtarget_selected_candidates.csv" \
  --out-dir "$OUT"

"${VST[@]}" simulate-pairs \
  --candidates "$OUT/candidates_ranked_post_human.csv" \
  --gff data/demo/reference.gff3 \
  --virus-alignment data/demo/virus_aligned.fasta \
  --reference-id HSV2_demo_ref \
  --out-dir "$OUT" \
  --min-distance 1 \
  --maximum-viral-occurrence-count 10 \
  --max-candidates 8

"${VST[@]}" report \
  --candidates "$OUT/candidates_ranked_post_human.csv" \
  --rejected "$OUT/candidates_rejected_pre_human.csv" \
  --pairs "$OUT/pair_hypotheses_same_gene.csv" \
  --predicted-hits "$OUT/predicted_human_hits.csv" \
  --out-dir "$OUT" \
  --title "ViralSafeTarget synthetic end-to-end verification"

python - <<'PY'
from pathlib import Path
import pandas as pd

out = Path("reports/synthetic_e2e")
candidates = pd.read_csv(out / "candidates_ranked_pre_human.csv")
assert candidates["pre_human_score"].nunique() > 1
assert candidates["candidate_id"].is_unique
assert (out / "run_manifest.json").is_file()
multi = pd.read_csv(out / "pair_hypotheses_multi_target.csv")
if not multi.empty:
    assert multi["deletion_length_bp"].isna().all()
print("Synthetic end-to-end verification passed")
PY
