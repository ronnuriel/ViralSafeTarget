#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

vst profiles validate \
  --virus-profile configs/viruses/hsv2.yaml \
  --host-profile configs/hosts/human_grch38.yaml \
  --nuclease-profile configs/nucleases/spcas9.yaml \
  --require-host-reference

vst analyze gene-function --out-dir reports/hsv2_gene_function

vst showcase build \
  --virus-profile configs/viruses/hsv2.yaml \
  --host-profile configs/hosts/human_grch38.yaml \
  --nuclease-profile configs/nucleases/spcas9.yaml \
  --out-dir reports/hsv2_showcase

echo "Showcase report: $ROOT/reports/hsv2_showcase/FINAL_REPORT.html"
