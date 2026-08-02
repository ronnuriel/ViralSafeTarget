#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"
WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/vst-release-validation.XXXXXX")"
DIST_DIR="${DIST_DIR:-$WORK_DIR/dist}"
trap 'rm -rf "$WORK_DIR"' EXIT

for module in build twine; do
  if ! "$PYTHON_BIN" -c "import $module" >/dev/null 2>&1; then
    echo "Missing Python module '$module'. Install release tools with:" >&2
    echo "  $PYTHON_BIN -m pip install --upgrade build twine" >&2
    exit 2
  fi
done

echo "Building distributions from $ROOT_DIR"
cd "$ROOT_DIR"
mkdir -p "$DIST_DIR"
"$PYTHON_BIN" -m build --outdir "$DIST_DIR"
"$PYTHON_BIN" -m twine check "$DIST_DIR"/*

echo "Installing the wheel outside the repository"
"$PYTHON_BIN" -m venv "$WORK_DIR/venv"
VENV_PYTHON="$WORK_DIR/venv/bin/python"
VST="$WORK_DIR/venv/bin/vst"
"$VENV_PYTHON" -m pip install --upgrade pip
"$VENV_PYTHON" -m pip install "$DIST_DIR"/*.whl

cd "$WORK_DIR"
"$VST" --version
"$VST" doctor --json > doctor.json
"$VST" tools setup --json > tool_setup.json
"$VST" quickstart --out demo-project
"$VST" plan demo-project/project.yaml --json > plan.json
"$VST" run demo-project/project.yaml > run.json
"$VST" status demo-project/project.yaml > status.json
"$VST" open demo-project/results --no-browser
"$VST" export demo-project/project.yaml

"$VENV_PYTHON" - <<'PY'
import json
from pathlib import Path

root = Path("demo-project/results")
required = {
    "START_HERE.html",
    "summary.json",
    "top_guides.csv",
    "top_genes.csv",
    "research_shortlist.csv",
    "stage_timings.json",
    "run_manifest.json",
    "export.zip",
}
missing = sorted(name for name in required if not (root / name).is_file())
if missing:
    raise SystemExit(f"Release smoke test is missing result files: {missing}")
summary = json.loads((root / "summary.json").read_text(encoding="utf-8"))
if summary["host_risk_status"] not in {"completed", "external_required"}:
    raise SystemExit(f"Unexpected host-risk status: {summary['host_risk_status']}")
print(f"Release smoke test passed with host-risk status: {summary['host_risk_status']}")
PY

echo "Release candidate validation passed. Temporary files were isolated in $WORK_DIR."
