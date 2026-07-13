import json
from pathlib import Path

from viral_safe_target import write_run_manifest

ROOT = Path(__file__).resolve().parents[1]


def test_manifest_contains_checksum(tmp_path):
    input_path = tmp_path / "input.txt"
    input_path.write_text("ACGT\n", encoding="utf-8")
    output_path = tmp_path / "manifest.json"
    write_run_manifest(output_path, [input_path], {"threshold": 0.95}, project_root=ROOT)
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    assert payload["inputs"][0]["sha256"]
    assert payload["parameters"]["threshold"] == 0.95
