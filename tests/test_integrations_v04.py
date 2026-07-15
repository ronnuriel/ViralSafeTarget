from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from viral_safe_target.integrations import CasOffinderAdapter, CrispritzAdapter
from viral_safe_target.integrations.base import (
    AdapterError,
    ToolAvailability,
    detect_executable,
)

ROOT = Path(__file__).resolve().parents[1]


def _candidates() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "candidate_id": "c1",
                "guide_sequence": "A" * 20,
                "pam": "AGG",
                "gene_name": "UL19",
            },
            {
                "candidate_id": "c2",
                "guide_sequence": "C" * 20,
                "pam": "TGG",
                "gene_name": "UL30",
            },
            {
                "candidate_id": "c3",
                "guide_sequence": "G" * 20,
                "pam": "CGG",
                "gene_name": "UL30",
            },
        ]
    )


def test_adapter_detection_is_actionable_for_missing_executable(monkeypatch):
    monkeypatch.setattr("shutil.which", lambda executable: None)
    status = detect_executable("missing-tool", "missing-tool", ("--version",))
    assert not status.available
    assert "not found" in status.message


def test_crispritz_input_parse_bulges_variants_and_provenance(tmp_path):
    adapter = CrispritzAdapter()
    manifest = adapter.build_input(
        _candidates(),
        {
            "genome_or_assembly": "GRCh38.p14",
            "mismatches": 3,
            "dna_bulges": 1,
            "rna_bulges": 1,
            "variant_aware": True,
        },
        tmp_path,
    )
    settings = json.loads(manifest.read_text())
    assert settings["candidate_count"] == 3
    assert settings["variant_aware"] is True
    assert (tmp_path / "guides.tsv").is_file()

    source = ROOT / "tests/fixtures/crispritz_results.tsv"
    parsed = adapter.parse(source, tmp_path / "guides.tsv")
    normalized = adapter.normalize(
        parsed,
        candidates=_candidates(),
        source_file=source,
        version="synthetic",
        assembly="GRCh38.p14+variants",
    )
    c1 = normalized[normalized["candidate_id"].eq("c1")].iloc[0]
    counts = json.loads(c1["explanation"])
    assert counts["mismatch_2"] == 1
    assert counts["dna_bulge_hits"] == 1
    assert counts["variant_enriched_hits"] == 1
    assert counts["annotation_counts"] == {"exon": 1, "intron": 1}
    assert counts["populations"] == ["EUR"]
    assert counts["samples"] == ["S1"]
    c3 = normalized[normalized["candidate_id"].eq("c3")].iloc[0]
    assert c3["raw_value"] == 0
    assert c3["source_file_sha256"]


def test_missing_external_output_is_not_zero_risk(tmp_path):
    with pytest.raises(AdapterError, match="not a zero-hit result"):
        CasOffinderAdapter().parse(tmp_path / "missing.tsv")


def test_native_cas_offinder_column_order_regression():
    parsed = CasOffinderAdapter().parse(ROOT / "tests/fixtures/cas_offinder_output.tsv")
    assert parsed.iloc[0]["query"].endswith("NNN")
    assert parsed.iloc[0]["chromosome"] == "chrSynthetic"
    assert parsed.iloc[0]["location_0based"] == 100


def test_crispritz_dry_run_records_command(monkeypatch, tmp_path):
    adapter = CrispritzAdapter()
    manifest = adapter.build_input(
        _candidates(),
        {"reference_genome": "genome", "pam_file": "pam.txt"},
        tmp_path,
    )
    monkeypatch.setattr(
        adapter,
        "detect",
        lambda: ToolAvailability("crispritz", True, "synthetic", "/usr/bin/crispritz.py", "native"),
    )
    execution = adapter.run(manifest, tmp_path, dry_run=True)
    assert execution.command[0] == "/usr/bin/crispritz.py"
    assert "search" in execution.command
    assert execution.returncode == 0
