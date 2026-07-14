from __future__ import annotations

import json

from viral_safe_target.cache import stage_is_current, write_stage_stamp
from viral_safe_target.provenance import write_run_manifest


def test_cache_invalidates_when_input_changes(tmp_path):
    source = tmp_path / "source.txt"
    output = tmp_path / "output.txt"
    stamp = tmp_path / "stamp.json"
    source.write_text("A", encoding="utf-8")
    output.write_text("result", encoding="utf-8")
    write_stage_stamp(stamp, [source], {"mode": "test"})
    assert stage_is_current(stamp, [output], [source], {"mode": "test"})
    source.write_text("B", encoding="utf-8")
    assert not stage_is_current(stamp, [output], [source], {"mode": "test"})


def test_manifest_contains_research_provenance(tmp_path):
    source = tmp_path / "source.txt"
    result = tmp_path / "result.txt"
    config = tmp_path / "config.yaml"
    source.write_text("ACGT", encoding="utf-8")
    result.write_text("result", encoding="utf-8")
    config.write_text("schema_version: '0.3'", encoding="utf-8")
    manifest = tmp_path / "manifest.json"
    write_run_manifest(
        manifest,
        [source],
        {"threshold": 0.95},
        config_path=config,
        editor_profile={"name": "SpCas9"},
        accepted_accessions=["A"],
        rejected_accessions=[{"accession": "B", "rejection_reason": "QC"}],
        human_assembly_identifier="GRCh38.p14",
        random_seed=7,
        output_paths=[result],
    )
    payload = json.loads(manifest.read_text(encoding="utf-8"))
    assert payload["configuration"]["sha256"]
    assert payload["editor_profile"]["name"] == "SpCas9"
    assert payload["accepted_accessions"] == ["A"]
    assert payload["rejected_accessions"][0]["rejection_reason"] == "QC"
    assert payload["outputs"][0]["sha256"]
    assert "dirty_worktree" in payload
