from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
import yaml

from viral_safe_target.provenance import sha256_file
from viral_safe_target.tool_benchmark import (
    freeze_panel,
    parse_crispritz_profile,
    run_ablation,
    run_tool_benchmark,
)


def _candidates() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "candidate_id": ["g1", "g2", "g3"],
            "guide_sequence": ["A" * 20, "C" * 20, "G" * 20],
            "pam": ["AGG", "CGG", "TGG"],
            "gene_name": ["a", "b", "c"],
            "pre_human_score": [0.9, 0.8, 0.7],
            "post_human_score": [0.9, 0.4, 0.7],
            "human_total_predicted_hits": [0, 2, 0],
            "conservation_score": [1.0, 0.8, 0.9],
            "viral_uniqueness_score": [1.0, 1.0, 0.5],
            "low_complexity_penalty": [0.0, 0.0, 0.0],
        }
    )


def test_freeze_panel_is_deterministic_and_rejects_duplicate_guides() -> None:
    candidates = _candidates().iloc[::-1]
    frozen = freeze_panel(candidates, expected_count=3)
    assert frozen["candidate_id"].tolist() == ["g1", "g2", "g3"]
    duplicated = candidates.copy()
    duplicated.loc[duplicated.index[1], "guide_sequence"] = duplicated.iloc[0]["guide_sequence"]
    with pytest.raises(ValueError, match="guide sequences must be unique"):
        freeze_panel(duplicated)


def test_parse_crispritz_profile_maps_counts_and_preserves_missing(tmp_path: Path) -> None:
    profile = tmp_path / "result.profile.xls"
    profile.write_text(
        "GUIDE\tONT\tOFFT\t0MM\t1MM\t2MM\t3MM\n"
        f"{'A' * 20}NNN\t0\t3\t0\t1\t2\t0\n"
        f"{'C' * 20}NNN\t0\t1\t0\t0\t0\t1\n",
        encoding="utf-8",
    )
    parsed = parse_crispritz_profile(profile, _candidates())
    assert parsed["predicted_offtarget_burden"].tolist()[:2] == [3.0, 1.0]
    assert pd.isna(parsed.loc[2, "predicted_offtarget_burden"])


def test_ablation_reports_rank_sensitivity_without_therapeutic_score() -> None:
    detail, summary = run_ablation(
        _candidates(),
        [
            {"name": "conservation", "column": "conservation_score", "weight": 0.7},
            {"name": "uniqueness", "column": "viral_uniqueness_score", "weight": 0.3},
        ],
        ["low_complexity_penalty"],
        [2],
    )
    assert set(detail["variant"]) == {
        "all_components",
        "without_conservation",
        "without_uniqueness",
    }
    assert "maximum_absolute_rank_shift" in summary
    assert "therapeutic_score" not in detail


def test_end_to_end_benchmark_keeps_unavailable_tools_pending(tmp_path: Path) -> None:
    candidates = _candidates()
    candidates.to_csv(tmp_path / "candidates.csv", index=False)
    ranking = {"ranking": {"weights": {"conservation": 0.7, "viral_uniqueness": 0.3}}}
    (tmp_path / "ranking.yaml").write_text(yaml.safe_dump(ranking), encoding="utf-8")
    (tmp_path / "capabilities.tsv").write_text(
        "tool_name\tcapability\tstatus\tnote\tofficial_source\taccessed_date\n"
        "ViralSafeTarget\tguide_design\tyes\ttest\thttps://example.org\t2026-07-15\n",
        encoding="utf-8",
    )
    config = {
        "schema_version": "1.0",
        "benchmark_id": "synthetic",
        "project_root": ".",
        "candidate_table": "candidates.csv",
        "expected_candidate_count": 3,
        "output_dir": "output",
        "capability_evidence_table": "capabilities.tsv",
        "top_k": [2],
        "tools": [
            {
                "id": "viral_safe_target_pre_human",
                "version": "test",
                "mode": "internal",
                "status": "completed",
                "official_source": "https://example.org",
            },
            {
                "id": "cas-offinder",
                "version": "test",
                "mode": "internal",
                "status": "completed",
                "official_source": "https://example.org",
            },
            {
                "id": "crispor",
                "version": "not run",
                "mode": "export",
                "status": "export_required",
                "official_source": "https://example.org",
            },
        ],
        "ablation": {
            "ranking_config": "ranking.yaml",
            "components": {
                "conservation": "conservation_score",
                "viral_uniqueness": "viral_uniqueness_score",
            },
            "penalty_columns": ["low_complexity_penalty"],
        },
    }
    config_path = tmp_path / "benchmark.yaml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
    result = run_tool_benchmark(config_path)
    status = pd.read_csv(tmp_path / "output" / "tool_execution_status.csv").set_index("tool_name")
    assert result["candidate_count"] == 3
    assert status.loc["crispor", "status"] == "export_required"
    assert status.loc["crispor", "reported_candidates"] == 0
    assert (tmp_path / "output" / "multitool_benchmark_report.html").is_file()


def test_public_hsv2_benchmark_snapshot_matches_committed_raw_output() -> None:
    root = Path(__file__).resolve().parents[1]
    output = root / "reports" / "hsv2_tool_benchmark"
    status = pd.read_csv(output / "tool_execution_status.csv").set_index("tool_name")
    agreement = pd.read_csv(output / "rank_agreement.csv")
    metadata = json.loads(
        (output / "raw" / "crispritz" / "run_metadata.json").read_text(encoding="utf-8")
    )
    manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    profile = output / "raw" / "crispritz" / "hsv2_257_grch38p14.profile.xls"

    assert manifest["candidate_count"] == 257
    assert manifest["unique_guide_count"] == 257
    assert status.loc["crispritz", "status"] == "completed"
    assert status.loc["crispritz", "reported_candidates"] == 257
    assert status.loc["crispor", "status"] == "export_required"
    assert status.loc["chopchop", "status"] == "export_required"
    assert status.loc["guidescan2", "status"] == "export_required"
    assert sum(1 for _ in profile.open(encoding="utf-8")) == 258
    assert sha256_file(profile) == metadata["profile_sha256"]
    row = agreement[
        agreement["tool_a"].eq("cas-offinder") & agreement["tool_b"].eq("crispritz")
    ].iloc[0]
    assert row["shared_candidates"] == 257
    assert row["spearman_rank_correlation"] == pytest.approx(0.8804128188049333)
