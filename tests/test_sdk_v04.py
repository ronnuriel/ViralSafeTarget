from __future__ import annotations

import json

import pandas as pd

import viral_safe_target as vst
from viral_safe_target.tables import TOOL_RESULT_COLUMNS


def test_public_sdk_imports_and_run_loading(tmp_path):
    candidates = pd.DataFrame(
        [{"candidate_id": "c1", "guide_sequence": "A" * 20, "gene_name": "UL19"}]
    )
    candidates.to_csv(tmp_path / "candidates_ranked_post_human.csv", index=False)
    pd.DataFrame([{"candidate_id": "c1", "mismatches": 2}]).to_csv(
        tmp_path / "predicted_human_hits.csv", index=False
    )
    pd.DataFrame([{"candidate_a": "c1", "candidate_b": "c2"}]).to_csv(
        tmp_path / "pair_hypotheses_same_gene.csv", index=False
    )
    (tmp_path / "run_manifest.json").write_text(json.dumps({"git_commit": "abc"}))

    run = vst.load_run(tmp_path)
    assert run.candidates.iloc[0]["candidate_id"] == "c1"
    assert len(run.human_hits) == 1
    assert len(run.same_gene_pairs) == 1
    assert len(run.candidate_table) == 1
    assert vst.__version__ == "0.6.0"
    for public_name in [
        "CandidateTable",
        "ToolResultTable",
        "ToolAdapter",
        "CandidateScorer",
        "compare_tools",
        "build_consensus",
        "load_external_results",
    ]:
        assert hasattr(vst, public_name)


def test_normalized_tool_schema_preserves_missing_values():
    frame = pd.DataFrame([{"candidate_id": "c1", "tool_name": "x"}])
    table = vst.ToolResultTable.from_frame(frame)
    assert list(table.dataframe.columns) == TOOL_RESULT_COLUMNS
    assert pd.isna(table.dataframe.iloc[0]["raw_value"])
    assert pd.isna(table.dataframe.iloc[0]["normalized_value"])
