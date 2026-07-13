from __future__ import annotations

import pandas as pd

from viral_safe_target import ExampleRuleScorer, build_consensus, compare_tools
from viral_safe_target.consensus import candidate_metrics_as_tool_results
from viral_safe_target.tables import ToolResultTable


def _candidates() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "candidate_id": f"c{index}",
                "guide_sequence": str(index) * 20,
                "gene_name": "UL19" if index < 3 else "UL30",
                "pre_human_score": 1 - index / 10,
                "post_human_score": 1 - index / 10,
                "human_total_predicted_hits": 0,
                "conservation_score": 1.0,
                "sequence_complexity_score": 0.5 + index / 20,
            }
            for index in range(1, 7)
        ]
    )


def _reverse_tool(candidates: pd.DataFrame) -> ToolResultTable:
    rows = []
    for index, candidate in candidates.iterrows():
        rows.append(
            {
                "candidate_id": candidate["candidate_id"],
                "guide_sequence": candidate["guide_sequence"],
                "gene_name": candidate["gene_name"],
                "tool_name": "reverse_tool",
                "metric_name": "predicted_efficiency",
                "raw_value": index,
                "normalized_value": index / (len(candidates) - 1),
                "rank": len(candidates) - index,
                "percentile_rank": index / (len(candidates) - 1),
                "status": "completed",
            }
        )
    return ToolResultTable.from_frame(pd.DataFrame(rows))


def test_rank_consensus_missing_tools_disagreement_and_determinism():
    candidates = _candidates()
    baseline = candidate_metrics_as_tool_results(candidates)
    reverse = _reverse_tool(candidates.iloc[:-1])
    first = compare_tools(
        candidates,
        [baseline, reverse],
        expected_tools=[
            "viral_safe_target_pre_human",
            "viral_safe_target_post_human",
            "cas-offinder",
            "reverse_tool",
            "pending_tool",
        ],
    )
    second = compare_tools(
        candidates,
        [baseline, reverse],
        expected_tools=[
            "viral_safe_target_pre_human",
            "viral_safe_target_post_human",
            "cas-offinder",
            "reverse_tool",
            "pending_tool",
        ],
    )
    pd.testing.assert_frame_equal(first.consensus_candidates, second.consensus_candidates)
    assert first.consensus_candidates["tools_missing"].min() >= 1
    assert (
        first.tool_coverage.query("tool_name == 'pending_tool'")["reported_candidates"].iloc[0] == 0
    )
    assert not first.model_agreement.empty
    assert {5, 10, 20} == set(first.model_agreement["top_k"])
    assert not first.disagreement_report.empty
    assert first.consensus_candidates["consensus_score"].notna().all()
    assert build_consensus(candidates, [baseline]).iloc[0]["consensus_rank"] == 1


def test_weighted_borda_uses_percentiles_not_raw_values():
    candidates = _candidates().iloc[:2]
    baseline = candidate_metrics_as_tool_results(candidates)
    baseline.dataframe.loc[:, "raw_value"] = [10**12] * len(baseline.dataframe)
    result = compare_tools(candidates, [baseline])
    assert result.consensus_candidates["consensus_score"].between(0, 1).all()


def test_custom_rule_scorer_contract():
    scored = ExampleRuleScorer().score(_candidates())
    assert list(scored.columns) == [
        "candidate_id",
        "scorer_name",
        "raw_score",
        "confidence",
        "explanation",
    ]
    assert scored["raw_score"].between(0, 1).all()
