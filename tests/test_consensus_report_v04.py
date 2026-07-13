from __future__ import annotations

import pandas as pd

from viral_safe_target.consensus import candidate_metrics_as_tool_results, compare_tools
from viral_safe_target.consensus_reporting import write_consensus_report


def test_partial_tool_report_is_created_with_pending_and_experimental_sections(tmp_path):
    candidates = pd.DataFrame(
        [
            {
                "candidate_id": "c1",
                "guide_sequence": "A" * 20,
                "gene_name": "UL19",
                "pre_human_score": 0.9,
                "post_human_score": 0.8,
                "human_total_predicted_hits": 0,
            }
        ]
    )
    result = compare_tools(
        candidates,
        [candidate_metrics_as_tool_results(candidates)],
        expected_tools=[
            "viral_safe_target_pre_human",
            "viral_safe_target_post_human",
            "cas-offinder",
            "crispritz",
        ],
    )
    availability = pd.DataFrame(
        [
            {
                "tool_name": "crispritz",
                "available": False,
                "version": pd.NA,
                "message": "pending",
            }
        ]
    )
    output = write_consensus_report(
        candidates,
        result,
        tmp_path / "report.html",
        tool_availability=availability,
    )
    text = output.read_text(encoding="utf-8")
    assert "Pending tools or missing exports" in text
    assert "Experimental results" in text
    assert "crispritz" in text
    assert "Consensus is prioritization, not proof" in text
