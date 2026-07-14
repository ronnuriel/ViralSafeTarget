from __future__ import annotations

import pandas as pd

from viral_safe_target.discovery_comparison import compare_discovery_modes


def test_discovery_comparison_reports_mode_sensitivity_without_safety_claim() -> None:
    balanced_candidates = pd.DataFrame(
        [
            {
                "candidate_id": "shared",
                "mapped_gene_names": "A",
                "post_human_rank": 2,
                "post_human_score": 0.8,
                "human_total_predicted_hits": 0,
                "decision": "retain",
                "screening_status": "completed",
            }
        ]
    )
    exhaustive_candidates = pd.concat(
        [
            balanced_candidates.assign(post_human_rank=3),
            pd.DataFrame(
                [
                    {
                        "candidate_id": "new",
                        "mapped_gene_names": "B",
                        "post_human_rank": 1,
                        "post_human_score": 0.9,
                        "human_total_predicted_hits": 0,
                        "decision": "retain",
                        "screening_status": "completed",
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    balanced_genes = pd.DataFrame(
        [{"gene_name": "A", "targetability_rank": 1, "targetability_score": 0.8}]
    )
    exhaustive_genes = pd.DataFrame(
        [
            {"gene_name": "B", "targetability_rank": 1, "targetability_score": 0.9},
            {"gene_name": "A", "targetability_rank": 2, "targetability_score": 0.8},
        ]
    )
    candidates, genes, summary = compare_discovery_modes(
        balanced_candidates, exhaustive_candidates, balanced_genes, exhaustive_genes
    )
    assert summary["exhaustive_top_gene"] == "B"
    assert summary["candidate_overlap_count"] == 1
    assert candidates.set_index("candidate_id").loc["new", "_merge"] == "right_only"
    assert genes.set_index("gene_name").loc["A", "targetability_rank_exhaustive"] == 2
    assert "not proof of safety" in summary["interpretation"]
