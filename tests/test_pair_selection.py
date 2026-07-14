from __future__ import annotations

import pandas as pd

from viral_safe_target.disruption import (
    rank_candidate_pairs,
    select_pair_candidates,
    simulate_candidate_pair,
)


def _candidates() -> pd.DataFrame:
    rows = []
    for index, (gene, coordinate, score) in enumerate(
        [
            ("G2", 500, 0.95),
            ("G1", 100, 0.90),
            ("G2", 600, 0.85),
            ("G1", 200, 0.80),
            ("G1", 300, 0.70),
        ]
    ):
        rows.append(
            {
                "candidate_id": f"c{index}",
                "reference_accession": "ref",
                "reference_start_1based": coordinate,
                "reference_end_1based": coordinate + 19,
                "site_start_1based": coordinate,
                "site_end_1based": coordinate + 22,
                "reference_site_plus_strand": "A" * 20 + "AGG",
                "strand": "+",
                "guide_sequence": ("ACGT" * 4) + f"{index:04d}".replace("0", "A")[-4:],
                "pam": "AGG",
                "gene_name": gene,
                "feature_type": "CDS",
                "virus_site_coverage": 1.0,
                "exact_strain_coverage": 1.0,
                "gc_fraction": 0.5,
                "reference_viral_occurrence_count": 1,
                "pre_human_score": score,
                "rejection_reasons": "",
            }
        )
    return pd.DataFrame(rows)


def test_pair_candidate_selection_is_stratified_and_row_order_independent():
    candidates = _candidates()
    first = select_pair_candidates(candidates, maximum_total_candidates=4)
    second = select_pair_candidates(
        candidates.sample(frac=1, random_state=42), maximum_total_candidates=4
    )
    assert first["candidate_id"].tolist() == second["candidate_id"].tolist()
    assert set(first.iloc[:2]["gene_name"]) == {"G1", "G2"}


def test_same_gene_deletion_calculation_and_cross_gene_prevention():
    candidates = _candidates()
    same = simulate_candidate_pair(candidates.iloc[1], candidates.iloc[3])
    assert same["hypothesis_type"] == "same_gene_deletion_hypothesis"
    assert same["deletion_length_bp"] == same["distance_bp"]
    assert same["coding_frame_disruption_heuristic"]
    cross = simulate_candidate_pair(candidates.iloc[0], candidates.iloc[1])
    assert cross["hypothesis_type"] == "multi_target_hypothesis"
    assert pd.isna(cross["deletion_length_bp"])
    assert "no single intervening deletion" in cross["interpretation"]


def test_pair_ranking_does_not_use_first_rows():
    candidates = _candidates()
    pairs = rank_candidate_pairs(
        candidates,
        same_feature_only=True,
        min_distance_bp=1,
        max_candidates=4,
        maximum_viral_occurrence_count=1,
    )
    selected_ids = set(pairs["candidate_a"]) | set(pairs["candidate_b"])
    assert "c0" in selected_ids
    assert "c1" in selected_ids
