from __future__ import annotations

import pandas as pd

from viral_safe_target.showcase import (
    build_comparison_sets,
    build_evidence_aware_candidates,
    build_research_findings,
    pareto_fronts,
    select_balanced_deep_panel,
)


def test_pareto_fronts_do_not_collapse_separate_objectives() -> None:
    frame = pd.DataFrame(
        {
            "targetability": [1.0, 0.2, 0.5, 0.1],
            "disruption": [0.2, 1.0, 0.5, 0.1],
        }
    )
    fronts = pareto_fronts(frame, ["targetability", "disruption"])
    assert fronts.tolist() == [1, 1, 1, 2]


def _inputs() -> tuple[pd.DataFrame, ...]:
    mapping = pd.DataFrame(
        [
            {
                "candidate_id": "a1",
                "gene_name": "A",
                "post_human_score": 0.9,
                "post_human_rank": 1,
                "exact_strain_coverage": 1.0,
                "human_total_predicted_hits": 0,
                "cut_domain_accessions": "IPR_A",
            },
            {
                "candidate_id": "a2",
                "gene_name": "A",
                "post_human_score": 0.8,
                "post_human_rank": 2,
                "exact_strain_coverage": 1.0,
                "human_total_predicted_hits": 0,
                "cut_domain_accessions": "",
            },
            {
                "candidate_id": "b1",
                "gene_name": "B",
                "post_human_score": 0.7,
                "post_human_rank": 3,
                "exact_strain_coverage": 0.9,
                "human_total_predicted_hits": 1,
                "cut_domain_accessions": "",
            },
        ]
    )
    outcomes = pd.DataFrame(
        [
            {
                "candidate_id": candidate,
                "event_class": "single_guide_indel",
                "indel_size_bp": size,
                "frameshift": size % 3 != 0,
                "retained_protein_fraction": retained,
                "premature_stop_position_aa": stop,
            }
            for candidate, retained, stop in (
                ("a1", 0.2, 10),
                ("a2", 0.8, pd.NA),
                ("b1", 0.5, 20),
            )
            for size in (-1, 1)
        ]
    )
    gene_scores = pd.DataFrame(
        [
            {
                "gene_name": "A",
                "sequence_targetability_score": 0.9,
                "evidence_based_essentiality_score": pd.NA,
                "hsv1_ortholog_essentiality_score": 1.0,
                "predicted_protein_disruption_score": 0.7,
                "evidence_coverage_score": 0.8,
            },
            {
                "gene_name": "B",
                "sequence_targetability_score": 0.7,
                "evidence_based_essentiality_score": pd.NA,
                "hsv1_ortholog_essentiality_score": pd.NA,
                "predicted_protein_disruption_score": 0.5,
                "evidence_coverage_score": 0.4,
            },
        ]
    )
    evidence = pd.DataFrame(
        [
            {
                "gene_name": "A",
                "virus_type": "HSV-1",
                "evidence_strength": "direct",
                "evidence_category": "null_mutant",
            },
            {
                "gene_name": "B",
                "virus_type": "HSV-2",
                "evidence_strength": "direct",
                "evidence_category": "knockdown",
            },
        ]
    )
    categories = pd.DataFrame(
        [
            {"gene_name": "A", "primary_category": "DNA_replication"},
            {"gene_name": "B", "primary_category": "structural_egress"},
        ]
    )
    return mapping, outcomes, gene_scores, evidence, categories


def test_hsv1_evidence_stays_separate_from_direct_hsv2_phenotype() -> None:
    candidates = build_evidence_aware_candidates(*_inputs()).set_index("candidate_id")
    assert candidates.loc["a1", "evidence_tier"] == "HSV1_ortholog_only"
    assert not bool(candidates.loc["a1", "direct_hsv2_phenotype_evidence"])
    assert candidates.loc["b1", "evidence_tier"] == "direct_HSV2_phenotype"
    assert pd.isna(candidates.loc["b1", "evidence_based_essentiality_score"])
    assert "therapeutic" in candidates.loc["a1", "priority_interpretation"]


def test_direct_hsv2_variability_is_not_mislabeled_as_phenotype() -> None:
    mapping, outcomes, gene_scores, evidence, categories = _inputs()
    evidence = pd.concat(
        [
            evidence,
            pd.DataFrame(
                [
                    {
                        "gene_name": "A",
                        "virus_type": "HSV-2",
                        "evidence_strength": "direct",
                        "evidence_category": "genetic_variability",
                    }
                ]
            ),
        ],
        ignore_index=True,
    )
    candidates = build_evidence_aware_candidates(
        mapping, outcomes, gene_scores, evidence, categories
    ).set_index("candidate_id")
    assert not bool(candidates.loc["a1", "direct_hsv2_phenotype_evidence"])
    assert candidates.loc["a1", "evidence_tier"] == "HSV1_ortholog_only"


def test_deep_panel_uses_balanced_gene_quota_and_comparison_sets_are_labeled() -> None:
    candidates = build_evidence_aware_candidates(*_inputs())
    panel = select_balanced_deep_panel(candidates, per_gene=1)
    assert panel.groupby("gene_name").size().to_dict() == {"A": 1, "B": 1}
    members, summary = build_comparison_sets(candidates)
    assert set(summary["strategy"]) == {
        "ranking_only",
        "predicted_disruption",
        "evidence_anchored",
        "mechanistically_diverse",
    }
    assert members["strategy_rationale"].str.len().gt(0).all()
    assert summary["interpretation"].str.contains("no joint editing").all()


def test_research_findings_keep_observations_and_limitations_together() -> None:
    candidates = build_evidence_aware_candidates(*_inputs())
    population = pd.DataFrame(
        [
            {
                "candidate_id": "a1",
                "population_validation_status": "exact_in_all_observable_records",
                "locus_observable_record_count": 10,
                "exact_target_in_observable_locus_count": 10,
                "observable_locus_exact_target_coverage": 1.0,
            }
        ]
    )
    gene_rankings = pd.DataFrame(
        [
            {"gene_name": "UL3", "targetability_rank": 1},
            {"gene_name": "UL30", "targetability_rank": 9},
        ]
    )
    gene_stability = pd.DataFrame([{"gene_name": "UL3", "top_10_stability": True}])
    gene_scores = pd.DataFrame(
        [
            {
                "gene_name": "UL3",
                "sequence_targetability_score": 0.9,
                "hsv2_evidence_based_essentiality_score": pd.NA,
                "hsv1_ortholog_essentiality_score": 0.0,
                "hsv1_essentiality_status": "nonessential_in_tested_cell_culture",
                "predicted_protein_disruption_score": 0.6,
            },
            {
                "gene_name": "A",
                "sequence_targetability_score": 0.8,
                "hsv2_evidence_based_essentiality_score": pd.NA,
                "hsv1_ortholog_essentiality_score": 1.0,
                "hsv1_essentiality_status": "supported_essential",
                "predicted_protein_disruption_score": 0.7,
            },
        ]
    )
    evolution = pd.DataFrame(
        [
            {"gene_name": "UL3", "mean_amino_acid_conservation": 0.99},
            {"gene_name": "A", "mean_amino_acid_conservation": 1.0},
        ]
    )
    findings = build_research_findings(
        gene_rankings,
        gene_stability,
        gene_scores,
        candidates,
        population,
        evolution,
    )
    assert {
        "genome_wide_reprioritization",
        "targetability_evidence_divergence",
        "quota_stable_gene_set",
        "direct_hsv2_essentiality_gap",
    }.issubset(set(findings["finding_id"]))
    assert findings["observation"].str.len().gt(0).all()
    assert findings["key_limitation"].str.len().gt(0).all()
