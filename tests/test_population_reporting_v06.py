from __future__ import annotations

import pandas as pd

from viral_safe_target.cli import build_parser
from viral_safe_target.population_reporting import build_population_comparison


def test_population_comparison_keeps_population_evidence_out_of_targetability() -> None:
    candidates = pd.DataFrame(
        [
            {
                "candidate_id": "c1",
                "mapped_gene_names": "UL30",
                "reference_viral_occurrence_count": 1,
                "post_human_rank": 1,
            },
            {
                "candidate_id": "c2",
                "mapped_gene_names": "UL30;UL31",
                "reference_viral_occurrence_count": 2,
                "post_human_rank": 2,
            },
        ]
    )
    validation = pd.DataFrame(
        [
            {
                "candidate_id": "c1",
                "population_record_count": 3,
                "locus_observable_record_count": 2,
                "exact_target_in_observable_locus_count": 2,
                "observable_locus_exact_target_coverage": 1.0,
                "observable_locus_without_exact_target_count": 0,
                "locus_unresolved_record_count": 1,
                "exact_target_anywhere_count": 2,
                "locus_validation_interpretation": "test",
            },
            {
                "candidate_id": "c2",
                "population_record_count": 3,
                "locus_observable_record_count": 2,
                "exact_target_in_observable_locus_count": 1,
                "observable_locus_exact_target_coverage": 0.5,
                "observable_locus_without_exact_target_count": 1,
                "locus_unresolved_record_count": 1,
                "exact_target_anywhere_count": 1,
                "locus_validation_interpretation": "test",
            },
        ]
    )
    comparison, genes = build_population_comparison(candidates, validation)
    indexed = comparison.set_index("candidate_id")
    assert indexed.loc["c1", "population_validation_status"] == ("exact_in_all_observable_records")
    assert indexed.loc["c2", "population_validation_status"] == ("multi_locus_attribution_limited")
    assert not comparison["population_validation_used_in_targetability_score"].any()
    assert set(genes["gene_name"]) == {"UL30", "UL31"}


def test_population_analysis_is_exposed_through_the_research_cli() -> None:
    args = build_parser().parse_args(
        [
            "analyze",
            "population",
            "--population-fasta",
            "population.fasta",
            "--reference-fasta",
            "reference.fasta",
            "--candidates",
            "candidates.csv",
            "--out-dir",
            "report",
        ]
    )
    assert args.func.__name__ == "_analyze_population"
