from __future__ import annotations

from collections import OrderedDict

import pandas as pd
import pytest

from viral_safe_target import (
    EditorProfile,
    rank_pre_human_candidates,
    scan_spcas9_candidates,
)
from viral_safe_target.scoring import gc_score, sequence_complexity_score


def test_stable_candidate_ids_do_not_depend_on_record_order():
    first = OrderedDict(
        [
            ("ref", "ACGTACGTACGTACGTACGTAGGTTTT"),
            ("strain", "ACGTACGTACGTACGTACGTAGGTTTT"),
        ]
    )
    second = OrderedDict(reversed(list(first.items())))
    a = scan_spcas9_candidates(first, "ref")
    b = scan_spcas9_candidates(second, "ref")
    assert set(a["candidate_id"]) == set(b["candidate_id"])


def test_duplicate_guides_are_retained_and_counted_across_genomes():
    site = "ACGTACGTACGTACGTACGTAGG"
    records = OrderedDict([("ref", site + "TTTT" + site), ("strain", site + "TTTT" + site)])
    candidates = scan_spcas9_candidates(records, "ref")
    repeated = candidates[candidates["guide_sequence"] == "ACGTACGTACGTACGTACGT"]
    assert len(repeated) == 2
    assert (repeated["reference_viral_occurrence_count"] == 2).all()
    assert (repeated["all_viral_occurrence_count"] == 4).all()
    assert (repeated["guide_genome_presence_count"] == 2).all()
    assert repeated["several_coordinates_share_guide"].all()
    assert not repeated["is_guide_unique_reference"].any()
    assert repeated["duplicate_handling"].eq("retained_all_coordinates").all()


def test_gc_and_complexity_components_are_conservative_and_configured():
    settings = {
        "accepted_min": 0.25,
        "preferred_min": 0.40,
        "preferred_max": 0.65,
        "accepted_max": 0.80,
    }
    assert gc_score(0.5, settings) == 1.0
    assert gc_score(0.1, settings) == 0.0
    assert gc_score(0.9, settings) == 0.0
    assert sequence_complexity_score("ACGT" * 5) == pytest.approx(1.0)
    assert sequence_complexity_score("A" * 20) == 0.0


def _candidate_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "candidate_id": "b",
                "reference_accession": "ref",
                "guide_sequence": "ACGT" * 5,
                "pam": "AGG",
                "virus_site_coverage": 1.0,
                "reference_viral_occurrence_count": 1,
                "feature_type": "CDS",
                "gene_name": "G1",
                "product": "",
            },
            {
                "candidate_id": "a",
                "reference_accession": "ref",
                "guide_sequence": "A" * 20,
                "pam": "AGG",
                "virus_site_coverage": 0.8,
                "reference_viral_occurrence_count": 2,
                "feature_type": "intergenic_or_unannotated",
                "gene_name": "",
                "product": "",
            },
        ]
    )


def test_missing_evidence_remains_missing_ranking_is_deterministic_and_reasons_are_explicit():
    first = rank_pre_human_candidates(_candidate_frame())
    second = rank_pre_human_candidates(_candidate_frame().sample(frac=1, random_state=9))
    assert first["gene_evidence_score"].isna().all()
    pd.testing.assert_series_equal(
        first.set_index("candidate_id")["pre_human_score"].sort_index(),
        second.set_index("candidate_id")["pre_human_score"].sort_index(),
    )
    rejected = first[first["candidate_id"] == "a"].iloc[0]
    assert "sequence_complexity_below_configured_minimum" in rejected["rejection_reasons"]
    assert "curated gene evidence unavailable" in first.iloc[0]["rank_explanation"]


def test_editor_profile_validation_and_reverse_strand_coordinates():
    with pytest.raises(ValueError, match="pam_orientation"):
        EditorProfile.from_mapping(
            {
                "name": "bad",
                "protospacer_length": 20,
                "pam_pattern": "NGG",
                "pam_orientation": "sideways",
                "cut_offset": 3,
                "mismatch_search_threshold": 3,
            }
        )
    records = OrderedDict([("ref", "CCG" + "ACGT" * 5), ("strain", "CCG" + "ACGT" * 5)])
    minus = scan_spcas9_candidates(records, "ref")
    assert (minus["strand"] == "-").any()
    row = minus[minus["strand"] == "-"].iloc[0]
    assert row["reference_start_1based"] == 4
    assert row["reference_end_1based"] == 23
