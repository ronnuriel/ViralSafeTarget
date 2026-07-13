from pathlib import Path

from viral_safe_target import (
    annotate_candidates,
    rank_candidate_pairs,
    read_fasta,
    read_gff3,
    scan_spcas9_candidates,
    spcas9_cut_after_1based,
)

ROOT = Path(__file__).resolve().parents[1]


def _demo_candidates():
    records = read_fasta(ROOT / "data/demo/virus_aligned.fasta")
    candidates = scan_spcas9_candidates(records, "HSV2_demo_ref", min_site_coverage=0.0)
    features = read_gff3(ROOT / "data/demo/reference.gff3")
    return records, features, annotate_candidates(candidates, features, "HSV2_demo_ref")


def test_cut_coordinate_is_inside_protospacer():
    _, _, candidates = _demo_candidates()
    for _, row in candidates.iterrows():
        cut = spcas9_cut_after_1based(row)
        assert int(row["reference_start_1based"]) <= cut
        assert cut < int(row["reference_end_1based"])


def test_pair_simulation_is_ordered_and_bounded():
    records, features, candidates = _demo_candidates()
    pairs = rank_candidate_pairs(
        candidates,
        features=features,
        aligned_records=records,
        reference_id="HSV2_demo_ref",
        same_feature_only=False,
        min_distance_bp=1,
        max_distance_bp=10_000,
    )
    assert not pairs.empty
    assert (pairs["deletion_length_bp"] > 0).all()
    assert pairs["exact_pair_coverage"].between(0, 1).all()
