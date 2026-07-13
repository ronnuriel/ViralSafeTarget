"""Sequence-level disruption simulations for candidate CRISPR sites.

These functions describe *what would happen to the reference sequence* if an
idealized cut or deletion occurred at the predicted coordinates. They do not
predict delivery, editing efficiency, DNA repair frequencies, viral viability,
latency clearance, toxicity, or clinical efficacy.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from itertools import combinations

import pandas as pd

from .io_utils import require_aligned


def spcas9_cut_after_1based(candidate: Mapping[str, object]) -> int:
    """Return the canonical SpCas9 cut boundary as a 1-based 'after base' coordinate.

    SpCas9 commonly cleaves about three nucleotides upstream of the PAM. For a
    plus-strand protospacer this is after protospacer base 17; for a minus-strand
    protospacer it is the corresponding boundary near the protospacer start.
    This is an idealized coordinate, not an experimental outcome prediction.
    """
    strand = str(candidate["strand"])
    start = int(candidate["reference_start_1based"])
    end = int(candidate["reference_end_1based"])
    if strand == "+":
        return end - 3
    if strand == "-":
        return start + 2
    raise ValueError(f"Unsupported strand: {strand!r}")


def hypothetical_indel_consequences(
    candidate: Mapping[str, object],
    deletion_sizes: Sequence[int] = tuple(range(1, 11)),
    insertion_sizes: Sequence[int] = tuple(range(1, 6)),
) -> pd.DataFrame:
    """List frame consequences for hypothetical indel sizes.

    This is arithmetic only. It does not estimate how frequently any indel is
    produced in a cell.
    """
    rows: list[dict[str, object]] = []
    candidate_id = str(candidate.get("candidate_id", "candidate"))
    feature_type = str(candidate.get("feature_type", ""))
    coding = feature_type == "CDS"
    for size in deletion_sizes:
        rows.append(
            {
                "candidate_id": candidate_id,
                "event": "deletion",
                "size_bp": int(size),
                "coding_feature": coding,
                "frame_consequence": (
                    "frameshift" if coding and size % 3 else "in_frame_or_not_applicable"
                ),
            }
        )
    for size in insertion_sizes:
        rows.append(
            {
                "candidate_id": candidate_id,
                "event": "insertion",
                "size_bp": int(size),
                "coding_feature": coding,
                "frame_consequence": (
                    "frameshift" if coding and size % 3 else "in_frame_or_not_applicable"
                ),
            }
        )
    return pd.DataFrame(rows)


def _reference_position_to_alignment_index(reference_aligned: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    reference_position = 0
    for alignment_index, base in enumerate(reference_aligned):
        if base != "-":
            reference_position += 1
            mapping[reference_position] = alignment_index
    return mapping


def _exact_site_present(
    aligned_sequence: str,
    ref_to_alignment: Mapping[int, int],
    start_1based: int,
    end_1based: int,
    expected_plus_strand_site: str,
) -> bool:
    observed = "".join(
        aligned_sequence[ref_to_alignment[position]]
        for position in range(start_1based, end_1based + 1)
    )
    return observed == expected_plus_strand_site


def exact_pair_coverage(
    first: Mapping[str, object],
    second: Mapping[str, object],
    aligned_records: Mapping[str, str],
    reference_id: str,
) -> tuple[float, int]:
    """Return fraction/count of genomes containing both exact 23-nt sites."""
    require_aligned(aligned_records)
    if reference_id not in aligned_records:
        raise KeyError(f"Reference id {reference_id!r} is not present in the alignment")
    ref_to_alignment = _reference_position_to_alignment_index(aligned_records[reference_id])
    exact_count = 0
    for sequence in aligned_records.values():
        first_ok = _exact_site_present(
            sequence,
            ref_to_alignment,
            int(first["site_start_1based"]),
            int(first["site_end_1based"]),
            str(first["reference_site_plus_strand"]),
        )
        second_ok = _exact_site_present(
            sequence,
            ref_to_alignment,
            int(second["site_start_1based"]),
            int(second["site_end_1based"]),
            str(second["reference_site_plus_strand"]),
        )
        exact_count += int(first_ok and second_ok)
    return exact_count / len(aligned_records), exact_count


def _overlap_summary(
    deletion_start_1based: int,
    deletion_end_1based: int,
    features: pd.DataFrame | None,
) -> dict[str, object]:
    if features is None or features.empty:
        return {
            "overlapping_feature_count": 0,
            "overlapping_features": "",
            "max_feature_fraction_removed": None,
        }
    overlaps = features[
        (features["start"] <= deletion_end_1based)
        & (features["end"] >= deletion_start_1based)
    ].copy()
    if overlaps.empty:
        return {
            "overlapping_feature_count": 0,
            "overlapping_features": "",
            "max_feature_fraction_removed": 0.0,
        }
    labels: list[str] = []
    fractions: list[float] = []
    for _, feature in overlaps.iterrows():
        overlap_start = max(deletion_start_1based, int(feature["start"]))
        overlap_end = min(deletion_end_1based, int(feature["end"]))
        overlap_bp = max(0, overlap_end - overlap_start + 1)
        feature_bp = int(feature["end"]) - int(feature["start"]) + 1
        fractions.append(overlap_bp / feature_bp if feature_bp else 0.0)
        label = str(feature.get("name") or feature.get("feature_id") or feature["feature_type"])
        labels.append(f"{label}:{feature['feature_type']}:{overlap_bp}bp")
    return {
        "overlapping_feature_count": len(overlaps),
        "overlapping_features": "; ".join(labels),
        "max_feature_fraction_removed": max(fractions) if fractions else 0.0,
    }


def simulate_candidate_pair(
    first: Mapping[str, object],
    second: Mapping[str, object],
    features: pd.DataFrame | None = None,
    aligned_records: Mapping[str, str] | None = None,
    reference_id: str | None = None,
) -> dict[str, object]:
    """Describe an idealized deletion between two canonical SpCas9 cut sites."""
    first_cut = spcas9_cut_after_1based(first)
    second_cut = spcas9_cut_after_1based(second)
    left_cut, right_cut = sorted((first_cut, second_cut))
    deletion_start = left_cut + 1
    deletion_end = right_cut
    deletion_length = max(0, right_cut - left_cut)

    result: dict[str, object] = {
        "candidate_a": str(first.get("candidate_id", "A")),
        "candidate_b": str(second.get("candidate_id", "B")),
        "cut_a_after_1based": first_cut,
        "cut_b_after_1based": second_cut,
        "deletion_start_1based": deletion_start,
        "deletion_end_1based": deletion_end,
        "deletion_length_bp": deletion_length,
        "deletion_length_mod_3": deletion_length % 3,
        "pair_same_annotated_feature": (
            bool(first.get("feature_id"))
            and str(first.get("feature_id")) == str(second.get("feature_id"))
        ),
        "individual_coverage_lower_bound": min(
            float(first.get("virus_site_coverage", 0.0)),
            float(second.get("virus_site_coverage", 0.0)),
        ),
    }
    result.update(_overlap_summary(deletion_start, deletion_end, features))

    if aligned_records is not None:
        if not reference_id:
            raise ValueError("reference_id is required when aligned_records are provided")
        coverage, count = exact_pair_coverage(
            first, second, aligned_records=aligned_records, reference_id=reference_id
        )
        result["exact_pair_coverage"] = coverage
        result["exact_pair_genome_count"] = count
        result["genome_count"] = len(aligned_records)
    return result


def rank_candidate_pairs(
    candidates: pd.DataFrame,
    features: pd.DataFrame | None = None,
    aligned_records: Mapping[str, str] | None = None,
    reference_id: str | None = None,
    same_feature_only: bool = True,
    min_distance_bp: int = 20,
    max_distance_bp: int = 10_000,
    max_candidates: int = 250,
) -> pd.DataFrame:
    """Enumerate and rank idealized two-cut deletions.

    The ranking is intentionally transparent and not biologically validated.
    It prioritizes exact strain coverage, overlap with an annotated feature, and
    a non-trivial deletion size. Researchers must validate all biological claims.
    """
    if candidates.empty:
        return pd.DataFrame()
    working = candidates.head(max_candidates).copy()
    rows: list[dict[str, object]] = []
    for (_, first), (_, second) in combinations(working.iterrows(), 2):
        if same_feature_only:
            first_feature = str(first.get("feature_id", ""))
            second_feature = str(second.get("feature_id", ""))
            if not first_feature or first_feature != second_feature:
                continue
        distance = abs(spcas9_cut_after_1based(first) - spcas9_cut_after_1based(second))
        if distance < min_distance_bp or distance > max_distance_bp:
            continue
        rows.append(
            simulate_candidate_pair(
                first,
                second,
                features=features,
                aligned_records=aligned_records,
                reference_id=reference_id,
            )
        )
    if not rows:
        return pd.DataFrame()
    pairs = pd.DataFrame(rows)
    coverage_column = (
        "exact_pair_coverage"
        if "exact_pair_coverage" in pairs.columns
        else "individual_coverage_lower_bound"
    )
    feature_fraction = pairs["max_feature_fraction_removed"].fillna(0.0).clip(0.0, 1.0)
    length_component = (pairs["deletion_length_bp"].clip(upper=1_000) / 1_000).astype(float)
    pairs["sequence_disruption_score"] = (
        0.65 * pairs[coverage_column].astype(float)
        + 0.25 * feature_fraction.astype(float)
        + 0.10 * length_component
    )
    pairs["interpretation"] = (
        "sequence_level_candidate_only_not_viability_prediction"
    )
    return pairs.sort_values(
        ["sequence_disruption_score", coverage_column, "deletion_length_bp"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
