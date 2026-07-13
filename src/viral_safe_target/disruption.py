"""Deterministic sequence-level pair hypotheses.

These calculations do not predict delivery, editing efficiency, repair
frequencies, viral viability, toxicity, or clinical efficacy.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from itertools import combinations
from pathlib import Path
from typing import Any

import pandas as pd

from .config import EditorProfile, get_editor, load_config
from .io_utils import require_aligned


def cut_after_1based(candidate: Mapping[str, object], editor: EditorProfile) -> int:
    strand = str(candidate["strand"])
    start = int(candidate["reference_start_1based"])
    end = int(candidate["reference_end_1based"])
    if strand == "+":
        return end - editor.cut_offset
    if strand == "-":
        return start + editor.cut_offset - 1
    raise ValueError(f"Unsupported strand: {strand!r}")


def spcas9_cut_after_1based(candidate: Mapping[str, object]) -> int:
    """Backward-compatible canonical SpCas9 cut-boundary calculation."""
    return cut_after_1based(candidate, get_editor(load_config()))


def hypothetical_indel_consequences(
    candidate: Mapping[str, object],
    deletion_sizes: Sequence[int] = tuple(range(1, 11)),
    insertion_sizes: Sequence[int] = tuple(range(1, 6)),
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    candidate_id = str(candidate.get("candidate_id", "candidate"))
    coding = str(candidate.get("feature_type", "")) == "CDS"
    for event, sizes in (("deletion", deletion_sizes), ("insertion", insertion_sizes)):
        for size in sizes:
            rows.append(
                {
                    "candidate_id": candidate_id,
                    "event": event,
                    "size_bp": int(size),
                    "coding_feature": coding,
                    "frame_consequence": (
                        "frameshift_heuristic"
                        if coding and size % 3
                        else "in_frame_or_not_applicable"
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
    require_aligned(aligned_records)
    if reference_id not in aligned_records:
        raise KeyError(f"Reference id {reference_id!r} is not present in the alignment")
    ref_to_alignment = _reference_position_to_alignment_index(aligned_records[reference_id])
    exact_count = 0
    for sequence in aligned_records.values():
        exact_count += int(
            _exact_site_present(
                sequence,
                ref_to_alignment,
                int(first["site_start_1based"]),
                int(first["site_end_1based"]),
                str(first["reference_site_plus_strand"]),
            )
            and _exact_site_present(
                sequence,
                ref_to_alignment,
                int(second["site_start_1based"]),
                int(second["site_end_1based"]),
                str(second["reference_site_plus_strand"]),
            )
        )
    return exact_count / len(aligned_records), exact_count


def _configured(config: dict[str, Any] | str | Path | None) -> dict[str, Any]:
    return config if isinstance(config, dict) else load_config(config)


def select_pair_candidates(
    candidates: pd.DataFrame,
    config: dict[str, Any] | str | Path | None = None,
    *,
    genes: Sequence[str] | None = None,
    feature_types: Sequence[str] | None = None,
    minimum_conservation: float | None = None,
    gc_range: tuple[float, float] | None = None,
    maximum_viral_occurrence_count: int | None = None,
    maximum_candidates_per_gene: int | None = None,
    maximum_total_candidates: int | None = None,
    stratify_by_gene: bool | None = None,
) -> pd.DataFrame:
    """Filter, rank, deduplicate by guide, then deterministically stratify."""
    if candidates.empty:
        return candidates.copy()
    settings = _configured(config)["pair_selection"]
    genes = list(genes if genes is not None else settings["genes"])
    feature_types = list(feature_types if feature_types is not None else settings["feature_types"])
    minimum_conservation = float(
        settings["minimum_conservation"] if minimum_conservation is None else minimum_conservation
    )
    gc_range = gc_range or (float(settings["gc_min"]), float(settings["gc_max"]))
    maximum_viral_occurrence_count = int(
        settings["maximum_viral_occurrence_count"]
        if maximum_viral_occurrence_count is None
        else maximum_viral_occurrence_count
    )
    maximum_candidates_per_gene = int(
        settings["maximum_candidates_per_gene"]
        if maximum_candidates_per_gene is None
        else maximum_candidates_per_gene
    )
    maximum_total_candidates = int(
        settings["maximum_total_candidates"]
        if maximum_total_candidates is None
        else maximum_total_candidates
    )
    stratify_by_gene = (
        bool(settings["stratify_by_gene"]) if stratify_by_gene is None else stratify_by_gene
    )

    working = candidates.copy()
    coverage_column = (
        "exact_strain_coverage" if "exact_strain_coverage" in working else "virus_site_coverage"
    )
    score_column = "post_human_score" if "post_human_score" in working else "pre_human_score"
    if score_column not in working:
        working[score_column] = pd.to_numeric(working[coverage_column], errors="coerce")
    if "gc_fraction" not in working:
        working["gc_fraction"] = (
            working["guide_sequence"].str.count("[GC]") / working["guide_sequence"].str.len()
        )
    if "reference_viral_occurrence_count" not in working:
        working["reference_viral_occurrence_count"] = 1
    if "gene_name" not in working:
        working["gene_name"] = ""
    if "feature_type" not in working:
        working["feature_type"] = "intergenic_or_unannotated"
    working = working[
        pd.to_numeric(working[coverage_column], errors="coerce").fillna(0) >= minimum_conservation
    ]
    working = working[working["gc_fraction"].between(*gc_range)]
    working = working[
        pd.to_numeric(working["reference_viral_occurrence_count"], errors="coerce").fillna(10**9)
        <= maximum_viral_occurrence_count
    ]
    if genes:
        working = working[working["gene_name"].fillna("").isin(genes)]
    if feature_types:
        working = working[working["feature_type"].isin(feature_types)]
    if "rejection_reasons" in working:
        working = working[working["rejection_reasons"].fillna("").eq("")]
    working = working.sort_values(
        [score_column, coverage_column, "candidate_id"],
        ascending=[False, False, True],
        kind="mergesort",
    )
    working = working.drop_duplicates("candidate_id", keep="first")
    working = working.drop_duplicates("guide_sequence", keep="first")
    working["selection_stratum"] = (
        working["gene_name"].fillna("").replace("", "intergenic_or_unannotated")
    )
    working["within_stratum_rank"] = working.groupby("selection_stratum").cumcount() + 1
    working = working[working["within_stratum_rank"] <= maximum_candidates_per_gene]
    if stratify_by_gene:
        working = working.sort_values(
            ["within_stratum_rank", score_column, "selection_stratum", "candidate_id"],
            ascending=[True, False, True, True],
            kind="mergesort",
        )
    return working.iloc[:maximum_total_candidates].reset_index(drop=True)


def _overlap_summary(
    deletion_start: int,
    deletion_end: int,
    features: pd.DataFrame | None,
) -> dict[str, object]:
    if features is None or features.empty:
        return {
            "overlapping_feature_count": 0,
            "overlapping_features": "",
            "fraction_of_feature_removed": pd.NA,
        }
    overlaps = features[(features["start"] <= deletion_end) & (features["end"] >= deletion_start)]
    labels: list[str] = []
    fractions: list[float] = []
    for _, feature in overlaps.iterrows():
        overlap_bp = max(
            0,
            min(deletion_end, int(feature["end"])) - max(deletion_start, int(feature["start"])) + 1,
        )
        feature_bp = int(feature["end"]) - int(feature["start"]) + 1
        fractions.append(overlap_bp / feature_bp if feature_bp else 0.0)
        label = str(feature.get("name") or feature.get("feature_id") or feature["feature_type"])
        labels.append(f"{label}:{feature['feature_type']}:{overlap_bp}bp")
    return {
        "overlapping_feature_count": len(overlaps),
        "overlapping_features": "; ".join(labels),
        "fraction_of_feature_removed": max(fractions) if fractions else 0.0,
    }


def simulate_candidate_pair(
    first: Mapping[str, object],
    second: Mapping[str, object],
    features: pd.DataFrame | None = None,
    aligned_records: Mapping[str, str] | None = None,
    reference_id: str | None = None,
    config: dict[str, Any] | str | Path | None = None,
) -> dict[str, object]:
    """Describe a deletion or multi-target computational hypothesis."""
    settings = _configured(config)
    editor = get_editor(settings)
    first_cut = cut_after_1based(first, editor)
    second_cut = cut_after_1based(second, editor)
    distance = abs(first_cut - second_cut)
    first_gene = str(first.get("gene_name", "") or "")
    second_gene = str(second.get("gene_name", "") or "")
    first_reference = str(first.get("reference_accession", reference_id or ""))
    second_reference = str(second.get("reference_accession", reference_id or ""))
    same_molecule = bool(first_reference and first_reference == second_reference)
    same_gene = bool(first_gene and first_gene == second_gene)
    maximum_deletion = int(settings["pair_selection"]["maximum_deletion_distance_bp"])
    same_region_limit = int(settings["pair_selection"]["same_region_distance_bp"])

    if same_molecule and same_gene and distance <= maximum_deletion:
        hypothesis_type = "same_gene_deletion_hypothesis"
    elif same_molecule and not first_gene and not second_gene and distance <= same_region_limit:
        hypothesis_type = "same_region_deletion_hypothesis"
    else:
        hypothesis_type = "multi_target_hypothesis"
    deletion_applicable = hypothesis_type.endswith("deletion_hypothesis")
    result: dict[str, object] = {
        "candidate_a": str(first.get("candidate_id", "A")),
        "candidate_b": str(second.get("candidate_id", "B")),
        "gene_a": first_gene,
        "gene_b": second_gene,
        "cut_coordinate_a": first_cut,
        "cut_coordinate_b": second_cut,
        "distance_bp": distance,
        "deletion_length_bp": distance if deletion_applicable else pd.NA,
        "hypothesis_type": hypothesis_type,
        "candidate_a_score": float(first.get("post_human_score", first.get("pre_human_score", 0))),
        "candidate_b_score": float(
            second.get("post_human_score", second.get("pre_human_score", 0))
        ),
        "same_reference_molecule": same_molecule,
    }
    if deletion_applicable:
        left_cut, right_cut = sorted((first_cut, second_cut))
        result.update(_overlap_summary(left_cut + 1, right_cut, features))
        coding = (
            str(first.get("feature_type", "")) == "CDS"
            or str(second.get("feature_type", "")) == "CDS"
        )
        result["coding_frame_disruption_heuristic"] = (
            "frameshift_possible_from_interval_modulo_3"
            if coding and distance % 3
            else "in_frame_or_not_applicable"
        )
        result["interpretation"] = (
            "intervening reference DNA deletion hypothesis if both predicted cuts occur"
        )
    else:
        result.update(
            {
                "overlapping_feature_count": pd.NA,
                "overlapping_features": "",
                "fraction_of_feature_removed": pd.NA,
                "coding_frame_disruption_heuristic": "not_applicable_to_multi_target_hypothesis",
                "interpretation": (
                    "two independent target hypotheses; no single intervening deletion inferred"
                ),
            }
        )
    result["limitations"] = (
        "sequence-level computational hypothesis; does not predict editing, repair, viral "
        "inactivation, delivery, toxicity, or clinical outcome"
    )
    if aligned_records is not None:
        if not reference_id:
            raise ValueError("reference_id is required when aligned_records are provided")
        coverage, count = exact_pair_coverage(first, second, aligned_records, reference_id)
        result["joint_strain_coverage"] = coverage
        result["exact_pair_genome_count"] = count
        result["genome_count"] = len(aligned_records)
    else:
        result["joint_strain_coverage"] = min(
            float(first.get("exact_strain_coverage", first.get("virus_site_coverage", 0))),
            float(second.get("exact_strain_coverage", second.get("virus_site_coverage", 0))),
        )
    return result


def rank_candidate_pairs(
    candidates: pd.DataFrame,
    features: pd.DataFrame | None = None,
    aligned_records: Mapping[str, str] | None = None,
    reference_id: str | None = None,
    same_feature_only: bool = True,
    min_distance_bp: int | None = None,
    max_distance_bp: int | None = None,
    max_candidates: int | None = None,
    config: dict[str, Any] | str | Path | None = None,
    genes: Sequence[str] | None = None,
    feature_types: Sequence[str] | None = None,
    minimum_conservation: float | None = None,
    gc_range: tuple[float, float] | None = None,
    maximum_viral_occurrence_count: int | None = None,
    maximum_candidates_per_gene: int | None = None,
    stratify_by_gene: bool | None = None,
) -> pd.DataFrame:
    """Select candidates without row-order dependence and rank pair hypotheses."""
    if candidates.empty:
        return pd.DataFrame()
    settings = _configured(config)
    selection = settings["pair_selection"]
    min_distance = int(
        selection["minimum_distance_bp"] if min_distance_bp is None else min_distance_bp
    )
    max_distance = int(
        selection["maximum_deletion_distance_bp"] if max_distance_bp is None else max_distance_bp
    )
    selected = select_pair_candidates(
        candidates,
        settings,
        genes=genes,
        feature_types=feature_types,
        minimum_conservation=minimum_conservation,
        gc_range=gc_range,
        maximum_viral_occurrence_count=maximum_viral_occurrence_count,
        maximum_candidates_per_gene=maximum_candidates_per_gene,
        maximum_total_candidates=max_candidates,
        stratify_by_gene=stratify_by_gene,
    )
    rows: list[dict[str, object]] = []
    for (_, first), (_, second) in combinations(selected.iterrows(), 2):
        same_gene = bool(first.get("gene_name")) and str(first.get("gene_name")) == str(
            second.get("gene_name")
        )
        if same_feature_only and not same_gene:
            continue
        distance = abs(
            cut_after_1based(first, get_editor(settings))
            - cut_after_1based(second, get_editor(settings))
        )
        if same_gene and not min_distance <= distance <= max_distance:
            continue
        if not same_gene and distance < min_distance:
            continue
        rows.append(
            simulate_candidate_pair(
                first,
                second,
                features=features,
                aligned_records=aligned_records,
                reference_id=reference_id,
                config=settings,
            )
        )
    if not rows:
        return pd.DataFrame()
    pairs = pd.DataFrame(rows)
    weights = {key: float(value) for key, value in settings["pair_scoring"]["weights"].items()}
    individual = (pairs["candidate_a_score"] + pairs["candidate_b_score"]) / 2
    fraction = pd.to_numeric(pairs["fraction_of_feature_removed"], errors="coerce").fillna(0)
    saturation = float(settings["pair_scoring"]["distance_saturation_bp"])
    distance_component = (
        pd.to_numeric(pairs["deletion_length_bp"], errors="coerce").fillna(0) / saturation
    ).clip(upper=1.0)
    pairs["pair_score"] = (
        weights["joint_strain_coverage"] * pairs["joint_strain_coverage"].astype(float)
        + weights["individual_candidate_scores"] * individual
        + weights["feature_fraction_removed"] * fraction
        + weights["distance_component"] * distance_component
    ) / sum(weights.values())
    pairs["pair_score_components"] = [
        json.dumps(
            {
                "joint_strain_coverage": round(float(coverage), 6),
                "individual_candidate_scores": round(float(candidate_score), 6),
                "feature_fraction_removed": round(float(feature_score), 6),
                "distance_component": round(float(distance_score), 6),
            },
            sort_keys=True,
        )
        for coverage, candidate_score, feature_score, distance_score in zip(
            pairs["joint_strain_coverage"], individual, fraction, distance_component, strict=True
        )
    ]
    # Compatibility aliases retained for existing downstream notebooks.
    pairs["cut_a_after_1based"] = pairs["cut_coordinate_a"]
    pairs["cut_b_after_1based"] = pairs["cut_coordinate_b"]
    pairs["max_feature_fraction_removed"] = pairs["fraction_of_feature_removed"]
    pairs["sequence_disruption_score"] = pairs["pair_score"]
    return pairs.sort_values(
        ["pair_score", "joint_strain_coverage", "candidate_a", "candidate_b"],
        ascending=[False, False, True, True],
        kind="mergesort",
    ).reset_index(drop=True)
