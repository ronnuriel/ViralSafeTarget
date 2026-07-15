"""Explainable pre-human and post-human candidate ranking."""

from __future__ import annotations

import math
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd

from .config import load_config
from .evidence import merge_gene_evidence, read_gene_evidence


def gc_fraction(sequence: str) -> float:
    sequence = sequence.upper()
    valid = [base for base in sequence if base in "ACGT"]
    return (sum(base in "GC" for base in valid) / len(valid)) if valid else 0.0


def gc_score(value: float, settings: dict[str, float]) -> float:
    """Piecewise conservative GC score configured entirely in YAML."""
    low = float(settings["accepted_min"])
    preferred_low = float(settings["preferred_min"])
    preferred_high = float(settings["preferred_max"])
    high = float(settings["accepted_max"])
    if value < low or value > high:
        return 0.0
    if preferred_low <= value <= preferred_high:
        return 1.0
    if value < preferred_low:
        return (value - low) / (preferred_low - low)
    return (high - value) / (high - preferred_high)


def sequence_complexity_score(sequence: str) -> float:
    """Normalized mono-nucleotide Shannon entropy in the range 0..1."""
    sequence = sequence.upper()
    valid = [base for base in sequence if base in "ACGT"]
    if not valid:
        return 0.0
    counts = Counter(valid)
    entropy = -sum(
        (count / len(valid)) * math.log2(count / len(valid)) for count in counts.values()
    )
    return entropy / 2.0


def _longest_homopolymer(sequence: str) -> int:
    longest = current = 0
    previous = ""
    for base in sequence.upper():
        current = current + 1 if base == previous else 1
        longest = max(longest, current)
        previous = base
    return longest


def _configured(config: dict[str, Any] | str | Path | None) -> dict[str, Any]:
    return config if isinstance(config, dict) else load_config(config)


def rank_pre_human_candidates(
    candidates: pd.DataFrame,
    config: dict[str, Any] | str | Path | None = None,
    evidence_path: str | Path | None = None,
) -> pd.DataFrame:
    """Calculate visible, deterministic sequence and evidence components."""
    if candidates.empty:
        return candidates.copy()
    settings = _configured(config)
    ranking = settings["ranking"]
    output = merge_gene_evidence(candidates, read_gene_evidence(evidence_path))

    output["exact_strain_coverage"] = output["virus_site_coverage"].astype(float)
    output["conservation_score"] = output["exact_strain_coverage"].clip(0.0, 1.0)
    occurrences = (
        output.get("reference_viral_occurrence_count", pd.Series(1, index=output.index))
        .fillna(1)
        .astype(float)
        .clip(lower=1)
    )
    output["viral_uniqueness_score"] = 1.0 / occurrences
    output["gc_fraction"] = output["guide_sequence"].astype(str).map(gc_fraction)
    output["gc_score"] = output["gc_fraction"].map(lambda value: gc_score(value, ranking["gc"]))
    output["sequence_complexity_score"] = (
        output["guide_sequence"].astype(str).map(sequence_complexity_score)
    )
    minimum_complexity = float(ranking["complexity"]["minimum_score"])
    output["low_complexity_penalty"] = (
        (minimum_complexity - output["sequence_complexity_score"]) / minimum_complexity
    ).clip(lower=0.0, upper=1.0)
    maximum_homopolymer = int(ranking["complexity"]["maximum_homopolymer"])
    runs = output["guide_sequence"].astype(str).map(_longest_homopolymer)
    output["homopolymer_length"] = runs
    output["homopolymer_penalty"] = (
        (runs - maximum_homopolymer).clip(lower=0) / output["guide_sequence"].str.len()
    ).clip(upper=1.0)
    ambiguous = (
        output["guide_sequence"]
        .astype(str)
        .map(lambda sequence: sum(base not in "ACGT" for base in sequence.upper()))
    )
    output["ambiguous_base_count"] = ambiguous
    output["ambiguous_base_penalty"] = ambiguous / output["guide_sequence"].str.len()
    feature_types = output.get(
        "feature_type", pd.Series("intergenic_or_unannotated", index=output.index)
    ).fillna("intergenic_or_unannotated")
    output["feature_type"] = feature_types
    output["annotation_score"] = feature_types.map(ranking["annotation_scores"]).fillna(0.0)
    evidence_levels = output.get("evidence_level", pd.Series(pd.NA, index=output.index))
    evidence_scores = ranking["evidence_scores"]
    output["gene_evidence_score"] = evidence_levels.map(evidence_scores).astype("Float64")

    weights = {key: float(value) for key, value in ranking["weights"].items()}
    components = {
        "conservation": output["conservation_score"],
        "viral_uniqueness": output["viral_uniqueness_score"],
        "gc": output["gc_score"],
        "sequence_complexity": output["sequence_complexity_score"],
        "annotation": output["annotation_score"],
        "gene_evidence": output["gene_evidence_score"],
    }
    denominator = sum(weights.values())
    score = pd.Series(0.0, index=output.index)
    available_weight = pd.Series(0.0, index=output.index)
    for name, component in components.items():
        numeric = pd.to_numeric(component, errors="coerce")
        score += numeric.fillna(0.0) * weights[name]
        available_weight += numeric.notna().astype(float) * weights[name]
    penalties = (
        output["low_complexity_penalty"]
        + output["homopolymer_penalty"]
        + output["ambiguous_base_penalty"]
    ).clip(upper=1.0)
    output["scoring_weight_coverage"] = available_weight / denominator
    output["pre_human_score"] = ((score / denominator) * (1.0 - penalties)).clip(0.0, 1.0)

    rejection = ranking["rejection"]
    gc_settings = ranking["gc"]

    def rejection_reasons(row: pd.Series) -> str:
        reasons: list[str] = []
        if float(row["exact_strain_coverage"]) < float(rejection["minimum_conservation"]):
            reasons.append("conservation_below_configured_minimum")
        if float(row["gc_fraction"]) < float(gc_settings["accepted_min"]):
            reasons.append("gc_below_configured_range")
        if float(row["gc_fraction"]) > float(gc_settings["accepted_max"]):
            reasons.append("gc_above_configured_range")
        if float(row["sequence_complexity_score"]) < minimum_complexity:
            reasons.append("sequence_complexity_below_configured_minimum")
        if int(row["ambiguous_base_count"]) > int(rejection["maximum_ambiguous_bases"]):
            reasons.append("ambiguous_bases_exceed_configured_maximum")
        return ";".join(reasons)

    output["rejection_reasons"] = output.apply(rejection_reasons, axis=1)
    output["pre_human_decision"] = output["rejection_reasons"].map(
        lambda reasons: "reject_pre_human" if reasons else "retain_for_human_off_target_screen"
    )
    output["decision"] = output["pre_human_decision"]

    def explanation(row: pd.Series) -> str:
        evidence_text = (
            "curated gene evidence unavailable; no positive evidence contribution"
            if pd.isna(row["gene_evidence_score"])
            else f"curated gene evidence component={float(row['gene_evidence_score']):.3f}"
        )
        return (
            f"conservation={row['conservation_score']:.3f}; "
            f"reference occurrences={int(row.get('reference_viral_occurrence_count', 1))}; "
            f"GC={row['gc_fraction']:.3f} (component={row['gc_score']:.3f}); "
            f"complexity={row['sequence_complexity_score']:.3f}; "
            f"annotation={row['feature_type']} (component={row['annotation_score']:.3f}); "
            f"{evidence_text}; penalties={float(penalties.loc[row.name]):.3f}"
        )

    output["rank_explanation"] = output.apply(explanation, axis=1)
    return output.sort_values(
        ["pre_human_score", "exact_strain_coverage", "candidate_id"],
        ascending=[False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)


def rank_post_human_candidates(
    candidates: pd.DataFrame,
    config: dict[str, Any] | str | Path | None = None,
) -> pd.DataFrame:
    """Apply a separate, transparent post-human predicted-risk layer."""
    if candidates.empty:
        return candidates.copy()
    settings = _configured(config)
    output = candidates.copy()
    risk_weights = settings["off_target"]["risk_weights"]
    hit_columns = {
        "exact": "human_exact_hit_count",
        "one_mismatch": "human_one_mismatch_hit_count",
        "two_mismatch": "human_two_mismatch_hit_count",
        "three_mismatch": "human_three_mismatch_hit_count",
    }
    risk = pd.Series(0.0, index=output.index)
    for label, column in hit_columns.items():
        counts = pd.to_numeric(output.get(column, 0), errors="coerce").fillna(0)
        risk += float(risk_weights[label]) * counts.clip(upper=1)
    output["predicted_offtarget_risk"] = risk.clip(0.0, 1.0)
    output["post_human_score"] = (
        output["pre_human_score"].astype(float) * (1.0 - output["predicted_offtarget_risk"])
    ).clip(0.0, 1.0)

    threshold = int(settings["editor"]["mismatch_search_threshold"])

    def decision(row: pd.Series) -> tuple[str, str]:
        rejection_reasons = row.get("rejection_reasons", "")
        if pd.notna(rejection_reasons) and str(rejection_reasons).strip():
            return "exclude_pre_human", str(rejection_reasons)
        if int(row.get("human_exact_hit_count", 0)) > 0:
            return (
                "exclude_or_expert_review",
                "predicted exact human off-target requires exclusion or expert review",
            )
        if int(row.get("human_total_predicted_hits", 0)) > 0:
            return (
                "expert_review_required",
                "predicted human off-target hit within configured mismatch threshold",
            )
        assembly = settings["off_target"]["human_assembly"]
        reason = (
            f"no predicted human hit up to {threshold} mismatches in {assembly} "
            "under this editor model"
        )
        return (
            "retain_computational_candidate",
            reason,
        )

    decisions = output.apply(decision, axis=1)
    output["decision"] = decisions.map(lambda value: value[0])
    output["decision_reason"] = decisions.map(lambda value: value[1])
    return output.sort_values(
        ["post_human_score", "pre_human_score", "candidate_id"],
        ascending=[False, False, True],
        kind="mergesort",
    ).reset_index(drop=True)


def rank_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    """Backward-compatible ranking entry point used by the demo and UI."""
    ranked = rank_pre_human_candidates(candidates)
    if "human_exact_hit_count" in ranked or "host_exact_matches" in ranked:
        if "human_exact_hit_count" not in ranked:
            ranked["human_exact_hit_count"] = ranked["host_exact_matches"]
            ranked["human_total_predicted_hits"] = ranked.filter(
                regex=r"host_matches_le_\d+_mismatches"
            ).max(axis=1)
            ranked["human_one_mismatch_hit_count"] = 0
            ranked["human_two_mismatch_hit_count"] = 0
            ranked["human_three_mismatch_hit_count"] = 0
        ranked = rank_post_human_candidates(ranked)
    ranked["demo_score"] = ranked.get("post_human_score", ranked["pre_human_score"])
    return ranked
