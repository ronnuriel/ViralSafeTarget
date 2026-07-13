from __future__ import annotations

import pandas as pd


def _human_safety(row: pd.Series) -> float:
    exact = row.get("host_exact_matches")
    minimum = row.get("host_min_mismatches")
    if pd.notna(exact) and float(exact) > 0:
        return 0.0
    if pd.isna(minimum):
        return 1.0
    # Demonstration-only monotonic transform; not a validated biological risk model.
    return min(1.0, max(0.0, float(minimum) / 4.0))


def rank_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    """Create an explainable DEMO ranking. Keep components visible; do not treat as validated."""
    if candidates.empty:
        return candidates.copy()
    ranked = candidates.copy()
    ranked["annotation_component"] = (
        ranked.get("feature_type", "intergenic_or_unannotated") != "intergenic_or_unannotated"
    ).astype(float)
    ranked["human_safety_component"] = ranked.apply(_human_safety, axis=1)
    ranked["demo_score"] = (
        0.60 * ranked["virus_site_coverage"].astype(float)
        + 0.15 * ranked["annotation_component"]
        + 0.25 * ranked["human_safety_component"]
    )
    ranked["decision"] = "review"
    if "host_exact_matches" in ranked:
        ranked.loc[
            ranked["host_exact_matches"].fillna(0) > 0,
            "decision",
        ] = "reject_demo_exact_host_match"
    ranked.loc[
        (ranked["virus_site_coverage"] >= 0.95)
        & (ranked["human_safety_component"] >= 0.75)
        & (ranked["annotation_component"] > 0),
        "decision",
    ] = "prioritize_for_literature_review"
    return ranked.sort_values(
        ["demo_score", "virus_site_coverage"], ascending=[False, False]
    ).reset_index(drop=True)
