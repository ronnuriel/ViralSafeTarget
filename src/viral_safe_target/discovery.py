"""Genome-wide candidate mapping, balanced selection, and gene-level ranking."""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .annotations import read_gff3
from .config import get_editor, load_config
from .disruption import cut_after_1based
from .provenance import sha256_file

FEATURE_MAP_COLUMNS = [
    "candidate_id",
    "guide_sequence",
    "reference_accession",
    "candidate_start",
    "candidate_end",
    "strand",
    "cut_position",
    "feature_id",
    "feature_name",
    "gene_name",
    "feature_type",
    "feature_start",
    "feature_end",
    "overlap_bp",
    "overlap_fraction",
    "mapping_status",
    "mapping_explanation",
    "annotation_source",
    "annotation_version",
]


@dataclass(frozen=True)
class DiscoverySelection:
    panel: pd.DataFrame
    audit: pd.DataFrame
    genes_without_candidates: pd.DataFrame


def build_candidate_feature_map(
    candidates: pd.DataFrame,
    features: pd.DataFrame | str | Path,
    *,
    config: dict[str, Any] | str | Path | None = None,
    annotation_source: str | Path = "",
) -> pd.DataFrame:
    """Return every inclusive 1-based candidate/feature overlap, retaining intergenic rows."""
    feature_frame = read_gff3(features) if isinstance(features, (str, Path)) else features.copy()
    settings = config if isinstance(config, dict) else load_config(config)
    editor = get_editor(settings)
    source = Path(annotation_source or features) if isinstance(features, (str, Path)) else None
    version = sha256_file(source)[:12] if source and source.is_file() else "unversioned"
    rows: list[dict[str, object]] = []
    for _, candidate in candidates.iterrows():
        start = int(candidate["reference_start_1based"])
        end = int(candidate["reference_end_1based"])
        strand = str(candidate["strand"])
        if start < 1 or end < start or strand not in {"+", "-"}:
            raise ValueError(
                f"Invalid 1-based inclusive candidate coordinates for {candidate['candidate_id']}: "
                f"{start}-{end} ({strand})"
            )
        relevant = feature_frame
        accession = str(candidate.get("reference_accession", ""))
        if accession and "seqid" in relevant:
            relevant = relevant[relevant["seqid"].astype(str).eq(accession)]
        overlaps = relevant[(relevant["start"] <= end) & (relevant["end"] >= start)]
        common = {
            "candidate_id": candidate["candidate_id"],
            "guide_sequence": candidate["guide_sequence"],
            "reference_accession": accession,
            "candidate_start": start,
            "candidate_end": end,
            "strand": strand,
            "cut_position": cut_after_1based(candidate, editor),
            "annotation_source": str(source.resolve()) if source and source.is_file() else "",
            "annotation_version": version,
        }
        if overlaps.empty:
            rows.append(
                {
                    **common,
                    "feature_id": "",
                    "feature_name": "",
                    "gene_name": "",
                    "feature_type": "intergenic_or_unannotated",
                    "feature_start": pd.NA,
                    "feature_end": pd.NA,
                    "overlap_bp": 0,
                    "overlap_fraction": 0.0,
                    "mapping_status": "intergenic",
                    "mapping_explanation": (
                        "No annotated feature overlaps this inclusive guide interval."
                    ),
                }
            )
            continue
        for _, feature in overlaps.sort_values(
            ["start", "end", "feature_type", "feature_id"], kind="mergesort"
        ).iterrows():
            overlap = min(end, int(feature["end"])) - max(start, int(feature["start"])) + 1
            name = str(feature.get("name", "") or "")
            status = "mapped" if feature.get("feature_id") else "mapped_incomplete_annotation"
            rows.append(
                {
                    **common,
                    "feature_id": feature.get("feature_id", ""),
                    "feature_name": name,
                    "gene_name": name,
                    "feature_type": feature["feature_type"],
                    "feature_start": int(feature["start"]),
                    "feature_end": int(feature["end"]),
                    "overlap_bp": overlap,
                    "overlap_fraction": overlap / (end - start + 1),
                    "mapping_status": status,
                    "mapping_explanation": (
                        "Inclusive 1-based candidate interval overlaps this feature; "
                        "one-to-many mappings are retained."
                    ),
                }
            )
    return pd.DataFrame(rows, columns=FEATURE_MAP_COLUMNS)


def _ranked_candidates(candidates: pd.DataFrame) -> pd.DataFrame:
    output = candidates.copy()
    if "rejection_reasons" in output:
        output = output[output["rejection_reasons"].fillna("").eq("")]
    output = output.sort_values(
        ["pre_human_score", "candidate_id"], ascending=[False, True], kind="mergesort"
    ).drop_duplicates("candidate_id")
    output["pre_human_rank"] = range(1, len(output) + 1)
    return output


def select_balanced_discovery_panel(
    candidates: pd.DataFrame,
    feature_map: pd.DataFrame,
    annotated_features: pd.DataFrame,
    *,
    top_per_gene: int = 50,
    global_top: int = 500,
    exhaustive: bool = False,
    confirm_exhaustive: bool = False,
) -> DiscoverySelection:
    """Select a deterministic pre-human panel without consulting human-screen results."""
    eligible = _ranked_candidates(candidates)
    if exhaustive and not confirm_exhaustive:
        raise ValueError(
            f"Exhaustive mode would include {len(eligible):,} candidates. "
            "Pass --confirm-exhaustive explicitly."
        )
    gene_features = annotated_features[annotated_features["feature_type"].eq("gene")].copy()
    annotated_genes = sorted(gene_features["name"].fillna("").replace("", pd.NA).dropna().unique())
    memberships = feature_map[
        feature_map["gene_name"].isin(annotated_genes)
        & feature_map["feature_type"].isin(["gene", "CDS", "ncRNA"])
    ][["candidate_id", "gene_name"]].drop_duplicates()
    audit: dict[str, dict[str, object]] = {}

    def record(candidate_id: str, reason: str, gene: str = "") -> None:
        item = audit.setdefault(
            candidate_id,
            {"candidate_id": candidate_id, "per_gene_quota_genes": set(), "global_top": False},
        )
        if reason == "per_gene_quota":
            item["per_gene_quota_genes"].add(gene)
        else:
            item["global_top"] = True

    if exhaustive:
        for candidate_id in eligible["candidate_id"]:
            record(str(candidate_id), "global_top")
    else:
        ranked_lookup = eligible.set_index("candidate_id")
        for gene in annotated_genes:
            ids = memberships.loc[memberships["gene_name"].eq(gene), "candidate_id"]
            gene_candidates = ranked_lookup.loc[ranked_lookup.index.intersection(ids)].sort_values(
                ["pre_human_score", "candidate_id"],
                ascending=[False, True],
                kind="mergesort",
            )
            for candidate_id in gene_candidates.head(top_per_gene).index:
                record(str(candidate_id), "per_gene_quota", gene)
        for candidate_id in eligible.head(global_top)["candidate_id"]:
            record(str(candidate_id), "global_top")

    audit_rows = []
    for candidate_id, item in sorted(audit.items()):
        genes = sorted(item["per_gene_quota_genes"])
        global_entry = bool(item["global_top"])
        reason = "both" if genes and global_entry else "per_gene_quota" if genes else "global_top"
        audit_rows.append(
            {
                "candidate_id": candidate_id,
                "selection_reason": reason if not exhaustive else "exhaustive",
                "per_gene_quota_genes": ";".join(genes),
                "global_top": global_entry,
                "top_per_gene": top_per_gene,
                "global_top_limit": global_top,
            }
        )
    audit_frame = pd.DataFrame(audit_rows)
    panel = eligible[eligible["candidate_id"].isin(audit)].merge(
        audit_frame, on="candidate_id", how="left", validate="one_to_one"
    )
    aggregated = (
        feature_map.groupby("candidate_id", sort=True)
        .agg(
            mapped_feature_ids=(
                "feature_id",
                lambda values: ";".join(sorted(set(filter(None, values.astype(str))))),
            ),
            mapped_gene_names=(
                "gene_name",
                lambda values: ";".join(sorted(set(filter(None, values.astype(str))))),
            ),
            mapped_feature_types=(
                "feature_type",
                lambda values: ";".join(sorted(set(values.astype(str)))),
            ),
            feature_mapping_count=("feature_id", "size"),
        )
        .reset_index()
    )
    panel = panel.merge(aggregated, on="candidate_id", how="left", validate="one_to_one")
    panel = panel.sort_values(["pre_human_rank", "candidate_id"], kind="mergesort").reset_index(
        drop=True
    )
    genes_with_candidates = set(memberships["gene_name"])
    no_candidates = gene_features[~gene_features["name"].isin(genes_with_candidates)][
        ["feature_id", "name", "feature_type", "start", "end"]
    ].rename(columns={"name": "gene_name", "start": "feature_start", "end": "feature_end"})
    no_candidates["reason"] = "no eligible pre-human candidate overlaps this annotated gene"
    return DiscoverySelection(panel, audit_frame, no_candidates.reset_index(drop=True))


def wilson_lower_bound(successes: int, total: int, z: float = 1.959963984540054) -> float:
    if total <= 0:
        return math.nan
    proportion = successes / total
    denominator = 1 + z**2 / total
    centre = proportion + z**2 / (2 * total)
    margin = z * math.sqrt((proportion * (1 - proportion) + z**2 / (4 * total)) / total)
    return (centre - margin) / denominator


def _desirability_percentile(values: pd.Series, *, higher_is_better: bool) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    count = int(numeric.notna().sum())
    ranks = numeric.rank(method="average", ascending=not higher_is_better, na_option="keep")
    if count == 0:
        return pd.Series(float("nan"), index=values.index)
    if count == 1:
        return numeric.notna().astype(float).replace(0.0, float("nan"))
    return 1 - (ranks - 1) / (count - 1)


def rank_genes(
    candidates: pd.DataFrame,
    feature_map: pd.DataFrame,
    annotated_features: pd.DataFrame,
    *,
    eligible_candidates: pd.DataFrame | None = None,
    evidence: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Compute separate targetability and evidence metrics for every annotated gene feature."""
    eligible_pool = candidates if eligible_candidates is None else eligible_candidates
    gene_features = annotated_features[annotated_features["feature_type"].eq("gene")].copy()
    rows = []
    for _, gene in gene_features.sort_values(["start", "end", "feature_id"]).iterrows():
        gene_name = str(gene["name"])
        feature_id = str(gene["feature_id"])
        member_ids = feature_map.loc[
            feature_map["feature_id"].eq(feature_id), "candidate_id"
        ].drop_duplicates()
        gene_candidates = eligible_pool[eligible_pool["candidate_id"].isin(member_ids)].copy()
        screened = candidates[
            candidates["candidate_id"].isin(member_ids)
            & candidates["screening_status"].eq("completed")
        ]
        ranked = screened[screened["post_human_rank"].notna()].sort_values(
            ["post_human_rank", "candidate_id"], kind="mergesort"
        )
        clean = screened[
            pd.to_numeric(screened["human_total_predicted_hits"], errors="coerce").eq(0)
        ]
        eligible_count = len(gene_candidates)
        screened_count = len(screened)
        length = int(gene["end"]) - int(gene["start"]) + 1
        evidence_rows = (
            evidence[evidence["gene_name"].eq(gene_name)]
            if evidence is not None
            else pd.DataFrame()
        )
        top5 = ranked.head(5)
        top10 = ranked.head(10)
        rows.append(
            {
                "gene_name": gene_name,
                "feature_id": feature_id,
                "feature_type": "gene",
                "gene_length_bp": length,
                "eligible_candidate_count": eligible_count,
                "screened_candidate_count": screened_count,
                "screening_fraction": screened_count / eligible_count
                if eligible_count
                else math.nan,
                "clean_candidate_count": len(clean),
                "expert_review_candidate_count": int(
                    screened["decision"].astype(str).eq("expert_review_required").sum()
                ),
                "failed_or_missing_count": eligible_count - screened_count,
                "clean_fraction": len(clean) / screened_count if screened_count else math.nan,
                "clean_fraction_wilson_lower_95": wilson_lower_bound(len(clean), screened_count),
                "candidate_density_per_kb": eligible_count / (length / 1000),
                "conserved_candidate_count": int(
                    pd.to_numeric(gene_candidates["exact_strain_coverage"], errors="coerce")
                    .ge(0.95)
                    .sum()
                ),
                "conserved_candidate_fraction": (
                    pd.to_numeric(gene_candidates["exact_strain_coverage"], errors="coerce")
                    .ge(0.95)
                    .mean()
                    if eligible_count
                    else math.nan
                ),
                "exact_all_accepted_genomes_count": int(
                    gene_candidates["exact_genome_count"].eq(gene_candidates["genome_count"]).sum()
                ),
                "best_candidate_id": ranked.iloc[0]["candidate_id"] if not ranked.empty else pd.NA,
                "best_candidate_rank": ranked.iloc[0]["post_human_rank"]
                if not ranked.empty
                else pd.NA,
                "best_candidate_score": ranked.iloc[0]["post_human_score"]
                if not ranked.empty
                else pd.NA,
                "median_top_5_candidate_rank": top5["post_human_rank"].median()
                if not top5.empty
                else pd.NA,
                "median_top_10_candidate_rank": top10["post_human_rank"].median()
                if not top10.empty
                else pd.NA,
                "worst_top_5_rank": top5["post_human_rank"].max() if not top5.empty else pd.NA,
                "median_human_hit_count": pd.to_numeric(
                    screened["human_total_predicted_hits"], errors="coerce"
                ).median(),
                "one_mismatch_candidate_count": int(
                    screened["human_one_mismatch_hit_count"].gt(0).sum()
                ),
                "two_mismatch_candidate_count": int(
                    screened["human_two_mismatch_hit_count"].gt(0).sum()
                ),
                "three_mismatch_candidate_count": int(
                    screened["human_three_mismatch_hit_count"].gt(0).sum()
                ),
                "tool_coverage": screened_count / eligible_count if eligible_count else math.nan,
                "rank_variance": pd.to_numeric(top10["post_human_rank"], errors="coerce").var(
                    ddof=0
                ),
                "evidence_coverage": (1.0 if not evidence_rows.empty else 0.0)
                if evidence is not None
                else pd.NA,
                "biological_evidence_status": "supplied"
                if not evidence_rows.empty
                else "biological evidence not supplied",
                "confidence_level": (
                    "high"
                    if screened_count >= 20
                    else "moderate"
                    if screened_count >= 10
                    else "low"
                ),
            }
        )
    output = pd.DataFrame(rows)
    components = {
        "best_candidate_component": ("best_candidate_rank", False),
        "top5_robustness_component": ("median_top_5_candidate_rank", False),
        "clean_fraction_component": ("clean_fraction_wilson_lower_95", True),
        "human_burden_component": ("median_human_hit_count", False),
        # A fraction, not a raw count, prevents long genes from receiving an
        # automatic advantage merely because they contain more guide sites.
        "conservation_component": ("conserved_candidate_fraction", True),
    }
    for component, (source, higher) in components.items():
        output[component] = _desirability_percentile(output[source], higher_is_better=higher)
    component_columns = list(components)
    output["targetability_score"] = output[component_columns].mean(axis=1, skipna=True)
    output.loc[output[component_columns].notna().sum(axis=1).lt(2), "targetability_score"] = pd.NA
    support_factor = output["screened_candidate_count"].clip(upper=10) / 10
    output["targetability_score"] = output["targetability_score"] * support_factor
    output["best_single_candidate_rank"] = output["best_candidate_rank"].rank(
        method="min", ascending=True, na_option="keep"
    )
    output["top5_robustness_rank"] = output["median_top_5_candidate_rank"].rank(
        method="min", ascending=True, na_option="keep"
    )
    output["clean_fraction_rank"] = output["clean_fraction_wilson_lower_95"].rank(
        method="min", ascending=False, na_option="keep"
    )
    output["targetability_rank"] = output["targetability_score"].rank(
        method="min", ascending=False, na_option="keep"
    )
    output["targetability_explanation"] = output.apply(
        lambda row: (
            "Rank-aggregated computational targetability from within-view percentiles; "
            f"screened={int(row['screened_candidate_count'])}, "
            f"eligible={int(row['eligible_candidate_count'])}, "
            f"confidence={row['confidence_level']}. "
            "Biological evidence is reported separately and is not added to this score."
        ),
        axis=1,
    )
    return output.sort_values(["targetability_rank", "gene_name", "feature_id"], kind="mergesort")


def gene_rank_stability(
    screened_candidates: pd.DataFrame,
    feature_map: pd.DataFrame,
    annotated_features: pd.DataFrame,
    quotas: Sequence[int] = (10, 25, 50),
) -> pd.DataFrame:
    """Recompute nested gene rankings from one completed maximum-quota screen."""
    rankings: dict[int, pd.Series] = {}
    for quota in quotas:
        allowed_ids: set[str] = set()
        for _, group in feature_map[
            feature_map["feature_type"].isin(["gene", "CDS", "ncRNA"])
            & feature_map["gene_name"].ne("")
        ].groupby("gene_name"):
            ids = set(group["candidate_id"])
            subset = screened_candidates[screened_candidates["candidate_id"].isin(ids)].sort_values(
                ["pre_human_score", "candidate_id"], ascending=[False, True], kind="mergesort"
            )
            allowed_ids.update(subset.head(quota)["candidate_id"].astype(str))
        limited = screened_candidates[screened_candidates["candidate_id"].isin(allowed_ids)]
        ranked = rank_genes(limited, feature_map, annotated_features)
        rankings[quota] = ranked.groupby("gene_name")["targetability_rank"].min()
    genes = sorted(set().union(*(set(series.index) for series in rankings.values())))
    rows = []
    for gene in genes:
        values = [float(rankings[quota].get(gene, math.nan)) for quota in quotas]
        finite = [value for value in values if not math.isnan(value)]
        row = {"gene_name": gene}
        for quota, value in zip(quotas, values, strict=True):
            row[f"rank_at_k{quota}"] = value
        row.update(
            {
                "mean_rank": np.mean(finite) if finite else math.nan,
                "rank_range": max(finite) - min(finite) if finite else math.nan,
                "rank_std": np.std(finite) if finite else math.nan,
                "top_5_stability": bool(finite and all(value <= 5 for value in finite)),
                "top_10_stability": bool(finite and all(value <= 10 for value in finite)),
                "stability_warning": (
                    "quota-sensitive ranking" if finite and max(finite) - min(finite) > 10 else ""
                ),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["mean_rank", "gene_name"], kind="mergesort")


def build_deep_screening_panel(
    candidates: pd.DataFrame,
    feature_map: pd.DataFrame,
    gene_rankings: pd.DataFrame,
    *,
    top_genes: int = 10,
    top_per_gene: int = 20,
    global_top: int = 100,
) -> pd.DataFrame:
    ranked_candidates = candidates[candidates["post_human_rank"].notna()].sort_values(
        ["post_human_rank", "candidate_id"], kind="mergesort"
    )
    selected: dict[str, set[str]] = {}
    genes = gene_rankings.head(top_genes)["gene_name"].astype(str)
    for gene in genes:
        ids = set(feature_map.loc[feature_map["gene_name"].eq(gene), "candidate_id"])
        for candidate_id in ranked_candidates[ranked_candidates["candidate_id"].isin(ids)].head(
            top_per_gene
        )["candidate_id"]:
            selected.setdefault(str(candidate_id), set()).add(f"top_gene:{gene}")
    for candidate_id in ranked_candidates.head(global_top)["candidate_id"]:
        selected.setdefault(str(candidate_id), set()).add("global_top")
    output = ranked_candidates[ranked_candidates["candidate_id"].isin(selected)].copy()
    output["deep_panel_reason"] = output["candidate_id"].map(
        lambda value: ";".join(sorted(selected[str(value)]))
    )
    return output.sort_values(["post_human_rank", "candidate_id"], kind="mergesort")


def build_bounded_pair_hypotheses(
    candidates: pd.DataFrame,
    feature_map: pd.DataFrame,
    gene_rankings: pd.DataFrame,
    *,
    top_genes: int = 10,
    candidates_per_gene: int = 10,
    maximum_pairs: int = 5000,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build bounded same-gene and cross-gene hypotheses after candidate ranking."""
    from itertools import combinations

    top_gene_names = gene_rankings.head(top_genes)["gene_name"].astype(str).tolist()
    candidate_lookup = candidates.set_index("candidate_id")
    per_gene: dict[str, list[str]] = {}
    for gene in top_gene_names:
        ids = set(feature_map.loc[feature_map["gene_name"].eq(gene), "candidate_id"])
        ranked = candidates[candidates["candidate_id"].isin(ids)].sort_values(
            ["post_human_rank", "candidate_id"], kind="mergesort"
        )
        per_gene[gene] = ranked.head(candidates_per_gene)["candidate_id"].astype(str).tolist()

    def coverage(first: pd.Series, second: pd.Series) -> float:
        left = set(str(first.get("exact_site_accessions", "")).split(";")) - {""}
        right = set(str(second.get("exact_site_accessions", "")).split(";")) - {""}
        denominator = max(int(first.get("genome_count", 0)), int(second.get("genome_count", 0)), 1)
        return len(left & right) / denominator

    same_rows = []
    for gene, ids in per_gene.items():
        for left_id, right_id in combinations(ids, 2):
            left, right = candidate_lookup.loc[left_id], candidate_lookup.loc[right_id]
            cut_left, cut_right = int(left["cut_position"]), int(right["cut_position"])
            joint = coverage(left, right)
            same_rows.append(
                {
                    "candidate_a": left_id,
                    "candidate_b": right_id,
                    "gene_a": gene,
                    "gene_b": gene,
                    "distance_bp": abs(cut_right - cut_left),
                    "deletion_length_bp": abs(cut_right - cut_left),
                    "joint_strain_coverage": joint,
                    "pair_score": np.mean(
                        [left["post_human_score"], right["post_human_score"], joint]
                    ),
                    "hypothesis_type": "same_gene_deletion_hypothesis",
                    "limitations": (
                        "theoretical cut-to-cut interval; requires independent and "
                        "experimental validation"
                    ),
                }
            )
    multi_rows = []
    for left_gene, right_gene in combinations(top_gene_names, 2):
        for left_id in per_gene[left_gene][:3]:
            for right_id in per_gene[right_gene][:3]:
                left, right = candidate_lookup.loc[left_id], candidate_lookup.loc[right_id]
                joint = coverage(left, right)
                multi_rows.append(
                    {
                        "candidate_a": left_id,
                        "candidate_b": right_id,
                        "gene_a": left_gene,
                        "gene_b": right_gene,
                        "distance_bp": pd.NA,
                        "deletion_length_bp": pd.NA,
                        "joint_strain_coverage": joint,
                        "pair_score": np.mean(
                            [left["post_human_score"], right["post_human_score"], joint]
                        ),
                        "hypothesis_type": "multi_target_hypothesis",
                        "limitations": (
                            "two separate target sites; not one physical deletion; "
                            "requires independent validation"
                        ),
                    }
                )
    columns = [
        "candidate_a",
        "candidate_b",
        "gene_a",
        "gene_b",
        "distance_bp",
        "deletion_length_bp",
        "joint_strain_coverage",
        "pair_score",
        "hypothesis_type",
        "limitations",
    ]
    same = (
        pd.DataFrame(same_rows, columns=columns)
        .sort_values(["pair_score", "candidate_a", "candidate_b"], ascending=[False, True, True])
        .head(maximum_pairs)
    )
    multi = (
        pd.DataFrame(multi_rows, columns=columns)
        .sort_values(["pair_score", "candidate_a", "candidate_b"], ascending=[False, True, True])
        .head(maximum_pairs)
    )
    summary = pd.DataFrame(columns=["gene_name", "pair_count", "best_pair_score"])
    if not same.empty:
        summary = (
            same.groupby("gene_a")
            .agg(pair_count=("candidate_a", "count"), best_pair_score=("pair_score", "max"))
            .reset_index()
            .rename(columns={"gene_a": "gene_name"})
        )
    return same, multi, summary
