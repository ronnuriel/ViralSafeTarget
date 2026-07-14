"""Multi-objective shortlist and comparison-set construction for research showcases."""

# Research finding prose is intentionally readable as complete strings.
# ruff: noqa: E501

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd

DIRECT_PHENOTYPE_CATEGORIES = {"knockdown", "null_mutant", "functional_mutagenesis"}


def _minmax(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    minimum = numeric.min()
    maximum = numeric.max()
    if pd.isna(minimum) or pd.isna(maximum):
        return pd.Series(np.nan, index=values.index, dtype=float)
    if maximum == minimum:
        return pd.Series(0.5, index=values.index, dtype=float)
    return (numeric - minimum) / (maximum - minimum)


def pareto_fronts(frame: pd.DataFrame, metrics: Iterable[str]) -> pd.Series:
    """Assign deterministic non-dominated fronts while maximizing all metrics."""
    columns = list(metrics)
    values = frame[columns].apply(pd.to_numeric, errors="coerce").fillna(-np.inf).to_numpy()
    remaining = list(range(len(frame)))
    fronts = np.zeros(len(frame), dtype=int)
    front = 1
    while remaining:
        selected: list[int] = []
        for index in remaining:
            vector = values[index]
            dominated = any(
                other != index
                and np.all(values[other] >= vector)
                and np.any(values[other] > vector)
                for other in remaining
            )
            if not dominated:
                selected.append(index)
        if not selected:
            selected = [remaining[0]]
        fronts[selected] = front
        chosen = set(selected)
        remaining = [index for index in remaining if index not in chosen]
        front += 1
    return pd.Series(fronts, index=frame.index, dtype=int)


def build_evidence_aware_candidates(
    mapping: pd.DataFrame,
    outcomes: pd.DataFrame,
    gene_scores: pd.DataFrame,
    evidence: pd.DataFrame,
    categories: pd.DataFrame,
) -> pd.DataFrame:
    single = outcomes[
        outcomes["event_class"].eq("single_guide_indel") & outcomes["indel_size_bp"].ne(0)
    ].copy()
    disruption = (
        single.groupby("candidate_id", sort=True)
        .agg(
            frameshift_size_fraction=("frameshift", "mean"),
            median_retained_protein_fraction=("retained_protein_fraction", "median"),
            predicted_premature_stop_fraction=(
                "premature_stop_position_aa",
                lambda x: x.notna().mean(),
            ),
        )
        .reset_index()
    )
    candidates = mapping.merge(disruption, on="candidate_id", how="left", validate="one_to_one")
    candidates = candidates.merge(gene_scores, on="gene_name", how="left", validate="many_to_one")
    candidates = candidates.merge(categories, on="gene_name", how="left", validate="many_to_one")
    # The mapping table is intentionally wide. Consolidate its blocks before adding
    # derived showcase columns to avoid pandas fragmentation and noisy warnings.
    candidates = candidates.copy()
    candidates["cut_inside_interpro_entry"] = candidates["cut_domain_accessions"].fillna("").ne("")
    candidates["candidate_predicted_disruption_score"] = (
        0.5 * candidates["frameshift_size_fraction"]
        + 0.3 * (1 - candidates["median_retained_protein_fraction"])
        + 0.2 * candidates["cut_inside_interpro_entry"].astype(float)
    )
    direct_hsv2 = evidence[
        evidence["virus_type"].eq("HSV-2")
        & evidence["evidence_strength"].eq("direct")
        & evidence["evidence_category"].isin(DIRECT_PHENOTYPE_CATEGORIES)
    ]
    direct_hsv2_genes = set(direct_hsv2["gene_name"])
    candidates["direct_hsv2_phenotype_evidence"] = candidates["gene_name"].isin(direct_hsv2_genes)
    candidates["evidence_tier"] = np.select(
        [
            candidates["evidence_based_essentiality_score"].notna(),
            candidates["direct_hsv2_phenotype_evidence"],
            candidates["hsv1_ortholog_essentiality_score"].notna(),
        ],
        ["direct_HSV2_essentiality", "direct_HSV2_phenotype", "HSV1_ortholog_only"],
        default="function_only_or_unknown",
    )
    candidates["targetability_percentile"] = _minmax(candidates["post_human_score"])
    candidates["disruption_percentile"] = _minmax(
        candidates["candidate_predicted_disruption_score"]
    )
    candidates["evidence_coverage_percentile"] = _minmax(candidates["evidence_coverage_score"])
    candidates["pareto_front"] = pareto_fronts(
        candidates,
        [
            "targetability_percentile",
            "disruption_percentile",
            "evidence_coverage_percentile",
        ],
    )
    candidates["priority_interpretation"] = (
        "Multi-objective research priority only; not a therapeutic or safety ranking."
    )
    return candidates.sort_values(
        ["pareto_front", "post_human_rank", "candidate_id"], kind="mergesort"
    ).reset_index(drop=True)


def select_balanced_deep_panel(
    candidates: pd.DataFrame,
    *,
    per_gene: int = 4,
) -> pd.DataFrame:
    working = candidates.sort_values(
        [
            "gene_name",
            "pareto_front",
            "post_human_rank",
            "candidate_predicted_disruption_score",
            "candidate_id",
        ],
        ascending=[True, True, True, False, True],
        kind="mergesort",
    ).copy()
    working["deep_panel_rank_within_gene"] = working.groupby("gene_name").cumcount() + 1
    panel = working[working["deep_panel_rank_within_gene"].le(per_gene)].copy()
    panel["deep_panel_selection_reason"] = (
        "balanced gene quota; multi-objective front; coordinate-level protein context retained"
    )
    return panel.sort_values(
        ["primary_category", "gene_name", "deep_panel_rank_within_gene"], kind="mergesort"
    ).reset_index(drop=True)


def _unique_gene_top(frame: pd.DataFrame, count: int) -> pd.DataFrame:
    return frame.drop_duplicates("gene_name", keep="first").head(count).copy()


def build_comparison_sets(candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    sets: list[pd.DataFrame] = []

    ranking = _unique_gene_top(
        candidates.sort_values(["post_human_rank", "candidate_id"], kind="mergesort"), 3
    )
    ranking["strategy"] = "ranking_only"
    ranking["strategy_rationale"] = "Top coordinate-level ranks without biological balancing"
    sets.append(ranking)

    disruption = _unique_gene_top(
        candidates.sort_values(
            ["candidate_predicted_disruption_score", "post_human_rank", "candidate_id"],
            ascending=[False, True, True],
            kind="mergesort",
        ),
        3,
    )
    disruption["strategy"] = "predicted_disruption"
    disruption["strategy_rationale"] = "Highest size-model disruption scores across unique genes"
    sets.append(disruption)

    evidence = candidates[candidates["evidence_tier"].ne("function_only_or_unknown")]
    evidence = _unique_gene_top(
        evidence.sort_values(
            [
                "direct_hsv2_phenotype_evidence",
                "hsv1_ortholog_essentiality_score",
                "post_human_rank",
            ],
            ascending=[False, False, True],
            na_position="last",
            kind="mergesort",
        ),
        3,
    )
    evidence["strategy"] = "evidence_anchored"
    evidence["strategy_rationale"] = (
        "Source-linked evidence prioritized while HSV-1 and HSV-2 scopes remain explicit"
    )
    sets.append(evidence)

    diverse = (
        candidates.sort_values(
            ["primary_category", "pareto_front", "post_human_rank", "candidate_id"],
            kind="mergesort",
        )
        .drop_duplicates("primary_category", keep="first")
        .head(4)
        .copy()
    )
    diverse["strategy"] = "mechanistically_diverse"
    diverse["strategy_rationale"] = "One computational representative per available category"
    sets.append(diverse)

    members = pd.concat(sets, ignore_index=True)
    members["strategy_member_order"] = members.groupby("strategy").cumcount() + 1
    summary = (
        members.groupby("strategy", sort=False)
        .agg(
            candidate_count=("candidate_id", "size"),
            unique_gene_count=("gene_name", "nunique"),
            mechanistic_category_count=("primary_category", "nunique"),
            minimum_exact_strain_coverage=("exact_strain_coverage", "min"),
            mean_sequence_targetability=("post_human_score", "mean"),
            mean_predicted_disruption=("candidate_predicted_disruption_score", "mean"),
            mean_evidence_coverage=("evidence_coverage_score", "mean"),
            direct_hsv2_phenotype_gene_count=("direct_hsv2_phenotype_evidence", "sum"),
        )
        .reset_index()
    )
    summary["interpretation"] = (
        "Comparison set only; no joint editing, delivery, efficacy, or safety prediction."
    )
    return members, summary


def build_research_findings(
    gene_rankings: pd.DataFrame,
    gene_stability: pd.DataFrame,
    gene_scores: pd.DataFrame,
    candidates: pd.DataFrame,
    population_candidates: pd.DataFrame,
    evolution: pd.DataFrame,
) -> pd.DataFrame:
    """Build auditable candidate findings without converting hypotheses into claims."""
    rows: list[dict[str, str]] = []

    def add(
        finding_id: str,
        classification: str,
        observation: str,
        research_value: str,
        support: str,
        limitation: str,
    ) -> None:
        rows.append(
            {
                "finding_id": finding_id,
                "classification": classification,
                "observation": observation,
                "potential_research_value": research_value,
                "computational_support": support,
                "key_limitation": limitation,
            }
        )

    ranked = gene_rankings.copy()
    if not ranked.empty and "targetability_rank" in ranked:
        ranked["targetability_rank"] = pd.to_numeric(ranked["targetability_rank"], errors="coerce")
        top = ranked.sort_values(["targetability_rank", "gene_name"], kind="mergesort").iloc[0]
        benchmark = ranked[ranked["gene_name"].eq("UL30")]
        benchmark_text = (
            f"; UL30 ranks {int(benchmark.iloc[0]['targetability_rank'])}"
            if not benchmark.empty and pd.notna(benchmark.iloc[0]["targetability_rank"])
            else ""
        )
        add(
            "genome_wide_reprioritization",
            "computational_observation",
            f"{top['gene_name']} ranks first for sequence targetability{benchmark_text}.",
            "A virus-wide search can surface technically tractable genes that a benchmark-only analysis would miss.",
            "Balanced genome-wide gene ranking after the completed GRCh38.p14 screen.",
            "The rank measures model-bounded targetability, not gene importance, editing, or viral inhibition; exhaustive sensitivity analysis is pending.",
        )

    score_index = gene_scores.set_index("gene_name", drop=False) if not gene_scores.empty else None
    if score_index is not None and "UL3" in score_index.index:
        ul3 = score_index.loc["UL3"]
        add(
            "targetability_evidence_divergence",
            "evidence_gap",
            "UL3 is technically highly targetable, while direct HSV-2 essentiality remains unknown and the curated HSV-1 ortholog evidence reports nonessentiality in the tested cell-culture context.",
            "UL3 is a useful example of why guide quality and target biology must remain separate research questions.",
            f"Sequence targetability={float(ul3['sequence_targetability_score']):.3f}; HSV-2 essentiality=unknown; HSV-1 status={ul3.get('hsv1_essentiality_status', 'unknown')}.",
            "HSV-1 context cannot establish HSV-2 function, and nonessentiality in one culture system does not exclude other phenotypes.",
        )

    population = population_candidates.copy()
    if not population.empty:
        exact_ids = set(
            population.loc[
                population["population_validation_status"].eq("exact_in_all_observable_records"),
                "candidate_id",
            ].astype(str)
        )
        convergent: list[str] = []
        for gene in ("UL52", "UL30", "UL18"):
            gene_candidates = candidates[
                candidates["gene_name"].eq(gene)
                & candidates["candidate_id"].astype(str).isin(exact_ids)
            ]
            if gene_candidates.empty or score_index is None or gene not in score_index.index:
                continue
            if pd.notna(score_index.loc[gene, "hsv1_ortholog_essentiality_score"]):
                convergent.append(gene)
        if convergent:
            add(
                "multi_axis_convergence",
                "research_hypothesis",
                f"{', '.join(convergent)} each combine at least one population-supported exact target with source-linked HSV-1 ortholog evidence and a mapped protein-disruption model.",
                "These genes form a defensible evidence-aware shortlist for independent mechanistic assessment without claiming a therapeutic rank.",
                "Candidate-level host-screen, held-out locus validation, protein mapping, and separately scoped ortholog evidence.",
                "No direct HSV-2 null-essentiality result is present for these genes, and population records are mostly partial.",
            )

        ranked_candidates = candidates.sort_values(
            ["post_human_rank", "candidate_id"], kind="mergesort"
        )
        if not ranked_candidates.empty:
            leading = ranked_candidates.iloc[0]
            pop_row = population[population["candidate_id"].eq(leading["candidate_id"])]
            if not pop_row.empty and pd.notna(
                pop_row.iloc[0].get("observable_locus_exact_target_coverage")
            ):
                observed = int(pop_row.iloc[0]["locus_observable_record_count"])
                exact = int(pop_row.iloc[0]["exact_target_in_observable_locus_count"])
                add(
                    "candidate_specific_population_variation",
                    "computational_observation",
                    f"The leading coordinate-level candidate {leading['candidate_id']} ({leading['gene_name']}) retained the exact target in {exact}/{observed} observable held-out loci.",
                    "Population support should be evaluated per coordinate rather than inferred from a gene-level conservation label.",
                    "Discovery-excluded, locus-aware exact guide/PAM comparison.",
                    "Exact sequence retention is not evidence of editor activity, host safety, delivery, or efficacy; unresolved partial records remain.",
                )

    if not evolution.empty and "predicted_protein_disruption_score" in gene_scores:
        merged = gene_scores.merge(evolution, on="gene_name", how="left", validate="one_to_one")
        strongest = merged.sort_values(
            ["predicted_protein_disruption_score", "gene_name"],
            ascending=[False, True],
            kind="mergesort",
        ).iloc[0]
        add(
            "conserved_disruption_signal",
            "research_hypothesis",
            f"{strongest['gene_name']} has the highest predicted protein-disruption score ({float(strongest['predicted_protein_disruption_score']):.3f}) among the nine mapped genes and mean amino-acid conservation {float(strongest['mean_amino_acid_conservation']):.3f} in the discovery alignment.",
            "A conserved protein with cuts mapping to disruptive sequence contexts is a useful hypothesis for deeper functional prioritization.",
            "Coding-coordinate mapping, InterPro overlap, deterministic -10..+10 bp outcome enumeration, and cross-strain translation.",
            "Size-enumerated outcomes are not repair-frequency predictions, and the 14-genome discovery alignment is not population representative.",
        )

    if not gene_stability.empty and "top_10_stability" in gene_stability:
        stable = gene_stability[gene_stability["top_10_stability"].fillna(False)]
        genes = stable.sort_values("gene_name")["gene_name"].astype(str).tolist()
        if genes:
            add(
                "quota_stable_gene_set",
                "robustness_observation",
                f"{', '.join(genes)} remain in the top 10 at K=10, 25, and 50.",
                "Quota-stable genes are less likely to be artifacts of one arbitrary per-gene sampling depth.",
                "Rank sensitivity analysis across three candidate quotas.",
                "Stability is internal to the current scoring model and balanced panel; it is not biological validation.",
            )

    direct_count = int(
        gene_scores["hsv2_evidence_based_essentiality_score"].notna().sum()
        if "hsv2_evidence_based_essentiality_score" in gene_scores
        else 0
    )
    add(
        "direct_hsv2_essentiality_gap",
        "evidence_gap",
        f"Direct, scored HSV-2 essentiality evidence is available for {direct_count}/{len(gene_scores)} deeply analyzed genes.",
        "The missing-evidence map identifies where literature curation or new functional datasets would most improve prioritization.",
        "Source-linked evidence table with HSV-1 and HSV-2 scopes kept separate.",
        "Absence of curated evidence is not evidence of nonessentiality and may reflect incomplete literature coverage.",
    )
    return pd.DataFrame(rows)
