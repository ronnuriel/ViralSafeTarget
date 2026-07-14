"""Researcher-facing HTML for the genome-wide discovery workflow."""

# HTML prose and column declarations are intentionally kept readable as complete strings.
# ruff: noqa: E501

from __future__ import annotations

import html
from pathlib import Path
from typing import Any

import pandas as pd

OLD_PILOT_CANDIDATE = "VST-240e20eb666f9c85"


def _table(frame: pd.DataFrame, columns: list[str], *, rows: int = 20) -> str:
    available = [column for column in columns if column in frame]
    if frame.empty or not available:
        return '<p class="pending">No completed result is available for this table.</p>'
    return frame[available].head(rows).to_html(index=False, border=0, classes="data")


def _gene_min(genes: pd.DataFrame, gene_name: str, column: str) -> float | None:
    values = pd.to_numeric(
        genes.loc[genes["gene_name"].eq(gene_name), column], errors="coerce"
    ).dropna()
    return float(values.min()) if not values.empty else None


def write_discovery_report(
    output_path: str | Path,
    *,
    candidates: pd.DataFrame,
    feature_map: pd.DataFrame,
    genes: pd.DataFrame,
    stability: pd.DataFrame,
    deep_panel: pd.DataFrame,
    top_per_gene_candidates: pd.DataFrame,
    same_pairs: pd.DataFrame,
    multi_pairs: pd.DataFrame,
    genes_without_candidates: pd.DataFrame,
    qc: pd.DataFrame,
    provenance: dict[str, Any],
    initial_candidate_count: int,
    eligible_candidate_count: int,
) -> dict[str, Any]:
    """Write all 25 required sections and return explicit benchmark answers."""
    completed = candidates[candidates["screening_status"].eq("completed")]
    ranked = candidates[candidates["post_human_rank"].notna()].sort_values("post_human_rank")
    ul30_rank = _gene_min(genes, "UL30", "targetability_rank")
    best_gene_rank = pd.to_numeric(genes["targetability_rank"], errors="coerce").min()
    another_gene_above_ul30 = (
        bool(best_gene_rank < ul30_rank)
        if ul30_rank is not None and pd.notna(best_gene_rank)
        else None
    )
    old = ranked[ranked["candidate_id"].eq(OLD_PILOT_CANDIDATE)]
    another_candidate_above_old = (
        bool(float(old.iloc[0]["post_human_rank"]) > 1) if not old.empty else None
    )
    stable_top = stability[stability.get("top_10_stability", False).fillna(False)]
    accepted = (
        int(qc.get("decision", pd.Series(dtype=str)).astype(str).eq("accepted").sum())
        if not qc.empty
        else 0
    )
    rejected = len(qc) - accepted if not qc.empty else 0
    annotated_gene_count = int(
        (feature_map["feature_type"] == "gene").groupby(feature_map["feature_id"]).any().sum()
    )
    total_features = int(
        feature_map.loc[feature_map["mapping_status"].ne("intergenic"), "feature_id"].nunique()
    )
    screening_fraction = len(completed) / len(candidates) if len(candidates) else 0.0
    selection_mode = str(provenance.get("selection_mode", "balanced"))
    panel_label = "exhaustive panel" if selection_mode == "exhaustive" else "balanced panel"

    if another_gene_above_ul30 is None:
        ul30_answer = "Not yet determinable: UL30 lacks a completed comparable gene rank."
    elif another_gene_above_ul30:
        ul30_answer = (
            "Yes. At least one other gene has a better computational targetability rank than UL30."
        )
    else:
        ul30_answer = (
            "No gene currently has a better completed computational targetability rank than UL30."
        )
    if another_candidate_above_old is None:
        candidate_answer = f"Not yet determinable: {OLD_PILOT_CANDIDATE} does not have a completed rank in this run."
    elif another_candidate_above_old:
        candidate_answer = f"Yes. At least one candidate ranks above {OLD_PILOT_CANDIDATE}."
    else:
        candidate_answer = f"No. {OLD_PILOT_CANDIDATE} is currently the top completed candidate."

    benchmark = genes[genes["gene_name"].isin(["UL19", "UL30"])].sort_values(
        ["gene_name", "feature_id"], kind="mergesort"
    )
    newly_surfaced = genes[~genes["gene_name"].isin(["UL19", "UL30"])].sort_values(
        ["targetability_rank", "gene_name"], na_position="last", kind="mergesort"
    )
    mismatch_counts = pd.DataFrame(
        {
            "mismatch_level": [0, 1, 2, 3],
            "candidate_count_with_at_least_one_hit": [
                int(completed.get(column, pd.Series(dtype=float)).gt(0).sum())
                for column in [
                    "human_exact_hit_count",
                    "human_one_mismatch_hit_count",
                    "human_two_mismatch_hit_count",
                    "human_three_mismatch_hit_count",
                ]
            ],
        }
    )
    status_counts = (
        candidates["screening_status"]
        .value_counts()
        .rename_axis("status")
        .reset_index(name="candidate_count")
    )
    funnel = pd.DataFrame(
        {
            "stage": [
                "initial",
                "eligible pre-human",
                panel_label,
                "human screen completed",
                "post-human ranked",
            ],
            "candidate_count": [
                initial_candidate_count,
                eligible_candidate_count,
                len(candidates),
                len(completed),
                len(ranked),
            ],
        }
    )
    evidence_message = (
        "Biological evidence not supplied. Missing evidence remains missing and contributes no negative score."
        if genes.get("evidence_coverage", pd.Series(dtype=float)).isna().all()
        or pd.to_numeric(genes.get("evidence_coverage", 0), errors="coerce").fillna(0).eq(0).all()
        else "Optional biological evidence was supplied and remains separate from computational targetability."
    )
    high_targetability_low_evidence = genes.sort_values(
        ["targetability_rank", "gene_name"], na_position="last", kind="mergesort"
    )
    high_targetability_low_evidence = high_targetability_low_evidence[
        pd.to_numeric(high_targetability_low_evidence["evidence_coverage"], errors="coerce")
        .fillna(0)
        .eq(0)
    ]
    high_evidence_low_targetability = genes[
        pd.to_numeric(genes["evidence_coverage"], errors="coerce").fillna(0).gt(0)
    ].sort_values(
        ["targetability_rank", "gene_name"],
        ascending=[False, True],
        na_position="last",
        kind="mergesort",
    )
    answer_box = (
        f"<p><strong>Did another gene rank above UL30?</strong> {html.escape(ul30_answer)}</p>"
        f"<p><strong>Did another candidate rank above the old pilot candidate?</strong> {html.escape(candidate_answer)}</p>"
        f"<p><strong>Genome coverage:</strong> {len(completed):,}/{len(candidates):,} panel candidates "
        f"({screening_fraction:.1%}) have completed Cas-OFFinder results.</p>"
        f"<p><strong>Quota stability:</strong> {len(stable_top):,} genes remain in the top 10 at K=10, 25, and 50.</p>"
    )
    sections = [
        (
            "Research question",
            "Which annotated HSV-2 genes and candidate sites have the strongest computational targetability under the declared genome, editor, conservation, and predicted-human-off-target model?"
            + answer_box,
        ),
        (
            "Input genomes and QC",
            f"<p>Human reference: {html.escape(str(provenance['human_assembly']))} ({html.escape(str(provenance['human_assembly_accession']))}). Viral accession QC records: {len(qc):,}.</p>",
        ),
        (
            "Accepted and rejected genome summary",
            f"<p>Accepted: {accepted:,}; rejected: {rejected:,}. Rejection remains visible rather than silently deleting genomes.</p>"
            + _table(qc, ["accession", "accepted", "reason", "sequence_length"], rows=50),
        ),
        (
            "Complete annotation coverage",
            f"<p>The normalized map contains {len(feature_map):,} one-to-many candidate-feature rows, including explicit intergenic mappings.</p>",
        ),
        (
            "Number of genes/features evaluated",
            f"<p>Annotated genes represented in the map: {annotated_gene_count:,}; distinct mapped feature IDs: {total_features:,}; gene ranking rows: {len(genes):,}.</p>",
        ),
        ("Candidate funnel", _table(funnel, ["stage", "candidate_count"])),
        (
            "Candidate selection method",
            (
                "<p>The exhaustive panel contains every retained pre-human candidate; "
                "stable candidate IDs and unique Cas-OFFinder queries preserve coordinate "
                "mapping without redundant searches.</p>"
                if selection_mode == "exhaustive"
                else "<p>The balanced panel is the deterministic union of the configured "
                "top candidates per annotated gene and global pre-human leaders. Selection "
                "uses no human-screen outcome; stable candidate IDs deduplicate the union, "
                "while one-to-many feature memberships remain in the mapping table.</p>"
            ),
        ),
        (
            "Cas-OFFinder execution and completeness",
            _table(status_counts, ["status", "candidate_count"])
            + "<p>A completed zero-hit search is distinct from a pending, failed, or missing batch. Zero predicted hits is not proof of safety.</p>",
        ),
        (
            "Ranking of all genes",
            _table(
                genes,
                [
                    "targetability_rank",
                    "gene_name",
                    "feature_id",
                    "eligible_candidate_count",
                    "screened_candidate_count",
                    "clean_fraction_wilson_lower_95",
                    "best_candidate_rank",
                    "median_top_5_candidate_rank",
                    "targetability_score",
                    "confidence_level",
                ],
                rows=100,
            ),
        ),
        (
            "Top gene by best candidate",
            _table(
                genes.sort_values("best_single_candidate_rank", na_position="last"),
                [
                    "best_single_candidate_rank",
                    "gene_name",
                    "best_candidate_id",
                    "best_candidate_rank",
                    "best_candidate_score",
                ],
            ),
        ),
        (
            "Top gene by robust top-5 performance",
            _table(
                genes.sort_values("top5_robustness_rank", na_position="last"),
                [
                    "top5_robustness_rank",
                    "gene_name",
                    "median_top_5_candidate_rank",
                    "worst_top_5_rank",
                    "confidence_level",
                ],
            ),
        ),
        (
            "Top gene by clean fraction",
            _table(
                genes.sort_values("clean_fraction_rank", na_position="last"),
                [
                    "clean_fraction_rank",
                    "gene_name",
                    "clean_fraction",
                    "clean_fraction_wilson_lower_95",
                    "screened_candidate_count",
                ],
            ),
        ),
        (
            "Targetability versus biological-evidence matrix",
            f"<p>{html.escape(evidence_message)}</p>"
            + _table(
                genes,
                [
                    "gene_name",
                    "targetability_rank",
                    "targetability_score",
                    "evidence_coverage",
                    "biological_evidence_status",
                    "confidence_level",
                ],
                rows=100,
            )
            + "<h3>High-targetability, low-evidence discovery view</h3>"
            + _table(
                high_targetability_low_evidence,
                ["gene_name", "targetability_rank", "evidence_coverage", "confidence_level"],
                rows=30,
            )
            + "<h3>High-evidence, low-targetability view</h3>"
            + _table(
                high_evidence_low_targetability,
                ["gene_name", "targetability_rank", "evidence_coverage", "confidence_level"],
                rows=30,
            ),
        ),
        (
            "Gene ranking sensitivity at K=10, K=25 and K=50",
            _table(
                stability,
                [
                    "gene_name",
                    "rank_at_k10",
                    "rank_at_k25",
                    "rank_at_k50",
                    "rank_range",
                    "rank_std",
                    "top_5_stability",
                    "top_10_stability",
                    "stability_warning",
                ],
                rows=100,
            ),
        ),
        (
            "Global top candidates",
            _table(
                ranked,
                [
                    "post_human_rank",
                    "candidate_id",
                    "guide_sequence",
                    "pam",
                    "mapped_gene_names",
                    "pre_human_score",
                    "post_human_score",
                    "human_total_predicted_hits",
                    "decision",
                ],
                rows=100,
            ),
        ),
        (
            "Top candidates per gene",
            _table(
                top_per_gene_candidates,
                [
                    "mapped_gene_for_view",
                    "post_human_rank",
                    "candidate_id",
                    "mapped_gene_names",
                    "post_human_score",
                    "decision",
                ],
                rows=750,
            ),
        ),
        (
            "UL19 and UL30 benchmark comparison",
            "<p>UL19 and UL30 are shown as benchmarks only and receive no scoring advantage.</p>"
            + _table(
                benchmark,
                [
                    "gene_name",
                    "targetability_rank",
                    "best_candidate_rank",
                    "median_top_5_candidate_rank",
                    "clean_fraction_wilson_lower_95",
                    "confidence_level",
                ],
            ),
        ),
        (
            "Newly surfaced genes not present in the old pilot",
            _table(
                newly_surfaced,
                [
                    "targetability_rank",
                    "gene_name",
                    "best_candidate_id",
                    "best_candidate_rank",
                    "median_top_5_candidate_rank",
                    "confidence_level",
                ],
                rows=30,
            ),
        ),
        (
            "Human-hit mismatch distribution",
            _table(mismatch_counts, ["mismatch_level", "candidate_count_with_at_least_one_hit"]),
        ),
        (
            "Deep-screening panel",
            f"<p>{len(deep_panel):,} deduplicated candidates are prepared for optional expensive or external tools.</p>"
            + _table(
                deep_panel,
                ["post_human_rank", "candidate_id", "mapped_gene_names", "deep_panel_reason"],
                rows=50,
            ),
        ),
        (
            "Pair hypotheses",
            f"<p>Same-gene hypotheses: {len(same_pairs):,}; multi-target hypotheses: {len(multi_pairs):,}. Cross-gene rows describe two separate target sites, not one physical deletion.</p>"
            + _table(
                same_pairs,
                [
                    "candidate_a",
                    "candidate_b",
                    "gene_a",
                    "distance_bp",
                    "joint_strain_coverage",
                    "pair_score",
                    "limitations",
                ],
            ),
        ),
        (
            "Pending external tools",
            "<p>CRISPRitz, CRISPOR, CHOPCHOP, and GuideScan2 inputs/import templates are optional and pending unless separately executed. Their absence does not masquerade as a favorable result.</p>",
        ),
        ("Provenance", f"<pre>{html.escape(str(provenance))}</pre>"),
        (
            "Limitations",
            "<p>This is a computational prioritization. Rankings depend on accepted viral genomes, the reference annotation, SpCas9 model, GRCh38.p14 reference search, mismatch threshold, quota, and available tool/evidence coverage. Variant-aware, bulge-aware, chromatin, delivery, phenotype, and experimental efficacy evidence are not established here.</p>",
        ),
        (
            "Non-clinical and non-experimental warning",
            "<p><strong>Research use only.</strong> These are computational candidates for further expert review, not validated interventions, clinical recommendations, safety conclusions, wet-lab protocols, or claims of cure.</p>",
        ),
    ]
    rendered = "".join(
        f'<section id="section-{index}"><h2>{index}. {html.escape(title)}</h2>{body}</section>'
        for index, (title, body) in enumerate(sections, start=1)
    )
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>HSV-2 genome-wide computational discovery</title>
<style>body{{font:16px/1.5 system-ui,sans-serif;max-width:1440px;margin:auto;padding:2rem;color:#15202b}}h1,h2{{color:#183b56}}section{{margin:2.2rem 0}}table{{border-collapse:collapse;width:100%;font-size:.85rem;display:block;overflow:auto}}th,td{{padding:.45rem;border:1px solid #d9e2ec;text-align:left;white-space:nowrap}}th{{background:#eef4f8}}.warning{{background:#fff4e5;padding:1rem;border-left:5px solid #d97706}}.pending{{color:#8a4b08}}pre{{white-space:pre-wrap;overflow-wrap:anywhere;background:#f4f7f9;padding:1rem}}</style>
</head><body><h1>HSV-2 genome-wide computational target discovery</h1>
<p class="warning">Computational research output only. A zero predicted human hit is model-bounded and is not evidence of safety.</p>
{rendered}</body></html>"""
    output = Path(output_path)
    output.write_text(document, encoding="utf-8")
    return {
        "another_gene_above_ul30": another_gene_above_ul30,
        "another_candidate_above_old_pilot": another_candidate_above_old,
        "ul30_answer": ul30_answer,
        "old_candidate_answer": candidate_answer,
        "screening_fraction": screening_fraction,
        "stable_top_10_gene_count": len(stable_top),
        "genes_without_eligible_candidates": len(genes_without_candidates),
    }
