"""Researcher-facing HTML, methods, and limitations outputs."""
# ruff: noqa: E501

from __future__ import annotations

import html
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

NOTICE = (
    "Computational research only. Candidates and predicted off-target risks require "
    "experimental and expert validation; these outputs do not establish safety, viral "
    "inactivation, treatment efficacy, or a cure."
)


def _table(frame: pd.DataFrame, columns: list[str], rows: int = 50) -> str:
    selected = [column for column in columns if column in frame]
    if frame.empty or not selected:
        return "<p>No data available for this section.</p>"
    return (
        frame[selected]
        .head(rows)
        .to_html(index=False, escape=True, float_format=lambda value: f"{value:.3f}")
    )


def write_html_report(
    candidates: pd.DataFrame,
    output_path: str | Path,
    title: str = "ViralSafeTarget report",
    *,
    rejected: pd.DataFrame | None = None,
    pairs: pd.DataFrame | None = None,
    predicted_hits: pd.DataFrame | None = None,
    output_links: list[str] | None = None,
) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    rejected = rejected if rejected is not None else pd.DataFrame()
    pairs = pairs if pairs is not None else pd.DataFrame()
    predicted_hits = predicted_hits if predicted_hits is not None else pd.DataFrame()
    gene_counts = (
        candidates.get("gene_name", pd.Series(dtype=str))
        .fillna("unannotated")
        .replace("", "unannotated")
        .value_counts()
        .rename_axis("gene")
        .reset_index(name="candidate_count")
    )
    components = [
        column
        for column in [
            "conservation_score",
            "viral_uniqueness_score",
            "gc_fraction",
            "gc_score",
            "sequence_complexity_score",
            "annotation_score",
            "gene_evidence_score",
            "pre_human_score",
            "post_human_score",
        ]
        if column in candidates
    ]
    distribution = (
        candidates[components].describe().transpose().reset_index(names="component")
        if components
        else pd.DataFrame()
    )
    uniqueness = (
        candidates.get("reference_viral_occurrence_count", pd.Series(dtype=float))
        .value_counts()
        .sort_index()
        .rename_axis("reference_occurrence_count")
        .reset_index(name="candidate_count")
    )
    links = "".join(
        f'<li><a href="{html.escape(link)}">{html.escape(link)}</a></li>'
        for link in (output_links or [])
    )
    candidate_columns = [
        "candidate_id",
        "gene_name",
        "feature_type",
        "guide_sequence",
        "pam",
        "exact_strain_coverage",
        "reference_viral_occurrence_count",
        "gc_fraction",
        "sequence_complexity_score",
        "gene_evidence_score",
        "pre_human_score",
        "post_human_score",
        "decision",
        "rank_explanation",
    ]
    pair_columns = [
        "candidate_a",
        "candidate_b",
        "gene_a",
        "gene_b",
        "hypothesis_type",
        "distance_bp",
        "deletion_length_bp",
        "joint_strain_coverage",
        "pair_score",
        "interpretation",
    ]
    html_text = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>{html.escape(title)}</title><style>
body {{ font-family: Arial,sans-serif; max-width: 1500px; margin: 28px auto; line-height: 1.45; }}
table {{ border-collapse: collapse; width: 100%; display: block; overflow-x: auto; }}
th,td {{ border:1px solid #ddd; padding:6px; text-align:left; vertical-align:top; }}
th {{ background:#f2f2f2; }} .notice {{ padding:12px; background:#fff4df; border:1px solid #b5843c; }}
</style></head><body>
<h1>{html.escape(title)}</h1><p>Generated: {datetime.now(timezone.utc).isoformat()}</p>
<div class="notice"><strong>Scope:</strong> {html.escape(NOTICE)}</div>
<h2>1. Dataset summary</h2><p>Retained candidates: {len(candidates):,}; rejected candidates: {len(rejected):,}.</p>
<h2>2. Quality-control summary</h2><p>See the run manifest for input hashes, accepted accessions, rejected accessions, and configuration.</p>
<h2>3. Candidate filtering funnel</h2>{_table(pd.DataFrame([{"scanned_or_ranked": len(candidates) + len(rejected), "retained": len(candidates), "rejected": len(rejected)}]), ["scanned_or_ranked", "retained", "rejected"])}
<h2>4. Ranking component distributions</h2>{_table(distribution, list(distribution.columns), 100)}
<h2>5. Candidates by gene</h2>{_table(gene_counts, ["gene", "candidate_count"], 100)}
<h2>6. Viral uniqueness distribution</h2>{_table(uniqueness, ["reference_occurrence_count", "candidate_count"], 100)}
<h2>7. GC and complexity distributions</h2>{_table(distribution[distribution.get("component", pd.Series(dtype=str)).isin(["gc_fraction", "sequence_complexity_score"])], list(distribution.columns), 10)}
<h2>8. Human off-target summary</h2><p>Predicted hits reported: {len(predicted_hits):,}. Absence of a predicted hit within the configured threshold does not establish safety.</p>{_table(predicted_hits, ["candidate_id", "chromosome", "human_coordinate_1based", "direction", "mismatches", "off_target_sequence", "human_annotation"], 50)}
<h2>9. Top candidates with explanations</h2>{_table(candidates, candidate_columns, 50)}
<h2>10. Pair hypotheses</h2>{_table(pairs, pair_columns, 50)}
<h2>11. Limitations and non-clinical-use warning</h2><p>{html.escape(NOTICE)} The model excludes delivery, accessibility, editing efficiency, repair outcomes, toxicity, latent infection biology, and reactivation.</p>
<h2>12. Machine-readable outputs</h2><ul>{links}</ul>
</body></html>"""
    output.write_text(html_text, encoding="utf-8")
    return output


def write_methods_and_limitations(output_directory: str | Path) -> tuple[Path, Path]:
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    methods = output / "methods.md"
    limitations = output / "limitations.md"
    methods.write_text(
        "# Methods\n\nCandidates are scanned on both strands under the configured editor profile, "
        "assigned content-derived IDs, annotated against GFF3, and ranked from visible sequence, "
        "conservation, uniqueness, annotation, and curated-evidence components. Human hits are "
        "enumerated externally and summarized separately. Pair hypotheses use deterministic "
        "ranked and gene-stratified selection. Exact settings and hashes are in `run_manifest.json`.\n",
        encoding="utf-8",
    )
    limitations.write_text(
        "# Limitations\n\n" + NOTICE + "\n\nThe model does not include latent chromatin "
        "accessibility, delivery, editing efficiency, repair distributions, toxicity, immune "
        "effects, reactivation, or clinical outcomes. Sequence disruption does not prove viral "
        "inactivation, and a predicted off-target scan does not establish safety.\n",
        encoding="utf-8",
    )
    return methods, limitations
