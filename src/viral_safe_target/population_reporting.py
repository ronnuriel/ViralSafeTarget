"""Reporting for held-out, locus-aware viral population validation."""

from __future__ import annotations

import html
import json
from pathlib import Path

import pandas as pd

FOCUS_GENES = ("UL3", "UL10", "UL18", "UL20", "UL36", "UL52", "UL53", "UL19", "UL30")


def _split_genes(value: object) -> list[str]:
    genes = [item.strip() for item in str(value or "").split(";") if item.strip()]
    return genes or ["unmapped"]


def build_population_comparison(
    candidates: pd.DataFrame,
    locus_validation: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Join discovery and held-out validation without altering targetability scores."""
    validation_columns = [
        "candidate_id",
        "population_record_count",
        "locus_observable_record_count",
        "exact_target_in_observable_locus_count",
        "observable_locus_exact_target_coverage",
        "observable_locus_without_exact_target_count",
        "locus_unresolved_record_count",
        "exact_target_anywhere_count",
        "locus_validation_interpretation",
    ]
    missing = sorted(set(validation_columns) - set(locus_validation.columns))
    if missing:
        raise ValueError(f"Locus validation table lacks columns: {missing}")
    comparison = candidates.merge(
        locus_validation[validation_columns], on="candidate_id", how="left", validate="one_to_one"
    )
    fallback_genes = comparison.get("gene_name", pd.Series("unmapped", index=comparison.index))
    comparison["mapped_gene_names"] = comparison.get("mapped_gene_names", fallback_genes).fillna(
        "unmapped"
    )
    comparison["reference_unique_target"] = pd.to_numeric(
        comparison.get("reference_viral_occurrence_count", 1), errors="coerce"
    ).eq(1)
    observable = pd.to_numeric(comparison["locus_observable_record_count"], errors="coerce").fillna(
        0
    )
    exact = pd.to_numeric(
        comparison["exact_target_in_observable_locus_count"], errors="coerce"
    ).fillna(0)
    comparison["population_validation_status"] = "not_evaluable"
    comparison.loc[observable.gt(0) & exact.eq(observable), "population_validation_status"] = (
        "exact_in_all_observable_records"
    )
    comparison.loc[observable.gt(0) & exact.lt(observable), "population_validation_status"] = (
        "exact_target_not_seen_in_all_observable_records"
    )
    comparison.loc[~comparison["reference_unique_target"], "population_validation_status"] = (
        "multi_locus_attribution_limited"
    )
    comparison["population_validation_used_in_targetability_score"] = False
    comparison["population_validation_limit"] = (
        "Held-out exact sequence/PAM evidence is reported separately; it is not editing, "
        "efficacy, delivery, host-safety, or therapeutic evidence."
    )

    exploded = comparison.assign(
        gene_name_for_summary=comparison["mapped_gene_names"].map(_split_genes)
    ).explode("gene_name_for_summary")
    rows: list[dict[str, object]] = []
    for gene_name, group in exploded.groupby("gene_name_for_summary", sort=True):
        evaluable = group[
            group["locus_observable_record_count"].fillna(0).gt(0)
            & group["reference_unique_target"]
        ].copy()
        supported = evaluable[
            evaluable["population_validation_status"].eq("exact_in_all_observable_records")
        ]
        ranked = group.sort_values(
            ["observable_locus_exact_target_coverage", "post_human_rank"],
            ascending=[False, True],
            na_position="last",
            kind="mergesort",
        )
        best = ranked.iloc[0]
        rows.append(
            {
                "gene_name": gene_name,
                "candidate_count": len(group),
                "reference_unique_candidate_count": int(group["reference_unique_target"].sum()),
                "locus_evaluable_unique_candidate_count": len(evaluable),
                "exact_in_all_observable_records_candidate_count": len(supported),
                "median_observable_record_count": evaluable[
                    "locus_observable_record_count"
                ].median()
                if len(evaluable)
                else pd.NA,
                "median_observable_locus_exact_target_coverage": evaluable[
                    "observable_locus_exact_target_coverage"
                ].median()
                if len(evaluable)
                else pd.NA,
                "minimum_observable_locus_exact_target_coverage": evaluable[
                    "observable_locus_exact_target_coverage"
                ].min()
                if len(evaluable)
                else pd.NA,
                "best_population_supported_candidate_id": best["candidate_id"],
                "best_population_supported_candidate_coverage": best[
                    "observable_locus_exact_target_coverage"
                ],
                "best_population_supported_candidate_post_human_rank": best.get(
                    "post_human_rank", pd.NA
                ),
                "population_evidence_status": "held_out_locus_aware_exact_sequence",
                "population_evidence_used_in_targetability_score": False,
            }
        )
    genes = pd.DataFrame(rows).sort_values(
        ["median_observable_locus_exact_target_coverage", "gene_name"],
        ascending=[False, True],
        na_position="last",
        kind="mergesort",
    )
    genes.insert(0, "population_support_rank", range(1, len(genes) + 1))
    return comparison, genes


def _table(frame: pd.DataFrame, columns: list[str], rows: int = 50) -> str:
    available = [column for column in columns if column in frame]
    if not available or frame.empty:
        return "<p>No rows available.</p>"
    return frame[available].head(rows).to_html(index=False, escape=True, border=0)


def write_population_report(
    comparison: pd.DataFrame,
    genes: pd.DataFrame,
    output_dir: str | Path,
    provenance: dict[str, object],
) -> None:
    """Write CSV, Markdown, HTML, and provenance artifacts."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(output / "candidate_population_comparison.csv", index=False)
    genes.to_csv(output / "gene_population_summary.csv", index=False)
    focus = genes[genes["gene_name"].isin(FOCUS_GENES)].copy()
    evaluable = comparison[
        comparison["locus_observable_record_count"].fillna(0).gt(0)
        & comparison["reference_unique_target"]
    ]
    fully_supported = evaluable[
        evaluable["population_validation_status"].eq("exact_in_all_observable_records")
    ]
    findings_lines = [
        "# Held-out HSV-2 population validation",
        "",
        (
            f"- Population records: {provenance['population_record_count']:,}; "
            "discovery genomes excluded."
        ),
        f"- Candidates evaluated: {len(comparison):,}.",
        (
            "- Reference-unique candidates with an observable locus in at least one "
            f"held-out record: {len(evaluable):,}."
        ),
        (
            "- Candidates exact in every record where their reference locus was observable: "
            f"{len(fully_supported):,}."
        ),
        "",
        (
            "Absence from a partial record is not counted as a variant unless a high-quality "
            "reference alignment covers the candidate interval."
        ),
        (
            "Population validation remains separate from sequence targetability and does not "
            "establish editing, efficacy, delivery, host safety, or therapeutic benefit."
        ),
    ]
    findings = "\n".join(findings_lines) + "\n"
    (output / "population_validation_findings.md").write_text(findings, encoding="utf-8")
    sections = [
        ("Study boundary", f"<pre>{html.escape(findings)}</pre>"),
        (
            "Focus-gene summary",
            _table(
                focus,
                [
                    "population_support_rank",
                    "gene_name",
                    "candidate_count",
                    "locus_evaluable_unique_candidate_count",
                    "exact_in_all_observable_records_candidate_count",
                    "median_observable_record_count",
                    "median_observable_locus_exact_target_coverage",
                    "minimum_observable_locus_exact_target_coverage",
                    "best_population_supported_candidate_id",
                ],
            ),
        ),
        (
            "All genes",
            _table(
                genes,
                [
                    "population_support_rank",
                    "gene_name",
                    "candidate_count",
                    "locus_evaluable_unique_candidate_count",
                    "exact_in_all_observable_records_candidate_count",
                    "median_observable_locus_exact_target_coverage",
                    "best_population_supported_candidate_id",
                ],
                rows=100,
            ),
        ),
        (
            "Candidate audit",
            _table(
                comparison.sort_values("post_human_rank", na_position="last"),
                [
                    "post_human_rank",
                    "candidate_id",
                    "mapped_gene_names",
                    "reference_unique_target",
                    "locus_observable_record_count",
                    "exact_target_in_observable_locus_count",
                    "observable_locus_exact_target_coverage",
                    "locus_unresolved_record_count",
                    "population_validation_status",
                ],
                rows=100,
            ),
        ),
        ("Provenance", f"<pre>{html.escape(json.dumps(provenance, indent=2))}</pre>"),
    ]
    body = "".join(
        f"<section><h2>{html.escape(title)}</h2>{content}</section>" for title, content in sections
    )
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Held-out HSV-2 population validation</title>
<style>
body{{font:16px/1.5 system-ui,sans-serif;max-width:1440px;margin:auto;padding:2rem;color:#172b4d}}
h1,h2{{color:#163a5f}}section{{margin:2rem 0}}table{{border-collapse:collapse;width:100%;
display:block;overflow:auto;font-size:.84rem}}th,td{{border:1px solid #d9e2ec;padding:.45rem;
text-align:left;white-space:nowrap}}th{{background:#eef4f8}}pre{{white-space:pre-wrap;
background:#f4f7f9;padding:1rem}}</style></head><body>
<h1>Held-out, locus-aware HSV-2 population validation</h1>
<p><strong>Interpretation:</strong> population-genomic evidence only. No wet-lab, safety,
or therapeutic claim is made.</p>{body}</body></html>
"""
    (output / "population_validation_report.html").write_text(document, encoding="utf-8")
    (output / "population_report_manifest.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
