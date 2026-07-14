"""Balanced-versus-exhaustive discovery sensitivity comparison."""

from __future__ import annotations

import html
import json
from pathlib import Path

import pandas as pd


def compare_discovery_modes(
    balanced_candidates: pd.DataFrame,
    exhaustive_candidates: pd.DataFrame,
    balanced_genes: pd.DataFrame,
    exhaustive_genes: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    """Compare completed outputs while preserving each mode's declared rank."""
    candidate_columns = [
        "candidate_id",
        "mapped_gene_names",
        "post_human_rank",
        "post_human_score",
        "human_total_predicted_hits",
        "decision",
        "screening_status",
    ]
    gene_columns = [
        "gene_name",
        "targetability_rank",
        "targetability_score",
        "eligible_candidate_count",
        "screened_candidate_count",
        "clean_candidate_count",
        "clean_fraction_wilson_lower_95",
        "best_candidate_id",
        "best_candidate_rank",
        "confidence_level",
    ]
    balanced_candidate_view = balanced_candidates[
        [column for column in candidate_columns if column in balanced_candidates]
    ]
    exhaustive_candidate_view = exhaustive_candidates[
        [column for column in candidate_columns if column in exhaustive_candidates]
    ]
    candidates = balanced_candidate_view.merge(
        exhaustive_candidate_view,
        on="candidate_id",
        how="outer",
        suffixes=("_balanced", "_exhaustive"),
        indicator=True,
        validate="one_to_one",
    )
    candidates["post_human_rank_change_balanced_minus_exhaustive"] = pd.to_numeric(
        candidates.get("post_human_rank_balanced"), errors="coerce"
    ) - pd.to_numeric(candidates.get("post_human_rank_exhaustive"), errors="coerce")
    candidates["comparison_interpretation"] = (
        "Rank sensitivity only; neither discovery mode establishes safety, editing, "
        "essentiality, efficacy, delivery, or therapeutic value."
    )

    balanced_gene_view = balanced_genes[
        [column for column in gene_columns if column in balanced_genes]
    ]
    exhaustive_gene_view = exhaustive_genes[
        [column for column in gene_columns if column in exhaustive_genes]
    ]
    genes = balanced_gene_view.merge(
        exhaustive_gene_view,
        on="gene_name",
        how="outer",
        suffixes=("_balanced", "_exhaustive"),
        indicator=True,
        validate="one_to_one",
    )
    genes["targetability_rank_change_balanced_minus_exhaustive"] = pd.to_numeric(
        genes.get("targetability_rank_balanced"), errors="coerce"
    ) - pd.to_numeric(genes.get("targetability_rank_exhaustive"), errors="coerce")
    genes["comparison_interpretation"] = (
        "Balanced-panel sensitivity analysis versus complete eligible-candidate host screen; "
        "biological evidence remains separate."
    )
    genes = genes.sort_values(
        ["targetability_rank_exhaustive", "targetability_rank_balanced", "gene_name"],
        na_position="last",
        kind="mergesort",
    )
    summary = {
        "balanced_candidate_count": len(balanced_candidates),
        "exhaustive_candidate_count": len(exhaustive_candidates),
        "candidate_overlap_count": int(candidates["_merge"].eq("both").sum()),
        "balanced_zero_predicted_hit_count": int(
            pd.to_numeric(balanced_candidates.get("human_total_predicted_hits"), errors="coerce")
            .eq(0)
            .sum()
        ),
        "exhaustive_zero_predicted_hit_count": int(
            pd.to_numeric(exhaustive_candidates.get("human_total_predicted_hits"), errors="coerce")
            .eq(0)
            .sum()
        ),
        "balanced_top_gene": (
            balanced_genes.sort_values("targetability_rank").iloc[0]["gene_name"]
            if not balanced_genes.empty
            else None
        ),
        "exhaustive_top_gene": (
            exhaustive_genes.sort_values("targetability_rank").iloc[0]["gene_name"]
            if not exhaustive_genes.empty
            else None
        ),
        "interpretation": (
            "Computational selection-mode sensitivity only; zero predicted hits remains "
            "model-bounded and is not proof of safety."
        ),
    }
    return candidates, genes, summary


def _table(frame: pd.DataFrame, columns: list[str], rows: int = 100) -> str:
    available = [column for column in columns if column in frame]
    return frame[available].head(rows).to_html(index=False, escape=True, border=0)


def write_discovery_comparison(
    candidates: pd.DataFrame,
    genes: pd.DataFrame,
    summary: dict[str, object],
    output_dir: str | Path,
) -> None:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(output / "candidate_rank_comparison.csv", index=False)
    genes.to_csv(output / "gene_rank_comparison.csv", index=False)
    (output / "comparison_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    summary_text = html.escape(json.dumps(summary, indent=2, sort_keys=True))
    gene_table = _table(
        genes,
        [
            "gene_name",
            "targetability_rank_balanced",
            "targetability_rank_exhaustive",
            "targetability_rank_change_balanced_minus_exhaustive",
            "screened_candidate_count_balanced",
            "screened_candidate_count_exhaustive",
            "clean_candidate_count_balanced",
            "clean_candidate_count_exhaustive",
            "best_candidate_id_balanced",
            "best_candidate_id_exhaustive",
        ],
    )
    candidate_table = _table(
        candidates.sort_values("post_human_rank_exhaustive", na_position="last"),
        [
            "candidate_id",
            "mapped_gene_names_exhaustive",
            "post_human_rank_balanced",
            "post_human_rank_exhaustive",
            "post_human_rank_change_balanced_minus_exhaustive",
            "human_total_predicted_hits_balanced",
            "human_total_predicted_hits_exhaustive",
            "_merge",
        ],
    )
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Balanced versus exhaustive discovery</title><style>
body{{font:16px/1.5 system-ui,sans-serif;max-width:1440px;margin:auto;padding:2rem;color:#172b4d}}
h1,h2{{color:#163a5f}}table{{border-collapse:collapse;width:100%;display:block;overflow:auto;
font-size:.83rem}}th,td{{border:1px solid #d9e2ec;padding:.45rem;text-align:left;
white-space:nowrap}}th{{background:#eef4f8}}pre{{background:#f4f7f9;padding:1rem}}</style>
</head><body><h1>Balanced versus exhaustive genome-wide discovery</h1>
<p><strong>Scope:</strong> computational selection-mode sensitivity only. A zero-hit
prediction is not proof of safety, and biological evidence is not inferred here.</p>
<h2>Summary</h2><pre>{summary_text}</pre><h2>Gene ranks</h2>{gene_table}
<h2>Candidate ranks</h2>{candidate_table}</body></html>"""
    (output / "discovery_mode_comparison.html").write_text(document, encoding="utf-8")
