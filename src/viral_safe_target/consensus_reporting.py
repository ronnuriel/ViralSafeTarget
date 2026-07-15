"""Combined explainable HTML report for partial or complete multi-tool comparisons."""
# ruff: noqa: E501

from __future__ import annotations

import html
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .consensus import ComparisonResult
from .reporting import NOTICE


def _table(frame: pd.DataFrame, rows: int = 50) -> str:
    if frame.empty:
        return "<p>No data are available; this stage is pending.</p>"
    return frame.head(rows).to_html(
        index=False, escape=True, float_format=lambda value: f"{value:.3f}"
    )


def write_consensus_report(
    candidates: pd.DataFrame,
    comparison: ComparisonResult,
    output_path: str | Path,
    *,
    tool_availability: pd.DataFrame | None = None,
    experimental_results: pd.DataFrame | None = None,
    provenance: dict | None = None,
) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    availability = tool_availability if tool_availability is not None else pd.DataFrame()
    experimental = experimental_results if experimental_results is not None else pd.DataFrame()
    gene_summary = (
        candidates.groupby("gene_name", dropna=False)
        .agg(
            candidate_count=("candidate_id", "count"),
            median_post_human_score=("post_human_score", "median"),
        )
        .reset_index()
    )
    matrix_columns = [
        column
        for column in comparison.candidate_tool_matrix.columns
        if column == "candidate_id" or column == "gene_name" or column.startswith("percentile__")
    ]
    top = comparison.consensus_candidates.head(20)
    pending = comparison.tool_coverage[comparison.tool_coverage["reported_candidates"].eq(0)].copy()
    if not availability.empty and not pending.empty:
        pending = pending.merge(availability, on="tool_name", how="left")
    provenance_text = html.escape(json.dumps(provenance or {}, indent=2, sort_keys=True))
    text = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width">
<title>HSV-2 multi-tool consensus</title><style>
body {{font-family:Arial,sans-serif;max-width:1500px;margin:28px auto;line-height:1.45}}
table {{border-collapse:collapse;width:100%;display:block;overflow-x:auto}}
th,td {{border:1px solid #ddd;padding:6px;text-align:left;vertical-align:top}}
th {{background:#f1f3f5}} .notice {{padding:12px;background:#fff4df;border:1px solid #b5843c}}
pre {{background:#f6f8fa;padding:12px;overflow:auto}}
</style></head><body>
<h1>HSV-2 multi-tool consensus report</h1>
<p>Generated {datetime.now(timezone.utc).isoformat()}</p>
<div class="notice"><strong>Research scope:</strong> {html.escape(NOTICE)} Consensus is prioritization, not proof.</div>
<h2>1. Current HSV-2 evidence summary</h2><p>The input is the 32 computational candidates with no predicted human hit through three mismatches in the completed GRCh38.p14 Cas-OFFinder model.</p>
<h2>2. Input candidate set</h2><p>{len(candidates)} candidates; missing external results remain missing.</p>{_table(candidates[["candidate_id", "gene_name", "guide_sequence", "post_human_score"]])}
<h2>3. Tool availability and versions</h2>{_table(availability)}
<h2>4. Per-tool result coverage</h2>{_table(comparison.tool_coverage)}
<h2>5. Candidate-tool score/rank matrix</h2><p>Cells are within-tool desirability percentiles, not raw scores.</p>{_table(comparison.candidate_tool_matrix[matrix_columns], 32)}
<h2>6. Pairwise tool agreement</h2>{_table(comparison.model_agreement)}
<h2>7. Top-k overlap</h2>{_table(comparison.model_agreement[[column for column in comparison.model_agreement if column in {"tool_a", "tool_b", "top_k", "top_k_overlap", "jaccard_overlap"}]])}
<h2>8. Consensus ranking</h2>{_table(top)}
<h2>9. High-disagreement candidates</h2>{_table(comparison.disagreement_report)}
<h2>10. UL19 versus UL30 descriptive comparison</h2>{_table(gene_summary)}
<h2>11. Pending tools or missing exports</h2>{_table(pending)}
<h2>12. Experimental results</h2><p>Predictions and measured CRISPResso2 metrics are kept separate.</p>{_table(experimental)}
<h2>13. Provenance</h2><pre>{provenance_text}</pre>
<h2>14. Limitations and non-clinical-use warning</h2><p>{html.escape(NOTICE)} Rank consensus does not identify a biologically correct model. No predicted hit within a configured threshold is not proof of safety. All candidates require expert review and experimental validation.</p>
</body></html>"""
    output.write_text(text, encoding="utf-8")
    return output
