"""Rank-based multi-tool comparison with explicit missingness and disagreement."""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .provenance import sha256_file
from .tables import TOOL_RESULT_COLUMNS, CandidateTable, ToolResultTable, as_dataframe

DEFAULT_METRIC_PRIORITY = [
    "post_human_score",
    "pre_human_score",
    "predicted_specificity",
    "predicted_efficiency",
    "predicted_offtarget_burden",
]


@dataclass(frozen=True)
class ComparisonResult:
    tool_results_long: pd.DataFrame
    candidate_tool_matrix: pd.DataFrame
    consensus_candidates: pd.DataFrame
    tool_coverage: pd.DataFrame
    model_agreement: pd.DataFrame
    disagreement_report: pd.DataFrame
    unmatched_external_rows: pd.DataFrame

    def write(self, output_directory: str | Path) -> dict[str, Path]:
        output = Path(output_directory)
        output.mkdir(parents=True, exist_ok=True)
        frames = {
            "tool_results_long.csv": self.tool_results_long,
            "candidate_tool_matrix.csv": self.candidate_tool_matrix,
            "consensus_candidates.csv": self.consensus_candidates,
            "tool_coverage.csv": self.tool_coverage,
            "model_agreement.csv": self.model_agreement,
            "disagreement_report.csv": self.disagreement_report,
            "unmatched_external_rows.csv": self.unmatched_external_rows,
        }
        paths = {}
        for name, frame in frames.items():
            path = output / name
            frame.to_csv(path, index=False)
            paths[name] = path
        return paths


def _percentile_from_rank(rank: pd.Series, count: int) -> pd.Series:
    if count <= 1:
        return pd.Series(1.0, index=rank.index)
    return 1.0 - (rank - 1.0) / (count - 1.0)


def candidate_metrics_as_tool_results(
    candidates: pd.DataFrame | CandidateTable,
    *,
    source_file: str | Path | None = None,
) -> ToolResultTable:
    """Represent VST ranks and existing Cas summaries without mixing their raw scales."""
    frame = as_dataframe(candidates)
    source = Path(source_file) if source_file else None
    timestamp = datetime.now(timezone.utc).isoformat()
    definitions = [
        ("viral_safe_target_pre_human", "pre_human_score", "candidate prioritization"),
        ("viral_safe_target_post_human", "post_human_score", "candidate prioritization"),
        ("cas-offinder", "human_total_predicted_hits", "reference-genome search"),
    ]
    rows: list[dict[str, object]] = []
    for tool, metric, mode in definitions:
        if metric not in frame:
            continue
        values = pd.to_numeric(frame[metric], errors="coerce")
        ascending = metric == "human_total_predicted_hits"
        ranks = values.rank(method="min", ascending=ascending, na_option="keep")
        percentiles = _percentile_from_rank(ranks, int(values.notna().sum()))
        normalized_metric = "predicted_offtarget_burden" if tool == "cas-offinder" else metric
        for index, candidate in frame.iterrows():
            available = pd.notna(values.loc[index])
            rows.append(
                {
                    "candidate_id": candidate["candidate_id"],
                    "guide_sequence": candidate.get("guide_sequence", ""),
                    "gene_name": candidate.get("gene_name", ""),
                    "tool_name": tool,
                    "tool_version": "ViralSafeTarget 0.4" if tool.startswith("viral_") else "2.4.1",
                    "tool_mode": mode,
                    "genome_or_assembly": candidate.get("human_assembly", "")
                    if tool == "cas-offinder"
                    else candidate.get("reference_accession", ""),
                    "editor": candidate.get("editor", ""),
                    "metric_name": normalized_metric,
                    "raw_value": values.loc[index] if available else pd.NA,
                    "normalized_value": percentiles.loc[index] if available else pd.NA,
                    "rank": ranks.loc[index] if available else pd.NA,
                    "percentile_rank": percentiles.loc[index] if available else pd.NA,
                    "decision": candidate.get("decision", ""),
                    "explanation": (
                        "Rank-normalized within this tool/metric; raw values are never "
                        "averaged across tools."
                    ),
                    "source_file": str(source.resolve()) if source else "",
                    "source_file_sha256": sha256_file(source)
                    if source and source.is_file()
                    else "",
                    "command_or_import_method": "loaded from completed ViralSafeTarget run",
                    "timestamp": timestamp,
                    "status": "completed" if available else "missing",
                    "error_message": "" if available else f"{metric} is missing",
                }
            )
    return ToolResultTable.from_frame(pd.DataFrame(rows))


def _select_primary_metrics(
    results: pd.DataFrame, metric_selection: dict[str, str] | None
) -> pd.DataFrame:
    selected_rows = []
    for tool, group in results.groupby("tool_name", sort=True):
        available_metrics = set(group["metric_name"].dropna().astype(str))
        requested = (metric_selection or {}).get(tool)
        if requested:
            metric = requested
        else:
            metric = next(
                (item for item in DEFAULT_METRIC_PRIORITY if item in available_metrics), None
            )
        if metric is None:
            documented = group[group["percentile_rank"].notna()]
            if documented.empty:
                continue
            metric = sorted(documented["metric_name"].astype(str).unique())[0]
        selected_rows.append(group[group["metric_name"].eq(metric)])
    if not selected_rows:
        return results.iloc[0:0].copy()
    selected = pd.concat(selected_rows, ignore_index=True)
    return selected.sort_values(["tool_name", "candidate_id"], kind="mergesort")


def _kendall_tau(left: pd.Series, right: pd.Series) -> float:
    concordant = 0
    discordant = 0
    values = list(zip(left, right, strict=True))
    for index, (left_a, right_a) in enumerate(values):
        for left_b, right_b in values[index + 1 :]:
            product = (left_a - left_b) * (right_a - right_b)
            concordant += product > 0
            discordant += product < 0
    denominator = concordant + discordant
    return (concordant - discordant) / denominator if denominator else np.nan


def _agreement(primary: pd.DataFrame, candidate_count: int) -> pd.DataFrame:
    matrix = primary.pivot_table(
        index="candidate_id", columns="tool_name", values="rank", aggfunc="first"
    )
    rows = []
    tools = sorted(matrix.columns)
    for index, left_tool in enumerate(tools):
        for right_tool in tools[index + 1 :]:
            shared = matrix[[left_tool, right_tool]].dropna()
            left_ranks = shared[left_tool].rank()
            right_ranks = shared[right_tool].rank()
            spearman = (
                left_ranks.corr(right_ranks)
                if len(shared) >= 3 and left_ranks.nunique() > 1 and right_ranks.nunique() > 1
                else np.nan
            )
            kendall = (
                _kendall_tau(shared[left_tool], shared[right_tool]) if len(shared) >= 3 else np.nan
            )
            for top_k in (5, 10, 20):
                effective_k = min(top_k, candidate_count)
                left_top = set(matrix[left_tool].dropna().nsmallest(effective_k).index.astype(str))
                right_top = set(
                    matrix[right_tool].dropna().nsmallest(effective_k).index.astype(str)
                )
                union = left_top | right_top
                intersection = left_top & right_top
                rows.append(
                    {
                        "tool_a": left_tool,
                        "tool_b": right_tool,
                        "shared_candidates": len(shared),
                        "spearman_rank_correlation": spearman,
                        "kendall_rank_correlation": kendall,
                        "top_k": top_k,
                        "top_k_overlap": len(intersection),
                        "jaccard_overlap": len(intersection) / len(union) if union else np.nan,
                    }
                )
    return pd.DataFrame(rows)


def compare_tools(
    candidates: pd.DataFrame | CandidateTable,
    tool_results: list[pd.DataFrame | ToolResultTable],
    *,
    method: str = "weighted_borda",
    weights: dict[str, float] | None = None,
    metric_selection: dict[str, str] | None = None,
    unmatched_external_rows: pd.DataFrame | None = None,
    expected_tools: list[str] | None = None,
) -> ComparisonResult:
    """Compare rank-normalized metrics while retaining missingness and disagreement."""
    candidate_frame = as_dataframe(candidates)
    if method not in {"weighted_borda", "median_rank", "worst_case_rank"}:
        raise ValueError("method must be weighted_borda, median_rank, or worst_case_rank")
    frames = [as_dataframe(item) for item in tool_results]
    long = pd.concat(frames, ignore_index=True) if frames else ToolResultTable.empty().dataframe
    ToolResultTable.from_frame(long)
    primary = _select_primary_metrics(long, metric_selection)
    tools = sorted(set(primary["tool_name"].dropna().unique()) | set(expected_tools or []))
    total_tools = len(tools)
    percentile_matrix = primary.pivot_table(
        index="candidate_id", columns="tool_name", values="percentile_rank", aggfunc="first"
    ).reindex(candidate_frame["candidate_id"])
    rank_matrix = primary.pivot_table(
        index="candidate_id", columns="tool_name", values="rank", aggfunc="first"
    ).reindex(candidate_frame["candidate_id"])
    percentile_matrix = percentile_matrix.reindex(columns=tools)
    rank_matrix = rank_matrix.reindex(columns=tools)
    tool_weights = pd.Series({tool: float((weights or {}).get(tool, 1.0)) for tool in tools})
    rows = []
    for _, candidate in candidate_frame.iterrows():
        candidate_id = candidate["candidate_id"]
        percentiles = percentile_matrix.loc[candidate_id] if tools else pd.Series(dtype=float)
        ranks = rank_matrix.loc[candidate_id] if tools else pd.Series(dtype=float)
        reported = percentiles.dropna()
        reporting_count = len(reported)
        coverage = reporting_count / total_tools if total_tools else 0.0
        if reporting_count:
            if method == "weighted_borda":
                active_weights = tool_weights.loc[reported.index]
                aggregate = float((reported * active_weights).sum() / active_weights.sum())
                score = aggregate * coverage
            elif method == "median_rank":
                median = float(ranks.dropna().median())
                aggregate = median
                score = (1.0 / (1.0 + median)) * coverage
            else:
                worst = float(ranks.dropna().max())
                aggregate = worst
                score = (1.0 / (1.0 + worst)) * coverage
            variance = float(reported.var(ddof=0))
            disagreement = float((reported.max() - reported.min() + np.sqrt(variance)) / 2.0)
        else:
            aggregate = np.nan
            score = np.nan
            variance = np.nan
            disagreement = np.nan
        rows.append(
            {
                "candidate_id": candidate_id,
                "guide_sequence": candidate.get("guide_sequence", ""),
                "gene_name": candidate.get("gene_name", ""),
                "tools_reporting": reporting_count,
                "tools_missing": total_tools - reporting_count,
                "tool_ranks": json.dumps(
                    {
                        tool: (None if pd.isna(value) else float(value))
                        for tool, value in ranks.items()
                    },
                    sort_keys=True,
                ),
                "aggregation_method": method,
                "aggregation_value": aggregate,
                "consensus_score": score,
                "rank_variance": variance,
                "disagreement_score": disagreement,
                "confidence_coverage_level": (
                    "high" if coverage >= 0.8 else "moderate" if coverage >= 0.5 else "low"
                ),
                "tool_coverage_fraction": coverage,
                "explanation": (
                    f"{reporting_count}/{total_tools} tools reported; ranks/percentiles were "
                    "aggregated, never raw incomparable scores."
                ),
            }
        )
    consensus = pd.DataFrame(rows).sort_values(
        ["consensus_score", "tool_coverage_fraction", "candidate_id"],
        ascending=[False, False, True],
        na_position="last",
        kind="mergesort",
    )
    consensus.insert(1, "consensus_rank", range(1, len(consensus) + 1))
    matrix = candidate_frame[["candidate_id", "gene_name"]].merge(
        percentile_matrix.add_prefix("percentile__").reset_index(), on="candidate_id", how="left"
    )
    rank_export = rank_matrix.add_prefix("rank__").reset_index()
    matrix = matrix.merge(rank_export, on="candidate_id", how="left")
    coverage_rows = []
    for tool in tools:
        reported = int(percentile_matrix[tool].notna().sum())
        coverage_rows.append(
            {
                "tool_name": tool,
                "candidate_count": len(candidate_frame),
                "reported_candidates": reported,
                "missing_candidates": len(candidate_frame) - reported,
                "coverage_fraction": reported / len(candidate_frame) if len(candidate_frame) else 0,
                "status": (
                    "complete"
                    if reported == len(candidate_frame)
                    else "pending"
                    if reported == 0
                    else "partial"
                ),
            }
        )
    coverage_frame = pd.DataFrame(coverage_rows)
    disagreement = consensus[
        consensus["disagreement_score"].fillna(0).ge(0.20) | consensus["tools_missing"].gt(0)
    ].copy()
    disagreement["disagreement_reason"] = disagreement.apply(
        lambda row: "; ".join(
            part
            for part in [
                "high rank spread" if row["disagreement_score"] >= 0.20 else "",
                "one or more tools missing" if row["tools_missing"] else "",
            ]
            if part
        ),
        axis=1,
    )
    return ComparisonResult(
        tool_results_long=long[TOOL_RESULT_COLUMNS].copy(),
        candidate_tool_matrix=matrix,
        consensus_candidates=consensus,
        tool_coverage=coverage_frame,
        model_agreement=_agreement(primary, len(candidate_frame)),
        disagreement_report=disagreement,
        unmatched_external_rows=(
            unmatched_external_rows.copy()
            if unmatched_external_rows is not None
            else pd.DataFrame(columns=["tool_name", "mapping_status"])
        ),
    )


def build_consensus(
    candidates: pd.DataFrame | CandidateTable,
    tool_results: list[pd.DataFrame | ToolResultTable],
    **kwargs: Any,
) -> pd.DataFrame:
    """Return the consensus candidate table from :func:`compare_tools`."""
    return compare_tools(candidates, tool_results, **kwargs).consensus_candidates
