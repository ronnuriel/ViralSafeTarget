"""Publication-facing multi-tool benchmark with explicit missingness.

The benchmark compares only commensurable within-tool ranks and counts.  A tool
without a committed raw export remains pending; it is never assigned zero hits or
an inferred rank.
"""

from __future__ import annotations

import json
import shutil
from datetime import datetime, timezone
from html import escape
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from .provenance import sha256_file

RESULT_COLUMNS = [
    "candidate_id",
    "guide_sequence",
    "gene_name",
    "tool_name",
    "tool_version",
    "metric_name",
    "metric_direction",
    "raw_value",
    "rank",
    "percentile_rank",
    "status",
    "source_file",
    "source_sha256",
    "notes",
]

TOOL_DISPLAY_NAMES = {
    "viral_safe_target_pre_human": "ViralSafeTarget\npre-host",
    "viral_safe_target_post_human": "ViralSafeTarget\npost-host",
    "cas-offinder": "Cas-OFFinder",
    "crispritz": "CRISPRitz",
    "crispor": "CRISPOR",
    "chopchop": "CHOPCHOP",
    "guidescan2": "GuideScan2",
}


def _resolve(root: Path, value: str | Path | None) -> Path | None:
    if value in (None, ""):
        return None
    path = Path(value)
    return path if path.is_absolute() else root / path


def _portable_path(path: str | Path, root: Path) -> str:
    """Prefer repository-relative provenance while retaining external paths."""
    resolved = Path(path).resolve()
    try:
        return str(resolved.relative_to(root.resolve()))
    except ValueError:
        return str(resolved)


def _rank(values: pd.Series, direction: str) -> tuple[pd.Series, pd.Series]:
    numeric = pd.to_numeric(values, errors="coerce")
    ranks = numeric.rank(method="min", ascending=direction == "lower", na_option="keep")
    count = int(numeric.notna().sum())
    if count <= 1:
        percentile = pd.Series(np.where(numeric.notna(), 1.0, np.nan), index=values.index)
    else:
        percentile = 1.0 - (ranks - 1.0) / (count - 1.0)
    return ranks, percentile


def _metric_rows(
    candidates: pd.DataFrame,
    *,
    tool_name: str,
    tool_version: str,
    metric_name: str,
    metric_direction: str,
    values: pd.Series,
    source_file: Path,
    notes: str,
) -> pd.DataFrame:
    ranks, percentiles = _rank(values, metric_direction)
    output = pd.DataFrame(
        {
            "candidate_id": candidates["candidate_id"].astype(str),
            "guide_sequence": candidates["guide_sequence"].astype(str),
            "gene_name": candidates.get("gene_name", pd.Series("", index=candidates.index)),
            "tool_name": tool_name,
            "tool_version": tool_version,
            "metric_name": metric_name,
            "metric_direction": metric_direction,
            "raw_value": pd.to_numeric(values, errors="coerce"),
            "rank": ranks,
            "percentile_rank": percentiles,
            "status": np.where(
                pd.to_numeric(values, errors="coerce").notna(), "completed", "missing"
            ),
            "source_file": str(source_file),
            "source_sha256": sha256_file(source_file),
            "notes": notes,
        }
    )
    return output[RESULT_COLUMNS]


def parse_crispritz_profile(
    profile_path: str | Path, candidates: pd.DataFrame
) -> pd.DataFrame:
    """Parse CRISPRitz ``*.profile.xls`` counts and map them by guide sequence."""
    path = Path(profile_path)
    profile = pd.read_csv(path, sep="\t")
    guide_column = next((column for column in profile if str(column).upper() == "GUIDE"), None)
    if guide_column is None:
        raise ValueError("CRISPRitz profile has no GUIDE column")
    mismatch_columns = [column for column in profile if str(column).endswith("MM")]
    if not mismatch_columns:
        raise ValueError("CRISPRitz profile has no mismatch-count columns")
    profile = profile.copy()
    profile["guide_sequence"] = (
        profile[guide_column].astype(str).str.upper().str.replace("-", "", regex=False).str[:20]
    )
    for column in mismatch_columns:
        profile[column] = pd.to_numeric(profile[column], errors="coerce")
    profile["predicted_offtarget_burden"] = profile[mismatch_columns].sum(axis=1, min_count=1)
    if profile["guide_sequence"].duplicated().any():
        raise ValueError("CRISPRitz profile contains duplicate guide rows")
    mapping = profile.set_index("guide_sequence")["predicted_offtarget_burden"]
    output = candidates[["candidate_id", "guide_sequence"]].copy()
    output["predicted_offtarget_burden"] = output["guide_sequence"].str.upper().map(mapping)
    return output


def freeze_panel(candidates: pd.DataFrame, expected_count: int | None = None) -> pd.DataFrame:
    required = {"candidate_id", "guide_sequence"}
    missing = sorted(required - set(candidates.columns))
    if missing:
        raise ValueError(f"Candidate panel is missing required columns: {missing}")
    if candidates["candidate_id"].duplicated().any():
        raise ValueError("Frozen benchmark candidate IDs must be unique")
    if candidates["guide_sequence"].astype(str).str.upper().duplicated().any():
        raise ValueError("Frozen benchmark guide sequences must be unique")
    frozen = candidates.sort_values("candidate_id", kind="mergesort").reset_index(drop=True)
    if expected_count is not None and len(frozen) != expected_count:
        raise ValueError(f"Expected {expected_count} frozen candidates, observed {len(frozen)}")
    return frozen


def write_tool_inputs(candidates: pd.DataFrame, output_directory: str | Path) -> dict[str, Path]:
    """Write documented, tool-neutral input bundles without invoking web services."""
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    identity = [
        column
        for column in (
            "candidate_id",
            "guide_sequence",
            "pam",
            "gene_name",
            "reference_start_1based",
            "strand",
        )
        if column in candidates
    ]
    paths: dict[str, Path] = {}
    paths["common"] = output / "frozen_panel.csv"
    candidates.to_csv(paths["common"], index=False)
    paths["identity"] = output / "candidate_identity.tsv"
    candidates[identity].to_csv(paths["identity"], sep="\t", index=False)
    paths["crispritz"] = output / "crispritz_guides.txt"
    paths["crispritz"].write_text(
        "".join(f"{sequence.upper()}NNN\n" for sequence in candidates["guide_sequence"]),
        encoding="utf-8",
    )
    paths["guidescan2"] = output / "guidescan2_guides.txt"
    paths["guidescan2"].write_text(
        "".join(f"{sequence.upper()}\n" for sequence in candidates["guide_sequence"]),
        encoding="utf-8",
    )
    paths["crispor"] = output / "crispor_batch.tsv"
    candidates[["candidate_id", "guide_sequence"]].to_csv(
        paths["crispor"], sep="\t", index=False
    )
    paths["chopchop"] = output / "chopchop_targets.fasta"
    paths["chopchop"].write_text(
        "".join(
            f">{row.candidate_id}\n"
            f"{str(row.guide_sequence).upper()}{str(getattr(row, 'pam', ''))}\n"
            for row in candidates.itertuples(index=False)
        ),
        encoding="utf-8",
    )
    manifest = {
        "schema_version": "1.0",
        "candidate_count": len(candidates),
        "guide_sequence_length": [
            int(value)
            for value in sorted(candidates["guide_sequence"].astype(str).str.len().unique())
        ],
        "files": {
            label: {"path": path.name, "sha256": sha256_file(path)}
            for label, path in paths.items()
        },
        "interpretation": (
            "Input preparation only. A tool is completed only after its raw output is committed "
            "and normalized; missing output remains pending."
        ),
    }
    paths["manifest"] = output / "input_manifest.json"
    paths["manifest"].write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return paths


def _rank_agreement(results: pd.DataFrame, top_k: list[int]) -> tuple[pd.DataFrame, pd.DataFrame]:
    primary = results[results["status"].eq("completed")].copy()
    matrix = primary.pivot_table(
        index="candidate_id", columns="tool_name", values="rank", aggfunc="first"
    )
    agreement: list[dict[str, object]] = []
    overlaps: list[dict[str, object]] = []
    tools = sorted(matrix.columns)
    for index, left in enumerate(tools):
        for right in tools[index + 1 :]:
            shared = matrix[[left, right]].dropna()
            correlation = (
                shared[left].rank().corr(shared[right].rank(), method="pearson")
                if len(shared) >= 3
                and shared[left].nunique() > 1
                and shared[right].nunique() > 1
                else np.nan
            )
            agreement.append(
                {
                    "tool_a": left,
                    "tool_b": right,
                    "shared_candidates": len(shared),
                    "spearman_rank_correlation": correlation,
                    "status": "completed" if len(shared) else "not_comparable",
                }
            )
            for k in top_k:
                left_count = int(matrix[left].notna().sum())
                left_top = set(matrix[left].dropna().nsmallest(min(k, left_count)).index)
                right_top = set(
                    matrix[right].dropna().nsmallest(min(k, matrix[right].notna().sum())).index
                )
                union = left_top | right_top
                overlaps.append(
                    {
                        "tool_a": left,
                        "tool_b": right,
                        "top_k": k,
                        "overlap_count": len(left_top & right_top),
                        "jaccard_overlap": (
                            len(left_top & right_top) / len(union) if union else np.nan
                        ),
                    }
                )
    return pd.DataFrame(agreement), pd.DataFrame(overlaps)


def run_ablation(
    candidates: pd.DataFrame,
    components: list[dict[str, object]],
    penalty_columns: list[str],
    top_k: list[int],
    rank_precision_decimals: int = 12,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Leave-one-component-out ranking; this is not model retraining."""
    component_names = [str(item["name"]) for item in components]
    variants = [("all_components", None), *[(f"without_{name}", name) for name in component_names]]
    penalties = sum(
        (
            pd.to_numeric(candidates.get(column, 0), errors="coerce").fillna(0)
            for column in penalty_columns
        ),
        start=pd.Series(0.0, index=candidates.index),
    ).clip(0, 1)
    detail_rows: list[pd.DataFrame] = []
    for variant, omitted in variants:
        active = [item for item in components if item["name"] != omitted]
        denominator = sum(float(item["weight"]) for item in active)
        score = pd.Series(0.0, index=candidates.index)
        observed_weight = pd.Series(0.0, index=candidates.index)
        for item in active:
            numeric = pd.to_numeric(candidates[str(item["column"])], errors="coerce")
            weight = float(item["weight"])
            score += numeric.fillna(0) * weight
            observed_weight += numeric.notna().astype(float) * weight
        score = ((score / denominator) * (1 - penalties)).clip(0, 1)
        order = pd.DataFrame(
            {
                "candidate_id": candidates["candidate_id"].astype(str),
                "score": score.round(rank_precision_decimals),
                "_index": candidates.index,
            }
        ).sort_values(
            ["score", "candidate_id"], ascending=[False, True], kind="mergesort"
        )
        rank = pd.Series(index=candidates.index, dtype=float)
        rank.loc[order["_index"]] = range(1, len(order) + 1)
        detail_rows.append(
            pd.DataFrame(
                {
                    "candidate_id": candidates["candidate_id"],
                    "variant": variant,
                    "omitted_component": omitted or "none",
                    "score": score,
                    "rank": rank,
                    "observed_weight_fraction": observed_weight / denominator,
                }
            )
        )
    detail = pd.concat(detail_rows, ignore_index=True)
    baseline = detail[detail["variant"].eq("all_components")].set_index("candidate_id")["rank"]
    summaries: list[dict[str, object]] = []
    for variant, group in detail.groupby("variant", sort=False):
        ranks = group.set_index("candidate_id")["rank"]
        shifts = (ranks - baseline).abs()
        row: dict[str, object] = {
            "variant": variant,
            "spearman_vs_all_components": ranks.corr(baseline, method="pearson"),
            "median_absolute_rank_shift": float(shifts.median()),
            "maximum_absolute_rank_shift": float(shifts.max()),
        }
        for k in top_k:
            base_top = set(baseline.nsmallest(k).index)
            variant_top = set(ranks.nsmallest(k).index)
            row[f"top_{k}_overlap"] = len(base_top & variant_top)
        summaries.append(row)
    return detail, pd.DataFrame(summaries)


def _status_table(
    tools: list[dict[str, object]], results: pd.DataFrame, candidate_count: int
) -> pd.DataFrame:
    rows = []
    for tool in tools:
        name = str(tool["id"])
        tool_results = results[results["tool_name"].eq(name)]
        completed = tool_results["status"].eq("completed")
        reported = int(tool_results.loc[completed, "candidate_id"].nunique())
        configured = str(tool.get("status", "pending"))
        if reported == candidate_count:
            status = "completed"
        elif reported:
            status = "partial"
        else:
            status = configured
        rows.append(
            {
                "tool_name": name,
                "tool_version": str(tool.get("version", "not assessed")),
                "execution_mode": str(tool.get("mode", "not assessed")),
                "status": status,
                "candidate_count": candidate_count,
                "reported_candidates": reported,
                "missing_candidates": candidate_count - reported,
                "coverage_fraction": reported / candidate_count if candidate_count else np.nan,
                "runtime_seconds": tool.get("runtime_seconds", pd.NA),
                "runtime_scope": str(tool.get("runtime_scope", "not assessed")),
                "reason": str(tool.get("reason", "")),
                "official_source": str(tool.get("official_source", "")),
            }
        )
    return pd.DataFrame(rows)


def _write_figures(
    status: pd.DataFrame, agreement: pd.DataFrame, ablation: pd.DataFrame, output: Path
) -> None:
    figures = output / "figures"
    figures.mkdir(exist_ok=True)
    colors = status["status"].map(
        {
            "completed": "#2a9d8f",
            "partial": "#e9c46a",
            "export_required": "#f4a261",
            "pending": "#b0b0b0",
        }
    ).fillna("#b0b0b0")
    labels = [TOOL_DISPLAY_NAMES.get(name, name) for name in status["tool_name"]]
    fig, axis = plt.subplots(figsize=(10, 5), constrained_layout=True)
    axis.bar(status["tool_name"], status["coverage_fraction"].fillna(0), color=colors)
    axis.set_ylim(0, 1.05)
    axis.set_ylabel("Reported-metric coverage of frozen panel")
    axis.set_title("Executable benchmark coverage; missing output is not zero")
    axis.set_xticks(range(len(labels)), labels, rotation=25, ha="right")
    fig.savefig(figures / "tool_result_coverage.png", dpi=180)
    plt.close(fig)

    if not agreement.empty:
        tools = sorted(set(agreement["tool_a"]) | set(agreement["tool_b"]))
        matrix = pd.DataFrame(np.nan, index=tools, columns=tools, dtype=float)
        for tool in tools:
            matrix.loc[tool, tool] = 1.0
        for row in agreement.itertuples(index=False):
            matrix.loc[row.tool_a, row.tool_b] = row.spearman_rank_correlation
            matrix.loc[row.tool_b, row.tool_a] = row.spearman_rank_correlation
        fig, axis = plt.subplots(figsize=(7.5, 6), constrained_layout=True)
        colour_map = plt.get_cmap("coolwarm").copy()
        colour_map.set_bad("#d9d9d9")
        image = axis.imshow(
            np.ma.masked_invalid(matrix.to_numpy(dtype=float)),
            vmin=-1,
            vmax=1,
            cmap=colour_map,
        )
        display_labels = [TOOL_DISPLAY_NAMES.get(name, name) for name in tools]
        axis.set_xticks(range(len(matrix.columns)), display_labels, rotation=30, ha="right")
        axis.set_yticks(range(len(matrix.index)), display_labels)
        axis.set_title("Pairwise within-tool rank agreement")
        fig.colorbar(image, ax=axis, label="Spearman rank correlation")
        fig.savefig(figures / "rank_agreement.png", dpi=180)
        plt.close(fig)

    plot = ablation[~ablation["variant"].eq("all_components")].copy()
    fig, axis = plt.subplots(figsize=(9, 4.5))
    labels = plot["variant"].str.replace("without_", "", regex=False)
    axis.bar(labels, plot["maximum_absolute_rank_shift"])
    axis.set_ylabel("Maximum absolute rank shift")
    axis.set_title("Leave-one-component-out sensitivity on the frozen panel")
    axis.tick_params(axis="x", rotation=30)
    fig.tight_layout()
    fig.savefig(figures / "ablation_rank_shift.png", dpi=180)
    plt.close(fig)


def _html_table(frame: pd.DataFrame, limit: int | None = None) -> str:
    visible = frame.head(limit) if limit else frame
    return visible.to_html(index=False, border=0, na_rep="unknown", escape=True)


def _write_report(
    output: Path,
    status: pd.DataFrame,
    agreement: pd.DataFrame,
    overlap: pd.DataFrame,
    ablation: pd.DataFrame,
    capabilities: pd.DataFrame,
    candidate_count: int,
) -> None:
    completed = status[status["status"].eq("completed")]["tool_name"].tolist()
    pending = status[~status["status"].eq("completed")]["tool_name"].tolist()
    css = """
body {font-family: Arial, sans-serif; max-width: 1180px; margin: 2rem auto;
      line-height: 1.5; color: #17202a}
table {border-collapse: collapse; width: 100%; font-size: .88rem}
th, td {border: 1px solid #ddd; padding: .45rem; text-align: left}
th {background: #eef3f7}
.warning {background: #fff4d6; padding: 1rem; border-left: 5px solid #e9a23b}
img {max-width: 100%}
"""
    completed_text = escape(", ".join(completed) or "none")
    pending_text = escape(", ".join(pending) or "none")
    sections = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        "<title>ViralSafeTarget multi-tool benchmark</title>",
        f"<style>{css}</style></head><body>",
        f"<h1>Multi-tool benchmark on a frozen {candidate_count}-guide panel</h1>",
        "<div class='warning'><strong>Interpretation boundary.</strong> Tool scores "
        "are not placed on a common biological scale. Missing output remains pending, "
        "not zero. This benchmark does not establish editing efficacy, safety, viral "
        "inactivation, treatment, or cure.</div>",
        "<h2>Execution status</h2>",
        f"<p>Completed: {completed_text}. Pending/export-required: {pending_text}.</p>",
        _html_table(status),
        "<img src='figures/tool_result_coverage.png' alt='Tool result coverage'>",
        "<h2>Rank agreement</h2>",
        _html_table(agreement),
        "<img src='figures/rank_agreement.png' alt='Rank agreement'>",
        "<h2>Top-K overlap</h2>",
        _html_table(overlap, 60),
        "<h2>ViralSafeTarget ablation</h2>",
        "<p>This is sensitivity analysis, not model retraining or biological validation.</p>",
        _html_table(ablation),
        "<img src='figures/ablation_rank_shift.png' alt='Ablation rank shifts'>",
        "<h2>Capability evidence matrix</h2>",
        "<p>Documented scope is separate from executable performance.</p>",
        _html_table(capabilities),
        "<h2>What this benchmark can support</h2>",
        "<p>It supports computational coverage, agreement, missingness, and reproducible "
        "integration claims. Predictive superiority requires experimental ground truth.</p>",
        "</body></html>",
    ]
    body = "\n".join(sections)
    (output / "multitool_benchmark_report.html").write_text(body, encoding="utf-8")


def _write_findings(
    output: Path,
    status: pd.DataFrame,
    agreement: pd.DataFrame,
    overlap: pd.DataFrame,
    ablation: pd.DataFrame,
    candidate_count: int,
) -> None:
    completed = status.loc[status["status"].eq("completed"), "tool_name"].tolist()
    incomplete = status.loc[~status["status"].eq("completed"), "tool_name"].tolist()
    strongest = ablation.loc[
        ablation["maximum_absolute_rank_shift"].idxmax()
    ].to_dict()
    correlation_lines = []
    for row in agreement.itertuples(index=False):
        value = (
            "not comparable"
            if pd.isna(row.spearman_rank_correlation)
            else f"{float(row.spearman_rank_correlation):.6f}"
        )
        correlation_lines.append(f"- `{row.tool_a}` versus `{row.tool_b}`: {value}.")
    top_lines = []
    for row in overlap[overlap["top_k"].isin([10, 25, 50])].itertuples(index=False):
        top_lines.append(
            f"- `{row.tool_a}` versus `{row.tool_b}`, K={row.top_k}: "
            f"{row.overlap_count} shared guides (Jaccard {row.jaccard_overlap:.6f})."
        )
    text = "\n".join(
        [
            "# Frozen-panel benchmark findings",
            "",
            "## Computational observations",
            "",
            f"- The identity-frozen benchmark contains {candidate_count} unique guides.",
            f"- Completed primary metrics: {', '.join(completed)}.",
            f"- Incomplete or export-required tools: {', '.join(incomplete)}.",
            "- Missing results remain missing; no unavailable output was interpreted as zero.",
            *correlation_lines,
            *top_lines,
            (
                "- The largest leave-one-component-out shift was observed for "
                f"`{strongest['variant']}`: median absolute shift "
                f"{strongest['median_absolute_rank_shift']:.1f}, maximum "
                f"{strongest['maximum_absolute_rank_shift']:.1f}."
            ),
            "",
            "## Research hypothesis",
            "",
            "The compared systems expose complementary rather than interchangeable axes. "
            "A virus-first workflow may add value by retaining viral-population support, "
            "host-risk results, gene/protein context, evidence provenance, and multiplex "
            "escape analysis in one auditable record.",
            "",
            "## Evidence gaps",
            "",
            "- CRISPOR, CHOPCHOP, and GuideScan2 require committed raw exports before "
            "quantitative rank comparison.",
            "- A completed second-virus benchmark is still required for generalization.",
            "- Independent experimental ground truth is required to test predictive superiority.",
            "",
            "## Limitations",
            "",
            "- Tool raw scores are not on a common biological scale and were not averaged.",
            "- Runtime is not compared across different guide counts, assemblies, hardware, or "
            "search models.",
            "- Ablation is conditioned on an enriched deep-screening panel and is not model "
            "retraining.",
            "- Capability evidence describes documented scope, not executable performance.",
            "- No editing, safety, efficacy, viral inactivation, treatment, or cure claim is made.",
            "",
        ]
    )
    (output / "FINDINGS.md").write_text(text, encoding="utf-8")


def run_tool_benchmark(config_path: str | Path) -> dict[str, object]:
    """Run a generic, config-driven benchmark and write auditable artifacts."""
    config_file = Path(config_path).resolve()
    settings = yaml.safe_load(config_file.read_text(encoding="utf-8")) or {}
    root = _resolve(config_file.parent, settings.get("project_root", "."))
    assert root is not None
    root = root.resolve()
    candidate_path = _resolve(root, settings["candidate_table"])
    output = _resolve(root, settings["output_dir"])
    capability_path = _resolve(root, settings["capability_evidence_table"])
    ranking_config_path = _resolve(root, settings["ablation"]["ranking_config"])
    assert candidate_path and output and capability_path and ranking_config_path
    if not candidate_path.is_file():
        raise FileNotFoundError(candidate_path)
    output.mkdir(parents=True, exist_ok=True)
    candidates = freeze_panel(
        pd.read_csv(candidate_path), int(settings["expected_candidate_count"])
    )
    inputs = write_tool_inputs(candidates, output / "inputs")

    result_frames: list[pd.DataFrame] = []
    for tool in settings["tools"]:
        tool_id = str(tool["id"])
        source = _resolve(root, tool.get("result_file")) or candidate_path
        if tool_id == "viral_safe_target_pre_human":
            result_frames.append(
                _metric_rows(
                    candidates,
                    tool_name=tool_id,
                    tool_version=str(tool["version"]),
                    metric_name="pre_human_score",
                    metric_direction="higher",
                    values=candidates["pre_human_score"],
                    source_file=candidate_path,
                    notes="Internal sequence/annotation score before host-risk adjustment.",
                )
            )
        elif tool_id == "viral_safe_target_post_human":
            result_frames.append(
                _metric_rows(
                    candidates,
                    tool_name=tool_id,
                    tool_version=str(tool["version"]),
                    metric_name="post_human_score",
                    metric_direction="higher",
                    values=candidates["post_human_score"],
                    source_file=candidate_path,
                    notes="Internal score after the configured Cas-OFFinder risk layer.",
                )
            )
        elif tool_id == "cas-offinder":
            result_frames.append(
                _metric_rows(
                    candidates,
                    tool_name=tool_id,
                    tool_version=str(tool["version"]),
                    metric_name="predicted_offtarget_burden",
                    metric_direction="lower",
                    values=candidates["human_total_predicted_hits"],
                    source_file=candidate_path,
                    notes=(
                        "Reference-genome similarities through the configured mismatch threshold."
                    ),
                )
            )
        elif tool_id == "crispritz" and source.is_file():
            parsed = parse_crispritz_profile(source, candidates)
            result_frames.append(
                _metric_rows(
                    candidates,
                    tool_name=tool_id,
                    tool_version=str(tool["version"]),
                    metric_name="predicted_offtarget_burden",
                    metric_direction="lower",
                    values=parsed["predicted_offtarget_burden"],
                    source_file=source,
                    notes=(
                        "CRISPRitz reference-genome mismatch profile; bulges/variants only if "
                        "configured."
                    ),
                )
            )
        elif source.is_file() and source != candidate_path:
            imported = pd.read_csv(source)
            missing = sorted(set(RESULT_COLUMNS) - set(imported.columns))
            if missing:
                raise ValueError(f"Normalized results for {tool_id} are missing columns: {missing}")
            result_frames.append(imported[RESULT_COLUMNS])
    results = (
        pd.concat(result_frames, ignore_index=True)
        if result_frames
        else pd.DataFrame(columns=RESULT_COLUMNS)
    )
    results["source_file"] = results["source_file"].map(
        lambda value: _portable_path(value, root)
    )
    results.to_csv(output / "tool_results_long.csv", index=False)
    rank_matrix = results.pivot_table(
        index="candidate_id", columns="tool_name", values="rank", aggfunc="first"
    ).reindex(candidates["candidate_id"])
    rank_matrix.reset_index().to_csv(output / "candidate_rank_matrix.csv", index=False)
    metric_summary = (
        results.groupby(["tool_name", "metric_name"], dropna=False)
        .agg(
            reported_candidates=("candidate_id", "nunique"),
            numeric_values=("raw_value", "count"),
            unique_numeric_values=("raw_value", "nunique"),
            minimum_raw_value=("raw_value", "min"),
            maximum_raw_value=("raw_value", "max"),
        )
        .reset_index()
    )
    metric_summary.to_csv(output / "tool_metric_summary.csv", index=False)

    top_k = [int(value) for value in settings.get("top_k", [10, 25, 50])]
    agreement, overlap = _rank_agreement(results, top_k)
    agreement.to_csv(output / "rank_agreement.csv", index=False)
    overlap.to_csv(output / "top_k_overlap.csv", index=False)
    status = _status_table(settings["tools"], results, len(candidates))
    status.to_csv(output / "tool_execution_status.csv", index=False)

    ranking_settings = yaml.safe_load(ranking_config_path.read_text(encoding="utf-8"))
    weights = ranking_settings["ranking"]["weights"]
    components = [
        {"name": name, "column": column, "weight": weights[name]}
        for name, column in settings["ablation"]["components"].items()
    ]
    ablation_detail, ablation_summary = run_ablation(
        candidates,
        components,
        list(settings["ablation"].get("penalty_columns", [])),
        top_k,
        int(settings["ablation"].get("rank_precision_decimals", 12)),
    )
    baseline_column = settings["ablation"].get("baseline_score_column")
    if baseline_column:
        baseline = ablation_detail[ablation_detail["variant"].eq("all_components")].set_index(
            "candidate_id"
        )["score"]
        expected = candidates.set_index("candidate_id")[str(baseline_column)]
        maximum_difference = float((baseline - expected).abs().max())
        tolerance = float(settings["ablation"].get("baseline_tolerance", 1e-12))
        if maximum_difference > tolerance:
            raise ValueError(
                "Ablation baseline does not reproduce the configured source score: "
                f"maximum difference {maximum_difference} exceeds {tolerance}"
            )
    ablation_detail.to_csv(output / "ablation_candidate_ranks.csv", index=False)
    ablation_summary.to_csv(output / "ablation_summary.csv", index=False)

    capabilities = pd.read_csv(capability_path, sep="\t")
    allowed = {"yes", "partial", "no", "not_assessed"}
    invalid = sorted(set(capabilities["status"]) - allowed)
    if invalid:
        raise ValueError(f"Capability table contains invalid statuses: {invalid}")
    capabilities.to_csv(output / "tool_capability_evidence.csv", index=False)
    capabilities.pivot(
        index="capability", columns="tool_name", values="status"
    ).reset_index().to_csv(output / "capability_matrix.csv", index=False)
    _write_figures(status, agreement, ablation_summary, output)
    _write_report(
        output, status, agreement, overlap, ablation_summary, capabilities, len(candidates)
    )
    _write_findings(
        output, status, agreement, overlap, ablation_summary, len(candidates)
    )

    source_files = [candidate_path, config_file, capability_path, ranking_config_path]
    source_files.extend(
        path
        for tool in settings["tools"]
        if (path := _resolve(root, tool.get("result_file"))) is not None and path.is_file()
    )
    manifest = {
        "schema_version": "1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "benchmark_id": settings["benchmark_id"],
        "candidate_count": len(candidates),
        "unique_guide_count": int(candidates["guide_sequence"].nunique()),
        "top_k": top_k,
        "source_files": [
            {"path": _portable_path(path, root), "sha256": sha256_file(path)}
            for path in dict.fromkeys(source_files)
        ],
        "input_manifest_sha256": sha256_file(inputs["manifest"]),
        "completed_tools": status.loc[status["status"].eq("completed"), "tool_name"].tolist(),
        "incomplete_tools": status.loc[~status["status"].eq("completed"), "tool_name"].tolist(),
        "raw_scores_directly_averaged": False,
        "ablation_baseline_maximum_absolute_difference": (
            maximum_difference if baseline_column else None
        ),
        "ablation_rank_precision_decimals": int(
            settings["ablation"].get("rank_precision_decimals", 12)
        ),
        "missing_output_interpretation": "pending_or_export_required; never zero",
        "claims_boundary": (
            "Computational coverage/agreement/sensitivity only; no safety, efficacy, viral "
            "inactivation, treatment, cure, or repair-frequency claim."
        ),
    }
    (output / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    shutil.copy2(config_file, output / "benchmark_config.yaml")
    return {
        "output_directory": output,
        "report": output / "multitool_benchmark_report.html",
        "candidate_count": len(candidates),
        "completed_tools": manifest["completed_tools"],
        "incomplete_tools": manifest["incomplete_tools"],
    }
