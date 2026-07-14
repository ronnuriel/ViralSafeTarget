"""Project-facing orchestration for virtual knockout and escape analyses."""

from __future__ import annotations

import html
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .escape import (
    select_strategy_panels,
    summarize_guide_escape,
    summarize_multiplex_strategies,
)
from .project_workflow import ProjectContext, load_project
from .provenance import sha256_file
from .virtual_knockout import (
    build_cds_models,
    enumerate_indel_hypotheses,
    load_annotation_regions,
    map_guides_to_cds,
    summarize_virtual_knockout,
)

ANALYSIS_VERSION = "1.0"
OUTPUT_NAMES = {
    "mapping": "guide_cds_mapping.csv",
    "hypotheses": "indel_sequence_hypotheses.csv",
    "guide_virtual": "guide_virtual_knockout.csv",
    "gene_virtual": "gene_virtual_knockout.csv",
    "counterfactuals": "single_nt_escape_counterfactuals.csv",
    "guide_escape": "guide_escape_robustness.csv",
    "members": "multiplex_panel_members.csv",
    "multiplex": "multiplex_escape_robustness.csv",
    "comparison": "strategy_comparison.csv",
    "report": "virtual_knockout_escape_report.html",
    "findings": "FINDINGS.md",
    "manifest": "run_manifest.json",
}


def _analysis(context: ProjectContext) -> dict[str, Any]:
    return dict(context.values.get("analysis") or {})


def _resolve(context: ProjectContext, value: str | None) -> Path | None:
    return context.resolve(value) if value else None


def _display_path(context: ProjectContext, path: Path) -> str:
    try:
        return str(path.resolve().relative_to(context.root))
    except ValueError:
        return str(path.resolve())


def analysis_output_dir(context: ProjectContext) -> Path:
    configured = str(_analysis(context).get("output_dir", "virtual_knockout_escape"))
    path = Path(configured)
    return path if path.is_absolute() else context.output_root / path


def _candidate_path(context: ProjectContext) -> Path:
    configured = _resolve(context, _analysis(context).get("candidate_table"))
    if configured:
        return configured
    post_host = context.output_root / "host_screen" / "candidates_ranked_post_host.csv"
    return (
        post_host if post_host.is_file() else context.output_root / "discovery/discovery_panel.csv"
    )


def _read_optional(path: Path | None, *, sep: str = ",") -> pd.DataFrame:
    return pd.read_csv(path, sep=sep) if path and path.is_file() else pd.DataFrame()


def _required_inputs(context: ProjectContext) -> tuple[Path, Path, Path]:
    reference = context.profiles.resolve(context.profiles.virus.get("reference_fasta"))
    gff = context.profiles.resolve(context.profiles.virus.get("annotation_gff"))
    candidates = _candidate_path(context)
    missing = [str(path) for path in (reference, gff, candidates) if not path or not path.is_file()]
    if missing:
        raise FileNotFoundError("Virtual analysis inputs are missing: " + ", ".join(missing))
    assert reference is not None and gff is not None
    return reference, gff, candidates


def run_virtual_knockout_analysis(project: str | Path | ProjectContext) -> dict[str, Any]:
    context = load_project(project) if not isinstance(project, ProjectContext) else project
    settings = _analysis(context)
    reference, gff, candidate_path = _required_inputs(context)
    output = analysis_output_dir(context)
    output.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(candidate_path, low_memory=False)
    minimum = int(settings.get("indel_min_bp", -10))
    maximum = int(settings.get("indel_max_bp", 10))
    if minimum > 0 or maximum < 0 or maximum - minimum > 200:
        raise ValueError("Indel bounds must span zero and contain at most 201 integer sizes")
    domains = context.profiles.resolve(context.profiles.virus.get("domain_table"))
    disorder = context.profiles.resolve(context.profiles.virus.get("disorder_table"))
    conserved = context.profiles.resolve(context.profiles.virus.get("conserved_region_table"))
    annotations = load_annotation_regions(domains, disorder, conserved)
    cds_models = build_cds_models(
        reference,
        gff,
        str(context.profiles.virus["reference_accession"]),
    )
    mapped = map_guides_to_cds(candidates, cds_models, context.profiles.editor, annotations)
    hypotheses = enumerate_indel_hypotheses(mapped, range(minimum, maximum + 1))
    guide_summary, gene_summary = summarize_virtual_knockout(mapped, hypotheses)
    mapped.drop(columns=["_cds_sequence"], errors="ignore").to_csv(
        output / OUTPUT_NAMES["mapping"], index=False
    )
    hypotheses.to_csv(output / OUTPUT_NAMES["hypotheses"], index=False)
    guide_summary.to_csv(output / OUTPUT_NAMES["guide_virtual"], index=False)
    gene_summary.to_csv(output / OUTPUT_NAMES["gene_virtual"], index=False)
    return {
        "output_dir": output,
        "candidate_count": candidates["candidate_id"].nunique(),
        "mapping_row_count": len(mapped),
        "mapped_candidate_count": guide_summary[
            guide_summary["mapping_status"].eq("mapped_to_cds")
        ]["candidate_id"].nunique(),
        "hypothesis_count": len(hypotheses),
        "gene_count": gene_summary["gene_name"].nunique() if not gene_summary.empty else 0,
    }


def run_escape_analysis(project: str | Path | ProjectContext) -> dict[str, Any]:
    context = load_project(project) if not isinstance(project, ProjectContext) else project
    _, _, candidate_path = _required_inputs(context)
    output = analysis_output_dir(context)
    output.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(candidate_path, low_memory=False)
    heldout_path = _resolve(context, _analysis(context).get("heldout_population_table"))
    heldout = _read_optional(heldout_path)
    guide_escape, counterfactuals = summarize_guide_escape(
        candidates, context.profiles.editor, heldout
    )
    guide_escape.to_csv(output / OUTPUT_NAMES["guide_escape"], index=False)
    counterfactuals.to_csv(output / OUTPUT_NAMES["counterfactuals"], index=False)
    return {
        "output_dir": output,
        "guide_count": len(guide_escape),
        "counterfactual_count": len(counterfactuals),
        "heldout_matched_count": int(
            pd.to_numeric(guide_escape["heldout_exact_target_coverage"], errors="coerce")
            .notna()
            .sum()
        ),
    }


def run_multiplex_analysis(project: str | Path | ProjectContext) -> dict[str, Any]:
    context = load_project(project) if not isinstance(project, ProjectContext) else project
    _, _, candidate_path = _required_inputs(context)
    output = analysis_output_dir(context)
    settings = _analysis(context)
    candidates = pd.read_csv(candidate_path, low_memory=False)
    categories_path = context.profiles.resolve(context.profiles.virus.get("gene_category_table"))
    categories = _read_optional(categories_path, sep="\t")
    definitions = list(settings.get("strategies") or [])
    if not definitions:
        definitions = [
            {
                "id": "top-ranking-only",
                "size": int(settings.get("default_panel_size", 3)),
                "unique_genes": False,
                "rationale": "Highest configured candidate rank only.",
            }
        ]
    members = select_strategy_panels(candidates, definitions, categories)
    guide_escape = pd.read_csv(output / OUTPUT_NAMES["guide_escape"])
    guide_virtual = _read_optional(output / OUTPUT_NAMES["guide_virtual"])
    multiplex, comparison = summarize_multiplex_strategies(
        members, context.profiles.editor, guide_escape, guide_virtual
    )
    members.to_csv(output / OUTPUT_NAMES["members"], index=False)
    multiplex.to_csv(output / OUTPUT_NAMES["multiplex"], index=False)
    comparison.to_csv(output / OUTPUT_NAMES["comparison"], index=False)
    return {
        "output_dir": output,
        "strategy_count": len(multiplex),
        "panel_member_count": len(members),
    }


def _validate_assertions(context: ProjectContext) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for assertion in _analysis(context).get("source_assertions", []):
        path = _resolve(context, str(assertion["path"]))
        if not path or not path.is_file():
            raise FileNotFoundError(f"Assertion source is missing: {path}")
        if path.suffix == ".json":
            payload = json.loads(path.read_text(encoding="utf-8"))
            observed = payload[str(assertion["json_key"])]
        else:
            frame = pd.read_csv(path, low_memory=False)
            metric = str(assertion["metric"])
            if metric == "rows":
                observed = len(frame)
            elif metric == "unique":
                observed = frame[str(assertion["column"])].nunique()
            elif metric == "zero_value_rows":
                observed = int(
                    pd.to_numeric(frame[str(assertion["column"])], errors="coerce").eq(0).sum()
                )
            else:
                raise ValueError(f"Unsupported source assertion metric: {metric}")
        expected = assertion["expected"]
        if observed != expected:
            raise ValueError(
                f"Source assertion failed for {path}: observed {observed}, expected {expected}"
            )
        rows.append(
            {
                "label": assertion.get("label", path.name),
                "source": _display_path(context, path),
                "observed": observed,
                "expected": expected,
                "status": "pass",
            }
        )
    return rows


def _write_figures(output: Path) -> list[Path]:
    figure_dir = output / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    figures: list[Path] = []
    gene = pd.read_csv(output / OUTPUT_NAMES["gene_virtual"])
    if not gene.empty:
        shown = gene.sort_values("mapped_guide_count", ascending=False).head(15)
        fig, axis = plt.subplots(figsize=(10, 5))
        axis.bar(shown["gene_name"], shown["mapped_guide_count"], color="#3b82f6")
        axis.set_ylabel("Mapped guides in analyzed panel")
        axis.set_title("Virtual-knockout annotation coverage by gene")
        axis.tick_params(axis="x", rotation=55)
        fig.tight_layout()
        path = figure_dir / "gene_annotation_coverage.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        figures.append(path)
    comparison = pd.read_csv(output / OUTPUT_NAMES["comparison"])
    if not comparison.empty:
        fig, axis = plt.subplots(figsize=(9, 5))
        axis.bar(
            comparison["strategy"],
            comparison["sequence_escape_barrier"],
            color="#0f766e",
        )
        axis.set_ylabel("Minimum target-disrupting substitutions")
        axis.set_title("Exact-target sequence escape barrier by configured strategy")
        axis.tick_params(axis="x", rotation=25)
        fig.tight_layout()
        path = figure_dir / "multiplex_sequence_escape_barrier.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        figures.append(path)
    guide = pd.read_csv(output / OUTPUT_NAMES["guide_escape"])
    heldout = pd.to_numeric(guide["heldout_exact_target_coverage"], errors="coerce")
    discovery = pd.to_numeric(guide["discovery_exact_target_coverage"], errors="coerce")
    if discovery.notna().any():
        fig, axis = plt.subplots(figsize=(6, 6))
        axis.scatter(discovery, heldout, alpha=0.6, s=22, color="#7c3aed")
        axis.set_xlabel("Discovery exact-target coverage")
        axis.set_ylabel("Held-out exact-target coverage")
        axis.set_title("Observed target support across viral populations")
        fig.tight_layout()
        path = figure_dir / "discovery_vs_heldout_coverage.png"
        fig.savefig(path, dpi=160)
        plt.close(fig)
        figures.append(path)
    return figures


def _write_findings(
    context: ProjectContext, output: Path, validation: list[dict[str, Any]]
) -> None:
    guide = pd.read_csv(output / OUTPUT_NAMES["guide_virtual"])
    escape = pd.read_csv(output / OUTPUT_NAMES["guide_escape"])
    multiplex = pd.read_csv(output / OUTPUT_NAMES["multiplex"])
    focus = [str(value) for value in _analysis(context).get("focus_genes", [])]
    focus_rows = guide[guide["gene_name"].astype(str).isin(focus)] if focus else guide
    notes = list(_analysis(context).get("case_study_notes") or [])
    heldout_count = (
        pd.to_numeric(escape["heldout_exact_target_coverage"], errors="coerce").notna().sum()
    )
    lines = [
        "# Virtual knockout and escape findings",
        "",
        "## Computational observations",
        "",
        (
            "- The analyzed candidate table contained "
            f"{escape['candidate_id'].nunique():,} unique guides."
        ),
        (
            f"- {focus_rows['candidate_id'].nunique():,} guide-to-CDS mappings were "
            "available for the configured focus set."
        ),
        (
            f"- {len(multiplex):,} configured multiplex strategies were compared on "
            "separate sequence, host-risk, disruption, and evidence axes."
        ),
        (
            "- Held-out exact-target coverage was available for "
            f"{heldout_count:,} of {escape['candidate_id'].nunique():,} analyzed guides; "
            "missing values remain unknown."
        ),
    ]
    lines.extend(f"- {note}" for note in notes)
    lines.extend(
        [
            "",
            "## Research hypotheses",
            "",
            (
                "- Guides combining observed population support, interpretable coding "
                "disruption hypotheses, and independent multiplex targets may be useful "
                "candidates for independent experimental review."
            ),
            (
                "- Differences between single-guide rank, gene-level targetability, and "
                "sequence escape barrier should be evaluated rather than collapsed into "
                "one therapeutic score."
            ),
            "",
            "## Evidence gaps",
            "",
            (
                "- Sequence hypotheses do not establish that an edit occurs or that a gene "
                "disruption changes viral phenotype."
            ),
            (
                "- Biological evidence remains source- and virus-specific and requires "
                "explicit human review before use."
            ),
            (
                "- Marginal per-guide population coverage does not establish joint panel "
                "coverage in the same genomes."
            ),
            "",
            "## Limitations",
            "",
            (
                "- Indel sizes are an equally weighted bounded grid, not biological repair "
                "probabilities."
            ),
            (
                "- The escape barrier is a minimum exact-target sequence-change count, not "
                "an evolutionary probability."
            ),
            "- Predicted host matches are model-bounded and do not establish safety.",
            (
                "- No editing, viral inactivation, treatment, efficacy, delivery, or cure "
                "claim is made."
            ),
            "",
            "## Source validation",
            "",
        ]
    )
    lines.extend(
        f"- {row['label']}: {row['observed']} (expected {row['expected']}) — {row['status']}"
        for row in validation
    )
    (output / OUTPUT_NAMES["findings"]).write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_report(context: ProjectContext, output: Path) -> None:
    comparison = pd.read_csv(output / OUTPUT_NAMES["comparison"])
    gene = pd.read_csv(output / OUTPUT_NAMES["gene_virtual"])
    links = "".join(
        f"<li><a href='{html.escape(name)}'>{html.escape(name)}</a></li>"
        for name in OUTPUT_NAMES.values()
        if (output / name).is_file() and not name.endswith(".html")
    )
    text = f"""<!doctype html><html lang='en'><head><meta charset='utf-8'>
<title>Virtual knockout and escape analysis</title><style>
body{{font-family:Arial,sans-serif;max-width:1100px;margin:2rem auto;line-height:1.5}}
table{{border-collapse:collapse;width:100%;font-size:.88rem}}
th,td{{border:1px solid #ddd;padding:.4rem}}
.warning{{background:#fff7ed;border-left:5px solid #ea580c;padding:1rem}}
</style></head><body><h1>{html.escape(str(context.values.get("display_name")))}</h1>
<h2>Virtual knockout and exact-target escape analysis</h2>
<div class='warning'><strong>Computational scope only.</strong> Indels are size-defined sequence
hypotheses, and escape barriers are exact-target sequence-change counts. No repair frequency,
editing, safety, efficacy, viral inactivation, treatment, delivery, or cure is inferred.</div>
<h2>Strategy comparison</h2>{comparison.to_html(index=False, na_rep="unknown")}
<h2>Gene-level descriptive summaries</h2>{gene.to_html(index=False, na_rep="unknown")}
<h2>Auditable outputs</h2><ul>{links}</ul>
<h2>Interpretation boundary</h2><p>Host-risk, biological evidence, predicted disruption, and
escape robustness remain separate axes. See FINDINGS.md and run_manifest.json.</p>
</body></html>"""
    (output / OUTPUT_NAMES["report"]).write_text(text, encoding="utf-8")


def run_full_virtual_analysis(project: str | Path | ProjectContext) -> dict[str, Any]:
    """Run all stages and write a publication-facing, auditable bundle."""
    context = load_project(project) if not isinstance(project, ProjectContext) else project
    virtual = run_virtual_knockout_analysis(context)
    escape = run_escape_analysis(context)
    multiplex = run_multiplex_analysis(context)
    output = analysis_output_dir(context)
    validation = _validate_assertions(context)
    figures = _write_figures(output)
    _write_findings(context, output, validation)
    _write_report(context, output)
    reference, gff, candidates = _required_inputs(context)
    input_paths = [reference, gff, candidates, context.source, *context.profiles.source_paths]
    for key in (
        "heldout_population_table",
        "gene_ranking_table",
    ):
        path = _resolve(context, _analysis(context).get(key))
        if path and path.is_file():
            input_paths.append(path)
    manifest = {
        "schema_version": ANALYSIS_VERSION,
        "analysis": "generic virtual knockout and exact-target escape robustness",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "project_id": context.values["id"],
        "parameters": {
            key: value
            for key, value in _analysis(context).items()
            if key not in {"case_study_notes"}
        },
        "inputs": [
            {"path": _display_path(context, path), "sha256": sha256_file(path)}
            for path in dict.fromkeys(input_paths)
            if path.is_file()
        ],
        "source_assertions": validation,
        "outputs": {
            "virtual_knockout": {
                key: _display_path(context, value) if isinstance(value, Path) else value
                for key, value in virtual.items()
            },
            "escape": {
                key: _display_path(context, value) if isinstance(value, Path) else value
                for key, value in escape.items()
            },
            "multiplex": {
                key: _display_path(context, value) if isinstance(value, Path) else value
                for key, value in multiplex.items()
            },
            "figures": [_display_path(context, path) for path in figures],
        },
        "interpretation": (
            "Computational sequence hypotheses only; no repair-frequency, editing, safety, "
            "efficacy, viral inactivation, treatment, delivery, or cure claim."
        ),
    }
    (output / OUTPUT_NAMES["manifest"]).write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return {
        "output_dir": str(output.resolve()),
        "report": str((output / OUTPUT_NAMES["report"]).resolve()),
        "virtual_knockout": virtual,
        "escape": escape,
        "multiplex": multiplex,
        "source_assertion_count": len(validation),
        "figure_count": len(figures),
    }
