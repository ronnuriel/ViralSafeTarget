"""Figures and reports for a presentation-ready ViralSafeTarget case study."""

# Report prose and table declarations are intentionally readable as complete strings.
# ruff: noqa: E501

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd

plt.switch_backend("Agg")


CATEGORY_COLORS = {
    "DNA_replication": "#3569a8",
    "capsid_assembly": "#7b61a8",
    "structural_egress": "#2f8f83",
    "exploratory_accessory": "#d17a22",
}


def _table(frame: pd.DataFrame, columns: list[str], rows: int = 40) -> str:
    available = [column for column in columns if column in frame]
    if frame.empty or not available:
        return '<p class="unknown">No completed rows are available.</p>'
    return (
        frame[available]
        .head(rows)
        .to_html(index=False, border=0, classes="data", render_links=True)
    )


def _save_figure(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def create_showcase_figures(
    figures_dir: Path,
    *,
    provenance: dict[str, Any],
    genome_candidates: pd.DataFrame,
    gene_scores: pd.DataFrame,
    candidates: pd.DataFrame,
    deep_panel: pd.DataFrame,
    strategy_summary: pd.DataFrame,
    population_genes: pd.DataFrame,
) -> list[Path]:
    figures_dir.mkdir(parents=True, exist_ok=True)
    created: list[Path] = []

    stages = [
        "Initial",
        "Pre-human eligible",
        "Human-screen panel",
        "Zero predicted hits",
        "Deep panel",
    ]
    counts = [
        int(provenance["initial_candidate_count"]),
        int(provenance["eligible_candidate_count"]),
        int(provenance["screening_panel_candidate_count"]),
        int(genome_candidates["human_total_predicted_hits"].eq(0).sum()),
        len(deep_panel),
    ]
    fig, ax = plt.subplots(figsize=(10, 4.8))
    bars = ax.barh(stages[::-1], counts[::-1], color="#3569a8")
    ax.set_xlabel("Candidate coordinates")
    ax.set_title("HSV-2 candidate funnel")
    ax.spines[["top", "right"]].set_visible(False)
    for bar, value in zip(bars, counts[::-1], strict=True):
        ax.text(bar.get_width(), bar.get_y() + bar.get_height() / 2, f"  {value:,}", va="center")
    path = figures_dir / "pipeline_funnel.png"
    _save_figure(fig, path)
    created.append(path)

    if not population_genes.empty:
        focus_genes = ["UL3", "UL10", "UL18", "UL20", "UL36", "UL52", "UL53", "UL19", "UL30"]
        population = population_genes[population_genes["gene_name"].isin(focus_genes)].copy()
        population["fully_supported_fraction"] = (
            population["exact_in_all_observable_records_candidate_count"]
            / population["locus_evaluable_unique_candidate_count"]
        )
        population = population.sort_values("fully_supported_fraction", ascending=True)
        fig, ax = plt.subplots(figsize=(9, 5.5))
        ax.barh(population["gene_name"], population["fully_supported_fraction"], color="#2f8f83")
        ax.set_xlim(0, 1)
        ax.set_xlabel("Fraction of evaluable guides exact in every observable held-out record")
        ax.set_title("Held-out population support remains separate from targetability")
        ax.grid(axis="x", alpha=0.2)
        ax.spines[["top", "right"]].set_visible(False)
        path = figures_dir / "heldout_population_support.png"
        _save_figure(fig, path)
        created.append(path)

    landscape = gene_scores.merge(
        candidates[["gene_name", "primary_category"]].drop_duplicates(), on="gene_name"
    )
    fig, ax = plt.subplots(figsize=(8.5, 6))
    for category, group in landscape.groupby("primary_category"):
        ax.scatter(
            group["sequence_targetability_score"],
            group["predicted_protein_disruption_score"],
            s=250 * group["evidence_coverage_score"],
            color=CATEGORY_COLORS.get(category, "#777777"),
            label=category.replace("_", " "),
            alpha=0.85,
            edgecolor="white",
        )
        for _, row in group.iterrows():
            ax.annotate(
                row["gene_name"],
                (row["sequence_targetability_score"], row["predicted_protein_disruption_score"]),
                xytext=(5, 5),
                textcoords="offset points",
                fontsize=9,
            )
    ax.set_xlabel("Sequence targetability")
    ax.set_ylabel("Predicted protein disruption")
    ax.set_title("Targetability and biological context are separate axes")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.2)
    ax.spines[["top", "right"]].set_visible(False)
    path = figures_dir / "gene_landscape.png"
    _save_figure(fig, path)
    created.append(path)

    heatmap = deep_panel.sort_values(
        ["primary_category", "gene_name", "deep_panel_rank_within_gene"]
    )
    metrics = [
        "targetability_percentile",
        "disruption_percentile",
        "evidence_coverage_percentile",
        "exact_strain_coverage",
    ]
    matrix = heatmap[metrics].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy()
    labels = [
        f"{row.gene_name} · {int(row.deep_panel_rank_within_gene)}" for row in heatmap.itertuples()
    ]
    fig, ax = plt.subplots(figsize=(8.5, max(7, len(labels) * 0.22)))
    image = ax.imshow(matrix, aspect="auto", cmap="Blues", vmin=0, vmax=1)
    ax.set_yticks(range(len(labels)), labels=labels, fontsize=7)
    ax.set_xticks(
        range(len(metrics)),
        labels=[item.replace("_", " ") for item in metrics],
        rotation=25,
        ha="right",
    )
    ax.set_title("Balanced deep panel: normalized computational evidence")
    fig.colorbar(image, ax=ax, fraction=0.03, pad=0.03)
    path = figures_dir / "deep_panel_heatmap.png"
    _save_figure(fig, path)
    created.append(path)

    plot = strategy_summary.set_index("strategy")[
        ["mean_sequence_targetability", "mean_predicted_disruption", "mean_evidence_coverage"]
    ]
    fig, ax = plt.subplots(figsize=(9.5, 5))
    plot.plot(kind="bar", ax=ax, color=["#3569a8", "#7b61a8", "#2f8f83"])
    ax.set_ylim(0, 1)
    ax.set_ylabel("Mean normalized or 0–1 score")
    ax.set_xlabel("")
    ax.set_title("Computational comparison sets emphasize different objectives")
    ax.legend(["Targetability", "Predicted disruption", "Evidence coverage"], frameon=False)
    ax.tick_params(axis="x", rotation=15)
    ax.spines[["top", "right"]].set_visible(False)
    path = figures_dir / "strategy_comparison.png"
    _save_figure(fig, path)
    created.append(path)
    return created


def write_showcase_documents(
    output: Path,
    *,
    profile_summary: dict[str, str],
    provenance: dict[str, Any],
    profile_checks: pd.DataFrame,
    genome_candidates: pd.DataFrame,
    gene_scores: pd.DataFrame,
    candidates: pd.DataFrame,
    deep_panel: pd.DataFrame,
    strategy_members: pd.DataFrame,
    strategy_summary: pd.DataFrame,
    external_validation: pd.DataFrame,
    population_candidates: pd.DataFrame,
    population_genes: pd.DataFrame,
    research_findings: pd.DataFrame,
) -> None:
    output.mkdir(parents=True, exist_ok=True)
    full_zero_hits = int(genome_candidates["human_total_predicted_hits"].eq(0).sum())
    mapped_zero_hits = int(candidates["human_total_predicted_hits"].eq(0).sum())
    direct_hsv2_essential_genes = int(
        candidates.loc[
            candidates["evidence_based_essentiality_score"].notna(), "gene_name"
        ].nunique()
    )
    population_evaluable = (
        int(population_candidates["locus_observable_record_count"].fillna(0).gt(0).sum())
        if not population_candidates.empty
        else 0
    )
    population_fully_supported = (
        int(
            population_candidates["population_validation_status"]
            .eq("exact_in_all_observable_records")
            .sum()
        )
        if not population_candidates.empty
        else 0
    )
    findings = (
        f"""# Findings

## Research question

Can a virus-first pipeline distinguish sequence targetability from biological target evidence and produce an auditable, mechanistically balanced shortlist?

## Data funnel

- {int(provenance["initial_candidate_count"]):,} initial candidate coordinates.
- {int(provenance["eligible_candidate_count"]):,} passed pre-human filters.
- {int(provenance["screening_panel_candidate_count"]):,} entered the balanced host screen.
- {full_zero_hits:,} candidates have zero predicted human hits under the declared model.
- {len(candidates):,} top candidates were mapped into the protein-disruption analysis; {mapped_zero_hits:,} of those have zero predicted human hits.
- {len(deep_panel):,} candidates form the balanced deep panel.

## Main finding

The highest sequence-targetability genes are not automatically the best-supported biological targets. The project therefore keeps targetability, direct essentiality evidence, predicted protein disruption, and evidence coverage as separate axes. No combined therapeutic score is reported.

The current curated set contains {direct_hsv2_essential_genes:,} genes with a direct HSV-2 essentiality score. Direct HSV-2 knockdown phenotypes remain phenotype evidence rather than null-essentiality claims. HSV-1 ortholog evidence is displayed separately.

## Potentially novel computational observations

The analysis produced {len(research_findings):,} auditable observations, hypotheses, robustness results, or evidence gaps in `research_findings.csv`. These rows are intended to help researchers decide what deserves independent investigation. They are not novelty claims against the complete literature and are not evidence of editing or treatment efficacy.

"""
        + "\n".join(
            f"- **{row.finding_id} ({row.classification}):** {row.observation} "
            f"Potential value: {row.potential_research_value} Limitation: {row.key_limitation}"
            for row in research_findings.itertuples(index=False)
        )
        + f"""

## External population-genomics context

The convenience strain alignment is checked against {len(external_validation):,} source-linked population or gene-variability findings. These studies support generally low HSV-2 divergence but identify lineage- and locus-specific exceptions, including UL30 and UL53. External studies contextualize the analysis; they do not replace held-out sequence validation.

In a separate held-out panel, {population_evaluable:,} candidates had an observable reference locus and {population_fully_supported:,} retained the exact guide/PAM in every observable record. Discovery genomes were excluded. Partial records remain unresolved where the locus is not aligned, and these results do not alter targetability scores.

## Presentation claim

ViralSafeTarget demonstrates a reproducible framework for discovering and comparing viral target hypotheses. It does not demonstrate editing, viral eradication, delivery, safety, efficacy, or a cure.
"""
    )
    methods = f"""# Methods

- Virus profile: `{profile_summary["virus_id"]}` ({profile_summary["virus_name"]}; {profile_summary["reference_accession"]}).
- Host profile: `{profile_summary["host_id"]}` ({profile_summary["host_assembly"]}).
- Nuclease profile: `{profile_summary["nuclease_id"]}` ({profile_summary["nuclease_name"]}).
- Host search: {provenance.get("cas_offinder_version", "Cas-OFFinder version unavailable")} through {profile_summary["mismatch_threshold"]} mismatches.
- Candidate shortlist: non-dominated fronts over targetability, predicted protein disruption, and evidence coverage, followed by a fixed per-gene quota.
- Biological evidence: source-linked rows only; HSV-1 and HSV-2 scopes are never merged.
- External validation: population-genomics findings are reported separately from essentiality evidence and computational scores.
- Held-out validation: discovery accessions are excluded; exact guide/PAM retention is evaluated only where a high-quality reference mapping covers the locus.
- Repair outcomes: deterministic size-defined sequence hypotheses, not repair-frequency predictions.
- Profile validation and input checksums are recorded in `run_manifest.json`.
"""
    limitations = """# Limitations

- No wet-lab validation, editing-efficiency measurement, delivery model, animal result, or clinical evidence is included.
- Zero predicted host hits is bounded by the assembly, editor, mismatch threshold, and search model; it is not proof of safety.
- The HSV-2 strain set is a convenience panel and is not a representative clinical population sample.
- Most held-out public HSV-2 records are partial; locus-specific denominators reduce but do not remove sampling and assembly bias.
- Whole-genome alignment, repeats, genome isomers, and annotation uncertainty can affect conservation and mapping.
- HSV-1 ortholog evidence cannot establish HSV-2 essentiality.
- Size-only indel models cannot predict repair frequencies or unspecified inserted bases.
- Comparison sets are computational scenarios, not experimental or treatment recommendations.
- Latent-neuron delivery and access to viral episomes remain outside this software's evidence.
"""
    (output / "FINDINGS.md").write_text(findings, encoding="utf-8")
    (output / "METHODS.md").write_text(methods, encoding="utf-8")
    (output / "LIMITATIONS.md").write_text(limitations, encoding="utf-8")

    sections = [
        ("Profile validation", _table(profile_checks, ["component", "status", "detail"], 30)),
        (
            "Potentially novel computational findings",
            "<p>These are auditable observations and research hypotheses generated from the completed balanced analysis. They may help prioritize independent work, but they are not claims of literature novelty, biological validation, safety, efficacy, or cure.</p>"
            + _table(
                research_findings,
                [
                    "finding_id",
                    "classification",
                    "observation",
                    "potential_research_value",
                    "computational_support",
                    "key_limitation",
                ],
                rows=20,
            ),
        ),
        (
            "Gene landscape",
            '<img src="figures/gene_landscape.png" alt="Scatter plot separating sequence targetability from predicted protein disruption">'
            + _table(
                gene_scores,
                [
                    "gene_name",
                    "sequence_targetability_score",
                    "evidence_based_essentiality_score",
                    "hsv1_ortholog_essentiality_score",
                    "predicted_protein_disruption_score",
                    "evidence_coverage_score",
                ],
            ),
        ),
        (
            "Discovery funnel",
            '<img src="figures/pipeline_funnel.png" alt="HSV-2 candidate funnel">',
        ),
        (
            "Balanced deep panel",
            '<img src="figures/deep_panel_heatmap.png" alt="Heatmap of the balanced deep panel">'
            + _table(
                deep_panel,
                [
                    "primary_category",
                    "gene_name",
                    "candidate_id",
                    "post_human_rank",
                    "pareto_front",
                    "candidate_predicted_disruption_score",
                    "evidence_tier",
                    "deep_panel_selection_reason",
                ],
                50,
            ),
        ),
        (
            "Computational comparison sets",
            '<img src="figures/strategy_comparison.png" alt="Comparison of computational strategy sets">'
            + _table(
                strategy_summary,
                [
                    "strategy",
                    "candidate_count",
                    "unique_gene_count",
                    "mechanistic_category_count",
                    "mean_sequence_targetability",
                    "mean_predicted_disruption",
                    "mean_evidence_coverage",
                    "interpretation",
                ],
            )
            + _table(
                strategy_members,
                [
                    "strategy",
                    "strategy_member_order",
                    "gene_name",
                    "candidate_id",
                    "primary_category",
                    "evidence_tier",
                    "strategy_rationale",
                ],
                30,
            ),
        ),
        (
            "Held-out population validation",
            '<img src="figures/heldout_population_support.png" alt="Held-out exact target support by focus gene">'
            + _table(
                population_genes[
                    population_genes["gene_name"].isin(
                        ["UL3", "UL10", "UL18", "UL20", "UL36", "UL52", "UL53", "UL19", "UL30"]
                    )
                ],
                [
                    "gene_name",
                    "locus_evaluable_unique_candidate_count",
                    "exact_in_all_observable_records_candidate_count",
                    "median_observable_locus_exact_target_coverage",
                    "minimum_observable_locus_exact_target_coverage",
                    "best_population_supported_candidate_id",
                ],
            )
            if not population_genes.empty
            else '<p class="unknown">Held-out population validation is unavailable.</p>',
        ),
        (
            "External population-genomics validation",
            "<p>Independent studies provide context for strain-panel generalizability. "
            "They are not used as essentiality evidence or combined into a therapeutic score.</p>"
            + _table(
                external_validation,
                [
                    "scope",
                    "gene_name",
                    "sample_context",
                    "sample_count",
                    "finding",
                    "implication_for_project",
                    "source_identifier",
                    "source_url",
                ],
                rows=20,
            ),
        ),
        (
            "Interpretation boundaries",
            "<p>Targetability is not essentiality. Ortholog evidence is not direct evidence. Predicted disruption is not viral inhibition. A comparison set is not a protocol or treatment recommendation.</p>",
        ),
    ]
    body = "".join(
        f"<section><h2>{html.escape(title)}</h2>{content}</section>" for title, content in sections
    )
    document = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>ViralSafeTarget HSV-2 showcase</title>
<style>body{{font:16px/1.5 system-ui,sans-serif;max-width:1400px;margin:auto;padding:2rem;color:#15202b}}h1,h2{{color:#183b56}}section{{margin:2.4rem 0}}img{{max-width:100%;height:auto;border:1px solid #d9e2ec}}table{{border-collapse:collapse;width:100%;font-size:.82rem;display:block;overflow:auto}}th,td{{padding:.45rem;border:1px solid #d9e2ec;text-align:left;vertical-align:top;white-space:normal;min-width:9rem;max-width:34rem}}th{{background:#eef4f8}}.warning{{background:#fff4e5;padding:1rem;border-left:5px solid #d97706}}.unknown{{color:#8a4b08}}</style></head><body>
<h1>ViralSafeTarget: HSV-2 evidence-aware case study</h1><div class="warning"><strong>Computational research showcase.</strong> No wet-lab protocol, safety conclusion, efficacy claim, or claim of cure is provided.</div>{body}
<section><h2>Reproducibility</h2><pre>{html.escape(json.dumps(profile_summary, indent=2, sort_keys=True))}</pre><p>See FINDINGS.md, METHODS.md, LIMITATIONS.md, and run_manifest.json in this directory.</p></section></body></html>"""
    (output / "FINAL_REPORT.html").write_text(document, encoding="utf-8")
