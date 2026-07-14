"""Configuration-driven presentation workflow for standardized discovery outputs."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from .profiles import load_profile_bundle, validate_profile_bundle
from .provenance import sha256_file
from .showcase import (
    build_comparison_sets,
    build_evidence_aware_candidates,
    build_research_findings,
    select_balanced_deep_panel,
)
from .showcase_reporting import create_showcase_figures, write_showcase_documents


def run_showcase(
    *,
    virus_profile: str | Path,
    host_profile: str | Path,
    nuclease_profile: str | Path,
    out_dir: str | Path,
    project_root: str | Path | None = None,
    per_gene: int = 4,
) -> dict[str, Any]:
    root = Path(project_root).resolve() if project_root else Path.cwd().resolve()
    bundle = load_profile_bundle(virus_profile, host_profile, nuclease_profile, project_root=root)
    checks = validate_profile_bundle(bundle)
    failures = checks[checks["status"].eq("fail")]
    if not failures.empty:
        raise ValueError("Profile validation failed:\n" + failures.to_string(index=False))
    case_outputs = bundle.virus.get("case_study_outputs", {})
    genome_dir = bundle.resolve(case_outputs.get("genome_wide"))
    function_dir = bundle.resolve(case_outputs.get("gene_function"))
    category_path = bundle.resolve(bundle.virus.get("gene_category_table"))
    if genome_dir is None or function_dir is None or category_path is None:
        raise ValueError("The case-study profile must configure standardized output directories")
    paths = {
        "genome_candidates": genome_dir / "genome_wide_candidates_post_human.csv",
        "gene_rankings": genome_dir / "gene_rankings.csv",
        "gene_stability": genome_dir / "gene_rank_stability.csv",
        "provenance": genome_dir / "provenance.json",
        "mapping": function_dir / "candidate_protein_mapping.csv",
        "outcomes": function_dir / "predicted_disruption_outcomes.csv",
        "gene_scores": function_dir / "gene_scores.csv",
        "evolution": function_dir / "gene_evolution.csv",
        "evidence": function_dir / "gene_evidence.tsv",
        "categories": category_path,
    }
    external_validation_path = bundle.resolve(bundle.virus.get("external_validation_table"))
    if external_validation_path is not None:
        paths["external_validation"] = external_validation_path
    population_candidates_path = bundle.resolve(
        bundle.virus.get("population_validation_candidates")
    )
    population_genes_path = bundle.resolve(bundle.virus.get("population_validation_genes"))
    if population_candidates_path is not None:
        paths["population_candidates"] = population_candidates_path
    if population_genes_path is not None:
        paths["population_genes"] = population_genes_path
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Standardized showcase inputs are missing: " + ", ".join(missing))
    genome_candidates = pd.read_csv(paths["genome_candidates"])
    gene_rankings = pd.read_csv(paths["gene_rankings"])
    gene_stability = pd.read_csv(paths["gene_stability"])
    mapping = pd.read_csv(paths["mapping"])
    outcomes = pd.read_csv(paths["outcomes"])
    gene_scores = pd.read_csv(paths["gene_scores"])
    evolution = pd.read_csv(paths["evolution"])
    evidence = pd.read_csv(paths["evidence"], sep="\t")
    categories = pd.read_csv(paths["categories"], sep="\t")
    external_validation = (
        pd.read_csv(paths["external_validation"], sep="\t")
        if "external_validation" in paths
        else pd.DataFrame()
    )
    population_candidates = (
        pd.read_csv(paths["population_candidates"])
        if "population_candidates" in paths
        else pd.DataFrame()
    )
    population_genes = (
        pd.read_csv(paths["population_genes"]) if "population_genes" in paths else pd.DataFrame()
    )
    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))

    candidates = build_evidence_aware_candidates(
        mapping, outcomes, gene_scores, evidence, categories
    )
    deep_panel = select_balanced_deep_panel(candidates, per_gene=per_gene)
    strategy_members, strategy_summary = build_comparison_sets(candidates)
    if not population_candidates.empty:
        population_columns = [
            "candidate_id",
            "locus_observable_record_count",
            "exact_target_in_observable_locus_count",
            "observable_locus_exact_target_coverage",
            "locus_unresolved_record_count",
            "population_validation_status",
            "population_validation_used_in_targetability_score",
        ]
        available = [column for column in population_columns if column in population_candidates]
        candidates = candidates.merge(
            population_candidates[available], on="candidate_id", how="left", validate="one_to_one"
        )
        deep_panel = deep_panel.merge(
            population_candidates[available], on="candidate_id", how="left", validate="one_to_one"
        )
    research_findings = build_research_findings(
        gene_rankings,
        gene_stability,
        gene_scores,
        candidates,
        population_candidates,
        evolution,
    )
    output = Path(out_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    candidates.to_csv(output / "candidates_evidence_aware.csv", index=False)
    deep_panel.to_csv(output / "deep_screening_panel.csv", index=False)
    strategy_members.to_csv(output / "comparison_set_members.csv", index=False)
    strategy_summary.to_csv(output / "strategy_comparison.csv", index=False)
    research_findings.to_csv(output / "research_findings.csv", index=False)
    checks.to_csv(output / "profile_validation.csv", index=False)
    external_validation.to_csv(output / "external_validation.csv", index=False)
    population_candidates.to_csv(output / "population_validation_candidates.csv", index=False)
    population_genes.to_csv(output / "population_validation_genes.csv", index=False)
    figures = create_showcase_figures(
        output / "figures",
        provenance=provenance,
        genome_candidates=genome_candidates,
        gene_scores=gene_scores,
        candidates=candidates,
        deep_panel=deep_panel,
        strategy_summary=strategy_summary,
        population_genes=population_genes,
    )
    profile_summary = {
        "virus_id": str(bundle.virus["id"]),
        "virus_name": str(bundle.virus["display_name"]),
        "reference_accession": str(bundle.virus["reference_accession"]),
        "host_id": str(bundle.host["id"]),
        "host_assembly": str(bundle.host["assembly_name"]),
        "nuclease_id": str(bundle.nuclease["id"]),
        "nuclease_name": str(bundle.editor.name),
        "mismatch_threshold": str(bundle.editor.mismatch_search_threshold),
    }
    write_showcase_documents(
        output,
        profile_summary=profile_summary,
        provenance=provenance,
        profile_checks=checks,
        genome_candidates=genome_candidates,
        gene_scores=gene_scores,
        candidates=candidates,
        deep_panel=deep_panel,
        strategy_members=strategy_members,
        strategy_summary=strategy_summary,
        external_validation=external_validation,
        population_candidates=population_candidates,
        population_genes=population_genes,
        research_findings=research_findings,
    )
    manifest = {
        "schema_version": "1.0",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "profile_summary": profile_summary,
        "profile_sources": [
            {"path": str(path), "sha256": sha256_file(path)} for path in bundle.source_paths
        ],
        "analysis_inputs": [
            {"name": name, "path": str(path), "sha256": sha256_file(path)}
            for name, path in paths.items()
        ],
        "parameters": {"deep_panel_per_gene": per_gene},
        "outputs": {
            "candidate_count": len(candidates),
            "deep_panel_count": len(deep_panel),
            "comparison_set_member_count": len(strategy_members),
            "research_finding_count": len(research_findings),
            "figure_count": len(figures),
        },
        "interpretation": (
            "Computational research showcase only; no wet-lab protocol, safety, efficacy, "
            "delivery, or cure claim."
        ),
    }
    (output / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    return {
        "output_dir": output,
        "candidate_count": len(candidates),
        "deep_panel_count": len(deep_panel),
        "strategy_count": strategy_summary["strategy"].nunique(),
        "research_finding_count": len(research_findings),
        "figure_count": len(figures),
    }
