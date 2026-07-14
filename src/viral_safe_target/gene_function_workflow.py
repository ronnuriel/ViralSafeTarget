"""End-to-end orchestration for HSV-2 gene function and disruption analysis."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from itertools import combinations
from pathlib import Path
from typing import Any

import pandas as pd

from .gene_function import (
    TARGET_GENES,
    build_domain_overlap,
    compute_gene_evolution,
    map_candidates_to_protein,
    read_target_cds,
    score_genes,
    select_top_candidates,
    simulate_indels,
    simulate_paired_deletions,
)
from .gene_function_reporting import write_gene_function_report
from .io_utils import read_fasta
from .provenance import sha256_file


def _read_tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).replace({"": pd.NA})


def run_gene_function_analysis(
    *,
    genome_wide_dir: str | Path = "reports/hsv2_genome_wide",
    hsv2_genbank: str | Path = "data/raw/hsv2_reference.gb",
    hsv1_genbank: str | Path = "data/raw/hsv1_reference.gb",
    virus_alignment: str | Path = "data/processed/hsv2_aligned_25.fasta",
    domain_table: str | Path = "data/curated/hsv2_target_domains.tsv",
    disorder_table: str | Path = "data/curated/hsv2_target_disorder.tsv",
    evidence_table: str | Path = "data/curated/hsv_gene_function_evidence.tsv",
    out_dir: str | Path = "reports/hsv2_gene_function",
    top_per_gene: int = 10,
) -> dict[str, Any]:
    source = Path(genome_wide_dir)
    output = Path(out_dir)
    output.mkdir(parents=True, exist_ok=True)
    candidates_path = source / "genome_wide_candidates_post_human.csv"
    feature_map_path = source / "candidate_feature_map.csv"
    pairs_path = source / "pair_hypotheses_same_gene.csv"
    input_paths = [
        candidates_path,
        feature_map_path,
        pairs_path,
        Path(hsv2_genbank),
        Path(hsv1_genbank),
        Path(virus_alignment),
        Path(domain_table),
        Path(disorder_table),
        Path(evidence_table),
    ]
    missing = [str(path) for path in input_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required analysis inputs are missing: " + ", ".join(missing))

    candidates = pd.read_csv(candidates_path)
    feature_map = pd.read_csv(feature_map_path)
    legacy_pairs = pd.read_csv(pairs_path)
    if legacy_pairs.empty:
        raise ValueError("The v0.5 pair-hypothesis artifact is unexpectedly empty")
    hsv2_cds = read_target_cds(hsv2_genbank)
    hsv1_cds = read_target_cds(hsv1_genbank)
    selected = select_top_candidates(candidates, feature_map, hsv2_cds, top_per_gene=top_per_gene)
    aligned_records = read_fasta(virus_alignment)
    reference_id = next(
        (accession for accession in aligned_records if accession.startswith("NC_001798")),
        next(iter(aligned_records)),
    )
    domains = _read_tsv(Path(domain_table))
    disorder = _read_tsv(Path(disorder_table))
    evidence = _read_tsv(Path(evidence_table))
    for frame, numeric_columns in (
        (domains, ["protein_start_1based", "protein_end_1based"]),
        (disorder, ["protein_start_1based", "protein_end_1based"]),
        (evidence, ["essentiality_score"]),
    ):
        for column in numeric_columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")

    evolution, aligned_cds = compute_gene_evolution(
        aligned_records, reference_id, hsv2_cds, hsv1_cds
    )
    mapping = map_candidates_to_protein(
        selected,
        hsv2_cds,
        aligned_records,
        aligned_cds,
        domains,
        disorder,
        evolution,
    )
    single_outcomes = simulate_indels(mapping, hsv2_cds, domains)
    pair_rows: list[dict[str, str]] = []
    for gene, group in mapping.groupby("gene_name", sort=True):
        ids = group.sort_values("within_gene_rank")["candidate_id"].astype(str).tolist()
        pair_rows.extend(
            {
                "candidate_a": first,
                "candidate_b": second,
                "gene_a": gene,
                "gene_b": gene,
            }
            for first, second in combinations(ids, 2)
        )
    pairs = pd.DataFrame(pair_rows)
    pair_outcomes = simulate_paired_deletions(pairs, candidates, feature_map, hsv2_cds, domains)
    outcomes = pd.concat([single_outcomes, pair_outcomes], ignore_index=True)
    overlap = build_domain_overlap(mapping, domains)
    gene_scores = score_genes(mapping, outcomes, evidence, domains, disorder, evolution)
    evidence_output = evidence.merge(
        gene_scores, on="gene_name", how="left", validate="many_to_one"
    )

    evidence_output.to_csv(output / "gene_evidence.tsv", sep="\t", index=False)
    mapping.to_csv(output / "candidate_protein_mapping.csv", index=False)
    outcomes.to_csv(output / "predicted_disruption_outcomes.csv", index=False)
    overlap.to_csv(output / "domain_overlap.csv", index=False)
    gene_scores.to_csv(output / "gene_scores.csv", index=False)
    evolution.to_csv(output / "gene_evolution.csv", index=False)
    write_gene_function_report(
        output / "gene_evidence_report.html",
        gene_scores=gene_scores,
        evidence=evidence_output,
        mapping=mapping,
        outcomes=outcomes,
        domains=domains,
        disorder=disorder,
        domain_overlap=overlap,
        evolution=evolution,
    )
    manifest = {
        "analysis": "HSV-2 gene function and predicted disruption",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "target_genes": list(TARGET_GENES),
        "top_per_gene": top_per_gene,
        "reference_id": reference_id,
        "inputs": [
            {"path": str(path.resolve()), "sha256": sha256_file(path)} for path in input_paths
        ],
        "definitions": {
            "small_indels": (
                "-10 through +10 bp; downstream-anchored deletions; unspecified N insertions"
            ),
            "paired_deletions": "theoretical cut-to-cut sequence deletion only",
            "essentiality": "cited evidence only; HSV-1 and HSV-2 kept separate",
        },
        "limitations": [
            "computational hypotheses only",
            "no wet-lab protocol",
            "no safety or efficacy conclusion",
            "absence of evidence remains unknown",
        ],
    }
    (output / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    return {
        "output_dir": output.resolve(),
        "candidate_count": len(mapping),
        "single_outcome_count": len(single_outcomes),
        "pair_outcome_count": len(pair_outcomes),
        "gene_count": mapping["gene_name"].nunique(),
    }
