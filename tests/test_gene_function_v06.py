from __future__ import annotations

from pathlib import Path

import pandas as pd

from viral_safe_target.gene_function import (
    TARGET_GENES,
    CdsRecord,
    build_domain_overlap,
    map_candidates_to_protein,
    score_genes,
    simulate_indels,
)
from viral_safe_target.gene_function_reporting import write_gene_function_report


def _cds() -> dict[str, CdsRecord]:
    return {
        "UL3": CdsRecord(
            "UL3", "ref", 1, 30, 1, "ATGAAACCCGGGTTTAAACCCGGGTTTTAA", "MKPGFKPGF", "p+", "plus"
        ),
        "UL18": CdsRecord(
            "UL18",
            "ref",
            101,
            130,
            -1,
            "ATGAAACCCGGGTTTAAACCCGGGTTTTAA",
            "MKPGFKPGF",
            "p-",
            "minus",
        ),
    }


def _domains() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "gene_name": "UL3",
                "interpro_accession": "IPR_TEST",
                "domain_name": "test domain",
                "protein_start_1based": 2,
                "protein_end_1based": 8,
                "source_url": "https://example.test/domain",
            },
            {
                "gene_name": "UL18",
                "interpro_accession": "IPR_MINUS",
                "domain_name": "minus domain",
                "protein_start_1based": 1,
                "protein_end_1based": 9,
                "source_url": "https://example.test/minus",
            },
        ]
    )


def _mapping() -> pd.DataFrame:
    selected = pd.DataFrame(
        [
            {
                "candidate_id": "plus",
                "mapped_gene_for_analysis": "UL3",
                "reference_start_1based": 4,
                "reference_end_1based": 23,
                "cut_position": 19,
                "within_gene_rank": 1,
                "post_human_score": 0.9,
            },
            {
                "candidate_id": "minus",
                "mapped_gene_for_analysis": "UL18",
                "reference_start_1based": 106,
                "reference_end_1based": 125,
                "cut_position": 120,
                "within_gene_rank": 1,
                "post_human_score": 0.8,
            },
        ]
    )
    alignment = {"ref": "A" * 130, "strain": "A" * 130}
    aligned_cds = {
        gene: {"ref": record.nucleotide_sequence, "strain": record.nucleotide_sequence}
        for gene, record in _cds().items()
    }
    evolution = pd.DataFrame(
        [
            {
                "gene_name": "UL3",
                "median_pairwise_dN_dS": 0.1,
                "hsv1_ortholog_protein_identity": 0.8,
            },
            {
                "gene_name": "UL18",
                "median_pairwise_dN_dS": 0.2,
                "hsv1_ortholog_protein_identity": 0.9,
            },
        ]
    )
    disorder = pd.DataFrame(columns=["gene_name", "protein_start_1based", "protein_end_1based"])
    return map_candidates_to_protein(
        selected, _cds(), alignment, aligned_cds, _domains(), disorder, evolution
    )


def test_plus_and_minus_cut_boundaries_map_to_coding_coordinates() -> None:
    mapped = _mapping().set_index("candidate_id")
    assert mapped.loc["plus", "cut_cds_offset_0based"] == 19
    assert mapped.loc["plus", "amino_acid_position_1based"] == 7
    assert mapped.loc["minus", "cut_cds_offset_0based"] == 10
    assert mapped.loc["minus", "amino_acid_position_1based"] == 4
    assert mapped.loc["minus", "cds_strand"] == "-"


def test_indel_grid_is_complete_and_frame_probability_is_size_based() -> None:
    outcomes = simulate_indels(_mapping().iloc[[0]], _cds(), _domains())
    assert outcomes["indel_size_bp"].tolist() == list(range(-10, 11))
    assert outcomes[outcomes["indel_size_bp"].ne(0)]["frameshift"].mean() == 0.7
    assert (
        outcomes.loc[outcomes["indel_size_bp"].gt(0), "inserted_sequence_status"]
        .eq("unspecified_N_placeholder")
        .all()
    )
    assert outcomes["limitations"].str.contains("not an efficacy prediction").all()


def test_domain_overlap_and_hsv1_evidence_do_not_fill_hsv2_score() -> None:
    mapped = _mapping()
    overlap = build_domain_overlap(mapped, _domains())
    assert set(overlap["relation_to_cut"]) == {"cut_inside_domain"}
    evidence = pd.DataFrame(
        [
            {
                "gene_name": "UL3",
                "virus_type": "HSV-1",
                "evidence_strength": "direct",
                "essentiality_score": 1.0,
                "essentiality_call": "supported_essential",
            },
            {
                "gene_name": "UL3",
                "virus_type": "HSV-1",
                "evidence_strength": "indirect",
                "essentiality_score": pd.NA,
                "essentiality_call": "unknown",
            },
            {
                "gene_name": "UL18",
                "virus_type": "HSV-2",
                "evidence_strength": "direct",
                "essentiality_score": pd.NA,
                "essentiality_call": "unknown",
            },
        ]
    )
    evolution = pd.DataFrame(
        [
            {"gene_name": gene, "median_pairwise_dN_dS": 0.1, "hsv1_ortholog_protein_identity": 0.8}
            for gene in TARGET_GENES
        ]
    )
    scores = score_genes(
        mapped,
        simulate_indels(mapped, _cds(), _domains()),
        evidence,
        _domains(),
        pd.DataFrame(columns=["gene_name"]),
        evolution,
    ).set_index("gene_name")
    assert pd.isna(scores.loc["UL3", "evidence_based_essentiality_score"])
    assert scores.loc["UL3", "hsv1_ortholog_essentiality_score"] == 1.0
    assert scores.loc["UL3", "hsv1_essentiality_status"] == "supported_essential"
    assert pd.isna(scores.loc["UL18", "evidence_based_essentiality_score"])


def test_report_labels_computational_scope_and_excludes_protocol(tmp_path: Path) -> None:
    mapped = _mapping()
    outcomes = simulate_indels(mapped, _cds(), _domains())
    output = tmp_path / "report.html"
    write_gene_function_report(
        output,
        gene_scores=pd.DataFrame([{"gene_name": "UL3", "sequence_targetability_score": 0.9}]),
        evidence=pd.DataFrame(columns=["virus_type"]),
        mapping=mapped,
        outcomes=outcomes,
        domains=_domains(),
        disorder=pd.DataFrame(),
        domain_overlap=build_domain_overlap(mapped, _domains()),
        evolution=pd.DataFrame([{"gene_name": "UL3"}]),
    )
    text = output.read_text(encoding="utf-8")
    assert "computational hypotheses" in text
    assert "No wet-lab protocol" in text
