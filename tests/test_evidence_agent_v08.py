from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from viral_safe_target.evidence import read_gene_evidence
from viral_safe_target.evidence_agent import (
    PROPOSAL_COLUMNS,
    apply_reviewed_evidence,
    build_gene_catalog,
    build_search_queries,
    extract_proposals,
)


def _profile() -> dict[str, object]:
    return {
        "id": "test-virus",
        "display_name": "Test virus 2",
        "scientific_name": "Test alphavirus 2",
        "tax_id": 1234,
        "reference_accession": "TEST_REF.1",
        "literature_search_names": ["Test alphavirus 2", "TV-2"],
        "ortholog_search_names": ["Test alphavirus 1", "TV-1"],
    }


def _gff(path: Path) -> Path:
    path.write_text(
        "##gff-version 3\n"
        "TEST_REF.1\tRefSeq\tgene\t1\t300\t.\t+\t.\t"
        "ID=geneA;locus_tag=LOC_A;Name=UL52;gene=UL52\n"
        "TEST_REF.1\tRefSeq\tCDS\t1\t297\t.\t+\t0\t"
        "ID=cdsA;locus_tag=LOC_A;Name=UL52;gene=UL52;"
        "protein_id=PROT_A;product=helicase-primase%20subunit\n",
        encoding="utf-8",
    )
    return path


def test_gene_catalog_and_queries_are_alias_aware_and_scope_separated(tmp_path: Path) -> None:
    catalog = build_gene_catalog(_gff(tmp_path / "virus.gff3"), _profile())
    assert catalog["gene_name"].tolist() == ["UL52"]
    aliases = set(catalog.iloc[0]["aliases"].split(";"))
    assert {"UL52", "LOC_A", "PROT_A"}.issubset(aliases)
    assert catalog.iloc[0]["product"] == "helicase-primase subunit"

    queries = build_search_queries(catalog, _profile())
    assert len(queries) == 8
    assert set(queries["evidence_scope"]) == {"direct_target_virus", "ortholog"}
    direct = queries[queries["evidence_scope"].eq("direct_target_virus")]
    assert direct["query_text"].str.contains("Test alphavirus 2", regex=False).all()
    assert direct["query_text"].str.contains("UL52", regex=False).all()


def test_expression_is_not_promoted_to_essentiality_and_all_proposals_are_pending(
    tmp_path: Path,
) -> None:
    catalog = build_gene_catalog(_gff(tmp_path / "virus.gff3"), _profile())
    queries = build_search_queries(catalog, _profile())
    query_id = queries.iloc[0]["query_id"]
    records = pd.DataFrame(
        [
            {
                "source_database": "PubMed",
                "source_identifier": "PMID:1",
                "source_url": "https://pubmed.ncbi.nlm.nih.gov/1/",
                "title": "UL52 expression in Test alphavirus 2 infection",
                "abstract_or_annotation": (
                    "UL52 expression increased during Test alphavirus 2 infection "
                    "in cultured cells."
                ),
                "publication_year": "2024",
                "doi": "",
                "query_ids": query_id,
                "source_record_sha256": "abc",
            }
        ]
    )
    proposals = extract_proposals(records, catalog, queries, _profile())
    assert len(proposals) == 1
    proposal = proposals.iloc[0]
    assert proposal["experiment_type"] == "expression"
    assert proposal["proposed_essentiality_call"] == "unknown"
    assert proposal["evidence_direction"] == "association_or_expression_only"
    assert proposal["review_status"] == "pending"
    assert len(proposal["quoted_evidence_span"].split()) <= 26


def test_null_mutant_interpretation_is_still_only_a_pending_proposal(tmp_path: Path) -> None:
    catalog = build_gene_catalog(_gff(tmp_path / "virus.gff3"), _profile())
    queries = build_search_queries(catalog, _profile())
    query_id = queries.iloc[0]["query_id"]
    records = pd.DataFrame(
        [
            {
                "source_database": "PubMed",
                "source_identifier": "PMID:2",
                "source_url": "https://pubmed.ncbi.nlm.nih.gov/2/",
                "title": "UL52 null mutant in Test alphavirus 2",
                "abstract_or_annotation": (
                    "A UL52 null mutant was unable to replicate in cultured cells."
                ),
                "publication_year": "2025",
                "doi": "",
                "query_ids": query_id,
                "source_record_sha256": "def",
            }
        ]
    )
    proposal = extract_proposals(records, catalog, queries, _profile()).iloc[0]
    assert proposal["experiment_type"] == "null_mutant"
    assert proposal["proposed_essentiality_call"] == "supported_required_in_reported_context"
    assert proposal["proposed_evidence_strength"] == "direct"
    assert proposal["review_status"] == "pending"


def test_another_virus_named_in_title_is_not_labeled_direct_evidence(tmp_path: Path) -> None:
    catalog = build_gene_catalog(_gff(tmp_path / "virus.gff3"), _profile())
    queries = build_search_queries(catalog, _profile())
    records = pd.DataFrame(
        [
            {
                "source_database": "PubMed",
                "source_identifier": "PMID:other",
                "source_url": "https://pubmed.ncbi.nlm.nih.gov/4/",
                "title": "Bovine herpesvirus UL52 null mutant",
                "abstract_or_annotation": (
                    "The discussion compares this result with Test alphavirus 2."
                ),
                "publication_year": "2025",
                "doi": "",
                "query_ids": queries.iloc[0]["query_id"],
                "source_record_sha256": "ghi",
            }
        ]
    )
    proposal = extract_proposals(records, catalog, queries, _profile()).iloc[0]
    assert proposal["evidence_scope"] == "mixed_or_other_virus"
    assert proposal["proposed_evidence_strength"] == "unknown"


def _review_queue(tmp_path: Path) -> Path:
    values = {column: "" for column in PROPOSAL_COLUMNS}
    values.update(
        {
            "proposal_id": "EP-approved",
            "gene_name": "UL52",
            "target_virus": "Test alphavirus 2",
            "evidence_virus": "Test alphavirus 1",
            "evidence_scope": "ortholog",
            "source_database": "PubMed",
            "source_identifier": "PMID:2",
            "source_url": "https://pubmed.ncbi.nlm.nih.gov/2/",
            "source_title": "A reviewed source",
            "evidence_category": "essentiality_or_replication",
            "experiment_type": "null_mutant",
            "model_system": "cultured_cells",
            "finding_summary": "A reviewed, context-bounded finding.",
            "proposed_essentiality_call": "supported_required_in_reported_context",
            "proposed_essentiality_score": "1.0",
            "proposed_evidence_strength": "indirect",
            "review_status": "approved",
            "reviewer": "Researcher A",
            "review_date": "2026-07-14",
            "review_notes": "Ortholog evidence only.",
        }
    )
    pending = dict(values)
    pending.update(
        {
            "proposal_id": "EP-pending",
            "source_identifier": "PMID:3",
            "review_status": "pending",
            "reviewer": "",
            "review_date": "",
        }
    )
    path = tmp_path / "review_queue.tsv"
    pd.DataFrame([values, pending], columns=PROPOSAL_COLUMNS).to_csv(path, sep="\t", index=False)
    return path


def test_apply_exports_only_explicitly_approved_rows_and_normalizes_for_scoring(
    tmp_path: Path,
) -> None:
    output = tmp_path / "gene_evidence.tsv"
    summary = apply_reviewed_evidence(
        _review_queue(tmp_path), output, reference_accession="TEST_REF.1"
    )
    assert summary["approved_count"] == 1
    assert summary["pending_count"] == 1
    curated = pd.read_csv(output, sep="\t")
    assert curated["source_identifier"].tolist() == ["PMID:2"]
    assert curated.iloc[0]["reference_accession"] == "TEST_REF.1"
    assert curated.iloc[0]["virus_type"] == "Test alphavirus 1"
    assert curated.iloc[0]["reviewer"] == "Researcher A"
    assert curated.iloc[0]["review_date"] == "2026-07-14"

    normalized = read_gene_evidence(output)
    assert normalized.iloc[0]["evidence_level"] == "limited"
    assert normalized.iloc[0]["essentiality_status"] == "supported"


def test_approved_rows_require_reviewer_and_review_date(tmp_path: Path) -> None:
    queue = _review_queue(tmp_path)
    frame = pd.read_csv(queue, sep="\t", dtype=str).fillna("")
    frame.loc[frame["proposal_id"].eq("EP-approved"), "reviewer"] = ""
    frame.to_csv(queue, sep="\t", index=False)
    with pytest.raises(ValueError, match="require reviewer"):
        apply_reviewed_evidence(queue, tmp_path / "gene_evidence.tsv")
