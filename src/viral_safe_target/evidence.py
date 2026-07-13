"""Curated viral gene-evidence loading without implicit biological claims."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

EVIDENCE_COLUMNS = [
    "virus",
    "reference_accession",
    "gene_name",
    "product",
    "biological_process",
    "essentiality_status",
    "latency_relevance",
    "evidence_level",
    "evidence_type",
    "source_identifier",
    "source_url_or_doi",
    "curator",
    "curation_date",
    "notes",
]


def read_gene_evidence(path: str | Path | None) -> pd.DataFrame:
    """Read a curated table; an absent table means evidence remains unavailable."""
    if path is None:
        return pd.DataFrame(columns=EVIDENCE_COLUMNS)
    frame = pd.read_csv(path, dtype=str).fillna("")
    missing = [column for column in EVIDENCE_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"Gene-evidence table is missing columns: {', '.join(missing)}")
    allowed_status = {"supported", "suggested", "unknown", "conflicting", ""}
    for column in ("essentiality_status", "latency_relevance"):
        invalid = sorted(set(frame[column]) - allowed_status)
        if invalid:
            raise ValueError(f"Invalid {column} values: {', '.join(invalid)}")
    return frame[EVIDENCE_COLUMNS].copy()


def merge_gene_evidence(candidates: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    """Attach at most one deterministic curated record to each candidate."""
    output = candidates.copy()
    if evidence.empty:
        for column in EVIDENCE_COLUMNS:
            if column not in {"reference_accession", "gene_name", "product"}:
                output[f"evidence_{column}"] = pd.NA
        output["evidence_level"] = pd.NA
        return output
    curated = evidence.sort_values(
        ["reference_accession", "gene_name", "source_identifier"], kind="mergesort"
    ).drop_duplicates(["reference_accession", "gene_name"], keep="first")
    rename = {
        column: f"evidence_{column}"
        for column in EVIDENCE_COLUMNS
        if column not in {"reference_accession", "gene_name", "product"}
    }
    curated = curated.rename(columns=rename).rename(
        columns={"evidence_evidence_level": "evidence_level"}
    )
    keep = ["reference_accession", "gene_name", *rename.values()]
    keep = ["evidence_level" if item == "evidence_evidence_level" else item for item in keep]
    return output.merge(curated[keep], on=["reference_accession", "gene_name"], how="left")
