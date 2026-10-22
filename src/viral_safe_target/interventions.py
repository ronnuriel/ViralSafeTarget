"""Source-linked intervention evidence with mandatory human review.

This module does not recommend treatments. It validates user-supplied evidence,
keeps reviewed and pending records separate, and exposes missingness explicitly.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from urllib.parse import quote_plus

import pandas as pd

INTERVENTION_EVIDENCE_COLUMNS = [
    "record_id",
    "gene_name",
    "target_virus",
    "evidence_virus",
    "intervention_name",
    "intervention_type",
    "target_role",
    "development_status",
    "evidence_category",
    "evidence_scope",
    "source_database",
    "source_identifier",
    "source_url",
    "finding",
    "directness_notes",
    "review_status",
    "reviewer",
    "review_date",
    "review_notes",
    "retrieved_utc",
    "source_record_sha256",
]

ENUMS = {
    "intervention_type": {"small_molecule", "biologic", "genome_editor", "other"},
    "target_role": {
        "direct_target",
        "prodrug_activation",
        "resistance_marker",
        "complex_component",
        "pathway_association",
        "unknown",
    },
    "development_status": {
        "approved_for_target_virus",
        "approved_other_indication",
        "clinical_investigation",
        "preclinical",
        "tool_compound",
        "unknown",
    },
    "evidence_category": {
        "regulatory_label",
        "clinical_trial",
        "biochemical",
        "phenotypic",
        "resistance",
        "structural",
        "database_record",
        "other",
    },
    "evidence_scope": {"direct_target_virus", "ortholog", "related_virus", "unknown"},
    "review_status": {"pending", "approved", "rejected"},
}

REQUIRED_NONEMPTY = {
    "record_id",
    "gene_name",
    "target_virus",
    "evidence_virus",
    "intervention_name",
    "intervention_type",
    "target_role",
    "development_status",
    "evidence_category",
    "evidence_scope",
    "source_database",
    "source_identifier",
    "source_url",
    "finding",
    "review_status",
}

INTERVENTION_LANDSCAPE_COLUMNS = [
    "gene_name",
    "evidence_records",
    "approved_records",
    "pending_records",
    "rejected_records",
    "direct_target_records",
    "interventions",
    "development_statuses",
    "evidence_scopes",
    "source_count",
    "review_marker",
    "interpretation",
]


def empty_intervention_evidence() -> pd.DataFrame:
    """Return an empty table with the stable public schema."""
    return pd.DataFrame(columns=INTERVENTION_EVIDENCE_COLUMNS)


def _record_hash(row: pd.Series) -> str:
    values = {column: str(row.get(column, "")) for column in INTERVENTION_EVIDENCE_COLUMNS}
    values["source_record_sha256"] = ""
    payload = json.dumps(values, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def validate_intervention_evidence(table: pd.DataFrame) -> pd.DataFrame:
    """Validate and normalize an intervention table without approving any record."""
    missing = sorted(set(INTERVENTION_EVIDENCE_COLUMNS) - set(table.columns))
    if missing:
        raise ValueError("Intervention evidence is missing columns: " + ", ".join(missing))
    normalized = table[INTERVENTION_EVIDENCE_COLUMNS].fillna("").astype(str).copy()
    for column in REQUIRED_NONEMPTY:
        empty_rows = normalized.index[normalized[column].str.strip().eq("")].tolist()
        if empty_rows:
            raise ValueError(f"Intervention evidence has empty {column!r} at rows {empty_rows}")
    if normalized["record_id"].duplicated().any():
        duplicate = normalized.loc[normalized["record_id"].duplicated(), "record_id"].iloc[0]
        raise ValueError(f"Duplicate intervention record_id: {duplicate}")
    for column, allowed in ENUMS.items():
        invalid = sorted(set(normalized[column]) - allowed)
        if invalid:
            raise ValueError(f"Invalid {column} values: {invalid}")
    invalid_urls = normalized.loc[
        ~normalized["source_url"].str.startswith("https://"), "source_url"
    ]
    if not invalid_urls.empty:
        raise ValueError(f"Intervention source_url must use HTTPS: {invalid_urls.iloc[0]}")
    approved = normalized["review_status"].eq("approved")
    incomplete_review = approved & (
        normalized["reviewer"].str.strip().eq("")
        | normalized["review_date"].str.strip().eq("")
    )
    if incomplete_review.any():
        records = normalized.loc[incomplete_review, "record_id"].tolist()
        raise ValueError(
            "Approved intervention records require reviewer and review_date: "
            + ", ".join(records)
        )
    for index, row in normalized.iterrows():
        digest = row["source_record_sha256"].strip()
        invalid_digest = len(digest) != 64 or any(
            character not in "0123456789abcdef" for character in digest
        )
        if digest and invalid_digest:
            raise ValueError(f"Invalid source_record_sha256 for {row['record_id']}")
        if not digest:
            normalized.at[index, "source_record_sha256"] = _record_hash(row)
    return normalized


def read_intervention_evidence(path: str | Path) -> pd.DataFrame:
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"Intervention evidence table not found: {source}")
    return validate_intervention_evidence(pd.read_csv(source, sep="\t", dtype=str))


def build_intervention_landscape(
    evidence: pd.DataFrame, *, genes: list[str] | None = None
) -> pd.DataFrame:
    """Summarize intervention evidence per gene without producing a therapeutic score."""
    validated = validate_intervention_evidence(evidence)
    gene_names = sorted(set(genes or []) | set(validated["gene_name"]))
    rows: list[dict[str, object]] = []
    for gene in gene_names:
        group = validated[validated["gene_name"].eq(gene)]
        pending = group[group["review_status"].eq("pending")]
        approved = group[group["review_status"].eq("approved")]
        active = group[~group["review_status"].eq("rejected")]
        names = sorted(set(active["intervention_name"]))
        statuses = sorted(set(active["development_status"]))
        scopes = sorted(set(active["evidence_scope"]))
        rows.append(
            {
                "gene_name": gene,
                "evidence_records": len(group),
                "approved_records": len(approved),
                "pending_records": len(pending),
                "rejected_records": int(group["review_status"].eq("rejected").sum()),
                "direct_target_records": int(active["target_role"].eq("direct_target").sum()),
                "interventions": "; ".join(names) if names else "not located in supplied evidence",
                "development_statuses": "; ".join(statuses) if statuses else "unknown",
                "evidence_scopes": "; ".join(scopes) if scopes else "unknown",
                "source_count": int(active["source_identifier"].nunique()),
                "review_marker": "*" if len(pending) else "",
                "interpretation": (
                    "pending human review"
                    if len(pending)
                    else "human-reviewed evidence available"
                    if len(approved)
                    else "not located in supplied evidence"
                ),
            }
        )
    return pd.DataFrame(rows, columns=INTERVENTION_LANDSCAPE_COLUMNS)


def build_intervention_queries(genes: list[str], virus_names: list[str]) -> pd.DataFrame:
    """Create reproducible database searches; the URLs are leads, not evidence."""
    virus_term = " OR ".join(f'"{name}"' for name in virus_names if name) or '"virus"'
    rows: list[dict[str, str]] = []
    for gene in sorted(set(genes)):
        pubmed_query = (
            f"({virus_term}) AND ({gene}) AND (drug OR inhibitor OR resistance OR antiviral)"
        )
        clinical_query = f"{virus_names[0] if virus_names else 'virus'} {gene} antiviral"
        chembl_query = f"{virus_names[0] if virus_names else 'virus'} {gene}"
        for database, query, url in (
            (
                "PubMed",
                pubmed_query,
                "https://pubmed.ncbi.nlm.nih.gov/?term=" + quote_plus(pubmed_query),
            ),
            (
                "ChEMBL",
                chembl_query,
                "https://www.ebi.ac.uk/chembl/explore/target?query=" + quote_plus(chembl_query),
            ),
            (
                "ClinicalTrials.gov",
                clinical_query,
                "https://clinicaltrials.gov/search?term=" + quote_plus(clinical_query),
            ),
        ):
            rows.append({"gene_name": gene, "database": database, "query": query, "url": url})
    return pd.DataFrame(rows)


def write_intervention_outputs(
    evidence_path: str | Path | pd.DataFrame,
    out_dir: str | Path,
    *,
    genes: list[str] | None = None,
    virus_names: list[str] | None = None,
) -> dict[str, Path]:
    """Validate evidence and write researcher-facing tables and an HTML report."""
    output = Path(out_dir)
    output.mkdir(parents=True, exist_ok=True)
    evidence = (
        validate_intervention_evidence(evidence_path)
        if isinstance(evidence_path, pd.DataFrame)
        else read_intervention_evidence(evidence_path)
    )
    landscape = build_intervention_landscape(evidence, genes=genes)
    review_queue = evidence[evidence["review_status"].eq("pending")].copy()
    queries = build_intervention_queries(
        landscape["gene_name"].tolist(), virus_names or ["virus"]
    )
    evidence_output = output / "intervention_evidence_validated.tsv"
    landscape_output = output / "intervention_landscape.tsv"
    review_output = output / "intervention_review_queue.tsv"
    query_output = output / "intervention_queries.tsv"
    report_output = output / "intervention_landscape.html"
    evidence.to_csv(evidence_output, sep="\t", index=False)
    landscape.to_csv(landscape_output, sep="\t", index=False)
    review_queue.to_csv(review_output, sep="\t", index=False)
    queries.to_csv(query_output, sep="\t", index=False)
    report_output.write_text(
        "<!doctype html><html lang='en'><meta charset='utf-8'>"
        "<title>Intervention evidence landscape</title>"
        "<style>body{font:16px system-ui;max-width:1180px;margin:2rem auto;padding:0 1rem;"
        "line-height:1.5;color:#172033}table{border-collapse:collapse;width:100%;display:block;"
        "overflow-x:auto}td,th{border:1px solid #d7dce5;padding:.5rem;text-align:left}"
        "th{background:#f3f6fa}.warning{background:#fff4d6;border-left:5px solid #d99b20;"
        "padding:1rem}code{background:#f5f6f8;padding:.15rem .3rem}</style>"
        "<h1>Intervention evidence landscape</h1>"
        "<div class='warning'><b>Research boundary.</b> Records marked <b>*</b> are pending "
        "human review. This report does not recommend a drug, genome-editing treatment, or "
        "clinical action. An absent record means only that no record was located in the supplied "
        "evidence table.</div>"
        "<h2>Target-to-intervention landscape</h2>"
        + landscape.to_html(index=False, escape=True)
        + "<h2>Pending human review</h2>"
        + review_queue.to_html(index=False, escape=True)
        + "<h2>Reproducible search leads</h2><p>These database queries are discovery aids, not "
        "evidence. Source records must be reviewed before use.</p>"
        + queries.to_html(index=False, escape=True, render_links=True)
        + "<h2>Interpretation</h2><ul><li>Genetic-editing targetability and pharmacological "
        "tractability are separate axes.</li><li>Complex-level evidence does not establish that "
        "every subunit is a direct binding target.</li><li>Ortholog evidence is not direct "
        "evidence "
        "for the target virus.</li><li>Unknown values are not converted to favorable values.</li>"
        "</ul></html>",
        encoding="utf-8",
    )
    return {
        "evidence": evidence_output,
        "landscape": landscape_output,
        "review_queue": review_output,
        "queries": query_output,
        "report": report_output,
    }
