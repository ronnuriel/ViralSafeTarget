"""Source-linked evidence discovery with a mandatory human-review boundary.

The agent may propose interpretations, but proposals are never treated as curated
evidence. Only rows explicitly marked ``approved`` by a researcher can be exported
to the project's ``gene_evidence.tsv`` table.
"""

from __future__ import annotations

import hashlib
import html
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from .annotations import read_gff3

GENE_CATALOG_COLUMNS = [
    "gene_name",
    "feature_id",
    "locus_tag",
    "product",
    "protein_id",
    "aliases",
    "reference_accession",
    "target_virus",
    "virus_tax_id",
]

QUERY_COLUMNS = [
    "query_id",
    "gene_name",
    "query_family",
    "evidence_scope",
    "query_text",
]

SOURCE_RECORD_COLUMNS = [
    "source_database",
    "source_identifier",
    "source_url",
    "title",
    "abstract_or_annotation",
    "publication_year",
    "doi",
    "query_ids",
    "source_record_sha256",
]

PROPOSAL_COLUMNS = [
    "proposal_id",
    "gene_name",
    "matched_alias",
    "target_virus",
    "evidence_virus",
    "evidence_scope",
    "source_database",
    "source_identifier",
    "source_url",
    "source_title",
    "publication_year",
    "query_ids",
    "evidence_category",
    "experiment_type",
    "model_system",
    "measurement",
    "finding_summary",
    "evidence_direction",
    "proposed_essentiality_call",
    "proposed_essentiality_score",
    "proposed_evidence_strength",
    "quoted_evidence_span",
    "extraction_method",
    "confidence",
    "review_status",
    "reviewer",
    "review_date",
    "review_notes",
    "source_record_sha256",
    "created_utc",
]

APPROVED_EVIDENCE_COLUMNS = [
    "proposal_id",
    "gene_name",
    "virus_type",
    "reference_accession",
    "evidence_category",
    "essentiality_call",
    "essentiality_score",
    "evidence_strength",
    "experimental_system",
    "finding",
    "source_identifier",
    "source_title",
    "source_url",
    "quoted_evidence_span",
    "reviewer",
    "review_date",
    "review_notes",
    "directness_notes",
]

_EXPERIMENT_PATTERNS: list[tuple[str, tuple[str, ...]]] = [
    ("null_mutant", ("null mutant", "deletion mutant", "gene deletion", "knockout")),
    ("knockdown", ("knockdown", "sirna", "short interfering rna")),
    ("functional_mutagenesis", ("mutagenesis", "substitution mutant", "point mutant")),
    ("replication_phenotype", ("plaque assay", "viral titre", "viral titer", "replication")),
    ("protein_interaction", ("interacts with", "binds to", "interaction with")),
    ("expression", ("expression", "transcript", "upregulated", "downregulated")),
]

_MODEL_PATTERNS: list[tuple[str, tuple[str, ...]]] = [
    ("animal_model", ("mouse model", "mice", "in vivo", "animal model")),
    ("neuronal_or_latency_model", ("neuron", "ganglia", "latency model", "reactivation model")),
    ("cultured_cells", ("cell culture", "cultured cells", "vero cells", "cell line")),
    ("purified_or_biochemical", ("purified protein", "biochemical assay", "enzyme assay")),
]

_MEASUREMENT_PATTERNS: list[tuple[str, tuple[str, ...]]] = [
    ("viral_replication_or_titre", ("viral titre", "viral titer", "virus yield", "replication")),
    ("plaque_formation", ("plaque", "plaque-forming")),
    ("interferon_or_innate_immunity", ("interferon", "cgas", "sting", "innate immune")),
    ("latency_or_reactivation", ("latency", "latent", "reactivation")),
    ("protein_localization_or_abundance", ("localization", "expression", "protein level")),
]


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _normalized_space(value: Any) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()


def _words(value: str) -> list[str]:
    return re.findall(r"\S+", _normalized_space(value))


def _short_span(value: str, *, maximum_words: int = 25) -> str:
    words = _words(value)
    return " ".join(words[:maximum_words]) + (" …" if len(words) > maximum_words else "")


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _stable_id(prefix: str, *values: Any) -> str:
    material = "\x1f".join(_normalized_space(value) for value in values)
    return f"{prefix}-{_sha256_text(material)[:16]}"


def _all_text(element: ET.Element | None) -> str:
    return _normalized_space("".join(element.itertext())) if element is not None else ""


def _first_matching(text: str, patterns: Iterable[tuple[str, tuple[str, ...]]]) -> str:
    lowered = text.lower()
    return next(
        (label for label, terms in patterns if any(term in lowered for term in terms)), "unknown"
    )


def _profile_names(virus_profile: Mapping[str, Any], key: str) -> list[str]:
    value = virus_profile.get(key, [])
    if isinstance(value, str):
        value = [value]
    return sorted({_normalized_space(item) for item in value if _normalized_space(item)})


@dataclass
class HttpTransport:
    """Small cached HTTP transport with transparent provenance and retry behavior."""

    cache_dir: Path
    email: str = ""
    api_key: str = ""
    timeout_seconds: int = 45
    retries: int = 3
    offline: bool = False
    _last_request: float = 0.0

    def _cache_path(self, url: str) -> Path:
        return self.cache_dir / f"{_sha256_text(url)}.json"

    def get_json(self, base_url: str, parameters: Mapping[str, Any]) -> tuple[Any, dict[str, Any]]:
        values = {key: value for key, value in parameters.items() if value not in (None, "")}
        url = f"{base_url}?{urllib.parse.urlencode(values, doseq=True)}"
        cache_path = self._cache_path(url)
        if cache_path.is_file():
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            return cached["payload"], cached["provenance"]
        if self.offline:
            raise FileNotFoundError(f"Offline cache miss for {url}")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        minimum_interval = 0.11 if self.api_key else 0.35
        for attempt in range(self.retries):
            wait = minimum_interval - (time.monotonic() - self._last_request)
            if wait > 0:
                time.sleep(wait)
            request = urllib.request.Request(
                url,
                headers={
                    "Accept": "application/json",
                    "User-Agent": (
                        "ViralSafeTarget/0.8 evidence-agent "
                        + (self.email or "contact-not-configured")
                    ),
                },
            )
            try:
                with urllib.request.urlopen(request, timeout=self.timeout_seconds) as response:
                    body = response.read().decode("utf-8")
                    self._last_request = time.monotonic()
                    payload = json.loads(body)
                    provenance = {
                        "url": url,
                        "retrieved_utc": _utc_now(),
                        "http_status": int(response.status),
                        "content_sha256": _sha256_text(body),
                        "content_type": response.headers.get("Content-Type", ""),
                    }
                    cache_path.write_text(
                        json.dumps(
                            {"provenance": provenance, "payload": payload},
                            indent=2,
                            sort_keys=True,
                        ),
                        encoding="utf-8",
                    )
                    return payload, provenance
            except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
                if attempt + 1 == self.retries:
                    raise
                time.sleep(2**attempt)
        raise RuntimeError("Unreachable HTTP retry state")

    def get_xml(
        self, base_url: str, parameters: Mapping[str, Any]
    ) -> tuple[ET.Element, dict[str, Any]]:
        values = {key: value for key, value in parameters.items() if value not in (None, "")}
        url = f"{base_url}?{urllib.parse.urlencode(values, doseq=True)}"
        cache_path = self._cache_path(url)
        if cache_path.is_file():
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            return ET.fromstring(cached["payload"]), cached["provenance"]
        if self.offline:
            raise FileNotFoundError(f"Offline cache miss for {url}")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        minimum_interval = 0.11 if self.api_key else 0.35
        for attempt in range(self.retries):
            wait = minimum_interval - (time.monotonic() - self._last_request)
            if wait > 0:
                time.sleep(wait)
            request = urllib.request.Request(
                url,
                headers={
                    "User-Agent": f"ViralSafeTarget/0.8 {self.email or 'contact-not-configured'}"
                },
            )
            try:
                with urllib.request.urlopen(request, timeout=self.timeout_seconds) as response:
                    body = response.read().decode("utf-8")
                    self._last_request = time.monotonic()
                    root = ET.fromstring(body)
                    provenance = {
                        "url": url,
                        "retrieved_utc": _utc_now(),
                        "http_status": int(response.status),
                        "content_sha256": _sha256_text(body),
                        "content_type": response.headers.get("Content-Type", ""),
                    }
                    cache_path.write_text(
                        json.dumps(
                            {"provenance": provenance, "payload": body},
                            indent=2,
                            sort_keys=True,
                        ),
                        encoding="utf-8",
                    )
                    return root, provenance
            except (urllib.error.URLError, TimeoutError, ET.ParseError):
                if attempt + 1 == self.retries:
                    raise
                time.sleep(2**attempt)
        raise RuntimeError("Unreachable HTTP retry state")


def build_gene_catalog(gff_path: str | Path, virus_profile: Mapping[str, Any]) -> pd.DataFrame:
    """Build a generic alias-aware gene catalog from GFF attributes."""
    features = read_gff3(gff_path)
    if features.empty:
        return pd.DataFrame(columns=GENE_CATALOG_COLUMNS)
    rows: list[dict[str, Any]] = []
    for _, feature in features[features["feature_type"].isin(["gene", "CDS", "ncRNA"])].iterrows():
        attributes = feature.get("attributes") or {}
        gene_name = _normalized_space(
            attributes.get("gene")
            or attributes.get("Name")
            or attributes.get("locus_tag")
            or feature.get("name")
        )
        if not gene_name:
            continue
        aliases: set[str] = {gene_name}
        for key in (
            "gene",
            "Name",
            "locus_tag",
            "old_locus_tag",
            "protein_id",
            "ID",
            "gene_synonym",
            "Alias",
        ):
            raw = _normalized_space(attributes.get(key, ""))
            aliases.update(item.strip() for item in re.split(r"[,|]", raw) if item.strip())
        rows.append(
            {
                "gene_name": gene_name,
                "feature_id": _normalized_space(feature.get("feature_id")),
                "locus_tag": _normalized_space(attributes.get("locus_tag", "")),
                "product": _normalized_space(feature.get("product")),
                "protein_id": _normalized_space(attributes.get("protein_id", "")),
                "aliases": ";".join(sorted(aliases, key=str.casefold)),
                "reference_accession": _normalized_space(
                    virus_profile.get("reference_accession", feature.get("seqid", ""))
                ),
                "target_virus": _normalized_space(
                    virus_profile.get("scientific_name")
                    or virus_profile.get("display_name")
                    or virus_profile.get("id")
                ),
                "virus_tax_id": _normalized_space(virus_profile.get("tax_id", "")),
                "_priority": 0 if feature["feature_type"] == "CDS" else 1,
            }
        )
    catalog = pd.DataFrame(rows)
    if catalog.empty:
        return pd.DataFrame(columns=GENE_CATALOG_COLUMNS)
    catalog = catalog.sort_values(["gene_name", "_priority", "feature_id"], kind="mergesort")
    combined: list[dict[str, Any]] = []
    for gene_name, group in catalog.groupby("gene_name", sort=True):
        first = group.iloc[0]
        aliases = sorted(
            {alias for value in group["aliases"] for alias in str(value).split(";") if alias},
            key=str.casefold,
        )
        combined.append(
            {
                **{column: first.get(column, "") for column in GENE_CATALOG_COLUMNS},
                "gene_name": gene_name,
                "aliases": ";".join(aliases),
                "product": next((str(value) for value in group["product"] if str(value)), ""),
                "protein_id": next((str(value) for value in group["protein_id"] if str(value)), ""),
            }
        )
    return pd.DataFrame(combined, columns=GENE_CATALOG_COLUMNS)


def build_search_queries(
    catalog: pd.DataFrame,
    virus_profile: Mapping[str, Any],
) -> pd.DataFrame:
    """Generate transparent direct-virus and optional ortholog query plans."""
    direct_names = _profile_names(virus_profile, "literature_search_names")
    if not direct_names:
        direct_names = [
            _normalized_space(
                virus_profile.get("scientific_name")
                or virus_profile.get("display_name")
                or virus_profile.get("id")
            )
        ]
    ortholog_names = _profile_names(virus_profile, "ortholog_search_names")
    families = {
        "essentiality": (
            "deletion mutant",
            "null mutant",
            "knockout",
            "essential",
            "replication",
        ),
        "immune_evasion": ("interferon", "immune evasion", "cGAS", "STING"),
        "latency_reactivation": ("latency", "reactivation", "latent infection"),
        "gene_function": ("function", "protein", "mutagenesis", "interaction"),
    }
    rows: list[dict[str, str]] = []
    for _, gene in catalog.iterrows():
        gene_terms = [term for term in str(gene["aliases"]).split(";") if term]
        if gene.get("product"):
            gene_terms.append(str(gene["product"]))
        gene_clause = " OR ".join(f'"{term}"' for term in sorted(set(gene_terms)))
        for scope, names in (("direct_target_virus", direct_names), ("ortholog", ortholog_names)):
            if not names:
                continue
            virus_clause = " OR ".join(f'"{name}"' for name in names)
            for family, terms in families.items():
                evidence_clause = " OR ".join(f'"{term}"' for term in terms)
                query = f"({virus_clause}) AND ({gene_clause}) AND ({evidence_clause})"
                rows.append(
                    {
                        "query_id": _stable_id("Q", gene["gene_name"], family, scope, query),
                        "gene_name": str(gene["gene_name"]),
                        "query_family": family,
                        "evidence_scope": scope,
                        "query_text": query,
                    }
                )
    return pd.DataFrame(rows, columns=QUERY_COLUMNS)


def _pubmed_records(
    queries: pd.DataFrame,
    transport: HttpTransport,
    *,
    maximum_results_per_query: int,
) -> list[dict[str, Any]]:
    endpoint = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    ids_to_queries: dict[str, set[str]] = {}
    for _, row in queries.iterrows():
        payload, _ = transport.get_json(
            f"{endpoint}/esearch.fcgi",
            {
                "db": "pubmed",
                "term": row["query_text"],
                "retmode": "json",
                "retmax": maximum_results_per_query,
                "sort": "relevance",
                "tool": "ViralSafeTarget",
                "email": transport.email,
                "api_key": transport.api_key,
            },
        )
        for identifier in payload.get("esearchresult", {}).get("idlist", []):
            ids_to_queries.setdefault(str(identifier), set()).add(str(row["query_id"]))
    records: list[dict[str, Any]] = []
    identifiers = sorted(ids_to_queries, key=lambda value: int(value))
    for start in range(0, len(identifiers), 100):
        block = identifiers[start : start + 100]
        root, _ = transport.get_xml(
            f"{endpoint}/efetch.fcgi",
            {
                "db": "pubmed",
                "id": ",".join(block),
                "retmode": "xml",
                "tool": "ViralSafeTarget",
                "email": transport.email,
                "api_key": transport.api_key,
            },
        )
        for article in root.findall(".//PubmedArticle"):
            pmid = _all_text(article.find(".//PMID"))
            title = _all_text(article.find(".//ArticleTitle"))
            abstracts = [_all_text(item) for item in article.findall(".//Abstract/AbstractText")]
            abstract = _normalized_space(" ".join(item for item in abstracts if item))
            year = (
                _all_text(article.find(".//PubDate/Year"))
                or _all_text(article.find(".//PubDate/MedlineDate"))[:4]
            )
            doi = next(
                (
                    _all_text(item)
                    for item in article.findall(".//ArticleId")
                    if item.attrib.get("IdType") == "doi"
                ),
                "",
            )
            records.append(
                {
                    "source_database": "PubMed",
                    "source_identifier": f"PMID:{pmid}",
                    "source_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    "title": title,
                    "abstract_or_annotation": abstract,
                    "publication_year": year,
                    "doi": doi,
                    "query_ids": ";".join(sorted(ids_to_queries.get(pmid, set()))),
                }
            )
    return records


def _europe_pmc_records(
    queries: pd.DataFrame,
    transport: HttpTransport,
    *,
    maximum_results_per_query: int,
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    endpoint = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    for _, row in queries.iterrows():
        payload, _ = transport.get_json(
            endpoint,
            {
                "query": row["query_text"],
                "format": "json",
                "resultType": "core",
                "pageSize": maximum_results_per_query,
            },
        )
        for result in payload.get("resultList", {}).get("result", []):
            pmid = _normalized_space(result.get("pmid", ""))
            pmcid = _normalized_space(result.get("pmcid", ""))
            identifier = f"PMID:{pmid}" if pmid else f"PMCID:{pmcid}" if pmcid else ""
            if not identifier:
                continue
            url_id = pmid or pmcid
            records.append(
                {
                    "source_database": "Europe PMC",
                    "source_identifier": identifier,
                    "source_url": f"https://europepmc.org/article/MED/{url_id}",
                    "title": _normalized_space(result.get("title", "")),
                    "abstract_or_annotation": _normalized_space(result.get("abstractText", "")),
                    "publication_year": _normalized_space(result.get("pubYear", "")),
                    "doi": _normalized_space(result.get("doi", "")),
                    "query_ids": str(row["query_id"]),
                }
            )
    return records


def _uniprot_records(
    catalog: pd.DataFrame,
    transport: HttpTransport,
    *,
    maximum_results_per_gene: int,
) -> list[dict[str, Any]]:
    endpoint = "https://rest.uniprot.org/uniprotkb/search"
    records: list[dict[str, Any]] = []
    for _, gene in catalog.iterrows():
        aliases = [item for item in str(gene["aliases"]).split(";") if item]
        alias_clause = " OR ".join(f"gene_exact:{alias}" for alias in aliases[:8])
        query = f"({alias_clause})"
        if gene.get("virus_tax_id"):
            query += f" AND organism_id:{gene['virus_tax_id']}"
        payload, _ = transport.get_json(
            endpoint,
            {
                "query": query,
                "format": "json",
                "size": maximum_results_per_gene,
                "fields": (
                    "accession,id,protein_name,gene_names,organism_name,ft_domain,"
                    "cc_function,cc_subcellular_location,lit_pubmed_id"
                ),
            },
        )
        for result in payload.get("results", []):
            accession = _normalized_space(result.get("primaryAccession", ""))
            description = result.get("proteinDescription", {})
            recommended = description.get("recommendedName", {}).get("fullName", {}).get("value")
            submitted = description.get("submissionNames", [{}])[0].get("fullName", {}).get("value")
            protein_name = _normalized_space(recommended or submitted or result.get("uniProtkbId"))
            comments = []
            for comment in result.get("comments", []):
                for text_item in comment.get("texts", []):
                    if text_item.get("value"):
                        comments.append(str(text_item["value"]))
            features = []
            for feature in result.get("features", []):
                if feature.get("type") == "Domain":
                    features.append(_normalized_space(feature.get("description", "domain")))
            annotation = _normalized_space(
                " ".join(
                    [
                        *(f"Function: {item}" for item in comments),
                        *(f"Domain: {item}" for item in features),
                    ]
                )
            )
            records.append(
                {
                    "source_database": "UniProt",
                    "source_identifier": f"UniProt:{accession}",
                    "source_url": f"https://www.uniprot.org/uniprotkb/{accession}/entry",
                    "title": protein_name,
                    "abstract_or_annotation": annotation,
                    "publication_year": "",
                    "doi": "",
                    "query_ids": _stable_id("UQ", gene["gene_name"], query),
                    "_gene_name": str(gene["gene_name"]),
                }
            )
    return records


def _deduplicate_records(records: Iterable[dict[str, Any]]) -> pd.DataFrame:
    grouped: dict[str, dict[str, Any]] = {}
    priority = {"PubMed": 0, "Europe PMC": 1, "UniProt": 2}
    for record in records:
        identifier = str(record["source_identifier"])
        key = (
            identifier
            if identifier.startswith(("PMID:", "PMCID:"))
            else (f"{record['source_database']}:{identifier}")
        )
        existing = grouped.get(key)
        if existing is None or priority.get(str(record["source_database"]), 9) < priority.get(
            str(existing["source_database"]), 9
        ):
            query_ids = set(str(record.get("query_ids", "")).split(";"))
            selected = dict(record)
            selected["_query_id_set"] = {item for item in query_ids if item}
            grouped[key] = selected
        else:
            existing["_query_id_set"].update(
                item for item in str(record.get("query_ids", "")).split(";") if item
            )
    output: list[dict[str, Any]] = []
    for record in grouped.values():
        record["query_ids"] = ";".join(sorted(record.pop("_query_id_set", set())))
        normalized = {column: record.get(column, "") for column in SOURCE_RECORD_COLUMNS}
        source_material = json.dumps(normalized, sort_keys=True, ensure_ascii=False)
        normalized["source_record_sha256"] = _sha256_text(source_material)
        if record.get("_gene_name"):
            normalized["_gene_name"] = record["_gene_name"]
        output.append(normalized)
    return (
        pd.DataFrame(output)
        .sort_values(["source_database", "source_identifier"], kind="mergesort")
        .reset_index(drop=True)
    )


def _matched_alias(text: str, aliases: Iterable[str]) -> str:
    for alias in sorted(set(aliases), key=lambda value: (-len(value), value.casefold())):
        if re.search(rf"(?<![A-Za-z0-9]){re.escape(alias)}(?![A-Za-z0-9])", text, re.I):
            return alias
    return ""


def _name_match(text: str, names: Iterable[str]) -> str:
    return next((name for name in names if name.casefold() in text.casefold()), "")


def _evidence_virus(
    title: str,
    body: str,
    virus_profile: Mapping[str, Any],
) -> tuple[str, str]:
    direct_names = _profile_names(virus_profile, "literature_search_names")
    ortholog_names = _profile_names(virus_profile, "ortholog_search_names")
    target = _normalized_space(
        virus_profile.get("scientific_name")
        or virus_profile.get("display_name")
        or virus_profile.get("id")
    )
    title_direct = _name_match(title, direct_names)
    title_ortholog = _name_match(title, ortholog_names)
    if title_direct and not title_ortholog:
        return target, "direct_target_virus"
    if title_ortholog and not title_direct:
        return title_ortholog, "ortholog"
    if re.search(r"\b(?:virus|herpesvirus|alphaherpesvirus|cytomegalovirus)\b", title, re.I):
        return "other or unresolved virus named in title", "mixed_or_other_virus"
    body_direct = _name_match(body, direct_names)
    body_ortholog = _name_match(body, ortholog_names)
    if body_direct and not body_ortholog:
        return target, "direct_target_virus"
    if body_ortholog and not body_direct:
        return body_ortholog, "ortholog"
    if body_direct and body_ortholog:
        return "multiple viruses mentioned", "mixed_or_other_virus"
    return "unresolved", "unresolved"


def _direction_and_call(text: str, experiment_type: str) -> tuple[str, str, str]:
    lowered = text.lower()
    nonessential = any(
        phrase in lowered
        for phrase in (
            "not essential",
            "nonessential",
            "non-essential",
            "dispensable",
            "replicated normally",
        )
    )
    required = any(
        phrase in lowered
        for phrase in (
            "essential",
            "required for replication",
            "failed to replicate",
            "unable to replicate",
            "abolished replication",
        )
    )
    if experiment_type == "null_mutant" and nonessential:
        return (
            "supports_nonessentiality_in_reported_context",
            "nonessential_in_reported_context",
            "",
        )
    if experiment_type in {"null_mutant", "functional_mutagenesis"} and required:
        return (
            "supports_requirement_in_reported_context",
            "supported_required_in_reported_context",
            "",
        )
    if experiment_type == "expression":
        return "association_or_expression_only", "unknown", ""
    return "phenotype_or_function_reported", "unknown", ""


def extract_proposals(
    records: pd.DataFrame,
    catalog: pd.DataFrame,
    queries: pd.DataFrame,
    virus_profile: Mapping[str, Any],
) -> pd.DataFrame:
    """Create conservative, review-pending proposals from normalized records."""
    query_to_gene = dict(zip(queries["query_id"], queries["gene_name"], strict=False))
    catalog_by_gene = catalog.set_index("gene_name", drop=False)
    target_virus = _normalized_space(
        virus_profile.get("scientific_name")
        or virus_profile.get("display_name")
        or virus_profile.get("id")
    )
    proposals: list[dict[str, Any]] = []
    for _, record in records.iterrows():
        record_text = _normalized_space(
            f"{record.get('title', '')} {record.get('abstract_or_annotation', '')}"
        )
        candidate_genes = {
            query_to_gene[query_id]
            for query_id in str(record.get("query_ids", "")).split(";")
            if query_id in query_to_gene
        }
        if record.get("_gene_name"):
            candidate_genes.add(str(record["_gene_name"]))
        for gene_name in sorted(candidate_genes):
            if gene_name not in catalog_by_gene.index:
                continue
            gene = catalog_by_gene.loc[gene_name]
            aliases = [item for item in str(gene["aliases"]).split(";") if item]
            matched_alias = _matched_alias(record_text, aliases)
            if not matched_alias and record["source_database"] != "UniProt":
                continue
            evidence_virus, scope = _evidence_virus(
                str(record.get("title", "")),
                str(record.get("abstract_or_annotation", "")),
                virus_profile,
            )
            if record["source_database"] == "UniProt" and scope == "unresolved":
                evidence_virus, scope = target_virus, "direct_target_virus"
            experiment = (
                "database_annotation"
                if record["source_database"] == "UniProt"
                else _first_matching(record_text, _EXPERIMENT_PATTERNS)
            )
            model = _first_matching(record_text, _MODEL_PATTERNS)
            measurement = _first_matching(record_text, _MEASUREMENT_PATTERNS)
            category = (
                "protein_annotation"
                if record["source_database"] == "UniProt"
                else "immune_evasion"
                if measurement == "interferon_or_innate_immunity"
                else "latency_reactivation"
                if measurement == "latency_or_reactivation"
                else "essentiality_or_replication"
                if experiment
                in {
                    "null_mutant",
                    "knockdown",
                    "functional_mutagenesis",
                    "replication_phenotype",
                }
                else "gene_function"
            )
            direction, essentiality_call, score = _direction_and_call(record_text, experiment)
            strength = (
                "direct"
                if scope == "direct_target_virus"
                and experiment in {"null_mutant", "knockdown", "functional_mutagenesis"}
                else "supporting"
                if scope == "direct_target_virus"
                else "indirect"
                if scope == "ortholog"
                else "unknown"
            )
            confidence = (
                "medium"
                if matched_alias
                and scope != "unresolved"
                and experiment not in {"unknown", "expression"}
                else "low"
            )
            finding_source = str(record.get("abstract_or_annotation") or record.get("title"))
            quoted = _short_span(finding_source)
            proposal_id = _stable_id(
                "EP", gene_name, record["source_identifier"], category, experiment
            )
            proposals.append(
                {
                    "proposal_id": proposal_id,
                    "gene_name": gene_name,
                    "matched_alias": matched_alias,
                    "target_virus": target_virus,
                    "evidence_virus": evidence_virus,
                    "evidence_scope": scope,
                    "source_database": record["source_database"],
                    "source_identifier": record["source_identifier"],
                    "source_url": record["source_url"],
                    "source_title": record["title"],
                    "publication_year": record["publication_year"],
                    "query_ids": record["query_ids"],
                    "evidence_category": category,
                    "experiment_type": experiment,
                    "model_system": model,
                    "measurement": measurement,
                    "finding_summary": quoted,
                    "evidence_direction": direction,
                    "proposed_essentiality_call": essentiality_call,
                    "proposed_essentiality_score": score,
                    "proposed_evidence_strength": strength,
                    "quoted_evidence_span": quoted,
                    "extraction_method": "conservative_rules_v1",
                    "confidence": confidence,
                    "review_status": "pending",
                    "reviewer": "",
                    "review_date": "",
                    "review_notes": "",
                    "source_record_sha256": record["source_record_sha256"],
                    "created_utc": _utc_now(),
                }
            )
    if not proposals:
        return pd.DataFrame(columns=PROPOSAL_COLUMNS)
    return (
        pd.DataFrame(proposals, columns=PROPOSAL_COLUMNS)
        .drop_duplicates("proposal_id")
        .sort_values(["gene_name", "source_identifier", "proposal_id"], kind="mergesort")
        .reset_index(drop=True)
    )


def fetch_refseq_metadata(
    reference_accession: str,
    transport: HttpTransport,
) -> dict[str, Any]:
    endpoint = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    search, search_provenance = transport.get_json(
        f"{endpoint}/esearch.fcgi",
        {
            "db": "nuccore",
            "term": f"{reference_accession}[Accession]",
            "retmode": "json",
            "retmax": 1,
            "tool": "ViralSafeTarget",
            "email": transport.email,
            "api_key": transport.api_key,
        },
    )
    identifiers = search.get("esearchresult", {}).get("idlist", [])
    if not identifiers:
        return {
            "reference_accession": reference_accession,
            "status": "not_found",
            "search_provenance": search_provenance,
        }
    summary, summary_provenance = transport.get_json(
        f"{endpoint}/esummary.fcgi",
        {
            "db": "nuccore",
            "id": identifiers[0],
            "retmode": "json",
            "tool": "ViralSafeTarget",
            "email": transport.email,
            "api_key": transport.api_key,
        },
    )
    result = summary.get("result", {}).get(str(identifiers[0]), {})
    return {
        "reference_accession": reference_accession,
        "status": "found",
        "uid": str(identifiers[0]),
        "caption": result.get("caption", ""),
        "title": result.get("title", ""),
        "tax_id": result.get("taxid", ""),
        "sequence_length": result.get("slen", ""),
        "updated_date": result.get("updatedate", ""),
        "search_provenance": search_provenance,
        "summary_provenance": summary_provenance,
    }


def _write_jsonl(frame: pd.DataFrame, path: Path) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for record in frame.to_dict("records"):
            handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True) + "\n")


def _table(frame: pd.DataFrame, columns: list[str], rows: int = 200) -> str:
    available = [column for column in columns if column in frame]
    if frame.empty or not available:
        return "<p>No records.</p>"
    return frame[available].head(rows).to_html(index=False, escape=True, classes="data")


def write_evidence_report(
    proposals: pd.DataFrame,
    catalog: pd.DataFrame,
    queries: pd.DataFrame,
    output: str | Path,
) -> Path:
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    statuses = proposals["review_status"].value_counts().to_dict() if not proposals.empty else {}
    genes_with_proposals = proposals["gene_name"].nunique() if not proposals.empty else 0
    review_columns = [
        "proposal_id",
        "gene_name",
        "evidence_virus",
        "evidence_scope",
        "evidence_category",
        "experiment_type",
        "evidence_direction",
        "proposed_essentiality_call",
        "quoted_evidence_span",
        "confidence",
        "review_status",
        "source_identifier",
        "source_url",
    ]
    style = (
        "body{font:15px/1.5 system-ui,sans-serif;max-width:1500px;margin:auto;"
        "padding:2rem;color:#172b4d}table{border-collapse:collapse;width:100%;"
        "display:block;overflow:auto;font-size:.82rem}th,td{border:1px solid #d8e2ec;"
        "padding:.4rem;text-align:left;vertical-align:top}th{background:#eef4f8}"
        ".warning{background:#fff4e5;border-left:5px solid #d97706;padding:1rem}"
        "code{background:#f4f7f9;padding:.1rem .3rem}"
    )
    warning = (
        "<div class='warning'><strong>Human review is mandatory.</strong> Proposals do "
        "not alter gene scores or become curated evidence until a researcher marks them "
        "<code>approved</code> and runs the explicit apply command. Expression is not "
        "essentiality; ortholog evidence is not direct target-virus evidence; absence of "
        "located evidence is not non-essentiality.</div>"
    )
    summary = (
        f"<p>Genes in annotation: {len(catalog):,}; queries: {len(queries):,}; "
        f"proposals: {len(proposals):,}; genes with proposals: {genes_with_proposals:,}; "
        f"review status: {html.escape(json.dumps(statuses, sort_keys=True))}.</p>"
    )
    boundary = (
        "<h2>Interpretation boundary</h2><p>This report supports literature triage. "
        "It does not establish gene essentiality, editing, viral inhibition, host safety, "
        "delivery, or therapeutic efficacy, and it contains no wet-lab protocol.</p>"
    )
    document = "".join(
        [
            "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
            "<meta name='viewport' content='width=device-width,initial-scale=1'>",
            f"<title>Evidence review queue</title><style>{style}</style></head><body>",
            "<h1>Source-linked evidence proposals</h1>",
            warning,
            summary,
            "<h2>Review queue</h2>",
            _table(proposals, review_columns),
            "<h2>Gene catalog</h2>",
            _table(catalog, GENE_CATALOG_COLUMNS),
            "<h2>Generated query plan</h2>",
            _table(queries, QUERY_COLUMNS),
            boundary,
            "</body></html>",
        ]
    )
    destination.write_text(document, encoding="utf-8")
    return destination


def discover_evidence(
    *,
    gff_path: str | Path,
    virus_profile: Mapping[str, Any],
    out_dir: str | Path,
    sources: Iterable[str] = ("pubmed", "europepmc", "uniprot", "ncbi_refseq"),
    maximum_results_per_query: int = 5,
    email: str = "",
    api_key: str = "",
    offline: bool = False,
    genes: Iterable[str] | None = None,
) -> dict[str, Any]:
    """Run source discovery and emit review-pending evidence artifacts."""
    output = Path(out_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    selected_sources = {source.lower() for source in sources}
    transport = HttpTransport(
        output / "cache",
        email=email or os.environ.get("NCBI_EMAIL", ""),
        api_key=api_key or os.environ.get("NCBI_API_KEY", ""),
        offline=offline,
    )
    catalog = build_gene_catalog(gff_path, virus_profile)
    selected_genes = {str(gene) for gene in genes or []}
    if selected_genes:
        missing_genes = sorted(selected_genes - set(catalog["gene_name"].astype(str)))
        if missing_genes:
            raise ValueError("Genes are absent from the annotation: " + ", ".join(missing_genes))
        catalog = catalog[catalog["gene_name"].astype(str).isin(selected_genes)].reset_index(
            drop=True
        )
    queries = build_search_queries(catalog, virus_profile)
    literature_queries = queries.copy()
    records: list[dict[str, Any]] = []
    source_errors: list[dict[str, str]] = []
    for source, action in (
        (
            "pubmed",
            lambda: _pubmed_records(
                literature_queries,
                transport,
                maximum_results_per_query=maximum_results_per_query,
            ),
        ),
        (
            "europepmc",
            lambda: _europe_pmc_records(
                literature_queries,
                transport,
                maximum_results_per_query=maximum_results_per_query,
            ),
        ),
        (
            "uniprot",
            lambda: _uniprot_records(
                catalog,
                transport,
                maximum_results_per_gene=maximum_results_per_query,
            ),
        ),
    ):
        if source not in selected_sources:
            continue
        try:
            records.extend(action())
        except Exception as error:  # Preserve partial source coverage explicitly.
            source_errors.append({"source": source, "error": str(error)})
    normalized_records = (
        _deduplicate_records(records) if records else pd.DataFrame(columns=SOURCE_RECORD_COLUMNS)
    )
    proposals = extract_proposals(normalized_records, catalog, queries, virus_profile)
    metadata: dict[str, Any] = {}
    if "ncbi_refseq" in selected_sources:
        try:
            metadata = fetch_refseq_metadata(
                str(virus_profile.get("reference_accession", "")), transport
            )
        except Exception as error:
            source_errors.append({"source": "ncbi_refseq", "error": str(error)})
            metadata = {"status": "error", "error": str(error)}

    catalog.to_csv(output / "gene_catalog.tsv", sep="\t", index=False)
    queries.to_csv(output / "search_queries.tsv", sep="\t", index=False)
    _write_jsonl(normalized_records, output / "source_records.jsonl")
    proposals.to_csv(output / "evidence_proposals.tsv", sep="\t", index=False)
    proposals.to_csv(output / "review_queue.tsv", sep="\t", index=False)
    (output / "virus_metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True), encoding="utf-8"
    )
    manifest = {
        "schema_version": "1.0",
        "created_utc": _utc_now(),
        "sources_requested": sorted(selected_sources),
        "source_errors": source_errors,
        "gene_count": len(catalog),
        "query_count": len(queries),
        "source_record_count": len(normalized_records),
        "proposal_count": len(proposals),
        "approved_count": 0,
        "automatic_score_integration": False,
        "interpretation": (
            "All rows are review-pending proposals. No essentiality, function, or score "
            "is assigned without explicit researcher approval."
        ),
    }
    (output / "evidence_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    write_evidence_report(proposals, catalog, queries, output / "evidence_review_report.html")
    return {"output_dir": output, **manifest}


def apply_reviewed_evidence(
    review_queue: str | Path,
    output_table: str | Path,
    *,
    append: bool = False,
    reference_accession: str = "",
) -> dict[str, Any]:
    """Export only explicitly approved proposals to the curated evidence schema."""
    review_path = Path(review_queue)
    frame = pd.read_csv(review_path, sep="\t", dtype=str).fillna("")
    missing = [column for column in PROPOSAL_COLUMNS if column not in frame]
    if missing:
        raise ValueError("Review queue is missing columns: " + ", ".join(missing))
    allowed = {"pending", "approved", "rejected", "needs_revision"}
    invalid = sorted(set(frame["review_status"]) - allowed)
    if invalid:
        raise ValueError("Invalid review_status values: " + ", ".join(invalid))
    approved = frame[frame["review_status"].eq("approved")].copy()
    if not approved.empty:
        incomplete = approved[
            approved["reviewer"].str.strip().eq("")
            | approved["review_date"].str.strip().eq("")
            | approved["source_url"].str.strip().eq("")
        ]
        if not incomplete.empty:
            raise ValueError(
                "Approved proposals require reviewer, review_date, and source_url: "
                + ", ".join(incomplete["proposal_id"].astype(str))
            )
    rows: list[dict[str, Any]] = []
    for _, proposal in approved.iterrows():
        score = proposal["proposed_essentiality_score"].strip()
        rows.append(
            {
                "proposal_id": proposal["proposal_id"],
                "gene_name": proposal["gene_name"],
                "virus_type": proposal["evidence_virus"],
                "reference_accession": reference_accession,
                "evidence_category": proposal["evidence_category"],
                "essentiality_call": proposal["proposed_essentiality_call"] or "unknown",
                "essentiality_score": float(score) if score else pd.NA,
                "evidence_strength": proposal["proposed_evidence_strength"] or "unknown",
                "experimental_system": ";".join(
                    item
                    for item in (proposal["experiment_type"], proposal["model_system"])
                    if item and item != "unknown"
                ),
                "finding": proposal["finding_summary"],
                "source_identifier": proposal["source_identifier"],
                "source_title": proposal["source_title"],
                "source_url": proposal["source_url"],
                "quoted_evidence_span": proposal["quoted_evidence_span"],
                "reviewer": proposal["reviewer"],
                "review_date": proposal["review_date"],
                "review_notes": proposal["review_notes"],
                "directness_notes": (
                    f"scope={proposal['evidence_scope']}; reviewed by {proposal['reviewer']} "
                    f"on {proposal['review_date']}; {proposal['review_notes']}"
                ).strip(),
            }
        )
    destination = Path(output_table)
    destination.parent.mkdir(parents=True, exist_ok=True)
    output = pd.DataFrame(rows, columns=APPROVED_EVIDENCE_COLUMNS)
    if append and destination.is_file():
        existing = pd.read_csv(destination, sep="\t", dtype=str)
        output = pd.concat([existing, output], ignore_index=True)
    if not output.empty:
        output = output.drop_duplicates(
            ["gene_name", "virus_type", "source_identifier", "evidence_category"],
            keep="last",
        ).sort_values(
            ["gene_name", "virus_type", "source_identifier", "evidence_category"],
            kind="mergesort",
        )
    output.to_csv(destination, sep="\t", index=False)
    summary = {
        "reviewed_proposal_count": len(frame),
        "approved_count": len(approved),
        "rejected_count": int(frame["review_status"].eq("rejected").sum()),
        "pending_count": int(frame["review_status"].eq("pending").sum()),
        "needs_revision_count": int(frame["review_status"].eq("needs_revision").sum()),
        "output_table": str(destination.resolve()),
    }
    destination.with_suffix(".review_manifest.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8"
    )
    return summary
