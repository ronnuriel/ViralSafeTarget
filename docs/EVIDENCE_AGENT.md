# Evidence Agent with mandatory researcher review

ViralSafeTarget can discover source-linked literature and protein annotations for genes
in a virus project. Discovery is deliberately separated from curation:

```text
official APIs -> normalized source records -> evidence proposals -> researcher review
             -> approved rows only -> gene_evidence.tsv -> evidence-aware rerun
```

An automatically discovered record is never treated as biological truth. Expression is
not essentiality, an interaction is not necessarily immune inhibition, evidence from an
ortholog is not direct evidence for the target virus, and failure to locate a paper is
not evidence of non-essentiality.

## Official sources

- [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/home/develop/api/) for PubMed and
  RefSeq metadata.
- [Europe PMC REST service](https://europepmc.org/RestfulWebService) for literature
  metadata and abstracts when available.
- [UniProt REST API](https://www.uniprot.org/help/programmatic) for structured protein
  annotations and source links.
- [NCBI Datasets Virus](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/virus/genome/)
  for reference and virus dataset acquisition outside the literature stage.

API responses are cached by request URL and content checksum. Set `NCBI_EMAIL` and,
optionally, `NCBI_API_KEY` in the environment for regular NCBI use. Do not commit API
keys or raw literature caches.

## 1. Configure virus identity

The virus profile accepts generic identity and search fields:

```yaml
scientific_name: Human alphaherpesvirus 2
tax_id: 10310
literature_search_names:
  - Human alphaherpesvirus 2
  - Herpes simplex virus 2
  - HSV-2
ortholog_search_names:
  - Human alphaherpesvirus 1
  - Herpes simplex virus 1
  - HSV-1
```

Genes, locus tags, products, protein accessions and aliases are read from the project
GFF. No virus-specific gene name is hard-coded in the agent.

## 2. Discover proposals

```bash
vst evidence discover \
  --project project.yaml \
  --max-results-per-query 5
```

For a focused review:

```bash
vst evidence discover \
  --project project.yaml \
  --genes UL3 UL30 UL52 \
  --max-results-per-query 10
```

The default output directory is `results/evidence/`:

- `gene_catalog.tsv`: annotation-derived genes and aliases.
- `search_queries.tsv`: exact generated query plan and evidence scope.
- `source_records.jsonl`: normalized source cache with record checksums.
- `evidence_proposals.tsv`: conservative extraction proposals.
- `review_queue.tsv`: editable researcher review queue.
- `virus_metadata.json`: reference lookup and API provenance.
- `evidence_manifest.json`: source coverage, errors and counts.
- `evidence_review_report.html`: linked human-readable review report.

Partial API failure is explicit in the manifest; it is never converted to “no evidence.”
Use `--offline` to reproduce from an existing cache without network access.

## 3. Review

Open `review_queue.tsv`. For every row, inspect the linked source and its context. Valid
review states are:

- `pending`
- `approved`
- `rejected`
- `needs_revision`

An approved row must contain `reviewer`, `review_date` and `source_url`. The researcher
may correct the proposed category, context, finding, direction or essentiality call
before approval. The quoted span is intentionally short; it is a navigation aid, not a
replacement for reading the source.

## 4. Apply approved rows explicitly

```bash
vst evidence apply \
  --project project.yaml \
  --review-queue results/evidence/review_queue.tsv
```

Only `approved` rows are exported to the virus profile's `gene_evidence.tsv`. Pending,
rejected and revision-needed rows cannot enter the curated table. HSV-1/ortholog scope
is retained as ortholog evidence and does not become direct HSV-2 evidence.

Then rerun:

```bash
vst project resume --project project.yaml
```

The evidence table checksum invalidates the relevant cached ranking stage. Source-linked
evidence coverage can then change, while sequence targetability, host-risk predictions
and predicted protein disruption remain separately visible.

## Scientific boundary

The evidence agent performs literature triage, not autonomous scientific adjudication.
It does not establish essentiality, editing, viral inhibition, safety, delivery or
clinical efficacy. It emits no wet-lab protocol. All computational proposals remain
hypotheses until independently reviewed and experimentally tested.
