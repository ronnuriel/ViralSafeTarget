# Curated viral gene evidence

`viral_gene_evidence.schema.csv` defines the required columns.
`hsv2_gene_evidence.csv` is intentionally empty because the repository does not yet
contain verified gene-level evidence suitable for scoring.

To extend a table, verify the exact virus, reference accession, gene, claim, and
source. Use `supported`, `suggested`, `unknown`, or `conflicting` for essentiality and
latency status. Preserve uncertainty, identify the curator and date, and link a DOI,
URL, accession, or other stable source identifier. Then pass the file to
`vst scan --gene-evidence ...` or `vst rank --gene-evidence ...`.
