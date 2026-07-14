# Machine-readable schemas

Schemas define the generic project and human-review contracts:

- `project.schema.json` — top-level project configuration, including bounded virtual
  knockout and configured multiplex strategies.
- `virus_profile.schema.json` — virus identity, inputs, aliases, evidence paths, and
  optional domain/disorder/conserved-region/category tables.
- `host_profile.schema.json` — host assembly configuration.
- `nuclease_profile.schema.json` — editor, PAM, and cut model.
- `gene_evidence.schema.json` — curated source-linked gene evidence.
- `evidence_proposal.schema.json` — machine-proposed, review-pending evidence row.
- `evidence_review.schema.json` — researcher review decision and provenance.
- `tool_benchmark.schema.json` — frozen-panel tool comparison, execution status, and
  ablation configuration.

The Python workflow validates the profiles it consumes. Missing biological evidence
must remain explicit and must not be converted to zero or a positive score.
