# Generic research profiles

ViralSafeTarget separates reusable program logic from case-study data. A run is
declared by three versioned YAML profiles:

- a **virus profile** with references, annotations, strain alignment, curated
  evidence and standardized upstream outputs;
- a **host profile** with an assembly identity and the external reference path;
- a **nuclease profile** with the tested sequence model, PAM and mismatch limit.

The bundled HSV-2 case study uses:

```text
configs/viruses/hsv2.yaml
configs/hosts/human_grch38.yaml
configs/nucleases/spcas9.yaml
```

Validate them before a report run:

```bash
vst profiles validate \
  --virus-profile configs/viruses/hsv2.yaml \
  --host-profile configs/hosts/human_grch38.yaml \
  --nuclease-profile configs/nucleases/spcas9.yaml
```

Use `--require-host-reference` when a genome-scale host FASTA must be present.
Without it, a missing large external reference is reported as `external_pending`
instead of being silently accepted.

## Add a virus

1. Add reference FASTA, GFF3 and an aligned strain panel.
2. Add a virus YAML matching `schemas/virus_profile.schema.json`.
3. Add a source-linked evidence table matching
   `schemas/gene_evidence.schema.json`.
4. Keep organism-specific categories and evidence in data files, never in Python
   conditionals.
5. Validate the profile and record all checksums in the resulting manifest.

An optional held-out population panel should exclude every accession used during
discovery. Public partial genomes require a locus-specific observability denominator;
whole-record absence alone is not evidence of a variant target. Population support is
reported separately and must not silently alter targetability or essentiality scores.

Evidence absence must remain `unknown`. Evidence from a related virus is stored
with its own `virus_type` and cannot fill a direct-evidence field for the target
virus.

The profiles are an integration contract, not a claim that every virus has the
same biology or that every nuclease model has been experimentally validated.
