# Adding a new virus

Start from the generated project contract rather than copying an HSV-specific
script:

```bash
vst project init \
  --id my-virus \
  --display-name "My virus" \
  --reference-accession REF_ACCESSION \
  --out-dir projects/my-virus
```

1. Choose and freeze a reference accession and a defensible public/local genome set.
2. Apply documented QC and retain accession-level acceptance/rejection reasons.
3. Produce an aligned multi-FASTA. Every record must have the same alignment length.
4. Supply a reference-matched GFF3 if annotation is available.
5. Copy `configs/research_v0.3.yaml`, adjust only documented thresholds, and retain it
   with the run manifest.
6. Run `vst scan`, inspect the alignment and output funnel, then build a full-host
   off-target input with an assembly identifier.
7. Create a held-out benchmark template rather than tuning against known positives.

Then use the canonical workflow:

```bash
vst project validate --project projects/my-virus/project.yaml
vst project run --project projects/my-virus/project.yaml
vst project status --project projects/my-virus/project.yaml
```

See [`NEW_VIRUS_WORKFLOW.md`](NEW_VIRUS_WORKFLOW.md) for host-screen execution,
resume semantics, and the current boundary of generic protein analysis.

FASTA contains nucleotide sequences. GFF3 contains coordinate annotations and
normally does not contain the genome sequence. Their reference identifiers and
coordinate systems must match. Annotation indicates location; it does not establish
essentiality, accessibility, or effectiveness.
