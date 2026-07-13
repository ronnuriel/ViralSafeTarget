# Architecture

## Design goals

1. **Virus-first:** begin with many viral strains rather than one preselected locus.
2. **Auditable:** retain coordinates, score components and rejection reasons.
3. **Composable:** use established alignment and off-target engines instead of replacing them with an opaque score.
4. **Reproducible:** record accessions, reference versions, parameters, checksums and code commits.
5. **Safe interpretation:** keep computational prioritization separate from experimental and clinical claims.

## Current pipeline

```text
Public/local viral FASTA
        ↓
QC and multiple-sequence alignment
        ↓
SpCas9-compatible site scan on the reference
        ↓
Exact site coverage across aligned strains
        ↓
GFF3 overlap annotation
        ↓
Small-host demo screen or export to a genome-scale off-target tool
        ↓
Candidate table + run manifest
        ↓
Optional idealized two-cut sequence-deletion report
```

## Components

### Data acquisition and QC

The HSV-2 runners can download public NCBI records. QC removes exact duplicates, unexpected lengths and records with excessive ambiguous bases. Dataset selection remains a research-design decision, not merely a software default.

### Alignment

MAFFT is the default external aligner. HSV repeats and genome isomers require careful inspection and future normalization; whole-genome alignment should not be accepted blindly.

### Candidate scanner

The current implementation scans both strands for 20-nt SpCas9 protospacers adjacent to an `NGG` PAM. Candidate coordinates are reported relative to the ungapped reference sequence.

### Annotation

GFF3 features are mapped to candidate coordinates. Annotation indicates where a site lies; it does not prove essentiality, accessibility or therapeutic value.

### Host off-target layer

A guarded exhaustive matcher is bundled for small teaching FASTA files. GRCh38-scale analysis is exported to Cas-OFFinder or another validated engine.

### Sequence-disruption simulation

The simulator calculates canonical cut boundaries and an idealized deletion interval between two sites. It reports feature overlap and exact pair coverage across strains. It is not a repair, viability or cure simulator.
