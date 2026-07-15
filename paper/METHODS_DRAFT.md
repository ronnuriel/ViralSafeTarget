# Methods draft

## Project and provenance

All analyses were executed through versioned ViralSafeTarget project profiles. Inputs,
parameters, SHA-256 checksums, output paths, and source-table assertions are recorded in
machine-readable manifests. Missing annotations or evidence remain unknown.

## HSV-2 discovery source

The committed exhaustive source contains 28,578 initial candidates, 23,108 eligible
candidate-coordinate rows, and 21,654 unique guide sequences. Cas-OFFinder 2.4.1 was
completed in 109 batches against GRCh38.p14 under the configured mismatch/PAM model.
The source table retains 440,341 predicted human-match rows; 2,668 candidate-coordinate
rows have zero predicted matches under that model. Zero predicted matches are not a
safety conclusion.

## Generic CDS mapping

Reference-matched FASTA and GFF CDS records were converted to coding-orientation
coordinate models. GFF phase and strand were respected. Every CDS containing a guide
cut boundary was retained, allowing one guide to map to overlapping genes. Protein
coordinates were calculated from the number of coding bases preceding the cut.

## Virtual knockout

For every guide-to-CDS mapping, integer indel sizes from -10 through +10 bp were
enumerated. Deletions were anchored downstream in coding orientation. Insertions were
represented by size only and their nucleotide sequence remained unspecified.
Frameshifts were classified by size modulo three. Premature stops were reported only
when determinable from a deletion sequence; insertion-dependent stop and retained-
protein fields remained unknown. Domain, disorder, and conserved-region overlaps were
reported only when compatible optional annotations existed.

Grid fractions describe equally weighted size hypotheses and are not repair
probabilities.

## Escape robustness

Discovery and held-out exact protospacer/PAM coverage were reported separately. All
three alternate bases at every protospacer and PAM position were enumerated as
single-nucleotide counterfactuals. Protospacer substitutions remove an exact match;
PAM substitutions were classified using the configured IUPAC PAM pattern.

For each configured multiplex panel, target-disrupting substitutions were represented
as sets of affected guides. An exact set-cover calculation returned the minimum number
of distinct substitutions whose union removed all exact targets. This is a sequence-
level barrier, not an evolutionary probability.

## Strategy comparison

Four three-guide strategies were configured outside core code: top-ranking-only,
essential/replication-focused, targetability-focused, and mechanism-diverse. Host-match
counts, observed population support, bounded disruption summaries, biological context,
and escape barrier remained separate columns. No combined therapeutic score was
computed.

## Systematic multi-tool benchmark

We froze the 257 unique guides in the HSV-2 deep-screening panel before comparison.
ViralSafeTarget pre-host and post-host ranks and completed Cas-OFFinder hit counts were
normalized as within-tool ranks. CRISPRitz profiles were parsed only when a raw profile
was available. The committed CRISPRitz 2.6.6 run used the official Docker image
(`sha256:8a6c8212621ee6cc467e7a5bc7ff4405cbb83bfa7c0487c232977f4a09ff0273`),
GRCh38.p14 (`GCF_000001405.40`), an NGG PAM, up to three mismatches, eight threads,
and no bulge or population-variant options. Its 257-row profile, target list, extended
profile, command, runtime, and hashes were retained. CRISPOR, CHOPCHOP, and GuideScan2
were treated as documented-export
integrations and remained export-required when no raw export was committed. Missing
output was retained as missing and was never interpreted as zero predicted risk.

We calculated candidate coverage, pairwise Spearman rank agreement for non-constant
shared rankings, and top-10, top-25, and top-50 overlap. Raw scores from different
tools were not averaged. Runtime was considered comparable only when hardware,
assembly, guide count, and search parameters were equivalent.

For sensitivity analysis, each configured ViralSafeTarget score component was omitted
in turn and the transparent weighted score was recomputed on the frozen panel. This
leave-one-component-out procedure was not treated as model retraining or biological
validation.

## Reproducibility and researcher usability

The 0.10.0 source distribution and wheel were built from the tagged package metadata.
A temporary environment outside the repository installed the wheel non-editably and
executed the public CLI (`doctor`, `quickstart`, `plan`, `run`, `status`, `resume`,
`open --no-browser`, and `export`). Each project persisted stage signatures, explicit
external-required states, elapsed timings, provenance, and a portable result archive.

As a software generalization check, the same installed CLI retrieved BK polyomavirus
NC_001538.1 through NCBI E-utilities, derived a reference-matched GFF3 from the GenBank
record, and ran the generic workflow. No BK-specific gene name or rule was introduced
in core code. Because this check used a reference-only strain panel and no host
assembly, population conservation and host-risk outputs were explicitly unavailable.
