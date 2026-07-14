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
