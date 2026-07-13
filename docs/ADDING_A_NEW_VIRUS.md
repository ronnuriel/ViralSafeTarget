# Adding a new virus

1. Choose and freeze a reference accession and a defensible public/local genome set.
2. Apply documented QC and retain accession-level acceptance/rejection reasons.
3. Produce an aligned multi-FASTA. Every record must have the same alignment length.
4. Supply a reference-matched GFF3 if annotation is available.
5. Copy `configs/research_v0.3.yaml`, adjust only documented thresholds, and retain it
   with the run manifest.
6. Run `vst scan`, inspect the alignment and output funnel, then build a full-host
   off-target input with an assembly identifier.
7. Create a held-out benchmark template rather than tuning against known positives.

FASTA contains nucleotide sequences. GFF3 contains coordinate annotations and
normally does not contain the genome sequence. Their reference identifiers and
coordinate systems must match. Annotation indicates location; it does not establish
essentiality, accessibility, or effectiveness.
