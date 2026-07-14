# Long project summary

ViralSafeTarget addresses a practical gap between finding a technically attractive
guide and selecting a biologically meaningful viral target. The framework scans
multiple viral genomes, checks exact target conservation, maps candidates to genes and
coding coordinates, imports or runs model-bounded host searches, separates individual
guide rank from gene-level targetability, and attaches source-provenanced biological
evidence only after human review.

The new virtual analysis maps each cut to every overlapping CDS, enumerates a bounded
indel-size grid without assigning repair probabilities, reports available domain and
disorder context, measures discovery and held-out target support, and computes a
minimum exact-target sequence-change barrier for configured multiplex panels.

The exhaustive HSV-2 source contains 23,108 eligible candidate-coordinate rows and
21,654 unique guide sequences. The publication-facing 257-guide panel produced 271
guide-to-CDS mappings, 5,691 indel hypotheses, and 17,733 single-nucleotide
counterfactuals. Four configured three-guide strategies each required three distinct
substitutions to remove all exact targets under the declared model.

We seek independent code and methods review, evidence curation, cross-tool comparison,
second-virus validation, and experimental collaboration. The project supplies no
wet-lab protocol and does not claim editing, safety, efficacy, viral inactivation,
treatment, or cure.
