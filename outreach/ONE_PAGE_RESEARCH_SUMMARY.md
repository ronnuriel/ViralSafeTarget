# ViralSafeTarget HSV-2 computational research summary

## What we built

ViralSafeTarget is a virus-first, reproducible framework that separates guide sequence
quality, gene-level targetability, predicted host matches, bounded coding disruption,
observed viral-population support, biological evidence provenance, and exact-target
multiplex escape robustness.

## What the HSV-2 case study found

- 28,578 initial and 23,108 eligible candidate-coordinate rows.
- 21,654 unique guide sequences.
- 440,341 predicted human-match rows from 109 completed Cas-OFFinder batches.
- 2,668 candidate-coordinate rows with zero predicted matches under the configured
  model; this is not proof of safety.
- UL3, UL10, and UL52 lead the exhaustive gene targetability table, while the leading
  individual guide is in UL36 and UL30 ranks eighth.
- A 257-guide deep panel produced 271 guide-to-CDS mappings, 5,691 bounded indel
  hypotheses, and 17,733 single-nucleotide exact-target counterfactuals.
- Four configured three-guide strategies each had a sequence escape barrier of three
  substitutions under the exact-target model.
- A frozen 257-guide benchmark produced complete ViralSafeTarget, Cas-OFFinder, and
  CRISPRitz metrics. Cas-OFFinder and CRISPRitz host-search ranks correlated at 0.880413
  and shared 49 of their top 50 guides under closely matched reference-genome settings.
- CRISPOR, CHOPCHOP, and GuideScan2 are documented in the capability matrix but remain
  export-required; their quantitative results were not fabricated or treated as zero.
- The installed 0.10.0 wheel completed a reference-only BK polyomavirus workflow without
  virus-specific core-code changes; this demonstrates usability, not biological
  generalization.

## Why it matters

The results demonstrate that the best individual guide, the most targetable gene, and
the gene with the strongest biological rationale can differ. A reproducible framework
can make those differences visible before costly experimental selection.

## What we are asking for

We welcome independent computational review, comparison with other tools, advice on
biological evidence curation, and experimental collaboration using collaborator-
approved protocols. The current outputs are hypotheses, not claims of editing,
safety, efficacy, viral inactivation, treatment, or cure.
