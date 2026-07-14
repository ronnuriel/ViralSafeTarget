# Existing tools and how ViralSafeTarget relates to them

## MAFFT

Multiple-sequence alignment. ViralSafeTarget expects an aligned multi-FASTA and can call MAFFT in the HSV-2 pilot workflow.

## Cas-OFFinder

Fast genome-wide enumeration of potential CRISPR off-target sites for configurable PAMs and mismatch limits. ViralSafeTarget exports candidate guides for this stage.

## CRISPRitz

Off-target search and analysis that can include mismatches, bulges and human genetic variation. It is a stronger option when population variants matter.

## GuideScan2, CRISPOR and CHOPCHOP

General-purpose guide-design and specificity tools. They are useful baselines for on-target/off-target scoring once a target region has already been selected.

## CRISPResso2

Analysis of actual deep-sequencing data from genome-editing experiments. It measures observed editing outcomes; it is not a pre-experiment simulator.

## OffRisk

Functional annotation and risk labeling of predicted human off-target sites.

## What ViralSafeTarget adds

The project begins with **many viral strains**, not one already selected locus. It connects:

- viral sequence quality control;
- strain conservation;
- editor/PAM compatibility;
- reference annotation;
- host off-target export;
- pair coverage;
- explicit rejection reasons;
- provenance and reproducibility.

Its research value must be demonstrated by benchmarks against the mature tools above.
