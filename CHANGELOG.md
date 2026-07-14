# Changelog

## 0.8.0 — 2026-07-14

- Added the generic Evidence Agent with alias-aware PubMed, Europe PMC, UniProt, and
  NCBI source discovery; every proposal requires explicit human review before use.
- Published compact, non-clinical HSV-2 result snapshots for the exhaustive host
  screen, balanced showcase, and review-pending evidence run.
- Added clean-clone notebook fallbacks and renamed the Hebrew notebook files with
  language-neutral English identifiers.
- Reorganized documentation by researcher goal, added complete script/notebook/config
  indexes, and added regression tests for repository links and structure.
- Made `make clean` preserve checked-in public result snapshots.

## 0.7.0 — 2026-07-14

- Added the canonical `vst project init/validate/run/resume/status` workflow for a
  self-contained new-virus project.
- Added checksum-aware stage reuse and explicit `external_required` state so a
  missing host-screen result is never treated as zero predicted hits.
- Added `vst reproduce hsv2` plan and guarded `--execute` mode for the complete
  public-data case study.
- Froze the HSV-2 discovery accession cohort and added automatic acquisition and
  validation of the declared HSV-1 ortholog reference.
- Added dedicated reproducibility and new-virus guides plus an end-to-end synthetic
  regression test for the researcher-facing interface.

## 0.6.0 — 2026-07-14

- Added configuration-driven virus, host, and nuclease profiles plus generic schemas.
- Added coding/protein-coordinate mapping, InterPro and disorder context, cross-strain
  entropy, amino-acid conservation, ortholog comparison, and bounded dN/dS summaries.
- Added deterministic -10 to +10 bp indel and theoretical paired-deletion consequence
  tables without claiming repair frequencies or experimental outcomes.
- Kept targetability, direct HSV-2 essentiality, HSV-1 ortholog evidence, predicted
  disruption, evidence coverage, and held-out population support as separate axes.
- Added a discovery-excluded HSV-2 population panel with IUPAC-aware QC and
  reference-aware locus denominators for partial public records.
- Added presentation reports, figures, English notebooks 10–12, and reusable scripts.

## 0.5.0 — 2026-07-14

- Added balanced and guarded exhaustive genome-wide discovery modes.
- Added one-to-many feature mapping, quota sensitivity at K=10/25/50, gene rankings,
  resumable Cas-OFFinder batches, and partial-state reporting.

## 0.4.0 — 2026-07-13

- Added multi-tool adapters, normalized imports, consensus reporting, the public Python
  run-loading API, and separate measured-result imports.

## 0.3.0 — 2026-07-13

- Added stable content-derived candidate IDs, viral occurrence counts, and explicit
  duplicate-guide metadata.
- Added versioned editor, ranking, pair-selection, off-target, and report settings.
- Replaced the tied pre-human score with visible conservation, uniqueness, GC,
  complexity, annotation, evidence, and penalty components.
- Added curated evidence and known-target benchmark schemas without unverified rows.
- Removed genomic-order pair truncation and separated deletion from multi-target
  hypotheses.
- Added Cas-OFFinder build/manifest/summary commands and separate post-human ranking.
- Added resumable HSV-2 workflows, environment doctor, richer provenance, synthetic
  end-to-end verification, and researcher documentation.

## 0.2.0 — 2026-07-13

- Added GitHub-ready English and Hebrew documentation.
- Added a local Streamlit upload interface.
- Added sequence-level single-site and two-cut disruption simulations.
- Added run provenance and SHA-256 manifests.
- Added CLI entry points, CI, contribution templates and Docker support.
- Improved deterministic HSV-2 pilot sampling.

## 0.1.0

- Initial educational Jupyter pipeline for conserved SpCas9 site scanning.
