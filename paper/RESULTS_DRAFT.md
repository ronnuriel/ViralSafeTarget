# Results draft

## Exhaustive source validation

Six assertions were checked against committed source tables before the new snapshot was
written: 28,578 initial candidates, 23,108 eligible rows, 21,654 unique guide
sequences, 440,341 predicted human-match rows, 2,668 zero-hit rows under the configured
model, and 109 completed batches. All assertions passed.

## Deep-panel virtual analysis

The configured deep panel contained 257 unique guides. Generic mapping produced 271
guide-to-CDS rows because overlapping annotations were preserved; 250 guides had at
least one CDS mapping. The -10 through +10 bp grid produced 5,691 sequence hypotheses.
Optional protein annotations were not imputed when unavailable.

## Observed variation and counterfactuals

The analysis enumerated 17,733 single-nucleotide protospacer/PAM counterfactuals. Every
SpCas9 NGG target contributed 69 substitutions, of which 66 removed the exact target:
60 protospacer substitutions and six substitutions at the two constrained PAM
positions. Held-out exact-target coverage was available for 200 of 257 guides; 57
guides remained unknown because the committed held-out snapshot did not cover them.

## Configured strategy comparison

Each of the four configured three-guide panels required three distinct substitutions
to remove all exact targets under the exact set-cover model. All panel members had
discovery exact-target coverage of 1.0 and zero predicted host matches under the
configured source model. Minimum held-out coverage was 0.973913 for the top-ranking,
targetability-focused, and mechanism-diverse panels, and 0.982609 for the
essential/replication-focused panel. These marginal guide-level values do not establish
joint panel coverage or safety.

## Gene-context observations

UL3 contributed 20 analyzed guides, with 10 cutting an annotated UL3 domain; direct
HSV-2 biological evidence remains incomplete. UL52 contributed 11 guides and UL30 20,
with 11 and 16 respectively cutting configured domains, while their replication-focused
rationales remain separate evidence fields. All three UL18 guides cut its annotated
domain and were positioned in the earlier portion of the protein (median relative
position 0.192). UL36 contributed seven guides, including the leading individual source
candidate, while its exhaustive gene-level portfolio rank remained 21.

These observations prioritize questions for review; they do not establish viral
phenotype or therapeutic relevance.

## Multi-tool benchmark

The systematic benchmark froze 257 unique candidate identities before comparison.
Primary metrics were complete for the ViralSafeTarget pre-host and post-host rankings,
Cas-OFFinder, and an independently executed CRISPRitz 2.6.6 reference-genome search.
CRISPOR, CHOPCHOP, and GuideScan2 remained export-required because no raw output was
committed; they were not assigned zero values or inferred ranks.

The two host-search tools showed strong but incomplete agreement: Cas-OFFinder and
CRISPRitz ranks had a Spearman correlation of 0.880413, with 9, 24, and 49 shared
guides in their respective top 10, 25, and 50 lists. In contrast, their primary
off-target-burden ranks had low agreement with the composite ViralSafeTarget ranks
(Spearman 0.085008-0.296323, depending on the comparison). This does not identify a
superior tool: the host-search metric and the composite virus-first score represent
different decision axes. ViralSafeTarget pre-host and post-host ranks remained highly
correlated (Spearman 0.961975), with identical top 10, 25, and 50 membership on this
enriched panel.

The CRISPRitz run used the official version 2.6.6 Docker image against GRCh38.p14,
through three mismatches, without bulges or population variants. It reported all 257
guides in 709.78 wall-clock seconds. Runtime was not compared with the recorded
Cas-OFFinder source runtime because the latter covered 23,108 candidates rather than
the frozen 257-guide panel.

Leave-one-component-out analysis quantified the sensitivity of ViralSafeTarget ranks
within the enriched deep-screening panel. Removing sequence complexity caused the
largest change (median absolute rank shift 75; maximum 230; top-10 overlap 0), whereas
removing GC changed at most nine rank positions and preserved the top 10, 25, and 50.
Conservation, viral uniqueness, annotation, and gene-evidence components produced no
rank change in this selected panel because those columns were constant there. The
ablation therefore diagnoses the panel and scoring configuration; it does not establish
which ranking better predicts editing or viral phenotype.

## Researcher usability and second-virus proof

Version 0.10.0 was installed as a non-editable wheel outside the source checkout and
completed the synthetic quickstart, plan, run, status, resume, report, and export
workflow. The same installed CLI retrieved and analyzed the public BK polyomavirus
reference NC_001538.1 without virus-specific core-code changes. It estimated 543
editor-compatible sites before the configured 500-row discovery cap and completed
annotation mapping, virtual analyses, multiplex comparison, reporting, caching, and
export. Host screening remained explicitly `external_required`, and the reference-only
run was not interpreted as population conservation or a biological case study.
