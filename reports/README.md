# Public HSV-2 result snapshots

This directory contains a deliberately small, reviewable snapshot of completed
computational runs. Large host-hit tables, reference genomes, tool caches, and raw
batch outputs remain excluded from Git. They can be regenerated with the versioned
workflows and notebooks. Compact raw outputs required to audit the frozen-panel
benchmark are an explicit exception.

The compact [BK polyomavirus second-virus snapshot](second_virus_bk_polyomavirus/README.md)
uses the same installed CLI and schemas without virus-specific core-code changes. It is
a software usability proof, not a therapeutic case study.

> These are computational research outputs. They do not establish editing,
> safety, viral inhibition, latency clearance, treatment efficacy, or a cure.

## Which result is current?

The genome-wide analyses were run at two different sampling depths:

1. `hsv2_showcase/` uses the earlier balanced panel of 2,952 candidate coordinates.
   It is useful for the multi-axis presentation and held-out-population walkthrough.
2. `hsv2_genome_wide_exhaustive/` is the later exhaustive host screen of all 23,108
   eligible candidate coordinates. Use this directory for the current sequence-
   targetability ranks.

The different sampling depths explain rank changes such as UL30 moving from rank 26
in the balanced showcase to rank 8 in the exhaustive screen.

## Exhaustive genome-wide result

- 23,108 eligible candidate coordinates representing 21,654 unique guide sequences.
- 109/109 Cas-OFFinder batches completed against GRCh38.p14.
- 440,341 predicted human-hit rows were retained.
- 2,668 candidate-coordinate rows had no predicted hit under the declared search
  model. This is not proof of safety.
- The five leading genes by exhaustive computational targetability are UL3, UL10,
  UL52, UL47, and UL11.
- UL30 ranks 8, UL20 ranks 11, UL18 ranks 12, UL53 ranks 17, UL19 ranks 20, and UL36
  ranks 21.
- The leading individual candidate remains `VST-2e9f052157f9bf29`, mapped to UL36.
- UL3, UL10, UL36, and UL20 remain top-10 genes at K=10, 25, and 50.

Start with:

- [`hsv2_genome_wide_exhaustive/report.html`](hsv2_genome_wide_exhaustive/report.html)
- [`hsv2_genome_wide_exhaustive/gene_rankings.csv`](hsv2_genome_wide_exhaustive/gene_rankings.csv)
- [`hsv2_genome_wide_exhaustive/top_candidates_global.csv`](hsv2_genome_wide_exhaustive/top_candidates_global.csv)
- [`hsv2_genome_wide_exhaustive/gene_rank_stability.csv`](hsv2_genome_wide_exhaustive/gene_rank_stability.csv)

## Evidence Agent result

The first review-pending evidence run covered UL3, UL10, UL18, UL20, UL36, UL52,
UL53, UL19, and UL30:

- 72 generated queries.
- 86 normalized source records.
- 76 source-linked proposals: 23 possible direct HSV-2 records, 31 HSV-1 ortholog
  records, 21 mixed/other-virus records, and one unresolved record.
- 0 approved proposals and 0 automatic biological-score changes.

Every proposal remains pending until a researcher checks the paper, excerpt,
experimental system, virus species, and interpretation. Missing evidence remains
unknown.

Start with:

- [`hsv2_evidence_agent/evidence_review_report.html`](hsv2_evidence_agent/evidence_review_report.html)
- [`hsv2_evidence_agent/review_queue.tsv`](hsv2_evidence_agent/review_queue.tsv)
- [`hsv2_evidence_agent/evidence_manifest.json`](hsv2_evidence_agent/evidence_manifest.json)

## Presentation snapshot

The balanced-panel presentation materials remain available under
[`hsv2_showcase/`](hsv2_showcase/). They preserve the earlier multi-axis analysis,
explicit limitations, research hypotheses, and held-out-population summaries. They
must not be confused with the later exhaustive targetability rank.

## Virtual knockout and exact-target escape snapshot

[`hsv2_virtual_knockout_escape/`](hsv2_virtual_knockout_escape/) is the current
publication-facing analysis of the exhaustive deep-screening panel:

- 257 unique guides and 271 guide-to-CDS rows, preserving overlapping annotations;
- 5,691 size-defined indel hypotheses across the configured -10 to +10 bp grid;
- 17,733 single-nucleotide protospacer/PAM counterfactuals;
- held-out exact-target coverage available for 200 guides and unknown for 57;
- four configured three-guide strategies with a sequence-level escape barrier;
- six source-table assertions checked before the snapshot was written.

Start with:

- [`hsv2_virtual_knockout_escape/virtual_knockout_escape_report.html`](hsv2_virtual_knockout_escape/virtual_knockout_escape_report.html)
- [`hsv2_virtual_knockout_escape/FINDINGS.md`](hsv2_virtual_knockout_escape/FINDINGS.md)
- [`hsv2_virtual_knockout_escape/strategy_comparison.csv`](hsv2_virtual_knockout_escape/strategy_comparison.csv)
- [`virtual knockout notebook`](../notebooks/advanced/disruption-escape/14_VIRTUAL_KNOCKOUT_AND_ESCAPE_EN.ipynb)

## Systematic multi-tool benchmark

[`hsv2_tool_benchmark/`](hsv2_tool_benchmark/) freezes 257 candidate identities and
separates executable comparison from documented capability evidence. It contains raw
outputs when available, normalized within-tool ranks, explicit missingness, top-K
overlap, ablation, figures, and a checksum manifest. A tool without a committed raw
export remains `pending` or `export_required`.

Current executable coverage:

- complete ViralSafeTarget pre-host and post-host ranks;
- complete Cas-OFFinder source counts and a complete CRISPRitz 2.6.6 profile;
- Cas-OFFinder/CRISPRitz rank correlation 0.880413 and top-50 overlap 49/50;
- CRISPOR, CHOPCHOP, and GuideScan2 retained as `export_required`, not zero;
- largest ablation shift after removing sequence complexity (median 75, maximum 230).

Start with:

- [`hsv2_tool_benchmark/multitool_benchmark_report.html`](hsv2_tool_benchmark/multitool_benchmark_report.html)
- [`hsv2_tool_benchmark/tool_execution_status.csv`](hsv2_tool_benchmark/tool_execution_status.csv)
- [`hsv2_tool_benchmark/rank_agreement.csv`](hsv2_tool_benchmark/rank_agreement.csv)
- [`systematic benchmark notebook`](../notebooks/advanced/benchmarking/15_SYSTEMATIC_MULTITOOL_BENCHMARK_EN.ipynb)

## Reproduction notebooks

- [`genome-wide discovery`](../notebooks/advanced/genome-wide/09_HSV2_GENOME_WIDE_DISCOVERY_EN.ipynb)
- [`gene function and disruption`](../notebooks/advanced/disruption-escape/10_HSV2_GENE_FUNCTION_AND_DISRUPTION_EN.ipynb)
- [`research showcase`](../notebooks/advanced/genome-wide/11_HSV2_RESEARCH_SHOWCASE_EN.ipynb)
- [`held-out population validation`](../notebooks/advanced/genome-wide/12_HSV2_HELDOUT_POPULATION_VALIDATION_EN.ipynb)
- [`evidence review`](../notebooks/advanced/evidence/13_EVIDENCE_AGENT_HUMAN_REVIEW_EN.ipynb)
- [`virtual knockout and escape`](../notebooks/advanced/disruption-escape/14_VIRTUAL_KNOCKOUT_AND_ESCAPE_EN.ipynb)
- [`systematic benchmark`](../notebooks/advanced/benchmarking/15_SYSTEMATIC_MULTITOOL_BENCHMARK_EN.ipynb)

The notebooks and machine-readable tables are the reproducible record. The checked-in
HTML files are convenience views of those outputs.
