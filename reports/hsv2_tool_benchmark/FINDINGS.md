# Frozen-panel benchmark findings

## Computational observations

- The identity-frozen benchmark contains 257 unique guides.
- Completed primary metrics: viral_safe_target_pre_human, viral_safe_target_post_human, cas-offinder, crispritz.
- Incomplete or export-required tools: crispor, chopchop, guidescan2.
- Missing results remain missing; no unavailable output was interpreted as zero.
- `cas-offinder` versus `crispritz`: 0.880413.
- `cas-offinder` versus `viral_safe_target_post_human`: 0.289675.
- `cas-offinder` versus `viral_safe_target_pre_human`: 0.085008.
- `crispritz` versus `viral_safe_target_post_human`: 0.296323.
- `crispritz` versus `viral_safe_target_pre_human`: 0.117298.
- `viral_safe_target_post_human` versus `viral_safe_target_pre_human`: 0.961975.
- `cas-offinder` versus `crispritz`, K=10: 9 shared guides (Jaccard 0.818182).
- `cas-offinder` versus `crispritz`, K=25: 24 shared guides (Jaccard 0.923077).
- `cas-offinder` versus `crispritz`, K=50: 49 shared guides (Jaccard 0.960784).
- `cas-offinder` versus `viral_safe_target_post_human`, K=10: 0 shared guides (Jaccard 0.000000).
- `cas-offinder` versus `viral_safe_target_post_human`, K=25: 4 shared guides (Jaccard 0.086957).
- `cas-offinder` versus `viral_safe_target_post_human`, K=50: 10 shared guides (Jaccard 0.111111).
- `cas-offinder` versus `viral_safe_target_pre_human`, K=10: 0 shared guides (Jaccard 0.000000).
- `cas-offinder` versus `viral_safe_target_pre_human`, K=25: 4 shared guides (Jaccard 0.086957).
- `cas-offinder` versus `viral_safe_target_pre_human`, K=50: 10 shared guides (Jaccard 0.111111).
- `crispritz` versus `viral_safe_target_post_human`, K=10: 0 shared guides (Jaccard 0.000000).
- `crispritz` versus `viral_safe_target_post_human`, K=25: 4 shared guides (Jaccard 0.086957).
- `crispritz` versus `viral_safe_target_post_human`, K=50: 10 shared guides (Jaccard 0.111111).
- `crispritz` versus `viral_safe_target_pre_human`, K=10: 0 shared guides (Jaccard 0.000000).
- `crispritz` versus `viral_safe_target_pre_human`, K=25: 4 shared guides (Jaccard 0.086957).
- `crispritz` versus `viral_safe_target_pre_human`, K=50: 10 shared guides (Jaccard 0.111111).
- `viral_safe_target_post_human` versus `viral_safe_target_pre_human`, K=10: 10 shared guides (Jaccard 1.000000).
- `viral_safe_target_post_human` versus `viral_safe_target_pre_human`, K=25: 25 shared guides (Jaccard 1.000000).
- `viral_safe_target_post_human` versus `viral_safe_target_pre_human`, K=50: 50 shared guides (Jaccard 1.000000).
- The largest leave-one-component-out shift was observed for `without_sequence_complexity`: median absolute shift 75.0, maximum 230.0.

## Research hypothesis

The compared systems expose complementary rather than interchangeable axes. A virus-first workflow may add value by retaining viral-population support, host-risk results, gene/protein context, evidence provenance, and multiplex escape analysis in one auditable record.

## Evidence gaps

- CRISPOR, CHOPCHOP, and GuideScan2 require committed raw exports before quantitative rank comparison.
- A completed second-virus benchmark is still required for generalization.
- Independent experimental ground truth is required to test predictive superiority.

## Limitations

- Tool raw scores are not on a common biological scale and were not averaged.
- Runtime is not compared across different guide counts, assemblies, hardware, or search models.
- Ablation is conditioned on an enriched deep-screening panel and is not model retraining.
- Capability evidence describes documented scope, not executable performance.
- No editing, safety, efficacy, viral inactivation, treatment, or cure claim is made.
