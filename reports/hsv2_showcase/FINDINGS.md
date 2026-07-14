# Findings

## Research question

Can a virus-first pipeline distinguish sequence targetability from biological target evidence and produce an auditable, mechanistically balanced shortlist?

## Data funnel

- 28,578 initial candidate coordinates.
- 23,108 passed pre-human filters.
- 2,952 entered the balanced host screen.
- 595 candidates have zero predicted human hits under the declared model.
- 88 top candidates were mapped into the protein-disruption analysis; 65 of those have zero predicted human hits.
- 36 candidates form the balanced deep panel.

## Main finding

The highest sequence-targetability genes are not automatically the best-supported biological targets. The project therefore keeps targetability, direct essentiality evidence, predicted protein disruption, and evidence coverage as separate axes. No combined therapeutic score is reported.

The current curated set contains 0 genes with a direct HSV-2 essentiality score. Direct HSV-2 knockdown phenotypes remain phenotype evidence rather than null-essentiality claims. HSV-1 ortholog evidence is displayed separately.

## Potentially novel computational observations

The analysis produced 7 auditable observations, hypotheses, robustness results, or evidence gaps in `research_findings.csv`. These rows are intended to help researchers decide what deserves independent investigation. They are not novelty claims against the complete literature and are not evidence of editing or treatment efficacy.

- **genome_wide_reprioritization (computational_observation):** UL3 ranks first for sequence targetability; UL30 ranks 26. Potential value: A virus-wide search can surface technically tractable genes that a benchmark-only analysis would miss. Limitation: The rank measures model-bounded targetability, not gene importance, editing, or viral inhibition; exhaustive sensitivity analysis is pending.
- **targetability_evidence_divergence (evidence_gap):** UL3 is technically highly targetable, while direct HSV-2 essentiality remains unknown and the curated HSV-1 ortholog evidence reports nonessentiality in the tested cell-culture context. Potential value: UL3 is a useful example of why guide quality and target biology must remain separate research questions. Limitation: HSV-1 context cannot establish HSV-2 function, and nonessentiality in one culture system does not exclude other phenotypes.
- **multi_axis_convergence (research_hypothesis):** UL52, UL30, UL18 each combine at least one population-supported exact target with source-linked HSV-1 ortholog evidence and a mapped protein-disruption model. Potential value: These genes form a defensible evidence-aware shortlist for independent mechanistic assessment without claiming a therapeutic rank. Limitation: No direct HSV-2 null-essentiality result is present for these genes, and population records are mostly partial.
- **candidate_specific_population_variation (computational_observation):** The leading coordinate-level candidate VST-2e9f052157f9bf29 (UL36) retained the exact target in 112/113 observable held-out loci. Potential value: Population support should be evaluated per coordinate rather than inferred from a gene-level conservation label. Limitation: Exact sequence retention is not evidence of editor activity, host safety, delivery, or efficacy; unresolved partial records remain.
- **conserved_disruption_signal (research_hypothesis):** UL18 has the highest predicted protein-disruption score (0.757) among the nine mapped genes and mean amino-acid conservation 1.000 in the discovery alignment. Potential value: A conserved protein with cuts mapping to disruptive sequence contexts is a useful hypothesis for deeper functional prioritization. Limitation: Size-enumerated outcomes are not repair-frequency predictions, and the 14-genome discovery alignment is not population representative.
- **quota_stable_gene_set (robustness_observation):** UL10, UL18, UL20, UL3, UL36 remain in the top 10 at K=10, 25, and 50. Potential value: Quota-stable genes are less likely to be artifacts of one arbitrary per-gene sampling depth. Limitation: Stability is internal to the current scoring model and balanced panel; it is not biological validation.
- **direct_hsv2_essentiality_gap (evidence_gap):** Direct, scored HSV-2 essentiality evidence is available for 0/9 deeply analyzed genes. Potential value: The missing-evidence map identifies where literature curation or new functional datasets would most improve prioritization. Limitation: Absence of curated evidence is not evidence of nonessentiality and may reflect incomplete literature coverage.

## External population-genomics context

The convenience strain alignment is checked against 4 source-linked population or gene-variability findings. These studies support generally low HSV-2 divergence but identify lineage- and locus-specific exceptions, including UL30 and UL53. External studies contextualize the analysis; they do not replace held-out sequence validation.

In a separate held-out panel, 2,952 candidates had an observable reference locus and 1,764 retained the exact guide/PAM in every observable record. Discovery genomes were excluded. Partial records remain unresolved where the locus is not aligned, and these results do not alter targetability scores.

## Presentation claim

ViralSafeTarget demonstrates a reproducible framework for discovering and comparing viral target hypotheses. It does not demonstrate editing, viral eradication, delivery, safety, efficacy, or a cure.
