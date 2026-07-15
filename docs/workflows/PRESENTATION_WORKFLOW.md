# Presentation workflow

The showcase turns completed standardized outputs into a reproducible research
brief. It does not rerun the expensive host search and does not generate a wet-lab
protocol.

## Build

```bash
bash scripts/build_hsv2_showcase.sh
```

Or run the two commands separately:

```bash
vst profiles validate \
  --virus-profile configs/viruses/hsv2.yaml \
  --host-profile configs/hosts/human_grch38.yaml \
  --nuclease-profile configs/nucleases/spcas9.yaml \
  --require-host-reference

vst showcase build \
  --virus-profile configs/viruses/hsv2.yaml \
  --host-profile configs/hosts/human_grch38.yaml \
  --nuclease-profile configs/nucleases/spcas9.yaml \
  --out-dir reports/hsv2_showcase
```

## Outputs

- `candidates_evidence_aware.csv`: candidate-level multi-objective table.
- `deep_screening_panel.csv`: fixed-quota shortlist across the configured genes.
- `comparison_set_members.csv`: membership and rationale for each comparison set.
- `strategy_comparison.csv`: set-level summaries without an efficacy claim.
- `research_findings.csv`: auditable computational observations, research
  hypotheses, robustness findings and evidence gaps, with a limitation beside
  every row.
- `population_validation_candidates.csv` and `population_validation_genes.csv`:
  held-out exact-target evidence with locus-aware denominators.
- `figures/`: presentation-ready PNG figures.
- `FINDINGS.md`, `METHODS.md`, `LIMITATIONS.md`: manuscript-oriented text.
- `FINAL_REPORT.html`: self-contained navigation and linked figures/tables.
- `run_manifest.json`: profiles, checksums, parameters and output counts.

The central design choice is to keep five questions separate:

1. How targetable is the sequence under the declared model?
2. Is there cited direct evidence that the gene is essential?
3. How disruptive are the deterministic sequence-level outcome hypotheses?
4. How complete is the evidence available for interpretation?
5. Does the exact target remain present in independent viral records where that
   reference locus is observable?

No combined “therapeutic score” is produced. A high score on one axis does not
compensate for missing evidence on another.

## Suggested presentation order

1. Research problem and scope boundary.
2. Candidate funnel.
3. Genome-wide targetability result.
4. Protein mapping and disruption analysis.
5. Direct HSV-2 evidence versus HSV-1 ortholog evidence.
6. Multi-objective landscape and balanced deep panel.
7. Comparison-set trade-offs.
8. Potentially novel computational observations and their explicit boundaries.
9. Discovery-excluded population validation and balanced-versus-exhaustive sensitivity.
10. Limitations, validation plan and generic-platform roadmap.

Open `notebooks/advanced/genome-wide/11_HSV2_RESEARCH_SHOWCASE_EN.ipynb` for a guided English
walkthrough of the generated artifacts.
Open `notebooks/advanced/genome-wide/12_HSV2_HELDOUT_POPULATION_VALIDATION_EN.ipynb` for the independent
population-panel audit.
