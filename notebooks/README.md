# Notebook guide

Notebooks explain and reproduce workflows; the `vst` CLI remains the canonical
automation interface. Run Jupyter from the repository root so paths resolve reliably:

```bash
jupyter lab
```

## Hebrew learning path

| Order | Notebook | Purpose |
|---:|---|---|
| 0 | [`00_START_HERE.ipynb`](00_START_HERE.ipynb) | Concepts, project map, and learning order |
| 1 | [`01_DEMO_PIPELINE.ipynb`](01_DEMO_PIPELINE.ipynb) | Small synthetic end-to-end demo |
| 2 | [`02_REAL_DATA_SETUP.ipynb`](02_REAL_DATA_SETUP.ipynb) | Public data, references, and alignment setup |
| 3 | [`03_RESEARCH_PROTOCOL.ipynb`](03_RESEARCH_PROTOCOL.ipynb) | Computational research question and benchmark design |
| 4 | [`04_CAS_OFFINDER.ipynb`](04_CAS_OFFINDER.ipynb) | Host-screen concepts and output interpretation |
| 5 | [`05_SEQUENCE_DISRUPTION_SIMULATION.ipynb`](05_SEQUENCE_DISRUPTION_SIMULATION.ipynb) | Bounded sequence/deletion hypotheses |

The filenames are English and language-neutral; the notebook content remains Hebrew.

## English case-study and research notebooks

| Notebook | Purpose | Expected inputs |
|---|---|---|
| [`07_HSV2_MULTITOOL_COMPARISON_EN.ipynb`](07_HSV2_MULTITOOL_COMPARISON_EN.ipynb) | Compare/import external tool results | Pilot outputs or fixtures |
| [`08_RUN_FULL_PIPELINE_EN.ipynb`](08_RUN_FULL_PIPELINE_EN.ipynb) | End-to-end HSV-2 orchestration | Public/cached inputs; fixture fallback |
| [`09_HSV2_GENOME_WIDE_DISCOVERY_EN.ipynb`](09_HSV2_GENOME_WIDE_DISCOVERY_EN.ipynb) | Genome-wide discovery and stability | Completed or synthetic discovery outputs |
| [`10_HSV2_GENE_FUNCTION_AND_DISRUPTION_EN.ipynb`](10_HSV2_GENE_FUNCTION_AND_DISRUPTION_EN.ipynb) | Protein mapping and disruption hypotheses | Genome-wide outputs and annotations |
| [`11_HSV2_RESEARCH_SHOWCASE_EN.ipynb`](11_HSV2_RESEARCH_SHOWCASE_EN.ipynb) | Presentation-ready multi-axis analysis | Showcase outputs |
| [`12_HSV2_HELDOUT_POPULATION_VALIDATION_EN.ipynb`](12_HSV2_HELDOUT_POPULATION_VALIDATION_EN.ipynb) | Held-out sequence support | Population-validation outputs |
| [`13_EVIDENCE_AGENT_HUMAN_REVIEW_EN.ipynb`](13_EVIDENCE_AGENT_HUMAN_REVIEW_EN.ipynb) | Literature proposal and human-review workflow | GFF or checked-in public gene catalog |

## Reproducibility rules

- Network access is off by default in tested notebooks.
- Missing real inputs use an explicit fixture or public-snapshot path when supported.
- A missing external result remains pending; it is never converted to a favorable score.
- Notebook conclusions must preserve the non-clinical scope described in
  [`../DISCLAIMER.md`](../DISCLAIMER.md).
