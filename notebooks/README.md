# Notebook guide

External researchers should start with the single canonical English notebook:

- [`00_VIRALSAFETARGET_END_TO_END_EN.ipynb`](00_VIRALSAFETARGET_END_TO_END_EN.ipynb)

It calls the installed `vst` CLI and supports `demo`, `hsv2_snapshot`, and
`custom_project` modes. Expensive external and network operations are opt-in.

## Advanced scientific notebooks

- [`advanced/genome-wide/`](advanced/genome-wide/) — discovery, full orchestration,
  showcase, and held-out population support.
- [`advanced/disruption-escape/`](advanced/disruption-escape/) — protein mapping,
  virtual-knockout hypotheses, and exact-target escape counterfactuals.
- [`advanced/evidence/`](advanced/evidence/) — source-linked literature proposals with
  mandatory human review.
- [`advanced/benchmarking/`](advanced/benchmarking/) — external-tool comparison and
  systematic benchmark analyses.

## Learning archive

The earlier Hebrew teaching notebooks are preserved under [`learning/he/`](learning/he/).
They are not part of the publication-facing scientific workflow.

## Reproducibility rules

- Missing external output remains pending, never zero.
- Unknown biological evidence remains unknown.
- Virtual outcomes are sequence hypotheses, not biological repair probabilities.
- Notebook outputs are cleared in Git; CI executes the canonical notebook from a clean wheel.
