# Benchmarking

Copy `benchmarks/known_targets.template.csv` and populate only sequences verified
against primary sources and the stated reference. Do not infer a positive label from
the fact that a gene or intervention appeared in a paper.

Run:

```bash
vst benchmark \
  --candidates reports/run/candidates_ranked_pre_human.csv \
  --known-targets benchmarks/known_targets.curated.csv \
  --out-dir reports/run/benchmark
```

The command reports regeneration, pre- and post-human ranks, percentile, exclusion
reason, recovery rate, top-k recall, and gene-level recovery. Empty templates are
valid and preferable to fabricated evidence. Freeze the benchmark before score
tuning and report negative as well as positive examples.
