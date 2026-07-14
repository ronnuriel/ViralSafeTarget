# Python SDK

```python
import viral_safe_target as vst

run = vst.load_run("reports/hsv2_pilot")
candidates = run.candidates
hits = run.human_hits
pairs = run.same_gene_pairs

baseline = vst.candidate_metrics_as_tool_results(candidates)
comparison = vst.compare_tools(candidates, [baseline])
comparison.consensus_candidates.head()
```

Stable public objects are exported from `viral_safe_target`: `load_run`, `ResearchRun`,
`CandidateTable`, `ToolResultTable`, `ToolAdapter`, `CandidateScorer`, `compare_tools`,
`build_consensus`, `candidate_metrics_as_tool_results`, and `load_external_results`.

Typed wrappers expose their underlying pandas object as `.dataframe`. The normalized tool table
validates all required columns and preserves missing values.
