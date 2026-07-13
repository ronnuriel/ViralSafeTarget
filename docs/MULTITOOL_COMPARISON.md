# Multi-tool comparison

Raw efficiency scores, hit counts and sequence-prioritization scores have different meanings and
must not be averaged. ViralSafeTarget first ranks each documented metric within its own tool and
converts that rank to a desirability percentile. It then supports weighted Borda/percentile
aggregation, median rank and worst-case rank.

Every consensus row reports tool coverage, missing tools, per-tool ranks, rank variance,
disagreement and an explanation. Coverage scales the weighted consensus so a result supported by
one tool receives a low-coverage warning. Missing remains `NaN`; it never becomes a perfect score.

Agreement outputs include Spearman and Kendall rank association for sufficiently shared candidate
sets, plus top-5/top-10/top-20 overlap and Jaccard overlap. Agreement is descriptive and does not
show that the most popular ranking is biologically correct.
