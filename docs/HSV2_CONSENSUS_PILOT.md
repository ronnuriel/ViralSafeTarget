# HSV-2 v0.4 consensus pilot

This English tutorial uses the completed v0.3 run and does not rescan the full candidate set:

```bash
bash scripts/run_hsv2_consensus.sh
```

The workflow selects exactly 32 computational candidates with no predicted human hit through three
mismatches in the configured Cas-OFFinder/SpCas9 model: 23 UL30 and 9 UL19. It builds cached
CRISPRitz inputs, imports any exports placed in `reports/hsv2_consensus/external_exports/`, and writes
a report even when external stages are pending.

Outputs are `tool_results_long.csv`, `candidate_tool_matrix.csv`, `consensus_candidates.csv`,
`tool_coverage.csv`, `model_agreement.csv`, `disagreement_report.csv`,
`unmatched_external_rows.csv`, `report.html`, and `run_manifest.json` under
`reports/hsv2_consensus/`.

The current top position can remain unchanged when only ViralSafeTarget and tied zero-hit
Cas-OFFinder ranks are available, but coverage is explicitly low. Add independent results before
interpreting consensus. These are computational prioritization hypotheses; no candidate is shown to
be safe, effective, validated, or curative.
