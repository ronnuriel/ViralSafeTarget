# HSV-2 genome-wide computational discovery

Version 0.5 extends the earlier UL19/UL30 pilot to every annotated HSV-2 gene with at
least one eligible pre-human candidate. UL19 and UL30 remain visible benchmarks but
receive no selection or scoring advantage.

## Run the workflow

Activate the project environment and run:

```bash
bash scripts/run_hsv2_genome_wide.sh
```

The script runs the external-tool doctor and delegates all analysis to the public CLI:

```bash
vst discover genome-wide \
  --virus hsv2 \
  --config configs/hsv2_genome_wide.yaml \
  --top-per-gene 50 \
  --global-top 500 \
  --run-cas-offinder \
  --out-dir reports/hsv2_genome_wide
```

Use `--analysis-only` to regenerate tables and the report from validated completed
batches without launching an external tool. Missing, failed, and pending batches remain
explicitly incomplete. Use `--exhaustive --confirm-exhaustive` only when the workload of
screening every retained pre-human candidate is intentional.

For a checksum-resumable exhaustive case-study run in a separate output directory:

```bash
CAS_OFFINDER_PATH=/path/to/cas-offinder bash scripts/run_hsv2_exhaustive.sh
```

The exhaustive mode is intended for sensitivity analysis against the balanced design.
It does not make a zero-hit prediction a safety result.

## Balanced discovery design

The default pre-human panel is the union of up to 50 candidates per annotated gene and
the top 500 candidates globally. Stable candidate IDs deduplicate this union. Selection
does not consult human-screen outcomes. A normalized mapping table preserves every
candidate-feature overlap, including overlapping genes, incomplete annotation, and
intergenic candidates.

Cas-OFFinder batches unique guide sequences to avoid redundant searches. Each batch
manifest maps a query back to every stable coordinate-level candidate ID. A batch is
reused only when its status is completed, its input and candidate-manifest checksums
match, and its raw output exists.

## Interpret the rankings

`gene_rankings.csv` separates computational targetability from biological evidence.
Targetability uses normalized ranking views rather than raw candidate count: best
candidate, robust top-five performance, Wilson-lower-bound clean fraction, predicted
human-hit burden, and conserved-candidate fraction. A support factor and confidence
label prevent a gene with very little completed screening from looking certain.

Biological evidence is optional. When no `data/curated/gene_evidence.tsv` exists, the
workflow reports “biological evidence not supplied”; missing evidence is not converted
to zero and is not a negative score. The import header is available at
`data/curated/gene_evidence.schema.tsv`.

`gene_rank_stability.csv` recomputes nested K=10, K=25, and K=50 views from the single
completed maximum-quota screen. It does not rerun Cas-OFFinder.

## Independent population validation

Population validation remains separate from targetability and host off-target scores.
The preparation workflow can explicitly exclude the discovery FASTA, audits valid IUPAC
ambiguity, and rejects records above a declared ambiguity threshold. Because public viral
records are frequently partial, the optional reference-aware workflow maps each record to
the reference before treating absence of an exact guide/PAM as an observable locus
difference. Install that optional mapper with `pip install -e '.[population]'`.

After both balanced and exhaustive screens finish, `scripts/compare_discovery_modes.py`
produces candidate- and gene-rank sensitivity tables. Rank changes are interpreted as
selection-mode sensitivity, not as biological efficacy or safety evidence.

## Main outputs

- `candidate_feature_map.csv`: one-to-many inclusive coordinate mapping.
- `genome_wide_screening_panel.csv`: balanced or exhaustive panel and selection reasons.
- `genome_wide_human_hits.csv`: every parsed, candidate-expanded predicted hit.
- `genome_wide_candidates_post_human.csv`: completed and explicitly incomplete rows.
- `gene_rankings.csv` and `gene_rank_stability.csv`: gene views and sensitivity.
- `deep_screening_panel.csv`: bounded inputs for optional external tools.
- `pair_hypotheses_*.csv`: bounded theoretical same-gene and multi-target hypotheses.
- `report.html`: the 25-section researcher report.
- `provenance.json`, `combined_batch_manifest.json`, and `stage_timings.json`: audit data.

## Limits

A completed zero-hit result means no predicted hit was enumerated within the declared
GRCh38.p14, SpCas9, PAM, and mismatch model. It is not proof of safety. These outputs are
computational candidates for expert review, not experimental results, clinical guidance,
or evidence of viral inactivation or cure.
