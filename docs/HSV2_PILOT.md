# HSV-2 UL19/UL30 computational pilot

## Purpose

The pilot produces a reproducible, focused set of unique SpCas9 candidates in UL19
and UL30. These genes are used as named pilot strata, not as fabricated claims of
essentiality or therapeutic validation. The shipped HSV-2 evidence table is empty
until repository sources are curated at the exact gene/reference level.

## Prepare or resume the real-data stages

```bash
bash scripts/doctor.sh
bash scripts/run_real_hsv2.sh --sample-size 25
```

The real-data runner validates cached ZIP files, retries interrupted NCBI downloads,
and uses checksum stamps for preparation, MAFFT, and candidate generation. It resumes
valid stages and falls back to `python -m viral_safe_target` when no console script is
installed.

To obtain GRCh38 locally and build the full pilot input:

```bash
bash scripts/run_real_hsv2.sh --with-human --sample-size 25
bash scripts/run_hsv2_pilot.sh
```

Alternatively set `HUMAN_FASTA_DIRECTORY` to an existing Cas-OFFinder-compatible
GRCh38 FASTA directory. The focused script applies the YAML conservation, GC,
complexity, uniqueness, and CDS filters; deterministically selects at most 100
candidates per gene; and writes:

```text
reports/hsv2_pilot/cas_offinder_input.txt
reports/hsv2_pilot/offtarget_selected_candidates.csv
reports/hsv2_pilot/pair_hypotheses_same_gene.csv
reports/hsv2_pilot/pair_hypotheses_multi_target.csv
```

Run Cas-OFFinder externally:

```bash
cas-offinder reports/hsv2_pilot/cas_offinder_input.txt C \
  reports/hsv2_pilot/cas_offinder_output.tsv
bash scripts/run_hsv2_pilot.sh
```

On the second invocation, the script detects the result and generates
`candidates_ranked_post_human.csv`, `predicted_human_hits.csv`, and `report.html`.
No predicted hit within three mismatches in the configured GRCh38 assembly is a
search result, not proof of safety.

## Alignment warning

HSV genomes contain repeats and alternative isomer orientations. Whole-genome MAFFT
output requires expert inspection and, where appropriate, normalization before
conservation estimates are treated as research findings.
