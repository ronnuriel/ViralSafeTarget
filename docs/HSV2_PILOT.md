# HSV-2 pilot

## Goal

Produce a reproducible **pre-human-screen** list of conserved SpCas9-compatible sites across a modest public HSV-2 genome set.

## Recommended run

```bash
conda env create -f environment.yml
conda activate viral-safe-target
bash scripts/run_real_hsv2.sh --sample-size 25
```

A cross-platform Entrez alternative is available when the NCBI Datasets CLI is inconvenient:

```bash
export NCBI_EMAIL="you@example.org"
python scripts/run_hsv2_pilot.py --email "$NCBI_EMAIL" --max-genomes 25
```

## Stages

1. Download reference accession `NC_001798.2` and public high-length HSV-2 records.
2. Remove exact duplicates, high-`N` records and unexpected lengths.
3. Convert reference GenBank features to GFF3.
4. Align the selected genomes with MAFFT.
5. Scan both strands for `20 nt + NGG` sites.
6. Calculate exact 23-nt site coverage across the alignment.
7. Map candidates to the reference annotation.
8. Generate an idealized two-cut sequence report.

## Critical alignment caveat

HSV genomes contain repeats and alternative isomer orientations. The pilot is a starting point; researchers must inspect and normalize the alignment before treating conservation estimates as authoritative.

## Human off-target stage

Prepare a fixed human assembly and run a validated genome-scale engine. Keep the genome build, tool version and parameters in the report.

```bash
cas-offinder reports/real_hsv2/cas_offinder_input.txt C \
  reports/real_hsv2/cas_offinder_output.tsv

python scripts/summarize_cas_offinder.py \
  --candidates reports/real_hsv2/candidates_pre_human.csv \
  --cas-output reports/real_hsv2/cas_offinder_output.tsv \
  --out-dir reports/real_hsv2/final
```

## Required next validation

- compare with published anti-HSV targets and tested negatives;
- add variant-aware human off-target analysis;
- separate gene annotation from experimentally supported essentiality;
- assess accessibility in latency-relevant models;
- obtain experimental collaboration before claiming biological activity.
