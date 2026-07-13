#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/raw

# Reference HSV-2 nucleotide record, with genome and annotation.
datasets download virus genome accession NC_001798.2 \
  --include genome,cds,protein,annotation \
  --filename data/raw/hsv2_reference.zip

# Complete HSV-2 genomes. This may be a larger download.
datasets download virus genome taxon "Human alphaherpesvirus 2" \
  --complete-only --include genome,annotation \
  --filename data/raw/hsv2_complete.zip

# Human GRCh38.p14 reference genome.
datasets download genome accession GCF_000001405.40 \
  --include genome,gff3 \
  --filename data/raw/human_GRCh38.zip

echo "Downloads completed. Inspect and unzip packages before analysis."
