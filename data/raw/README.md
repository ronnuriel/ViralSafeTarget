# Raw data

Large biological datasets are intentionally excluded from Git.

Use `scripts/run_real_hsv2.sh` or the NCBI Datasets CLI to populate this
directory locally. Record accessions, versions, checksums and download dates in
the generated run manifest. Do not commit GRCh38, FASTQ/BAM files, or large
multi-genome FASTA files to the repository.
