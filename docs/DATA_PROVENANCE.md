# Data provenance

Every real run writes `run_manifest.json` with the UTC timestamp, Git commit and
dirty status, configuration checksum, editor profile, input SHA-256 hashes,
accepted and rejected accessions with reasons, human assembly identifier, command
line, package/tool versions, random seed, and output hashes.

The accession QC table is `accession_qc.csv`; the frozen selected list is
`accessions_used.txt`. A researcher should archive both beside the manifest.

Downloaded genomes, human references, local/private reports, temporary alignments,
and Cas-OFFinder result files are ignored by Git. Small synthetic fixtures remain
tracked so CI can reproduce behavior without network access.

Before publication, copy the configuration and manifest into the research archive,
record the exact external-tool binaries, inspect alignment quality, and retain the
original public-accession list. A manifest records what ran; it does not certify that
the sampling design or biological interpretation is valid.
