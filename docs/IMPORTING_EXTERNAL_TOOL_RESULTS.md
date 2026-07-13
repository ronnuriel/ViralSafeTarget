# Importing external tool results

CRISPOR, CHOPCHOP and GuideScan2 are supported through researcher-supplied CSV, TSV or JSON exports.
Public interfaces are not scraped.

```bash
vst tools import --tool crispor --input export.tsv \
  --mapping configs/import_maps/crispor.yaml \
  --candidates reports/hsv2_pilot/candidates_ranked_post_human.csv \
  --out-dir reports/hsv2_consensus/crispor
```

Mapping templates identify guide, PAM, strand and coordinate columns and document known metric
directions. Matching uses stable candidate IDs. Unmatched and ambiguous rows are retained. Unknown
columns become named raw metrics without invented normalization. Imported source rows and SHA-256
provenance remain available for audit.
