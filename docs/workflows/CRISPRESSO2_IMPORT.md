# CRISPResso2 measured-result import

CRISPResso2 output represents measured sequencing data, not another prediction score. Import an
existing result directory without any experimental procedure:

```bash
vst experimental import-crispresso2 \
  --input CRISPResso2_output_directory \
  --candidate-map candidate_map.csv \
  --out-dir reports/my_run/experimental
```

When present, the adapter extracts aligned reads, modified/inserted/deleted/substituted percentages,
frameshift and in-frame percentages, quantification windows, source hashes and run metadata. These
measurements are written under `experimental/` and are not automatically aggregated with predicted
metrics.
