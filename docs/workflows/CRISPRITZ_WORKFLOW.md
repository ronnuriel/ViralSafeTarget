# CRISPRitz workflow

Generate inputs for the existing 32-candidate HSV-2 set:

```bash
vst tools crispritz build-input \
  --run-dir reports/hsv2_pilot \
  --out-dir reports/hsv2_consensus/crispritz \
  --dna-bulges 1 --rna-bulges 1
```

Add `--reference-genome` and `--pam-file` to prepare execution. Inspect the command without running:

```bash
vst tools crispritz run \
  --input-dir reports/hsv2_consensus/crispritz \
  --out-dir reports/hsv2_consensus/crispritz --dry-run
```

Native execution is preferred when available. On macOS the adapter can construct a Docker command,
but it never installs Docker or silently pulls an image. Import an existing output with:

```bash
vst tools crispritz import --results results.tsv \
  --manifest reports/hsv2_consensus/crispritz/crispritz_manifest.json \
  --out-dir reports/hsv2_consensus/crispritz
```

The parser retains mismatch, DNA/RNA bulge, annotation, reference/variant, population and sample
fields when present. A missing result file is pending, not evidence of zero predicted hits.
