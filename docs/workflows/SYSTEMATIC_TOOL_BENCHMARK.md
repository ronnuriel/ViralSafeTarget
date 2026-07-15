# Systematic multi-tool benchmark

The publication benchmark uses one frozen candidate identity table and separates two questions:

1. **Executable comparison:** which candidates and primary metrics were actually reported by each tool?
2. **Capability evidence:** which functions are documented by official papers, documentation, or source repositories?

These questions must not be collapsed. A documented capability is not a completed result, and missing output is not zero predicted risk.

## HSV-2 benchmark

```bash
vst tools benchmark --config configs/benchmarks/hsv2_multitool.yaml
```

The frozen panel contains 257 unique candidate IDs and 257 unique guide sequences. The command writes normalized results, rank agreement, top-K overlap, execution status, leave-one-component-out sensitivity, figures, an HTML report, and a checksum manifest.

## Tool input bundles

The workflow writes documented input bundles under `inputs/`:

- `frozen_panel.csv`: complete frozen identity and source attributes;
- `candidate_identity.tsv`: compact ID/sequence/coordinate map;
- `crispritz_guides.txt`: protospacers with the configured PAM placeholder;
- `crispor_batch.tsv`: candidate IDs and guide sequences;
- `chopchop_targets.fasta`: named target sequences;
- `guidescan2_guides.txt`: one guide per line;
- `input_manifest.json`: checksums and candidate invariants.

Public web interfaces are not scraped. CRISPOR, CHOPCHOP, and GuideScan2 remain `export_required` until a researcher supplies an official export. The HSV-2 snapshot includes a completed CRISPRitz 2.6.6 Docker run on the frozen panel against GRCh38.p14 through three mismatches, without bulges or population variants. Each import must retain the raw file, version, command or import method, assembly, editor, date, and SHA-256 hash.

The completed executable comparison currently contains ViralSafeTarget pre-host and post-host rankings, Cas-OFFinder source counts, and CRISPRitz counts. The capability matrix includes all named tools but does not imply that an executable export exists.

## Comparable metrics

Raw values from different tools are never averaged. The report compares:

- candidate coverage and missingness;
- within-tool rank and percentile;
- pairwise Spearman rank agreement when both tools have non-constant ranks;
- top-10, top-25, and top-50 overlap;
- off-target burden only when search models are documented;
- runtime only when hardware, genome, guide count, and search parameters are comparable.

## Ablation

The benchmark removes one ViralSafeTarget scoring component at a time and recalculates the transparent weighted score. This quantifies rank sensitivity on the frozen panel. It is not model retraining, biological validation, or evidence that any component causes editing activity.

## Interpretation boundary

The benchmark can establish reproducibility, coverage, disagreement, and scope. It cannot establish editing efficiency, safety, viral phenotype, delivery, treatment, or cure. Predictive superiority requires independent experimental ground truth.
