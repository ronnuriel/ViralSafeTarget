# Script catalog

Prefer the `vst` CLI for new research projects. These scripts preserve auditable HSV-2
case-study runs, data preparation, and maintenance operations.

## Researcher-facing launchers

| Script | Purpose |
|---|---|
| `run_demo.sh` | Tiny synthetic smoke test |
| `run_synthetic_e2e.sh` | Synthetic end-to-end integration run |
| `run_real_hsv2.sh` | Download/prepare a sampled public HSV-2 run |
| `run_hsv2_pilot.sh` | Focused UL19/UL30 pilot |
| `run_hsv2_consensus.sh` | Multi-tool consensus case study |
| `run_hsv2_genome_wide.sh` | Balanced genome-wide screen |
| `run_hsv2_exhaustive.sh` | Confirmed exhaustive host screen |
| `run_hsv2_population_validation.sh` | Held-out population workflow |
| `build_hsv2_showcase.sh` | Presentation-ready result bundle |
| `finalize_hsv2_exhaustive.sh` | Resume/finalize the long exhaustive run |
| `run_hsv2_tool_benchmark.sh` | Frozen 257-guide systematic multi-tool benchmark |
| `doctor.sh` | Environment and external-tool diagnostics |

## Data and reference preparation helpers

- `download_ncbi_data.sh`
- `download_hsv2_pilot.py`
- `fetch_reference_genbank.py`
- `prepare_real_hsv2.py`
- `prepare_hsv2_population_panel.py`
- `generate_real_candidates.py`
- `run_mafft.sh`

## Workflow implementation helpers

- `run_pipeline.py`
- `run_hsv2_pilot.py`
- `run_hsv2_consensus.py`
- `run_population_validation.py`
- `run_locus_aware_population_validation.py`
- `build_population_validation_report.py`
- `resume_cached_cas_batches.py`

## Audit and maintenance helpers

- `cache_stage.py`
- `compare_discovery_modes.py`
- `summarize_cas_offinder.py`
- `validate_release_candidate.sh` — build, metadata-check, clean-install, and execute the
  complete five-minute researcher smoke flow outside the repository
- `write_hsv2_pilot_manifest.py`

Helpers may have stricter input assumptions than the `vst project` interface. Before
running a long external stage, use `vst project validate`, `vst doctor`, or
`bash scripts/doctor.sh` and inspect the generated command/input files.
