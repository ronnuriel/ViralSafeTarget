# Reproducing the HSV-2 computational case study

This document is the canonical entry point for the published HSV-2 case study.
Historical scripts and notebooks remain available for provenance, but researchers
should not need to infer their execution order.

The workflow produces computational hypotheses. It does not provide a wet-lab
protocol or establish editing, viral inactivation, safety, efficacy, delivery,
latency clearance, or clinical utility.

## 1. Install and audit the environment

```bash
git clone https://github.com/ronnuriel/ViralSafeTarget.git
cd ViralSafeTarget
conda env create -f environment.yml
conda activate viral-safe-target
python -m pip install -e '.[population,notebooks]'

vst doctor
vst tools doctor
pytest -q
ruff check .
```

The full host screen additionally requires a local Cas-OFFinder executable and
the GRCh38 download performed by the reproduction workflow.

```bash
export CAS_OFFINDER_PATH=/absolute/path/to/cas-offinder
```

## 2. Inspect the complete plan without downloading or computing

```bash
vst reproduce hsv2
```

The plan reports every stage, its command, purpose, expected output, and whether
a cached output is already present.

## 3. Execute or resume the complete case study

```bash
vst reproduce hsv2 --execute
```

The workflow:

1. downloads the frozen HSV-2 discovery accessions listed in
   `data/curated/hsv2_discovery_accessions.txt`;
2. validates and aligns the viral records;
3. performs the balanced genome-wide host screen;
4. downloads the declared HSV-1 ortholog reference;
5. validates a discovery-excluded population panel;
6. maps selected sites to proteins and computes size-defined disruption
   hypotheses;
7. builds the final multi-objective report and a reproduction manifest.

Existing checksum-valid stages and Cas-OFFinder batches are reused. Re-running
the command therefore resumes rather than silently restarting long computations.

Use `--skip-population` only for a partial reproduction. The resulting run must
be described as partial.

## 4. Primary outputs

```text
reports/hsv2_genome_wide/
reports/hsv2_population_report_balanced/
reports/hsv2_gene_function/
reports/hsv2_showcase/FINAL_REPORT.html
reports/hsv2_reproduction/reproduction_manifest.json
```

Generated reports and large public references are intentionally not committed to
Git. Manifests record input and output checksums. A changing external database or
tool version may still change a result; the frozen accession list prevents silent
changes to the discovery cohort.

## 5. Notebook trail

- `09_HSV2_GENOME_WIDE_DISCOVERY_EN.ipynb`: genome-wide discovery logic.
- `10_HSV2_GENE_FUNCTION_AND_DISRUPTION_EN.ipynb`: protein mapping and disruption.
- `11_HSV2_RESEARCH_SHOWCASE_EN.ipynb`: presentation-ready analysis.
- `12_HSV2_HELDOUT_POPULATION_VALIDATION_EN.ipynb`: held-out population validation.

The CLI workflow is authoritative for execution. Notebooks explain and inspect
the generated artifacts; they are not hidden state required to produce them.
