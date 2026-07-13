# ViralSafeTarget

**A reproducible, virus-first research toolkit for conserved genome target prioritization, host off-target triage, and sequence-level disruption simulation.**

ViralSafeTarget is designed for computational researchers who want to start with a collection of viral genomes and produce an auditable table of candidate target sites. The first case study is HSV-2, but the core pipeline is intended to be reusable for other DNA viruses.

> **Scope:** this repository produces computational hypotheses. It does not prove editing, viral inactivation, safety, latency clearance, or clinical efficacy. It contains no wet-lab protocol.

## What the pipeline does

```mermaid
flowchart LR
    A[Viral genomes<br/>aligned FASTA] --> B[Conservation scan]
    B --> C[Editor-compatible sites<br/>guide + PAM]
    C --> D[Reference annotation<br/>GFF3]
    D --> E[Host off-target screen]
    E --> F[Ranked candidates]
    F --> G[Idealized cut/deletion<br/>sequence simulation]
    G --> H[Experimental validation<br/>outside this repository]
```

It answers questions such as:

- Is this target sequence present across many viral strains?
- Which annotated viral feature contains the target?
- Does a selected editor have a compatible PAM at that location?
- Are there similar host sites that require off-target review?
- If two idealized cuts occurred, which reference interval would be deleted?

It does **not** answer whether delivery succeeds, whether chromatin is accessible, how often a repair outcome occurs, whether the virus remains viable, or whether a treatment is safe.

## Quick start

### Conda (recommended)

```bash
conda env create -f environment.yml
conda activate viral-safe-target
python -m pip install -e .
pytest -q
ruff check .
bash scripts/run_demo.sh
```

Open the generated files under `reports/demo/`.

### Local upload interface

```bash
streamlit run app.py
```

Upload an **already aligned** multi-FASTA and, optionally, a matching GFF3 annotation. The browser interface is intended for local research use.

### Jupyter

```bash
jupyter lab
```

Start with `notebooks/00_START_HERE_HE.ipynb` or the English README and documentation.

## Run the HSV-2 pilot on public data

```bash
bash scripts/run_real_hsv2.sh --sample-size 25
```

This validates or downloads public NCBI records, resumes checksum-cached stages,
selects a deterministic pilot sample, aligns it with MAFFT, and produces explainable
pre-human rankings plus corrected pair hypotheses.

To also download GRCh38 and prepare a Cas-OFFinder input file:

```bash
bash scripts/run_real_hsv2.sh --with-human --sample-size 25
```

The human genome is intentionally not committed to Git.

For the focused UL19/UL30 workflow, including deterministic Cas-OFFinder selection:

```bash
bash scripts/run_hsv2_pilot.sh
```

See [`docs/HSV2_PILOT.md`](docs/HSV2_PILOT.md) for the two-stage off-target run.

## Run the v0.4 multi-tool consensus

The consensus pilot reuses the completed v0.3 outputs and selects only the 32 candidates with no
predicted human hit within the configured Cas-OFFinder model and three-mismatch threshold:

```bash
bash scripts/run_hsv2_consensus.sh
```

It generates CRISPRitz inputs, imports researcher-supplied CRISPOR/CHOPCHOP/GuideScan2 exports when
present, preserves missing stages, and writes `reports/hsv2_consensus/report.html`. See
[`docs/HSV2_CONSENSUS_PILOT.md`](docs/HSV2_CONSENSUS_PILOT.md),
[`docs/PYTHON_SDK.md`](docs/PYTHON_SDK.md), and
[`docs/INTEGRATIONS.md`](docs/INTEGRATIONS.md).

## Main outputs

- `candidates_ranked_pre_human.csv`: retained candidates with visible score components.
- `candidates_rejected_pre_human.csv`: rejected candidates and explicit reasons.
- `candidates_ranked_post_human.csv`: separate post-human score when results exist.
- `predicted_human_hits.csv`: parsed predicted hits when results exist.
- `pair_hypotheses_same_gene.csv`: bounded same-gene deletion hypotheses.
- `pair_hypotheses_multi_target.csv`: independent cross-gene target hypotheses.
- `report.html`: human-readable summary.
- `run_manifest.json`: input checksums, parameters, environment and Git commit.

## What “simulation” means here

The included simulator calculates canonical SpCas9 cut coordinates and the sequence interval that would be removed **if** two selected cuts occurred. It can report overlap with annotated features and exact pair coverage across the input viral alignment.

This is a sequence transformation, not a cell or organism simulator. To measure real editing, researchers need sequencing data and tools such as CRISPResso2. To evaluate potential host off-target sites at genome scale, use a validated engine such as Cas-OFFinder or CRISPRitz. To establish viral inactivation, appropriate virology experiments are required.

See [`docs/SIMULATION_LIMITS.md`](docs/SIMULATION_LIMITS.md) and [`docs/KNOWN_LIMITATIONS.md`](docs/KNOWN_LIMITATIONS.md).

## Curated evidence

The schema is `data/curated/viral_gene_evidence.schema.csv`. Add only source-linked,
reviewed rows to a virus-specific table and pass it with `--gene-evidence`. Supported
status vocabulary includes `supported`, `suggested`, `unknown`, and `conflicting`.
Unknown or missing evidence remains missing and never silently receives a positive
score. See [`docs/OUTPUT_INTERPRETATION.md`](docs/OUTPUT_INTERPRETATION.md).

## Repository structure

```text
ViralSafeTarget/
├── app.py                    # local upload UI
├── configs/                  # example research configurations
├── data/demo/                # synthetic smoke-test data
├── data/knowledge/           # evidence templates, not biological truth tables
├── docs/                     # concepts, formats, validation and tool map
├── notebooks/                # guided research notebooks
├── scripts/                  # demo and real-data workflows
├── src/viral_safe_target/    # reusable Python package
├── tests/                    # unit tests
└── .github/                  # CI and contribution templates
```

## Research contribution target

A useful paper would not claim “we found a cure” from computational scores. A credible contribution would demonstrate that ViralSafeTarget:

1. recovers published viral targets without using them to tune the result;
2. improves strain coverage or host-specificity over existing baselines;
3. produces reproducible decisions with explicit rejection reasons;
4. generalizes to held-out viruses or newly published genomes; and
5. ideally prioritizes candidates later validated by an independent laboratory.

## Existing tools we complement

ViralSafeTarget is an orchestration and provenance layer, not a replacement for mature tools:

- MAFFT — multiple sequence alignment
- Cas-OFFinder — genome-scale off-target enumeration
- CRISPRitz — variant-aware off-target analysis
- GuideScan2 / CRISPOR / CHOPCHOP — guide design and specificity scoring
- CRISPResso2 — analysis of measured genome-editing sequencing data

See [`docs/EXISTING_TOOLS.md`](docs/EXISTING_TOOLS.md).

## Safety and responsible use

Read [`DISCLAIMER.md`](DISCLAIMER.md) and [`SECURITY.md`](SECURITY.md). The project is for legitimate computational research and education. Do not interpret rankings as treatment recommendations or use them as a substitute for institutional biosafety, ethics, clinical, or regulatory review.

## Citation

Use [`CITATION.cff`](CITATION.cff). Before public release, replace the placeholder repository owner and author metadata.

## License

MIT.
