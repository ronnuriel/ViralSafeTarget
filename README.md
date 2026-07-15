# ViralSafeTarget

ViralSafeTarget is a virus-first computational research toolkit for conserved guide
discovery, explicit host-search status, virtual-knockout sequence hypotheses,
exact-target escape robustness, multiplex comparison, and source-linked evidence review.

> The outputs are computational research hypotheses. They do not establish editing,
> safety, viral inhibition, delivery, treatment, clearance, or cure. No wet-lab protocol
> is provided.

## Five-minute installation

```bash
python -m venv .venv
source .venv/bin/activate
pip install "viral-safe-target[notebooks]"
vst --version
vst doctor
```

Create and run the bundled demonstration:

```bash
vst quickstart --out demo-project
vst plan demo-project/project.yaml
vst run demo-project/project.yaml
vst open demo-project/results
```

See the tested [five-minute quick start](docs/getting-started/FIVE_MINUTE_QUICKSTART.md).

## New-virus quick start

```bash
vst init my-virus --interactive
vst plan my-virus/project.yaml
vst run my-virus/project.yaml
vst status my-virus/project.yaml
vst resume my-virus/project.yaml
vst open my-virus/results
vst export my-virus/project.yaml
```

The wizard supports an NCBI accession or local FASTA/GFF3 files, aligned or unaligned
strain panels, a host profile or local host FASTA, and explicit sequence-only mode.
MAFFT, Cas-OFFinder, and CRISPRitz remain external programs; missing results remain
`external_required`, never zero.

## Canonical notebook

External researchers should use one notebook:

[`notebooks/00_VIRALSAFETARGET_END_TO_END_EN.ipynb`](notebooks/00_VIRALSAFETARGET_END_TO_END_EN.ipynb)

It supports `demo`, `hsv2_snapshot`, and `custom_project` modes and calls the installed
public CLI rather than duplicating algorithms. Earlier scientific and teaching notebooks
are preserved under [`notebooks/advanced/`](notebooks/advanced/) and
[`notebooks/learning/he/`](notebooks/learning/he/).

## Completed HSV-2 results

The frozen exhaustive snapshot reports 28,578 initial sites, 23,108 eligible candidate
coordinates, 21,654 unique guide sequences, 109/109 Cas-OFFinder batches, 440,341
predicted human matches under the configured model, and 2,668 candidate rows without a
predicted match under that model. UL3, UL10, UL52, UL47, and UL11 lead gene-level
targetability; the leading individual guide is `VST-2e9f052157f9bf29` in UL36.

These are targetability observations, not biological or therapeutic rankings. Start at
the [results index](reports/README.md), the
[exhaustive report](reports/hsv2_genome_wide_exhaustive/report.html), the
[virtual-knockout/escape findings](reports/hsv2_virtual_knockout_escape/FINDINGS.md),
and the [systematic benchmark](docs/workflows/SYSTEMATIC_TOOL_BENCHMARK.md).

Reproduction is explicit and resumable:

```bash
vst reproduce hsv2        # plan only
vst reproduce hsv2 --execute
```

## What the workflow keeps separate

1. Guide sequence quality and population conservation.
2. Model-bounded host-search results and missing external output.
3. Gene-level targetability.
4. Size-defined protein-disruption hypotheses.
5. Exact-target sequence escape robustness.
6. Direct-virus, ortholog, pending, and unknown biological evidence.

No combined therapeutic score is produced. The benchmark compares compatible exported
results from CRISPOR, CHOPCHOP, CRISPRitz, GuideScan2 and other tools while preserving
missingness and provenance.

## Documentation and development

- [Documentation](docs/README.md)
- [New-virus workflow](docs/getting-started/NEW_VIRUS_WORKFLOW.md)
- [Notebook index](notebooks/README.md)
- [Scientific boundaries](DISCLAIMER.md)
- [Contributing](CONTRIBUTING.md)

```bash
python -m pip install -e ".[dev,notebooks]"
pytest -q
ruff check .
```

Use [`CITATION.cff`](CITATION.cff) for citation. The project is licensed under the
[MIT License](LICENSE).
