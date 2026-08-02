# Five-minute quick start

## A. Install

### Current reviewed source (before the first PyPI release)

```bash
python -m venv .venv
source .venv/bin/activate
pip install "viral-safe-target[notebooks] @ git+https://github.com/ronnuriel/ViralSafeTarget.git@main"
vst --version
vst doctor
```

### Stable PyPI release

Use this form only after the release is visible on the public PyPI project page:

```bash
python -m venv .venv
source .venv/bin/activate
pip install "viral-safe-target[notebooks]"
vst --version
vst doctor
```

MAFFT, Cas-OFFinder, and CRISPRitz are explicit external programs. A pip install does
not pretend they are present.

## B. Run the synthetic demonstration

```bash
vst quickstart --out demo-project
vst plan demo-project/project.yaml
vst run demo-project/project.yaml
vst open demo-project/results
```

Open `demo-project/results/START_HERE.html`. The demonstration is computational and
contains no clinical or wet-lab claim.

## C. Start from an NCBI accession

```bash
vst init my-virus \
  --virus-name "Virus scientific name" \
  --tax-id TAX_ID \
  --reference-accession ACCESSION
vst plan my-virus/project.yaml
```

ViralSafeTarget retrieves the reference and GenBank annotation through NCBI E-utilities,
derives a reference-matched GFF3, and records URLs, dates, versions, and SHA-256 hashes.
A reference-only run does not establish population conservation; add a strain panel for
that question.

## D. Start from local files

```bash
vst init my-virus \
  --virus-name "Virus scientific name" \
  --reference-fasta inputs/reference.fasta \
  --annotation-gff inputs/reference.gff3 \
  --strains-fasta inputs/strains.aligned.fasta \
  --host-fasta inputs/host.fasta
```

For an unaligned panel, add `--strains-unaligned`; MAFFT must be available. If no
annotation exists, add `--sequence-only`. Gene, protein, domain, and evidence analyses
then remain explicitly unavailable.

## E. Continue without external host-search tools

Run the project normally. The host stage will be labeled `external_required`, and all
sequence-only stages and the partial report will still be produced.

```bash
vst tools status
vst tools setup
vst tools setup --tool cas-offinder
vst status my-virus/project.yaml
```

After installing Cas-OFFinder or supplying its output, continue with:

```bash
vst resume my-virus/project.yaml --run-external
```

Missing output is never interpreted as zero predicted hits.

## F. Send results to collaborators

```bash
vst export my-virus/project.yaml
```

Send `results/export.zip` and ask collaborators to open `START_HERE.html`. The portable
archive excludes reference genomes and large raw host-search outputs by default. It
includes the research shortlist, stage status, evidence queue, methods, limitations,
timings, and provenance.
