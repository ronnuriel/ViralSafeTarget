# Running ViralSafeTarget on a new virus

The `vst project` interface is the canonical researcher workflow for a new DNA
virus. It separates project data from reusable program code and records the state
of each stage.

## 1. Create a project

```bash
vst project init \
  --id my-virus \
  --display-name "My virus" \
  --reference-accession REF_ACCESSION \
  --out-dir projects/my-virus
```

This creates:

```text
projects/my-virus/
├── project.yaml
├── profiles/
│   ├── virus.yaml
│   ├── host.yaml
│   ├── nuclease.yaml
│   └── ranking.yaml
├── data/
├── evidence/gene_evidence.tsv
├── external/host/
└── results/
```

## 2. Add declared inputs

Add the following files:

```text
data/reference.fasta
data/reference.gff3
data/strains.aligned.fasta
external/host/*.fasta
```

Requirements:

- every aligned viral sequence has the same alignment length;
- the declared reference is present in the alignment;
- the ungapped aligned reference matches the reference FASTA;
- GFF3 `seqid` values use the declared reference identifier and coordinate system;
- the host profile declares the exact assembly searched;
- evidence rows include a source URL and keep evidence from related viruses
  separate from direct evidence in the target virus.

An empty evidence table is valid. Missing evidence remains unknown.

## 3. Validate before computing

```bash
vst project validate --project projects/my-virus/project.yaml
```

Require the large host reference as well with:

```bash
vst project validate \
  --project projects/my-virus/project.yaml \
  --require-host-reference
```

## 4. Run the sequence-level workflow

```bash
vst project run --project projects/my-virus/project.yaml
```

This performs:

- conserved editor-compatible site discovery;
- reference GFF3 mapping, including one-to-many overlaps;
- transparent pre-host sequence ranking;
- balanced per-gene plus global panel selection;
- same-feature and multi-target pair hypotheses;
- Cas-OFFinder input preparation;
- a partial report and checksum-linked manifest.

If Cas-OFFinder is installed:

```bash
export CAS_OFFINDER_PATH=/absolute/path/to/cas-offinder
vst project run \
  --project projects/my-virus/project.yaml \
  --run-external
```

If it is run separately, place its native output at:

```text
projects/my-virus/results/host_screen/cas_offinder_output.tsv
```

Then resume:

```bash
vst project resume --project projects/my-virus/project.yaml
```

Inspect progress at any time:

```bash
vst project status --project projects/my-virus/project.yaml
```

The state vocabulary distinguishes `completed`, `failed`, `pending`, and
`external_required`. Missing external output is never interpreted as a zero-hit
result.

## 5. Current generic boundary

The project workflow is generic for sequence discovery, annotation, balanced
selection, host-screen orchestration, pair hypotheses, reporting, and provenance.

The advanced protein/domain/ortholog analysis used in the HSV-2 case study still
requires a compatible GenBank annotation and virus-specific evidence/domain inputs.
It is not silently executed for a new virus. This boundary is reported explicitly
rather than producing fabricated essentiality or ortholog conclusions.
