# Local clean-environment validation

ViralSafeTarget does not require a project server. Validate releases in a fresh local
virtual environment so that cached research files or an editable checkout cannot hide
packaging problems.

## Baseline checks

```bash
python3 -m venv .release-venv
.release-venv/bin/python -m pip install --upgrade pip build twine
.release-venv/bin/python -m build
.release-venv/bin/python -m twine check dist/*
```

Install the resulting wheel from a directory outside the checkout, then run:

```bash
vst --version
vst doctor
vst quickstart --out demo-project
vst plan demo-project/project.yaml
vst run demo-project/project.yaml
vst evidence interventions --project demo-project/project.yaml
vst open demo-project/results --no-browser
vst export demo-project/project.yaml
```

Execute the canonical notebook in both `demo` and verified `hsv2_snapshot` modes before
a release. The committed notebook must remain output-free.

## External tools

Run:

```bash
vst doctor
vst tools setup
```

MAFFT, Cas-OFFinder, CRISPRitz, Docker, GPU/OpenCL support, and large host references are
optional external capabilities. Missing components must remain explicit. A local partial
run is valid when unavailable stages remain visible and are never interpreted as zero or
favorable results.

Record the operating system, CPU/GPU mode, memory, disk, tool versions, container digests,
commands, measured stage timings, and output checksums in the run manifest.
