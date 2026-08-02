# Server validation

Use a fresh clone or worktree for server validation. Do not reuse a research checkout
that contains local notebooks, staging data, or experimental source files.

## Baseline checks

```bash
git clone https://github.com/ronnuriel/ViralSafeTarget.git ViralSafeTarget-release
cd ViralSafeTarget-release
python3 -m venv .venv
.venv/bin/python -m pip install --upgrade pip build twine
.venv/bin/python -m build
.venv/bin/python -m twine check dist/*
```

Install the wheel from a directory outside the checkout, then run the same clean-wheel
smoke flow used by GitHub Actions.

## External tools

Run:

```bash
vst doctor
vst tools setup
```

Missing MAFFT, Cas-OFFinder, CRISPRitz, Docker, GPU/OpenCL support, or a host reference
must remain explicit. Do not repair GPU drivers or install system packages as part of a
scientific workflow run; those are separate administrator actions and may require a
reboot.

Record the server hostname, operating system, CPU/GPU mode, memory, disk, tool versions,
container digests, commands, measured stage timings, and output checksums in the run
manifest. A partial run is acceptable when its missing stages remain visible.
