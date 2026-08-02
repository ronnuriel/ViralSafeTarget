# Release checklist

This checklist keeps a public ViralSafeTarget release deliberate and reproducible.
Publishing is never triggered by an ordinary branch push or by creating a version tag.

## 1. Freeze the release candidate

- Confirm that the release commit is on `main` and the worktree is clean.
- Confirm one version in `pyproject.toml`, `CITATION.cff`, `CHANGELOG.md`, and
  `viral_safe_target.__version__`.
- Confirm author metadata, license, project URLs, and the scientific boundary language.
- Confirm that frozen HSV-2 snapshots have not changed without a documented bug report.

## 2. Validate locally

```bash
python -m pip install --upgrade build twine
python -m build
python -m twine check dist/*
pytest -q
ruff check src tests scripts
```

Install the wheel in a new virtual environment outside the repository and run:

```bash
vst --version
vst doctor
vst quickstart --out demo-project
vst plan demo-project/project.yaml
vst run demo-project/project.yaml
vst open demo-project/results --no-browser
vst export demo-project/project.yaml
```

## 3. Build in GitHub Actions

Run **Build release artifacts** manually with `publish_target=none`. Download and inspect
the wheel and source distribution produced on Linux and smoke-tested on Linux and macOS.

## 4. TestPyPI

Configure a protected GitHub environment named `testpypi` and its Trusted Publisher.
Run the workflow manually with `publish_target=testpypi`. Test installation in a clean
environment, allowing dependencies to come from the production index when necessary.

Do not interpret a successful package installation as validation of external tools or
scientific outputs.

## 5. Production approval

- Obtain explicit approval from the project maintainers.
- Create an annotated `vX.Y.Z` tag on the approved commit.
- Select that tag in GitHub Actions and manually run **Build release artifacts** with
  `publish_target=pypi`.
- The production job is guarded: it runs only for a manual dispatch on a `v*` tag and
  through the protected `pypi` environment.
- Create the GitHub Release and Zenodo archive only after the uploaded artifacts and
  checksums have been verified.

## 6. Post-release

- Install `viral-safe-target` from public PyPI in a new environment.
- Re-run the five-minute demonstration.
- Update the README only after the public project page resolves.
- Record the DOI and immutable release URL in `CITATION.cff`.

No release step establishes editing efficacy, host safety, viral inhibition, delivery,
treatment, clearance, or cure.
