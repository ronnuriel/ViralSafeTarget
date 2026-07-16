# BMC Bioinformatics submission workspace

This directory contains the publication build for a **Software Article** describing
ViralSafeTarget 0.10.0 and its computational HSV-2 case study.

The package is deliberately labelled **DO NOT SUBMIT** until every row in
`final/HUMAN_REVIEW_REQUIRED.csv` has a real reviewer name, review date, checked source,
and an approved or edited decision. An asterisk in the working manuscript identifies a
claim that remains subject to that review.

The frozen analysis software commit is `de7868a83d3d1e30729323e53a735615d43fc231`.
Submission-artifact commits do not change the frozen HSV-2 results.

Build commands:

```bash
python scripts/build_bmc_figures.py
python scripts/build_bmc_documents.py
```

The second command requires `python-docx` and is run with the bundled document runtime
in the Codex workspace. Figures are rebuilt from committed CSV/JSON sources.

