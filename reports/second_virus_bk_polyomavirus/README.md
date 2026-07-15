# Second-virus usability proof: BK polyomavirus

This compact snapshot was generated with the installed researcher-facing CLI, without
virus-specific core-code changes:

```bash
vst init bk-polyomavirus \
  --virus-name "BK polyomavirus" \
  --tax-id 1891762 \
  --reference-accession NC_001538.1
vst plan bk-polyomavirus/project.yaml
vst run bk-polyomavirus/project.yaml
```

The accession was retrieved through NCBI E-utilities. `input_provenance.json` records
the URLs, retrieval time, observed accession, and checksums. This is a reference-only
usability run: it does **not** establish population conservation.

## Result

- 543 editor-compatible sites were estimated before the configured 500-row discovery cap.
- 500 candidate rows entered the public result bundle.
- Four annotated gene labels were represented in the bundle.
- Gene/CDS mapping, virtual-knockout hypotheses, escape counterfactuals, pair analysis,
  multiplex comparison, reporting, resume/cache, and export completed generically.
- Host screening remained `external_required` because a host assembly was not supplied.
- Evidence discovery was not run; no biological claim was inferred.

This is a software usability/generalization test, not a therapeutic case study. Open
`START_HERE.html` first.
