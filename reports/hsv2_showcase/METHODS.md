# Methods

- Virus profile: `hsv2` (Human alphaherpesvirus 2; NC_001798.2).
- Host profile: `human_grch38` (GRCh38.p14).
- Nuclease profile: `spcas9` (SpCas9).
- Host search: Cas-OFFinder v2.4.1 (Jul 13 2026) through 3 mismatches.
- Candidate shortlist: non-dominated fronts over targetability, predicted protein disruption, and evidence coverage, followed by a fixed per-gene quota.
- Biological evidence: source-linked rows only; HSV-1 and HSV-2 scopes are never merged.
- External validation: population-genomics findings are reported separately from essentiality evidence and computational scores.
- Held-out validation: discovery accessions are excluded; exact guide/PAM retention is evaluated only where a high-quality reference mapping covers the locus.
- Repair outcomes: deterministic size-defined sequence hypotheses, not repair-frequency predictions.
- Profile validation and input checksums are recorded in `run_manifest.json`.
