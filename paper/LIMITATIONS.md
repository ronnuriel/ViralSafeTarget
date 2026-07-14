# Limitations

- No guide was experimentally tested by this project.
- No editing, repair-frequency, viral-inactivation, safety, efficacy, delivery,
  treatment, latency-clearance, or cure conclusion is supported.
- Indels are a bounded, equally weighted size grid; insertion sequence is unspecified.
- Exact-target counterfactuals do not model mismatch-tolerant activity after mutation.
- The multiplex barrier is not an evolutionary probability and omits mutation rate,
  selection, viral fitness, linkage, and population dynamics.
- Held-out coverage is missing for 57 of 257 deep-panel guides.
- Marginal per-guide coverage does not prove joint panel coverage within genomes.
- Predicted host matches are limited by the configured assembly, editor, mismatch/PAM
  model, and software behavior; zero predicted matches do not establish safety.
- Protein-region annotations exist only for a subset of genes.
- Evidence Agent proposals remain pending until explicit human review.
- HSV-1 or other-virus evidence cannot be treated as direct HSV-2 evidence.
- Generic architecture has not yet been demonstrated with a completed second-virus
  public case study.
