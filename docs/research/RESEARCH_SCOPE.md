# Research scope

ViralSafeTarget is a virus-first computational target-discovery and prioritization
pipeline. It produces inspectable candidate tables, predicted human off-target
summaries, sequence-level pair hypotheses, benchmarks, and provenance records.

The unit of output is a **computational candidate**, not a validated guide. A high
rank indicates agreement with configured sequence and annotation criteria. It does
not prove cleavage, effectiveness, viral inactivation, safety, or clinical value.

## In scope

- exact target-site conservation across aligned viral genomes;
- reference occurrence counting and duplicate-guide disclosure;
- editor/PAM compatibility under a versioned profile;
- reference GFF3 overlap;
- explicitly curated, source-linked gene evidence when available;
- Cas-OFFinder input generation and predicted-hit summarization;
- idealized deletion and independent multi-target hypotheses;
- deterministic ranking, benchmarking, caching, and provenance.

## Outside the model

Delivery, latent chromatin accessibility, editing efficiency, repair outcomes,
toxicity, immune effects, viral reactivation, and clinical outcomes remain outside
the current model. No output is intended for clinical use or as a wet-lab protocol.
