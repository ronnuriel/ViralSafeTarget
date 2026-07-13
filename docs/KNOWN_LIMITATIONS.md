# Known limitations

- SpCas9/NGG is the only tested editor profile. The schema is extensible, but other
  profiles are not claimed as supported.
- Exact protospacer-plus-PAM coverage is intentionally strict, sensitive to viral
  sampling, and dependent on alignment quality.
- HSV repeats and genome isomers can create misleading whole-genome alignments.
- Viral occurrence counting is exact-string based and does not model functional
  equivalence or mismatch-tolerant viral editing.
- GFF3 overlap indicates location, not essentiality, accessibility, or effectiveness.
- Curated evidence is absent unless a researcher supplies source-linked rows; unknown
  evidence remains null and supplies no positive score.
- Cas-OFFinder enumerates sequence matches under a configured model. It does not
  establish cleavage probability or safety, and a reference-only human search omits
  population variation unless a separate variant-aware design is used.
- Pair deletion coordinates assume two canonical cuts and an idealized intervening
  deletion. Repair distributions and cleavage co-occurrence are not modeled.
- The transparent scores are not experimentally calibrated probabilities.
- Delivery, latent infection, chromatin accessibility, editing efficiency, repair
  outcomes, toxicity, immune effects, viral reactivation, and clinical outcomes are
  outside the model.
