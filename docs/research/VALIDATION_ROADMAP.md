# Validation roadmap

## Phase 0 — software correctness

- Unit tests for both strands, PAM matching and coordinate mapping.
- Deterministic demo output.
- CI on multiple Python versions.
- Explicit software and data versions.

## Phase 1 — retrospective computational benchmark

Build a frozen table of published viral targets and negatives. Do not use the held-out outcomes to tune the final ranking. Measure:

- recall of experimentally active targets;
- precision among the top-ranked sites;
- viral strain coverage;
- host off-target burden;
- agreement and disagreement with CRISPOR/CHOPCHOP/CRISPRitz;
- runtime and reproducibility.

## Phase 2 — temporal validation

Use only information published before a chosen date to rank candidates, then test whether later publications preferentially evaluated or validated high-ranked sites.

## Phase 3 — external laboratory validation

A collaborating virology laboratory selects a blinded, preregistered candidate set and measures editing, intact viral DNA, infectious replication, cell viability and experimentally detected off-targets.

## Phase 4 — latency-relevant validation

Evaluate delivery and editing in a biologically relevant latency/reactivation model. A simple lytic cell culture result is not evidence of latent-reservoir elimination.

## Publication threshold

A software/tool paper may be justified once the pipeline is reproducible and clearly outperforms or complements existing tools on a public benchmark. A therapeutic claim requires experimental evidence well beyond the repository.
