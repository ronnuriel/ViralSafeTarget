# What the simulator means — and what it does not

## Included simulation

ViralSafeTarget performs a deterministic sequence-level calculation:

1. infer the canonical SpCas9 cut boundary from guide orientation and PAM location;
2. choose two candidate sites;
3. calculate the reference interval between the two cut boundaries;
4. report which annotated features overlap that interval;
5. measure whether both exact target sites are present across the input viral strains.

## Not simulated

The project does not model:

- delivery to infected neurons;
- latent viral chromatin accessibility;
- on-target cleavage probability;
- DNA-repair outcome distributions;
- off-target cleavage frequencies;
- cell viability, immune responses or toxicity;
- viral replication, reactivation or clinical shedding.

## How real evidence is added

- **Before experiments:** use guide-efficiency and off-target predictors as separate evidence layers.
- **After editing experiments:** analyze amplicon sequencing with CRISPResso2 or an equivalent validated workflow.
- **For antiviral effect:** measure appropriate virological endpoints in a qualified laboratory.
- **For latency/cure claims:** use relevant neuronal and in-vivo models; a sequence score is never sufficient.

A good report should label every field as observed, predicted, inferred, or experimentally measured.
