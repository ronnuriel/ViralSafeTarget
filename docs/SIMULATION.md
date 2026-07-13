# Sequence-level disruption simulation

## Question answered

> If two canonical SpCas9 cuts occurred at the selected sites, which interval of the reference sequence would lie between them, and which annotations would overlap that interval?

## Calculation

For each candidate, ViralSafeTarget derives an idealized cut boundary approximately three bases upstream of the PAM. For a pair:

1. order the two cut boundaries;
2. calculate the intervening reference interval;
3. report deletion length and modulo-three arithmetic;
4. intersect the interval with GFF3 features;
5. calculate the fraction of input strains containing both exact 23-nt target sites.

## What the score means

`sequence_disruption_score` is a transparent convenience ranking based on sequence coverage, annotation overlap and deletion size. It is not experimentally calibrated and must not be interpreted as a probability of editing or viral inactivation.

## Missing biology

The calculation does not model delivery, chromatin access, cleavage efficiency, DNA repair, toxicity, immune effects, viral replication, latency or reactivation. Those require separate prediction tools and experimental evidence.
