# Output interpretation

## Candidate tables

`candidates_ranked_pre_human.csv` ranks retained candidates before any human search.
Every component is visible. `gene_evidence_score` remains null when no curated
record exists; it is never silently replaced by `1.0`.

`candidates_rejected_pre_human.csv` contains candidates that crossed a configured
threshold, together with semicolon-separated `rejection_reasons`. Rejection is a
computational filtering decision, not proof that a sequence cannot work.

`candidates_ranked_post_human.csv` keeps `pre_human_score` unchanged and adds a
separate `post_human_score`. `predicted_offtarget_risk` summarizes enumerated hits
under the configured editor and mismatch limit. “No predicted hit” means only that
the selected search did not report one; it does not establish safety.

## Conservation and effectiveness are different

`exact_strain_coverage` is the fraction of supplied aligned genomes containing the
exact protospacer-plus-PAM site. It depends on sampling and alignment quality. A
conserved sequence can still be inaccessible or inefficient, while a less conserved
sequence may work in a subset of strains. Conservation is not experimental activity.

## Pair hypotheses

Same-gene and same-region rows may contain an idealized intervening deletion only
when both sites share a reference molecule and satisfy the distance threshold.
Cross-gene rows are `multi_target_hypothesis` entries and do not claim one physical
deletion. Even a correctly calculated sequence deletion does not prove viral
inactivation or predict a repair outcome.
