# Results draft

## Exhaustive source validation

Six assertions were checked against committed source tables before the new snapshot was
written: 28,578 initial candidates, 23,108 eligible rows, 21,654 unique guide
sequences, 440,341 predicted human-match rows, 2,668 zero-hit rows under the configured
model, and 109 completed batches. All assertions passed.

## Deep-panel virtual analysis

The configured deep panel contained 257 unique guides. Generic mapping produced 271
guide-to-CDS rows because overlapping annotations were preserved; 250 guides had at
least one CDS mapping. The -10 through +10 bp grid produced 5,691 sequence hypotheses.
Optional protein annotations were not imputed when unavailable.

## Observed variation and counterfactuals

The analysis enumerated 17,733 single-nucleotide protospacer/PAM counterfactuals. Every
SpCas9 NGG target contributed 69 substitutions, of which 66 removed the exact target:
60 protospacer substitutions and six substitutions at the two constrained PAM
positions. Held-out exact-target coverage was available for 200 of 257 guides; 57
guides remained unknown because the committed held-out snapshot did not cover them.

## Configured strategy comparison

Each of the four configured three-guide panels required three distinct substitutions
to remove all exact targets under the exact set-cover model. All panel members had
discovery exact-target coverage of 1.0 and zero predicted host matches under the
configured source model. Minimum held-out coverage was 0.973913 for the top-ranking,
targetability-focused, and mechanism-diverse panels, and 0.982609 for the
essential/replication-focused panel. These marginal guide-level values do not establish
joint panel coverage or safety.

## Gene-context observations

UL3 contributed 20 analyzed guides, with 10 cutting an annotated UL3 domain; direct
HSV-2 biological evidence remains incomplete. UL52 contributed 11 guides and UL30 20,
with 11 and 16 respectively cutting configured domains, while their replication-focused
rationales remain separate evidence fields. All three UL18 guides cut its annotated
domain and were positioned in the earlier portion of the protein (median relative
position 0.192). UL36 contributed seven guides, including the leading individual source
candidate, while its exhaustive gene-level portfolio rank remained 21.

These observations prioritize questions for review; they do not establish viral
phenotype or therapeutic relevance.
