# Data contract

## Required input

### Aligned viral multi-FASTA

- At least two records.
- Unique record IDs.
- Equal aligned lengths.
- One record ID selected as the reference.
- DNA alphabet should primarily contain `A`, `C`, `G`, `T`, `N`, and `-`.

## Optional input

### Reference GFF3

The `seqid` must match the selected reference FASTA record. Core columns used:

- `seqid`
- `feature_type`
- `start`, `end`
- `strand`
- attributes such as `ID`, `Name`, `gene`, `product`

### Host FASTA

The built-in exhaustive screen is intentionally limited to a small FASTA of at most 5 Mb. Use Cas-OFFinder, CRISPRitz or another validated genome-scale engine for GRCh38.

## Candidate output

Core columns include:

- `candidate_id`
- `guide_sequence`
- `pam`
- `strand`
- `reference_start_1based`, `reference_end_1based`
- `site_start_1based`, `site_end_1based`
- `virus_site_coverage`
- `gene_name`, `product`, `feature_type` when GFF3 is supplied
- host-screen fields when available
- `decision` and transparent component scores

## Pair simulation output

- candidate IDs
- idealized canonical cut boundaries
- idealized deletion interval and length
- exact pair strain coverage when alignment is supplied
- annotation overlap summary
- transparent, unvalidated sequence-disruption score
