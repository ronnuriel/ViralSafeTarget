# Target-to-intervention evidence

ViralSafeTarget answers two different research questions without merging them into a
therapeutic score:

1. **Genetic-editing route:** which conserved viral sites and genes are technically
   targetable under the configured sequence and host-search model?
2. **Pharmacological route:** which source-linked drugs, investigational compounds,
   resistance observations, or pathway associations have been reported for those genes?

## Generic workflow

Every project created with `vst init` contains:

```text
evidence/intervention_evidence.tsv
```

Researchers add source-linked records and run:

```bash
vst evidence interventions --project project.yaml
vst run project.yaml
```

The command validates controlled fields and source URLs, refuses to accept an approved
record without a named reviewer and review date, and generates:

- `intervention_evidence_validated.tsv`;
- `intervention_landscape.tsv`;
- `intervention_review_queue.tsv`;
- `intervention_queries.tsv`;
- `intervention_landscape.html`.

The final `START_HERE.html` report links these files and adds the same evidence to each
shortlisted guide explanation.

## What a record means

Each row identifies a virus, gene, intervention, relationship, development status,
evidence category, evidence scope, source, finding, and review state. Important
relationships include:

- `direct_target`: evidence indicates that the intervention acts on the gene product;
- `prodrug_activation`: the gene product activates an intervention but is not necessarily
  its inhibitory binding target;
- `resistance_marker`: variation is associated with or evaluated for resistance;
- `complex_component`: the protein belongs to a targeted complex, without claiming that
  every subunit is directly bound;
- `pathway_association`: indirect pathway evidence only.

## Mandatory human review

New records start as `pending` and are displayed with `*`. A record can be `approved`
only when `reviewer` and `review_date` are populated. Approval means the source-to-record
extraction was checked; it does not mean that the intervention is safe, effective, or
appropriate for a patient.

Direct evidence in the target virus remains distinct from ortholog or related-virus
evidence. Missingness is preserved: `not located in supplied evidence` does not mean that
no intervention exists.

## HSV-2 example

The case-study table illustrates three different situations:

- **UL30:** a classical pharmacological route linked to activated acyclovir/valacyclovir
  through the viral DNA polymerase mechanism;
- **UL5/UL52:** an investigational helicase-primase route linked to pritelivir, with
  clinical, resistance, and HSV-1 structural records kept separate;
- **UL3:** high computational targetability, but no intervention record located in the
  supplied table.

All case-study records remain pending until a named virology, pharmacology, or clinical
reviewer checks the primary source. The example does not recommend a drug combination or
genome-editing treatment.

## Route toward biological or therapeutic validation

The software can prioritize hypotheses and expose evidence gaps. Moving toward a treatment
still requires independent work outside the repository: expert evidence review,
retrospective validation, target-engagement studies, antiviral phenotype measurements,
host-toxicity assessment, delivery research for editing approaches, resistance analysis,
and appropriately governed preclinical and clinical studies.
