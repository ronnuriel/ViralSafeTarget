# Current status and translational path

## Where the project stands

ViralSafeTarget is ready to generate an auditable computational research shortlist for
a new DNA virus from a reference, annotation, strain panel, host genome, and nuclease
profile. It can preserve sequence conservation, host-search status, gene mapping,
predicted coding disruption, exact-target escape counterfactuals, multiplex comparisons,
biological evidence, and intervention evidence as separate axes.

The HSV-2 case study has completed the exhaustive sequence and host-search workflow and
has public result snapshots. It has not demonstrated genome editing, target engagement,
viral inhibition, host safety, delivery, reservoir clearance, or clinical benefit.

## What the current HSV-2 evidence suggests

| Route | Example genes | Current interpretation | Main gap |
|---|---|---|---|
| Established pharmacology | UL23, UL30 | Acyclovir/valacyclovir provides a known activation/polymerase mechanism | It treats infection but does not establish eradication of latent HSV-2 |
| Investigational pharmacology | UL5, UL52, UL8 | Pritelivir supports the helicase-primase complex as a clinically investigated route | Complex membership, direct binding, resistance, and target-virus transfer require expert review |
| Genome-editing hypothesis | UL3, UL10, UL36 and others | The screen identifies technically strong or unusual targetability patterns | No demonstrated delivery, editing, antiviral phenotype, or safety |
| Mechanistically balanced panel | UL30/UL52 plus a complementary target | May reduce dependence on one sequence or biological mechanism | Requires preregistered experimental comparison; the software does not prescribe a combination |

Every biological and intervention claim must be traced to a primary source and reviewed.
The current HSV-2 intervention table contains pending records, not approved conclusions.

## What remains before researchers can pursue a treatment claim

1. **Named human evidence review.** Virologists, pharmacologists, and clinicians must
   approve, correct, or reject the pending source extractions.
2. **Retrospective validation.** Freeze known positive and negative viral targets and
   measure whether the ranking recovers them without tuning on the answers.
3. **Independent laboratory evaluation.** Test a blinded shortlist and suitable controls
   for editing or target engagement, antiviral phenotype, infectious output, cell
   viability, and measured host effects.
4. **Resistance and combination logic.** Compare single and multiplex/combination
   hypotheses under transparent assumptions. Sequence-level escape counts are not
   evolutionary probabilities.
5. **Delivery for editing approaches.** Demonstrate access to the biologically relevant
   infected and latent cell populations. A lytic cell-culture result is not evidence of
   reservoir elimination.
6. **Preclinical and clinical development.** Use appropriately governed studies before
   any statement about efficacy, safety, treatment, or cure.

## What an external researcher receives

After `vst run`, the researcher opens `results/START_HERE.html`. The report provides the
shortlist, gene targetability, host-search state, disruption/escape hypotheses,
biological evidence queue, target-to-intervention landscape, missing stages, provenance,
and exact commands. `vst export` creates a portable collaborator bundle.

This makes the project useful today as a decision-support and evidence-gap tool. Its
largest potential contribution will be demonstrated only when external groups show that
the structured shortlist improves experimental selection, reproducibility, or research
efficiency relative to guide-only workflows.
