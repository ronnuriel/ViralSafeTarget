# ViralSafeTarget: evidence-aware, virus-first CRISPR target prioritization with an exhaustive HSV-2 case study

**Article type:** Software Article
**Journal:** BMC Bioinformatics
**Submission status:** **DO NOT SUBMIT — HUMAN REVIEW PENDING**
**Authors:** Ron Nuriel¹* (ORCID: 0009-0008-3970-2591); Sarel Cohen¹ (ORCID: 0000-0003-4578-1245)
**Affiliation:** ¹[AFFILIATION 1 REQUIRED]
**Corresponding author:** ron.nuriel01@post.runi.ac.il

## Abstract

### Background

Guide-design software commonly emphasizes sequence-level activity features and
predicted host off-targets. Viral target selection poses a broader decision problem:
candidate sites should also be evaluated across viral strains, mapped to gene and
protein context, distinguished from gene-level targetability, linked to biological
evidence with species provenance, and considered as multiplex panels with explicit
sequence-level escape assumptions.

### Results

We developed ViralSafeTarget, an open and reproducible framework that keeps these
decision axes separate. In an exhaustive herpes simplex virus 2 (HSV-2) case study,
28,578 initial candidate coordinates yielded 23,108 eligible rows and 21,654 unique
guide sequences. A completed 109-batch Cas-OFFinder search against GRCh38.p14 retained
440,341 predicted human-match rows; 2,668 candidate-coordinate rows had no predicted
human hit within the declared assembly, PAM and mismatch model. UL3, UL10 and UL52 led
the gene-targetability ranking, whereas the highest-ranked individual guide mapped to
UL36, whose gene-level portfolio ranked 21st. A frozen 257-guide panel produced 271
guide-to-CDS mappings, 5,691 bounded indel hypotheses and 17,733 single-nucleotide
exact-target counterfactuals. Four configured three-guide panels each required at least
three distinct substitutions to remove every exact guide/PAM target; this is a
sequence-level barrier, not an evolutionary probability. On the same frozen panel,
Cas-OFFinder and CRISPRitz host-search ranks showed Spearman correlation 0.880 and
top-50 overlap 49/50. ViralSafeTarget pre- and post-host ranks correlated at 0.962.

### Conclusions

ViralSafeTarget provides auditable decision support from viral genomes to research
shortlists while preserving missing data and evidence provenance. The HSV-2 analysis
shows that top individual guides, targetable genes and biologically established targets
are distinct concepts. Results are computational hypotheses for independent review and
do not establish editing, host safety, viral inhibition, delivery, treatment or cure.

**Keywords:** CRISPR; viral genomics; HSV-2; guide RNA; off-target screening;
multi-strain conservation; evidence provenance; reproducible bioinformatics

## Background

CRISPR-associated nucleases are programmable by short guide sequences, making target
selection a central computational step in genome editing [1]. Mature tools such as
CRISPOR and CHOPCHOP support guide design, activity scoring and predicted off-target
assessment [2,3], while Cas-OFFinder, CRISPRitz and GuideScan2 emphasize scalable
specificity or off-target analyses [4–6]. These capabilities are complementary rather
than interchangeable.

For a virus, a technically attractive guide may be absent from circulating strains,
may fall in a gene of uncertain relevance, or may produce limited predicted coding
disruption. Conversely, a classically important viral gene may offer fewer conserved
and host-distinct sites. A single composite “therapeutic” score would obscure these
differences and overstate what sequence analysis can establish.

ViralSafeTarget was designed as a virus-first framework that integrates existing
alignment and host-search engines with reproducible viral-population analysis,
one-to-many annotation mapping, separate guide- and gene-level rankings, bounded
coding-disruption hypotheses, source-provenanced literature evidence, and configured
multiplex comparisons. Its intended output is an auditable research shortlist, not a
clinical recommendation.

We report the software architecture and an exhaustive HSV-2 computational case study.
The case study tests three methodological propositions: (i) guide rank and gene-level
targetability can diverge; (ii) model-bounded host-search agreement is an axis distinct
from virus-first composite ranking; and (iii) exact-target robustness of a multiplex
panel can be represented reproducibly without claiming evolutionary probability.

## Implementation

### Software architecture and project contract

ViralSafeTarget 0.10.0 is implemented in Python and distributed as the
`viral-safe-target` package with the `vst` command-line interface. A project is defined
by versioned virus, host and nuclease profiles plus a project configuration. The same
public commands support initialization, validation, planning, execution, resumption,
status reporting, opening a primary HTML report, and exporting a portable result bundle.

Each stage records input checksums, parameters, tool versions, explicit missing or
external-required states, measured timings and output paths. Cached stages are reused
only when their signatures remain valid. Missing external results are never converted
to zero predicted hits.

### Viral inputs, quality control and conservation

The HSV-2 reference was RefSeq accession NC_001798.2. Discovery used a frozen set of 14
accepted HSV-2 genomes after explicit quality-control rejection of incomplete,
length-incompatible or ambiguity-heavy records. Viral sequences were aligned with
MAFFT [7]. Candidate support was measured as exact protospacer-plus-PAM presence across
accepted genomes; a separate held-out sequence collection was retained for population
support analysis when loci were observable.

### Candidate enumeration and annotation mapping

The reference genome was scanned for configured SpCas9 targets consisting of a
20-nucleotide protospacer adjacent to an NGG PAM. Sequence features included GC balance,
complexity, viral occurrence counts, exact strain coverage and annotation context.
Coordinates were mapped to GFF3 gene and CDS features. Overlapping annotations were
preserved as one-to-many mappings rather than collapsed to a single label. CDS strand
and phase were used to map cut coordinates to coding and protein positions.

### Host-genome screening

Eligible candidate sequences were screened against the GRCh38.p14 reference assembly
(GCF_000001405.40) with Cas-OFFinder 2.4.1 [4], using the configured NGG PAM and up to
three mismatches. The exhaustive workflow was divided into 109 checksum-resumable
batches. A completed search with zero returned sites was distinguished from a pending,
failed or missing batch. “Zero predicted hit” therefore always refers to this bounded
reference-assembly search and is not evidence of biological safety.

### Guide and gene targetability

Guide-level ranks retained conservation, viral uniqueness, GC, sequence complexity,
annotation and model-bounded host-search components with transparent explanations.
Gene-level targetability was computed separately from the portfolio of eligible and
screened candidates per annotated gene. It included the best candidate, robustness of
leading candidates, the lower confidence bound on the clean fraction, host-hit burden
and conservation. Biological evidence was reported separately and was not added to the
exhaustive gene-targetability score.

### Virtual coding-disruption hypotheses

For each guide-to-CDS mapping, integer indel sizes from −10 through +10 bp were
enumerated. Indel class, reading-frame remainder, determinable premature-stop position,
retained protein fraction and overlap with optional domain or disorder annotations were
reported. Insertions were represented by size only when their sequence was unknown.
Fractions over the equally weighted size grid are descriptive hypotheses and are not
repair probabilities.

### Exact-target escape counterfactuals and multiplex panels

All three alternative bases at every protospacer and PAM position were enumerated for
the frozen panel. A substitution was classified by whether it removed the exact
protospacer/PAM target under the configured IUPAC PAM pattern. For each configured
multiplex panel, an exact set-cover calculation identified the minimum number of
distinct substitutions required to remove all exact sites. This quantity omits mutation
rates, mismatch-tolerant cutting, viral fitness, selection, linkage and population
dynamics; it is not an evolutionary escape probability.

### Evidence Agent and mandatory human review

The Evidence Agent generates gene- and virus-aware queries, retrieves structured
records from PubMed, Europe PMC, UniProt and NCBI, and proposes evidence rows with
tested virus, experiment, phenotype, source and directness. AI-generated proposals do
not affect biological scores until approved by a named domain-qualified reviewer.
Evidence from HSV-1 orthologs remains separate from direct HSV-2 evidence. Biological
interpretations marked with an asterisk in this working manuscript remain pending
human verification.*

### Systematic multi-tool benchmark

A 257-guide panel was frozen before comparison. ViralSafeTarget pre-host and post-host
ranks, Cas-OFFinder hit-burden ranks, and CRISPRitz 2.6.6 profile ranks were normalized
within tool. The independent CRISPRitz run used its official container image digest,
GRCh38.p14, NGG, up to three mismatches, eight threads, and no bulge or
population-variant options [5]. Raw scores from different systems were not averaged.
CRISPOR, CHOPCHOP and GuideScan2 were included in a source-based capability matrix but
remained export-required in the quantitative benchmark because no raw export was
committed.

Rank agreement used Spearman correlation on shared candidates. Top-k overlap was
reported for k = 10, 25 and 50. A leave-one-component-out analysis recomputed the
transparent ViralSafeTarget score after omitting each configured component; it is a
sensitivity diagnostic, not model training or biological validation.

### Reproducibility, portability and generative-AI disclosure

The installed wheel was tested outside the repository checkout. A reference-only BK
polyomavirus workflow (NC_001538.1) used the same CLI and schemas without a
virus-specific branch in core code. Unavailable host and population stages remained
explicitly unavailable. This demonstrates software portability, not biological
generalization.

Generative AI assisted with code organization, documentation, literature-query
preparation and manuscript language. All numerical outputs reported here were generated
by versioned software from committed machine-readable sources and are checked by build
assertions. AI is not an author. Biological claims and source interpretations marked
with an asterisk remain pending author or domain-expert verification.*

## Results

### Exhaustive HSV-2 source validation

The publication build revalidated every primary count before generating tables or
figures (Table 1; Fig. 2). From 28,578 initial candidate coordinates, 23,108 were
eligible after sequence and annotation filters. These corresponded to 21,654 unique
guide sequences. All 109 Cas-OFFinder batches completed. The committed host-search
table contained 440,341 predicted human-match rows, and 2,668 candidate-coordinate rows
had no predicted match under the configured model.

### Guide rank and gene targetability diverged

UL3, UL10, UL52, UL47 and UL11 occupied the first five exhaustive gene-targetability
ranks (Table 2; Fig. 3). The highest-ranked individual candidate,
VST-2e9f052157f9bf29, mapped to UL36. UL36 ranked 21st at gene level, whereas UL3 led
the gene portfolio despite not containing the top individual guide. Thus the best
single candidate and the strongest gene-level portfolio were empirically distinct.

This result does not imply that UL3 is biologically or therapeutically preferable to
UL30 or UL52. Direct HSV-2 essentiality evidence for the focus genes remains
incomplete.* HSV-1 functional studies provide replication-focused context for UL30 and
UL52, but those observations remain ortholog evidence rather than direct HSV-2 null
evidence.*

### The deep panel added coding and sequence-robustness context

The frozen deep panel contained 257 unique guides. Generic mapping produced 271
guide-to-CDS rows because overlapping CDS records were retained; 250 guides mapped to
at least one CDS (Fig. 4). The −10 to +10 bp grid generated 5,691 bounded indel
hypotheses. The analysis enumerated 17,733 single-nucleotide protospacer/PAM
counterfactuals. Held-out exact-target coverage was available for 200 guides and
remained unknown for 57.

For the configured focus genes, all three analyzed UL18 guides intersected the supplied
domain annotation and occurred relatively early in the protein, giving UL18 a strong
predicted-disruption context.* This is annotation-dependent and does not demonstrate an
HSV-2 phenotype. UL36 contributed seven deep-panel guides, including the leading
individual candidate, but its gene-level portfolio remained lower ranked.

### Multiplex panels had a three-substitution exact-target barrier

Four three-guide strategies were evaluated: top-ranking-only,
essential/replication-focused, targetability-focused and mechanism-diverse (Fig. 5).
Each required three distinct substitutions to remove all exact guide/PAM sites under
the set-cover model. Minimum available held-out support ranged from 0.973913 to
0.982609. These are marginal per-guide values and do not establish within-genome joint
coverage or reduced biological escape.

### Host-search engines agreed while composite ranks addressed a different axis

Cas-OFFinder and CRISPRitz ranks correlated at 0.880413 over all 257 frozen guides
(Fig. 6). Their top-k overlaps were 9/10, 24/25 and 49/50. ViralSafeTarget pre-host and
post-host ranks correlated at 0.961975 and retained identical membership at top 10, 25
and 50. In contrast, host-burden ranks showed low correlation with ViralSafeTarget
composite ranks (0.085008–0.296323), as expected because the latter also incorporated
viral conservation and annotation context.

The leave-one-component-out analysis showed greatest sensitivity to sequence
complexity within this enriched panel: omitting it caused median absolute rank shift 75,
maximum shift 230 and top-10 overlap 0 (Supplementary Fig. S1). Conservation, viral
uniqueness, annotation and gene-evidence columns were constant in the selected panel,
so their omission produced no rank change. This diagnoses the frozen panel and scoring
configuration, not predictive correctness.

### The installed workflow generalized operationally to a second virus

The 0.10.0 wheel completed initialization, planning, execution, resumption, reporting
and export outside the source checkout. The same installed interface retrieved and
processed BK polyomavirus NC_001538.1 without core-code modification. It estimated 543
editor-compatible sites before a configured 500-row cap and completed annotation,
virtual analyses, multiplex comparison and export. Host screening and population
conservation remained unavailable, as no host assembly or multi-strain panel was
supplied.

## Discussion

ViralSafeTarget addresses a decision problem that extends beyond guide generation. Its
principal contribution is the auditable integration—and deliberate separation—of viral
population support, host-search status, guide rank, gene-level targetability, predicted
coding disruption, evidence provenance and exact-target multiplex robustness.

The HSV-2 case study illustrates why these axes should not be collapsed. UL36 contained
the leading individual guide but provided a weaker gene portfolio. UL3 led gene
targetability, while its direct HSV-2 biological importance remains unresolved.* UL30
and UL52 retained clearer HSV-1 replication-focused rationales*, yet did not occupy the
first two targetability ranks. These differences define testable research questions;
they do not establish therapeutic priority.

The independent host-search comparison increases confidence that the reported
reference-genome burden is not unique to one search engine under matched settings.
However, agreement cannot establish safety, and divergence from the composite rank is
not evidence that ViralSafeTarget is more accurate. The systems answer overlapping but
different questions.

To our knowledge, we did not identify an existing platform that jointly integrates
multi-strain viral conservation, model-bounded host-genome risk, separate guide- and
gene-level targetability, source-provenanced biological evidence, bounded coding
disruption, and configured exact-target multiplex escape analysis in one reproducible
project workflow. Individual capabilities are available in established tools, and the
claim concerns integration rather than invention of each component.

The intended use is early-stage computational prioritization and transparent handoff to
domain experts. Future work should complete human evidence review, obtain quantitative
exports from additional guide-design systems, evaluate a fully multi-strain second
virus, recruit external users, and compare computational ranks with experimental ground
truth.

## Conclusions

ViralSafeTarget converts viral reference, population, annotation and host inputs into
an auditable research shortlist with explicit uncertainty and provenance. Its HSV-2
case study demonstrates that guide quality, gene targetability, host-search burden,
predicted disruption, evidence coverage and exact-target robustness are distinct
properties. The software and frozen computational outputs support independent review
and reproducibility, but not claims of editing, safety, viral inhibition, treatment or
cure.

## Availability and requirements

- **Project name:** ViralSafeTarget
- **Project home page:** https://github.com/ronnuriel/ViralSafeTarget
- **Archived version:** [DOI PENDING]
- **Operating systems:** Platform-independent Python; external engines may have
  platform-specific requirements
- **Programming language:** Python 3.10 or later
- **Other requirements:** Optional MAFFT, Cas-OFFinder and CRISPRitz for corresponding
  stages; missing tools remain explicit
- **License:** MIT
- **Restrictions for non-academic use:** None under the MIT license

## List of abbreviations

CDS: coding sequence; CI: confidence interval; CRISPR: clustered regularly interspaced
short palindromic repeats; GFF3: General Feature Format version 3; HSV-1: herpes simplex
virus 1; HSV-2: herpes simplex virus 2; PAM: protospacer-adjacent motif; SNV:
single-nucleotide variant; VST: ViralSafeTarget.

## Declarations

### Ethics approval and consent to participate

Not applicable. This computational study used publicly available viral and reference
genome data and did not involve human participants, human tissue or animals.

### Consent for publication

Not applicable.

### Availability of data and materials

Source code and compact result snapshots are available at
https://github.com/ronnuriel/ViralSafeTarget. The archived software and supporting data
release will be cited here after DOI minting: [DOI PENDING]. Public sequence accessions,
configuration files, checksums and provenance manifests are included in the repository
and additional files.

### Competing interests

The authors declare no competing interests. [AUTHOR CONFIRMATION REQUIRED]

### Funding

No external funding was reported for this work. [AUTHOR CONFIRMATION REQUIRED]

### Authors' contributions

RN conceived the project, implemented and evaluated the software, interpreted the
computational findings and prepared the manuscript. SC's contribution statement is
[SAREL COHEN CONTRIBUTION REQUIRED]. Both authors reviewed and approved the final
manuscript. [AUTHOR CONFIRMATION REQUIRED]

### Acknowledgements

Generative AI systems assisted with code organization, documentation, query preparation
and language editing under author supervision. AI systems are not authors. Biological
claims marked with an asterisk remain pending human verification. [AUTHOR CONFIRMATION
REQUIRED]

### Authors' information

[OPTIONAL AUTHOR INFORMATION]

## References

1. Jinek M, Chylinski K, Fonfara I, Hauer M, Doudna JA, Charpentier E. A programmable dual-RNA-guided DNA endonuclease in adaptive bacterial immunity. Science. 2012;337:816–821. doi:10.1126/science.1225829.
2. Concordet JP, Haeussler M. CRISPOR: intuitive guide selection for CRISPR/Cas9 genome editing experiments and screens. Nucleic Acids Res. 2018;46:W242–W245. doi:10.1093/nar/gky354.
3. Labun K, Montague TG, Krause M, Torres Cleuren YN, Tjeldnes H, Valen E. CHOPCHOP v3: expanding the CRISPR web toolbox beyond genome editing. Nucleic Acids Res. 2019;47:W171–W174. doi:10.1093/nar/gkz365.
4. Bae S, Park J, Kim JS. Cas-OFFinder: a fast and versatile algorithm that searches for potential off-target sites of Cas9 RNA-guided endonucleases. Bioinformatics. 2014;30:1473–1475. doi:10.1093/bioinformatics/btu048.
5. Cancellieri S, Canver MC, Bombieri N, Giugno R, Pinello L. CRISPRitz: rapid, high-throughput and variant-aware in silico off-target site identification for CRISPR genome editing. Bioinformatics. 2020;36:2001–2008. doi:10.1093/bioinformatics/btz867.
6. Schmidt H, Zhang M, Mourelatos H, Sánchez-Rivera FJ, Lowe SW, Ventura A, et al. GuideScan2 for improved CRISPR guide RNA design. Nat Biotechnol. 2022;40:789–790. doi:10.1038/s41587-022-01278-8.
7. Katoh K, Standley DM. MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Mol Biol Evol. 2013;30:772–780. doi:10.1093/molbev/mst010.
8. National Center for Biotechnology Information. Human alphaherpesvirus 2 reference genome NC_001798.2. NCBI Nucleotide. https://www.ncbi.nlm.nih.gov/nuccore/NC_001798.2. Accessed 16 Jul 2026.
9. National Center for Biotechnology Information. GRCh38.p14 assembly GCF_000001405.40. NCBI Assembly. https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/. Accessed 16 Jul 2026.
10. Springer Nature. BMC Bioinformatics Software article submission guidelines. https://link.springer.com/journal/12859/submission-guidelines/software-article. Accessed 16 Jul 2026.
11. Markovitz NS, Filatov F, Roizman B. The UL3 protein of herpes simplex virus 1 is translated predominantly from the second in-frame methionine codon and is subject to at least two posttranslational modifications. J Virol. 1999;73:8010–8018. doi:10.1128/JVI.73.10.8010-8018.1999. PMID:10482549.*
12. MacLean CA, Robertson LM, Jamieson FE. Characterization of the UL10 gene product of herpes simplex virus type 1 and investigation of its role in vivo. J Gen Virol. 1993;74:975–983. doi:10.1099/0022-1317-74-6-975. PMID:8389812.*
13. Fulmer PA, Melancon JM, Baines JD, Kousoulas KG. UL20 protein functions precede and are required for the UL11 functions of herpes simplex virus type 1 cytoplasmic virion envelopment. J Virol. 2007;81:3097–3108. doi:10.1128/JVI.02201-06. PMID:17215291.*
14. Newcomb WW, Trus BL, Cheng N, Steven AC, Sheaffer AK, Tenney DJ, et al. Isolation of herpes simplex virus procapsids from cells infected with a protease-deficient mutant virus. J Virol. 2000;74:1663–1673. PMID:10644336.*
15. Desai PJ. A null mutation in the UL36 gene of herpes simplex virus type 1 results in accumulation of unenveloped DNA-filled capsids in the cytoplasm of infected cells. J Virol. 2000;74:11608–11618. doi:10.1128/JVI.74.24.11608-11618.2000. PMID:11090159.*
16. Klinedinst DK, Challberg MD. Helicase-primase complex of herpes simplex virus type 1: a mutation in the UL52 subunit abolishes primase activity. J Virol. 1994;68:3693–3701. doi:10.1128/JVI.68.6.3693-3701.1994. PMID:8189507.*
17. Stow ND. Sequences at the C-terminus of the herpes simplex virus type 1 UL30 protein are dispensable for DNA polymerase activity but not for viral origin-dependent DNA replication. Nucleic Acids Res. 1993;21:87–92. doi:10.1093/nar/21.1.87. PMID:8382792.*
18. Jayachandra S, Baghian A, Kousoulas KG. Herpes simplex virus type 1 glycoprotein K is not essential for infectious virus production in actively replicating cells but is required for efficient envelopment and translocation of infectious virions from the cytoplasm to the extracellular space. J Virol. 1997;71:5012–5024. doi:10.1128/JVI.71.7.5012-5024.1997. PMID:9188566.*
19. Jin F, Li S, Zheng K, Zhuo C, Ma K, Chen M, et al. Silencing herpes simplex virus type 1 capsid protein encoding genes by siRNA: a promising antiviral therapeutic approach. PLoS One. 2014;9:e96623. doi:10.1371/journal.pone.0096623. PMID:24794394.*
20. ViralSafeTarget contributors. ViralSafeTarget source repository. https://github.com/ronnuriel/ViralSafeTarget. Accessed 16 Jul 2026.

## Figure legends

**Figure 1. ViralSafeTarget separates complementary viral target decision layers.**
Versioned inputs pass through viral conservation, annotation, host screening, separate
guide and gene ranking, coding-disruption hypotheses, exact-target robustness and
human-reviewed evidence. Outputs remain research candidates rather than treatment
recommendations.

**Figure 2. Exhaustive HSV-2 screen and bounded host-search funnel.**
Counts are verified from committed source tables. The 2,668 zero-hit rows had no
predicted GRCh38.p14 match through three mismatches under the configured SpCas9 model;
this is not evidence of safety.

**Figure 3. Individual-guide rank diverges from gene targetability rank.**
Gene targetability ranks use candidate portfolios. The leading individual guide maps to
UL36, while UL36 ranks 21st as a gene and UL3 leads the gene-level table.

**Figure 4. Frozen-panel coding mappings preserve overlapping annotations.**
The 257-guide panel produced 271 guide-to-CDS mappings; 250 guides mapped to at least
one CDS. Values describe annotation coverage, not editing outcomes.

**Figure 5. Configured panels show exact-target sequence barriers.**
All four three-guide panels required three distinct substitutions to remove every exact
guide/PAM target. Held-out values are marginal per-guide observations. Barriers are not
evolutionary probabilities.

**Figure 6. Host-search rankings agree more closely than composite rankings.**
Cas-OFFinder and CRISPRitz were compared on the same 257 guides. ViralSafeTarget ranks
include viral and annotation features and therefore address a broader prioritization
axis; divergence is not evidence of superiority.

## Additional files

- **Additional file 1:** Supplementary Figure S1 (PDF). Leave-one-component-out rank
  sensitivity on the frozen deep panel.
- **Additional file 2:** Machine-readable source tables (ZIP). Tables underlying the
  reported counts, rankings, virtual analyses and benchmark figures.
- **Additional file 3:** ViralSafeTarget 0.10.0 source archive (ZIP) with SHA-256
  checksum, frozen from commit `de7868a83d3d1e30729323e53a735615d43fc231`.
