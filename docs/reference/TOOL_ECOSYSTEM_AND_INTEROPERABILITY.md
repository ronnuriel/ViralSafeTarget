# Tool ecosystem and interoperability

Last audited: 2026-10-22

ViralSafeTarget is an orchestration and decision-support layer. It reuses established
alignment, guide-design, host-search, experimental-analysis, and evidence resources while
keeping every source, parameter, missing result, and interpretation boundary visible. It
does not claim that one combined score replaces the underlying tools or biological review.

## What researchers use and how ViralSafeTarget connects

| System | Established role | Exchange route | ViralSafeTarget status |
|---|---|---|---|
| MAFFT | Multiple-sequence alignment | FASTA in; aligned FASTA out | Native executable adapter |
| Cas-OFFinder | Fast mismatch-bounded genome search | Native input and tabular hit output | Native execution and parser |
| CRISPRitz | Mismatch, DNA/RNA bulge, annotation, and variant-aware analysis | Native/container bundle and result files | Native/container execution and import |
| CRISPOR | Guide design, activity/specificity scores, off-target context, and primer helpers | Researcher-downloaded CSV/TSV/JSON | Versioned import map; public site is not scraped |
| CHOPCHOP | Interactive multi-editor guide design, ranking, visualization, and downloadable results | Researcher-downloaded text/TSV/FASTA/GenBank output | Versioned import map; public site is not scraped |
| GuideScan2 | Genome-wide guide databases and specificity analysis | Downloaded table or CLI/database result | Versioned import map; CLI execution is not yet native |
| CRISPResso2 | Quantification and visualization of measured genome-editing outcomes | Existing CRISPResso2 result directory | Measured-result import, kept separate from predictions |

The normalized import schema records stable candidate identity, raw and normalized metric,
metric direction, tool/version/mode, genome, source hash, command/import method, timestamp,
status, and error. Incomparable scores are not averaged. Unmatched and ambiguous rows are
retained for inspection.

### Current commands

```bash
vst tools status --project my-virus/project.yaml
vst tools import --tool crispor --input crispor.tsv \
  --mapping configs/import_maps/crispor.yaml \
  --candidates my-virus/results/discovery/discovery_panel.csv \
  --out-dir my-virus/results/external/crispor
vst tools import --tool chopchop --input chopchop.tsv \
  --mapping configs/import_maps/chopchop.yaml \
  --candidates my-virus/results/discovery/discovery_panel.csv \
  --out-dir my-virus/results/external/chopchop
vst tools import --tool guidescan2 --input guidescan2.tsv \
  --mapping configs/import_maps/guidescan2.yaml \
  --candidates my-virus/results/discovery/discovery_panel.csv \
  --out-dir my-virus/results/external/guidescan2
vst experimental import-crispresso2 --input existing_crispresso_results \
  --candidate-map candidate_amplicon_map.csv \
  --out-dir my-virus/results/experimental/crispresso2
```

The same CRISPOR, CHOPCHOP, and GuideScan2 export validation is available in the local UI
under **Tool & evidence hub**.

## Evidence Agent: accepted sources and human control

The current agent uses documented APIs rather than scraping journal or database pages:

| Source | Agent use | Boundary |
|---|---|---|
| NCBI E-utilities / PubMed | Search and retrieve PMID-linked article metadata and abstracts | API limits and copyright rules apply |
| Europe PMC REST API | Retrieve life-science publication metadata, abstracts, OA links, and preprint labels | Preprints and peer-reviewed articles remain distinguishable |
| UniProt REST API | Retrieve organism-specific protein names, functional annotation, domains, and publication links | Curated and computational annotations must not be conflated |

The agent creates queries, retrieves records, deduplicates identifiers, hashes source
material, links a proposed claim to a source, and places it in a review queue. Only a
named human reviewer may approve, reject, or qualify a biological claim. No result means
`unknown`, not evidence that a gene is non-essential or irrelevant.

## Audited next connectors

These resources are useful for a future translational layer, but they are **not connected
in the current released workflow**:

| Candidate connector | Potential contribution | Required safeguard before implementation |
|---|---|---|
| InterPro / InterProScan | Independent protein-family, domain, and functional-site annotation | Pin database/version and preserve sequence-coordinate mapping |
| ChEMBL Web Services | Curated molecule, assay, target, and bioactivity records | Keep assay evidence separate from viral-genome target ranking |
| PubChem PUG REST | Compound, assay, gene/protein, structure, and identifier retrieval | Do not translate database association into efficacy |
| Open Targets Platform API | Target-disease evidence and tractability context | Viral proteins and host targets require separate semantics |

Adding these connectors could help researchers move from a viral target hypothesis toward
an auditable intervention-evidence map. It would not predict a drug, recommend treatment,
or replace medicinal chemistry, pharmacology, safety assessment, or experiments.

## Why this integration is useful

Each established tool answers a narrower question well. ViralSafeTarget keeps their
answers separate and brings them into one traceable workflow:

```text
viral population + annotation
        -> candidate discovery and conservation
        -> independent host-search/design tools
        -> normalized, provenance-preserving comparison
        -> sequence/protein disruption and exact-target escape hypotheses
        -> source-linked evidence proposals
        -> mandatory human review
        -> research shortlist and portable result bundle
```

The intended contribution is not “AI decides the cure.” It is a reproducible research
workspace that reduces manual stitching, exposes disagreements and missing data, and makes
the next biological questions easier to audit and hand to collaborators.

## Official sources

- [CRISPOR manual](https://github.com/maximilianh/crisporWebsite/blob/master/manual/manual.md)
- [CHOPCHOP instructions](https://chopchop.cbu.uib.no/instructions)
- [CHOPCHOP v2 paper and export description](https://doi.org/10.1093/nar/gkw398)
- [CRISPRitz official repository](https://github.com/pinellolab/CRISPRitz)
- [GuideScan2 web documentation](https://guidescan.com/about)
- [GuideScan2 CLI repository](https://github.com/pritykinlab/guidescan-cli)
- [CRISPResso2 output documentation](https://crispresso.pinellolab.org/help)
- [NCBI developer APIs and E-utilities](https://www.ncbi.nlm.nih.gov/home/develop/api/)
- [Europe PMC developer resources](https://europepmc.org/developers)
- [UniProt REST API documentation](https://www.uniprot.org/help/api)
- [ChEMBL Data Web Services](https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services)
- [PubChem PUG REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest)
- [Open Targets technical overview](https://opentargets.org/technical-knowledge-base/platform/technical-overview/)

