# Virtual knockout and exact-target escape workflow

This workflow extends a ViralSafeTarget project after candidate discovery. It is
generic: core code reads the project virus, nuclease, reference, annotation, optional
protein-region tables, population table, and configured strategies. No viral gene name
is embedded in the algorithms.

## Commands

```bash
vst analyze virtual-knockout --project projects/my-virus/project.yaml
vst analyze escape --project projects/my-virus/project.yaml
vst analyze multiplex --project projects/my-virus/project.yaml
```

`vst project run` and `vst project resume` execute the same stages with
checksum/configuration-aware caching. `vst project status` reports
`virtual_knockout`, `escape`, and `multiplex` separately.

## Virtual-knockout model

For every candidate, the configured editor determines a cut boundary. The boundary is
mapped to every CDS containing it, so overlapping genes and reverse strands are
preserved. The configured integer range, by default -10 through +10 bp, is enumerated
as an equally weighted hypothesis grid.

Deletions are anchored downstream in coding orientation. Insertion length is known,
but sequence is unspecified; therefore insertion-dependent premature stops and
retained protein fractions remain unknown. Frameshift classification follows indel
size modulo three. Domain, disorder, and conserved-region labels are used only when a
profile supplies compatible protein-coordinate tables.

The fraction of grid rows classified as frameshifts is a property of the configured
grid. It is not a repair probability.

## Escape model

Observed support and counterfactual robustness are separate:

1. Discovery and held-out columns report observed exact protospacer plus compatible
   PAM presence when source tables are available.
2. Each possible single-nucleotide substitution in the protospacer/PAM is classified
   under the configured exact-target and IUPAC PAM model.
3. A multiplex barrier is solved as an exact set-cover problem: the minimum number of
   distinct genomic substitutions whose union removes every exact panel target.

The barrier is not an evolutionary probability. It does not model mutation rate,
selection, fitness, population dynamics, editing, or joint target occupancy.

## Configuration

The project `analysis` section controls:

- candidate and held-out tables;
- output directory;
- indel bounds;
- focus genes used only for case-study summaries;
- source-count assertions;
- strategy definitions with optional genes or biological categories.

Virus-specific gene lists belong in project data/configuration, never in core code.
If no strategy is configured, the project creates a top-ranking-only panel.

## Outputs

- `guide_cds_mapping.csv`
- `indel_sequence_hypotheses.csv`
- `guide_virtual_knockout.csv`
- `gene_virtual_knockout.csv`
- `single_nt_escape_counterfactuals.csv`
- `guide_escape_robustness.csv`
- `multiplex_panel_members.csv`
- `multiplex_escape_robustness.csv`
- `strategy_comparison.csv`
- `virtual_knockout_escape_report.html`
- `FINDINGS.md`
- `run_manifest.json`
- `figures/`

Unknown annotation, held-out support, or biological evidence remains missing. The
workflow does not create a combined therapeutic score.

## Interpretation boundary

These outputs are computational sequence hypotheses. They do not establish repair
frequency, successful editing, viral inactivation, safety, efficacy, delivery,
treatment, or cure. Independent biological review and experimental validation are
required for phenotype claims.
