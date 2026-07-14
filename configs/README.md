# Configuration guide

Configuration is split between reusable profiles and versioned case-study settings.

## Reusable profiles

- `viruses/hsv2.yaml` — virus identity, reference/annotation paths, aliases, and
  evidence sources.
- `hosts/human_grch38.yaml` — host assembly identity and reference location.
- `nucleases/spcas9.yaml` — editor/PAM/cut model.
- `projects/hsv2_case_study.yaml` — binds virus, host, nuclease, ranking, output,
  virtual-knockout bounds, source assertions, and multiplex strategies into one
  researcher-facing project.

Add a new virus by creating a project with `vst project init`; do not add gene-name
conditionals to Python code.

## Case-study and compatibility settings

- `hsv2_pilot.yaml`
- `hsv2_consensus.yaml`
- `hsv2_genome_wide.yaml`
- `research_v0.3.yaml`

These preserve historical HSV-2 analyses and should not be copied as the primary way
to start a new virus project.

`examples/research_v0_3.example.yaml` documents the older monolithic configuration
contract for compatibility and migration; new work should use the profile bundle.

## External import maps

`import_maps/` contains column mappings for CRISPOR, CHOPCHOP, and GuideScan2 exports.
Missing imported fields remain missing rather than receiving synthetic values.

Machine-readable contracts are under [`../schemas/`](../schemas/README.md).
