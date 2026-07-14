# Benchmarks

`published_studies.csv` is a literature index, not a ground-truth guide dataset. It intentionally omits exact editor sequences until each sequence is independently verified against the primary paper and reference genome.

A future machine-readable benchmark should distinguish:

- tested and active targets;
- tested and inactive targets;
- untested targets;
- lytic versus latent models;
- editor and delivery system;
- assay endpoint;
- host species and cell type;
- exact viral strain and reference coordinates.

Do not treat “targeted in a paper” as a positive biological label without reading the experimental result.

Use `known_targets.template.csv` as the empty curation contract and
`known_targets.schema.csv` as the machine-readable header. Run `vst benchmark` only
after exact sequences, PAMs, references, systems, and expected statuses have been
verified. See `docs/BENCHMARKING.md`.
