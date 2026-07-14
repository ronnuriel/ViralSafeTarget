# Adding an editor

Editor profiles live in versioned YAML and require:

- `name`
- `protospacer_length`
- `pam_pattern` (IUPAC)
- `pam_orientation` (`3prime` or `5prime`)
- `cut_offset`
- `mismatch_search_threshold`
- optional `notes` and a truthful `tested` flag

SpCas9/NGG is the only scanner profile tested in v0.3. The configuration model can
represent additional editors, but a profile must not be marked supported until PAM
matching, both strands, boundary coordinates, Cas-OFFinder formatting, and synthetic
end-to-end behavior have dedicated tests. Thresholds describe a search model, not a
claim of equivalent biochemical behavior across editors.
