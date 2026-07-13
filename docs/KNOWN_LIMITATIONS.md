# Known limitations

- The current core supports SpCas9/NGG only. Editor profiles should become configurable before claiming generality.
- Exact 23-nt coverage is intentionally strict and does not model mismatch tolerance at the viral target.
- The built-in host screen is a small-FASTA teaching implementation, not a GRCh38 engine.
- Whole HSV genomes contain repeats and alternative isomer orientations. Input genomes should be normalized and alignment quality inspected; a naive whole-genome alignment may create false conservation or false variation.
- GFF annotation indicates location and known feature names, not gene essentiality or target accessibility.
- The pair simulator assumes canonical cut coordinates and an idealized deletion. It does not predict repair distributions.
- The demo ranking is transparent but not experimentally calibrated.
- Population-representative sampling, human pangenome variation and clinically relevant diversity require dedicated study designs.
