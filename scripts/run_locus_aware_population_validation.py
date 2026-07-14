#!/usr/bin/env python3
"""Run reference-aware exact-target validation on a held-out viral population panel."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from viral_safe_target.config import get_editor, load_config
from viral_safe_target.io_utils import read_fasta
from viral_safe_target.population_validation import (
    exact_guide_presence_by_accession,
    map_population_to_reference,
    summarize_locus_aware_population_validation,
)
from viral_safe_target.provenance import sha256_file


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--population-fasta", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=Path("configs/hsv2_pilot.yaml"))
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--minimum-mapq", type=int, default=20)
    parser.add_argument("--minimum-identity", type=float, default=0.9)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(args.candidates)
    records = read_fasta(args.population_fasta)
    editor = get_editor(load_config(args.config))
    presence = exact_guide_presence_by_accession(candidates, records, editor)
    alignments = map_population_to_reference(
        records,
        args.reference_fasta,
        minimum_mapq=args.minimum_mapq,
        minimum_identity=args.minimum_identity,
    )
    validation = summarize_locus_aware_population_validation(
        candidates, presence, alignments, list(records)
    )
    alignments.to_csv(args.out_dir / "population_reference_alignments.csv", index=False)
    validation.to_csv(args.out_dir / "candidate_locus_population_validation.csv", index=False)
    manifest = {
        "schema_version": "1.0",
        "population_record_count": len(records),
        "candidate_count": len(validation),
        "reference_alignment_count": len(alignments),
        "records_with_reference_alignment": int(alignments["accession"].nunique())
        if not alignments.empty
        else 0,
        "minimum_mapq": args.minimum_mapq,
        "minimum_identity": args.minimum_identity,
        "population_fasta": str(args.population_fasta.resolve()),
        "population_fasta_sha256": sha256_file(args.population_fasta),
        "reference_fasta": str(args.reference_fasta.resolve()),
        "reference_fasta_sha256": sha256_file(args.reference_fasta),
        "candidate_source": str(args.candidates.resolve()),
        "candidate_source_sha256": sha256_file(args.candidates),
        "interpretation": (
            "Reference-aware exact target validation only; no editing efficacy, safety, "
            "delivery, or therapeutic claim."
        ),
    }
    (args.out_dir / "locus_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
