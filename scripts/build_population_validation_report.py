#!/usr/bin/env python3
"""Build a held-out HSV-2 population-validation comparison report."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from viral_safe_target.population_reporting import (
    build_population_comparison,
    write_population_report,
)
from viral_safe_target.provenance import sha256_file


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--locus-validation", type=Path, required=True)
    parser.add_argument("--locus-manifest", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    candidates = pd.read_csv(args.candidates)
    validation = pd.read_csv(args.locus_validation)
    comparison, genes = build_population_comparison(candidates, validation)
    locus_manifest = json.loads(args.locus_manifest.read_text(encoding="utf-8"))
    provenance = {
        **locus_manifest,
        "candidate_source": str(args.candidates.resolve()),
        "candidate_source_sha256": sha256_file(args.candidates),
        "locus_validation_source": str(args.locus_validation.resolve()),
        "locus_validation_source_sha256": sha256_file(args.locus_validation),
        "score_integration": "none; population evidence is reported separately",
    }
    write_population_report(comparison, genes, args.out_dir, provenance)
    print(json.dumps(provenance, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
