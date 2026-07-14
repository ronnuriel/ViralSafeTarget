#!/usr/bin/env python3
"""Validate candidate guide/PAM presence against an independent viral FASTA panel."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from viral_safe_target.config import get_editor, load_config
from viral_safe_target.io_utils import read_fasta
from viral_safe_target.population_validation import candidate_population_validation
from viral_safe_target.provenance import sha256_file


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--population-fasta", type=Path, required=True)
    parser.add_argument(
        "--population-qc",
        type=Path,
        help="Optional QC table used to stratify exact coverage by completeness group.",
    )
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=Path("configs/hsv2_pilot.yaml"))
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(args.candidates)
    records = read_fasta(args.population_fasta)
    record_groups = None
    group_counts: dict[str, int] = {}
    if args.population_qc:
        qc = pd.read_csv(args.population_qc)
        qc = qc[qc["decision"].eq("accepted")].copy()
        record_groups = dict(
            zip(
                qc["accession"].astype(str),
                qc["submitter_completeness"].fillna("unknown").astype(str),
                strict=False,
            )
        )
        group_counts = pd.Series(record_groups).value_counts().sort_index().astype(int).to_dict()
    result = candidate_population_validation(
        candidates,
        records,
        get_editor(load_config(args.config)),
        record_groups=record_groups,
    )
    result.to_csv(args.out_dir / "candidate_population_validation.csv", index=False)
    summary = {
        "schema_version": "1.0",
        "candidate_count": len(result),
        "unique_guide_count": int(candidates["guide_sequence"].nunique()),
        "population_genome_count": len(records),
        "population_fasta": str(args.population_fasta.resolve()),
        "population_fasta_sha256": sha256_file(args.population_fasta),
        "candidate_source": str(args.candidates.resolve()),
        "candidate_source_sha256": sha256_file(args.candidates),
        "population_qc": str(args.population_qc.resolve()) if args.population_qc else None,
        "population_qc_sha256": sha256_file(args.population_qc) if args.population_qc else None,
        "population_group_counts": group_counts,
        "fully_conserved_candidate_count": int(
            result["population_exact_pam_compatible_coverage"].eq(1).sum()
        ),
        "interpretation": (
            "Exact sequence/PAM held-out validation only; no editing, efficacy, safety, "
            "or locus-completeness claim."
        ),
    }
    (args.out_dir / "population_validation_manifest.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
