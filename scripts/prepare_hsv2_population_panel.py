#!/usr/bin/env python3
"""Prepare a reproducible length-bounded HSV-2 population panel from NCBI outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from viral_safe_target.io_utils import read_fasta, write_fasta
from viral_safe_target.population_validation import (
    qc_population_records,
    select_population_accessions,
)
from viral_safe_target.provenance import sha256_file


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary-jsonl", type=Path, required=True)
    parser.add_argument("--source-fasta", type=Path)
    parser.add_argument(
        "--exclude-fasta",
        type=Path,
        help="Optional discovery/training FASTA whose accessions are excluded from validation.",
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--tax-id", type=int, default=10310)
    parser.add_argument("--minimum-length", type=int, default=120_000)
    parser.add_argument("--maximum-length", type=int, default=170_000)
    parser.add_argument(
        "--maximum-n-fraction",
        type=float,
        default=0.01,
        help=(
            "Maximum fraction of all IUPAC-ambiguous bases. The legacy option name "
            "is retained for command-line compatibility."
        ),
    )
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    selected = select_population_accessions(
        args.summary_jsonl,
        tax_id=args.tax_id,
        minimum_length=args.minimum_length,
        maximum_length=args.maximum_length,
    )
    selected.to_csv(args.out_dir / "accession_metadata.csv", index=False)
    (args.out_dir / "accessions.txt").write_text(
        "\n".join(selected["accession"].astype(str)) + "\n", encoding="utf-8"
    )
    manifest: dict[str, object] = {
        "schema_version": "1.0",
        "summary_jsonl": str(args.summary_jsonl.resolve()),
        "summary_sha256": sha256_file(args.summary_jsonl),
        "tax_id": args.tax_id,
        "minimum_length": args.minimum_length,
        "maximum_length": args.maximum_length,
        "selected_accession_count": len(selected),
    }
    if args.source_fasta:
        records = read_fasta(args.source_fasta)
        excluded_accessions = set(read_fasta(args.exclude_fasta)) if args.exclude_fasta else set()
        accepted, audit = qc_population_records(
            records,
            selected,
            maximum_n_fraction=args.maximum_n_fraction,
            excluded_accessions=excluded_accessions,
        )
        audit.to_csv(args.out_dir / "population_qc.csv", index=False)
        write_fasta(accepted, args.out_dir / "population_unique.fasta")
        manifest.update(
            {
                "source_fasta": str(args.source_fasta.resolve()),
                "source_fasta_sha256": sha256_file(args.source_fasta),
                "maximum_n_fraction": args.maximum_n_fraction,
                "maximum_ambiguous_base_fraction": args.maximum_n_fraction,
                "exclude_fasta": str(args.exclude_fasta.resolve()) if args.exclude_fasta else None,
                "exclude_fasta_sha256": sha256_file(args.exclude_fasta)
                if args.exclude_fasta
                else None,
                "excluded_accession_count": len(excluded_accessions),
                "accepted_unique_count": len(accepted),
                "rejected_count": len(audit) - len(accepted),
            }
        )
    (args.out_dir / "population_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
