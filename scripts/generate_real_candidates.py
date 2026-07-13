#!/usr/bin/env python3
"""Generate conserved SpCas9 candidate sites from an aligned HSV-2 pilot dataset.

The output is a computational prioritization table only. It does not establish
biological efficacy or safety.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from viral_safe_target import (  # noqa: E402
    annotate_candidates,
    rank_pre_human_candidates,
    read_fasta,
    read_gff3,
    scan_spcas9_candidates,
    write_cas_offinder_input,
    write_html_report,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--gff", type=Path, required=True)
    parser.add_argument("--reference-id", default="NC_001798.2")
    parser.add_argument("--min-coverage", type=float, default=0.95)
    parser.add_argument("--out-dir", type=Path, default=Path("reports/real_hsv2"))
    parser.add_argument("--human-fasta-directory", type=Path)
    parser.add_argument("--max-mismatches", type=int, default=3)
    parser.add_argument("--config", type=Path, default=Path("configs/research_v0.3.yaml"))
    parser.add_argument("--gene-evidence", type=Path)
    args = parser.parse_args()

    records = read_fasta(args.alignment)
    candidates = scan_spcas9_candidates(records, args.reference_id, args.min_coverage)
    features = read_gff3(args.gff)
    candidates = annotate_candidates(candidates, features, seqid=args.reference_id)

    candidates = rank_pre_human_candidates(candidates, args.config, args.gene_evidence)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    retained = candidates[candidates["rejection_reasons"].fillna("").eq("")].copy()
    rejected = candidates[candidates["rejection_reasons"].fillna("").ne("")].copy()
    csv_path = args.out_dir / "candidates_ranked_pre_human.csv"
    legacy_csv_path = args.out_dir / "candidates_pre_human.csv"
    html_path = args.out_dir / "report.html"
    retained.to_csv(csv_path, index=False)
    retained.to_csv(legacy_csv_path, index=False)
    rejected.to_csv(args.out_dir / "candidates_rejected_pre_human.csv", index=False)

    write_html_report(
        retained,
        html_path,
        title="ViralSafeTarget — HSV-2 real-data pilot (before human screening)",
    )

    print(f"Generated {len(retained)} retained candidates; rejected {len(rejected)}")
    print(f"CSV:  {csv_path}")
    print(f"HTML: {html_path}")

    if args.human_fasta_directory:
        cas_input = args.out_dir / "cas_offinder_input.txt"
        write_cas_offinder_input(
            retained,
            human_fasta_directory=args.human_fasta_directory,
            output_path=cas_input,
            max_mismatches=args.max_mismatches,
        )
        print(f"Cas-OFFinder input: {cas_input}")


if __name__ == "__main__":
    main()
