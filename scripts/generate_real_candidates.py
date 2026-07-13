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
    args = parser.parse_args()

    records = read_fasta(args.alignment)
    candidates = scan_spcas9_candidates(records, args.reference_id, args.min_coverage)
    features = read_gff3(args.gff)
    candidates = annotate_candidates(candidates, features, seqid=args.reference_id)

    # Pre-human prioritization: conservation dominates; annotation is only context.
    if not candidates.empty:
        candidates["annotation_component"] = (
            candidates["feature_type"] != "intergenic_or_unannotated"
        ).astype(float)
        candidates["pre_human_score"] = (
            0.85 * candidates["virus_site_coverage"].astype(float)
            + 0.15 * candidates["annotation_component"]
        )
        candidates["decision"] = "human_off_target_screen_pending"
        candidates = candidates.sort_values(
            ["pre_human_score", "virus_site_coverage"], ascending=False
        ).reset_index(drop=True)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / "candidates_pre_human.csv"
    html_path = args.out_dir / "report_pre_human.html"
    candidates.to_csv(csv_path, index=False)

    report_df = candidates.rename(columns={"pre_human_score": "demo_score"})
    write_html_report(
        report_df,
        html_path,
        title="ViralSafeTarget — HSV-2 real-data pilot (before human screening)",
    )

    print(f"Generated {len(candidates)} candidates")
    print(f"CSV:  {csv_path}")
    print(f"HTML: {html_path}")

    if args.human_fasta_directory:
        cas_input = args.out_dir / "cas_offinder_input.txt"
        write_cas_offinder_input(
            candidates,
            human_fasta_directory=args.human_fasta_directory,
            output_path=cas_input,
            max_mismatches=args.max_mismatches,
        )
        print(f"Cas-OFFinder input: {cas_input}")


if __name__ == "__main__":
    main()
