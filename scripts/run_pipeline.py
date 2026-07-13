from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from viral_safe_target import (  # noqa: E402
    annotate_candidates,
    rank_candidates,
    read_fasta,
    read_gff3,
    scan_spcas9_candidates,
    screen_against_small_fasta,
    write_html_report,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the ViralSafeTarget demo pipeline")
    parser.add_argument("--virus-alignment", required=True)
    parser.add_argument("--reference-id", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--human-fasta", required=True, help="Small/demo FASTA only")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--min-coverage", type=float, default=0.0)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    virus = read_fasta(args.virus_alignment)
    candidates = scan_spcas9_candidates(virus, args.reference_id, args.min_coverage)
    gff = read_gff3(args.gff)
    candidates = annotate_candidates(candidates, gff, seqid=args.reference_id)
    human = read_fasta(args.human_fasta)
    candidates = screen_against_small_fasta(candidates, human)
    candidates = rank_candidates(candidates)

    csv_path = out_dir / "candidates.csv"
    html_path = out_dir / "report.html"
    candidates.to_csv(csv_path, index=False)
    write_html_report(candidates, html_path, title="ViralSafeTarget demo")
    print(f"Saved {len(candidates)} candidates to {csv_path}")
    print(f"Saved report to {html_path}")


if __name__ == "__main__":
    main()
