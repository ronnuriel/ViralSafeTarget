"""Command-line interface for ViralSafeTarget."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .annotations import annotate_candidates, read_gff3
from .crispr import scan_spcas9_candidates
from .disruption import rank_candidate_pairs
from .io_utils import read_fasta
from .offtarget import screen_against_small_fasta
from .provenance import write_run_manifest
from .reporting import write_html_report
from .scoring import rank_candidates


def _scan(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    virus = read_fasta(args.virus_alignment)
    candidates = scan_spcas9_candidates(virus, args.reference_id, args.min_coverage)
    features = None
    if args.gff:
        features = read_gff3(args.gff)
        candidates = annotate_candidates(candidates, features, seqid=args.reference_id)
    if args.small_host_fasta:
        host = read_fasta(args.small_host_fasta)
        candidates = screen_against_small_fasta(
            candidates, host, max_mismatches=args.max_mismatches
        )
        candidates = rank_candidates(candidates)
    else:
        candidates["decision"] = "full_host_off_target_screen_pending"
    candidates.to_csv(out_dir / "candidates.csv", index=False)
    write_html_report(candidates, out_dir / "report.html", title="ViralSafeTarget scan")
    input_paths = [args.virus_alignment]
    if args.gff:
        input_paths.append(args.gff)
    if args.small_host_fasta:
        input_paths.append(args.small_host_fasta)
    write_run_manifest(
        out_dir / "run_manifest.json",
        input_paths,
        {
            "reference_id": args.reference_id,
            "min_coverage": args.min_coverage,
            "max_mismatches": args.max_mismatches,
        },
    )
    print(f"Wrote {len(candidates)} candidates to {out_dir}")


def _simulate_pairs(args: argparse.Namespace) -> None:
    candidates = pd.read_csv(args.candidates)
    features = read_gff3(args.gff) if args.gff else None
    alignment = read_fasta(args.virus_alignment) if args.virus_alignment else None
    pairs = rank_candidate_pairs(
        candidates,
        features=features,
        aligned_records=alignment,
        reference_id=args.reference_id,
        same_feature_only=not args.allow_cross_feature,
        min_distance_bp=args.min_distance,
        max_distance_bp=args.max_distance,
        max_candidates=args.max_candidates,
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(output, index=False)
    print(f"Wrote {len(pairs)} idealized pair simulations to {output}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="viral-safe-target",
        description="Conserved viral target prioritization for computational research",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    scan = subparsers.add_parser("scan", help="scan an aligned viral genome collection")
    scan.add_argument("--virus-alignment", required=True)
    scan.add_argument("--reference-id", required=True)
    scan.add_argument("--gff")
    scan.add_argument("--small-host-fasta", help="demo only; guarded at 5 Mb")
    scan.add_argument("--out-dir", required=True)
    scan.add_argument("--min-coverage", type=float, default=0.95)
    scan.add_argument("--max-mismatches", type=int, default=3)
    scan.set_defaults(func=_scan)

    pairs = subparsers.add_parser(
        "simulate-pairs", help="simulate idealized sequence deletions between candidate pairs"
    )
    pairs.add_argument("--candidates", required=True)
    pairs.add_argument("--gff")
    pairs.add_argument("--virus-alignment")
    pairs.add_argument("--reference-id")
    pairs.add_argument("--output", required=True)
    pairs.add_argument("--min-distance", type=int, default=20)
    pairs.add_argument("--max-distance", type=int, default=10_000)
    pairs.add_argument("--max-candidates", type=int, default=250)
    pairs.add_argument("--allow-cross-feature", action="store_true")
    pairs.set_defaults(func=_simulate_pairs)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
