#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from viral_safe_target import (
    rank_candidates,
    read_cas_offinder_output,
    summarize_cas_offinder_hits,
    write_html_report,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Merge Cas-OFFinder output into a candidate report"
    )
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--cas-output", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--max-mismatches", type=int, default=3)
    args = parser.parse_args()

    candidates = pd.read_csv(args.candidates)
    hits = read_cas_offinder_output(args.cas_output)
    merged = summarize_cas_offinder_hits(
        candidates,
        hits,
        max_mismatches=args.max_mismatches,
    )
    ranked = rank_candidates(merged)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(args.out_dir / "candidates_with_human_screen.csv", index=False)
    write_html_report(
        ranked,
        args.out_dir / "report_with_human_screen.html",
        title="ViralSafeTarget — candidates after human off-target screen",
    )
    print(f"Saved {len(ranked)} candidates to {args.out_dir}")


if __name__ == "__main__":
    main()
