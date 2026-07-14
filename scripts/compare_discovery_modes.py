#!/usr/bin/env python3
"""Compare balanced and exhaustive completed genome-wide discovery outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from viral_safe_target.discovery_comparison import (
    compare_discovery_modes,
    write_discovery_comparison,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--balanced-dir", type=Path, required=True)
    parser.add_argument("--exhaustive-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    candidates, genes, summary = compare_discovery_modes(
        pd.read_csv(args.balanced_dir / "genome_wide_candidates_post_human.csv"),
        pd.read_csv(args.exhaustive_dir / "genome_wide_candidates_post_human.csv"),
        pd.read_csv(args.balanced_dir / "gene_rankings.csv"),
        pd.read_csv(args.exhaustive_dir / "gene_rankings.csv"),
    )
    write_discovery_comparison(candidates, genes, summary, args.out_dir)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
