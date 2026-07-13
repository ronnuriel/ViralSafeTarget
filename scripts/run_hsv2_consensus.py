#!/usr/bin/env python3
from __future__ import annotations

import argparse

from viral_safe_target.hsv2_consensus import run_hsv2_consensus


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="configs/hsv2_consensus.yaml")
    parser.add_argument("--out-dir", default="reports/hsv2_consensus")
    args = parser.parse_args()
    run_hsv2_consensus(args.config, args.out_dir)


if __name__ == "__main__":
    main()
