#!/usr/bin/env python3
"""Resume previously prepared Cas-OFFinder batches without regenerating inputs."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from viral_safe_target.discovery_workflow import run_batches


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("batch_root", type=Path)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--device", choices=("C", "G", "A"), default="C")
    return parser


def main() -> None:
    args = _parser().parse_args()
    batches: list[dict[str, object]] = []
    for batch_dir in sorted(args.batch_root.glob("batch_*")):
        status_path = batch_dir / "manifest.json"
        if not status_path.is_file():
            continue
        record = json.loads(status_path.read_text(encoding="utf-8"))
        batches.append(
            {
                "batch_number": int(record["batch_number"]),
                "batch_dir": batch_dir,
                "input_path": batch_dir / "input.txt",
                "manifest_path": batch_dir / "candidate_manifest.csv",
                "raw_path": batch_dir / "raw_output.tsv",
                "status_path": status_path,
                "input_sha256": str(record["input_sha256"]),
                "candidate_manifest_sha256": str(record["candidate_manifest_sha256"]),
                "unique_query_count": int(record["unique_query_count"]),
                "candidate_count": int(record["candidate_count"]),
                "status": str(record["status"]),
                "previous": record,
            }
        )
    if not batches:
        raise SystemExit(f"No prepared batches found under {args.batch_root}")
    for batch in batches:
        if batch["status"] != "completed":
            Path(batch["raw_path"]).unlink(missing_ok=True)
    os.environ["CAS_OFFINDER_DEVICE"] = args.device
    run_batches(batches, args.executable.resolve(), execute=True)


if __name__ == "__main__":
    main()
