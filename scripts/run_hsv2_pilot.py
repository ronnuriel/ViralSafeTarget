#!/usr/bin/env python3
"""Cross-platform HSV-2 pilot runner."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


def run(command: list[str]) -> None:
    print("+", " ".join(command))
    subprocess.run(command, check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL"))
    parser.add_argument("--max-genomes", type=int, default=50)
    parser.add_argument("--min-coverage", type=float, default=0.95)
    parser.add_argument("--skip-download", action="store_true")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    raw = root / "data/raw/hsv2_entrez"
    processed = root / "data/processed"
    reports = root / "reports/real_hsv2"
    processed.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)

    if not args.skip_download:
        if not args.email:
            raise SystemExit("Pass --email or set NCBI_EMAIL.")
        run(
            [
                sys.executable,
                str(root / "scripts/download_hsv2_pilot.py"),
                "--email",
                args.email,
                "--max-genomes",
                str(args.max_genomes),
                "--out-dir",
                str(raw),
            ]
        )

    required = [
        raw / "hsv2_reference.fasta",
        raw / "hsv2_reference.gb",
        raw / "hsv2_complete_genomes.fasta",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("Missing data files: " + ", ".join(missing))
    if shutil.which("mafft") is None:
        raise SystemExit(
            "MAFFT is not installed. Create the Conda environment from environment.yml."
        )

    sample_fasta = processed / "hsv2_sample.fasta"
    gff = processed / "hsv2_reference.gff3"
    alignment = processed / "hsv2_aligned.fasta"
    run(
        [
            sys.executable,
            str(root / "scripts/prepare_real_hsv2.py"),
            "--reference-fasta",
            str(raw / "hsv2_reference.fasta"),
            "--reference-genbank",
            str(raw / "hsv2_reference.gb"),
            "--all-genomes-fasta",
            str(raw / "hsv2_complete_genomes.fasta"),
            "--output-fasta",
            str(sample_fasta),
            "--output-gff",
            str(gff),
            "--sample-size",
            str(args.max_genomes),
        ]
    )
    print("+ mafft --auto --thread -1", sample_fasta, ">", alignment)
    with alignment.open("w", encoding="utf-8") as handle:
        subprocess.run(
            ["mafft", "--auto", "--thread", "-1", str(sample_fasta)],
            check=True,
            stdout=handle,
        )
    run(
        [
            sys.executable,
            str(root / "scripts/generate_real_candidates.py"),
            "--alignment",
            str(alignment),
            "--gff",
            str(gff),
            "--reference-id",
            "NC_001798.2",
            "--min-coverage",
            str(args.min_coverage),
            "--out-dir",
            str(reports),
        ]
    )
    print(f"Pilot complete. Open {reports / 'report_pre_human.html'}")


if __name__ == "__main__":
    main()
