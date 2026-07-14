"""Auditable case-study reproduction entry points."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .provenance import sha256_file


def hsv2_reproduction_plan(
    project_root: str | Path = ".", *, skip_population: bool = False
) -> list[dict[str, Any]]:
    root = Path(project_root).resolve()
    stages = [
        {
            "stage": "environment",
            "command": "vst doctor && vst tools doctor",
            "purpose": "Record Python and external-tool availability.",
            "expected_output": "console audit",
        },
        {
            "stage": "public_data_and_candidates",
            "command": (
                "bash scripts/run_real_hsv2.sh --sample-size 25 --with-human "
                "--accessions-file data/curated/hsv2_discovery_accessions.txt"
            ),
            "purpose": "Fetch frozen public accessions, QC, align, scan, and prepare GRCh38.",
            "expected_output": "reports/real_hsv2/run_manifest.json",
        },
        {
            "stage": "ortholog_reference",
            "command": (
                f"{sys.executable} scripts/fetch_reference_genbank.py "
                "--accession NC_001806.2 --output data/raw/hsv1_reference.gb"
            ),
            "purpose": "Fetch the declared HSV-1 ortholog reference used by protein analysis.",
            "expected_output": "data/raw/hsv1_reference.gb",
        },
        {
            "stage": "genome_wide_host_screen",
            "command": "bash scripts/run_hsv2_genome_wide.sh",
            "purpose": "Run or resume the balanced Cas-OFFinder screen.",
            "expected_output": "reports/hsv2_genome_wide/report.html",
        },
    ]
    if not skip_population:
        stages.append(
            {
                "stage": "heldout_population",
                "command": "bash scripts/run_hsv2_population_validation.sh",
                "purpose": "Validate candidate observability in a discovery-excluded population.",
                "expected_output": (
                    "reports/hsv2_population_report_balanced/population_validation_report.html"
                ),
            }
        )
    stages.append(
        {
            "stage": "gene_function_and_showcase",
            "command": "bash scripts/build_hsv2_showcase.sh",
            "purpose": "Build protein-disruption analysis and presentation report.",
            "expected_output": "reports/hsv2_showcase/FINAL_REPORT.html",
        }
    )
    for stage in stages:
        output = root / str(stage["expected_output"])
        stage["current_status"] = "cached_output_present" if output.exists() else "pending"
    return stages


def _execute(command: list[str], root: Path, environment: dict[str, str]) -> None:
    completed = subprocess.run(command, cwd=root, env=environment, check=False)
    if completed.returncode:
        raise RuntimeError(
            f"Reproduction command failed with exit code {completed.returncode}: "
            + " ".join(command)
        )


def reproduce_hsv2(
    project_root: str | Path = ".",
    *,
    execute: bool = False,
    skip_population: bool = False,
) -> dict[str, Any]:
    root = Path(project_root).resolve()
    plan = hsv2_reproduction_plan(root, skip_population=skip_population)
    if not execute:
        return {"mode": "plan", "project_root": str(root), "stages": plan}

    required = ["bash", "mafft", "datasets", "unzip"]
    missing = [name for name in required if shutil.which(name) is None]
    cas = os.environ.get("CAS_OFFINDER_PATH") or shutil.which("cas-offinder")
    if not cas or not Path(cas).is_file():
        missing.append("cas-offinder (or CAS_OFFINDER_PATH)")
    if missing:
        raise RuntimeError("Missing reproduction dependencies: " + ", ".join(missing))
    environment = os.environ.copy()
    environment["CAS_OFFINDER_PATH"] = str(cas)
    stages_run: list[dict[str, str]] = []

    commands = [
        ("environment", [sys.executable, "-m", "viral_safe_target", "doctor"]),
        (
            "public_data_and_candidates",
            [
                "bash",
                "scripts/run_real_hsv2.sh",
                "--sample-size",
                "25",
                "--with-human",
                "--accessions-file",
                "data/curated/hsv2_discovery_accessions.txt",
            ],
        ),
        (
            "ortholog_reference",
            [
                sys.executable,
                "scripts/fetch_reference_genbank.py",
                "--accession",
                "NC_001806.2",
                "--output",
                "data/raw/hsv1_reference.gb",
            ],
        ),
        ("genome_wide_host_screen", ["bash", "scripts/run_hsv2_genome_wide.sh"]),
    ]
    if not skip_population:
        commands.append(
            ("heldout_population", ["bash", "scripts/run_hsv2_population_validation.sh"])
        )
    commands.append(("gene_function_and_showcase", ["bash", "scripts/build_hsv2_showcase.sh"]))
    for stage, command in commands:
        _execute(command, root, environment)
        stages_run.append({"stage": stage, "status": "completed", "command": " ".join(command)})

    report = root / "reports" / "hsv2_showcase" / "FINAL_REPORT.html"
    manifest_path = root / "reports" / "hsv2_reproduction" / "reproduction_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema_version": "1.0",
        "case_study": "hsv2",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": str(root),
        "stages": stages_run,
        "frozen_accession_file": str(
            (root / "data/curated/hsv2_discovery_accessions.txt").resolve()
        ),
        "frozen_accession_file_sha256": sha256_file(
            root / "data/curated/hsv2_discovery_accessions.txt"
        ),
        "final_report": str(report.resolve()),
        "final_report_sha256": sha256_file(report) if report.is_file() else None,
        "interpretation": (
            "Computational reproduction only; no editing, safety, efficacy, delivery, "
            "wet-lab, or therapeutic claim."
        ),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "mode": "executed",
        "manifest": str(manifest_path),
        "final_report": str(report),
        "stages": stages_run,
    }
