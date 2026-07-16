#!/usr/bin/env python3
"""Assemble the rendered BMC working package and DOI-preparation archives."""

from __future__ import annotations

import hashlib
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT / "paper/bmc_bioinformatics/final"
RELEASE = FINAL / "release"
PACKAGE = ROOT / "paper/bmc_bioinformatics/ViralSafeTarget_BMC_Submission_WORKING.zip"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_zip(output: Path, files: list[tuple[Path, str]]) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED) as archive:
        for source, destination in files:
            if not source.exists():
                raise FileNotFoundError(f"Required publication artifact is missing: {source}")
            archive.write(source, destination)


def build_release_archives() -> None:
    source = FINAL / "additional_files/Additional_file_3_ViralSafeTarget_0.10.0_Source.zip"
    source_copy = RELEASE / "ViralSafeTarget-0.10.0-source.zip"
    source_copy.parent.mkdir(parents=True, exist_ok=True)
    source_copy.write_bytes(source.read_bytes())
    (RELEASE / "ViralSafeTarget-0.10.0-source.sha256").write_text(
        f"{_sha256(source_copy)}  {source_copy.name}\n"
    )

    result_files = [
        ROOT / "reports/hsv2_genome_wide_exhaustive/gene_rankings.csv",
        ROOT / "reports/hsv2_genome_wide_exhaustive/top_candidates_global.csv",
        ROOT / "reports/hsv2_genome_wide_exhaustive/deep_screening_panel.csv",
        ROOT / "reports/hsv2_virtual_knockout_escape/gene_virtual_knockout.csv",
        ROOT / "reports/hsv2_virtual_knockout_escape/strategy_comparison.csv",
        ROOT / "reports/hsv2_tool_benchmark/rank_agreement.csv",
        ROOT / "reports/hsv2_tool_benchmark/top_k_overlap.csv",
        ROOT / "reports/hsv2_tool_benchmark/ablation_summary.csv",
        FINAL / "verified_statistics.json",
        FINAL / "build_metadata.json",
        ROOT / ".zenodo.json",
    ]
    data_zip = RELEASE / "ViralSafeTarget-0.10.0-compact-data.zip"
    _write_zip(
        data_zip,
        [(path, str(path.relative_to(ROOT))) for path in result_files],
    )
    (RELEASE / "ViralSafeTarget-0.10.0-compact-data.sha256").write_text(
        f"{_sha256(data_zip)}  {data_zip.name}\n"
    )


def build_submission_package() -> None:
    included: list[tuple[Path, str]] = []
    for path in sorted(FINAL.rglob("*")):
        if path.is_file():
            included.append(
                (
                    path,
                    str(Path("ViralSafeTarget_BMC_Submission_WORKING") / path.relative_to(FINAL)),
                )
            )
    included.extend(
        [
            (
                ROOT / "paper/bmc_bioinformatics/README.md",
                "ViralSafeTarget_BMC_Submission_WORKING/README.md",
            ),
            (
                ROOT / "paper/bmc_bioinformatics/MANUSCRIPT_SOURCE.md",
                "ViralSafeTarget_BMC_Submission_WORKING/MANUSCRIPT_SOURCE.md",
            ),
            (
                ROOT / ".zenodo.json",
                "ViralSafeTarget_BMC_Submission_WORKING/.zenodo.json",
            ),
        ]
    )
    _write_zip(PACKAGE, included)
    (PACKAGE.with_suffix(".sha256")).write_text(f"{_sha256(PACKAGE)}  {PACKAGE.name}\n")


def main() -> None:
    build_release_archives()
    build_submission_package()
    print(f"Built {PACKAGE}")


if __name__ == "__main__":
    main()
