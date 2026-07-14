from __future__ import annotations

import json
import shutil
from pathlib import Path

from viral_safe_target.project_workflow import (
    initialize_project,
    project_status,
    run_project,
    validate_project,
)
from viral_safe_target.reproduction import hsv2_reproduction_plan, reproduce_hsv2

ROOT = Path(__file__).resolve().parents[1]


def _project(tmp_path: Path, reference: str = "HSV2_demo_ref") -> Path:
    project = initialize_project(
        tmp_path / "virus-project",
        project_id="demo-virus",
        display_name="Demo virus",
        reference_accession=reference,
    )
    root = project.parent
    shutil.copyfile(ROOT / "data/demo/virus_aligned.fasta", root / "data/reference.fasta")
    shutil.copyfile(ROOT / "data/demo/virus_aligned.fasta", root / "data/strains.aligned.fasta")
    shutil.copyfile(ROOT / "data/demo/reference.gff3", root / "data/reference.gff3")
    shutil.copyfile(ROOT / "data/demo/human_mini.fasta", root / "external/host/host.fasta")
    return project


def test_new_virus_project_runs_and_resumes_without_inventing_host_results(tmp_path: Path) -> None:
    project = _project(tmp_path)
    checks = validate_project(project)
    assert not checks["status"].eq("fail").any()

    first = run_project(project)
    stages = {row["stage"]: row for row in first["stages"]}
    assert stages["discover"]["status"] == "completed"
    assert stages["host_screen"]["status"] == "external_required"
    assert stages["pairs"]["status"] == "completed"
    assert stages["report"]["status"] == "completed"

    output = project.parent / "results"
    discovery = output / "discovery/discovery_panel.csv"
    report = output / "report/report.html"
    manifest = output / "run_manifest.json"
    assert discovery.is_file() and report.is_file() and manifest.is_file()
    assert not (output / "host_screen/candidates_ranked_post_host.csv").exists()
    first_mtime = discovery.stat().st_mtime_ns

    second = run_project(project)
    assert discovery.stat().st_mtime_ns == first_mtime
    assert second == project_status(project)
    state = json.loads((output / "workflow_state.json").read_text(encoding="utf-8"))
    assert state["stages"]["host_screen"]["status"] == "external_required"


def test_project_validation_rejects_reference_identity_mismatch(tmp_path: Path) -> None:
    project = _project(tmp_path, reference="NOT_IN_ALIGNMENT")
    checks = validate_project(project)
    failed = checks[checks["status"].eq("fail")]
    assert "aligned viral panel" in set(failed["component"])
    assert "annotation/reference identity" in set(failed["component"])


def test_hsv2_reproduction_defaults_to_a_non_mutating_plan() -> None:
    plan = hsv2_reproduction_plan(ROOT)
    stages = [row["stage"] for row in plan]
    assert stages == [
        "environment",
        "public_data_and_candidates",
        "ortholog_reference",
        "genome_wide_host_screen",
        "heldout_population",
        "gene_function_and_showcase",
    ]
    assert "hsv2_discovery_accessions.txt" in plan[1]["command"]
    result = reproduce_hsv2(ROOT)
    assert result["mode"] == "plan"
    assert len(result["stages"]) == 6
