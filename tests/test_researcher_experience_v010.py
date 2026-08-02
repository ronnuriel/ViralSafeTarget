from __future__ import annotations

import json
import zipfile
from pathlib import Path

import viral_safe_target.researcher as researcher
from viral_safe_target.project_workflow import run_project
from viral_safe_target.researcher import (
    create_demo_project,
    create_project,
    doctor_report,
    export_project,
    open_results,
    plan_project,
    tool_setup_report,
)


def test_demo_plan_run_bundle_resume_and_export(tmp_path: Path) -> None:
    project = create_demo_project(tmp_path / "demo")
    plan = plan_project(project)
    assert plan["candidate_count_estimate"] > 0
    assert plan["warning"].startswith("Runtime estimates")

    first = run_project(project)
    results = Path(first["output_root"])
    required = {
        "START_HERE.html",
        "SUMMARY.md",
        "summary.json",
        "top_guides.csv",
        "top_genes.csv",
        "research_shortlist.csv",
        "multiplex_panels.csv",
        "evidence_review_queue.tsv",
        "stage_timings.json",
        "run_manifest.json",
        "export.zip",
    }
    assert required <= {path.name for path in results.iterdir()}
    assert open_results(results, no_browser=True) == results / "START_HERE.html"
    start_here = (results / "START_HERE.html").read_text(encoding="utf-8")
    assert "What to do next" in start_here
    assert "guide_explanations/" in start_here
    assert "vst tools setup --tool cas-offinder" in start_here

    second = run_project(project)
    assert all(row["status"] in {"completed", "external_required"} for row in second["stages"])
    timings = json.loads((results / "stage_timings.json").read_text(encoding="utf-8"))
    assert any(row["cache_reused"] for row in timings["stages"])

    archive = export_project(project, output=tmp_path / "portable.zip")
    with zipfile.ZipFile(archive) as handle:
        names = set(handle.namelist())
    assert "results/START_HERE.html" in names
    assert not any("reference.fasta" in name for name in names)


def test_doctor_preserves_external_tool_missingness() -> None:
    report = doctor_report()
    assert report["viral_safe_target_version"] == "0.10.0"
    assert set(report["tools"]) >= {"mafft", "cas_offinder", "docker", "podman"}
    for tool in report["tools"].values():
        assert isinstance(tool["available"], bool)


def test_doctor_uses_supported_crispritz_version_probe(monkeypatch) -> None:
    calls: dict[str, list[str]] = {}

    def fake_program(command: str, version_args: list[str]) -> dict[str, object]:
        calls[command] = version_args
        return {"available": True, "path": f"/tools/{command}", "version": "test"}

    monkeypatch.setattr(researcher, "_program", fake_program)
    doctor_report()

    assert calls["crispritz.py"] == []


def test_tool_setup_guidance_is_scoped_actionable_and_non_mutating() -> None:
    report = tool_setup_report("cas-offinder")
    assert report["safety"].endswith("never executed by vst tools setup.")
    assert report["next_check"] == "vst tools status"
    assert [tool["id"] for tool in report["tools"]] == ["cas-offinder"]
    tool = report["tools"][0]
    assert tool["official_url"].startswith("https://github.com/")
    assert "source_build" in tool["commands"]


def test_sequence_only_mode_preserves_unavailable_annotation_stages(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    project = create_project(
        tmp_path / "sequence-only",
        project_name="sequence-only",
        virus_name="Sequence-only test",
        reference_fasta=root / "src/viral_safe_target/resources/demo/reference.fasta",
        sequence_only=True,
    )
    result = run_project(project)
    statuses = {row["stage"]: row["status"] for row in result["stages"]}
    assert statuses["discover"] == "completed"
    assert statuses["pairs"] == "unavailable"
    assert statuses["virtual_knockout"] == "disabled"
    assert statuses["host_screen"] == "external_required"


def test_canonical_notebook_is_english_and_output_free() -> None:
    root = Path(__file__).resolve().parents[1]
    path = root / "notebooks" / "00_VIRALSAFETARGET_END_TO_END_EN.ipynb"
    text = path.read_text(encoding="utf-8")
    assert not any("\u0590" <= character <= "\u05ff" for character in text)
    notebook = json.loads(text)
    assert all(
        cell.get("execution_count") is None and cell.get("outputs") == []
        for cell in notebook["cells"]
        if cell["cell_type"] == "code"
    )
