from __future__ import annotations

from pathlib import Path

from viral_safe_target.cli import build_parser
from viral_safe_target.researcher import resource_path
from viral_safe_target.ui_app import (
    _format_seconds,
    _remaining_estimate,
    _resolve_project,
    _status_badge,
    _workflow_cards,
)


def test_packaged_ui_helpers_and_cli_entrypoint(tmp_path: Path) -> None:
    project = tmp_path / "project.yaml"
    project.write_text("schema_version: '1.0'\n", encoding="utf-8")
    assert _resolve_project(str(project)) == project.resolve()
    assert _resolve_project(str(tmp_path / "missing.yaml")) is None
    assert _format_seconds(None) == "Unavailable"
    assert _format_seconds(90) == "1.5 min"
    assert "vst-completed" in _status_badge("completed")
    assert "vst-running" in _status_badge("running")
    cards = _workflow_cards(
        {"validate": "completed", "discover": "running"},
        estimates={"discover": 90.0},
    )
    assert "Guide discovery &amp; conservation" in cards
    assert "Estimated: 1.5 min" in cards
    statuses = {stage: "completed" for stage in (
        "validate",
        "discover",
        "host_screen",
        "pairs",
        "virtual_knockout",
        "escape",
        "multiplex",
        "report",
        "bundle",
    )}
    statuses["discover"] = "pending"
    assert _remaining_estimate(statuses, {"discover": 4.0}) == 4.0
    assert resource_path("import_maps", "crispor.yaml").is_file()
    args = build_parser().parse_args(["ui", "--port", "8765", "--headless"])
    assert args.port == 8765
    assert args.headless is True


def test_readme_explains_pipeline_and_visual_workflow() -> None:
    root = Path(__file__).resolve().parents[1]
    readme = (root / "README.md").read_text(encoding="utf-8")
    assert "## How the research pipeline works" in readme
    assert "```mermaid" in readme
    assert "external_required" in readme
    assert "not collapsed into a therapeutic score" in readme
    assert "vst ui" in readme


def test_public_ui_is_english_and_preserves_scientific_boundaries() -> None:
    root = Path(__file__).resolve().parents[1]
    source = (root / "src" / "viral_safe_target" / "ui_app.py").read_text(encoding="utf-8")
    assert not any("\u0590" <= character <= "\u05ff" for character in source)
    for phrase in (
        "Research boundary",
        "external_required",
        "Human review queue",
        "not biological evidence",
        "do not establish editing",
        "Tool & evidence hub",
        "Current stage",
    ):
        assert phrase in source


def test_interoperability_document_uses_explicit_status_and_official_sources() -> None:
    root = Path(__file__).resolve().parents[1]
    text = (
        root / "docs/reference/TOOL_ECOSYSTEM_AND_INTEROPERABILITY.md"
    ).read_text(encoding="utf-8")
    for phrase in (
        "Versioned import map",
        "named human reviewer",
        "not connected",
        "CRISPResso2",
        "NCBI E-utilities",
        "Official sources",
    ):
        assert phrase in text
