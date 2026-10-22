from __future__ import annotations

from pathlib import Path

from viral_safe_target.cli import build_parser
from viral_safe_target.ui_app import _format_seconds, _resolve_project, _status_badge


def test_packaged_ui_helpers_and_cli_entrypoint(tmp_path: Path) -> None:
    project = tmp_path / "project.yaml"
    project.write_text("schema_version: '1.0'\n", encoding="utf-8")
    assert _resolve_project(str(project)) == project.resolve()
    assert _resolve_project(str(tmp_path / "missing.yaml")) is None
    assert _format_seconds(None) == "Unavailable"
    assert _format_seconds(90) == "1.5 min"
    assert "vst-completed" in _status_badge("completed")
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
    ):
        assert phrase in source
