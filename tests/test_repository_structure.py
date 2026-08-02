from __future__ import annotations

import re
from pathlib import Path
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
MARKDOWN_LINK = re.compile(r"\[[^\]]*\]\(([^)]+)\)")


def _markdown_files() -> list[Path]:
    roots = [
        ROOT / "README.md",
        ROOT / "README_HE.md",
        ROOT / "CONTRIBUTING.md",
        ROOT / "DISCLAIMER.md",
        ROOT / "ROADMAP.md",
        ROOT / "SECURITY.md",
        ROOT / "REFERENCES.md",
    ]
    for directory in ("docs", "notebooks", "scripts", "configs", "schemas"):
        roots.extend(sorted((ROOT / directory).rglob("*.md")))
    roots.append(ROOT / "reports" / "README.md")
    return [path for path in roots if path.is_file()]


def _local_target(source: Path, raw_target: str) -> Path | None:
    target = raw_target.strip().split(maxsplit=1)[0].strip("<>")
    if not target or target.startswith(("#", "mailto:")) or "://" in target:
        return None
    without_anchor = unquote(target.split("#", maxsplit=1)[0])
    if not without_anchor:
        return None
    return (source.parent / without_anchor).resolve()


def test_repository_markdown_links_resolve() -> None:
    broken: list[str] = []
    for source in _markdown_files():
        text = source.read_text(encoding="utf-8")
        for raw_target in MARKDOWN_LINK.findall(text):
            target = _local_target(source, raw_target)
            if target is not None and not target.exists():
                broken.append(f"{source.relative_to(ROOT)} -> {raw_target}")
    assert not broken, "Broken repository links:\n" + "\n".join(broken)


def test_documentation_has_a_single_indexed_root() -> None:
    root_documents = sorted(path.name for path in (ROOT / "docs").glob("*.md"))
    assert root_documents == ["README.md"]
    assert {path.name for path in (ROOT / "docs").iterdir() if path.is_dir()} == {
        "getting-started",
        "maintenance",
        "reference",
        "research",
        "workflows",
    }


def test_every_script_and_notebook_is_listed() -> None:
    script_index = (ROOT / "scripts" / "README.md").read_text(encoding="utf-8")
    scripts = {
        path.name
        for path in (ROOT / "scripts").iterdir()
        if path.is_file() and path.name != "README.md"
    }
    assert not sorted(name for name in scripts if f"`{name}`" not in script_index)

    notebook_index = (ROOT / "notebooks" / "README.md").read_text(encoding="utf-8")
    notebooks = {path.name for path in (ROOT / "notebooks").glob("*.ipynb")}
    assert not sorted(name for name in notebooks if f"]({name})" not in notebook_index)
    assert not list((ROOT / "notebooks").glob("*_HE.ipynb"))


def test_public_reports_are_not_deleted_by_default_clean_target() -> None:
    makefile = (ROOT / "Makefile").read_text(encoding="utf-8")
    assert "rm -rf reports/*" not in makefile
    for directory in (
        "reports/hsv2_showcase",
        "reports/hsv2_genome_wide_exhaustive",
        "reports/hsv2_evidence_agent",
        "reports/hsv2_virtual_knockout_escape",
        "reports/hsv2_tool_benchmark",
    ):
        assert directory not in makefile.split("clean-generated:", maxsplit=1)[-1]
