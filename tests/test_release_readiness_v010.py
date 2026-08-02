from __future__ import annotations

from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def test_release_workflow_requires_manual_production_publication() -> None:
    workflow = (ROOT / ".github/workflows/release.yml").read_text(encoding="utf-8")
    assert "publish_target:" in workflow
    assert "inputs.publish_target == 'testpypi'" in workflow
    assert "inputs.publish_target == 'pypi'" in workflow
    assert "startsWith(github.ref, 'refs/tags/v')" in workflow
    assert "environment: testpypi" in workflow
    assert "environment: pypi" in workflow
    assert "canonical-notebook:" in workflow
    assert "VST_NOTEBOOK_MODE: demo" in workflow
    assert "VST_NOTEBOOK_MODE: hsv2_snapshot" in workflow


def test_public_author_metadata_is_complete() -> None:
    citation = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    authors = citation["authors"]
    assert [(row["given-names"], row["family-names"]) for row in authors] == [
        ("Ron", "Nuriel"),
        ("Sarel", "Cohen"),
    ]
    assert authors[0]["orcid"].endswith("0009-0008-3970-2591")
    assert authors[1]["orcid"].endswith("0000-0003-4578-1245")


def test_installation_docs_distinguish_source_and_public_release() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    quickstart = (
        ROOT / "docs/getting-started/FIVE_MINUTE_QUICKSTART.md"
    ).read_text(encoding="utf-8")
    for text in (readme, quickstart):
        normalized = " ".join(text.lower().split())
        assert "before the first pypi release" in normalized
        assert "git+https://github.com/ronnuriel/ViralSafeTarget.git@main" in text
        assert 'pip install "viral-safe-target[notebooks]"' in text
