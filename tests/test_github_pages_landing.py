import re
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
INDEX = DOCS / "index.html"


class _ReferenceParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.references: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        values = dict(attrs)
        for attribute in ("href", "src"):
            value = values.get(attribute)
            if value:
                self.references.append(value)


def _normalized_text() -> str:
    return " ".join(INDEX.read_text(encoding="utf-8").split())


def test_landing_page_assets_and_local_links_exist() -> None:
    assert INDEX.is_file()
    assert (DOCS / ".nojekyll").is_file()
    assert (DOCS / "assets/landing/styles.css").is_file()
    assert (DOCS / "assets/landing/script.js").is_file()
    assert (DOCS / "assets/landing/favicon.svg").is_file()

    parser = _ReferenceParser()
    parser.feed(INDEX.read_text(encoding="utf-8"))
    for reference in parser.references:
        if reference.startswith(("#", "http://", "https://", "mailto:", "data:")):
            continue
        clean_reference = reference.split("?", 1)[0].split("#", 1)[0]
        assert (DOCS / clean_reference).resolve().is_file(), reference


def test_landing_page_preserves_frozen_hsv2_claims_and_boundaries() -> None:
    text = _normalized_text()
    for expected in (
        "28,578",
        "23,108",
        "21,654",
        "109 / 109",
        "2,668",
        "0.880",
        "49 / 50",
        "UL3",
        "UL30",
        "UL36",
    ):
        assert expected in text

    assert "No predicted hit under a configured search model is not proof of safety" in text
    assert "No claim of editing, safety, viral inhibition, delivery, treatment" in text
    assert "they are not therapeutic scores" in text


def test_landing_page_is_english_only_and_does_not_claim_a_pypi_release() -> None:
    text = INDEX.read_text(encoding="utf-8")
    assert re.search(r"[\u0590-\u05FF]", text) is None
    assert 'pip install "viral-safe-target' not in text
    assert 'pip install ".[ui,notebooks]"' in text
    assert "public package release remains a separate approval step" in text
