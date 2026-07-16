import importlib.util
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
WORK = ROOT / "paper/bmc_bioinformatics"
FINAL = WORK / "final"


def _load_figure_builder():
    script = ROOT / "scripts/build_bmc_figures.py"
    spec = importlib.util.spec_from_file_location("build_bmc_figures", script)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_bmc_source_statistics_are_frozen() -> None:
    observed = _load_figure_builder().validate_sources()
    assert observed["initial_candidates"] == 28_578
    assert observed["eligible_candidates"] == 23_108
    assert observed["unique_guides"] == 21_654
    assert observed["human_matches"] == 440_341
    assert observed["zero_hit_rows"] == 2_668


def test_bmc_manuscript_has_required_software_article_sections() -> None:
    source = (WORK / "MANUSCRIPT_SOURCE.md").read_text()
    required = [
        "## Abstract",
        "## Background",
        "## Implementation",
        "## Results",
        "## Discussion",
        "## Conclusions",
        "## Availability and requirements",
        "## List of abbreviations",
        "## Declarations",
        "## References",
        "## Figure legends",
        "## Additional files",
    ]
    for heading in required:
        assert heading in source
    assert "DO NOT SUBMIT" in source
    assert "*" in source


def test_human_evidence_review_is_unresolved() -> None:
    review = pd.read_csv(FINAL / "HUMAN_REVIEW_REQUIRED.csv")
    assert not review.empty
    assert set(review["status"]) == {"pending"}
    assert review["reviewer_name"].isna().all()
    assert review["decision"].isna().all()
    critical = review[review.source_identifier == "PMID:24794394"]
    assert len(critical) >= 1
    assert critical.review_issue.str.contains("CRITICAL").all()


def test_bmc_figures_are_separate_and_below_upload_limit() -> None:
    for number in range(1, 7):
        for suffix in ("pdf", "png"):
            figure = FINAL / "figures" / f"Figure_{number}.{suffix}"
            assert figure.exists()
            assert figure.stat().st_size < 10 * 1024 * 1024


def test_bmc_working_files_remain_blocked() -> None:
    required = [
        "Main_Manuscript_BMC_Bioinformatics.docx",
        "Main_Manuscript_BMC_Bioinformatics.pdf",
        "Cover_Letter_BMC_Bioinformatics.docx",
        "BMC_PreSubmission_Checklist.md",
        "BMC_Submission_Metadata.json",
    ]
    for name in required:
        assert (FINAL / name).exists()
    assert "DO NOT SUBMIT" in (FINAL / "BMC_PreSubmission_Checklist.md").read_text()
