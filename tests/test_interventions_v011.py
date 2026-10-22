from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import yaml

from viral_safe_target.interventions import (
    INTERVENTION_EVIDENCE_COLUMNS,
    build_intervention_landscape,
    validate_intervention_evidence,
    write_intervention_outputs,
)
from viral_safe_target.project_workflow import initialize_project, run_project
from viral_safe_target.researcher import create_demo_project


def _record(**updates: str) -> dict[str, str]:
    record = {column: "" for column in INTERVENTION_EVIDENCE_COLUMNS}
    record.update(
        {
            "record_id": "IE-1",
            "gene_name": "UL30",
            "target_virus": "HSV-2",
            "evidence_virus": "HSV-2",
            "intervention_name": "example inhibitor",
            "intervention_type": "small_molecule",
            "target_role": "direct_target",
            "development_status": "clinical_investigation",
            "evidence_category": "clinical_trial",
            "evidence_scope": "direct_target_virus",
            "source_database": "PubMed",
            "source_identifier": "PMID:1",
            "source_url": "https://pubmed.ncbi.nlm.nih.gov/1/",
            "finding": "Example source-linked finding.",
            "review_status": "pending",
        }
    )
    record.update(updates)
    return record


def test_project_init_creates_intervention_evidence_template(tmp_path: Path) -> None:
    project = initialize_project(
        tmp_path / "project",
        project_id="virus-test",
        display_name="Virus test",
    )
    profile = yaml.safe_load((project.parent / "profiles" / "virus.yaml").read_text())
    assert profile["intervention_evidence_table"] == "evidence/intervention_evidence.tsv"
    table = pd.read_csv(project.parent / "evidence" / "intervention_evidence.tsv", sep="\t")
    assert list(table.columns) == INTERVENTION_EVIDENCE_COLUMNS


def test_approved_intervention_requires_named_dated_review() -> None:
    table = pd.DataFrame([_record(review_status="approved")])
    with pytest.raises(ValueError, match="require reviewer and review_date"):
        validate_intervention_evidence(table)


def test_pending_records_are_marked_and_missingness_is_explicit(tmp_path: Path) -> None:
    evidence = pd.DataFrame([_record()])
    landscape = build_intervention_landscape(evidence, genes=["UL3", "UL30"])
    ul30 = landscape.loc[landscape["gene_name"].eq("UL30")].iloc[0]
    ul3 = landscape.loc[landscape["gene_name"].eq("UL3")].iloc[0]
    assert ul30["review_marker"] == "*"
    assert ul30["interpretation"] == "pending human review"
    assert ul3["interventions"] == "not located in supplied evidence"

    outputs = write_intervention_outputs(
        evidence,
        tmp_path,
        genes=["UL3", "UL30"],
        virus_names=["HSV-2"],
    )
    report = outputs["report"].read_text(encoding="utf-8")
    assert "pending human review" in report
    assert "does not recommend a drug" in report
    assert outputs["review_queue"].is_file()


def test_result_bundle_contains_intervention_landscape(tmp_path: Path) -> None:
    project = create_demo_project(tmp_path / "demo")
    run_project(project)
    results = project.parent / "results"
    assert (results / "intervention_landscape.html").is_file()
    assert (results / "intervention_landscape.tsv").is_file()
    assert (results / "intervention_review_queue.tsv").is_file()
    assert "Target-to-intervention evidence" in (results / "START_HERE.html").read_text()
