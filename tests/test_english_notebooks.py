from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = [
    ROOT / "notebooks/advanced/benchmarking/07_HSV2_MULTITOOL_COMPARISON_EN.ipynb",
    ROOT / "notebooks/advanced/genome-wide/08_RUN_FULL_PIPELINE_EN.ipynb",
    ROOT / "notebooks/advanced/genome-wide/09_HSV2_GENOME_WIDE_DISCOVERY_EN.ipynb",
    ROOT / "notebooks/advanced/disruption-escape/10_HSV2_GENE_FUNCTION_AND_DISRUPTION_EN.ipynb",
    ROOT / "notebooks/advanced/genome-wide/11_HSV2_RESEARCH_SHOWCASE_EN.ipynb",
    ROOT / "notebooks/advanced/genome-wide/12_HSV2_HELDOUT_POPULATION_VALIDATION_EN.ipynb",
    ROOT / "notebooks/advanced/evidence/13_EVIDENCE_AGENT_HUMAN_REVIEW_EN.ipynb",
    ROOT / "notebooks/advanced/disruption-escape/14_VIRTUAL_KNOCKOUT_AND_ESCAPE_EN.ipynb",
    ROOT / "notebooks/advanced/benchmarking/15_SYSTEMATIC_MULTITOOL_BENCHMARK_EN.ipynb",
]


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_v04_notebooks_are_english_and_store_no_generated_outputs():
    assert not list((ROOT / "notebooks").glob("*_HE.ipynb"))
    for path in NOTEBOOKS:
        notebook = _load(path)
        assert not re.search(r"[\u0590-\u05ff]", path.read_text(encoding="utf-8"))
        for cell in notebook["cells"]:
            if cell["cell_type"] == "code":
                assert cell["execution_count"] is None
                assert cell["outputs"] == []


def test_end_to_end_notebook_executes_all_cells_in_synthetic_mode(monkeypatch):
    matplotlib.use("Agg")
    monkeypatch.chdir(ROOT)
    notebook = _load(NOTEBOOKS[1])
    namespace: dict[str, object] = {}
    executed = 0
    for number, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] != "code":
            continue
        source = "".join(cell["source"])
        exec(compile(source, f"{NOTEBOOKS[1]}:cell-{number}", "exec"), namespace)
        executed += 1
    assert executed == 12
    assert namespace["SYNTHETIC_MODE"] is True
    assert len(namespace["post_human_candidates"]) > 0
    assert namespace["tool_coverage"]["status"].eq("pending").any()


def test_genome_wide_notebook_executes_in_english_synthetic_mode(monkeypatch):
    matplotlib.use("Agg")
    # Jupyter starts kernels in the notebook directory, not necessarily repository root.
    monkeypatch.chdir(ROOT / "notebooks")
    notebook = _load(NOTEBOOKS[2])
    namespace: dict[str, object] = {}
    for number, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] == "code":
            source = "".join(cell["source"])
            exec(compile(source, f"{NOTEBOOKS[2]}:cell-{number}", "exec"), namespace)
    assert namespace["SYNTHETIC_MODE"] is True
    assert len(namespace["screening_panel"]) > 0
    assert namespace["summary"]["external_tool_status"] == "pending Cas-OFFinder completion"


def test_evidence_agent_notebook_executes_without_network_or_automatic_approval(monkeypatch):
    monkeypatch.chdir(ROOT / "notebooks")
    path = ROOT / "notebooks/advanced/evidence/13_EVIDENCE_AGENT_HUMAN_REVIEW_EN.ipynb"
    notebook = _load(path)
    namespace: dict[str, object] = {}
    for number, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] == "code":
            source = "".join(cell["source"])
            exec(compile(source, f"{path}:cell-{number}", "exec"), namespace)
    assert namespace["RUN_NETWORK"] is False
    assert len(namespace["gene_catalog"]) > 0
    assert namespace["review_checks"]["automatic_approval"] is False


def test_virtual_knockout_escape_notebook_executes_from_clean_synthetic_inputs(monkeypatch):
    monkeypatch.chdir(ROOT / "notebooks")
    path = ROOT / "notebooks/advanced/disruption-escape/14_VIRTUAL_KNOCKOUT_AND_ESCAPE_EN.ipynb"
    notebook = _load(path)
    namespace: dict[str, object] = {}
    for number, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] == "code":
            source = "".join(cell["source"])
            exec(compile(source, f"{path}:cell-{number}", "exec"), namespace)
    assert namespace["MODE"] == "synthetic"
    assert namespace["summary"]["guide_count"] > 0
    assert namespace["strategy_comparison"]["combined_therapeutic_score"].isna().all()


def test_systematic_benchmark_notebook_executes_in_synthetic_mode(monkeypatch):
    monkeypatch.chdir(ROOT / "notebooks")
    path = ROOT / "notebooks/advanced/benchmarking/15_SYSTEMATIC_MULTITOOL_BENCHMARK_EN.ipynb"
    notebook = _load(path)
    namespace: dict[str, object] = {}
    for number, cell in enumerate(notebook["cells"]):
        if cell["cell_type"] == "code":
            source = "".join(cell["source"])
            exec(compile(source, f"{path}:cell-{number}", "exec"), namespace)
    assert namespace["MODE"] == "synthetic"
    assert set(namespace["ablation"]["variant"]) == {
        "all_components",
        "without_conservation",
        "without_uniqueness",
    }
