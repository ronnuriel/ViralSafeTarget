from __future__ import annotations

import os
import stat
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from viral_safe_target.notebook_helpers import (
    clear_cache_stamps,
    detect_cas_offinder,
    find_project_root,
    load_notebook_run,
    result_funnel,
    run_streaming,
    safe_read_csv,
    valid_cas_offinder_output,
)

ROOT = Path(__file__).resolve().parents[1]


def test_project_root_detection_from_nested_directory():
    assert find_project_root(ROOT / "notebooks") == ROOT


def test_streaming_runner_uses_argument_list_and_surfaces_failure(tmp_path, capsys):
    result = run_streaming([sys.executable, "-c", "print('streamed')"], cwd=tmp_path)
    assert result.returncode == 0
    assert "streamed" in capsys.readouterr().out
    with pytest.raises(subprocess.CalledProcessError):
        run_streaming([sys.executable, "-c", "raise SystemExit(3)"], cwd=tmp_path)


def test_cas_offinder_resolution_order_and_output_validation(tmp_path, monkeypatch):
    configured = tmp_path / "configured-cas"
    configured.write_text("#!/bin/sh\n", encoding="utf-8")
    configured.chmod(configured.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", os.environ.get("PATH", ""))
    assert detect_cas_offinder(tmp_path, configured) == configured.resolve()
    output = tmp_path / "cas.tsv"
    assert not valid_cas_offinder_output(output)
    output.write_text("query\tchr1\t1\tobserved\t+\t2\n", encoding="utf-8")
    assert valid_cas_offinder_output(output)


def test_safe_csv_requires_expected_schema(tmp_path):
    path = tmp_path / "table.csv"
    pd.DataFrame([{"candidate_id": "c1"}]).to_csv(path, index=False)
    assert len(safe_read_csv(path, ["candidate_id"])) == 1
    with pytest.raises(ValueError, match="guide_sequence"):
        safe_read_csv(path, ["candidate_id", "guide_sequence"])


def test_force_rerun_helper_only_removes_repository_cache_stamps(tmp_path):
    cache = tmp_path / "reports/run/.cache"
    cache.mkdir(parents=True)
    stamp = cache / "stage.json"
    stamp.write_text("{}", encoding="utf-8")
    retained = cache / "notes.txt"
    retained.write_text("keep", encoding="utf-8")
    assert clear_cache_stamps(tmp_path, ["reports/run/.cache"]) == [stamp]
    assert retained.is_file()
    with pytest.raises(ValueError, match="Refusing"):
        clear_cache_stamps(tmp_path, ["reports/run"])


def test_synthetic_notebook_data_and_funnel_use_bundled_outputs():
    data = load_notebook_run(ROOT, synthetic=True)
    assert not data["candidates"].empty
    assert {"candidate_id", "guide_sequence", "gene_name"} <= set(data["candidates"].columns)
    funnel = result_funnel(data)
    assert funnel["stage"].str.contains("Selected pilot").any()
    assert funnel["count"].ge(0).all()
