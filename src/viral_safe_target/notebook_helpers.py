"""Small, testable helpers for researcher-facing orchestration notebooks."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

import pandas as pd

from . import __version__
from .sdk import load_run


def find_project_root(start: str | Path | None = None) -> Path:
    """Find the repository root from the root itself or a nested notebook directory."""
    current = Path(start or Path.cwd()).expanduser().resolve()
    if current.is_file():
        current = current.parent
    for candidate in [current, *current.parents]:
        if (candidate / "pyproject.toml").is_file() and (candidate / "scripts").is_dir():
            return candidate
    raise FileNotFoundError(
        f"Could not find the ViralSafeTarget root from {current}; "
        "open the notebook inside the repository."
    )


def run_streaming(
    command: Sequence[str | Path],
    *,
    cwd: str | Path,
    check: bool = True,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run an argument list without a shell and stream combined stdout/stderr."""
    arguments = [str(item) for item in command]
    print("$", " ".join(arguments), flush=True)
    process = subprocess.Popen(
        arguments,
        cwd=Path(cwd),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    captured: list[str] = []
    assert process.stdout is not None
    for line in process.stdout:
        print(line, end="", flush=True)
        captured.append(line)
    returncode = process.wait()
    result = subprocess.CompletedProcess(arguments, returncode, "".join(captured), "")
    if check and returncode:
        raise subprocess.CalledProcessError(returncode, arguments, output=result.stdout)
    return result


def _command_summary(executable: str, arguments: Sequence[str]) -> tuple[str, str]:
    path = shutil.which(executable)
    if not path:
        return "missing", f"{executable} was not found on PATH"
    try:
        completed = subprocess.run(
            [path, *arguments], capture_output=True, text=True, timeout=10, check=False
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        return "warning", str(error)
    lines = (completed.stdout or completed.stderr).strip().splitlines()
    return "ready", f"{path} — {lines[0] if lines else 'available'}"


def detect_cas_offinder(
    project_root: str | Path,
    configured_path: str | Path | None = None,
) -> Path | None:
    """Resolve Cas-OFFinder in configured, repository-local, then PATH order."""
    candidates: list[Path] = []
    if configured_path:
        candidates.append(Path(configured_path).expanduser())
    candidates.append(Path(project_root) / "tools/bin/cas-offinder")
    on_path = shutil.which("cas-offinder")
    if on_path:
        candidates.append(Path(on_path))
    for candidate in candidates:
        resolved = candidate.resolve()
        if resolved.is_file() and os.access(resolved, os.X_OK):
            return resolved
    return None


def valid_cas_offinder_output(path: str | Path) -> bool:
    """Check that an output contains at least one supported tab-delimited result row."""
    output = Path(path)
    if not output.is_file() or output.stat().st_size == 0:
        return False
    try:
        with output.open(encoding="utf-8") as handle:
            for line in handle:
                if line.strip() and not line.startswith("#"):
                    return len(line.rstrip("\n").split("\t")) >= 6
    except (OSError, UnicodeDecodeError):
        return False
    return False


def safe_read_csv(
    path: str | Path,
    required_columns: Iterable[str] = (),
    *,
    label: str = "CSV",
) -> pd.DataFrame:
    """Read a CSV only after existence and expected-schema checks."""
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(f"{label} is missing: {source}")
    frame = pd.read_csv(source)
    missing = sorted(set(required_columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")
    return frame


def clear_cache_stamps(
    project_root: str | Path, relative_cache_directories: Iterable[str | Path]
) -> list[Path]:
    """Remove only JSON cache stamps under explicitly named repository .cache directories."""
    root = Path(project_root).resolve()
    removed: list[Path] = []
    for relative in relative_cache_directories:
        directory = (root / relative).resolve()
        if root not in directory.parents or directory.name != ".cache":
            raise ValueError(f"Refusing to clear a non-repository .cache directory: {directory}")
        if not directory.is_dir():
            continue
        for stamp in directory.glob("*.json"):
            stamp.unlink()
            removed.append(stamp)
    return removed


def environment_status(project_root: str | Path) -> pd.DataFrame:
    """Return a concise doctor table with actionable missing-stage guidance."""
    root = Path(project_root)
    git_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=root, capture_output=True, text=True, check=False
    ).stdout.strip()
    git_dirty = bool(
        subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=root,
            capture_output=True,
            text=True,
            check=False,
        ).stdout.strip()
    )
    rows = [
        {
            "component": "Python",
            "status": "ready",
            "details": f"{sys.executable} — {sys.version.split()[0]}",
            "action": "",
        },
        {
            "component": "ViralSafeTarget",
            "status": "ready",
            "details": __version__,
            "action": "python -m pip install -e .",
        },
        {
            "component": "Git",
            "status": "warning" if git_dirty else "ready",
            "details": f"{git_commit[:12]} — {'dirty' if git_dirty else 'clean'}",
            "action": "Review git status before recording final provenance" if git_dirty else "",
        },
    ]
    checks = [
        ("MAFFT", "mafft", ("--version",), "Install MAFFT in the active environment"),
        ("NCBI Datasets", "datasets", ("version",), "Install ncbi-datasets-cli"),
        (
            "Cas-OFFinder",
            "cas-offinder",
            ("--help",),
            "Set CAS_OFFINDER_PATH or install Cas-OFFinder",
        ),
        ("OpenCL", "clinfo", ("--version",), "Install/configure an OpenCL runtime"),
        (
            "CRISPRitz",
            "crispritz.py",
            ("--help",),
            "Optional: install CRISPRitz or use its documented Docker/import workflow",
        ),
    ]
    for label, executable, arguments, action in checks:
        status, details = _command_summary(executable, arguments)
        if label == "CRISPRitz" and status == "missing" and shutil.which("docker"):
            status = "pending"
            details = "Native CLI missing; Docker is present, image availability not assumed"
        rows.append(
            {
                "component": label,
                "status": status,
                "details": details,
                "action": "" if status == "ready" else action,
            }
        )
    hsv_cached = (root / "data/processed/hsv2_aligned_25.fasta").is_file()
    human_cached = any((root / "data/raw/human_GRCh38").glob("**/*_genomic.fna"))
    free_gib = shutil.disk_usage(root).free / (1024**3)
    rows.extend(
        [
            {
                "component": "Cached HSV-2 data",
                "status": "ready" if hsv_cached else "pending",
                "details": "alignment present" if hsv_cached else "alignment absent",
                "action": "bash scripts/run_real_hsv2.sh --sample-size 25"
                if not hsv_cached
                else "",
            },
            {
                "component": "Cached GRCh38 data",
                "status": "ready" if human_cached else "pending",
                "details": "genomic FASTA present" if human_cached else "genomic FASTA absent",
                "action": "bash scripts/run_real_hsv2.sh --with-human --sample-size 25"
                if not human_cached
                else "",
            },
            {
                "component": "Free disk",
                "status": "ready" if free_gib >= 20 else "warning",
                "details": f"{free_gib:.2f} GiB",
                "action": "Free disk space before downloading genome data" if free_gib < 20 else "",
            },
        ]
    )
    return pd.DataFrame(rows, columns=["component", "status", "details", "action"])


def load_notebook_run(project_root: str | Path, *, synthetic: bool) -> dict[str, Any]:
    """Load consistent real or bundled synthetic tables for notebook analysis."""
    root = Path(project_root)
    run_directory = root / ("reports/demo" if synthetic else "reports/hsv2_pilot")
    if synthetic and not (run_directory / "candidates.csv").is_file():
        run_directory = root / "examples/demo_output"
    run = load_run(run_directory)
    candidates = run.candidates.copy()
    required = {"candidate_id", "guide_sequence", "gene_name", "decision"}
    missing = sorted(required - set(candidates.columns))
    if missing:
        raise ValueError(f"Candidate table is missing notebook fields: {missing}")
    pre_path = run_directory / "candidates_ranked_pre_human.csv"
    rejected_path = run_directory / "candidates_rejected_pre_human.csv"
    pre = pd.read_csv(pre_path) if pre_path.is_file() else candidates.copy()
    rejected = pd.read_csv(rejected_path) if rejected_path.is_file() else candidates.iloc[0:0]
    selected_path = run_directory / "offtarget_selected_candidates.csv"
    selected = pd.read_csv(selected_path) if selected_path.is_file() else candidates.copy()
    qc_path = root / "reports/real_hsv2/accession_qc.csv"
    if not synthetic and qc_path.is_file():
        qc = pd.read_csv(qc_path)
    else:
        qc = pd.DataFrame(
            {
                "accession": run.manifest.get("accepted_accessions", []),
                "decision": "accepted",
                "rejection_reason": "synthetic fixture",
            }
        )
    same_pairs = run.same_gene_pairs.copy()
    if synthetic and same_pairs.empty:
        pair_path = run_directory / "simulated_pairs.csv"
        same_pairs = pd.read_csv(pair_path) if pair_path.is_file() else pd.DataFrame()
    return {
        "run_directory": run_directory,
        "run": run,
        "candidates": candidates,
        "pre_human": pre,
        "rejected_pre_human": rejected,
        "selected": selected,
        "human_hits": run.human_hits.copy(),
        "same_gene_pairs": same_pairs,
        "multi_target_pairs": run.multi_target_pairs.copy(),
        "qc": qc,
        "manifest": run.manifest,
    }


def result_funnel(data: dict[str, Any]) -> pd.DataFrame:
    """Build the notebook's deterministic result-funnel table."""
    qc = data["qc"]
    candidates = data["candidates"]
    pre = data["pre_human"]
    rejected = data["rejected_pre_human"]
    selected = data["selected"]
    no_hit = (
        int(pd.to_numeric(candidates["human_total_predicted_hits"], errors="coerce").eq(0).sum())
        if "human_total_predicted_hits" in candidates
        else 0
    )
    expert = (
        int(candidates["decision"].astype(str).eq("expert_review_required").sum())
        if "decision" in candidates
        else 0
    )
    values = [
        ("Downloaded/QC genome records", len(qc)),
        ("Accepted genomes", int(qc.get("decision", pd.Series(dtype=str)).eq("accepted").sum())),
        ("Rejected genomes", int(qc.get("decision", pd.Series(dtype=str)).eq("rejected").sum())),
        ("Initial ranked candidates", len(pre) + len(rejected)),
        ("Selected pilot candidates", len(selected)),
        ("Expert-review candidates", expert),
        ("No predicted hit within configured threshold", no_hit),
    ]
    return pd.DataFrame(values, columns=["stage", "count"])
