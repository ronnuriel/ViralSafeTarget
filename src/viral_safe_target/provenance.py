"""Research-run provenance and checksums."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import shlex
import subprocess
import sys
from collections.abc import Iterable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _command_output(command: list[str], root: Path) -> str | None:
    try:
        return subprocess.check_output(
            command, cwd=root, text=True, stderr=subprocess.STDOUT
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _version(command: list[str], root: Path) -> str | None:
    try:
        completed = subprocess.run(
            command,
            cwd=root,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    lines = completed.stdout.strip().splitlines()
    return lines[0] if lines else None


def _package_versions() -> dict[str, str]:
    packages = ["viral-safe-target", "biopython", "pandas", "numpy", "PyYAML"]
    versions: dict[str, str] = {}
    for package in packages:
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not-installed"
    return versions


def write_run_manifest(
    output_path: str | Path,
    input_paths: Iterable[str | Path],
    parameters: dict[str, object],
    project_root: str | Path | None = None,
    *,
    config_path: str | Path | None = None,
    editor_profile: dict[str, Any] | None = None,
    accepted_accessions: Iterable[str] | None = None,
    rejected_accessions: Iterable[dict[str, str]] | None = None,
    human_assembly_identifier: str | None = None,
    command_line: Iterable[str] | None = None,
    random_seed: int | None = None,
    output_paths: Iterable[str | Path] | None = None,
) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    root = Path(project_root or Path.cwd()).resolve()

    def file_record(path_value: str | Path) -> dict[str, object]:
        path = Path(path_value)
        return {
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }

    config_record = file_record(config_path) if config_path else None
    outputs = [
        file_record(path)
        for path in (output_paths or [])
        if Path(path).is_file() and Path(path).resolve() != output.resolve()
    ]
    git_status = _command_output(["git", "status", "--porcelain"], root)
    manifest = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": _command_output(["git", "rev-parse", "HEAD"], root),
        "dirty_worktree": bool(git_status),
        "configuration": config_record,
        "editor_profile": editor_profile,
        "inputs": [file_record(path) for path in input_paths],
        "accepted_accessions": sorted(accepted_accessions or []),
        "rejected_accessions": list(rejected_accessions or []),
        "human_assembly_identifier": human_assembly_identifier,
        "command_line": shlex.join(command_line or sys.argv),
        "parameters": parameters,
        "package_versions": _package_versions(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "mafft_version": _version(["mafft", "--version"], root),
        "cas_offinder_version": _version(["cas-offinder", "--help"], root),
        "random_seed": random_seed,
        "outputs": outputs,
    }
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return output
