"""Common executable/import adapter contracts and provenance helpers."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Protocol, runtime_checkable

import pandas as pd


class AdapterError(RuntimeError):
    """An actionable external-tool integration failure."""


@dataclass(frozen=True)
class ToolAvailability:
    name: str
    available: bool
    version: str | None = None
    executable: str | None = None
    execution_mode: str | None = None
    message: str = ""


@dataclass(frozen=True)
class ToolExecution:
    command: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str
    output_path: Path | None = None


def detect_executable(
    name: str,
    executable: str,
    version_arguments: tuple[str, ...],
) -> ToolAvailability:
    path = shutil.which(executable)
    if not path:
        return ToolAvailability(
            name=name,
            available=False,
            message=(
                f"{name} was not found on PATH; use build-input/import mode or "
                "install it separately."
            ),
        )
    try:
        completed = subprocess.run(
            [path, *version_arguments],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
        lines = (completed.stdout or completed.stderr).strip().splitlines()
        version = lines[0] if lines else "available (version not reported)"
    except (OSError, subprocess.TimeoutExpired) as error:
        return ToolAvailability(
            name=name,
            available=False,
            executable=path,
            message=f"{name} version detection failed: {error}",
        )
    return ToolAvailability(name, True, version, path, "native")


@runtime_checkable
class ToolAdapter(Protocol):
    name: str

    def detect(self) -> ToolAvailability: ...

    def build_input(
        self, candidates: pd.DataFrame, config: dict, output_dir: str | Path
    ) -> Path: ...

    def run(
        self, input_path: str | Path, output_dir: str | Path, *, dry_run: bool = False
    ) -> ToolExecution: ...

    def parse(
        self, output_path: str | Path, manifest: str | Path | None = None
    ) -> pd.DataFrame: ...

    def normalize(self, parsed_results: pd.DataFrame) -> pd.DataFrame: ...
