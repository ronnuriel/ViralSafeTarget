"""MAFFT command adapter; alignment remains delegated to MAFFT."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pandas as pd

from .base import AdapterError, ToolAvailability, ToolExecution, detect_executable


class MafftAdapter:
    name = "mafft"

    def detect(self) -> ToolAvailability:
        return detect_executable(self.name, "mafft", ("--version",))

    def build_input(self, candidates: pd.DataFrame, config: dict, output_dir: str | Path) -> Path:
        del candidates
        source = Path(config.get("input_fasta", ""))
        if not source.is_file():
            raise AdapterError("MAFFT input_fasta must point to an existing FASTA file.")
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        manifest = output / "mafft_input_manifest.json"
        manifest.write_text(
            json.dumps({"input_fasta": str(source.resolve())}, indent=2) + "\n",
            encoding="utf-8",
        )
        return manifest

    def run(
        self, input_path: str | Path, output_dir: str | Path, *, dry_run: bool = False
    ) -> ToolExecution:
        availability = self.detect()
        manifest = json.loads(Path(input_path).read_text(encoding="utf-8"))
        output = Path(output_dir) / "alignment.fasta"
        command = (availability.executable or "mafft", "--auto", manifest["input_fasta"])
        if dry_run:
            return ToolExecution(command, 0, "", "", output)
        if not availability.available:
            raise AdapterError(availability.message)
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        if completed.returncode:
            raise AdapterError(f"MAFFT failed ({completed.returncode}): {completed.stderr.strip()}")
        output.write_text(completed.stdout, encoding="utf-8")
        return ToolExecution(command, completed.returncode, "", completed.stderr, output)

    def parse(self, output_path: str | Path, manifest: str | Path | None = None) -> pd.DataFrame:
        del manifest
        path = Path(output_path)
        if not path.is_file():
            raise AdapterError(f"MAFFT output is missing: {path}")
        return pd.DataFrame([{"alignment_file": str(path.resolve()), "status": "completed"}])

    def normalize(self, parsed_results: pd.DataFrame) -> pd.DataFrame:
        return parsed_results.copy()
