"""Cas-OFFinder adapter preserving native parser semantics and provenance."""

from __future__ import annotations

import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from ..offtarget import (
    build_cas_offinder_input,
    read_cas_offinder_output,
    summarize_cas_offinder_hits,
)
from ..provenance import sha256_file
from ..tables import TOOL_RESULT_COLUMNS
from .base import AdapterError, ToolAvailability, ToolExecution, detect_executable


class CasOffinderAdapter:
    name = "cas-offinder"

    def detect(self) -> ToolAvailability:
        return detect_executable(self.name, "cas-offinder", ("--help",))

    def build_input(self, candidates: pd.DataFrame, config: dict, output_dir: str | Path) -> Path:
        human_directory = config.get("human_fasta_directory")
        if not human_directory:
            raise AdapterError("Cas-OFFinder requires config['human_fasta_directory'].")
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        input_path = output / "cas_offinder_input.txt"
        manifest_path = output / "cas_offinder_candidates.csv"
        build_cas_offinder_input(
            candidates,
            human_directory,
            input_path,
            manifest_path,
            maximum_candidates=int(config.get("maximum_candidates", len(candidates))),
            stratify_by_gene=bool(config.get("stratify_by_gene", False)),
            config=config.get("vst_config"),
        )
        return input_path

    def run(
        self, input_path: str | Path, output_dir: str | Path, *, dry_run: bool = False
    ) -> ToolExecution:
        availability = self.detect()
        output = Path(output_dir) / "cas_offinder_output.tsv"
        command = (availability.executable or "cas-offinder", str(input_path), "C", str(output))
        if dry_run:
            return ToolExecution(command, 0, "", "", output)
        if not availability.available:
            raise AdapterError(availability.message)
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        if completed.returncode:
            raise AdapterError(
                f"Cas-OFFinder failed ({completed.returncode}): "
                f"{(completed.stderr or completed.stdout).strip()}"
            )
        return ToolExecution(
            command, completed.returncode, completed.stdout, completed.stderr, output
        )

    def parse(self, output_path: str | Path, manifest: str | Path | None = None) -> pd.DataFrame:
        del manifest
        path = Path(output_path)
        if not path.is_file():
            raise AdapterError(
                f"Cas-OFFinder output is missing: {path}. Missing output is not a zero-hit result."
            )
        return read_cas_offinder_output(path)

    def normalize(
        self,
        parsed_results: pd.DataFrame,
        *,
        candidates: pd.DataFrame | None = None,
        manifest: str | Path | None = None,
        source_file: str | Path | None = None,
        version: str | None = None,
        assembly: str = "",
        editor: str = "",
        command: str = "imported native output",
    ) -> pd.DataFrame:
        if candidates is None:
            raise AdapterError("Cas-OFFinder normalization requires the candidate table.")
        summarized = summarize_cas_offinder_hits(
            candidates, parsed_results, selected_manifest=manifest
        )
        source = Path(source_file) if source_file else None
        now = datetime.now(timezone.utc).isoformat()
        rows = []
        for _, row in summarized.iterrows():
            count = int(row.get("human_total_predicted_hits", 0))
            rows.append(
                {
                    "candidate_id": row["candidate_id"],
                    "guide_sequence": row.get("guide_sequence", ""),
                    "gene_name": row.get("gene_name", ""),
                    "tool_name": self.name,
                    "tool_version": version or self.detect().version or "unknown",
                    "tool_mode": "reference-genome search",
                    "genome_or_assembly": assembly,
                    "editor": editor,
                    "metric_name": "predicted_offtarget_burden",
                    "raw_value": count,
                    "normalized_value": pd.NA,
                    "rank": pd.NA,
                    "percentile_rank": pd.NA,
                    "decision": (
                        "no predicted hit within configured model and threshold"
                        if count == 0
                        else "requires expert review"
                    ),
                    "explanation": json.dumps(
                        {
                            "exact": int(row.get("human_exact_hit_count", 0)),
                            "one_mismatch": int(row.get("human_one_mismatch_hit_count", 0)),
                            "two_mismatch": int(row.get("human_two_mismatch_hit_count", 0)),
                            "three_mismatch": int(row.get("human_three_mismatch_hit_count", 0)),
                        },
                        sort_keys=True,
                    ),
                    "source_file": str(source.resolve()) if source else "",
                    "source_file_sha256": sha256_file(source)
                    if source and source.is_file()
                    else "",
                    "command_or_import_method": command,
                    "timestamp": now,
                    "status": "completed",
                    "error_message": "",
                }
            )
        output = pd.DataFrame(rows, columns=TOOL_RESULT_COLUMNS)
        if output.empty:
            return output
        numeric = pd.to_numeric(output["raw_value"], errors="coerce")
        output["rank"] = numeric.rank(method="min", ascending=True)
        n = len(output)
        output["percentile_rank"] = 1.0 if n == 1 else 1.0 - (output["rank"] - 1.0) / (n - 1.0)
        output["normalized_value"] = output["percentile_rank"]
        return output
