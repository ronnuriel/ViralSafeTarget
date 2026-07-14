"""Documented-export import adapters for CRISPOR, CHOPCHOP, and GuideScan2."""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from ..provenance import sha256_file
from ..tables import TOOL_RESULT_COLUMNS, ToolResultTable
from .base import AdapterError, ToolAvailability, ToolExecution

SUPPORTED_IMPORT_TOOLS = {"crispor", "chopchop", "guidescan2"}


@dataclass(frozen=True)
class ExternalImportResult:
    results: ToolResultTable
    unmatched_rows: pd.DataFrame
    ambiguous_rows: pd.DataFrame
    raw_rows: pd.DataFrame


def _read_export(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload, dict):
            payload = payload.get("results", payload.get("data", [payload]))
        return pd.json_normalize(payload)
    return pd.read_csv(path, sep=None, engine="python")


def _normalize_sequence(value: object) -> str:
    return str(value or "").upper().replace(" ", "").replace("-", "")


def _match_candidates(
    row: pd.Series, candidates: pd.DataFrame, columns: dict[str, str]
) -> pd.DataFrame:
    matches = candidates
    comparisons = {
        "guide_sequence": lambda series, value: series.map(_normalize_sequence).eq(
            _normalize_sequence(value)
        ),
        "pam": lambda series, value: series.astype(str).str.upper().eq(str(value).upper()),
        "strand": lambda series, value: series.astype(str).eq(str(value)),
        "coordinate": lambda series, value: pd.to_numeric(series, errors="coerce").eq(
            pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
        ),
    }
    candidate_columns = {
        "guide_sequence": "guide_sequence",
        "pam": "pam",
        "strand": "strand",
        "coordinate": "reference_start_1based",
    }
    for identity, imported_column in columns.items():
        candidate_column = candidate_columns.get(identity)
        if (
            identity not in comparisons
            or imported_column not in row.index
            or candidate_column not in matches
            or pd.isna(row[imported_column])
            or str(row[imported_column]).strip() == ""
        ):
            continue
        matches = matches[comparisons[identity](matches[candidate_column], row[imported_column])]
    return matches


def _metric_configuration(mapping: dict[str, Any], raw_columns: list[str]) -> list[dict[str, str]]:
    configured = []
    for source, value in mapping.get("metrics", {}).items():
        if isinstance(value, str):
            configured.append({"source": source, "name": value, "direction": "higher"})
        else:
            configured.append(
                {
                    "source": source,
                    "name": value.get("name", source),
                    "direction": value.get("direction", "higher"),
                }
            )
    identity_sources = set(mapping.get("columns", {}).values())
    known_sources = {item["source"] for item in configured}
    for source in raw_columns:
        if source not in identity_sources and source not in known_sources:
            configured.append({"source": source, "name": f"raw.{source}", "direction": "none"})
    return configured


def load_external_results(
    tool: str,
    input_path: str | Path,
    mapping_path: str | Path,
    candidates: pd.DataFrame,
) -> ExternalImportResult:
    """Map a researcher-supplied export into the normalized long schema."""
    tool = tool.lower()
    if tool not in SUPPORTED_IMPORT_TOOLS:
        raise ValueError(
            f"Unsupported import tool {tool!r}; choose {sorted(SUPPORTED_IMPORT_TOOLS)}"
        )
    source = Path(input_path)
    mapping_file = Path(mapping_path)
    if not source.is_file():
        raise FileNotFoundError(f"External result export is missing: {source}")
    mapping = yaml.safe_load(mapping_file.read_text(encoding="utf-8")) or {}
    raw = _read_export(source)
    raw = raw.copy()
    raw["_source_row"] = range(1, len(raw) + 1)
    columns = mapping.get("columns", {})
    metrics = _metric_configuration(
        mapping, [str(column) for column in raw.columns if column != "_source_row"]
    )
    now = datetime.now(timezone.utc).isoformat()
    normalized_rows: list[dict[str, object]] = []
    unmatched: list[dict[str, object]] = []
    ambiguous: list[dict[str, object]] = []
    for _, row in raw.iterrows():
        matches = _match_candidates(row, candidates, columns)
        if len(matches) != 1:
            record = row.to_dict()
            record["tool_name"] = tool
            record["mapping_status"] = "unmatched" if matches.empty else "ambiguous"
            record["candidate_matches"] = ";".join(
                matches.get("candidate_id", pd.Series(dtype=str))
            )
            (unmatched if matches.empty else ambiguous).append(record)
            continue
        candidate = matches.iloc[0]
        for metric in metrics:
            source_column = metric["source"]
            if source_column not in row:
                continue
            raw_value = row[source_column]
            numeric = pd.to_numeric(pd.Series([raw_value]), errors="coerce").iloc[0]
            normalized_rows.append(
                {
                    "candidate_id": candidate["candidate_id"],
                    "guide_sequence": candidate.get("guide_sequence", ""),
                    "gene_name": candidate.get("gene_name", ""),
                    "tool_name": tool,
                    "tool_version": str(mapping.get("tool_version", "researcher supplied")),
                    "tool_mode": "documented export import",
                    "genome_or_assembly": str(mapping.get("genome_or_assembly", "")),
                    "editor": str(mapping.get("editor", "")),
                    "metric_name": metric["name"],
                    "raw_value": raw_value,
                    "normalized_value": pd.NA,
                    "rank": pd.NA,
                    "percentile_rank": pd.NA,
                    "decision": "imported; interpretation depends on documented metric",
                    "explanation": f"Imported {source_column}; direction={metric['direction']}",
                    "source_file": str(source.resolve()),
                    "source_file_sha256": sha256_file(source),
                    "command_or_import_method": f"mapping={mapping_file.resolve()}",
                    "timestamp": now,
                    "status": "completed",
                    "error_message": "" if pd.notna(raw_value) else "metric missing in source row",
                    "_direction": metric["direction"],
                    "_numeric": numeric,
                }
            )
    results = pd.DataFrame(normalized_rows)
    if results.empty:
        normalized = pd.DataFrame(columns=TOOL_RESULT_COLUMNS)
    else:
        for (_, _metric), indexes in results.groupby(["tool_name", "metric_name"]).groups.items():
            direction = str(results.loc[indexes[0], "_direction"])
            numeric = pd.to_numeric(results.loc[indexes, "_numeric"], errors="coerce")
            if direction not in {"higher", "lower"} or numeric.notna().sum() == 0:
                continue
            ranks = numeric.rank(method="min", ascending=direction == "lower", na_option="keep")
            count = numeric.notna().sum()
            percentile = 1.0 if count == 1 else 1.0 - (ranks - 1.0) / (count - 1.0)
            results.loc[indexes, "rank"] = ranks
            results.loc[indexes, "percentile_rank"] = percentile
            results.loc[indexes, "normalized_value"] = percentile
        normalized = results[TOOL_RESULT_COLUMNS].copy()
    return ExternalImportResult(
        results=ToolResultTable.from_frame(normalized),
        unmatched_rows=pd.DataFrame(unmatched),
        ambiguous_rows=pd.DataFrame(ambiguous),
        raw_rows=raw,
    )


class ExternalImportAdapter:
    """Import-only adapter; public web interfaces are intentionally not scraped."""

    def __init__(self, name: str) -> None:
        if name.lower() not in SUPPORTED_IMPORT_TOOLS:
            raise ValueError(f"Unsupported external import adapter: {name}")
        self.name = name.lower()

    def detect(self) -> ToolAvailability:
        return ToolAvailability(
            self.name,
            True,
            "import adapter 0.4",
            execution_mode="import-only",
            message="Researcher-supplied CSV/TSV/JSON exports are supported.",
        )

    def build_input(self, candidates: pd.DataFrame, config: dict, output_dir: str | Path) -> Path:
        del config
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        path = output / f"{self.name}_candidate_reference.csv"
        candidates.to_csv(path, index=False)
        return path

    def run(
        self, input_path: str | Path, output_dir: str | Path, *, dry_run: bool = False
    ) -> ToolExecution:
        del output_dir, dry_run
        raise AdapterError(
            f"{self.name} is import-only. Export results through the tool's documented interface "
            f"and import them; candidate reference: {input_path}"
        )

    def parse(self, output_path: str | Path, manifest: str | Path | None = None) -> pd.DataFrame:
        del manifest
        return _read_export(Path(output_path))

    def normalize(self, parsed_results: pd.DataFrame) -> pd.DataFrame:
        raise AdapterError("Use load_external_results with a mapping YAML for normalization.")
