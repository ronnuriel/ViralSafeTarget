"""Typed DataFrame wrappers and the normalized external-tool schema."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from typing import Any

import pandas as pd

TOOL_RESULT_COLUMNS = [
    "candidate_id",
    "guide_sequence",
    "gene_name",
    "tool_name",
    "tool_version",
    "tool_mode",
    "genome_or_assembly",
    "editor",
    "metric_name",
    "raw_value",
    "normalized_value",
    "rank",
    "percentile_rank",
    "decision",
    "explanation",
    "source_file",
    "source_file_sha256",
    "command_or_import_method",
    "timestamp",
    "status",
    "error_message",
]


@dataclass(frozen=True)
class CandidateTable:
    """A validated candidate table that still behaves like a pandas container."""

    dataframe: pd.DataFrame

    def __post_init__(self) -> None:
        required = {"candidate_id", "guide_sequence"}
        missing = sorted(required - set(self.dataframe.columns))
        if missing:
            raise ValueError(f"Candidate table is missing required columns: {missing}")
        if self.dataframe["candidate_id"].duplicated().any():
            raise ValueError("Candidate IDs must be unique")

    def __len__(self) -> int:
        return len(self.dataframe)

    def __getitem__(self, key: Any) -> Any:
        return self.dataframe.__getitem__(key)

    def __iter__(self) -> Iterator[str]:
        return iter(self.dataframe)

    def __getattr__(self, name: str) -> Any:
        return getattr(self.dataframe, name)

    def copy(self) -> pd.DataFrame:
        return self.dataframe.copy()


@dataclass(frozen=True)
class ToolResultTable:
    """Long-form normalized tool results with strict missing-value semantics."""

    dataframe: pd.DataFrame

    def __post_init__(self) -> None:
        missing = sorted(set(TOOL_RESULT_COLUMNS) - set(self.dataframe.columns))
        if missing:
            raise ValueError(f"Tool result table is missing normalized columns: {missing}")

    @classmethod
    def empty(cls) -> ToolResultTable:
        return cls(pd.DataFrame(columns=TOOL_RESULT_COLUMNS))

    @classmethod
    def from_frame(cls, frame: pd.DataFrame) -> ToolResultTable:
        output = frame.copy()
        for column in TOOL_RESULT_COLUMNS:
            if column not in output:
                output[column] = pd.NA
        return cls(output[TOOL_RESULT_COLUMNS])

    def __len__(self) -> int:
        return len(self.dataframe)

    def __getitem__(self, key: Any) -> Any:
        return self.dataframe.__getitem__(key)

    def __getattr__(self, name: str) -> Any:
        return getattr(self.dataframe, name)

    def copy(self) -> pd.DataFrame:
        return self.dataframe.copy()


def as_dataframe(value: pd.DataFrame | CandidateTable | ToolResultTable) -> pd.DataFrame:
    """Return an isolated DataFrame from a supported public table type."""
    if isinstance(value, (CandidateTable, ToolResultTable)):
        return value.dataframe.copy()
    if isinstance(value, pd.DataFrame):
        return value.copy()
    raise TypeError(f"Expected a pandas DataFrame or ViralSafeTarget table, got {type(value)!r}")
