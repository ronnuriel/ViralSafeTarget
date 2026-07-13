"""Stable researcher-facing run loading API."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from .tables import CandidateTable


def _read_optional(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.is_file() else pd.DataFrame()


@dataclass(frozen=True)
class ResearchRun:
    """A completed run loaded from its machine-readable artifacts."""

    path: Path
    candidates: pd.DataFrame
    human_hits: pd.DataFrame
    same_gene_pairs: pd.DataFrame
    multi_target_pairs: pd.DataFrame
    manifest: dict[str, Any]

    @property
    def candidate_table(self) -> CandidateTable:
        return CandidateTable(self.candidates.copy())


def load_run(run_directory: str | Path) -> ResearchRun:
    """Load a v0.3+ run without mutating or regenerating its artifacts."""
    directory = Path(run_directory)
    if not directory.is_dir():
        raise FileNotFoundError(f"Run directory does not exist: {directory}")
    candidate_paths = [
        directory / "candidates_ranked_post_human.csv",
        directory / "candidates_ranked_pre_human.csv",
        directory / "candidates.csv",
    ]
    candidate_path = next((path for path in candidate_paths if path.is_file()), None)
    if candidate_path is None:
        raise FileNotFoundError(
            f"No candidate table found in {directory}; expected one of "
            + ", ".join(path.name for path in candidate_paths)
        )
    candidates = pd.read_csv(candidate_path)
    CandidateTable(candidates)
    manifest_path = directory / "run_manifest.json"
    manifest = (
        json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.is_file() else {}
    )
    return ResearchRun(
        path=directory.resolve(),
        candidates=candidates,
        human_hits=_read_optional(directory / "predicted_human_hits.csv"),
        same_gene_pairs=_read_optional(directory / "pair_hypotheses_same_gene.csv"),
        multi_target_pairs=_read_optional(directory / "pair_hypotheses_multi_target.csv"),
        manifest=manifest,
    )
