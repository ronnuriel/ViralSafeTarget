"""Researcher-defined scorer protocol, example rules, and optional PyTorch support."""

from __future__ import annotations

import json
from collections.abc import Callable
from pathlib import Path
from typing import Any, Protocol, runtime_checkable

import numpy as np
import pandas as pd

SCORER_COLUMNS = ["candidate_id", "scorer_name", "raw_score", "confidence", "explanation"]


@runtime_checkable
class CandidateScorer(Protocol):
    name: str
    version: str

    def score(self, candidates: pd.DataFrame) -> pd.DataFrame: ...


class ExampleRuleScorer:
    """Transparent example combining conservation and sequence complexity."""

    name = "example_rule_scorer"
    version = "1.0"

    def score(self, candidates: pd.DataFrame) -> pd.DataFrame:
        conservation = pd.to_numeric(
            candidates.get("conservation_score", pd.Series(0.0, index=candidates.index)),
            errors="coerce",
        )
        complexity = pd.to_numeric(
            candidates.get("sequence_complexity_score", pd.Series(0.0, index=candidates.index)),
            errors="coerce",
        )
        score = 0.7 * conservation + 0.3 * complexity
        return pd.DataFrame(
            {
                "candidate_id": candidates["candidate_id"],
                "scorer_name": self.name,
                "raw_score": score,
                "confidence": pd.NA,
                "explanation": "0.7 × conservation_score + 0.3 × sequence_complexity_score",
            },
            columns=SCORER_COLUMNS,
        )


class TorchScorer:
    """Optional, explicit-model PyTorch inference without automatic pickle loading."""

    def __init__(
        self,
        *,
        model: Any,
        feature_extractor: Callable[[pd.DataFrame], Any],
        name: str,
        version: str,
        batch_size: int = 64,
        deterministic: bool = True,
        feature_names: list[str] | None = None,
    ) -> None:
        self.model = model
        self.feature_extractor = feature_extractor
        self.name = name
        self.version = version
        self.batch_size = batch_size
        self.deterministic = deterministic
        self.feature_names = feature_names or []

    @classmethod
    def from_state_dict(
        cls,
        *,
        model_factory: Callable[[], Any],
        state_dict_path: str | Path,
        feature_extractor: Callable[[pd.DataFrame], Any],
        name: str,
        version: str,
        **kwargs: Any,
    ) -> TorchScorer:
        try:
            import torch
        except ImportError as error:
            raise ImportError('TorchScorer requires pip install -e ".[torch]"') from error
        model = model_factory()
        try:
            state = torch.load(state_dict_path, map_location="cpu", weights_only=True)
        except TypeError:
            state = torch.load(state_dict_path, map_location="cpu")
        model.load_state_dict(state)
        return cls(
            model=model,
            feature_extractor=feature_extractor,
            name=name,
            version=version,
            **kwargs,
        )

    def score(self, candidates: pd.DataFrame) -> pd.DataFrame:
        try:
            import torch
        except ImportError as error:
            raise ImportError('TorchScorer requires pip install -e ".[torch]"') from error
        if self.deterministic:
            torch.manual_seed(0)
            torch.use_deterministic_algorithms(True, warn_only=True)
        features = self.feature_extractor(candidates.copy())
        tensor = torch.as_tensor(np.array(features, copy=True), dtype=torch.float32, device="cpu")
        self.model.to("cpu")
        self.model.eval()
        predictions = []
        with torch.no_grad():
            for start in range(0, len(tensor), self.batch_size):
                batch = self.model(tensor[start : start + self.batch_size])
                predictions.extend(batch.detach().cpu().reshape(-1).tolist())
        manifest = json.dumps(
            {
                "model_name": self.name,
                "model_version": self.version,
                "features": self.feature_names,
                "device": "cpu",
                "deterministic": self.deterministic,
            },
            sort_keys=True,
        )
        return pd.DataFrame(
            {
                "candidate_id": candidates["candidate_id"],
                "scorer_name": self.name,
                "raw_score": predictions,
                "confidence": pd.NA,
                "explanation": manifest,
            },
            columns=SCORER_COLUMNS,
        )
