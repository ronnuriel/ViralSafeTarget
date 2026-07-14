from __future__ import annotations

import pandas as pd
import pytest

torch = pytest.importorskip("torch")

from viral_safe_target.scorers import TorchScorer  # noqa: E402


def test_optional_torch_scorer_cpu_deterministic():
    model = torch.nn.Linear(2, 1, bias=False)
    with torch.no_grad():
        model.weight[:] = torch.tensor([[1.0, 2.0]])
    candidates = pd.DataFrame(
        [{"candidate_id": "c1", "a": 1.0, "b": 2.0}, {"candidate_id": "c2", "a": 2.0, "b": 3.0}]
    )
    scorer = TorchScorer(
        model=model,
        feature_extractor=lambda frame: frame[["a", "b"]].to_numpy(),
        name="synthetic",
        version="1",
        feature_names=["a", "b"],
    )
    first = scorer.score(candidates)
    second = scorer.score(candidates)
    pd.testing.assert_frame_equal(first, second)
    assert first["raw_score"].tolist() == [5.0, 8.0]
