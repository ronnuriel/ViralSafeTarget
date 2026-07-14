# Optional PyTorch scorer

Install only when needed:

```bash
python -m pip install -e ".[torch]"
```

`TorchScorer` requires an explicit model object or a model factory plus a `state_dict`. It performs
CPU batch inference, can request deterministic algorithms, and records model/version/features in
the explanation field. It does not automatically load arbitrary pickled model objects.

```python
from viral_safe_target.scorers import TorchScorer

scorer = TorchScorer(
    model=model,
    feature_extractor=my_feature_extractor,
    name="lab_efficiency_model",
    version="1.0",
    feature_names=["feature_a", "feature_b"],
)
scored = scorer.score(run.candidates)
```

Core installation and CI do not require PyTorch.
