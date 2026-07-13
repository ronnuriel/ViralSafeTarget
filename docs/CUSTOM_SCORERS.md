# Custom candidate scorers

A scorer supplies `name`, `version`, and `score(candidates)`. It returns `candidate_id`,
`scorer_name`, `raw_score`, optional `confidence`, and optional `explanation`.

```python
class MyScorer:
    name = "lab_model"
    version = "1.0"

    def score(self, candidates):
        return candidates[["candidate_id"]].assign(
            scorer_name=self.name,
            raw_score=candidates["conservation_score"],
            confidence=None,
            explanation="Lab-defined prioritization feature",
        )
```

`ExampleRuleScorer` is a tested API example, not a validated biological model. Convert a custom
score into a documented within-model rank before adding it to consensus.
