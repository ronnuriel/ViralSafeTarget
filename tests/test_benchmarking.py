from __future__ import annotations

import pandas as pd

from viral_safe_target.benchmarking import BENCHMARK_COLUMNS, run_benchmark


def test_benchmark_reports_recovery_rank_and_missing_reason():
    candidates = pd.DataFrame(
        [
            {
                "candidate_id": "c1",
                "guide_sequence": "A" * 20,
                "pam": "AGG",
                "pre_human_score": 0.9,
                "rejection_reasons": "",
            }
        ]
    )
    known = pd.DataFrame(
        [
            dict.fromkeys(BENCHMARK_COLUMNS, "")
            | {
                "target_sequence": "A" * 20,
                "PAM": "AGG",
                "gene": "G1",
                "expected_status": "positive",
            },
            dict.fromkeys(BENCHMARK_COLUMNS, "")
            | {
                "target_sequence": "C" * 20,
                "PAM": "CGG",
                "gene": "G2",
                "expected_status": "positive",
            },
        ]
    )
    detail, summary = run_benchmark(candidates, known)
    assert detail["regenerated"].tolist() == [True, False]
    assert detail.iloc[0]["pre_human_rank"] == 1
    assert detail.iloc[1]["exclusion_reason"]
    assert summary["recovery_rate"] == 0.5
