"""Known-target recovery benchmark support."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

BENCHMARK_COLUMNS = [
    "virus",
    "virus_type",
    "reference_accession",
    "editor",
    "target_sequence",
    "PAM",
    "gene",
    "evidence_context",
    "experimental_system",
    "source_identifier",
    "evidence_level",
    "expected_status",
    "notes",
]


def run_benchmark(
    candidates: pd.DataFrame, known_targets: pd.DataFrame
) -> tuple[pd.DataFrame, dict[str, object]]:
    missing = [column for column in BENCHMARK_COLUMNS if column not in known_targets]
    if missing:
        raise ValueError(f"Benchmark table is missing columns: {', '.join(missing)}")
    rows: list[dict[str, object]] = []
    pre_rank = {candidate: rank + 1 for rank, candidate in enumerate(candidates["candidate_id"])}
    if "post_human_score" in candidates:
        post_order = candidates.sort_values(
            ["post_human_score", "candidate_id"], ascending=[False, True], kind="mergesort"
        )
        post_rank = {
            candidate: rank + 1 for rank, candidate in enumerate(post_order["candidate_id"])
        }
    else:
        post_rank = {}
    for _, target in known_targets.iterrows():
        target_sequence = str(target["target_sequence"]).upper()
        pam = str(target["PAM"]).upper()
        matches = candidates[
            (candidates["guide_sequence"].astype(str).str.upper() == target_sequence)
            & (candidates["pam"].astype(str).str.upper() == pam)
        ]
        if matches.empty:
            candidate_id = ""
            exclusion = "not_regenerated_or_filtered_before_candidate_output"
            regenerated = False
            first_rank = pd.NA
            final_rank = pd.NA
            percentile = pd.NA
        else:
            match = matches.sort_values(
                ["pre_human_score", "candidate_id"], ascending=[False, True]
            ).iloc[0]
            candidate_id = str(match["candidate_id"])
            regenerated = True
            first_rank = pre_rank[candidate_id]
            final_rank = post_rank.get(candidate_id, pd.NA)
            percentile = 1.0 - ((first_rank - 1) / max(1, len(candidates)))
            exclusion = str(match.get("rejection_reasons", ""))
        rows.append(
            {
                **target.to_dict(),
                "candidate_id": candidate_id,
                "regenerated": regenerated,
                "pre_human_rank": first_rank,
                "post_human_rank": final_rank,
                "pre_human_percentile": percentile,
                "exclusion_reason": exclusion,
            }
        )
    detail = pd.DataFrame(rows)
    expected = detail[
        detail["expected_status"].astype(str).str.lower().isin({"positive", "recover"})
    ]
    denominator = len(expected)
    summary: dict[str, object] = {
        "known_target_count": len(detail),
        "expected_positive_count": denominator,
        "recovered_positive_count": int(expected["regenerated"].sum()) if denominator else 0,
        "recovery_rate": float(expected["regenerated"].mean()) if denominator else None,
        "top_k_recall": {
            str(k): (
                float((pd.to_numeric(expected["pre_human_rank"], errors="coerce") <= k).mean())
                if denominator
                else None
            )
            for k in (10, 50, 100)
        },
        "gene_level_recovery": (
            expected.groupby("gene")["regenerated"].mean().dropna().to_dict() if denominator else {}
        ),
    }
    return detail, summary


def write_benchmark_outputs(
    detail: pd.DataFrame, summary: dict[str, object], output_directory: str | Path
) -> None:
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    detail.to_csv(output / "benchmark_results.csv", index=False)
    (output / "benchmark_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
