from __future__ import annotations

from collections import Counter
from collections.abc import Mapping

import pandas as pd

from .io_utils import require_aligned

VALID_BASES = set("ACGT")


def conservation_profile(records: Mapping[str, str]) -> pd.DataFrame:
    """Calculate per-column consensus, identity and gap/ambiguous fractions."""
    alignment_length = require_aligned(records)
    n = len(records)
    rows: list[dict] = []

    sequences = list(records.values())
    for idx in range(alignment_length):
        column = [seq[idx] for seq in sequences]
        base_counts = Counter(base for base in column if base in VALID_BASES)
        informative = sum(base_counts.values())
        consensus, top_count = (base_counts.most_common(1)[0] if base_counts else ("N", 0))
        rows.append({
            "alignment_position_1based": idx + 1,
            "consensus": consensus,
            "identity_all_sequences": top_count / n,
            "identity_informative": top_count / informative if informative else 0.0,
            "informative_count": informative,
            "gap_or_ambiguous_fraction": 1.0 - informative / n,
        })
    return pd.DataFrame(rows)


def find_conserved_runs(
    profile: pd.DataFrame,
    identity_threshold: float = 0.95,
    max_gap_fraction: float = 0.0,
    min_length: int = 20,
) -> pd.DataFrame:
    """Return contiguous alignment runs that satisfy conservation thresholds."""
    mask = (
        (profile["identity_all_sequences"] >= identity_threshold)
        & (profile["gap_or_ambiguous_fraction"] <= max_gap_fraction)
    ).tolist()

    runs: list[dict] = []
    start = None
    for i, good in enumerate(mask):
        if good and start is None:
            start = i
        if start is not None and (not good or i == len(mask) - 1):
            end = i if good and i == len(mask) - 1 else i - 1
            length = end - start + 1
            if length >= min_length:
                subset = profile.iloc[start : end + 1]
                runs.append({
                    "alignment_start_1based": start + 1,
                    "alignment_end_1based": end + 1,
                    "length": length,
                    "mean_identity": subset["identity_all_sequences"].mean(),
                })
            start = None
    return pd.DataFrame(runs)
