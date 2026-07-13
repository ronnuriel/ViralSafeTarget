from __future__ import annotations

import pandas as pd

from viral_safe_target.experimental import import_crispresso2_results


def test_crispresso2_synthetic_import_keeps_measurements_separate(tmp_path):
    result_directory = tmp_path / "CRISPResso_on_test"
    result_directory.mkdir()
    pd.DataFrame(
        [
            {
                "Amplicon": "amp1",
                "Reads_aligned": 1000,
                "Modified%": 12.5,
                "Insertions%": 2.0,
                "Deletions%": 9.0,
                "Substitutions%": 1.5,
                "Frameshift%": 7.0,
                "In-frame%": 5.5,
                "Quantification_window": "-3:3",
            }
        ]
    ).to_csv(
        result_directory / "CRISPResso_quantification_of_editing_frequency.txt",
        sep="\t",
        index=False,
    )
    candidate_map = tmp_path / "map.csv"
    pd.DataFrame(
        [{"candidate_id": "c1", "guide_sequence": "A" * 20, "amplicon_name": "amp1"}]
    ).to_csv(candidate_map, index=False)
    imported = import_crispresso2_results(result_directory, candidate_map)
    measurement = imported.measurements.iloc[0]
    assert measurement["modified_read_percentage"] == 12.5
    assert measurement["measurement_type"] == "measured experimental metric"
    assert "prediction" not in measurement["measurement_type"]
    assert imported.metadata["measurement_count"] == 1
