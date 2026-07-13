from __future__ import annotations

import pandas as pd

from viral_safe_target.integrations import load_external_results


def test_generic_import_mapping_unmatched_ambiguous_and_unknown_metrics(tmp_path):
    candidates = pd.DataFrame(
        [
            {
                "candidate_id": "c1",
                "guide_sequence": "A" * 20,
                "pam": "AGG",
                "strand": "+",
                "reference_start_1based": 10,
                "gene_name": "UL19",
            },
            {
                "candidate_id": "c2",
                "guide_sequence": "C" * 20,
                "pam": "TGG",
                "strand": "+",
                "reference_start_1based": 20,
                "gene_name": "UL30",
            },
            {
                "candidate_id": "c3",
                "guide_sequence": "C" * 20,
                "pam": "TGG",
                "strand": "-",
                "reference_start_1based": 30,
                "gene_name": "UL30",
            },
        ]
    )
    export = tmp_path / "export.tsv"
    pd.DataFrame(
        [
            {"guide": "A" * 20, "pam": "AGG", "score": 80, "note": "known"},
            {"guide": "C" * 20, "pam": "TGG", "score": 50, "note": "ambiguous"},
            {"guide": "G" * 20, "pam": "CGG", "score": 20, "note": "unmatched"},
        ]
    ).to_csv(export, sep="\t", index=False)
    mapping = tmp_path / "mapping.yaml"
    mapping.write_text(
        "columns:\n"
        "  guide_sequence: guide\n"
        "  pam: pam\n"
        "metrics:\n"
        "  score:\n"
        "    name: predicted_efficiency\n"
        "    direction: higher\n",
        encoding="utf-8",
    )
    imported = load_external_results("crispor", export, mapping, candidates)
    assert set(imported.results.dataframe["metric_name"]) == {
        "predicted_efficiency",
        "raw.note",
    }
    assert len(imported.ambiguous_rows) == 1
    assert len(imported.unmatched_rows) == 1
    known = imported.results.dataframe.query("metric_name == 'predicted_efficiency'").iloc[0]
    assert known["normalized_value"] == 1.0
    unknown = imported.results.dataframe.query("metric_name == 'raw.note'").iloc[0]
    assert pd.isna(unknown["normalized_value"])
