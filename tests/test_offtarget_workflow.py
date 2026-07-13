from __future__ import annotations

from pathlib import Path

import pandas as pd

from viral_safe_target import (
    build_cas_offinder_input,
    rank_post_human_candidates,
    read_cas_offinder_output,
    summarize_cas_offinder_hits,
)

ROOT = Path(__file__).resolve().parents[1]


def _candidates() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "candidate_id": "c1",
                "guide_sequence": "TACGATGCTAACCGGTTAAC",
                "gene_name": "G1",
                "pre_human_score": 0.9,
                "rejection_reasons": "",
            },
            {
                "candidate_id": "c2",
                "guide_sequence": "TTGCAACGTTGCAACGTTGC",
                "gene_name": "G2",
                "pre_human_score": 0.8,
                "rejection_reasons": "",
            },
        ]
    )


def test_cas_offinder_input_and_manifest_are_valid_and_deterministic(tmp_path):
    output = tmp_path / "input.txt"
    manifest = tmp_path / "manifest.csv"
    selected = build_cas_offinder_input(
        _candidates(),
        tmp_path,
        output,
        manifest,
        maximum_candidates=2,
    )
    lines = output.read_text(encoding="utf-8").splitlines()
    assert lines[1] == "N" * 20 + "NGG"
    assert lines[2].endswith("NNN 3")
    assert "c1" not in lines[2]
    assert selected["candidate_id"].tolist() == pd.read_csv(manifest)["candidate_id"].tolist()


def test_output_parsing_mismatch_counts_and_post_human_ranking(tmp_path):
    manifest_path = tmp_path / "manifest.csv"
    candidates = _candidates()
    selected = candidates.copy()
    selected["cas_offinder_query"] = selected["guide_sequence"] + "NNN"
    selected.to_csv(manifest_path, index=False)
    hits = read_cas_offinder_output(ROOT / "tests/fixtures/cas_offinder_output.tsv")
    summarized = summarize_cas_offinder_hits(
        candidates,
        hits,
        selected_manifest=manifest_path,
    )
    c1 = summarized[summarized["candidate_id"] == "c1"].iloc[0]
    c2 = summarized[summarized["candidate_id"] == "c2"].iloc[0]
    assert c1["human_exact_hit_count"] == 1
    assert c1["human_minimum_mismatch_count"] == 0
    assert c2["human_one_mismatch_hit_count"] == 1
    ranked = rank_post_human_candidates(summarized)
    assert "post_human_score" in ranked
    assert ranked.loc[ranked["candidate_id"] == "c1", "decision"].iloc[0] == (
        "exclude_or_expert_review"
    )


def test_empty_and_bulge_enabled_outputs_are_supported(tmp_path):
    empty_path = tmp_path / "empty.tsv"
    empty_path.write_text("", encoding="utf-8")
    empty = read_cas_offinder_output(empty_path)
    summarized = summarize_cas_offinder_hits(_candidates(), empty)
    assert summarized["human_total_predicted_hits"].eq(0).all()

    bulge_path = tmp_path / "bulge.tsv"
    bulge_path.write_text(
        "TACGATGCTAACCGGTTAACNNN\tDNA\tTACGATGCTAACCGGTTAACAGG\tchr1\t10\t+\t1\t1\n",
        encoding="utf-8",
    )
    bulge = read_cas_offinder_output(bulge_path)
    assert bulge.iloc[0]["bulge_type"] == "DNA"
    assert bulge.iloc[0]["bulge_size"] == 1
