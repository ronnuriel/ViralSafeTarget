from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from viral_safe_target.config import load_config
from viral_safe_target.discovery import (
    build_bounded_pair_hypotheses,
    build_candidate_feature_map,
    gene_rank_stability,
    rank_genes,
    select_balanced_discovery_panel,
    wilson_lower_bound,
)
from viral_safe_target.discovery_reporting import write_discovery_report
from viral_safe_target.discovery_workflow import merge_batch_results, prepare_batches


def _features() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "seqid": "REF",
                "feature_type": "gene",
                "start": 1,
                "end": 30,
                "feature_id": "gA",
                "name": "geneA",
            },
            {
                "seqid": "REF",
                "feature_type": "gene",
                "start": 20,
                "end": 50,
                "feature_id": "gB",
                "name": "geneB",
            },
            {
                "seqid": "REF",
                "feature_type": "gene",
                "start": 60,
                "end": 90,
                "feature_id": "gC",
                "name": "geneC",
            },
        ]
    )


def _candidates(count: int = 5) -> pd.DataFrame:
    starts = [1, 20, 40, 60, 100][:count]
    rows = []
    for index, start in enumerate(starts):
        rows.append(
            {
                "candidate_id": f"c{index}",
                "reference_accession": "REF",
                "reference_start_1based": start,
                "reference_end_1based": start + 19,
                "strand": "-" if index == 1 else "+",
                "guide_sequence": ("ACGT" * 5)[:-1] + "ACGT"[index % 4],
                "pam": "TGG",
                "pre_human_score": 1 - index / 10,
                "rejection_reasons": "",
                "exact_strain_coverage": 1.0,
                "exact_genome_count": 3,
                "genome_count": 3,
                "exact_site_accessions": "a;b;c",
            }
        )
    return pd.DataFrame(rows)


def _settings() -> dict:
    return load_config("configs/hsv2_pilot.yaml")


def test_mapping_retains_boundaries_overlap_opposite_strand_and_intergenic() -> None:
    mapping = build_candidate_feature_map(_candidates(), _features(), config=_settings())
    assert set(mapping.loc[mapping["candidate_id"].eq("c1"), "gene_name"]) == {"geneA", "geneB"}
    assert mapping.loc[mapping["candidate_id"].eq("c1"), "strand"].eq("-").all()
    boundary = mapping[mapping["candidate_id"].eq("c0")]
    assert set(boundary["gene_name"]) == {"geneA", "geneB"}
    assert boundary.loc[boundary["gene_name"].eq("geneB"), "overlap_bp"].eq(1).all()
    assert mapping.loc[mapping["candidate_id"].eq("c4"), "mapping_status"].eq("intergenic").all()
    assert mapping["overlap_bp"].ge(0).all()


def test_balanced_selection_is_stable_union_and_not_gene_hardcoded() -> None:
    candidates = _candidates()
    mapping = build_candidate_feature_map(candidates, _features(), config=_settings())
    first = select_balanced_discovery_panel(
        candidates, mapping, _features(), top_per_gene=1, global_top=2
    )
    second = select_balanced_discovery_panel(
        candidates, mapping, _features(), top_per_gene=1, global_top=2
    )
    assert first.panel["candidate_id"].tolist() == second.panel["candidate_id"].tolist()
    assert set(first.audit["selection_reason"]) <= {"per_gene_quota", "global_top", "both"}
    assert {"geneA", "geneB", "geneC"} <= set(
        mapping.loc[mapping["candidate_id"].isin(first.panel["candidate_id"]), "gene_name"]
    )
    assert not first.panel["candidate_id"].duplicated().any()
    assert "UL19" not in first.panel.to_csv(index=False)


def test_exhaustive_requires_explicit_confirmation() -> None:
    candidates = _candidates()
    mapping = build_candidate_feature_map(candidates, _features(), config=_settings())
    with pytest.raises(ValueError, match="confirm-exhaustive"):
        select_balanced_discovery_panel(candidates, mapping, _features(), exhaustive=True)


def _screened_candidates() -> tuple[pd.DataFrame, pd.DataFrame]:
    candidates = _candidates(4)
    mapping = build_candidate_feature_map(candidates, _features(), config=_settings())
    candidates["screening_status"] = "completed"
    candidates["post_human_rank"] = [1, 2, 3, 4]
    candidates["post_human_score"] = [0.9, 0.8, 0.7, 0.6]
    candidates["decision"] = "retain_computational_candidate"
    candidates["human_exact_hit_count"] = 0
    candidates["human_one_mismatch_hit_count"] = [0, 0, 1, 0]
    candidates["human_two_mismatch_hit_count"] = 0
    candidates["human_three_mismatch_hit_count"] = 0
    candidates["human_total_predicted_hits"] = [0, 0, 1, 0]
    candidates["cut_position"] = candidates["reference_start_1based"] + 16
    return candidates, mapping


def test_gene_ranking_uses_fraction_uncertainty_and_missing_evidence() -> None:
    candidates, mapping = _screened_candidates()
    genes = rank_genes(candidates, mapping, _features())
    assert "conserved_candidate_fraction" in genes
    assert genes["biological_evidence_status"].eq("biological evidence not supplied").all()
    assert genes["evidence_coverage"].isna().all()
    assert wilson_lower_bound(1, 1) < wilson_lower_bound(10, 10)
    assert genes.loc[genes["screened_candidate_count"].lt(10), "confidence_level"].eq("low").all()


def test_gene_metrics_use_full_eligible_denominator_not_only_screened_panel() -> None:
    candidates, mapping = _screened_candidates()
    screened_panel = candidates.iloc[:2].copy()
    genes = rank_genes(
        screened_panel,
        mapping,
        _features(),
        eligible_candidates=candidates,
    )
    gene_b = genes[genes["gene_name"].eq("geneB")].iloc[0]
    assert gene_b["eligible_candidate_count"] > gene_b["screened_candidate_count"]
    assert gene_b["screening_fraction"] < 1


def test_stability_fields_exist_without_rerunning_external_search() -> None:
    candidates, mapping = _screened_candidates()
    stability = gene_rank_stability(candidates, mapping, _features())
    assert {"rank_at_k10", "rank_at_k25", "rank_at_k50", "stability_warning"} <= set(stability)


def test_batches_resume_expand_duplicate_guides_and_keep_missing_incomplete(tmp_path: Path) -> None:
    candidates = _candidates(3)
    candidates.loc[1, "guide_sequence"] = candidates.loc[0, "guide_sequence"]
    batches = prepare_batches(candidates, tmp_path / "batches", tmp_path, _settings(), batch_size=1)
    assert len(batches) == 2
    first = batches[0]
    query = pd.read_csv(first["manifest_path"])["cas_offinder_query"].iloc[0]
    first["raw_path"].write_text(f"{query}\tchr1\t5\t{query}\t+\t1\n", encoding="utf-8")
    first["status"] = "completed"
    first["status_path"].write_text(
        json.dumps(
            {
                "status": "completed",
                "input_sha256": first["input_sha256"],
                "candidate_manifest_sha256": first["candidate_manifest_sha256"],
            }
        ),
        encoding="utf-8",
    )
    resumed = prepare_batches(candidates, tmp_path / "batches", tmp_path, _settings(), batch_size=1)
    assert resumed[0]["status"] == "completed"
    panel = candidates.copy()
    panel["pre_human_rank"] = range(1, len(panel) + 1)
    hits, post = merge_batch_results(panel, resumed, _settings())
    first_ids = set(pd.read_csv(resumed[0]["manifest_path"])["candidate_id"])
    assert set(hits["candidate_id"]) == first_ids
    assert post.loc[post["candidate_id"].isin(first_ids), "screening_status"].eq("completed").all()
    assert (
        post.loc[~post["candidate_id"].isin(first_ids), "decision"].eq("screening_incomplete").all()
    )


def test_pair_bounds_and_cross_gene_wording() -> None:
    candidates, mapping = _screened_candidates()
    genes = rank_genes(candidates, mapping, _features())
    same, multi, _ = build_bounded_pair_hypotheses(
        candidates, mapping, genes, top_genes=3, candidates_per_gene=3, maximum_pairs=1
    )
    assert len(same) <= 1 and len(multi) <= 1
    if not multi.empty:
        assert "not one physical deletion" in multi.iloc[0]["limitations"]


def test_feature_map_has_one_cut_position_per_candidate() -> None:
    candidates = _candidates()
    mapping = build_candidate_feature_map(candidates, _features(), config=_settings())
    cut_positions = mapping[["candidate_id", "cut_position"]].drop_duplicates()
    assert not cut_positions["candidate_id"].duplicated().any()
    joined = candidates.merge(cut_positions, on="candidate_id", validate="one_to_one")
    assert joined["cut_position"].notna().all()


def test_partial_report_states_answers_are_not_determinable(tmp_path: Path) -> None:
    candidates, mapping = _screened_candidates()
    candidates["screening_status"] = "pending"
    candidates["post_human_rank"] = pd.NA
    candidates["post_human_score"] = pd.NA
    genes = rank_genes(candidates, mapping, _features())
    stability = gene_rank_stability(candidates, mapping, _features())
    output = tmp_path / "report.html"
    answers = write_discovery_report(
        output,
        candidates=candidates,
        feature_map=mapping,
        genes=genes,
        stability=stability,
        deep_panel=candidates.head(0),
        top_per_gene_candidates=candidates.head(0),
        same_pairs=pd.DataFrame(),
        multi_pairs=pd.DataFrame(),
        genes_without_candidates=pd.DataFrame(),
        qc=pd.DataFrame(),
        provenance={"human_assembly": "synthetic", "human_assembly_accession": "SYN"},
        initial_candidate_count=4,
        eligible_candidate_count=4,
    )
    assert answers["another_gene_above_ul30"] is None
    assert "Research use only" in output.read_text(encoding="utf-8")
