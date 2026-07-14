"""Observed-variation and sequence-counterfactual escape robustness.

Escape barriers reported here count sequence changes needed to remove exact
protospacer/PAM targets under a configured editor model.  They are not
evolutionary probabilities and do not predict selection, fitness, editing,
viral inactivation, safety, efficacy, treatment, or cure.
"""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping, Sequence
from typing import Any

import pandas as pd

from .config import EditorProfile

IUPAC: dict[str, set[str]] = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},
    "Y": {"C", "T"},
    "S": {"G", "C"},
    "W": {"A", "T"},
    "K": {"G", "T"},
    "M": {"A", "C"},
    "B": {"C", "G", "T"},
    "D": {"A", "G", "T"},
    "H": {"A", "C", "T"},
    "V": {"A", "C", "G"},
    "N": {"A", "C", "G", "T"},
}
COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _target_parts(candidate: Mapping[str, Any], editor: EditorProfile) -> tuple[str, str]:
    guide = str(candidate["guide_sequence"]).upper()
    pam = str(candidate["pam"]).upper()
    if len(guide) != editor.protospacer_length:
        raise ValueError(
            f"Guide {candidate.get('candidate_id', '')!r} has length {len(guide)}; "
            f"expected {editor.protospacer_length}"
        )
    if len(pam) != len(editor.pam_pattern):
        raise ValueError(f"PAM {pam!r} has length {len(pam)}; expected {len(editor.pam_pattern)}")
    return guide, pam


def single_nucleotide_counterfactuals(
    candidate: Mapping[str, Any], editor: EditorProfile
) -> pd.DataFrame:
    """Enumerate all bounded single-nucleotide target-site substitutions."""
    guide, pam = _target_parts(candidate, editor)
    target = guide + pam if editor.pam_orientation == "3prime" else pam + guide
    guide_indexes = (
        set(range(len(guide)))
        if editor.pam_orientation == "3prime"
        else set(range(len(pam), len(target)))
    )
    pam_start = len(guide) if editor.pam_orientation == "3prime" else 0
    site_start = int(candidate["site_start_1based"])
    site_end = int(candidate["site_end_1based"])
    if site_end - site_start + 1 != len(target):
        raise ValueError("Candidate site coordinates do not match guide plus PAM length")
    strand = str(candidate["strand"])
    genomic_positions = (
        list(range(site_start, site_end + 1))
        if strand == "+"
        else list(range(site_end, site_start - 1, -1))
    )
    rows: list[dict[str, Any]] = []
    for index, reference_base in enumerate(target):
        component = "protospacer" if index in guide_indexes else "pam"
        pam_index = index - pam_start if component == "pam" else None
        for alternate in "ACGT":
            if alternate == reference_base:
                continue
            if component == "protospacer":
                disrupts = True
                reason = "exact_protospacer_match_removed"
            else:
                allowed = IUPAC.get(editor.pam_pattern[pam_index].upper(), set())
                disrupts = alternate not in allowed
                reason = (
                    "pam_incompatible"
                    if disrupts
                    else "pam_remains_compatible_under_editor_pattern"
                )
            plus_reference = (
                reference_base if strand == "+" else reference_base.translate(COMPLEMENT)
            )
            plus_alternate = alternate if strand == "+" else alternate.translate(COMPLEMENT)
            rows.append(
                {
                    "candidate_id": candidate.get("candidate_id"),
                    "component": component,
                    "target_position_1based": index + 1,
                    "component_position_1based": (
                        index + 1
                        if component == "protospacer" and editor.pam_orientation == "3prime"
                        else index - len(pam) + 1
                        if component == "protospacer"
                        else pam_index + 1
                    ),
                    "genomic_position_1based": genomic_positions[index],
                    "reference_base_target_strand": reference_base,
                    "alternate_base_target_strand": alternate,
                    "reference_base_plus_strand": plus_reference,
                    "alternate_base_plus_strand": plus_alternate,
                    "disrupts_exact_target": disrupts,
                    "classification": reason,
                    "counterfactual_scope": (
                        "Single-nucleotide sequence counterfactual; no mutation rate or "
                        "evolutionary probability is assigned."
                    ),
                }
            )
    return pd.DataFrame(rows)


def summarize_guide_escape(
    candidates: pd.DataFrame,
    editor: EditorProfile,
    heldout: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarize discovery/held-out support and exact-target counterfactuals."""
    heldout = heldout if heldout is not None else pd.DataFrame()
    heldout_lookup = (
        heldout.drop_duplicates("candidate_id").set_index("candidate_id").to_dict("index")
        if not heldout.empty and "candidate_id" in heldout
        else {}
    )
    detail_frames: list[pd.DataFrame] = []
    rows: list[dict[str, Any]] = []
    for _, candidate in candidates.drop_duplicates("candidate_id").iterrows():
        detail = single_nucleotide_counterfactuals(candidate, editor)
        detail_frames.append(detail)
        disrupting = detail[detail["disrupts_exact_target"]]
        disruptive_positions = (
            disrupting.groupby("component")["component_position_1based"]
            .apply(lambda values: ";".join(map(str, sorted(set(map(int, values))))))
            .to_dict()
        )
        held = heldout_lookup.get(str(candidate["candidate_id"]), {})
        held_coverage = held.get(
            "population_exact_pam_compatible_coverage",
            held.get("heldout_exact_pam_compatible_coverage", math.nan),
        )
        rows.append(
            {
                "candidate_id": candidate["candidate_id"],
                "guide_sequence": candidate["guide_sequence"],
                "gene_name": candidate.get("gene_name", pd.NA),
                "discovery_exact_target_coverage": candidate.get(
                    "exact_strain_coverage", candidate.get("virus_site_coverage", pd.NA)
                ),
                "heldout_exact_target_coverage": held_coverage,
                "heldout_population_genome_count": held.get("population_genome_count", pd.NA),
                "single_nt_counterfactual_count": len(detail),
                "exact_target_disrupting_counterfactual_count": len(disrupting),
                "protospacer_disruptive_positions": disruptive_positions.get("protospacer", ""),
                "pam_disruptive_positions": disruptive_positions.get("pam", ""),
                "target_disruptible_by_one_substitution": bool(len(disrupting)),
                "host_total_predicted_hits": candidate.get("human_total_predicted_hits", pd.NA),
                "post_human_score": candidate.get("post_human_score", pd.NA),
                "interpretation": (
                    "Observed exact-target coverage and single-nucleotide exact-match "
                    "counterfactuals are separate sequence-level axes; neither is an "
                    "evolutionary probability."
                ),
            }
        )
    detail_output = pd.concat(detail_frames, ignore_index=True) if detail_frames else pd.DataFrame()
    return pd.DataFrame(rows), detail_output


def _mutation_coverage(
    candidates: pd.DataFrame, editor: EditorProfile
) -> tuple[dict[tuple[int, str], set[str]], dict[str, pd.DataFrame]]:
    mutation_to_guides: dict[tuple[int, str], set[str]] = defaultdict(set)
    details: dict[str, pd.DataFrame] = {}
    for _, candidate in candidates.drop_duplicates("candidate_id").iterrows():
        candidate_id = str(candidate["candidate_id"])
        detail = single_nucleotide_counterfactuals(candidate, editor)
        details[candidate_id] = detail
        for _, mutation in detail[detail["disrupts_exact_target"]].iterrows():
            key = (
                int(mutation["genomic_position_1based"]),
                str(mutation["alternate_base_plus_strand"]),
            )
            mutation_to_guides[key].add(candidate_id)
    return mutation_to_guides, details


def multiplex_escape_barrier(candidates: pd.DataFrame, editor: EditorProfile) -> dict[str, Any]:
    """Return the exact minimum substitution cover for all panel targets."""
    members = candidates.drop_duplicates("candidate_id").copy()
    ids = members["candidate_id"].astype(str).tolist()
    if not ids:
        return {
            "panel_size": 0,
            "minimum_independent_target_disrupting_substitutions": pd.NA,
            "barrier_status": "unknown_empty_panel",
        }
    mutation_to_guides, _ = _mutation_coverage(members, editor)
    index = {candidate_id: position for position, candidate_id in enumerate(ids)}
    all_mask = (1 << len(ids)) - 1
    masks = {
        sum(1 << index[candidate_id] for candidate_id in guide_ids)
        for guide_ids in mutation_to_guides.values()
    }
    best: dict[int, int] = {0: 0}
    for mask in masks:
        for existing, count in list(best.items())[::-1]:
            combined = existing | mask
            best[combined] = min(best.get(combined, len(ids) + 1), count + 1)
    barrier = best.get(all_mask)
    return {
        "panel_size": len(ids),
        "unique_target_site_count": members[["site_start_1based", "site_end_1based", "strand"]]
        .drop_duplicates()
        .shape[0],
        "minimum_independent_target_disrupting_substitutions": (
            barrier if barrier is not None else pd.NA
        ),
        "barrier_status": "computed_exact_set_cover" if barrier is not None else "unknown",
        "barrier_interpretation": (
            "Minimum number of distinct single-nucleotide substitutions whose union removes "
            "all exact guide/PAM targets. This is a sequence-level barrier, not an "
            "evolutionary probability."
        ),
    }


def select_strategy_panels(
    candidates: pd.DataFrame,
    definitions: Sequence[Mapping[str, Any]],
    gene_categories: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Select configured strategy panels without virus-specific core logic."""
    categories = gene_categories if gene_categories is not None else pd.DataFrame()
    working = candidates.drop_duplicates("candidate_id").copy()
    score_column = "post_human_score" if "post_human_score" in working else "pre_human_score"
    rank_column = "post_human_rank" if "post_human_rank" in working else "pre_human_rank"
    if score_column not in working:
        working[score_column] = pd.NA
    if rank_column not in working:
        working[rank_column] = range(1, len(working) + 1)
    if not categories.empty and "gene_name" in categories:
        working = working.merge(categories, on="gene_name", how="left", validate="many_to_one")
    rows: list[pd.DataFrame] = []
    for definition in definitions:
        strategy = str(definition["id"])
        size = int(definition.get("size", 3))
        selected = working.copy()
        genes = [str(value) for value in definition.get("genes", [])]
        if genes:
            selected = selected[selected["gene_name"].astype(str).isin(genes)]
            gene_order = {gene: index for index, gene in enumerate(genes)}
            selected["_gene_order"] = selected["gene_name"].map(gene_order)
        else:
            selected["_gene_order"] = 0
        categories_wanted = [str(value) for value in definition.get("categories", [])]
        if categories_wanted and "primary_category" in selected:
            selected = selected[selected["primary_category"].isin(categories_wanted)]
        selected = selected.sort_values(
            ["_gene_order", rank_column, score_column, "candidate_id"],
            ascending=[True, True, False, True],
            kind="mergesort",
        )
        if bool(definition.get("unique_genes", True)):
            selected = selected.drop_duplicates("gene_name")
        if bool(definition.get("unique_categories", False)) and "primary_category" in selected:
            selected = selected.drop_duplicates("primary_category")
        selected = selected.head(size).copy()
        selected["strategy"] = strategy
        selected["strategy_member_order"] = range(1, len(selected) + 1)
        selected["strategy_rationale"] = str(definition.get("rationale", ""))
        public_columns = [
            "strategy",
            "strategy_member_order",
            "strategy_rationale",
            "candidate_id",
            "guide_sequence",
            "pam",
            "strand",
            "site_start_1based",
            "site_end_1based",
            "gene_name",
            "post_human_rank",
            "post_human_score",
            "human_total_predicted_hits",
            "exact_strain_coverage",
            "primary_category",
            "lifecycle_stage",
            "biological_rationale",
            "key_limitation",
        ]
        rows.append(selected[[column for column in public_columns if column in selected]])
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def summarize_multiplex_strategies(
    members: pd.DataFrame,
    editor: EditorProfile,
    guide_escape: pd.DataFrame,
    guide_virtual_knockout: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Summarize panels while keeping evidence/risk/disruption as separate axes."""
    virtual = (
        guide_virtual_knockout.drop_duplicates("candidate_id")
        if guide_virtual_knockout is not None and not guide_virtual_knockout.empty
        else pd.DataFrame()
    )
    escape_lookup = guide_escape.drop_duplicates("candidate_id")
    rows: list[dict[str, Any]] = []
    comparison: list[dict[str, Any]] = []

    def numeric_stat(frame: pd.DataFrame, column: str, operation: str) -> float:
        if column not in frame:
            return math.nan
        values = pd.to_numeric(frame[column], errors="coerce")
        if not values.notna().any():
            return math.nan
        return float(getattr(values, operation)())

    for strategy, group in members.groupby("strategy", sort=False):
        barrier = multiplex_escape_barrier(group, editor)
        merged = group.merge(
            escape_lookup,
            on="candidate_id",
            how="left",
            suffixes=("", "_escape"),
            validate="one_to_one",
        )
        if not virtual.empty:
            merged = merged.merge(
                virtual[
                    [
                        "candidate_id",
                        "frameshift_hypothesis_fraction",
                        "premature_stop_hypothesis_count",
                    ]
                ],
                on="candidate_id",
                how="left",
                validate="one_to_one",
            )
        row = {
            "strategy": strategy,
            "candidate_ids": ";".join(group["candidate_id"].astype(str)),
            "gene_names": ";".join(group["gene_name"].fillna("unknown").astype(str)),
            **barrier,
            "minimum_discovery_exact_target_coverage": numeric_stat(
                merged, "discovery_exact_target_coverage", "min"
            ),
            "minimum_heldout_exact_target_coverage": numeric_stat(
                merged, "heldout_exact_target_coverage", "min"
            ),
            "maximum_host_predicted_hits": numeric_stat(merged, "host_total_predicted_hits", "max"),
            "mean_frameshift_hypothesis_fraction": numeric_stat(
                merged, "frameshift_hypothesis_fraction", "mean"
            ),
            "evidence_axis": "not_combined_into_escape_barrier",
            "biological_evidence_context": "; ".join(
                sorted(set(group.get("biological_rationale", pd.Series(dtype=str)).dropna()))
            ),
            "biological_evidence_limitations": "; ".join(
                sorted(set(group.get("key_limitation", pd.Series(dtype=str)).dropna()))
            ),
            "host_risk_axis": "reported_separately_not_a_safety_conclusion",
            "disruption_axis": "hypothesis_grid_not_repair_probability",
        }
        rows.append(row)
        comparison.append(
            {
                "strategy": strategy,
                "panel_size": barrier["panel_size"],
                "sequence_escape_barrier": barrier[
                    "minimum_independent_target_disrupting_substitutions"
                ],
                "discovery_support_min": row["minimum_discovery_exact_target_coverage"],
                "heldout_support_min": row["minimum_heldout_exact_target_coverage"],
                "predicted_host_hits_max": row["maximum_host_predicted_hits"],
                "disruption_hypothesis_fraction_mean": row["mean_frameshift_hypothesis_fraction"],
                "combined_therapeutic_score": pd.NA,
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(comparison)
