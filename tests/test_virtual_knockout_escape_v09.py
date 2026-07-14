from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd
import yaml

from viral_safe_target.config import EditorProfile
from viral_safe_target.escape import (
    multiplex_escape_barrier,
    single_nucleotide_counterfactuals,
)
from viral_safe_target.project_workflow import initialize_project, run_project
from viral_safe_target.virtual_knockout import (
    build_cds_models,
    enumerate_indel_hypotheses,
    map_guides_to_cds,
    summarize_virtual_knockout,
)

ROOT = Path(__file__).resolve().parents[1]
EDITOR = EditorProfile("SpCas9", 20, "NGG", "3prime", 3, 3, tested=True)


def _inputs(tmp_path: Path) -> tuple[Path, Path]:
    fasta = tmp_path / "reference.fasta"
    fasta.write_text(">ref\n" + "ATG" * 50 + "\n", encoding="utf-8")
    gff = tmp_path / "reference.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "ref\ttest\tCDS\t1\t120\t.\t+\t0\tID=plus;Name=PLUS\n"
        "ref\ttest\tCDS\t31\t150\t.\t-\t0\tID=minus;Name=MINUS\n",
        encoding="utf-8",
    )
    return fasta, gff


def _candidate(identifier: str, *, strand: str = "+", start: int = 40) -> dict[str, object]:
    if strand == "+":
        site_start, site_end = start, start + 22
        protospacer_start, protospacer_end = start, start + 19
    else:
        site_start, site_end = start, start + 22
        protospacer_start, protospacer_end = start + 3, start + 22
    return {
        "candidate_id": identifier,
        "guide_sequence": "ACGTACGTACGTACGTACGT",
        "pam": "AGG",
        "strand": strand,
        "reference_start_1based": protospacer_start,
        "reference_end_1based": protospacer_end,
        "site_start_1based": site_start,
        "site_end_1based": site_end,
        "exact_strain_coverage": 1.0,
    }


def test_coordinate_handling_reverse_strand_and_overlapping_genes(tmp_path: Path) -> None:
    fasta, gff = _inputs(tmp_path)
    models = build_cds_models(fasta, gff, "ref")
    candidates = pd.DataFrame([_candidate("plus"), _candidate("minus", strand="-")])
    mapped = map_guides_to_cds(candidates, models, EDITOR)
    assert set(mapped[mapped["candidate_id"].eq("plus")]["mapped_gene_name"]) == {
        "PLUS",
        "MINUS",
    }
    minus = mapped[
        mapped["candidate_id"].eq("minus") & mapped["mapped_gene_name"].eq("MINUS")
    ].iloc[0]
    assert minus["cds_strand"] == "-"
    assert int(minus["cut_cds_offset_0based"]) == 105


def test_frameshift_domain_overlap_and_missing_annotation_are_explicit(tmp_path: Path) -> None:
    fasta, gff = _inputs(tmp_path)
    models = build_cds_models(fasta, gff, "ref")
    candidate = pd.DataFrame([_candidate("guide")])
    regions = pd.DataFrame(
        [
            {
                "gene_name": "PLUS",
                "protein_start_1based": 1,
                "protein_end_1based": 40,
                "region_name": "test domain",
                "region_kind": "domain",
            }
        ]
    )
    mapped = map_guides_to_cds(candidate, models, EDITOR, regions)
    hypotheses = enumerate_indel_hypotheses(mapped, range(-3, 4))
    guide, _ = summarize_virtual_knockout(mapped, hypotheses)
    plus = guide[guide["gene_name"].eq("PLUS")].iloc[0]
    minus = guide[guide["gene_name"].eq("MINUS")].iloc[0]
    assert plus["frameshift_hypothesis_fraction"] == 4 / 6
    assert plus["cut_domain_names"] == "test domain"
    assert minus["conserved_region_status"] == "unknown_no_annotation"
    insertions = hypotheses[hypotheses["indel_size_bp"].gt(0)]
    assert insertions["premature_stop_status"].eq("unknown_unspecified_insertion_sequence").all()
    assert insertions["protein_fraction_remaining"].isna().all()

    unmapped = map_guides_to_cds(pd.DataFrame([_candidate("outside", start=170)]), models, EDITOR)
    assert unmapped.iloc[0]["mapping_status"] == "unmapped_or_cut_outside_cds"
    assert pd.isna(unmapped.iloc[0]["amino_acid_position_1based"])


def test_single_guide_counterfactuals_respect_pam_and_reverse_coordinates() -> None:
    plus = single_nucleotide_counterfactuals(_candidate("plus"), EDITOR)
    assert len(plus) == 69
    assert int(plus["disrupts_exact_target"].sum()) == 66
    assert not plus[plus["component"].eq("pam") & plus["component_position_1based"].eq(1)][
        "disrupts_exact_target"
    ].any()
    reverse = single_nucleotide_counterfactuals(_candidate("minus", strand="-"), EDITOR)
    assert reverse.iloc[0]["genomic_position_1based"] == 62
    assert reverse.iloc[-1]["genomic_position_1based"] == 40


def test_multiplex_escape_barrier_accounts_for_shared_and_independent_sites() -> None:
    shared = pd.DataFrame([_candidate("a"), _candidate("b")])
    assert (
        multiplex_escape_barrier(shared, EDITOR)[
            "minimum_independent_target_disrupting_substitutions"
        ]
        == 1
    )
    independent = pd.DataFrame([_candidate("a", start=10), _candidate("b", start=70)])
    assert (
        multiplex_escape_barrier(independent, EDITOR)[
            "minimum_independent_target_disrupting_substitutions"
        ]
        == 2
    )


def test_virtual_analysis_cache_invalidates_when_indel_grid_changes(tmp_path: Path) -> None:
    project = initialize_project(
        tmp_path / "project",
        project_id="cache-virus",
        display_name="Cache virus",
        reference_accession="HSV2_demo_ref",
    )
    root = project.parent
    shutil.copyfile(ROOT / "data/demo/virus_aligned.fasta", root / "data/reference.fasta")
    shutil.copyfile(ROOT / "data/demo/virus_aligned.fasta", root / "data/strains.aligned.fasta")
    shutil.copyfile(ROOT / "data/demo/reference.gff3", root / "data/reference.gff3")
    shutil.copyfile(ROOT / "data/demo/human_mini.fasta", root / "external/host/host.fasta")
    run_project(project)
    state_path = root / "results/workflow_state.json"
    first = json.loads(state_path.read_text(encoding="utf-8"))
    first_signature = first["stages"]["virtual_knockout"]["signature"]
    values = yaml.safe_load(project.read_text(encoding="utf-8"))
    values["analysis"]["indel_min_bp"] = -2
    values["analysis"]["indel_max_bp"] = 2
    project.write_text(yaml.safe_dump(values, sort_keys=False), encoding="utf-8")
    run_project(project)
    second = json.loads(state_path.read_text(encoding="utf-8"))
    second_signature = second["stages"]["virtual_knockout"]["signature"]
    hypotheses = pd.read_csv(root / "results/virtual_knockout_escape/indel_sequence_hypotheses.csv")
    mappings = pd.read_csv(root / "results/virtual_knockout_escape/guide_cds_mapping.csv")
    assert first_signature != second_signature
    assert len(hypotheses) == len(mappings) * 5
