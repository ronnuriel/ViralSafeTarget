from __future__ import annotations

import json
from collections import OrderedDict
from pathlib import Path

import pandas as pd

from viral_safe_target.config import EditorProfile
from viral_safe_target.crispr import reverse_complement
from viral_safe_target.population_validation import (
    candidate_population_validation,
    qc_population_records,
    select_population_accessions,
    summarize_locus_aware_population_validation,
)


def _editor() -> EditorProfile:
    return EditorProfile("SpCas9", 20, "NGG", "3prime", 3, 3, tested=True)


def test_population_selection_and_qc_keep_unknown_evidence_explicit(tmp_path: Path) -> None:
    summary = tmp_path / "summary.jsonl"
    rows = [
        {
            "accession": "a",
            "length": 23,
            "completeness": "PARTIAL",
            "virus": {"lineage": [{"tax_id": 10310}]},
        },
        {
            "accession": "b",
            "length": 23,
            "completeness": "COMPLETE",
            "virus": {"lineage": [{"tax_id": 10310}]},
        },
    ]
    summary.write_text("".join(json.dumps(row) + "\n" for row in rows), encoding="utf-8")
    selected = select_population_accessions(
        summary, tax_id=10310, minimum_length=20, maximum_length=30
    )
    accepted, audit = qc_population_records(
        OrderedDict([("a", "A" * 23), ("b", "A" * 23)]), selected
    )
    assert list(accepted) == ["a"]
    assert audit.set_index("accession").loc["b", "reason"] == "exact_sequence_duplicate"


def test_population_qc_accepts_valid_iupac_within_declared_threshold() -> None:
    selected = pd.DataFrame(
        [
            {"accession": "iupac", "length": 100, "completeness": "PARTIAL"},
            {"accession": "too_ambiguous", "length": 100, "completeness": "PARTIAL"},
            {"accession": "invalid", "length": 100, "completeness": "PARTIAL"},
        ]
    )
    records = OrderedDict(
        [
            ("iupac", "A" * 99 + "R"),
            ("too_ambiguous", "A" * 98 + "RY"),
            ("invalid", "A" * 99 + "?"),
        ]
    )
    accepted, audit = qc_population_records(records, selected, maximum_n_fraction=0.01)
    indexed = audit.set_index("accession")
    assert list(accepted) == ["iupac"]
    assert indexed.loc["iupac", "ambiguous_base_fraction"] == 0.01
    assert indexed.loc["too_ambiguous", "reason"] == "ambiguous_base_fraction_above_threshold"
    assert indexed.loc["invalid", "reason"] == "invalid_sequence_characters"
    assert indexed.loc["invalid", "invalid_characters"] == "?"


def test_population_qc_explicitly_excludes_discovery_accessions() -> None:
    selected = pd.DataFrame([{"accession": "training", "length": 23, "completeness": "COMPLETE"}])
    accepted, audit = qc_population_records(
        {"training": "A" * 23}, selected, excluded_accessions={"training"}
    )
    assert not accepted
    assert audit.iloc[0]["reason"] == "excluded_discovery_genome"


def test_population_validation_detects_plus_and_reverse_pam_sites() -> None:
    guide = "ACGTACGTACGTACGTACGT"
    plus = guide + "TGG"
    reverse_site = reverse_complement(guide + "AGG")
    records = OrderedDict([("plus", plus), ("minus", reverse_site), ("absent", "A" * 23)])
    candidates = pd.DataFrame([{"candidate_id": "c1", "guide_sequence": guide}])
    result = candidate_population_validation(
        candidates,
        records,
        _editor(),
        record_groups={"plus": "COMPLETE", "minus": "PARTIAL", "absent": "PARTIAL"},
    ).iloc[0]
    assert result["population_exact_pam_compatible_genome_count"] == 2
    assert result["population_exact_pam_compatible_coverage"] == 2 / 3
    assert result["population_complete_exact_pam_compatible_coverage"] == 1
    assert result["population_partial_exact_pam_compatible_coverage"] == 0.5
    assert "partial records" in result["population_validation_interpretation"]


def test_locus_aware_validation_uses_only_observable_records_as_denominator() -> None:
    candidates = pd.DataFrame(
        [
            {
                "candidate_id": "c1",
                "guide_sequence": "A" * 20,
                "site_start_1based": 101,
                "site_end_1based": 123,
            }
        ]
    )
    alignments = pd.DataFrame(
        [
            {
                "accession": "exact",
                "reference_start_0based": 0,
                "reference_end_0based_exclusive": 200,
            },
            {
                "accession": "variant",
                "reference_start_0based": 50,
                "reference_end_0based_exclusive": 150,
            },
            {
                "accession": "partial_elsewhere",
                "reference_start_0based": 500,
                "reference_end_0based_exclusive": 700,
            },
        ]
    )
    result = summarize_locus_aware_population_validation(
        candidates,
        {"exact": {"A" * 20}, "variant": set(), "partial_elsewhere": set()},
        alignments,
        ["exact", "variant", "partial_elsewhere"],
    ).iloc[0]
    assert result["locus_observable_record_count"] == 2
    assert result["exact_target_in_observable_locus_count"] == 1
    assert result["observable_locus_exact_target_coverage"] == 0.5
    assert result["locus_unresolved_record_count"] == 1
