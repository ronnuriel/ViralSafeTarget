"""Held-out viral population-panel preparation and exact target validation."""

from __future__ import annotations

import json
import re
from collections import OrderedDict
from collections.abc import Mapping, Sequence
from itertools import product
from pathlib import Path

import pandas as pd

from .config import EditorProfile
from .crispr import IUPAC_BASES, reverse_complement

VALID_DNA_IUPAC = frozenset("ACGTNRYSWKMBDHV")
UNAMBIGUOUS_DNA = frozenset("ACGT")


def select_population_accessions(
    summary_jsonl: str | Path,
    *,
    tax_id: int,
    minimum_length: int,
    maximum_length: int,
) -> pd.DataFrame:
    """Select length-bounded records from an NCBI Datasets JSONL summary."""
    rows: list[dict[str, object]] = []
    with Path(summary_jsonl).open(encoding="utf-8") as handle:
        for line in handle:
            record = json.loads(line)
            lineage = record.get("virus", {}).get("lineage", [])
            tax_ids = {int(item["tax_id"]) for item in lineage if item.get("tax_id")}
            length = int(record.get("length", 0))
            accession = str(record.get("accession", ""))
            if tax_id not in tax_ids or not accession:
                continue
            if not minimum_length <= length <= maximum_length:
                continue
            rows.append(
                {
                    "accession": accession,
                    "length": length,
                    "completeness": str(record.get("completeness", "unknown")),
                    "is_annotated": bool(record.get("is_annotated", False)),
                    "isolate": str(record.get("isolate", {}).get("name", "")),
                    "release_date": str(record.get("release_date", "")),
                    "update_date": str(record.get("update_date", "")),
                    "sequence_hash": str(record.get("nucleotide", {}).get("sequence_hash", "")),
                }
            )
    return (
        pd.DataFrame(rows)
        .sort_values(["accession", "length"], ascending=[True, False], kind="mergesort")
        .drop_duplicates("accession", keep="first")
        .reset_index(drop=True)
    )


def qc_population_records(
    records: Mapping[str, str],
    selected: pd.DataFrame,
    *,
    maximum_n_fraction: float = 0.01,
    excluded_accessions: set[str] | None = None,
) -> tuple[OrderedDict[str, str], pd.DataFrame]:
    """Reject missing, overly ambiguous, malformed, or exact-duplicate records.

    Valid IUPAC ambiguity codes are retained when the total ambiguous-base fraction
    is within the declared threshold. They never count as an exact target match.
    """
    selected_lookup = selected.set_index("accession")
    seen_sequences: dict[str, str] = {}
    accepted: OrderedDict[str, str] = OrderedDict()
    audit: list[dict[str, object]] = []
    excluded_accessions = excluded_accessions or set()
    for accession in selected["accession"].astype(str):
        sequence = str(records.get(accession, "")).upper().replace("-", "")
        invalid_characters = sorted(set(sequence) - VALID_DNA_IUPAC)
        ambiguous_base_count = sum(base not in UNAMBIGUOUS_DNA for base in sequence)
        ambiguous_base_fraction = ambiguous_base_count / len(sequence) if sequence else pd.NA
        decision = "accepted"
        reason = ""
        duplicate_of = ""
        if accession in excluded_accessions:
            decision, reason = "rejected", "excluded_discovery_genome"
        elif not sequence:
            decision, reason = "rejected", "missing_from_downloaded_fasta"
        elif invalid_characters:
            decision, reason = "rejected", "invalid_sequence_characters"
        elif len(sequence) != int(selected_lookup.loc[accession, "length"]):
            decision, reason = "rejected", "downloaded_length_mismatch"
        elif ambiguous_base_fraction > maximum_n_fraction:
            decision, reason = "rejected", "ambiguous_base_fraction_above_threshold"
        elif sequence in seen_sequences:
            decision, reason = "rejected", "exact_sequence_duplicate"
            duplicate_of = seen_sequences[sequence]
        if decision == "accepted":
            seen_sequences[sequence] = accession
            accepted[accession] = sequence
        audit.append(
            {
                "accession": accession,
                "decision": decision,
                "reason": reason,
                "duplicate_of": duplicate_of,
                "sequence_length": len(sequence),
                "n_fraction": sequence.count("N") / len(sequence) if sequence else pd.NA,
                "ambiguous_base_fraction": ambiguous_base_fraction,
                "invalid_characters": "".join(invalid_characters),
                "submitter_completeness": selected_lookup.loc[accession, "completeness"],
            }
        )
    return accepted, pd.DataFrame(audit)


def exact_guide_presence_by_accession(
    candidates: pd.DataFrame,
    records: Mapping[str, str],
    editor: EditorProfile,
) -> dict[str, set[str]]:
    """Return exact protospacer-plus-PAM presence for each population accession."""
    if editor.pam_orientation != "3prime":
        raise NotImplementedError("Population validation currently supports 3-prime PAM editors")
    target_guides = set(candidates["guide_sequence"].astype(str).str.upper())
    pams = [
        "".join(bases)
        for bases in product(*(sorted(IUPAC_BASES[code]) for code in editor.pam_pattern))
    ]
    site_targets: dict[str, set[str]] = {}
    for guide in target_guides:
        for pam in pams:
            forward = guide + pam
            site_targets.setdefault(forward, set()).add(guide)
            site_targets.setdefault(reverse_complement(forward), set()).add(guide)
    site_length = editor.protospacer_length + len(editor.pam_pattern)
    presence: dict[str, set[str]] = {}
    for accession, raw_sequence in records.items():
        sequence = raw_sequence.upper().replace("-", "")
        present: set[str] = set()
        for index in range(max(0, len(sequence) - site_length + 1)):
            matched = site_targets.get(sequence[index : index + site_length])
            if matched:
                present.update(matched)
        presence[accession] = present
    return presence


def candidate_population_validation(
    candidates: pd.DataFrame,
    records: Mapping[str, str],
    editor: EditorProfile,
    *,
    record_groups: Mapping[str, str] | None = None,
) -> pd.DataFrame:
    """Count exact protospacer plus compatible-PAM presence in each held-out genome.

    Records are searched independently of coordinate alignments. For partial records,
    absence is a conservative lower-bound result rather than proof of a mutated locus.
    """
    if editor.pam_orientation != "3prime":
        raise NotImplementedError("Population validation currently supports 3-prime PAM editors")
    target_guides = set(candidates["guide_sequence"].astype(str).str.upper())
    presence_counts = dict.fromkeys(target_guides, 0)
    groups = {
        accession: str((record_groups or {}).get(accession, "all")).strip().lower()
        for accession in records
    }
    group_names = sorted(set(groups.values())) if record_groups else []
    group_sizes = {group: sum(value == group for value in groups.values()) for group in group_names}
    presence_counts_by_group = {group: dict.fromkeys(target_guides, 0) for group in group_names}
    presence_by_accession = exact_guide_presence_by_accession(candidates, records, editor)
    for accession, present in presence_by_accession.items():
        for guide in present:
            presence_counts[guide] += 1
            if record_groups:
                presence_counts_by_group[groups[accession]][guide] += 1
    output = candidates[["candidate_id", "guide_sequence"]].copy()
    output["population_exact_pam_compatible_genome_count"] = (
        output["guide_sequence"].astype(str).str.upper().map(presence_counts).fillna(0).astype(int)
    )
    output["population_genome_count"] = len(records)
    output["population_exact_pam_compatible_coverage"] = (
        output["population_exact_pam_compatible_genome_count"] / len(records) if records else pd.NA
    )
    for group in group_names:
        safe_group = re.sub(r"[^a-z0-9]+", "_", group).strip("_") or "unknown"
        count_column = f"population_{safe_group}_exact_pam_compatible_genome_count"
        total_column = f"population_{safe_group}_genome_count"
        coverage_column = f"population_{safe_group}_exact_pam_compatible_coverage"
        output[count_column] = (
            output["guide_sequence"]
            .astype(str)
            .str.upper()
            .map(presence_counts_by_group[group])
            .fillna(0)
            .astype(int)
        )
        output[total_column] = group_sizes[group]
        output[coverage_column] = (
            output[count_column] / group_sizes[group] if group_sizes[group] else pd.NA
        )
    output["population_validation_interpretation"] = (
        "Exact protospacer plus compatible PAM presence; absence in partial records is a "
        "lower-bound result, not proof of locus mutation."
    )
    return output


def map_population_to_reference(
    records: Mapping[str, str],
    reference_fasta: str | Path,
    *,
    minimum_mapq: int = 20,
    minimum_identity: float = 0.9,
) -> pd.DataFrame:
    """Map population records to a reference and retain high-confidence intervals."""
    try:
        import mappy as mp
    except ImportError as error:  # pragma: no cover - environment-dependent optional tool
        raise RuntimeError(
            "Reference-aware population validation requires the optional 'mappy' package."
        ) from error
    aligner = mp.Aligner(str(reference_fasta), preset="asm5", best_n=20)
    if not aligner:
        raise ValueError(f"Could not build a mappy reference index from {reference_fasta}")
    rows: list[dict[str, object]] = []
    for accession, sequence in records.items():
        for hit in aligner.map(sequence):
            identity = hit.mlen / hit.blen if hit.blen else 0.0
            if hit.mapq < minimum_mapq or identity < minimum_identity:
                continue
            rows.append(
                {
                    "accession": accession,
                    "reference_contig": hit.ctg,
                    "reference_start_0based": hit.r_st,
                    "reference_end_0based_exclusive": hit.r_en,
                    "query_start_0based": hit.q_st,
                    "query_end_0based_exclusive": hit.q_en,
                    "strand": "+" if hit.strand == 1 else "-",
                    "mapq": hit.mapq,
                    "alignment_identity": identity,
                    "alignment_block_length": hit.blen,
                    "is_primary": bool(hit.is_primary),
                }
            )
    return pd.DataFrame(rows)


def summarize_locus_aware_population_validation(
    candidates: pd.DataFrame,
    presence_by_accession: Mapping[str, set[str]],
    alignments: pd.DataFrame,
    accessions: Sequence[str],
) -> pd.DataFrame:
    """Distinguish exact target loss from loci absent in partial population records."""
    required = {"site_start_1based", "site_end_1based", "guide_sequence", "candidate_id"}
    missing = sorted(required - set(candidates.columns))
    if missing:
        raise ValueError(f"Candidate table lacks locus-validation columns: {missing}")
    intervals: dict[str, list[tuple[int, int]]] = {accession: [] for accession in accessions}
    if not alignments.empty:
        for row in alignments.itertuples(index=False):
            intervals.setdefault(str(row.accession), []).append(
                (
                    int(row.reference_start_0based),
                    int(row.reference_end_0based_exclusive),
                )
            )
    rows: list[dict[str, object]] = []
    for candidate in candidates.itertuples(index=False):
        start = int(candidate.site_start_1based) - 1
        end = int(candidate.site_end_1based)
        guide = str(candidate.guide_sequence).upper()
        observable = {
            accession
            for accession in accessions
            if any(left <= start and end <= right for left, right in intervals.get(accession, []))
        }
        exact_anywhere = {
            accession
            for accession in accessions
            if guide in presence_by_accession.get(accession, set())
        }
        exact_observable = observable & exact_anywhere
        denominator = len(observable)
        rows.append(
            {
                "candidate_id": candidate.candidate_id,
                "guide_sequence": guide,
                "population_record_count": len(accessions),
                "locus_observable_record_count": denominator,
                "exact_target_in_observable_locus_count": len(exact_observable),
                "observable_locus_exact_target_coverage": (
                    len(exact_observable) / denominator if denominator else pd.NA
                ),
                "observable_locus_without_exact_target_count": denominator - len(exact_observable),
                "locus_unresolved_record_count": len(accessions) - denominator,
                "exact_target_anywhere_count": len(exact_anywhere),
                "locus_validation_interpretation": (
                    "Reference interval observability is inferred from high-quality whole-record "
                    "mapping. Exact-target attribution is strongest for reference-unique guides; "
                    "this is population-genomic validation, not editing or safety evidence."
                ),
            }
        )
    return pd.DataFrame(rows)
