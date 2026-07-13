from __future__ import annotations

from collections.abc import Mapping

import pandas as pd

from .io_utils import require_aligned

DNA_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(DNA_COMPLEMENT)[::-1]


def _reference_maps(reference_aligned: str) -> tuple[dict[int, int], dict[int, int], str]:
    """Return ref_position->alignment_index, alignment_index->ref_position, ungapped reference."""
    ref_to_alignment: dict[int, int] = {}
    alignment_to_ref: dict[int, int] = {}
    ref_position = 0
    ungapped: list[str] = []
    for alignment_index, base in enumerate(reference_aligned):
        if base != "-":
            ref_position += 1
            ref_to_alignment[ref_position] = alignment_index
            alignment_to_ref[alignment_index] = ref_position
            ungapped.append(base)
    return ref_to_alignment, alignment_to_ref, "".join(ungapped)


def _extract_by_reference_positions(
    aligned_sequence: str,
    ref_to_alignment: dict[int, int],
    start_1based: int,
    end_1based: int,
) -> str:
    chars: list[str] = []
    for ref_pos in range(start_1based, end_1based + 1):
        alignment_idx = ref_to_alignment[ref_pos]
        chars.append(aligned_sequence[alignment_idx])
    return "".join(chars)


def _site_coverage(
    records: Mapping[str, str],
    ref_to_alignment: dict[int, int],
    start_1based: int,
    end_1based: int,
    reference_site: str,
) -> tuple[float, int]:
    exact = 0
    for sequence in records.values():
        observed = _extract_by_reference_positions(
            sequence, ref_to_alignment, start_1based, end_1based
        )
        if observed == reference_site:
            exact += 1
    return exact / len(records), exact


def scan_spcas9_candidates(
    aligned_records: Mapping[str, str],
    reference_id: str,
    min_site_coverage: float = 0.0,
) -> pd.DataFrame:
    """
    Scan an aligned viral genome collection for SpCas9 sites on both strands.

    Candidate coordinates are 1-based relative to the ungapped reference.
    Coverage is the fraction of input genomes with the exact 23-nt site.
    """
    require_aligned(aligned_records)
    if reference_id not in aligned_records:
        raise KeyError(
            f"Reference id '{reference_id}' not found. Available: {list(aligned_records)}"
        )

    ref_to_alignment, _, reference = _reference_maps(aligned_records[reference_id])
    candidates: list[dict] = []
    candidate_id = 0

    # Plus-strand sites: 20-nt protospacer followed by NGG.
    for i in range(0, len(reference) - 22):
        protospacer = reference[i : i + 20]
        pam = reference[i + 20 : i + 23]
        if set(protospacer + pam) <= set("ACGT") and pam[1:] == "GG":
            start, end = i + 1, i + 23
            site = reference[i : i + 23]
            coverage, exact_count = _site_coverage(
                aligned_records, ref_to_alignment, start, end, site
            )
            if coverage >= min_site_coverage:
                candidate_id += 1
                candidates.append({
                    "candidate_id": f"T{candidate_id:05d}",
                    "reference_start_1based": start,
                    "reference_end_1based": i + 20,
                    "site_start_1based": start,
                    "site_end_1based": end,
                    "strand": "+",
                    "guide_sequence": protospacer,
                    "pam": pam,
                    "reference_site_plus_strand": site,
                    "virus_site_coverage": coverage,
                    "exact_genome_count": exact_count,
                    "genome_count": len(aligned_records),
                })

    # Minus-strand sites appear as CCN followed by a 20-nt genomic segment.
    for i in range(0, len(reference) - 22):
        pam_on_plus = reference[i : i + 3]
        genomic_target = reference[i + 3 : i + 23]
        if set(pam_on_plus + genomic_target) <= set("ACGT") and pam_on_plus[:2] == "CC":
            site = reference[i : i + 23]
            coverage, exact_count = _site_coverage(
                aligned_records, ref_to_alignment, i + 1, i + 23, site
            )
            if coverage >= min_site_coverage:
                candidate_id += 1
                candidates.append({
                    "candidate_id": f"T{candidate_id:05d}",
                    "reference_start_1based": i + 4,
                    "reference_end_1based": i + 23,
                    "site_start_1based": i + 1,
                    "site_end_1based": i + 23,
                    "strand": "-",
                    "guide_sequence": reverse_complement(genomic_target),
                    "pam": reverse_complement(pam_on_plus),
                    "reference_site_plus_strand": site,
                    "virus_site_coverage": coverage,
                    "exact_genome_count": exact_count,
                    "genome_count": len(aligned_records),
                })

    columns = [
        "candidate_id", "reference_start_1based", "reference_end_1based",
        "site_start_1based", "site_end_1based", "strand", "guide_sequence",
        "pam", "reference_site_plus_strand", "virus_site_coverage",
        "exact_genome_count", "genome_count",
    ]
    return pd.DataFrame(candidates, columns=columns)
