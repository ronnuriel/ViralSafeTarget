from __future__ import annotations

import hashlib
from collections import Counter
from collections.abc import Mapping

import pandas as pd

from .config import EditorProfile, get_editor, load_config
from .io_utils import require_aligned

DNA_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
IUPAC_BASES = {
    "A": set("A"),
    "C": set("C"),
    "G": set("G"),
    "T": set("T"),
    "N": set("ACGT"),
    "R": set("AG"),
    "Y": set("CT"),
    "S": set("GC"),
    "W": set("AT"),
    "K": set("GT"),
    "M": set("AC"),
    "B": set("CGT"),
    "D": set("AGT"),
    "H": set("ACT"),
    "V": set("ACG"),
}


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(DNA_COMPLEMENT)[::-1]


def _matches_iupac(sequence: str, pattern: str) -> bool:
    return len(sequence) == len(pattern) and all(
        base in IUPAC_BASES[code] for base, code in zip(sequence, pattern, strict=True)
    )


def stable_candidate_id(
    reference_accession: str,
    editor_name: str,
    reference_start_1based: int,
    strand: str,
    guide_sequence: str,
    pam: str,
) -> str:
    """Return a stable content-derived identifier, independent of row ordering."""
    identity = "|".join(
        [
            reference_accession,
            editor_name,
            str(reference_start_1based),
            strand,
            guide_sequence.upper(),
            pam.upper(),
        ]
    )
    return "VST-" + hashlib.sha256(identity.encode("utf-8")).hexdigest()[:16]


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
) -> tuple[float, int, list[str]]:
    exact = 0
    accessions: list[str] = []
    for accession, sequence in records.items():
        observed = _extract_by_reference_positions(
            sequence, ref_to_alignment, start_1based, end_1based
        )
        if observed == reference_site:
            exact += 1
            accessions.append(accession)
    return exact / len(records), exact, accessions


def _window_counts(sequence: str, length: int) -> Counter[str]:
    ungapped = sequence.replace("-", "")
    return Counter(
        ungapped[index : index + length]
        for index in range(max(0, len(ungapped) - length + 1))
        if set(ungapped[index : index + length]) <= set("ACGT")
    )


def _both_strand_count(counts: Mapping[str, int], guide: str) -> int:
    reverse = reverse_complement(guide)
    return int(counts.get(guide, 0)) + (int(counts.get(reverse, 0)) if reverse != guide else 0)


def scan_editor_candidates(
    aligned_records: Mapping[str, str],
    reference_id: str,
    editor: EditorProfile,
    min_site_coverage: float = 0.0,
) -> pd.DataFrame:
    """Scan both strands using a validated editor profile.

    Version 0.3 tests the 3-prime-PAM SpCas9 path. Other profiles can be
    represented in configuration but are rejected until a scanner is tested.
    """
    require_aligned(aligned_records)
    if reference_id not in aligned_records:
        raise KeyError(
            f"Reference id '{reference_id}' not found. Available: {list(aligned_records)}"
        )
    if editor.pam_orientation != "3prime":
        raise NotImplementedError("Only tested 3-prime PAM profiles are scan-enabled in v0.3")

    ref_to_alignment, _, reference = _reference_maps(aligned_records[reference_id])
    guide_length = editor.protospacer_length
    pam_length = len(editor.pam_pattern)
    site_length = guide_length + pam_length
    accession_counts = {
        accession: _window_counts(sequence, guide_length)
        for accession, sequence in aligned_records.items()
    }
    reference_counts = accession_counts[reference_id]
    source_accessions = ";".join(aligned_records)
    candidates: list[dict[str, object]] = []

    for index in range(0, len(reference) - site_length + 1):
        protospacer = reference[index : index + guide_length]
        pam = reference[index + guide_length : index + site_length]
        if set(protospacer + pam) <= set("ACGT") and _matches_iupac(pam, editor.pam_pattern):
            start, end = index + 1, index + site_length
            site = reference[index : index + site_length]
            coverage, exact_count, exact_accessions = _site_coverage(
                aligned_records, ref_to_alignment, start, end, site
            )
            if coverage >= min_site_coverage:
                candidates.append(
                    _candidate_row(
                        reference_id,
                        editor,
                        start,
                        index + guide_length,
                        start,
                        end,
                        "+",
                        protospacer,
                        pam,
                        site,
                        coverage,
                        exact_count,
                        exact_accessions,
                        accession_counts,
                        reference_counts,
                        source_accessions,
                    )
                )

        pam_on_plus = reference[index : index + pam_length]
        genomic_target = reference[index + pam_length : index + site_length]
        pam_as_read_by_editor = reverse_complement(pam_on_plus)
        if set(pam_on_plus + genomic_target) <= set("ACGT") and _matches_iupac(
            pam_as_read_by_editor, editor.pam_pattern
        ):
            site = reference[index : index + site_length]
            coverage, exact_count, exact_accessions = _site_coverage(
                aligned_records, ref_to_alignment, index + 1, index + site_length, site
            )
            if coverage >= min_site_coverage:
                guide = reverse_complement(genomic_target)
                candidates.append(
                    _candidate_row(
                        reference_id,
                        editor,
                        index + pam_length + 1,
                        index + site_length,
                        index + 1,
                        index + site_length,
                        "-",
                        guide,
                        pam_as_read_by_editor,
                        site,
                        coverage,
                        exact_count,
                        exact_accessions,
                        accession_counts,
                        reference_counts,
                        source_accessions,
                    )
                )

    columns = [
        "candidate_id",
        "reference_accession",
        "editor",
        "reference_start_1based",
        "reference_end_1based",
        "site_start_1based",
        "site_end_1based",
        "strand",
        "guide_sequence",
        "pam",
        "reference_site_plus_strand",
        "virus_site_coverage",
        "exact_genome_count",
        "genome_count",
        "reference_viral_occurrence_count",
        "all_viral_occurrence_count",
        "guide_genome_presence_count",
        "is_guide_unique_reference",
        "exact_site_accessions",
        "source_accessions",
        "source_genome_count",
    ]
    frame = pd.DataFrame(candidates, columns=columns)
    if frame.empty:
        return frame
    coordinate_map = (
        frame.groupby("guide_sequence", sort=True)["reference_start_1based"]
        .apply(lambda values: ";".join(str(value) for value in sorted(set(values))))
        .to_dict()
    )
    group_sizes = frame.groupby("guide_sequence")["candidate_id"].transform("size")
    frame["duplicate_guide_group_size"] = group_sizes.astype(int)
    frame["several_coordinates_share_guide"] = group_sizes.gt(1)
    frame["duplicate_guide_coordinates"] = frame["guide_sequence"].map(coordinate_map)
    frame["duplicate_handling"] = "retained_all_coordinates"
    return frame.sort_values(
        ["reference_start_1based", "strand", "candidate_id"], kind="mergesort"
    ).reset_index(drop=True)


def _candidate_row(
    reference_id: str,
    editor: EditorProfile,
    reference_start: int,
    reference_end: int,
    site_start: int,
    site_end: int,
    strand: str,
    guide: str,
    pam: str,
    site: str,
    coverage: float,
    exact_count: int,
    exact_accessions: list[str],
    accession_counts: Mapping[str, Mapping[str, int]],
    reference_counts: Mapping[str, int],
    source_accessions: str,
) -> dict[str, object]:
    occurrences = {
        accession: _both_strand_count(counts, guide)
        for accession, counts in accession_counts.items()
    }
    reference_occurrences = _both_strand_count(reference_counts, guide)
    return {
        "candidate_id": stable_candidate_id(
            reference_id, editor.name, reference_start, strand, guide, pam
        ),
        "reference_accession": reference_id,
        "editor": editor.name,
        "reference_start_1based": reference_start,
        "reference_end_1based": reference_end,
        "site_start_1based": site_start,
        "site_end_1based": site_end,
        "strand": strand,
        "guide_sequence": guide,
        "pam": pam,
        "reference_site_plus_strand": site,
        "virus_site_coverage": coverage,
        "exact_genome_count": exact_count,
        "genome_count": len(accession_counts),
        "reference_viral_occurrence_count": reference_occurrences,
        "all_viral_occurrence_count": sum(occurrences.values()),
        "guide_genome_presence_count": sum(value > 0 for value in occurrences.values()),
        "is_guide_unique_reference": reference_occurrences == 1,
        "exact_site_accessions": ";".join(exact_accessions),
        "source_accessions": source_accessions,
        "source_genome_count": len(accession_counts),
    }


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
    editor = get_editor(load_config())
    return scan_editor_candidates(
        aligned_records, reference_id, editor, min_site_coverage=min_site_coverage
    )
