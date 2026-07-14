"""Generic, deterministic virtual-knockout sequence hypotheses.

The functions in this module describe bounded sequence counterfactuals.  They
do not model DNA-repair frequencies, editing efficiency, viral viability,
delivery, safety, efficacy, treatment, or cure.
"""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
from Bio.Seq import Seq

from .annotations import read_gff3
from .config import EditorProfile
from .disruption import cut_after_1based
from .io_utils import read_fasta


@dataclass(frozen=True)
class GenericCds:
    """A CDS represented as genomic positions in coding (5' to 3') order."""

    gene_name: str
    feature_id: str
    seqid: str
    strand: str
    genomic_positions: tuple[int, ...]
    nucleotide_sequence: str
    product: str = ""

    @property
    def coding_length(self) -> int:
        return len(self.nucleotide_sequence)

    @property
    def protein_length(self) -> int:
        translated = str(
            Seq(self.nucleotide_sequence[: self.coding_length - self.coding_length % 3]).translate(
                to_stop=False
            )
        )
        first_stop = translated.find("*")
        return first_stop if first_stop >= 0 else len(translated)

    def cut_offset(self, cut_after: int) -> int | None:
        """Return coding bases before a plus-strand genomic cut boundary."""
        if not self.genomic_positions:
            return None
        if self.strand == "+":
            return sum(position <= cut_after for position in self.genomic_positions)
        if self.strand == "-":
            return sum(position > cut_after for position in self.genomic_positions)
        return None


def _reference_sequence(path: str | Path, reference_id: str) -> str:
    records = read_fasta(path)
    sequence = next(
        (
            value
            for accession, value in records.items()
            if accession == reference_id or accession.startswith(reference_id)
        ),
        None,
    )
    if sequence is None:
        if len(records) != 1:
            raise KeyError(f"Reference {reference_id!r} is absent from {path}")
        sequence = next(iter(records.values()))
    return sequence.replace("-", "").upper()


def build_cds_models(
    reference_fasta: str | Path,
    annotation_gff: str | Path,
    reference_id: str,
) -> list[GenericCds]:
    """Build generic CDS models from reference-matched FASTA and GFF inputs.

    CDS rows sharing a non-empty gene/name/parent key are joined in coding
    order.  GFF phase is applied at the 5' edge of each segment.
    """
    reference = _reference_sequence(reference_fasta, reference_id)
    features = read_gff3(annotation_gff)
    features = features[
        features["feature_type"].eq("CDS") & features["seqid"].astype(str).eq(reference_id)
    ].copy()
    if features.empty:
        return []
    features["gene_key"] = features.apply(
        lambda row: str(row.get("name") or row.get("parent") or row.get("feature_id") or ""),
        axis=1,
    )
    models: list[GenericCds] = []
    for gene_key, group in features.groupby("gene_key", sort=True, dropna=False):
        if not gene_key:
            continue
        strands = set(group["strand"].astype(str))
        if len(strands) != 1 or next(iter(strands)) not in {"+", "-"}:
            continue
        strand = next(iter(strands))
        ordered = group.sort_values("start", ascending=strand == "+", kind="mergesort")
        positions: list[int] = []
        sequence_parts: list[str] = []
        for _, segment in ordered.iterrows():
            segment_positions = list(range(int(segment["start"]), int(segment["end"]) + 1))
            sequence = reference[int(segment["start"]) - 1 : int(segment["end"])]
            if strand == "-":
                segment_positions.reverse()
                sequence = str(Seq(sequence).reverse_complement())
            phase_text = str(segment.get("phase", "."))
            phase = int(phase_text) if phase_text in {"0", "1", "2"} else 0
            positions.extend(segment_positions[phase:])
            sequence_parts.append(sequence[phase:])
        models.append(
            GenericCds(
                gene_name=str(gene_key),
                feature_id=";".join(sorted(set(group["feature_id"].astype(str)))),
                seqid=reference_id,
                strand=strand,
                genomic_positions=tuple(positions),
                nucleotide_sequence="".join(sequence_parts),
                product=";".join(sorted(set(group["product"].dropna().astype(str)) - {""})),
            )
        )
    return models


def _optional_regions(path: str | Path | None, region_kind: str) -> pd.DataFrame:
    columns = [
        "gene_name",
        "protein_start_1based",
        "protein_end_1based",
        "region_name",
        "region_kind",
    ]
    if not path or not Path(path).is_file():
        return pd.DataFrame(columns=columns)
    frame = pd.read_csv(path, sep="\t")
    rename = {
        "domain_name": "region_name",
        "region_type": "region_name",
    }
    frame = frame.rename(columns=rename)
    required = {"gene_name", "protein_start_1based", "protein_end_1based"}
    if not required <= set(frame):
        return pd.DataFrame(columns=columns)
    if "region_name" not in frame:
        frame["region_name"] = region_kind
    frame["region_kind"] = region_kind
    return frame[columns].copy()


def load_annotation_regions(
    domain_table: str | Path | None = None,
    disorder_table: str | Path | None = None,
    conserved_region_table: str | Path | None = None,
) -> pd.DataFrame:
    """Load optional protein-coordinate annotations without filling gaps."""
    return pd.concat(
        [
            _optional_regions(domain_table, "domain"),
            _optional_regions(disorder_table, "disorder"),
            _optional_regions(conserved_region_table, "conserved_region"),
        ],
        ignore_index=True,
    )


def map_guides_to_cds(
    candidates: pd.DataFrame,
    cds_models: Iterable[GenericCds],
    editor: EditorProfile,
    annotations: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Map every candidate to every CDS containing its cut boundary.

    This deliberately preserves one-to-many mappings for overlapping genes.
    Candidates without a CDS mapping receive one row with unknown coordinates.
    """
    models = list(cds_models)
    annotations = annotations if annotations is not None else pd.DataFrame()
    rows: list[dict[str, Any]] = []
    for _, candidate in candidates.iterrows():
        cut = int(candidate.get("cut_position", cut_after_1based(candidate, editor)))
        mappings: list[dict[str, Any]] = []
        for cds in models:
            offset = cds.cut_offset(cut)
            if offset is None or offset <= 0 or offset >= cds.coding_length:
                continue
            amino_acid = offset // 3 + 1
            regions = (
                annotations[
                    annotations.get("gene_name", pd.Series(dtype=str)).astype(str).eq(cds.gene_name)
                    & pd.to_numeric(annotations.get("protein_start_1based"), errors="coerce").le(
                        amino_acid
                    )
                    & pd.to_numeric(annotations.get("protein_end_1based"), errors="coerce").ge(
                        amino_acid
                    )
                ]
                if not annotations.empty
                else pd.DataFrame()
            )

            def region_labels(kind: str, frame: pd.DataFrame = regions) -> str:
                if frame.empty or not {"region_kind", "region_name"} <= set(frame):
                    return ""
                return ";".join(
                    sorted(
                        set(
                            frame.loc[frame["region_kind"].eq(kind), "region_name"]
                            .dropna()
                            .astype(str)
                        )
                    )
                )

            mappings.append(
                {
                    **candidate.to_dict(),
                    "cut_position_1based_boundary": cut,
                    "mapped_gene_name": cds.gene_name,
                    "mapped_feature_id": cds.feature_id,
                    "cds_strand": cds.strand,
                    "cds_length_bp": cds.coding_length,
                    "protein_length_aa": cds.protein_length,
                    "cut_cds_offset_0based": offset,
                    "amino_acid_position_1based": amino_acid,
                    "relative_protein_position": (
                        amino_acid / cds.protein_length if cds.protein_length else math.nan
                    ),
                    "cut_domain_names": region_labels("domain"),
                    "cut_disorder_regions": region_labels("disorder"),
                    "cut_conserved_region_names": region_labels("conserved_region"),
                    "conserved_region_status": (
                        "annotated_overlap"
                        if region_labels("conserved_region")
                        else "annotated_no_overlap"
                        if not annotations[
                            annotations.get("region_kind", pd.Series(dtype=str)).eq(
                                "conserved_region"
                            )
                            & annotations.get("gene_name", pd.Series(dtype=str))
                            .astype(str)
                            .eq(cds.gene_name)
                        ].empty
                        else "unknown_no_annotation"
                    ),
                    "mapping_status": "mapped_to_cds",
                    "_cds_sequence": cds.nucleotide_sequence,
                }
            )
        if not mappings:
            mappings.append(
                {
                    **candidate.to_dict(),
                    "cut_position_1based_boundary": cut,
                    "mapped_gene_name": pd.NA,
                    "mapped_feature_id": pd.NA,
                    "cds_strand": pd.NA,
                    "cds_length_bp": pd.NA,
                    "protein_length_aa": pd.NA,
                    "cut_cds_offset_0based": pd.NA,
                    "amino_acid_position_1based": pd.NA,
                    "relative_protein_position": pd.NA,
                    "cut_domain_names": pd.NA,
                    "cut_disorder_regions": pd.NA,
                    "cut_conserved_region_names": pd.NA,
                    "conserved_region_status": "unknown_no_cds_mapping",
                    "mapping_status": "unmapped_or_cut_outside_cds",
                    "_cds_sequence": pd.NA,
                }
            )
        rows.extend(mappings)
    return pd.DataFrame(rows)


def _first_stop(sequence: str, *, after_amino_acid: int) -> int | None:
    translated = str(Seq(sequence[: len(sequence) - len(sequence) % 3]).translate(to_stop=False))
    for index, amino_acid in enumerate(translated, start=1):
        if amino_acid == "*" and index >= after_amino_acid:
            return index
    return None


def enumerate_indel_hypotheses(
    mapped: pd.DataFrame,
    sizes: Iterable[int] = range(-10, 11),
) -> pd.DataFrame:
    """Enumerate a bounded grid of size-defined sequence hypotheses."""
    size_grid = sorted(set(int(size) for size in sizes))
    if not size_grid or min(size_grid) > 0 or max(size_grid) < 0:
        raise ValueError("The indel grid must contain bounded deletion and insertion sizes")
    rows: list[dict[str, Any]] = []
    for _, guide in mapped.iterrows():
        mapped_status = guide.get("mapping_status") == "mapped_to_cds"
        sequence = str(guide.get("_cds_sequence", "")) if mapped_status else ""
        offset = int(guide["cut_cds_offset_0based"]) if mapped_status else 0
        original_length = len(sequence)
        cut_aa = int(guide["amino_acid_position_1based"]) if mapped_status else 0
        for size in size_grid:
            event = "reference" if size == 0 else "deletion" if size < 0 else "insertion"
            frameshift: object = (
                bool(abs(size) % 3)
                if mapped_status and size
                else (False if mapped_status else pd.NA)
            )
            edited_length: object = original_length + size if mapped_status else pd.NA
            coding_fraction: object = (
                max(edited_length, 0) / original_length
                if mapped_status and original_length
                else pd.NA
            )
            stop_position: object = pd.NA
            stop_status = "unknown_no_cds_mapping"
            protein_fraction: object = pd.NA
            if mapped_status and size <= 0:
                deleted_end = min(offset + abs(size), original_length)
                edited = sequence[:offset] + sequence[deleted_end:]
                stop = _first_stop(edited, after_amino_acid=cut_aa)
                if stop is not None and stop <= int(guide["protein_length_aa"]):
                    stop_position = stop
                    stop_status = "premature_stop_determinable"
                    protein_fraction = stop / int(guide["protein_length_aa"])
                elif frameshift:
                    stop_status = "no_premature_stop_in_bounded_translation"
                    protein_fraction = pd.NA
                else:
                    stop_status = "no_premature_stop_in_hypothesis"
                    protein_fraction = min(coding_fraction, 1.0)
            elif mapped_status and size > 0:
                stop_status = "unknown_unspecified_insertion_sequence"
            rows.append(
                {
                    "candidate_id": guide.get("candidate_id"),
                    "gene_name": guide.get("mapped_gene_name"),
                    "event": event,
                    "indel_size_bp": size,
                    "frameshift": frameshift,
                    "premature_stop_status": stop_status,
                    "premature_stop_position_aa": stop_position,
                    "coding_sequence_fraction_remaining": coding_fraction,
                    "protein_fraction_remaining": protein_fraction,
                    "affected_domain_names": guide.get("cut_domain_names", pd.NA),
                    "affected_disorder_regions": guide.get("cut_disorder_regions", pd.NA),
                    "conserved_region_status": guide.get("conserved_region_status", pd.NA),
                    "hypothesis_definition": (
                        "size-defined downstream deletion in coding orientation"
                        if size < 0
                        else "size-defined insertion with unspecified sequence"
                        if size > 0
                        else "unmodified reference sequence"
                    ),
                    "interpretation_limit": (
                        "Sequence hypothesis only; not a repair-frequency, viability, safety, "
                        "efficacy, treatment, or cure prediction."
                    ),
                }
            )
    return pd.DataFrame(rows)


def summarize_virtual_knockout(
    mapped: pd.DataFrame,
    hypotheses: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create guide- and gene-level summaries without a therapeutic score."""
    summary_rows: list[dict[str, Any]] = []
    keys = ["candidate_id", "mapped_gene_name"]
    for key, group in mapped.groupby(keys, dropna=False, sort=True):
        candidate_id, gene_name = key
        events = hypotheses[
            hypotheses["candidate_id"].astype(str).eq(str(candidate_id))
            & hypotheses["gene_name"]
            .fillna("")
            .astype(str)
            .eq("" if pd.isna(gene_name) else str(gene_name))
            & hypotheses["indel_size_bp"].ne(0)
        ]
        known_frame = events["frameshift"].dropna()
        deterministic_stop = events["premature_stop_status"].eq("premature_stop_determinable")
        first = group.iloc[0]
        protein_values = pd.to_numeric(events["protein_fraction_remaining"], errors="coerce")
        protein_median = protein_values.median() if protein_values.notna().any() else math.nan
        summary_rows.append(
            {
                "candidate_id": candidate_id,
                "guide_sequence": first.get("guide_sequence", pd.NA),
                "gene_name": gene_name,
                "mapping_status": first.get("mapping_status", pd.NA),
                "cut_position_1based_boundary": first.get("cut_position_1based_boundary", pd.NA),
                "cds_strand": first.get("cds_strand", pd.NA),
                "amino_acid_position_1based": first.get("amino_acid_position_1based", pd.NA),
                "relative_protein_position": first.get("relative_protein_position", pd.NA),
                "hypothesis_count_excluding_reference": len(events),
                "frameshift_hypothesis_fraction": (
                    float(pd.Series(known_frame, dtype=float).mean())
                    if len(known_frame)
                    else math.nan
                ),
                "premature_stop_hypothesis_count": int(deterministic_stop.sum()),
                "median_determinable_protein_fraction_remaining": protein_median,
                "cut_domain_names": first.get("cut_domain_names", pd.NA),
                "cut_disorder_regions": first.get("cut_disorder_regions", pd.NA),
                "conserved_region_status": first.get("conserved_region_status", pd.NA),
                "exact_strain_coverage": first.get("exact_strain_coverage", pd.NA),
                "human_total_predicted_hits": first.get("human_total_predicted_hits", pd.NA),
                "post_human_score": first.get("post_human_score", pd.NA),
                "interpretation": (
                    "Fractions summarize an equally weighted hypothesis grid; they are not "
                    "biological repair probabilities or therapeutic scores."
                ),
            }
        )
    guide_summary = pd.DataFrame(summary_rows)
    gene_rows: list[dict[str, Any]] = []
    mapped_guides = guide_summary[guide_summary["mapping_status"].eq("mapped_to_cds")]
    for gene, group in mapped_guides.groupby("gene_name", sort=True):
        gene_rows.append(
            {
                "gene_name": gene,
                "mapped_guide_count": group["candidate_id"].nunique(),
                "median_relative_protein_position": pd.to_numeric(
                    group["relative_protein_position"], errors="coerce"
                ).median(),
                "median_frameshift_hypothesis_fraction": pd.to_numeric(
                    group["frameshift_hypothesis_fraction"], errors="coerce"
                ).median(),
                "guides_with_determinable_premature_stop": int(
                    group["premature_stop_hypothesis_count"].gt(0).sum()
                ),
                "guides_cutting_annotated_domain": int(
                    group["cut_domain_names"].fillna("").ne("").sum()
                ),
                "guides_cutting_annotated_disorder": int(
                    group["cut_disorder_regions"].fillna("").ne("").sum()
                ),
                "guides_in_annotated_conserved_region": int(
                    group["conserved_region_status"].eq("annotated_overlap").sum()
                ),
                "summary_scope": (
                    "Descriptive aggregation only; no biological importance or therapeutic "
                    "priority is inferred."
                ),
            }
        )
    return guide_summary, pd.DataFrame(gene_rows)
