"""Gene-function and sequence-level protein-disruption analysis.

All repair outcomes are deterministic computational hypotheses.  They do not
predict editing efficiency, repair frequencies, viral viability, delivery,
toxicity, safety, or clinical efficacy.
"""

from __future__ import annotations

import math
from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Align, SeqIO
from Bio.Seq import Seq

TARGET_GENES = ("UL3", "UL10", "UL18", "UL20", "UL36", "UL52", "UL53", "UL19", "UL30")


@dataclass(frozen=True)
class CdsRecord:
    gene_name: str
    accession: str
    start_1based: int
    end_1based: int
    strand: int
    nucleotide_sequence: str
    protein_sequence: str
    protein_id: str
    product: str

    @property
    def protein_length(self) -> int:
        return len(self.protein_sequence)

    @property
    def coding_length(self) -> int:
        return len(self.nucleotide_sequence)

    def cut_offset_0based(self, cut_after_1based: int) -> int:
        """Return coding bases before a genomic cut boundary."""
        if self.strand == 1:
            return cut_after_1based - self.start_1based + 1
        return self.end_1based - cut_after_1based

    def genomic_to_cds_position(self, genomic_position_1based: int) -> int:
        if self.strand == 1:
            return genomic_position_1based - self.start_1based + 1
        return self.end_1based - genomic_position_1based + 1


def read_target_cds(
    genbank_path: str | Path,
    genes: Iterable[str] = TARGET_GENES,
) -> dict[str, CdsRecord]:
    record = SeqIO.read(genbank_path, "genbank")
    wanted = set(genes)
    output: dict[str, CdsRecord] = {}
    for feature in record.features:
        if feature.type != "CDS":
            continue
        gene = feature.qualifiers.get("gene", [""])[0]
        if gene not in wanted:
            continue
        nucleotide = str(feature.extract(record.seq)).upper()
        protein = feature.qualifiers.get("translation", [""])[0]
        output[gene] = CdsRecord(
            gene_name=gene,
            accession=record.id,
            start_1based=int(feature.location.start) + 1,
            end_1based=int(feature.location.end),
            strand=int(feature.location.strand or 1),
            nucleotide_sequence=nucleotide,
            protein_sequence=protein,
            protein_id=feature.qualifiers.get("protein_id", [""])[0],
            product=feature.qualifiers.get("product", [""])[0],
        )
    missing = sorted(wanted - set(output))
    if missing:
        raise ValueError(f"Target CDS records are absent from {genbank_path}: {', '.join(missing)}")
    return output


def select_top_candidates(
    candidates: pd.DataFrame,
    feature_map: pd.DataFrame,
    cds_records: Mapping[str, CdsRecord],
    *,
    top_per_gene: int = 10,
) -> pd.DataFrame:
    """Select ranked candidates that overlap each requested CDS."""
    cds_map = feature_map[
        feature_map["feature_type"].eq("CDS")
        & feature_map["gene_name"].isin(cds_records)
        & feature_map["overlap_bp"].gt(0)
    ][["candidate_id", "gene_name", "feature_id", "overlap_bp"]].copy()
    joined = cds_map.merge(candidates, on="candidate_id", how="inner", validate="many_to_one")
    joined = joined[joined["screening_status"].eq("completed")]
    joined = joined.sort_values(
        ["gene_name_x", "post_human_rank", "candidate_id"], kind="mergesort"
    )
    joined["mapped_gene_for_analysis"] = joined["gene_name_x"]
    joined = joined.drop(columns=["gene_name_x"]).rename(columns={"gene_name_y": "source_gene"})
    joined = joined.drop_duplicates(["mapped_gene_for_analysis", "candidate_id"])
    joined["within_gene_rank"] = joined.groupby("mapped_gene_for_analysis").cumcount() + 1
    return joined[joined["within_gene_rank"].le(top_per_gene)].reset_index(drop=True)


def _reference_to_alignment_columns(reference_aligned: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    position = 0
    for column, base in enumerate(reference_aligned):
        if base != "-":
            position += 1
            mapping[position] = column
    return mapping


def _entropy(values: Iterable[str]) -> float:
    clean = [str(value).upper() for value in values if str(value).upper() in "ACGT"]
    if not clean:
        return math.nan
    counts = Counter(clean)
    return float(
        -sum((count / len(clean)) * math.log2(count / len(clean)) for count in counts.values())
    )


def _categorical_entropy(values: Iterable[str]) -> float:
    clean = [str(value) for value in values if str(value) not in {"", "X", "-", "?"}]
    if not clean:
        return math.nan
    counts = Counter(clean)
    return float(
        -sum((count / len(clean)) * math.log2(count / len(clean)) for count in counts.values())
    )


def _aligned_cds_sequences(
    aligned_records: Mapping[str, str],
    reference_id: str,
    cds: CdsRecord,
) -> dict[str, str]:
    columns = _reference_to_alignment_columns(aligned_records[reference_id])
    genomic_positions = range(cds.start_1based, cds.end_1based + 1)
    if cds.strand == -1:
        genomic_positions = reversed(list(genomic_positions))
    positions = list(genomic_positions)
    output: dict[str, str] = {}
    for accession, sequence in aligned_records.items():
        bases = [sequence[columns[position]].upper() for position in positions]
        if cds.strand == -1:
            bases = [str(Seq(base).complement()) if base in "ACGT" else base for base in bases]
        output[accession] = "".join(bases)
    return output


def _protein_identity(reference: str, ortholog: str) -> tuple[float, float]:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(reference, ortholog)[0]
    matches = 0
    aligned = 0
    for (ref_start, ref_end), (ortho_start, ortho_end) in zip(
        alignment.aligned[0], alignment.aligned[1], strict=True
    ):
        length = min(ref_end - ref_start, ortho_end - ortho_start)
        matches += sum(
            reference[ref_start + offset] == ortholog[ortho_start + offset]
            for offset in range(length)
        )
        aligned += length
    identity = matches / aligned if aligned else math.nan
    coverage = aligned / len(reference) if reference else math.nan
    return identity, coverage


def _pairwise_dnds(reference_cds: str, comparison_cds: str) -> tuple[float, float, float, str]:
    try:
        from Bio.codonalign.codonseq import CodonSeq, cal_dn_ds

        ref_codons: list[str] = []
        other_codons: list[str] = []
        for index in range(0, min(len(reference_cds), len(comparison_cds)) - 2, 3):
            ref = reference_cds[index : index + 3]
            other = comparison_cds[index : index + 3]
            if set(ref + other) <= set("ACGT") and "*" not in (
                str(Seq(ref).translate()),
                str(Seq(other).translate()),
            ):
                ref_codons.append(ref)
                other_codons.append(other)
        if len(ref_codons) < 30:
            return math.nan, math.nan, math.nan, "insufficient_valid_codons"
        dn, ds = cal_dn_ds(CodonSeq("".join(ref_codons)), CodonSeq("".join(other_codons)))
        ratio = dn / ds if ds and np.isfinite(ds) else math.nan
        return float(dn), float(ds), float(ratio), "computed_NG86"
    except (ValueError, ZeroDivisionError, RuntimeError):
        return math.nan, math.nan, math.nan, "not_estimable"


def compute_gene_evolution(
    aligned_records: Mapping[str, str],
    reference_id: str,
    cds_records: Mapping[str, CdsRecord],
    hsv1_cds: Mapping[str, CdsRecord],
) -> tuple[pd.DataFrame, dict[str, dict[str, str]]]:
    rows: list[dict[str, object]] = []
    aligned_by_gene: dict[str, dict[str, str]] = {}
    for gene, cds in cds_records.items():
        aligned = _aligned_cds_sequences(aligned_records, reference_id, cds)
        aligned_by_gene[gene] = aligned
        ref = aligned[reference_id]
        nucleotide_entropies = [
            _entropy(sequence[index] for sequence in aligned.values()) for index in range(len(ref))
        ]
        amino_acid_conservation: list[float] = []
        for codon_index, reference_amino_acid in enumerate(cds.protein_sequence):
            start = codon_index * 3
            translated = [
                str(Seq(sequence[start : start + 3]).translate())
                for sequence in aligned.values()
                if len(sequence[start : start + 3]) == 3
                and set(sequence[start : start + 3]) <= set("ACGT")
            ]
            if translated:
                amino_acid_conservation.append(
                    sum(value == reference_amino_acid for value in translated) / len(translated)
                )
        dn_values: list[float] = []
        ds_values: list[float] = []
        ratio_values: list[float] = []
        statuses: Counter[str] = Counter()
        for accession, sequence in aligned.items():
            if accession == reference_id:
                continue
            dn, ds, ratio, status = _pairwise_dnds(ref, sequence)
            statuses[status] += 1
            if np.isfinite(dn):
                dn_values.append(dn)
            if np.isfinite(ds):
                ds_values.append(ds)
            if np.isfinite(ratio):
                ratio_values.append(ratio)
        identity, coverage = _protein_identity(
            cds.protein_sequence, hsv1_cds[gene].protein_sequence
        )
        rows.append(
            {
                "gene_name": gene,
                "hsv2_strain_count": len(aligned),
                "mean_cds_nucleotide_entropy_bits": float(np.nanmean(nucleotide_entropies)),
                "mean_amino_acid_conservation": (
                    float(np.mean(amino_acid_conservation)) if amino_acid_conservation else pd.NA
                ),
                "invariant_amino_acid_fraction": (
                    float(np.mean(np.asarray(amino_acid_conservation) == 1.0))
                    if amino_acid_conservation
                    else pd.NA
                ),
                "median_pairwise_dN": np.median(dn_values) if dn_values else pd.NA,
                "median_pairwise_dS": np.median(ds_values) if ds_values else pd.NA,
                "median_pairwise_dN_dS": np.median(ratio_values) if ratio_values else pd.NA,
                "dN_dS_estimable_comparison_count": len(ratio_values),
                "dN_dS_method": "NG86 pairwise median; Biopython experimental codonalign",
                "dN_dS_status_counts": ";".join(
                    f"{key}:{value}" for key, value in sorted(statuses.items())
                ),
                "hsv1_ortholog_accession": hsv1_cds[gene].protein_id,
                "hsv1_ortholog_protein_identity": identity,
                "hsv1_ortholog_reference_coverage": coverage,
                "ortholog_comparison_scope": "HSV-2 reference protein versus HSV-1 strain 17",
            }
        )
    return pd.DataFrame(rows), aligned_by_gene


def map_candidates_to_protein(
    selected: pd.DataFrame,
    cds_records: Mapping[str, CdsRecord],
    aligned_records: Mapping[str, str],
    aligned_cds: Mapping[str, Mapping[str, str]],
    domains: pd.DataFrame,
    disorder: pd.DataFrame,
    evolution: pd.DataFrame,
) -> pd.DataFrame:
    columns = _reference_to_alignment_columns(aligned_records[next(iter(aligned_records))])
    rows: list[dict[str, object]] = []
    for _, candidate in selected.iterrows():
        gene = str(candidate["mapped_gene_for_analysis"])
        cds = cds_records[gene]
        cut = int(candidate["cut_position"])
        offset = cds.cut_offset_0based(cut)
        if offset < 0 or offset > cds.coding_length:
            continue
        affected_aa = min(offset // 3 + 1, cds.protein_length)
        codon_start = (affected_aa - 1) * 3
        codon = cds.nucleotide_sequence[codon_start : codon_start + 3]
        amino_acid = cds.protein_sequence[affected_aa - 1]
        guide_start = cds.genomic_to_cds_position(int(candidate["reference_start_1based"]))
        guide_end = cds.genomic_to_cds_position(int(candidate["reference_end_1based"]))
        genomic_positions = range(
            int(candidate["reference_start_1based"]),
            int(candidate["reference_end_1based"]) + 1,
        )
        entropies = [
            _entropy(sequence[columns[position]] for sequence in aligned_records.values())
            for position in genomic_positions
        ]
        codon_variants = [
            sequence[codon_start : codon_start + 3] for sequence in aligned_cds[gene].values()
        ]
        aa_variants = [
            str(Seq(value).translate())
            for value in codon_variants
            if len(value) == 3 and set(value) <= set("ACGT")
        ]
        aa_conservation = (
            sum(value == amino_acid for value in aa_variants) / len(aa_variants)
            if aa_variants
            else math.nan
        )
        gene_domains = domains[domains["gene_name"].eq(gene)]
        cut_domains = gene_domains[
            gene_domains["protein_start_1based"].le(affected_aa)
            & gene_domains["protein_end_1based"].ge(affected_aa)
        ]
        gene_disorder = disorder[disorder["gene_name"].eq(gene)]
        in_disorder = bool(
            (
                gene_disorder["protein_start_1based"].le(affected_aa)
                & gene_disorder["protein_end_1based"].ge(affected_aa)
            ).any()
        )
        structure_status = (
            "predicted_disordered"
            if in_disorder
            else "predicted_structured_complement"
            if not gene_disorder.empty
            else "unknown_no_MobiDB_lite_region_set"
        )
        row = candidate.to_dict()
        row.update(
            {
                "gene_name": gene,
                "cds_start_1based": cds.start_1based,
                "cds_end_1based": cds.end_1based,
                "cds_strand": "+" if cds.strand == 1 else "-",
                "cds_position_start_1based": min(guide_start, guide_end),
                "cds_position_end_1based": max(guide_start, guide_end),
                "cut_cds_offset_0based": offset,
                "cut_phase": offset % 3,
                "codon_index_1based": affected_aa,
                "reference_codon": codon,
                "reference_amino_acid": amino_acid,
                "amino_acid_position_1based": affected_aa,
                "protein_length_aa": cds.protein_length,
                "relative_protein_position": affected_aa / cds.protein_length,
                "protein_id": cds.protein_id,
                "mean_guide_nucleotide_entropy_bits": float(np.nanmean(entropies)),
                "max_guide_nucleotide_entropy_bits": float(np.nanmax(entropies)),
                "cut_codon_amino_acid_conservation": aa_conservation,
                "cut_codon_amino_acid_entropy_bits": _categorical_entropy(aa_variants),
                "cut_domain_accessions": ";".join(cut_domains["interpro_accession"]),
                "cut_domain_names": ";".join(cut_domains["domain_name"]),
                "predicted_region_class": structure_status,
                "coordinate_convention": (
                    "cut_position is a genomic boundary after the stated plus-strand base; "
                    "cut_cds_offset counts coding bases before that boundary"
                ),
            }
        )
        rows.append(row)
    mapped = pd.DataFrame(rows)
    return mapped.merge(evolution, on="gene_name", how="left", validate="many_to_one")


def _translate_outcome(sequence: str, ambiguity_end_nt: int = 0) -> tuple[object, float, str]:
    complete_length = len(sequence) - (len(sequence) % 3)
    protein = str(Seq(sequence[:complete_length]).translate(to_stop=False))
    for index, amino_acid in enumerate(protein, start=1):
        if amino_acid == "*" and (index * 3) > ambiguity_end_nt:
            return index, float(index), "predictable_stop_found"
    return pd.NA, math.nan, "no_predictable_premature_stop"


def simulate_indels(
    mapped: pd.DataFrame,
    cds_records: Mapping[str, CdsRecord],
    domains: pd.DataFrame,
    sizes: Iterable[int] = range(-10, 11),
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, candidate in mapped.iterrows():
        gene = str(candidate["gene_name"])
        cds = cds_records[gene]
        offset = int(candidate["cut_cds_offset_0based"])
        affected_aa = int(candidate["amino_acid_position_1based"])
        for size in sizes:
            size = int(size)
            if size < 0:
                deletion_end = min(offset + abs(size), len(cds.nucleotide_sequence))
                edited = cds.nucleotide_sequence[:offset] + cds.nucleotide_sequence[deletion_end:]
                ambiguity_end = 0
                event = "deletion"
            elif size > 0:
                edited = (
                    cds.nucleotide_sequence[:offset]
                    + ("N" * size)
                    + cds.nucleotide_sequence[offset:]
                )
                ambiguity_end = offset + size + 3
                event = "insertion_unspecified_sequence"
            else:
                edited = cds.nucleotide_sequence
                ambiguity_end = 0
                event = "reference_no_indel"
            frameshift = bool(size % 3) if size else False
            stop, stop_value, stop_status = _translate_outcome(edited, ambiguity_end)
            premature = bool(np.isfinite(stop_value) and stop_value <= cds.protein_length)
            retained = min(stop_value / cds.protein_length, 1.0) if premature else 1.0
            gene_domains = domains[domains["gene_name"].eq(gene)]
            if frameshift:
                affected = gene_domains[gene_domains["protein_end_1based"].ge(affected_aa)]
            else:
                affected = gene_domains[
                    gene_domains["protein_start_1based"].le(affected_aa)
                    & gene_domains["protein_end_1based"].ge(affected_aa)
                ]
            rows.append(
                {
                    "event_class": "single_guide_indel",
                    "candidate_id": candidate["candidate_id"],
                    "candidate_a": candidate["candidate_id"],
                    "candidate_b": pd.NA,
                    "gene_name": gene,
                    "event": event,
                    "indel_size_bp": size,
                    "repair_model": "size-only; deletion anchored downstream in coding orientation",
                    "inserted_sequence_status": (
                        "unspecified_N_placeholder" if size > 0 else "not_applicable"
                    ),
                    "frameshift": frameshift,
                    "frameshift_probability_for_this_size": float(frameshift),
                    "premature_stop_position_aa": stop if premature else pd.NA,
                    "premature_stop_status": (
                        "sequence-dependent_near_insertion; " + stop_status
                        if size > 0
                        else stop_status
                    ),
                    "retained_protein_fraction": retained,
                    "affected_domain_accessions": ";".join(affected["interpro_accession"]),
                    "affected_domain_names": ";".join(affected["domain_name"]),
                    "deleted_coding_bp": abs(size) if size < 0 else 0,
                    "theoretical_deletion_start_cds_1based": offset + 1 if size < 0 else pd.NA,
                    "theoretical_deletion_end_cds_1based": (
                        min(offset + abs(size), cds.coding_length) if size < 0 else pd.NA
                    ),
                    "limitations": (
                        "No repair-frequency or insertion-sequence model; "
                        "this is not an efficacy prediction."
                    ),
                }
            )
    return pd.DataFrame(rows)


def simulate_paired_deletions(
    pairs: pd.DataFrame,
    candidates: pd.DataFrame,
    feature_map: pd.DataFrame,
    cds_records: Mapping[str, CdsRecord],
    domains: pd.DataFrame,
    genes: Iterable[str] = TARGET_GENES,
) -> pd.DataFrame:
    wanted = set(genes)
    candidate_lookup = candidates.set_index("candidate_id").to_dict("index")
    cds_membership = feature_map[
        feature_map["feature_type"].eq("CDS") & feature_map["gene_name"].isin(wanted)
    ][["candidate_id", "gene_name"]].drop_duplicates()
    memberships = set(map(tuple, cds_membership.itertuples(index=False, name=None)))
    rows: list[dict[str, object]] = []
    for _, pair in pairs[
        pairs["gene_a"].isin(wanted) & pairs["gene_a"].eq(pairs["gene_b"])
    ].iterrows():
        gene = str(pair["gene_a"])
        first_id = str(pair["candidate_a"])
        second_id = str(pair["candidate_b"])
        if (first_id, gene) not in memberships or (second_id, gene) not in memberships:
            continue
        cds = cds_records[gene]
        first_offset = cds.cut_offset_0based(int(candidate_lookup[first_id]["cut_position"]))
        second_offset = cds.cut_offset_0based(int(candidate_lookup[second_id]["cut_position"]))
        start, end = sorted((first_offset, second_offset))
        deleted = end - start
        edited = cds.nucleotide_sequence[:start] + cds.nucleotide_sequence[end:]
        frameshift = bool(deleted % 3)
        stop, stop_value, stop_status = _translate_outcome(edited)
        premature = bool(np.isfinite(stop_value) and stop_value <= cds.protein_length)
        retained = min(stop_value / cds.protein_length, 1.0) if premature else 1.0
        aa_start = start // 3 + 1
        aa_end = max(aa_start, math.ceil(end / 3))
        gene_domains = domains[domains["gene_name"].eq(gene)]
        affected = gene_domains[
            gene_domains["protein_start_1based"].le(aa_end)
            & gene_domains["protein_end_1based"].ge(aa_start)
        ]
        if frameshift:
            affected = gene_domains[gene_domains["protein_end_1based"].ge(aa_start)]
        rows.append(
            {
                "event_class": "paired_guide_theoretical_deletion",
                "candidate_id": f"{first_id}__{second_id}",
                "candidate_a": first_id,
                "candidate_b": second_id,
                "gene_name": gene,
                "event": "paired_deletion",
                "indel_size_bp": -deleted,
                "repair_model": "deterministic cut-to-cut deletion in coding orientation",
                "inserted_sequence_status": "not_applicable",
                "frameshift": frameshift,
                "frameshift_probability_for_this_size": float(frameshift),
                "premature_stop_position_aa": stop if premature else pd.NA,
                "premature_stop_status": stop_status,
                "retained_protein_fraction": retained,
                "affected_domain_accessions": ";".join(affected["interpro_accession"]),
                "affected_domain_names": ";".join(affected["domain_name"]),
                "deleted_coding_bp": deleted,
                "theoretical_deletion_start_cds_1based": start + 1,
                "theoretical_deletion_end_cds_1based": end,
                "limitations": (
                    "Theoretical cut-to-cut deletion only; no joint editing "
                    "or repair-frequency prediction."
                ),
            }
        )
    return pd.DataFrame(rows)


def build_domain_overlap(mapped: pd.DataFrame, domains: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, candidate in mapped.iterrows():
        aa = int(candidate["amino_acid_position_1based"])
        for _, domain in domains[domains["gene_name"].eq(candidate["gene_name"])].iterrows():
            if aa < int(domain["protein_start_1based"]):
                relation = "cut_upstream_of_domain"
            elif aa > int(domain["protein_end_1based"]):
                relation = "cut_downstream_of_domain"
            else:
                relation = "cut_inside_domain"
            rows.append(
                {
                    "candidate_id": candidate["candidate_id"],
                    "gene_name": candidate["gene_name"],
                    "amino_acid_position_1based": aa,
                    "interpro_accession": domain["interpro_accession"],
                    "interpro_entry_type": domain.get("entry_type", pd.NA),
                    "domain_name": domain["domain_name"],
                    "domain_start_1based": domain["protein_start_1based"],
                    "domain_end_1based": domain["protein_end_1based"],
                    "relation_to_cut": relation,
                    "affected_by_local_in_frame_indel": relation == "cut_inside_domain",
                    "potentially_affected_by_frameshift": aa <= int(domain["protein_end_1based"]),
                    "domain_source_url": domain["source_url"],
                }
            )
    return pd.DataFrame(rows)


def score_genes(
    mapped: pd.DataFrame,
    outcomes: pd.DataFrame,
    evidence: pd.DataFrame,
    domains: pd.DataFrame,
    disorder: pd.DataFrame,
    evolution: pd.DataFrame,
) -> pd.DataFrame:
    def evidence_status(frame: pd.DataFrame) -> str:
        if frame.empty:
            return "unknown"
        scored = frame[pd.to_numeric(frame["essentiality_score"], errors="coerce").notna()]
        source = scored if not scored.empty else frame
        calls = {
            str(value).strip()
            for value in source["essentiality_call"].dropna()
            if str(value).strip()
        }
        if len(calls) > 1:
            calls.discard("unknown")
        return ";".join(sorted(calls)) if calls else "unknown"

    rows: list[dict[str, object]] = []
    singles = outcomes[outcomes["event_class"].eq("single_guide_indel")]
    for gene in TARGET_GENES:
        candidates = mapped[mapped["gene_name"].eq(gene)]
        gene_outcomes = singles[singles["gene_name"].eq(gene) & singles["indel_size_bp"].ne(0)]
        hsv2_evidence = evidence[
            evidence["gene_name"].eq(gene) & evidence["virus_type"].eq("HSV-2")
        ]
        hsv1_evidence = evidence[
            evidence["gene_name"].eq(gene) & evidence["virus_type"].eq("HSV-1")
        ]
        direct_hsv2_scores = pd.to_numeric(
            hsv2_evidence.loc[
                hsv2_evidence["evidence_strength"].eq("direct"), "essentiality_score"
            ],
            errors="coerce",
        ).dropna()
        hsv1_scores = pd.to_numeric(hsv1_evidence["essentiality_score"], errors="coerce").dropna()
        frameshift_fraction = pd.to_numeric(
            gene_outcomes["frameshift_probability_for_this_size"], errors="coerce"
        ).mean()
        retained = pd.to_numeric(
            gene_outcomes.loc[gene_outcomes["frameshift"].eq(True), "retained_protein_fraction"],
            errors="coerce",
        ).median()
        domain_fraction = candidates["cut_domain_accessions"].fillna("").ne("").mean()
        disruption_score = (
            0.5 * frameshift_fraction + 0.3 * (1 - retained) + 0.2 * domain_fraction
            if np.isfinite(frameshift_fraction) and np.isfinite(retained)
            else math.nan
        )
        evo = evolution[evolution["gene_name"].eq(gene)].iloc[0]
        coverage_flags = {
            "function_evidence": not evidence[evidence["gene_name"].eq(gene)].empty,
            "direct_HSV2_essentiality": not direct_hsv2_scores.empty,
            "HSV1_ortholog_evidence": not hsv1_evidence.empty,
            "domain_annotation": not domains[domains["gene_name"].eq(gene)].empty,
            "disorder_prediction": not disorder[disorder["gene_name"].eq(gene)].empty,
            "ortholog_comparison": pd.notna(evo["hsv1_ortholog_protein_identity"]),
            "dN_dS": pd.notna(evo["median_pairwise_dN_dS"]),
        }
        rows.append(
            {
                "gene_name": gene,
                "top_candidate_count": len(candidates),
                "sequence_targetability_score": pd.to_numeric(
                    candidates["post_human_score"], errors="coerce"
                ).median(),
                "sequence_targetability_score_definition": (
                    "median post-human score of top mapped candidates"
                ),
                "hsv2_evidence_based_essentiality_score": (
                    direct_hsv2_scores.mean() if not direct_hsv2_scores.empty else pd.NA
                ),
                "evidence_based_essentiality_score": (
                    direct_hsv2_scores.mean() if not direct_hsv2_scores.empty else pd.NA
                ),
                "evidence_based_essentiality_scope": "direct HSV-2 evidence only",
                "hsv2_essentiality_status": (evidence_status(hsv2_evidence)),
                "hsv1_ortholog_essentiality_score": (
                    hsv1_scores.mean() if not hsv1_scores.empty else pd.NA
                ),
                "hsv1_essentiality_status": (evidence_status(hsv1_evidence)),
                "predicted_protein_disruption_score": disruption_score,
                "predicted_disruption_score_definition": (
                    "0.5*frameshift fraction + 0.3*(1-median retained fraction) + "
                    "0.2*fraction of cuts inside InterPro entries"
                ),
                "evidence_coverage_score": sum(coverage_flags.values()) / len(coverage_flags),
                "evidence_coverage_components": ";".join(
                    f"{key}={int(value)}" for key, value in coverage_flags.items()
                ),
            }
        )
    return pd.DataFrame(rows)
