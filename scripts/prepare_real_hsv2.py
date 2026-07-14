#!/usr/bin/env python3
"""Prepare a deterministic HSV-2 pilot dataset and reference GFF3.

This script is data-wrangling only. It does not design or validate a treatment.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path
from urllib.parse import quote

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def clean_description(value: str) -> str:
    return value.replace("\t", " ").replace("\n", " ").strip()


def feature_name(feature) -> tuple[str, str, str]:
    q = feature.qualifiers
    feature_id = (q.get("locus_tag") or q.get("gene") or q.get("protein_id") or [""])[0]
    name = (q.get("gene") or q.get("locus_tag") or q.get("protein_id") or [""])[0]
    product = (q.get("product") or [""])[0]
    return str(feature_id), str(name), str(product)


def genbank_to_gff3(genbank_path: Path, output_path: Path) -> str:
    records = list(SeqIO.parse(str(genbank_path), "genbank"))
    if not records:
        raise ValueError(f"No GenBank record found in {genbank_path}")
    record = records[0]
    seqid = record.id
    lines = ["##gff-version 3"]
    feature_counter = 0
    allowed = {"gene", "CDS", "mRNA", "ncRNA", "repeat_region", "regulatory", "misc_feature"}
    for feature in record.features:
        if feature.type not in allowed or feature.location is None:
            continue
        # Compound locations are represented by their outer span for this pilot.
        start = int(feature.location.start) + 1
        end = int(feature.location.end)
        if start < 1 or end < start:
            continue
        strand_value = getattr(feature.location, "strand", None)
        strand = "+" if strand_value == 1 else "-" if strand_value == -1 else "."
        feature_counter += 1
        feature_id, name, product = feature_name(feature)
        if not feature_id:
            feature_id = f"feature_{feature_counter}"
        attrs = [f"ID={quote(clean_description(feature_id))}"]
        if name:
            attrs.append(f"Name={quote(clean_description(name))}")
        if product:
            attrs.append(f"product={quote(clean_description(product))}")
        note = (feature.qualifiers.get("note") or [""])[0]
        if note:
            attrs.append(f"Note={quote(clean_description(str(note)))}")
        phase = "."
        if feature.type == "CDS":
            codon_start = (feature.qualifiers.get("codon_start") or ["1"])[0]
            try:
                phase = str((int(codon_start) - 1) % 3)
            except (TypeError, ValueError):
                phase = "0"
        lines.append(
            "\t".join(
                [
                    seqid,
                    "NCBI_GenBank",
                    feature.type,
                    str(start),
                    str(end),
                    ".",
                    strand,
                    phase,
                    ";".join(attrs),
                ]
            )
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return seqid


def load_reference(reference_fasta: Path, preferred_accession: str) -> SeqRecord:
    records = list(SeqIO.parse(str(reference_fasta), "fasta"))
    if not records:
        raise ValueError(f"No FASTA records in {reference_fasta}")
    for record in records:
        if record.id == preferred_accession or record.id.startswith(preferred_accession):
            record.id = preferred_accession
            record.name = preferred_accession
            record.description = "HSV-2 reference NC_001798.2"
            return record
    record = records[0]
    record.id = preferred_accession
    record.name = preferred_accession
    record.description = "HSV-2 reference NC_001798.2"
    return record


def quality_ok(sequence: str, min_length: int, max_length: int, max_n_fraction: float) -> bool:
    seq = sequence.upper()
    if not (min_length <= len(seq) <= max_length):
        return False
    if not seq:
        return False
    return seq.count("N") / len(seq) <= max_n_fraction


def build_sample(
    reference: SeqRecord,
    all_genomes_fasta: Path,
    output_path: Path,
    sample_size: int,
    min_length: int,
    max_length: int,
    max_n_fraction: float,
) -> tuple[int, int, list[dict[str, str]]]:
    seen_hashes = {hashlib.sha256(str(reference.seq).upper().encode()).hexdigest()}
    selected: list[SeqRecord] = [reference]
    considered = 0
    usable: list[SeqRecord] = []
    audit: list[dict[str, str]] = [
        {"accession": reference.id, "decision": "accepted", "rejection_reason": "reference"}
    ]

    for record in SeqIO.parse(str(all_genomes_fasta), "fasta"):
        considered += 1
        if record.id == reference.id or record.id.startswith(reference.id):
            audit.append(
                {
                    "accession": record.id,
                    "decision": "rejected",
                    "rejection_reason": "duplicate_reference_accession",
                }
            )
            continue
        sequence = str(record.seq).upper().replace("-", "")
        if not min_length <= len(sequence) <= max_length:
            audit.append(
                {
                    "accession": record.id,
                    "decision": "rejected",
                    "rejection_reason": "length_outside_configured_range",
                }
            )
            continue
        if not sequence or sequence.count("N") / len(sequence) > max_n_fraction:
            audit.append(
                {
                    "accession": record.id,
                    "decision": "rejected",
                    "rejection_reason": "ambiguous_fraction_above_configured_maximum",
                }
            )
            continue
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        if digest in seen_hashes:
            audit.append(
                {
                    "accession": record.id,
                    "decision": "rejected",
                    "rejection_reason": "exact_sequence_duplicate",
                }
            )
            continue
        seen_hashes.add(digest)
        record.seq = record.seq.__class__(sequence)
        record.id = record.id.split()[0]
        record.name = record.id
        record.description = ""
        usable.append(record)

    # Stable selection independent of NCBI package ordering.
    usable.sort(key=lambda item: item.id)
    selected_usable = usable[: max(0, sample_size - 1)]
    selected.extend(selected_usable)
    selected_ids = {record.id for record in selected_usable}
    for record in usable:
        audit.append(
            {
                "accession": record.id,
                "decision": "accepted" if record.id in selected_ids else "rejected",
                "rejection_reason": ""
                if record.id in selected_ids
                else "outside_deterministic_sample_limit",
            }
        )

    if len(selected) < 2:
        raise ValueError("Fewer than two usable HSV-2 genomes were found after filtering")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(selected, str(output_path), "fasta")
    return considered, len(selected), audit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--reference-genbank", type=Path, required=True)
    parser.add_argument("--all-genomes-fasta", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--output-gff", type=Path, required=True)
    parser.add_argument("--sample-size", type=int, default=100)
    parser.add_argument("--min-length", type=int, default=140_000)
    parser.add_argument("--max-length", type=int, default=170_000)
    parser.add_argument("--max-n-fraction", type=float, default=0.01)
    parser.add_argument("--qc-report", type=Path)
    parser.add_argument("--accessions-used", type=Path)
    args = parser.parse_args()

    reference = load_reference(args.reference_fasta, "NC_001798.2")
    gff_seqid = genbank_to_gff3(args.reference_genbank, args.output_gff)
    if gff_seqid != reference.id:
        print(
            f"WARNING: GenBank seqid={gff_seqid}, FASTA reference id={reference.id}. "
            "The generated GFF uses the GenBank seqid; verify consistency."
        )
    considered, selected, audit = build_sample(
        reference,
        args.all_genomes_fasta,
        args.output_fasta,
        args.sample_size,
        args.min_length,
        args.max_length,
        args.max_n_fraction,
    )
    if args.qc_report:
        args.qc_report.parent.mkdir(parents=True, exist_ok=True)
        with args.qc_report.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=["accession", "decision", "rejection_reason"]
            )
            writer.writeheader()
            writer.writerows(audit)
    if args.accessions_used:
        args.accessions_used.parent.mkdir(parents=True, exist_ok=True)
        args.accessions_used.write_text(
            "\n".join(row["accession"] for row in audit if row["decision"] == "accepted") + "\n",
            encoding="utf-8",
        )
    print(
        f"Considered {considered} downloaded records; "
        f"wrote {selected} genomes to {args.output_fasta}"
    )
    print(f"Wrote reference annotations to {args.output_gff}")


if __name__ == "__main__":
    main()
