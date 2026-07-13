#!/usr/bin/env python3
"""Download a small public HSV-2 pilot dataset from NCBI Entrez.

The script intentionally downloads a modest number of complete/high-length
records for reproducible prototyping. It is not a population-representative
sampling strategy.
"""
from __future__ import annotations

import argparse
import os
import time
from io import StringIO
from pathlib import Path

from Bio import Entrez, SeqIO


def _email(value: str | None) -> str:
    result = value or os.environ.get("NCBI_EMAIL", "")
    if not result or "@" not in result:
        raise SystemExit("Provide --email or set NCBI_EMAIL to comply with NCBI Entrez policy.")
    return result


def _fetch_text(db: str, ids: list[str] | str, rettype: str, retmode: str) -> str:
    handle = Entrez.efetch(db=db, id=ids, rettype=rettype, retmode=retmode)
    try:
        return handle.read()
    finally:
        handle.close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--email")
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY"))
    parser.add_argument("--max-genomes", type=int, default=50)
    parser.add_argument("--out-dir", type=Path, default=Path("data/raw/hsv2_entrez"))
    args = parser.parse_args()

    Entrez.email = _email(args.email)
    Entrez.tool = "ViralSafeTarget"
    if args.api_key:
        Entrez.api_key = args.api_key
    args.out_dir.mkdir(parents=True, exist_ok=True)

    reference_accession = "NC_001798.2"
    reference_fasta_text = _fetch_text("nuccore", reference_accession, "fasta", "text")
    reference_gb_text = _fetch_text("nuccore", reference_accession, "gbwithparts", "text")
    (args.out_dir / "hsv2_reference.fasta").write_text(reference_fasta_text, encoding="utf-8")
    (args.out_dir / "hsv2_reference.gb").write_text(reference_gb_text, encoding="utf-8")

    query = (
        '"Human alphaherpesvirus 2"[Organism] '
        "AND 140000:170000[SLEN] "
        'AND ("complete genome"[Title] OR "complete sequence"[Title])'
    )
    search_handle = Entrez.esearch(
        db="nuccore",
        term=query,
        retmax=max(args.max_genomes * 5, 100),
        retmode="xml",
    )
    try:
        search_result = Entrez.read(search_handle)
    finally:
        search_handle.close()
    ids = list(search_result.get("IdList", []))
    if not ids:
        raise SystemExit("NCBI returned no HSV-2 records for the pilot query.")

    records = []
    seen_sequences: set[str] = set()
    batch_size = 50
    delay = 0.11 if args.api_key else 0.34
    for start in range(0, len(ids), batch_size):
        text = _fetch_text("nuccore", ids[start : start + batch_size], "fasta", "text")
        for record in SeqIO.parse(StringIO(text), "fasta"):
            sequence = str(record.seq).upper().replace("-", "")
            if not 140_000 <= len(sequence) <= 170_000:
                continue
            if sequence.count("N") / len(sequence) > 0.01:
                continue
            if sequence in seen_sequences:
                continue
            seen_sequences.add(sequence)
            record.seq = record.seq.__class__(sequence)
            record.id = record.id.split()[0]
            record.name = record.id
            record.description = ""
            records.append(record)
            if len(records) >= args.max_genomes:
                break
        if len(records) >= args.max_genomes:
            break
        time.sleep(delay)

    if len(records) < 2:
        raise SystemExit("Fewer than two usable complete HSV-2 records were downloaded.")
    output = args.out_dir / "hsv2_complete_genomes.fasta"
    SeqIO.write(records, str(output), "fasta")
    (args.out_dir / "accessions.txt").write_text(
        "\n".join(record.id for record in records) + "\n",
        encoding="utf-8",
    )
    print(f"Downloaded {len(records)} public HSV-2 records to {output}")
    print(f"Reference files saved under {args.out_dir}")


if __name__ == "__main__":
    main()
