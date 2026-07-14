#!/usr/bin/env python3
"""Fetch and validate one public GenBank reference record."""

from __future__ import annotations

import argparse
import os
import tempfile
import urllib.parse
import urllib.request
from pathlib import Path

from Bio import SeqIO


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--accession", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""))
    args = parser.parse_args()
    if args.output.is_file():
        records = list(SeqIO.parse(str(args.output), "genbank"))
        if records and records[0].id.startswith(args.accession.split(".")[0]):
            print(f"Validated cached GenBank reference: {args.output}")
            return

    query = {
        "db": "nuccore",
        "id": args.accession,
        "rettype": "gbwithparts",
        "retmode": "text",
        "tool": "ViralSafeTarget",
    }
    if args.email:
        query["email"] = args.email
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urllib.parse.urlencode(
        query
    )
    request = urllib.request.Request(url, headers={"User-Agent": "ViralSafeTarget/0.7"})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(request, timeout=120) as response:
        payload = response.read()
    with tempfile.NamedTemporaryFile(
        mode="wb", dir=args.output.parent, prefix=args.output.name, delete=False
    ) as handle:
        handle.write(payload)
        temporary = Path(handle.name)
    try:
        records = list(SeqIO.parse(str(temporary), "genbank"))
        if len(records) != 1 or not records[0].id.startswith(args.accession.split(".")[0]):
            raise ValueError(f"NCBI response did not contain expected record {args.accession}")
        temporary.replace(args.output)
    finally:
        temporary.unlink(missing_ok=True)
    print(f"Downloaded {args.accession} to {args.output}")


if __name__ == "__main__":
    main()
