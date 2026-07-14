from __future__ import annotations

from collections import OrderedDict
from collections.abc import Mapping
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_fasta(path: str | Path) -> OrderedDict[str, str]:
    """Read a FASTA file and return ordered {record_id: uppercase_sequence}."""
    records: OrderedDict[str, str] = OrderedDict()
    for record in SeqIO.parse(str(path), "fasta"):
        if record.id in records:
            raise ValueError(f"Duplicate FASTA id: {record.id}")
        records[record.id] = str(record.seq).upper()
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def write_fasta(records: Mapping[str, str], path: str | Path) -> None:
    """Write a mapping of record IDs to sequences in FASTA format."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    seq_records = [
        SeqRecord(Seq(seq), id=record_id, description="") for record_id, seq in records.items()
    ]
    SeqIO.write(seq_records, str(path), "fasta")


def require_aligned(records: Mapping[str, str]) -> int:
    """Validate that all sequences have the same alignment length."""
    lengths = {len(seq) for seq in records.values()}
    if len(lengths) != 1:
        raise ValueError(f"Sequences are not aligned; observed lengths: {sorted(lengths)}")
    return next(iter(lengths))
