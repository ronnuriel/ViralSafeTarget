from pathlib import Path

from viral_safe_target import (
    annotate_candidates,
    read_fasta,
    read_gff3,
    scan_spcas9_candidates,
)

ROOT = Path(__file__).resolve().parents[1]


def test_annotation_mapping():
    records = read_fasta(ROOT / "data/demo/virus_aligned.fasta")
    candidates = scan_spcas9_candidates(records, "HSV2_demo_ref")
    features = read_gff3(ROOT / "data/demo/reference.gff3")
    annotated = annotate_candidates(candidates, features, seqid="HSV2_demo_ref")
    assert "gene_name" in annotated.columns
    assert annotated["gene_name"].ne("").any()
