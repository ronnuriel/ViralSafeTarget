from pathlib import Path

from viral_safe_target import read_fasta, scan_spcas9_candidates

ROOT = Path(__file__).resolve().parents[1]


def test_spcas9_candidates_are_valid():
    records = read_fasta(ROOT / "data/demo/virus_aligned.fasta")
    candidates = scan_spcas9_candidates(records, "HSV2_demo_ref")
    assert len(candidates) > 0
    assert (candidates["guide_sequence"].str.len() == 20).all()
    assert candidates["pam"].str.endswith("GG").all()
