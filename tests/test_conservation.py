from pathlib import Path

from viral_safe_target import conservation_profile, read_fasta

ROOT = Path(__file__).resolve().parents[1]


def test_profile_length():
    records = read_fasta(ROOT / "data/demo/virus_aligned.fasta")
    profile = conservation_profile(records)
    assert len(profile) == len(next(iter(records.values())))
    assert profile["identity_all_sequences"].between(0, 1).all()
