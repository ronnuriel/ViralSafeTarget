from __future__ import annotations

from collections.abc import Iterable, Mapping
from pathlib import Path

import pandas as pd

from .crispr import reverse_complement


def _hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError("Hamming distance requires equal-length strings")
    return sum(x != y for x, y in zip(a, b, strict=True))


def _enumerate_spcas9_sites(records: Mapping[str, str]) -> Iterable[dict]:
    for seqid, raw_sequence in records.items():
        sequence = raw_sequence.upper().replace("-", "")
        for i in range(0, len(sequence) - 22):
            protospacer = sequence[i : i + 20]
            pam = sequence[i + 20 : i + 23]
            if set(protospacer + pam) <= set("ACGT") and pam[1:] == "GG":
                yield {
                    "seqid": seqid,
                    "start_1based": i + 1,
                    "end_1based": i + 20,
                    "strand": "+",
                    "guide": protospacer,
                    "pam": pam,
                }

            pam_on_plus = sequence[i : i + 3]
            genomic_target = sequence[i + 3 : i + 23]
            if set(pam_on_plus + genomic_target) <= set("ACGT") and pam_on_plus[:2] == "CC":
                yield {
                    "seqid": seqid,
                    "start_1based": i + 4,
                    "end_1based": i + 23,
                    "strand": "-",
                    "guide": reverse_complement(genomic_target),
                    "pam": reverse_complement(pam_on_plus),
                }


def screen_against_small_fasta(
    candidates: pd.DataFrame,
    host_records: Mapping[str, str],
    max_mismatches: int = 3,
    max_total_host_bases: int = 5_000_000,
) -> pd.DataFrame:
    """
    Exhaustively screen candidates against a SMALL host FASTA.

    This educational implementation is intentionally guarded. For GRCh38 use
    Cas-OFFinder/CRISPRitz or another validated genome-scale engine.
    """
    total_bases = sum(len(seq.replace("-", "")) for seq in host_records.values())
    if total_bases > max_total_host_bases:
        raise ValueError(
            f"Host FASTA has {total_bases:,} bases; use Cas-OFFinder for genome-scale screening."
        )

    host_sites = list(_enumerate_spcas9_sites(host_records))
    summaries: list[dict] = []
    for _, candidate in candidates.iterrows():
        guide = str(candidate["guide_sequence"])
        hits = []
        for site in host_sites:
            mismatches = _hamming(guide, site["guide"])
            if mismatches <= max_mismatches:
                hits.append({**site, "mismatches": mismatches})
        hits.sort(key=lambda x: (x["mismatches"], x["seqid"], x["start_1based"]))
        best = hits[0] if hits else None
        summaries.append({
            "candidate_id": candidate["candidate_id"],
            "host_exact_matches": sum(hit["mismatches"] == 0 for hit in hits),
            f"host_matches_le_{max_mismatches}_mismatches": len(hits),
            "host_min_mismatches": best["mismatches"] if best else None,
            "host_best_location": (
                f"{best['seqid']}:{best['start_1based']}-{best['end_1based']}({best['strand']})"
                if best else ""
            ),
            "host_best_pam": best["pam"] if best else "",
        })
    return candidates.merge(pd.DataFrame(summaries), on="candidate_id", how="left")


def write_cas_offinder_input(
    candidates: pd.DataFrame,
    human_fasta_directory: str | Path,
    output_path: str | Path,
    max_mismatches: int = 3,
) -> Path:
    """Write a Cas-OFFinder input file for 20-nt SpCas9 guides and NGG PAM."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        str(Path(human_fasta_directory).resolve()),
        ("N" * 21) + "GG",  # 20 guide positions + N from the NGG PAM
    ]
    for _, row in candidates.iterrows():
        lines.append(f"{row['guide_sequence']}NNN {max_mismatches} {row['candidate_id']}")
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_path


def read_cas_offinder_output(path: str | Path) -> pd.DataFrame:
    """Parse Cas-OFFinder tab-separated output (v2/v3 common columns)."""
    columns = [
        "candidate_id", "bulge_type", "query", "off_target_sequence", "chromosome",
        "location_0based", "direction", "mismatches", "bulge_size",
    ]
    df = pd.read_csv(path, sep="\t", header=None, names=columns, comment="#")
    return df


def summarize_cas_offinder_hits(
    candidates: pd.DataFrame,
    hits: pd.DataFrame,
    *,
    max_mismatches: int = 3,
) -> pd.DataFrame:
    """Merge Cas-OFFinder hits into one transparent summary per candidate.

    The function assumes query IDs in the Cas-OFFinder input were the
    ``candidate_id`` values emitted by ViralSafeTarget.
    """
    if candidates.empty:
        return candidates.copy()
    if hits.empty:
        out = candidates.copy()
        out["host_exact_matches"] = 0
        out[f"host_matches_le_{max_mismatches}_mismatches"] = 0
        out["host_min_mismatches"] = pd.NA
        out["host_best_location"] = ""
        return out

    filtered = hits[hits["mismatches"].astype(int) <= max_mismatches].copy()
    rows: list[dict] = []
    for candidate_id, group in filtered.groupby("candidate_id", sort=False):
        group = group.sort_values(["mismatches", "chromosome", "location_0based"])
        best = group.iloc[0]
        rows.append({
            "candidate_id": str(candidate_id),
            "host_exact_matches": int((group["mismatches"].astype(int) == 0).sum()),
            f"host_matches_le_{max_mismatches}_mismatches": int(len(group)),
            "host_min_mismatches": int(best["mismatches"]),
            "host_best_location": (
                f"{best['chromosome']}:{int(best['location_0based']) + 1}({best['direction']})"
            ),
            "host_best_off_target_sequence": str(best["off_target_sequence"]),
            "host_best_bulge_type": str(best["bulge_type"]),
            "host_best_bulge_size": int(best["bulge_size"]),
        })
    summary = pd.DataFrame(rows)
    out = candidates.merge(summary, on="candidate_id", how="left")
    out["host_exact_matches"] = out["host_exact_matches"].fillna(0).astype(int)
    count_col = f"host_matches_le_{max_mismatches}_mismatches"
    out[count_col] = out[count_col].fillna(0).astype(int)
    out["host_best_location"] = out["host_best_location"].fillna("")
    return out
