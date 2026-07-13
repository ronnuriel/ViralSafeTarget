from __future__ import annotations

from pathlib import Path
from urllib.parse import unquote

import pandas as pd


def _parse_attributes(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = unquote(value)
        elif " " in item:
            key, value = item.split(" ", 1)
            attrs[key] = value.strip('"')
    return attrs


def read_gff3(path: str | Path) -> pd.DataFrame:
    """Read a GFF3/GFF-like file into a DataFrame."""
    rows: list[dict] = []
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
            attrs = _parse_attributes(attributes)
            rows.append({
                "seqid": seqid,
                "source": source,
                "feature_type": feature_type,
                "start": int(start),
                "end": int(end),
                "score": score,
                "strand": strand,
                "phase": phase,
                "attributes": attrs,
                "feature_id": attrs.get("ID") or attrs.get("gene") or attrs.get("Name") or "",
                "parent": attrs.get("Parent", ""),
                "name": attrs.get("Name") or attrs.get("gene") or attrs.get("locus_tag") or "",
                "product": attrs.get("product", ""),
                "note": attrs.get("Note", ""),
            })
    return pd.DataFrame(rows)


def annotate_candidates(
    candidates: pd.DataFrame,
    features: pd.DataFrame,
    seqid: str | None = None,
) -> pd.DataFrame:
    """Annotate candidates by overlapping GFF features, preferring CDS/gene features."""
    if candidates.empty:
        return candidates.copy()
    relevant = features.copy()
    if seqid is not None and not relevant.empty:
        relevant = relevant[relevant["seqid"] == seqid]

    priority = {"CDS": 0, "gene": 1, "mRNA": 2}
    output_rows = []
    for _, candidate in candidates.iterrows():
        overlaps = relevant[
            (relevant["start"] <= int(candidate["reference_end_1based"]))
            & (relevant["end"] >= int(candidate["reference_start_1based"]))
        ].copy()
        if overlaps.empty:
            selected = None
        else:
            overlaps["_priority"] = overlaps["feature_type"].map(priority).fillna(9)
            selected = overlaps.sort_values(["_priority", "start", "end"]).iloc[0]

        row = candidate.to_dict()
        if selected is None:
            row.update({
                "feature_type": "intergenic_or_unannotated",
                "feature_id": "",
                "gene_name": "",
                "product": "",
                "feature_start": None,
                "feature_end": None,
            })
        else:
            row.update({
                "feature_type": selected["feature_type"],
                "feature_id": selected["feature_id"],
                "gene_name": selected["name"],
                "product": selected["product"],
                "feature_start": int(selected["start"]),
                "feature_end": int(selected["end"]),
            })
        output_rows.append(row)
    return pd.DataFrame(output_rows)
