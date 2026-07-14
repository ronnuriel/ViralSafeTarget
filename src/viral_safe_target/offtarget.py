"""Small-host screening and Cas-OFFinder exchange workflows."""

from __future__ import annotations

import bisect
import gzip
import json
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path
from typing import Any
from urllib.parse import unquote

import pandas as pd

from .config import EditorProfile, get_editor, load_config
from .crispr import reverse_complement

HIT_COLUMNS = [
    "candidate_id",
    "query",
    "off_target_sequence",
    "chromosome",
    "location_0based",
    "direction",
    "mismatches",
    "bulge_type",
    "bulge_size",
]


def _hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError("Hamming distance requires equal-length strings")
    return sum(x != y for x, y in zip(a, b, strict=True))


def _enumerate_spcas9_sites(records: Mapping[str, str]) -> Iterable[dict[str, object]]:
    for seqid, raw_sequence in records.items():
        sequence = raw_sequence.upper().replace("-", "")
        for index in range(0, len(sequence) - 22):
            protospacer = sequence[index : index + 20]
            pam = sequence[index + 20 : index + 23]
            if set(protospacer + pam) <= set("ACGT") and pam[1:] == "GG":
                yield {
                    "seqid": seqid,
                    "start_1based": index + 1,
                    "end_1based": index + 20,
                    "strand": "+",
                    "guide": protospacer,
                    "pam": pam,
                }
            pam_on_plus = sequence[index : index + 3]
            genomic_target = sequence[index + 3 : index + 23]
            if set(pam_on_plus + genomic_target) <= set("ACGT") and pam_on_plus[:2] == "CC":
                yield {
                    "seqid": seqid,
                    "start_1based": index + 4,
                    "end_1based": index + 23,
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
    """Exhaustively screen a teaching-scale host FASTA, never a full human genome."""
    total_bases = sum(len(seq.replace("-", "")) for seq in host_records.values())
    if total_bases > max_total_host_bases:
        raise ValueError(
            f"Host FASTA has {total_bases:,} bases; use Cas-OFFinder for genome-scale screening."
        )
    host_sites = list(_enumerate_spcas9_sites(host_records))
    rows: list[dict[str, object]] = []
    for _, candidate in candidates.iterrows():
        hits = []
        for site in host_sites:
            mismatches = _hamming(str(candidate["guide_sequence"]), str(site["guide"]))
            if mismatches <= max_mismatches:
                hits.append({**site, "mismatches": mismatches})
        hits.sort(key=lambda item: (item["mismatches"], item["seqid"], item["start_1based"]))
        best = hits[0] if hits else None
        row: dict[str, object] = {
            "candidate_id": candidate["candidate_id"],
            "human_exact_hit_count": sum(hit["mismatches"] == 0 for hit in hits),
            "human_one_mismatch_hit_count": sum(hit["mismatches"] == 1 for hit in hits),
            "human_two_mismatch_hit_count": sum(hit["mismatches"] == 2 for hit in hits),
            "human_three_mismatch_hit_count": sum(hit["mismatches"] == 3 for hit in hits),
            "human_total_predicted_hits": len(hits),
            "human_minimum_mismatch_count": best["mismatches"] if best else pd.NA,
            "human_highest_risk_location": (
                f"{best['seqid']}:{best['start_1based']}({best['strand']})" if best else ""
            ),
            "host_exact_matches": sum(hit["mismatches"] == 0 for hit in hits),
            f"host_matches_le_{max_mismatches}_mismatches": len(hits),
            "host_min_mismatches": best["mismatches"] if best else pd.NA,
        }
        rows.append(row)
    return candidates.merge(pd.DataFrame(rows), on="candidate_id", how="left")


def _configured(config: dict[str, Any] | str | Path | None) -> dict[str, Any]:
    return config if isinstance(config, dict) else load_config(config)


def select_offtarget_candidates(
    candidates: pd.DataFrame,
    *,
    maximum_candidates: int,
    genes: Sequence[str] | None = None,
    stratify_by_gene: bool = True,
) -> pd.DataFrame:
    """Select ranked, unique guides deterministically with optional gene strata."""
    working = candidates.copy()
    if "rejection_reasons" in working:
        working = working[working["rejection_reasons"].fillna("").eq("")]
    if genes:
        working = working[working.get("gene_name", "").fillna("").isin(genes)]
    if "gene_name" not in working:
        working["gene_name"] = ""
    score = "post_human_score" if "post_human_score" in working else "pre_human_score"
    if score not in working:
        score = "virus_site_coverage"
    working = working.sort_values(
        [score, "candidate_id"], ascending=[False, True], kind="mergesort"
    ).drop_duplicates("guide_sequence", keep="first")
    working["selection_stratum"] = (
        working["gene_name"].fillna("").replace("", "intergenic_or_unannotated")
    )
    working["within_stratum_rank"] = working.groupby("selection_stratum").cumcount() + 1
    if stratify_by_gene:
        working = working.sort_values(
            ["within_stratum_rank", score, "selection_stratum", "candidate_id"],
            ascending=[True, False, True, True],
            kind="mergesort",
        )
    return working.iloc[:maximum_candidates].reset_index(drop=True)


def _cas_pattern(editor: EditorProfile) -> str:
    if editor.pam_orientation == "3prime":
        return "N" * editor.protospacer_length + editor.pam_pattern
    return editor.pam_pattern + "N" * editor.protospacer_length


def _cas_query(guide: str, editor: EditorProfile) -> str:
    pam_wildcards = "N" * len(editor.pam_pattern)
    return guide + pam_wildcards if editor.pam_orientation == "3prime" else pam_wildcards + guide


def build_cas_offinder_input(
    candidates: pd.DataFrame,
    human_fasta_directory: str | Path,
    output_path: str | Path,
    manifest_path: str | Path,
    *,
    maximum_candidates: int,
    genes: Sequence[str] | None = None,
    stratify_by_gene: bool = True,
    config: dict[str, Any] | str | Path | None = None,
) -> pd.DataFrame:
    """Write valid Cas-OFFinder input plus the required query-to-candidate manifest."""
    settings = _configured(config)
    editor = get_editor(settings)
    selected = select_offtarget_candidates(
        candidates,
        maximum_candidates=maximum_candidates,
        genes=genes,
        stratify_by_gene=stratify_by_gene,
    )
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(Path(human_fasta_directory).resolve()), _cas_pattern(editor)]
    selected = selected.copy()
    selected["cas_offinder_query"] = selected["guide_sequence"].map(
        lambda guide: _cas_query(str(guide), editor)
    )
    selected["human_assembly"] = settings["off_target"]["human_assembly"]
    selected["human_assembly_accession"] = settings["off_target"]["human_assembly_accession"]
    selected["editor_profile"] = editor.name
    selected["mismatch_search_threshold"] = editor.mismatch_search_threshold
    lines.extend(
        f"{query} {editor.mismatch_search_threshold}" for query in selected["cas_offinder_query"]
    )
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    manifest = Path(manifest_path)
    manifest.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(manifest, index=False)
    return selected


def write_cas_offinder_input(
    candidates: pd.DataFrame,
    human_fasta_directory: str | Path,
    output_path: str | Path,
    max_mismatches: int = 3,
) -> Path:
    """Backward-compatible standards-compliant Cas-OFFinder input writer."""
    settings = load_config()
    settings["editor"]["mismatch_search_threshold"] = max_mismatches
    output = Path(output_path)
    build_cas_offinder_input(
        candidates,
        human_fasta_directory,
        output,
        output.with_suffix(output.suffix + ".manifest.csv"),
        maximum_candidates=len(candidates),
        stratify_by_gene=False,
        config=settings,
    )
    return output


def read_cas_offinder_output(path: str | Path) -> pd.DataFrame:
    """Parse standard six-column Cas-OFFinder or common bulge-enabled output."""
    rows: list[dict[str, object]] = []
    with Path(path).open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) == 6:
                query, chromosome, location, observed, direction, mismatches = fields
                bulge_type, bulge_size = "none", 0
            elif len(fields) >= 8 and fields[1] in {"DNA", "RNA", "X", "none"}:
                (
                    query,
                    bulge_type,
                    observed,
                    chromosome,
                    location,
                    direction,
                    mismatches,
                    bulge_size,
                ) = fields[:8]
            elif len(fields) >= 8 and fields[0] in {"DNA", "RNA", "X", "none"}:
                (
                    bulge_type,
                    query,
                    observed,
                    chromosome,
                    location,
                    direction,
                    mismatches,
                    bulge_size,
                ) = fields[:8]
            elif len(fields) >= 7:
                # Legacy repository fixture with an explicit candidate ID prefix.
                candidate_id, query, observed, chromosome, location, direction, mismatches = fields[
                    :7
                ]
                bulge_type, bulge_size = "none", 0
            else:
                raise ValueError(f"Unsupported Cas-OFFinder row at line {line_number}: {line!r}")
            rows.append(
                {
                    "candidate_id": locals().get("candidate_id", ""),
                    "query": query.upper(),
                    "off_target_sequence": observed.upper(),
                    "chromosome": chromosome,
                    "location_0based": int(location),
                    "direction": direction,
                    "mismatches": int(mismatches),
                    "bulge_type": bulge_type,
                    "bulge_size": int(bulge_size),
                }
            )
            candidate_id = ""
    return pd.DataFrame(rows, columns=HIT_COLUMNS)


def annotate_human_hits(hits: pd.DataFrame, human_gff: str | Path | None) -> pd.DataFrame:
    """Annotate hit coordinates with a single streaming pass over a potentially large GFF."""
    output = hits.copy()
    output["human_annotation"] = ""
    if human_gff is None or output.empty:
        return output

    positions: dict[str, list[tuple[int, object]]] = {}
    for index, hit in output.iterrows():
        coordinate = int(hit["location_0based"]) + 1
        hit_seqid = str(hit["chromosome"]).split()[0]
        positions.setdefault(hit_seqid, []).append((coordinate, index))
    for values in positions.values():
        values.sort()
    coordinate_values = {
        seqid: [coordinate for coordinate, _ in values] for seqid, values in positions.items()
    }

    annotations: dict[object, set[str]] = {index: set() for index in output.index}
    path = Path(human_gff)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[0] not in positions:
                continue
            start, end = int(fields[3]), int(fields[4])
            coordinate_rows = positions[fields[0]]
            coordinates = coordinate_values[fields[0]]
            left = bisect.bisect_left(coordinates, start)
            right = bisect.bisect_right(coordinates, end)
            if left == right:
                continue
            attributes: dict[str, str] = {}
            for item in fields[8].strip().strip(";").split(";"):
                if "=" in item:
                    key, value = item.split("=", 1)
                    attributes[key] = unquote(value)
            label = (
                attributes.get("Name")
                or attributes.get("gene")
                or attributes.get("ID")
                or fields[2]
            )
            for _, row_index in coordinate_rows[left:right]:
                annotations[row_index].add(str(label))
    output["human_annotation"] = [
        ";".join(sorted(annotations[index])) for index in output.index
    ]
    return output


def summarize_cas_offinder_hits(
    candidates: pd.DataFrame,
    hits: pd.DataFrame,
    *,
    max_mismatches: int = 3,
    selected_manifest: pd.DataFrame | str | Path | None = None,
    human_gff: str | Path | None = None,
    config: dict[str, Any] | str | Path | None = None,
) -> pd.DataFrame:
    """Return one transparent per-candidate summary from predicted human hits."""
    settings = _configured(config)
    editor = get_editor(settings)
    if isinstance(selected_manifest, (str, Path)):
        selected_manifest = pd.read_csv(selected_manifest)
    mapped_hits = hits.copy()
    for column in HIT_COLUMNS:
        if column not in mapped_hits:
            mapped_hits[column] = pd.Series(dtype="object")
    if not mapped_hits.empty and selected_manifest is not None:
        query_map = selected_manifest.set_index("cas_offinder_query")["candidate_id"].to_dict()
        mapped_hits["candidate_id"] = (
            mapped_hits["query"].map(query_map).fillna(mapped_hits.get("candidate_id", ""))
        )
    if not mapped_hits.empty and mapped_hits["candidate_id"].fillna("").eq("").any():
        guide_map = (
            candidates.drop_duplicates("guide_sequence")
            .set_index("guide_sequence")["candidate_id"]
            .to_dict()
        )
        mapped_hits.loc[mapped_hits["candidate_id"].eq(""), "candidate_id"] = (
            mapped_hits.loc[mapped_hits["candidate_id"].eq(""), "query"]
            .str[: editor.protospacer_length]
            .map(guide_map)
        )
    mapped_hits = mapped_hits[mapped_hits["mismatches"] <= max_mismatches].copy()
    mapped_hits = annotate_human_hits(mapped_hits, human_gff)
    if not mapped_hits.empty:
        pam_length = len(editor.pam_pattern)
        mapped_hits["observed_human_sequence"] = mapped_hits["off_target_sequence"]
        mapped_hits["pam_compatibility"] = (
            mapped_hits["off_target_sequence"]
            .str[-pam_length:]
            .map(lambda pam: pam.endswith("GG") if editor.name == "SpCas9" else pd.NA)
        )
        mapped_hits["human_coordinate_1based"] = mapped_hits["location_0based"] + 1
    rows: list[dict[str, object]] = []
    for candidate_id in candidates["candidate_id"].astype(str):
        group = mapped_hits[mapped_hits["candidate_id"].astype(str) == candidate_id].sort_values(
            ["mismatches", "chromosome", "location_0based"], kind="mergesort"
        )
        best = group.iloc[0] if not group.empty else None
        rows.append(
            {
                "candidate_id": candidate_id,
                "human_exact_hit_count": int((group["mismatches"] == 0).sum()),
                "human_one_mismatch_hit_count": int((group["mismatches"] == 1).sum()),
                "human_two_mismatch_hit_count": int((group["mismatches"] == 2).sum()),
                "human_three_mismatch_hit_count": int((group["mismatches"] == 3).sum()),
                "human_minimum_mismatch_count": (
                    int(best["mismatches"]) if best is not None else pd.NA
                ),
                "human_total_predicted_hits": int(len(group)),
                "highest_risk_predicted_hit": (
                    str(best["off_target_sequence"]) if best is not None else ""
                ),
                "human_hit_chromosome_or_contig": (
                    str(best["chromosome"]) if best is not None else ""
                ),
                "human_hit_coordinate_1based": (
                    int(best["location_0based"]) + 1 if best is not None else pd.NA
                ),
                "human_hit_strand": str(best["direction"]) if best is not None else "",
                "observed_human_sequence": (
                    str(best["off_target_sequence"]) if best is not None else ""
                ),
                "pam_compatibility": best["pam_compatibility"] if best is not None else pd.NA,
                "human_annotation": str(best["human_annotation"]) if best is not None else "",
            }
        )
    output = candidates.merge(pd.DataFrame(rows), on="candidate_id", how="left")
    output.attrs["predicted_human_hits"] = mapped_hits
    return output


def write_offtarget_metadata(path: str | Path, settings: dict[str, Any]) -> Path:
    output = Path(path)
    output.write_text(
        json.dumps(
            {
                "human_assembly": settings["off_target"]["human_assembly"],
                "human_assembly_accession": settings["off_target"]["human_assembly_accession"],
                "editor": settings["editor"],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return output
