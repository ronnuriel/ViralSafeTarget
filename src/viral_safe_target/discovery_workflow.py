"""Resumable HSV-2 genome-wide discovery orchestration."""

from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from .annotations import read_gff3
from .config import get_editor, load_config
from .discovery import (
    build_bounded_pair_hypotheses,
    build_candidate_feature_map,
    build_deep_screening_panel,
    gene_rank_stability,
    rank_genes,
    select_balanced_discovery_panel,
)
from .offtarget import HIT_COLUMNS, read_cas_offinder_output, summarize_cas_offinder_hits
from .provenance import sha256_file
from .scoring import rank_post_human_candidates

HUMAN_COUNT_COLUMNS = [
    "human_exact_hit_count",
    "human_one_mismatch_hit_count",
    "human_two_mismatch_hit_count",
    "human_three_mismatch_hit_count",
    "human_total_predicted_hits",
]


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        values = yaml.safe_load(handle) or {}
    if str(values.get("schema_version")) != "0.5":
        raise ValueError("Genome-wide configuration schema_version must be '0.5'.")
    return values


def _resolve(root: Path, value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else root / path


def _json(path: Path, values: dict[str, Any]) -> None:
    path.write_text(json.dumps(values, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _cas_executable(root: Path) -> Path | None:
    choices = [
        os.environ.get("CAS_OFFINDER_PATH", ""),
        str(root / "tools/bin/cas-offinder"),
        shutil.which("cas-offinder") or "",
    ]
    for value in choices:
        path = Path(value).expanduser() if value else None
        if path and path.is_file() and os.access(path, os.X_OK):
            return path.resolve()
    return None


def _tool_version(executable: Path | None) -> str:
    if executable is None:
        return "unavailable"
    try:
        result = subprocess.run(
            [str(executable), "--help"], capture_output=True, text=True, timeout=10, check=False
        )
        line = (result.stdout or result.stderr).strip().splitlines()
        return line[0] if line else "available; version not reported"
    except (OSError, subprocess.TimeoutExpired) as error:
        return f"version check unavailable: {error}"


def _opencl_information() -> str:
    executable = shutil.which("clinfo")
    if not executable:
        return "clinfo unavailable; Cas-OFFinder device selector requested: C"
    result = subprocess.run(
        [executable, "--raw"], capture_output=True, text=True, timeout=15, check=False
    )
    lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    return " | ".join(lines[:8]) or "clinfo returned no device description"


def _find_human_directory(root: Path) -> tuple[Path, Path]:
    fasta_files = sorted(root.glob("**/*_genomic.fna"))
    if not fasta_files:
        raise FileNotFoundError(f"No cached GRCh38 genomic FASTA was found below {root}.")
    return fasta_files[0].parent.resolve(), fasta_files[0].resolve()


def _cas_query(guide: str, editor: Any) -> str:
    wildcards = "N" * len(editor.pam_pattern)
    return guide + wildcards if editor.pam_orientation == "3prime" else wildcards + guide


def _cas_pattern(editor: Any) -> str:
    guide = "N" * editor.protospacer_length
    return (
        guide + editor.pam_pattern
        if editor.pam_orientation == "3prime"
        else editor.pam_pattern + guide
    )


def prepare_batches(
    panel: pd.DataFrame,
    batch_root: Path,
    human_directory: Path,
    settings: dict[str, Any],
    *,
    batch_size: int,
) -> list[dict[str, Any]]:
    """Create deterministic unique-query batches while retaining all candidate mappings."""
    editor = get_editor(settings)
    batch_root.mkdir(parents=True, exist_ok=True)
    query_table = panel[["candidate_id", "guide_sequence"]].copy()
    query_table["cas_offinder_query"] = query_table["guide_sequence"].map(
        lambda value: _cas_query(str(value), editor)
    )
    unique = query_table.drop_duplicates("cas_offinder_query").sort_values(
        ["cas_offinder_query", "candidate_id"], kind="mergesort"
    )
    batches: list[dict[str, Any]] = []
    for start in range(0, len(unique), batch_size):
        number = start // batch_size + 1
        batch_dir = batch_root / f"batch_{number:04d}"
        batch_dir.mkdir(parents=True, exist_ok=True)
        queries = unique.iloc[start : start + batch_size]["cas_offinder_query"]
        manifest = query_table[query_table["cas_offinder_query"].isin(queries)].sort_values(
            ["cas_offinder_query", "candidate_id"], kind="mergesort"
        )
        manifest_path = batch_dir / "candidate_manifest.csv"
        input_path = batch_dir / "input.txt"
        raw_path = batch_dir / "raw_output.tsv"
        status_path = batch_dir / "manifest.json"
        previous = json.loads(status_path.read_text()) if status_path.is_file() else {}
        manifest.to_csv(manifest_path, index=False)
        lines = [str(human_directory), _cas_pattern(editor)]
        lines.extend(f"{query} {editor.mismatch_search_threshold}" for query in queries.astype(str))
        input_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        input_hash = sha256_file(input_path)
        candidate_hash = sha256_file(manifest_path)
        same_inputs = (
            previous.get("input_sha256") == input_hash
            and previous.get("candidate_manifest_sha256") == candidate_hash
        )
        reusable = previous.get("status") == "completed" and same_inputs and raw_path.is_file()
        status = (
            "completed"
            if reusable
            else "failed"
            if previous.get("status") == "failed" and same_inputs
            else "pending"
        )
        batch = {
            "batch_number": number,
            "batch_dir": batch_dir,
            "input_path": input_path,
            "manifest_path": manifest_path,
            "raw_path": raw_path,
            "status_path": status_path,
            "input_sha256": input_hash,
            "candidate_manifest_sha256": candidate_hash,
            "unique_query_count": len(queries),
            "candidate_count": len(manifest),
            "status": status,
            "previous": previous,
        }
        batches.append(batch)
        if not reusable:
            _json(
                status_path,
                {
                    **previous,
                    "schema_version": "0.5",
                    "batch_number": number,
                    "status": status,
                    "candidate_count": len(manifest),
                    "unique_query_count": len(queries),
                    "input_sha256": input_hash,
                    "candidate_manifest_sha256": candidate_hash,
                    "raw_output": str(raw_path.resolve()),
                },
            )
    return batches


def _write_batch_status(batch: dict[str, Any], **updates: Any) -> dict[str, Any]:
    record = {
        "schema_version": "0.5",
        "batch_number": batch["batch_number"],
        "status": batch["status"],
        "candidate_count": batch["candidate_count"],
        "unique_query_count": batch["unique_query_count"],
        "input_sha256": batch["input_sha256"],
        "candidate_manifest_sha256": batch["candidate_manifest_sha256"],
        "raw_output": str(batch["raw_path"].resolve()),
        **updates,
    }
    _json(batch["status_path"], record)
    batch["status"] = str(record["status"])
    return record


def run_batches(
    batches: list[dict[str, Any]],
    executable: Path | None,
    *,
    execute: bool,
) -> None:
    completed_at_start = sum(batch["status"] == "completed" for batch in batches)
    pending = len(batches) - completed_at_start
    total_queries = sum(batch["unique_query_count"] for batch in batches)
    print(f"Discovery panel unique guides: {total_queries:,}")
    print(f"Batches: {len(batches):,}; completed: {completed_at_start:,}; pending: {pending:,}")
    print(f"Relative workload versus a 200-guide pilot: {total_queries / 200:.2f}x")
    print("Previous pilot duration: unavailable in a machine-readable timing artifact")
    if not execute:
        print("Analysis-only mode: pending and failed batches will remain explicitly incomplete.")
        return
    if executable is None:
        print("Cas-OFFinder is unavailable; leaving pending batches incomplete.")
        return
    started = time.monotonic()
    total = len(batches)
    for index, batch in enumerate(batches, start=1):
        if batch["status"] == "completed":
            continue
        command = [str(executable), str(batch["input_path"]), "C", str(batch["raw_path"])]
        batch_started = time.monotonic()
        started_utc = _utc_now()
        print(f"Batch {index}/{total} ({100 * (index - 1) / max(total, 1):.1f}% complete)")
        result = subprocess.run(command, capture_output=True, text=True, check=False)
        elapsed = time.monotonic() - batch_started
        status = "completed" if result.returncode == 0 and batch["raw_path"].is_file() else "failed"
        _write_batch_status(
            batch,
            status=status,
            command=command,
            returncode=result.returncode,
            stdout=result.stdout[-4000:],
            stderr=result.stderr[-4000:],
            started_utc=started_utc,
            completed_utc=_utc_now(),
            elapsed_seconds=elapsed,
            executable_sha256=sha256_file(executable),
            tool_version=_tool_version(executable),
            device_information=_opencl_information(),
            raw_output_sha256=(
                sha256_file(batch["raw_path"]) if batch["raw_path"].is_file() else ""
            ),
        )
        print(f"  {status} in {elapsed:.1f}s; elapsed total {time.monotonic() - started:.1f}s")


def merge_batch_results(
    panel: pd.DataFrame,
    batches: list[dict[str, Any]],
    settings: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Merge only validated completed batches; pending is never interpreted as zero hits."""
    hit_frames: list[pd.DataFrame] = []
    completed_ids: set[str] = set()
    status_by_id: dict[str, str] = {}
    for batch in batches:
        manifest = pd.read_csv(batch["manifest_path"])
        candidate_ids = manifest["candidate_id"].astype(str)
        for candidate_id in candidate_ids:
            status_by_id[candidate_id] = batch["status"]
        if batch["status"] != "completed":
            continue
        completed_ids.update(candidate_ids)
        parsed = read_cas_offinder_output(batch["raw_path"])
        if parsed.empty:
            continue
        expanded = parsed.drop(columns=["candidate_id"], errors="ignore").merge(
            manifest[["candidate_id", "guide_sequence", "cas_offinder_query"]],
            left_on="query",
            right_on="cas_offinder_query",
            how="left",
            validate="many_to_many",
        )
        expanded["batch_number"] = batch["batch_number"]
        expanded["raw_output_source"] = str(batch["raw_path"].resolve())
        expanded["raw_output_sha256"] = sha256_file(batch["raw_path"])
        expanded["tool_name"] = "cas-offinder"
        expanded["human_assembly"] = settings["off_target"]["human_assembly"]
        hit_frames.append(expanded)
    hit_columns = [
        *HIT_COLUMNS,
        "guide_sequence",
        "cas_offinder_query",
        "batch_number",
        "raw_output_source",
        "raw_output_sha256",
        "tool_name",
        "human_assembly",
    ]
    hits = (
        pd.concat(hit_frames, ignore_index=True)
        if hit_frames
        else pd.DataFrame(columns=hit_columns)
    )
    if not hits.empty:
        hits = hits.sort_values(
            ["candidate_id", "mismatches", "chromosome", "location_0based", "direction"],
            kind="mergesort",
        ).reset_index(drop=True)
    completed = panel[panel["candidate_id"].astype(str).isin(completed_ids)].copy()
    if completed.empty:
        ranked_completed = completed
    else:
        summarized = summarize_cas_offinder_hits(
            completed,
            hits,
            max_mismatches=int(settings["editor"]["mismatch_search_threshold"]),
            config=settings,
        )
        ranked_completed = rank_post_human_candidates(summarized, settings)
        ranked_completed["post_human_rank"] = range(1, len(ranked_completed) + 1)
    post = panel.copy()
    result_columns = [
        *HUMAN_COUNT_COLUMNS,
        "human_minimum_mismatch_count",
        "highest_risk_predicted_hit",
        "human_hit_chromosome_or_contig",
        "human_hit_coordinate_1based",
        "human_hit_strand",
        "observed_human_sequence",
        "pam_compatibility",
        "human_annotation",
        "predicted_offtarget_risk",
        "post_human_score",
        "post_human_rank",
        "decision",
        "decision_reason",
    ]
    available = [column for column in result_columns if column in ranked_completed]
    if available:
        replacement = ranked_completed[["candidate_id", *available]]
        post = post.drop(columns=available, errors="ignore").merge(
            replacement, on="candidate_id", how="left", validate="one_to_one"
        )
    for column in result_columns:
        if column not in post:
            post[column] = pd.NA
    post["screening_status"] = (
        post["candidate_id"].astype(str).map(status_by_id).fillna("incomplete")
    )
    incomplete = post["screening_status"].ne("completed")
    post.loc[
        incomplete,
        HUMAN_COUNT_COLUMNS + ["predicted_offtarget_risk", "post_human_score", "post_human_rank"],
    ] = pd.NA
    post.loc[incomplete, "decision"] = "screening_incomplete"
    post.loc[incomplete, "decision_reason"] = (
        "Cas-OFFinder batch is pending, failed, or missing; absence of a result "
        "is not a zero-hit result."
    )
    post["tool_coverage"] = post["screening_status"].eq("completed").astype(float)
    post["human_screen_explanation"] = post.apply(
        lambda row: (
            f"Cas-OFFinder completed against {settings['off_target']['human_assembly']} through "
            f"{settings['editor']['mismatch_search_threshold']} mismatches; zero hits "
            "is a model-bounded prediction, not proof of safety."
            if row["screening_status"] == "completed"
            else str(row["decision_reason"])
        ),
        axis=1,
    )
    post = post.sort_values(
        ["post_human_rank", "pre_human_rank", "candidate_id"], na_position="last", kind="mergesort"
    ).reset_index(drop=True)
    return hits, post


def _aggregate_top_per_gene(
    candidates: pd.DataFrame, feature_map: pd.DataFrame, *, count: int = 10
) -> pd.DataFrame:
    rows = []
    ranked = candidates[candidates["post_human_rank"].notna()]
    for gene, mappings in feature_map[feature_map["gene_name"].ne("")].groupby("gene_name"):
        ids = set(mappings["candidate_id"])
        subset = ranked[ranked["candidate_id"].isin(ids)].head(count).copy()
        subset.insert(0, "mapped_gene_for_view", gene)
        rows.append(subset)
    return pd.concat(rows, ignore_index=True) if rows else ranked.head(0).copy()


def _external_templates(deep: pd.DataFrame, output: Path) -> None:
    columns = [
        column
        for column in ["candidate_id", "guide_sequence", "pam", "mapped_gene_names"]
        if column in deep
    ]
    for tool in ("crispritz", "crispor", "chopchop", "guidescan2"):
        directory = output / "external_tools" / tool
        directory.mkdir(parents=True, exist_ok=True)
        deep[columns].to_csv(directory / "candidate_input.csv", index=False)
        _json(
            directory / "status.json",
            {
                "tool": tool,
                "status": "pending_external_execution_or_import",
                "candidate_count": len(deep),
                "note": "Optional external results are not required for the genome-wide report.",
            },
        )


def run_genome_wide_discovery(
    *,
    virus: str,
    run_dir: str | Path | None,
    config_path: str | Path,
    out_dir: str | Path,
    top_per_gene: int | None = None,
    global_top: int | None = None,
    batch_size: int | None = None,
    run_cas_offinder: bool = False,
    analysis_only: bool = False,
    exhaustive: bool = False,
    confirm_exhaustive: bool = False,
) -> dict[str, Any]:
    if virus.lower() != "hsv2":
        raise ValueError(
            "v0.5 currently provides a validated genome-wide profile for virus='hsv2'."
        )
    if analysis_only and run_cas_offinder:
        raise ValueError("--analysis-only and --run-cas-offinder are mutually exclusive.")
    root = Path.cwd().resolve()
    discovery_config_path = _resolve(root, config_path).resolve()
    config = _read_yaml(discovery_config_path)
    source = _resolve(root, run_dir or config["source_run"]).resolve()
    output = _resolve(root, out_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    editor_config_path = _resolve(root, config["editor_config"]).resolve()
    settings = load_config(editor_config_path)
    candidate_path = source / "candidates_ranked_pre_human.csv"
    annotation_path = _resolve(root, config["annotation"]).resolve()
    qc_path = source / "accession_qc.csv"
    if not candidate_path.is_file() or not annotation_path.is_file():
        raise FileNotFoundError(
            "The cached pre-human candidate table or HSV-2 annotation is missing."
        )
    timings: dict[str, Any] = {"started_utc": _utc_now(), "stages": {}}

    stage = time.monotonic()
    candidates = pd.read_csv(candidate_path)
    eligible = candidates[candidates["rejection_reasons"].fillna("").eq("")].copy()
    eligible = eligible.sort_values(
        ["pre_human_score", "candidate_id"],
        ascending=[False, True],
        kind="mergesort",
    ).reset_index(drop=True)
    eligible["pre_human_rank"] = range(1, len(eligible) + 1)
    eligible["candidate_provenance"] = str(candidate_path)
    eligible["candidate_source_sha256"] = sha256_file(candidate_path)
    features = read_gff3(annotation_path)
    feature_map = build_candidate_feature_map(
        eligible, features, config=settings, annotation_source=annotation_path
    )
    feature_map.to_csv(output / "candidate_feature_map.csv", index=False)
    mapping_summary = (
        feature_map.groupby("candidate_id", sort=True)
        .agg(
            mapped_feature_ids=(
                "feature_id",
                lambda values: ";".join(sorted(set(filter(None, values.astype(str))))),
            ),
            mapped_gene_names=(
                "gene_name",
                lambda values: ";".join(sorted(set(filter(None, values.astype(str))))),
            ),
            mapped_feature_types=(
                "feature_type",
                lambda values: ";".join(sorted(set(values.astype(str)))),
            ),
            feature_mapping_count=("feature_id", "size"),
        )
        .reset_index()
    )
    eligible_export = eligible.merge(
        mapping_summary, on="candidate_id", how="left", validate="one_to_one"
    )
    eligible_export["mapping_provenance"] = str((output / "candidate_feature_map.csv").resolve())
    eligible_export["screening_status"] = "not_selected_or_pending"
    eligible_export["tool_coverage"] = 0.0
    eligible_export["post_human_score"] = pd.NA
    eligible_export["post_human_rank"] = pd.NA
    for column in HUMAN_COUNT_COLUMNS:
        eligible_export[column] = pd.NA
    eligible_export.to_csv(output / "genome_wide_candidates_pre_human.csv", index=False)
    timings["stages"]["mapping_seconds"] = time.monotonic() - stage

    selected = select_balanced_discovery_panel(
        eligible,
        feature_map,
        features,
        top_per_gene=top_per_gene or int(config["top_per_gene"]),
        global_top=global_top or int(config["global_top"]),
        exhaustive=exhaustive,
        confirm_exhaustive=confirm_exhaustive,
    )
    panel = selected.panel.copy()
    panel["mapping_provenance"] = str((output / "candidate_feature_map.csv").resolve())
    panel_export = panel.copy()
    panel_export["screening_status"] = "pending"
    panel_export["tool_coverage"] = 0.0
    panel_export["post_human_score"] = pd.NA
    panel_export["post_human_rank"] = pd.NA
    for column in HUMAN_COUNT_COLUMNS:
        panel_export[column] = pd.NA
    panel_export.to_csv(output / "genome_wide_screening_panel.csv", index=False)
    selected.audit.to_csv(output / "selection_audit.csv", index=False)
    selected.genes_without_candidates.to_csv(
        output / "genes_without_eligible_candidates.csv", index=False
    )
    print(f"Eligible pre-human candidates: {len(eligible):,}")
    print(f"Balanced discovery panel candidates: {len(panel):,}")
    print(f"Panel unique guides: {panel['guide_sequence'].nunique():,}")

    human_directory, human_fasta = _find_human_directory(_resolve(root, config["human_fasta_root"]))
    size_ratio = panel["guide_sequence"].nunique() / 200
    prior_raw = source / "cas_offinder_output.tsv"
    estimated_bytes = int(prior_raw.stat().st_size * size_ratio) if prior_raw.is_file() else None
    print(
        f"Estimated raw-output disk requirement: {estimated_bytes / 1024**2:.1f} MiB based on pilot"
        if estimated_bytes is not None
        else "Estimated raw-output disk requirement: unavailable (no cached pilot raw output)"
    )
    batches = prepare_batches(
        panel,
        output / "batches",
        human_directory,
        settings,
        batch_size=batch_size or int(config["batch_size"]),
    )
    stage = time.monotonic()
    executable = _cas_executable(root)
    run_batches(batches, executable, execute=run_cas_offinder and not analysis_only)
    timings["stages"]["cas_offinder_seconds"] = time.monotonic() - stage
    _json(
        output / "combined_batch_manifest.json",
        {
            "schema_version": "0.5",
            "candidate_count": len(panel),
            "unique_guide_count": int(panel["guide_sequence"].nunique()),
            "batches": [
                {
                    "batch_number": batch["batch_number"],
                    "status": batch["status"],
                    "manifest": str(batch["manifest_path"].resolve()),
                    "raw_output": str(batch["raw_path"].resolve()),
                }
                for batch in batches
            ],
        },
    )

    stage = time.monotonic()
    hits, post = merge_batch_results(panel, batches, settings)
    hits.to_csv(output / "genome_wide_human_hits.csv", index=False)
    post.to_csv(output / "genome_wide_candidates_post_human.csv", index=False)
    top_global = post[post["post_human_rank"].notna()].head(100)
    top_global.to_csv(output / "top_candidates_global.csv", index=False)
    top_per_gene_frame = _aggregate_top_per_gene(post, feature_map)
    top_per_gene_frame.to_csv(output / "top_candidates_per_gene.csv", index=False)

    evidence_path = _resolve(root, config.get("gene_evidence", "data/curated/gene_evidence.tsv"))
    evidence = pd.read_csv(evidence_path, sep="\t") if evidence_path.is_file() else None
    gene_rankings = rank_genes(post, feature_map, features, evidence=evidence)
    gene_rankings["ranking_provenance"] = str(
        (output / "genome_wide_candidates_post_human.csv").resolve()
    )
    gene_rankings.to_csv(output / "gene_rankings.csv", index=False)
    stability = gene_rank_stability(post, feature_map, features)
    stability.to_csv(output / "gene_rank_stability.csv", index=False)
    deep_config = config["deep_panel"]
    deep = build_deep_screening_panel(
        post,
        feature_map,
        gene_rankings,
        top_genes=int(deep_config["top_genes"]),
        top_per_gene=int(deep_config["top_per_gene"]),
        global_top=int(deep_config["global_top"]),
    )
    deep.to_csv(output / "deep_screening_panel.csv", index=False)
    _external_templates(deep, output)
    pair_config = config["pairs"]
    same, multi, pair_summary = build_bounded_pair_hypotheses(
        post[post["post_human_rank"].notna()],
        feature_map,
        gene_rankings,
        top_genes=int(pair_config["top_genes"]),
        candidates_per_gene=int(pair_config["candidates_per_gene"]),
        maximum_pairs=int(pair_config["maximum_pairs"]),
    )
    same.to_csv(output / "pair_hypotheses_same_gene.csv", index=False)
    multi.to_csv(output / "pair_hypotheses_multi_target.csv", index=False)
    pair_summary.to_csv(output / "pair_summary_by_gene.csv", index=False)
    timings["stages"]["analysis_seconds"] = time.monotonic() - stage

    qc = pd.read_csv(qc_path) if qc_path.is_file() else pd.DataFrame()
    provenance = {
        "schema_version": "0.5",
        "created_utc": _utc_now(),
        "command": sys.argv,
        "platform": platform.platform(),
        "discovery_config": str(discovery_config_path),
        "discovery_config_sha256": sha256_file(discovery_config_path),
        "editor_config": str(editor_config_path),
        "editor_config_sha256": sha256_file(editor_config_path),
        "candidate_source": str(candidate_path),
        "candidate_source_sha256": sha256_file(candidate_path),
        "annotation_source": str(annotation_path),
        "annotation_source_sha256": sha256_file(annotation_path),
        "human_fasta": str(human_fasta),
        "human_fasta_size_bytes": human_fasta.stat().st_size,
        "human_assembly": config["human_assembly"],
        "human_assembly_accession": config["human_assembly_accession"],
        "cas_offinder_executable": str(executable) if executable else "unavailable",
        "cas_offinder_version": _tool_version(executable),
        "opencl_device_information": _opencl_information(),
        "analysis_only": analysis_only,
        "completed_batches": sum(batch["status"] == "completed" for batch in batches),
        "total_batches": len(batches),
    }
    _json(output / "provenance.json", provenance)
    timings["completed_utc"] = _utc_now()
    timings["total_seconds"] = sum(float(value) for value in timings["stages"].values())
    _json(output / "stage_timings.json", timings)

    from .discovery_reporting import write_discovery_report

    answers = write_discovery_report(
        output / "report.html",
        candidates=post,
        feature_map=feature_map,
        genes=gene_rankings,
        stability=stability,
        deep_panel=deep,
        same_pairs=same,
        multi_pairs=multi,
        genes_without_candidates=selected.genes_without_candidates,
        qc=qc,
        provenance=provenance,
        initial_candidate_count=len(candidates),
        eligible_candidate_count=len(eligible),
    )
    print("Output paths:")
    for path in sorted(output.iterdir()):
        if path.is_file():
            print(f"- {path}")
    return {"output_dir": output, "answers": answers, "provenance": provenance}
