"""Cached HSV-2 32-candidate multi-tool consensus workflow."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from .cache import stage_is_current, write_stage_stamp
from .consensus import candidate_metrics_as_tool_results, compare_tools
from .consensus_reporting import write_consensus_report
from .integrations import (
    CasOffinderAdapter,
    CrispritzAdapter,
    ExternalImportAdapter,
    MafftAdapter,
    load_external_results,
)
from .provenance import sha256_file, write_run_manifest
from .sdk import load_run
from .tables import ToolResultTable

OUTPUT_NAMES = [
    "tool_results_long.csv",
    "candidate_tool_matrix.csv",
    "consensus_candidates.csv",
    "tool_coverage.csv",
    "model_agreement.csv",
    "disagreement_report.csv",
    "unmatched_external_rows.csv",
    "report.html",
    "run_manifest.json",
]


def _load_settings(path: str | Path) -> dict[str, Any]:
    return yaml.safe_load(Path(path).read_text(encoding="utf-8"))


def _verify_candidate_set(candidates: pd.DataFrame, settings: dict[str, Any]) -> pd.DataFrame:
    filter_settings = settings["candidate_filter"]
    column, expected_value = next(
        (key, value)
        for key, value in filter_settings.items()
        if key not in {"expected_total", "expected_genes"}
    )
    selected = candidates[
        pd.to_numeric(candidates[column], errors="coerce").eq(expected_value)
    ].copy()
    selected = selected.sort_values("candidate_id", kind="mergesort").reset_index(drop=True)
    expected_total = int(filter_settings["expected_total"])
    actual_genes = selected["gene_name"].value_counts().to_dict()
    expected_genes = filter_settings["expected_genes"]
    if len(selected) != expected_total or actual_genes != expected_genes:
        raise ValueError(
            f"HSV-2 consensus candidate invariant failed: expected {expected_total} "
            f"and genes {expected_genes}, observed {len(selected)} and {actual_genes}."
        )
    return selected


def _availability() -> pd.DataFrame:
    adapters = [MafftAdapter(), CasOffinderAdapter(), CrispritzAdapter()]
    adapters.extend(ExternalImportAdapter(tool) for tool in ("crispor", "chopchop", "guidescan2"))
    rows = []
    for adapter in adapters:
        status = adapter.detect()
        rows.append(
            {
                "tool_name": status.name,
                "available": status.available,
                "version": status.version,
                "execution_mode": status.execution_mode,
                "executable": status.executable,
                "message": status.message,
            }
        )
    return pd.DataFrame(rows)


def _find_export(directory: Path, tool: str) -> Path | None:
    for suffix in (".csv", ".tsv", ".json"):
        path = directory / f"{tool}{suffix}"
        if path.is_file():
            return path
    return None


def run_hsv2_consensus(
    config_path: str | Path = "configs/hsv2_consensus.yaml",
    output_directory: str | Path = "reports/hsv2_consensus",
) -> Path:
    settings = _load_settings(config_path)
    source_directory = Path(settings["source_run"])
    source_candidate_path = source_directory / "candidates_ranked_post_human.csv"
    run = load_run(source_directory)
    candidates = _verify_candidate_set(run.candidates, settings)
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    cache_directory = output / ".cache"
    cache_directory.mkdir(parents=True, exist_ok=True)
    candidate_path = output / "hsv2_candidates_32.csv"
    candidates.to_csv(candidate_path, index=False)

    crispritz_directory = output / "crispritz"
    crispritz_manifest = crispritz_directory / "crispritz_manifest.json"
    crispritz_guides = crispritz_directory / "guides.tsv"
    crispritz_parameters = settings["crispritz"]
    input_stamp = cache_directory / "crispritz_input.json"
    if stage_is_current(
        input_stamp,
        [crispritz_manifest, crispritz_guides],
        [source_candidate_path, config_path],
        crispritz_parameters,
    ):
        print("CRISPRitz input: cached")
    else:
        CrispritzAdapter(crispritz_parameters["docker_image"]).build_input(
            candidates, crispritz_parameters, crispritz_directory
        )
        write_stage_stamp(input_stamp, [source_candidate_path, config_path], crispritz_parameters)
        print("CRISPRitz input: generated")

    tool_result_frames = [
        candidate_metrics_as_tool_results(candidates, source_file=source_candidate_path)
    ]
    unmatched_frames = []
    optional_inputs: list[Path] = []
    crispritz_results = crispritz_directory / "tool_results_long.csv"
    if crispritz_results.is_file():
        tool_result_frames.append(ToolResultTable.from_frame(pd.read_csv(crispritz_results)))
        optional_inputs.append(crispritz_results)

    export_directory = Path(settings["external_exports_directory"])
    for tool in ("crispor", "chopchop", "guidescan2"):
        imported_path = output / tool / "tool_results_long.csv"
        unmatched_path = output / tool / "unmatched_external_rows.csv"
        export = _find_export(export_directory, tool)
        if export:
            mapping = Path("configs/import_maps") / f"{tool}.yaml"
            imported = load_external_results(tool, export, mapping, candidates)
            imported_path.parent.mkdir(parents=True, exist_ok=True)
            imported.results.dataframe.to_csv(imported_path, index=False)
            imported.raw_rows.to_csv(imported_path.parent / "raw_external_rows.csv", index=False)
            combined_unmatched = pd.concat(
                [imported.unmatched_rows, imported.ambiguous_rows], ignore_index=True
            )
            combined_unmatched.to_csv(unmatched_path, index=False)
            optional_inputs.extend([export, mapping])
        if imported_path.is_file():
            tool_result_frames.append(ToolResultTable.from_frame(pd.read_csv(imported_path)))
            optional_inputs.append(imported_path)
        if unmatched_path.is_file():
            unmatched_frames.append(pd.read_csv(unmatched_path))

    input_paths = [source_candidate_path, Path(config_path), *optional_inputs]
    final_outputs = [output / name for name in OUTPUT_NAMES]
    consensus_stamp = cache_directory / "consensus.json"
    signature_parameters = {
        "method": settings["consensus"]["method"],
        "weights": json.dumps(settings["consensus"]["weights"], sort_keys=True),
        "expected_tools": ",".join(settings["expected_tools"]),
    }
    if stage_is_current(consensus_stamp, final_outputs, input_paths, signature_parameters):
        print(f"HSV-2 consensus: cached ({output / 'report.html'})")
        return output / "report.html"

    unmatched = (
        pd.concat(unmatched_frames, ignore_index=True)
        if unmatched_frames
        else pd.DataFrame(columns=["tool_name", "mapping_status"])
    )
    comparison = compare_tools(
        candidates,
        tool_result_frames,
        method=settings["consensus"]["method"],
        weights=settings["consensus"]["weights"],
        unmatched_external_rows=unmatched,
        expected_tools=settings["expected_tools"],
    )
    comparison.write(output)
    availability = _availability()
    availability.to_csv(output / "tool_availability.csv", index=False)
    experimental_path = output / "experimental" / "crispresso2_measured_metrics.csv"
    experimental = pd.read_csv(experimental_path) if experimental_path.is_file() else pd.DataFrame()
    provenance = {
        "source_run": str(source_directory.resolve()),
        "source_candidate_sha256": sha256_file(source_candidate_path),
        "configuration": str(Path(config_path).resolve()),
        "candidate_count": len(candidates),
        "gene_counts": candidates["gene_name"].value_counts().to_dict(),
        "missing_results_are_not_zero_risk": True,
        "raw_scores_directly_averaged": False,
    }
    write_consensus_report(
        candidates,
        comparison,
        output / "report.html",
        tool_availability=availability,
        experimental_results=experimental,
        provenance=provenance,
    )
    write_run_manifest(
        output / "run_manifest.json",
        input_paths,
        signature_parameters,
        config_path=config_path,
        accepted_accessions=candidates["candidate_id"].astype(str).tolist(),
        human_assembly_identifier="GCF_000001405.40 / GRCh38.p14",
        command_line=["bash", "scripts/run_hsv2_consensus.sh"],
        random_seed=int(settings["random_seed"]),
        output_paths=[path for path in output.iterdir() if path.is_file()],
    )
    write_stage_stamp(consensus_stamp, input_paths, signature_parameters)
    print(f"HSV-2 consensus: wrote {output / 'report.html'}")
    return output / "report.html"
