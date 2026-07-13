"""CRISPRitz executable/import adapter with bulge and variant-aware summaries."""

from __future__ import annotations

import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from ..provenance import sha256_file
from ..tables import TOOL_RESULT_COLUMNS
from .base import AdapterError, ToolAvailability, ToolExecution, detect_executable


def _first(frame: pd.DataFrame, names: tuple[str, ...], default: object = "") -> pd.Series:
    lowered = {str(column).lower(): column for column in frame.columns}
    for name in names:
        if name.lower() in lowered:
            return frame[lowered[name.lower()]]
    return pd.Series([default] * len(frame), index=frame.index)


class CrispritzAdapter:
    name = "crispritz"

    def __init__(self, docker_image: str = "pinellolab/crispritz:latest") -> None:
        self.docker_image = docker_image

    def detect(self) -> ToolAvailability:
        for executable in ("crispritz.py", "crispritz"):
            result = detect_executable(self.name, executable, ("--help",))
            if result.available:
                return result
        docker = shutil.which("docker")
        if docker:
            return ToolAvailability(
                self.name,
                True,
                f"Docker image {self.docker_image} (not pulled automatically)",
                docker,
                "docker",
                "Docker is available; the configured image must already be present or be "
                "pulled by the researcher.",
            )
        return ToolAvailability(
            self.name,
            False,
            message=(
                "CRISPRitz was not found natively and Docker is unavailable. Input-only and import "
                "modes remain available; install CRISPRitz or Docker separately."
            ),
        )

    def build_input(self, candidates: pd.DataFrame, config: dict, output_dir: str | Path) -> Path:
        required = {"candidate_id", "guide_sequence"}
        missing = sorted(required - set(candidates.columns))
        if missing:
            raise AdapterError(f"CRISPRitz input is missing candidate columns: {missing}")
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        selected = candidates.copy().sort_values("candidate_id", kind="mergesort")
        guides = output / "guides.tsv"
        columns = [
            column
            for column in ["candidate_id", "guide_sequence", "pam", "gene_name"]
            if column in selected
        ]
        selected[columns].to_csv(guides, sep="\t", index=False)
        settings = {
            "schema_version": "0.4",
            "tool": self.name,
            "guides_file": guides.name,
            "candidate_count": len(selected),
            "reference_genome": config.get("reference_genome"),
            "genome_or_assembly": config.get("genome_or_assembly", ""),
            "pam_file": config.get("pam_file"),
            "mismatches": int(config.get("mismatches", 3)),
            "dna_bulges": int(config.get("dna_bulges", 0)),
            "rna_bulges": int(config.get("rna_bulges", 0)),
            "variant_aware": bool(config.get("variant_aware", False)),
            "variant_file": config.get("variant_file"),
            "docker_image": self.docker_image,
            "status": "input_ready",
        }
        manifest = output / "crispritz_manifest.json"
        manifest.write_text(json.dumps(settings, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return manifest

    def _command(self, manifest_path: Path, output_dir: Path) -> tuple[str, ...]:
        settings = json.loads(manifest_path.read_text(encoding="utf-8"))
        availability = self.detect()
        genome = settings.get("reference_genome")
        pam = settings.get("pam_file")
        if not genome or not pam:
            raise AdapterError(
                "CRISPRitz execution requires reference_genome and pam_file in the input manifest; "
                "the generated guides and manifest can still be exported/imported."
            )
        guides = manifest_path.parent / settings["guides_file"]
        arguments = (
            "search",
            str(genome),
            str(pam),
            str(guides),
            str(output_dir / "crispritz_results"),
            "-mm",
            str(settings["mismatches"]),
            "-bDNA",
            str(settings["dna_bulges"]),
            "-bRNA",
            str(settings["rna_bulges"]),
        )
        if settings.get("variant_aware"):
            arguments += ("--variants", str(settings.get("variant_file") or ""))
        if availability.execution_mode == "docker":
            root = manifest_path.parent.resolve()
            return (
                availability.executable or "docker",
                "run",
                "--rm",
                "-v",
                f"{root}:/work",
                self.docker_image,
                "crispritz.py",
                *arguments,
            )
        return (availability.executable or "crispritz.py", *arguments)

    def run(
        self, input_path: str | Path, output_dir: str | Path, *, dry_run: bool = False
    ) -> ToolExecution:
        manifest = Path(input_path)
        output = Path(output_dir)
        output.mkdir(parents=True, exist_ok=True)
        command = self._command(manifest, output)
        if dry_run:
            return ToolExecution(command, 0, "", "", output / "crispritz_results.targets.txt")
        availability = self.detect()
        if not availability.available:
            raise AdapterError(availability.message)
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        if completed.returncode:
            raise AdapterError(
                f"CRISPRitz failed ({completed.returncode}): "
                f"{(completed.stderr or completed.stdout).strip()}"
            )
        return ToolExecution(command, completed.returncode, completed.stdout, completed.stderr)

    def parse(self, output_path: str | Path, manifest: str | Path | None = None) -> pd.DataFrame:
        path = Path(output_path)
        if not path.is_file():
            raise AdapterError(
                f"CRISPRitz output is missing: {path}. Missing output is pending, not zero risk."
            )
        frame = pd.read_csv(path, sep=None, engine="python", comment="#")
        output = frame.copy()
        output["candidate_id"] = _first(frame, ("candidate_id", "Candidate", "ID"))
        output["guide_sequence"] = _first(frame, ("guide_sequence", "crRNA", "Guide"))
        output["mismatches"] = pd.to_numeric(
            _first(frame, ("mismatches", "MM", "Mismatch"), 0), errors="coerce"
        )
        output["bulge_type"] = _first(frame, ("bulge_type", "Bulge Type", "Bulge"), "none").astype(
            str
        )
        output["bulge_size"] = pd.to_numeric(
            _first(frame, ("bulge_size", "Bulge Size", "Bulge_Size"), 0), errors="coerce"
        ).fillna(0)
        output["annotation_type"] = _first(
            frame, ("annotation_type", "Annotation", "Gene Type"), ""
        )
        output["genome_mode"] = _first(
            frame, ("genome_mode", "Genome Mode", "Reference/Variant"), "reference-only"
        )
        output["population"] = _first(frame, ("population", "Population"), "")
        output["sample"] = _first(frame, ("sample", "Sample", "Sample ID"), "")
        if manifest and output["candidate_id"].fillna("").eq("").any():
            manifest_frame = pd.read_csv(manifest, sep=None, engine="python")
            mapping = manifest_frame.drop_duplicates("guide_sequence").set_index("guide_sequence")[
                "candidate_id"
            ]
            output.loc[output["candidate_id"].fillna("").eq(""), "candidate_id"] = output.loc[
                output["candidate_id"].fillna("").eq(""), "guide_sequence"
            ].map(mapping)
        return output

    def normalize(
        self,
        parsed_results: pd.DataFrame,
        *,
        candidates: pd.DataFrame | None = None,
        source_file: str | Path | None = None,
        version: str | None = None,
        assembly: str = "",
        editor: str = "",
        command: str = "imported CRISPRitz output",
    ) -> pd.DataFrame:
        if candidates is None:
            raise AdapterError("CRISPRitz normalization requires a candidate table.")
        source = Path(source_file) if source_file else None
        now = datetime.now(timezone.utc).isoformat()
        grouped = {
            key: group for key, group in parsed_results.groupby("candidate_id", dropna=False)
        }
        rows = []
        for _, candidate in candidates.iterrows():
            candidate_id = candidate["candidate_id"]
            group = grouped.get(candidate_id, parsed_results.iloc[0:0])
            counts = {
                f"mismatch_{value}": int(group["mismatches"].eq(value).sum())
                for value in sorted(group["mismatches"].dropna().astype(int).unique())
            }
            counts.update(
                {
                    "dna_bulge_hits": int(group["bulge_type"].str.upper().eq("DNA").sum())
                    if not group.empty
                    else 0,
                    "rna_bulge_hits": int(group["bulge_type"].str.upper().eq("RNA").sum())
                    if not group.empty
                    else 0,
                    "reference_hits": int(
                        group["genome_mode"].str.lower().str.contains("reference").sum()
                    )
                    if not group.empty
                    else 0,
                    "variant_enriched_hits": int(
                        group["genome_mode"].str.lower().str.contains("variant").sum()
                    )
                    if not group.empty
                    else 0,
                }
            )
            counts["annotation_counts"] = (
                group["annotation_type"]
                .fillna("unannotated")
                .replace("", "unannotated")
                .value_counts()
                .sort_index()
                .astype(int)
                .to_dict()
            )
            counts["populations"] = sorted(
                value for value in group["population"].dropna().astype(str).unique() if value
            )
            counts["samples"] = sorted(
                value for value in group["sample"].dropna().astype(str).unique() if value
            )
            burden = len(group)
            rows.append(
                {
                    "candidate_id": candidate_id,
                    "guide_sequence": candidate.get("guide_sequence", ""),
                    "gene_name": candidate.get("gene_name", ""),
                    "tool_name": self.name,
                    "tool_version": version or self.detect().version or "unknown",
                    "tool_mode": "variant-aware search"
                    if counts["variant_enriched_hits"]
                    else "reference-genome search",
                    "genome_or_assembly": assembly,
                    "editor": editor,
                    "metric_name": "predicted_offtarget_burden",
                    "raw_value": burden,
                    "normalized_value": pd.NA,
                    "rank": pd.NA,
                    "percentile_rank": pd.NA,
                    "decision": "requires expert review"
                    if burden
                    else "no predicted hit within imported CRISPRitz model and threshold",
                    "explanation": json.dumps(counts, sort_keys=True),
                    "source_file": str(source.resolve()) if source else "",
                    "source_file_sha256": sha256_file(source)
                    if source and source.is_file()
                    else "",
                    "command_or_import_method": command,
                    "timestamp": now,
                    "status": "completed",
                    "error_message": "",
                }
            )
        output = pd.DataFrame(rows, columns=TOOL_RESULT_COLUMNS)
        numeric = pd.to_numeric(output["raw_value"], errors="coerce")
        output["rank"] = numeric.rank(method="min", ascending=True)
        n = len(output)
        output["percentile_rank"] = 1.0 if n == 1 else 1.0 - (output["rank"] - 1) / (n - 1)
        output["normalized_value"] = output["percentile_rank"]
        return output
