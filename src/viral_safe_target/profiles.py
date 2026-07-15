"""Configuration-driven virus, host, and nuclease research profiles."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from .config import EditorProfile

PROFILE_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class ResearchProfileBundle:
    virus: dict[str, Any]
    host: dict[str, Any]
    nuclease: dict[str, Any]
    project_root: Path
    source_paths: tuple[Path, Path, Path]

    @property
    def editor(self) -> EditorProfile:
        return EditorProfile.from_mapping(self.nuclease)

    def resolve(self, value: str | None) -> Path | None:
        if not value:
            return None
        path = Path(value)
        return path if path.is_absolute() else (self.project_root / path).resolve()


def _read_profile(path: str | Path, expected_type: str) -> tuple[dict[str, Any], Path]:
    source = Path(path).resolve()
    with source.open(encoding="utf-8") as handle:
        values = yaml.safe_load(handle) or {}
    if values.get("schema_version") != PROFILE_SCHEMA_VERSION:
        raise ValueError(f"{source} must use profile schema_version {PROFILE_SCHEMA_VERSION!r}")
    if values.get("profile_type") != expected_type:
        raise ValueError(f"{source} is not a {expected_type!r} profile")
    identifier = str(values.get("id", ""))
    if not identifier or any(character.isspace() for character in identifier):
        raise ValueError(f"{source} has an invalid or missing profile id")
    return values, source


def load_profile_bundle(
    virus_profile: str | Path,
    host_profile: str | Path,
    nuclease_profile: str | Path,
    *,
    project_root: str | Path | None = None,
) -> ResearchProfileBundle:
    virus, virus_path = _read_profile(virus_profile, "virus")
    host, host_path = _read_profile(host_profile, "host")
    nuclease, nuclease_path = _read_profile(nuclease_profile, "nuclease")
    root = Path(project_root).resolve() if project_root else Path.cwd().resolve()
    bundle = ResearchProfileBundle(
        virus=virus,
        host=host,
        nuclease=nuclease,
        project_root=root,
        source_paths=(virus_path, host_path, nuclease_path),
    )
    bundle.editor.validate()
    return bundle


def validate_profile_bundle(
    bundle: ResearchProfileBundle,
    *,
    require_large_host_reference: bool = False,
    require_virus_inputs: bool = False,
) -> pd.DataFrame:
    checks: list[dict[str, object]] = []

    def display_path(path: Path) -> str:
        try:
            return str(path.relative_to(bundle.project_root))
        except ValueError:
            return str(path)

    def add(component: str, status: str, detail: str) -> None:
        checks.append({"component": component, "status": status, "detail": detail})

    required_virus = {
        "reference_accession",
        "reference_fasta",
        "annotation_gff",
        "strain_alignment",
    }
    missing = sorted(required_virus - bundle.virus.keys())
    add(
        "virus required fields",
        "pass" if not missing else "fail",
        "complete" if not missing else "missing: " + ", ".join(missing),
    )
    required_input_fields = {"reference_fasta", "annotation_gff", "strain_alignment"}
    for field in (
        "reference_fasta",
        "reference_genbank",
        "annotation_gff",
        "strain_alignment",
        "evidence_table",
        "domain_table",
        "disorder_table",
        "conserved_region_table",
        "gene_category_table",
        "external_validation_table",
        "population_validation_candidates",
        "population_validation_genes",
    ):
        path = bundle.resolve(bundle.virus.get(field))
        if path is None:
            add(f"virus path: {field}", "optional_missing", "not configured")
        else:
            exists = path.exists()
            if exists:
                status = "pass"
            elif field in required_input_fields:
                status = "fail" if require_virus_inputs else "input_pending"
            else:
                status = "optional_missing"
            add(
                f"virus path: {field}",
                status,
                display_path(path),
            )
    host_root = bundle.resolve(bundle.host.get("fasta_root"))
    host_exists = False
    if host_root and host_root.is_file():
        host_exists = host_root.stat().st_size > 0
    elif host_root and host_root.is_dir():
        host_exists = any(
            path.is_file() and path.stat().st_size > 0
            for path in host_root.rglob("*")
            if path.suffix.lower() in {".fa", ".fna", ".fasta", ".fas"}
        )
    add(
        "host reference",
        "pass" if host_exists else "fail" if require_large_host_reference else "external_pending",
        display_path(host_root) if host_root else "not configured",
    )
    add(
        "nuclease profile",
        "pass",
        f"{bundle.editor.name}; PAM={bundle.editor.pam_pattern}; tested={bundle.editor.tested}",
    )
    evidence = bundle.resolve(bundle.virus.get("evidence_table"))
    if evidence and evidence.is_file():
        table = pd.read_csv(evidence, sep="\t")
        essential = {
            "gene_name",
            "virus_type",
            "essentiality_call",
            "evidence_strength",
            "source_identifier",
            "source_url",
        }
        absent = sorted(essential - set(table))
        uncited = int(table.get("source_url", pd.Series(dtype=str)).fillna("").eq("").sum())
        add(
            "gene evidence schema",
            "pass" if not absent and not uncited else "fail",
            f"rows={len(table)}; missing_columns={absent}; uncited_rows={uncited}",
        )
    return pd.DataFrame(checks)
