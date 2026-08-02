"""Stable researcher-facing project experience.

This module orchestrates existing sequence analyses. It does not add a scientific
score or infer favorable values from missing data.
"""

from __future__ import annotations

import json
import os
import platform
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request
import webbrowser
import zipfile
from datetime import datetime, timezone
from html import escape
from pathlib import Path
from typing import Any

import pandas as pd
import yaml
from Bio import SeqIO

from . import __version__
from .crispr import scan_editor_candidates
from .io_utils import read_fasta
from .project_workflow import (
    STAGE_ORDER,
    ProjectContext,
    initialize_project,
    load_project,
    project_status,
)
from .provenance import sha256_file


def resource_path(*parts: str) -> Path:
    path = Path(__file__).resolve().parent / "resources"
    return path.joinpath(*parts)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _slug(value: str) -> str:
    cleaned = "".join(character.lower() if character.isalnum() else "-" for character in value)
    return "-".join(part for part in cleaned.split("-") if part) or "virus-project"


def _copy(source: str | Path, destination: Path) -> None:
    source_path = Path(source).expanduser().resolve()
    if not source_path.is_file():
        raise FileNotFoundError(source_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_path, destination)


def _download_text(url: str, destination: Path) -> None:
    request = urllib.request.Request(url, headers={"User-Agent": "ViralSafeTarget/0.10.0"})
    with urllib.request.urlopen(request, timeout=60) as response:
        content = response.read()
    if not content.strip():
        raise RuntimeError(f"NCBI returned an empty response: {url}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(content)


def _genbank_to_gff(genbank: Path, output: Path) -> int:
    """Convert NCBI GenBank annotations into a reference-matched compact GFF3."""
    rows = ["##gff-version 3"]
    feature_count = 0
    for record in SeqIO.parse(genbank, "genbank"):
        seqid = record.id
        for index, feature in enumerate(record.features, start=1):
            if feature.type not in {"gene", "CDS"}:
                continue
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            strand = "+" if feature.location.strand != -1 else "-"
            qualifiers = feature.qualifiers
            name = (qualifiers.get("gene") or qualifiers.get("locus_tag") or [""])[0]
            product = (qualifiers.get("product") or [""])[0]
            identifier = (qualifiers.get("protein_id") or [f"{feature.type.lower()}_{index}"])[0]
            attributes = {
                "ID": identifier,
                "Name": name or identifier,
                "product": product,
            }
            encoded = ";".join(
                f"{key}={urllib.parse.quote(str(value), safe='_.-')}"
                for key, value in attributes.items()
                if value
            )
            phase = "0" if feature.type == "CDS" else "."
            rows.append(
                "\t".join(
                    [
                        seqid,
                        "NCBI-GenBank",
                        feature.type,
                        str(start),
                        str(end),
                        ".",
                        strand,
                        phase,
                        encoded,
                    ]
                )
            )
            feature_count += 1
    if not feature_count:
        raise RuntimeError("No gene/CDS annotation was found in the NCBI GenBank record")
    output.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return feature_count


def acquire_ncbi_accession(accession: str, data_dir: Path) -> dict[str, Any]:
    encoded = urllib.parse.quote(accession)
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    fasta_url = f"{base}?db=nuccore&id={encoded}&rettype=fasta&retmode=text"
    genbank_url = f"{base}?db=nuccore&id={encoded}&rettype=gbwithparts&retmode=text"
    reference = data_dir / "reference.fasta"
    genbank = data_dir / "reference.gb"
    annotation = data_dir / "reference.gff3"
    _download_text(fasta_url, reference)
    _download_text(genbank_url, genbank)
    feature_count = _genbank_to_gff(genbank, annotation)
    records = read_fasta(reference)
    observed_accession = next(iter(records))
    strains = data_dir / "strains.aligned.fasta"
    shutil.copy2(reference, strains)
    return {
        "requested_accession": accession,
        "observed_accession": observed_accession,
        "retrieved_utc": _utc_now(),
        "sources": [
            {
                "role": "reference_fasta",
                "url": fasta_url,
                "path": str(reference),
                "sha256": sha256_file(reference),
            },
            {
                "role": "reference_genbank",
                "url": genbank_url,
                "path": str(genbank),
                "sha256": sha256_file(genbank),
            },
            {
                "role": "derived_gff3",
                "source": str(genbank),
                "path": str(annotation),
                "sha256": sha256_file(annotation),
            },
        ],
        "annotation_feature_count": feature_count,
        "strain_panel_note": (
            "Reference-only aligned panel; population conservation is not established."
        ),
    }


def _run_mafft(source: Path, destination: Path) -> None:
    executable = shutil.which("mafft")
    if not executable:
        raise RuntimeError(
            "Unaligned strain FASTA was supplied but MAFFT is unavailable. Install MAFFT "
            "and rerun, or provide --strains-aligned."
        )
    completed = subprocess.run(
        [executable, "--auto", str(source)], capture_output=True, text=True, check=False
    )
    if completed.returncode:
        raise RuntimeError(f"MAFFT failed: {completed.stderr.strip()}")
    destination.write_text(completed.stdout, encoding="utf-8")


def create_project(
    out_dir: str | Path,
    *,
    project_name: str,
    virus_name: str,
    tax_id: str | None = None,
    reference_accession: str | None = None,
    reference_fasta: str | Path | None = None,
    annotation_gff: str | Path | None = None,
    strains_fasta: str | Path | None = None,
    strains_aligned: bool = True,
    host_profile: str = "human_grch38",
    host_fasta: str | Path | None = None,
    nuclease_profile: str = "spcas9",
    evidence_enabled: bool = False,
    run_external: bool = False,
    sequence_only: bool = False,
    force: bool = False,
) -> Path:
    root = Path(out_dir).expanduser().resolve()
    identifier = _slug(project_name)
    project_file = initialize_project(
        root,
        project_id=identifier,
        display_name=virus_name,
        reference_accession=reference_accession or "LOCAL_REFERENCE",
        force=force,
    )
    provenance: dict[str, Any] = {"created_utc": _utc_now(), "inputs": []}
    data_dir = root / "data"
    if reference_accession and not reference_fasta:
        provenance["ncbi"] = acquire_ncbi_accession(reference_accession, data_dir)
    elif reference_fasta:
        _copy(reference_fasta, data_dir / "reference.fasta")
        provenance["inputs"].append(
            {
                "role": "reference_fasta",
                "source": str(Path(reference_fasta).expanduser()),
                "sha256": sha256_file(data_dir / "reference.fasta"),
            }
        )
    if annotation_gff:
        _copy(annotation_gff, data_dir / "reference.gff3")
    elif not (data_dir / "reference.gff3").is_file() and not sequence_only:
        raise ValueError(
            "Annotation is required unless --sequence-only is explicitly selected. "
            "Provide --annotation-gff or an accession with NCBI annotation."
        )
    if strains_fasta:
        source = Path(strains_fasta).expanduser().resolve()
        if strains_aligned:
            _copy(source, data_dir / "strains.aligned.fasta")
        else:
            _run_mafft(source, data_dir / "strains.aligned.fasta")
    elif (data_dir / "reference.fasta").is_file() and not (
        data_dir / "strains.aligned.fasta"
    ).is_file():
        shutil.copy2(data_dir / "reference.fasta", data_dir / "strains.aligned.fasta")
    if not (data_dir / "reference.fasta").is_file():
        raise ValueError("Provide either --reference-accession or --reference-fasta")

    project = yaml.safe_load(project_file.read_text(encoding="utf-8"))
    project["workflow"].update(
        {
            "sequence_only": bool(sequence_only),
            "evidence_discovery_enabled": bool(evidence_enabled),
            "run_external_by_default": bool(run_external),
        }
    )
    if sequence_only:
        project["analysis"]["enabled"] = False
    project_file.write_text(yaml.safe_dump(project, sort_keys=False), encoding="utf-8")

    virus_path = root / "profiles" / "virus.yaml"
    virus = yaml.safe_load(virus_path.read_text(encoding="utf-8"))
    records = read_fasta(data_dir / "reference.fasta")
    reference_id = next(iter(records))
    virus.update(
        {
            "scientific_name": virus_name,
            "tax_id": int(tax_id) if tax_id else None,
            "reference_accession": reference_id,
            "annotation_gff": None if sequence_only else "data/reference.gff3",
            "notes": "Sequence-only mode; gene/protein/evidence stages unavailable."
            if sequence_only
            else "Researcher-created project with explicit input provenance.",
        }
    )
    virus_path.write_text(yaml.safe_dump(virus, sort_keys=False), encoding="utf-8")

    host_path = root / "profiles" / "host.yaml"
    bundled_host = resource_path("profiles", "human_grch38.yaml")
    if host_profile == "human_grch38":
        shutil.copy2(bundled_host, host_path)
    host = yaml.safe_load(host_path.read_text(encoding="utf-8"))
    if host_fasta:
        _copy(host_fasta, root / "external" / "host" / "host.fasta")
        host["id"] = "local_host"
        host["display_name"] = "Researcher-supplied local host"
        host["assembly_name"] = "local_host"
        host["assembly_accession"] = "local"
    host_path.write_text(yaml.safe_dump(host, sort_keys=False), encoding="utf-8")
    if nuclease_profile != "spcas9":
        raise ValueError("Only the tested bundled spcas9 profile is currently selectable")
    shutil.copy2(resource_path("profiles", "spcas9.yaml"), root / "profiles" / "nuclease.yaml")
    (root / "input_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )
    return project_file


def create_demo_project(out_dir: str | Path, *, force: bool = False) -> Path:
    root = Path(out_dir).expanduser().resolve()
    return create_project(
        root,
        project_name="demo-project",
        virus_name="Synthetic demonstration virus",
        reference_fasta=resource_path("demo", "reference.fasta"),
        annotation_gff=resource_path("demo", "reference.gff3"),
        strains_fasta=resource_path("demo", "strains.aligned.fasta"),
        host_fasta=resource_path("demo", "host.fasta"),
        force=force,
    )


def _program(command: str, version_args: list[str]) -> dict[str, Any]:
    executable = shutil.which(command)
    if command == "cas-offinder":
        configured = os.environ.get("CAS_OFFINDER_PATH")
        if (
            configured
            and Path(configured).expanduser().is_file()
            and os.access(Path(configured).expanduser(), os.X_OK)
        ):
            executable = str(Path(configured).expanduser().resolve())
    if not executable:
        return {"available": False, "path": None, "version": "not found"}
    try:
        result = subprocess.run(
            [executable, *version_args], capture_output=True, text=True, timeout=10, check=False
        )
        lines = (result.stdout or result.stderr).strip().splitlines()
        version = lines[0] if lines else "available"
    except (OSError, subprocess.TimeoutExpired) as error:
        version = f"version unavailable: {error}"
    return {"available": True, "path": executable, "version": version}


def tool_setup_report(selected_tool: str = "all") -> dict[str, Any]:
    """Return actionable, non-mutating setup guidance for external programs."""
    system = platform.system().lower()
    managers = {
        "brew": bool(shutil.which("brew")),
        "apt": bool(shutil.which("apt-get")),
        "conda": bool(
            shutil.which("conda") or shutil.which("mamba") or shutil.which("micromamba")
        ),
        "docker": bool(shutil.which("docker")),
        "podman": bool(shutil.which("podman")),
    }
    status = doctor_report()["tools"]
    recipes: dict[str, dict[str, Any]] = {
        "mafft": {
            "name": "MAFFT",
            "required_for": "aligning an unaligned viral strain panel",
            "official_url": "https://mafft.cbrc.jp/alignment/software/",
            "commands": {
                "macos_homebrew": ["brew install mafft"],
                "ubuntu_debian": ["sudo apt-get update", "sudo apt-get install -y mafft"],
                "conda": ["conda install -c bioconda mafft"],
            },
        },
        "cas-offinder": {
            "name": "Cas-OFFinder",
            "required_for": "the configured genome-scale host-search stage",
            "official_url": "https://github.com/snugel/cas-offinder",
            "commands": {
                "source_build": [
                    "git clone https://github.com/snugel/cas-offinder.git",
                    "cmake -S cas-offinder -B cas-offinder/build",
                    "cmake --build cas-offinder/build --parallel",
                    "export CAS_OFFINDER_PATH=$PWD/cas-offinder/build/cas-offinder",
                ]
            },
            "note": (
                "OpenCL and a C++ build toolchain are required; follow the official "
                "repository for platform prerequisites."
            ),
        },
        "crispritz": {
            "name": "CRISPRitz",
            "required_for": (
                "optional independent host-search comparison, bulges, and "
                "variant-aware analyses"
            ),
            "official_url": "https://github.com/pinellolab/CRISPRitz",
            "commands": {
                "docker": ["docker pull pinellolab/crispritz:latest"],
                "native": [
                    "Follow the pinned release instructions in the official CRISPRitz repository"
                ],
            },
            "note": (
                "Input-only and import modes remain available; unavailable output remains "
                "pending, never zero."
            ),
        },
    }
    if selected_tool != "all" and selected_tool not in recipes:
        raise ValueError(f"Unknown tool: {selected_tool}")
    selected = recipes if selected_tool == "all" else {selected_tool: recipes[selected_tool]}
    tools: list[dict[str, Any]] = []
    for tool_id, recipe in selected.items():
        status_key = "cas_offinder" if tool_id == "cas-offinder" else tool_id
        detected = status[status_key]
        recommendation = "already available"
        if not detected["available"]:
            if tool_id == "mafft" and system == "darwin" and managers["brew"]:
                recommendation = "macos_homebrew"
            elif tool_id == "mafft" and managers["apt"]:
                recommendation = "ubuntu_debian"
            elif tool_id == "mafft" and managers["conda"]:
                recommendation = "conda"
            elif tool_id == "crispritz" and managers["docker"]:
                recommendation = "docker"
            elif tool_id == "crispritz":
                recommendation = "native"
            else:
                recommendation = "source_build"
        tools.append(
            {
                "id": tool_id,
                **recipe,
                "available": detected["available"],
                "detected_path": detected["path"],
                "detected_version": detected["version"],
                "recommended_recipe": recommendation,
            }
        )
    return {
        "platform": platform.platform(),
        "detected_package_managers": managers,
        "tools": tools,
        "next_check": "vst tools status",
        "safety": "Commands are printed for review and are never executed by vst tools setup.",
    }


def doctor_report(project: str | Path | None = None) -> dict[str, Any]:
    disk = shutil.disk_usage(Path.cwd())
    memory: int | None = None
    try:
        memory = int(os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES"))
    except (AttributeError, OSError, ValueError):
        pass
    tools = {
        "mafft": _program("mafft", ["--version"]),
        "cas_offinder": _program("cas-offinder", ["--help"]),
        "crispritz": _program("crispritz.py", ["--help"]),
        "docker": _program("docker", ["--version"]),
        "podman": _program("podman", ["--version"]),
        "ncbi_datasets": _program("datasets", ["version"]),
    }
    report: dict[str, Any] = {
        "python": {"executable": sys.executable, "version": platform.python_version()},
        "viral_safe_target_version": __version__,
        "platform": platform.platform(),
        "memory_gib": round(memory / 1024**3, 2) if memory else None,
        "disk_free_gib": round(disk.free / 1024**3, 2),
        "tools": tools,
        "scientific_boundary": "Missing external tools remain external_required, never zero hits.",
        "setup_help": "Run `vst tools setup` for platform-aware, non-mutating guidance.",
    }
    if project:
        context = load_project(project)
        host_root = context.profiles.resolve(context.profiles.host.get("fasta_root"))
        host_available = bool(
            host_root
            and (
                (host_root.is_file() and host_root.stat().st_size > 0)
                or (
                    host_root.is_dir()
                    and any(
                        path.is_file() and path.stat().st_size > 0
                        for path in host_root.rglob("*")
                        if path.suffix.lower() in {".fa", ".fna", ".fasta", ".fas"}
                    )
                )
            )
        )
        report["project"] = {
            "path": str(context.source),
            "host_reference_available": host_available,
            "can_run_sequence_stages": True,
            "can_run_external_host_stage": bool(
                host_available and tools["cas_offinder"]["available"]
            ),
        }
    return report


def plan_project(project: str | Path) -> dict[str, Any]:
    context = load_project(project)
    status = project_status(context)
    virus = context.profiles.virus
    alignment = context.profiles.resolve(virus.get("strain_alignment"))
    reference = context.profiles.resolve(virus.get("reference_fasta"))
    annotation = context.profiles.resolve(virus.get("annotation_gff"))
    host = context.profiles.resolve(context.profiles.host.get("fasta_root"))
    files = {
        "reference_fasta": reference,
        "annotation_gff": annotation,
        "strain_alignment": alignment,
        "host_reference": host,
    }
    validated = {
        name: {
            "path": str(path) if path else None,
            "status": "available" if path and path.exists() else "missing",
            "size_bytes": path.stat().st_size if path and path.is_file() else None,
        }
        for name, path in files.items()
    }
    candidate_estimate: int | None = None
    if alignment and alignment.is_file():
        records = read_fasta(alignment)
        reference_id = str(virus.get("reference_accession"))
        if reference_id in records:
            candidate_estimate = len(
                scan_editor_candidates(records, reference_id, context.profiles.editor, 0.0)
            )
    stage_rows = []
    timing_calibration = {
        "discovery_seconds_per_candidate": 0.015,
        "host_search_seconds_per_candidate_grch38_cpu": 0.814,
        "source": "committed synthetic and HSV-2 manifests; hardware-dependent",
    }
    total = 0.0
    unavailable_total = False
    existing = {row["stage"]: row for row in status["stages"]}
    for stage in STAGE_ORDER:
        cached = existing.get(stage, {}).get("status") == "completed"
        estimate: float | None
        confidence = "low"
        reason = "best-effort sequence-workflow calibration"
        if cached:
            estimate, confidence, reason = 0.0, "high", "cached output with valid state"
        elif stage == "host_screen":
            if not host or not host.exists():
                estimate, confidence, reason = None, "unavailable", "host reference unavailable"
            elif candidate_estimate is None:
                estimate, confidence, reason = None, "unavailable", "candidate estimate unavailable"
            else:
                estimate = (
                    candidate_estimate
                    * timing_calibration["host_search_seconds_per_candidate_grch38_cpu"]
                )
                confidence, reason = "low", "scaled from HSV-2 CPU source timing"
        elif candidate_estimate is None:
            estimate, confidence, reason = None, "unavailable", "no compatible calibration input"
        elif stage == "discover":
            estimate = max(
                1.0, candidate_estimate * timing_calibration["discovery_seconds_per_candidate"]
            )
        else:
            estimate = max(0.5, min(30.0, candidate_estimate * 0.004))
        if estimate is None:
            unavailable_total = True
        else:
            total += estimate
        stage_rows.append(
            {
                "stage": stage,
                "cached": cached,
                "estimated_seconds": round(estimate, 2) if estimate is not None else None,
                "confidence": confidence,
                "reason": reason,
            }
        )
    expected_host_workload = candidate_estimate
    estimated_disk = None if candidate_estimate is None else candidate_estimate * 25_000 + 5_000_000
    return {
        "schema_version": "1.0",
        "project": str(context.source),
        "validated_inputs": validated,
        "missing_inputs": [
            name for name, value in validated.items() if value["status"] == "missing"
        ],
        "candidate_count_estimate": candidate_estimate,
        "expected_host_search_queries": expected_host_workload,
        "estimated_disk_bytes": estimated_disk,
        "stages": stage_rows,
        "total_estimated_seconds": None if unavailable_total else round(total, 2),
        "total_estimate_confidence": "unavailable" if unavailable_total else "low",
        "assumptions": timing_calibration,
        "external_programs": doctor_report(project)["tools"],
        "execution_environment": {
            "cpu_count": os.cpu_count(),
            "cas_offinder_device": os.environ.get("CAS_OFFINDER_DEVICE", "C"),
            "container_mode": "available"
            if shutil.which("docker") or shutil.which("podman")
            else "unavailable",
        },
        "warning": "Runtime estimates are best effort and are not guarantees.",
    }


def _score_column(frame: pd.DataFrame) -> str | None:
    for column in ("post_human_score", "pre_human_score", "rank_score"):
        if column in frame:
            return column
    return None


def _portable_copy(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.is_file():
        shutil.copy2(source, destination)


def build_result_bundle(project: ProjectContext | str | Path) -> dict[str, Any]:
    context = load_project(project) if not isinstance(project, ProjectContext) else project
    root = context.output_root
    for directory in ("guide_explanations", "figures", "tables", "logs"):
        (root / directory).mkdir(parents=True, exist_ok=True)
    post = root / "host_screen" / "candidates_ranked_post_host.csv"
    pre = root / "discovery" / "discovery_panel.csv"
    candidates = pd.read_csv(post if post.is_file() else pre)
    score = _score_column(candidates)
    ordered = (
        candidates.sort_values(score, ascending=False, kind="mergesort") if score else candidates
    )
    top_guides = ordered.head(100).copy()
    shortlist = ordered.head(min(40, len(ordered))).copy()
    top_guides.to_csv(root / "top_guides.csv", index=False)
    shortlist.to_csv(root / "research_shortlist.csv", index=False)
    gene_column = (
        "gene_name"
        if "gene_name" in candidates
        else "feature_name"
        if "feature_name" in candidates
        else None
    )
    if gene_column:
        aggregations: dict[str, tuple[str, str]] = {"guide_count": ("candidate_id", "nunique")}
        if score:
            aggregations["best_guide_score"] = (score, "max")
        genes = (
            candidates.dropna(subset=[gene_column])
            .groupby(gene_column)
            .agg(**aggregations)
            .reset_index()
        )
        genes = genes.sort_values(list(aggregations), ascending=False, kind="mergesort")
    else:
        genes = pd.DataFrame(columns=["gene_name", "guide_count", "best_guide_score"])
    genes.to_csv(root / "top_genes.csv", index=False)
    virtual_dir = root / str(
        (context.values.get("analysis") or {}).get("output_dir", "virtual_knockout_escape")
    )
    multiplex_source = virtual_dir / "strategy_comparison.csv"
    multiplex = pd.read_csv(multiplex_source) if multiplex_source.is_file() else pd.DataFrame()
    multiplex.to_csv(root / "multiplex_panels.csv", index=False)
    evidence_source = root / "evidence" / "review_queue.tsv"
    if evidence_source.is_file():
        _portable_copy(evidence_source, root / "evidence_review_queue.tsv")
    else:
        pd.DataFrame(columns=["proposal_id", "gene_name", "review_status"]).to_csv(
            root / "evidence_review_queue.tsv", sep="\t", index=False
        )
    guide_virtual_path = virtual_dir / "guide_virtual_knockout.csv"
    guide_escape_path = virtual_dir / "guide_escape_robustness.csv"
    guide_virtual = (
        pd.read_csv(guide_virtual_path) if guide_virtual_path.is_file() else pd.DataFrame()
    )
    guide_escape = pd.read_csv(guide_escape_path) if guide_escape_path.is_file() else pd.DataFrame()
    for row in shortlist.itertuples(index=False):
        candidate_id = str(row.candidate_id)
        record = {column: getattr(row, column) for column in shortlist.columns}
        virtual = (
            guide_virtual[
                guide_virtual.get("candidate_id", pd.Series(dtype=str)).astype(str).eq(candidate_id)
            ]
            if not guide_virtual.empty
            else pd.DataFrame()
        )
        escape_rows = (
            guide_escape[
                guide_escape.get("candidate_id", pd.Series(dtype=str)).astype(str).eq(candidate_id)
            ]
            if not guide_escape.empty
            else pd.DataFrame()
        )
        explanation = {
            "candidate_id": candidate_id,
            "why_ranked": {
                key: value
                for key, value in record.items()
                if key.endswith("_score") or key.endswith("_penalty")
            },
            "conservation": record.get("conservation_score", "unknown"),
            "host_risk_status": "completed" if post.is_file() else "external_required",
            "predicted_host_hits": record.get("human_total_predicted_hits", "unknown"),
            "gene": record.get("gene_name", record.get("feature_name", "unknown")),
            "protein_mapping": virtual.to_dict("records"),
            "escape_summary": escape_rows.to_dict("records"),
            "direct_virus_evidence": record.get("direct_virus_evidence", "unknown"),
            "ortholog_evidence": record.get("ortholog_evidence", "unknown"),
            "limitations": "Computational explanation only; missing fields remain unknown.",
            "source_candidate_table": str(
                (post if post.is_file() else pre).relative_to(context.root)
            ),
        }
        explanation_path = root / "guide_explanations" / f"{candidate_id}.json"
        explanation_path.write_text(
            json.dumps(explanation, indent=2, default=str) + "\n", encoding="utf-8"
        )
        html_path = explanation_path.with_suffix(".html")
        html_path.write_text(
            "<!doctype html><html lang='en'><meta charset='utf-8'><title>Guide explanation</title>"
            f"<h1>{escape(candidate_id)}</h1><pre>"
            f"{escape(json.dumps(explanation, indent=2, default=str))}</pre>",
            encoding="utf-8",
        )
    state = project_status(context)
    summary = {
        "schema_version": "1.0",
        "project_id": context.values["id"],
        "created_utc": _utc_now(),
        "candidate_count": len(candidates),
        "shortlist_count": len(shortlist),
        "gene_count": len(genes),
        "host_risk_status": next(
            (row["status"] for row in state["stages"] if row["stage"] == "host_screen"), "pending"
        ),
        "stage_status": state["stages"],
        "scientific_boundary": (
            "Candidates for independent evaluation; not treatment recommendations."
        ),
    }
    (root / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (root / "SUMMARY.md").write_text(
        f"# {context.values.get('display_name', context.values['id'])} research summary\n\n"
        f"- Candidate rows: {len(candidates)}\n- Research shortlist: {len(shortlist)}\n"
        f"- Gene rows: {len(genes)}\n- Host-risk stage: {summary['host_risk_status']}\n\n"
        "These outputs are computational hypotheses for independent evaluation.\n",
        encoding="utf-8",
    )
    links = [
        "top_guides.csv",
        "top_genes.csv",
        "research_shortlist.csv",
        "multiplex_panels.csv",
        "evidence_review_queue.tsv",
        "stage_timings.json",
        "run_manifest.json",
        "export.zip",
    ]
    stage_html = "".join(
        f"<tr><td>{escape(row['stage'])}</td><td><span class='status {escape(row['status'])}'>"
        f"{escape(row['status'])}</span></td><td>{escape(row['message'])}</td></tr>"
        for row in state["stages"]
    )
    top_guide_columns = [
        column
        for column in ("candidate_id", "guide_sequence", "gene_name", score)
        if column and column in top_guides
    ]
    top_gene_columns = [column for column in genes.columns if column][:4]
    top_guide_display = top_guides[top_guide_columns].head(10).copy()
    if "candidate_id" in top_guide_display:
        top_guide_display["candidate_id"] = top_guide_display["candidate_id"].map(
            lambda value: (
                f"<a href='guide_explanations/{escape(str(value))}.html'>"
                f"{escape(str(value))}</a>"
            )
        )
    top_guide_html = top_guide_display.to_html(index=False, escape=False)
    top_gene_html = genes[top_gene_columns].head(10).to_html(index=False, escape=True)
    missing_stages = [row["stage"] for row in state["stages"] if row["status"] != "completed"]
    link_html = "".join(f"<li><a href='{name}'>{name}</a></li>" for name in links)
    next_actions: list[str] = []
    if summary["host_risk_status"] == "external_required":
        next_actions.append(
            "Install or configure Cas-OFFinder (`vst tools setup --tool cas-offinder`), "
            f"then run `vst resume {context.source} --run-external`."
        )
    if missing_stages:
        next_actions.append(
            "Review the stage table before interpreting rankings; partial or unavailable stages "
            "must not be treated as favorable results."
        )
    next_actions.append(
        "Open the evidence review queue and obtain named human review before biological evidence "
        "affects interpretation."
    )
    next_actions_html = "".join(f"<li>{escape(action)}</li>" for action in next_actions)
    (root / "START_HERE.html").write_text(
        "<!doctype html><html lang='en'><meta charset='utf-8'>"
        "<title>ViralSafeTarget results</title>"
        "<style>:root{color-scheme:light}body{font:16px system-ui,-apple-system,sans-serif;"
        "max-width:1120px;margin:2rem auto;padding:0 1rem;line-height:1.55;color:#172033}"
        "h1,h2{line-height:1.2}a{color:#0756a8}table{border-collapse:collapse;width:100%;"
        "display:block;overflow-x:auto}td,th{border:1px solid #d7dce5;padding:.55rem;"
        "text-align:left}th{background:#f3f6fa}.warning{background:#fff4d6;padding:1rem;"
        "border-left:5px solid #d99b20}.cards{display:grid;grid-template-columns:repeat(3,1fr);"
        "gap:1rem;margin:1.25rem 0}.card{background:#f3f6fa;border-radius:.5rem;padding:1rem}"
        ".card b{display:block;font-size:1.5rem}.status{border-radius:1rem;padding:.2rem .55rem;"
        "font-size:.85rem;background:#e8edf4}.status.completed{background:#dff4e6;color:#155d32}"
        ".status.external_required,.status.unavailable{background:#fff0cf;color:#704d00}"
        "code,pre{background:#f5f6f8;border-radius:.35rem;padding:.15rem .3rem}"
        "pre{padding:1rem;overflow:auto}@media(max-width:700px){.cards{grid-template-columns:1fr}}"
        "</style>"
        f"<h1>{escape(str(context.values.get('display_name', context.values['id'])))}</h1>"
        "<div class='warning'><b>Research boundary.</b> This is a computational research "
        "shortlist, "
        "not evidence of editing, safety, viral inhibition, treatment, or cure.</div>"
        "<div class='cards'>"
        f"<div class='card'><b>{len(candidates)}</b>candidate rows</div>"
        f"<div class='card'><b>{len(shortlist)}</b>research shortlist</div>"
        f"<div class='card'><b>{escape(summary['host_risk_status'])}</b>host-risk stage</div>"
        "</div>"
        "<h2>Research question and inputs</h2><p>Which conserved editor-compatible "
        "viral sites merit independent evaluation after preserving host-risk status, gene "
        "targetability, disruption hypotheses, escape robustness, and evidence as separate "
        "axes?</p>"
        f"<p><b>Project:</b> {escape(str(context.source))}<br><b>Reference:</b> "
        f"{escape(str(context.profiles.virus.get('reference_accession', 'unknown')))}</p>"
        "<h2>Candidate funnel</h2>"
        f"<p>{len(candidates)} rows were available at the final completed sequence stage; "
        f"{len(top_guides)} are shown in top_guides.csv and {len(shortlist)} form the research "
        "shortlist.</p>"
        f"<h2>Top guides</h2>{top_guide_html}<h2>Top genes</h2>{top_gene_html}"
        "<h2>Stage status</h2><table><tr><th>Stage</th><th>Status</th><th>Message</th></tr>"
        f"{stage_html}</table><h2>Start with these files</h2><ul>{link_html}</ul>"
        f"<h2>Missing or partial stages</h2><p>{escape(', '.join(missing_stages) or 'None')}</p>"
        f"<h2>What to do next</h2><ol>{next_actions_html}</ol>"
        "<h2>Detailed analyses</h2><p>The linked virtual-knockout/escape report and "
        "multiplex table preserve predicted disruption and exact-target escape as separate "
        "sequence-level hypotheses. The evidence queue requires explicit human approval.</p>"
        f"<h2>Reproduce</h2><pre>vst plan {escape(str(context.source))}\n"
        f"vst run {escape(str(context.source))}\n"
        f"vst export {escape(str(context.source))}</pre>"
        "<h2>Interpretation</h2><p>Guide quality, gene targetability, host-risk status, "
        "predicted disruption, evidence, and escape robustness remain separate axes. Unknowns "
        "are not converted to favorable values.</p></html>",
        encoding="utf-8",
    )
    export_project(context.source, output=root / "export.zip")
    return {"results": str(root), "start_here": str(root / "START_HERE.html"), **summary}


def export_project(
    project: str | Path,
    *,
    output: str | Path | None = None,
    include_large_raw: bool = False,
) -> Path:
    context = load_project(project)
    destination = Path(output).resolve() if output else context.output_root / "export.zip"
    destination.parent.mkdir(parents=True, exist_ok=True)
    excluded_parts = {"external", "raw", "reference.fasta", "reference.gb"}
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        sources = [context.source, *context.profiles.source_paths, context.ranking_config]
        for source in sources:
            if source.is_file():
                archive.write(source, f"project/{source.name}")
        for path in context.output_root.rglob("*"):
            if not path.is_file() or path.resolve() == destination.resolve():
                continue
            relative = path.relative_to(context.output_root)
            if not include_large_raw and (
                set(relative.parts) & excluded_parts or path.stat().st_size > 25_000_000
            ):
                continue
            archive.write(path, f"results/{relative}")
    return destination


def open_results(path: str | Path, *, no_browser: bool = False) -> Path:
    source = Path(path).expanduser().resolve()
    if source.name == "project.yaml":
        source = load_project(source).output_root
    report = source if source.is_file() else source / "START_HERE.html"
    if not report.is_file():
        raise FileNotFoundError(f"START_HERE.html not found: {report}")
    if not no_browser:
        webbrowser.open(report.as_uri())
    return report
