"""Single-entry, configuration-driven research project workflow.

The project workflow deliberately stops at computational hypotheses. It does
not provide wet-lab instructions or claim editing, safety, efficacy, viral
inactivation, latency clearance, or clinical utility.
"""

from __future__ import annotations

import copy
import hashlib
import json
import os
import shutil
import subprocess
from collections.abc import Callable
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from .annotations import annotate_candidates, read_gff3
from .config import DEFAULT_CONFIG_PATH, load_config
from .crispr import scan_editor_candidates
from .discovery import build_candidate_feature_map, select_balanced_discovery_panel
from .disruption import rank_candidate_pairs
from .io_utils import read_fasta, require_aligned
from .offtarget import (
    build_cas_offinder_input,
    read_cas_offinder_output,
    summarize_cas_offinder_hits,
)
from .profiles import ResearchProfileBundle, load_profile_bundle, validate_profile_bundle
from .provenance import sha256_file
from .reporting import write_html_report, write_methods_and_limitations
from .scoring import rank_post_human_candidates, rank_pre_human_candidates

PROJECT_SCHEMA_VERSION = "1.0"
PROJECT_TYPE = "viral_safe_target_project"
STAGE_ORDER = ("validate", "discover", "host_screen", "pairs", "report")

SOURCE_LINKED_EVIDENCE_COLUMNS = [
    "gene_name",
    "virus_type",
    "reference_accession",
    "evidence_category",
    "essentiality_call",
    "essentiality_score",
    "evidence_strength",
    "experimental_system",
    "finding",
    "source_identifier",
    "source_title",
    "source_url",
    "directness_notes",
]


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass(frozen=True)
class ProjectContext:
    source: Path
    root: Path
    values: dict[str, Any]
    profiles: ResearchProfileBundle
    ranking_config: Path
    output_root: Path

    def resolve(self, value: str | Path) -> Path:
        path = Path(value)
        return path if path.is_absolute() else (self.root / path).resolve()


def _read_yaml(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        values = yaml.safe_load(handle) or {}
    if not isinstance(values, dict):
        raise ValueError(f"YAML root must be a mapping: {path}")
    return values


def load_project(path: str | Path) -> ProjectContext:
    source = Path(path).resolve()
    values = _read_yaml(source)
    if values.get("schema_version") != PROJECT_SCHEMA_VERSION:
        raise ValueError(f"{source} must use project schema_version {PROJECT_SCHEMA_VERSION!r}")
    if values.get("project_type") != PROJECT_TYPE:
        raise ValueError(f"{source} is not a {PROJECT_TYPE!r} file")
    identifier = str(values.get("id", ""))
    if not identifier or any(character.isspace() for character in identifier):
        raise ValueError("Project id is missing or contains whitespace")
    profiles = values.get("profiles") or {}
    required_profiles = {"virus", "host", "nuclease"}
    missing = sorted(required_profiles - set(profiles))
    if missing:
        raise ValueError("Project is missing profile paths: " + ", ".join(missing))
    root_value = str(values.get("project_root", "."))
    root_path = Path(root_value)
    root = root_path.resolve() if root_path.is_absolute() else (source.parent / root_path).resolve()

    def project_path(value: str) -> Path:
        candidate = Path(value)
        return candidate if candidate.is_absolute() else (root / candidate).resolve()

    bundle = load_profile_bundle(
        project_path(str(profiles["virus"])),
        project_path(str(profiles["host"])),
        project_path(str(profiles["nuclease"])),
        project_root=root,
    )
    workflow = values.get("workflow") or {}
    ranking_config = project_path(
        str(workflow.get("ranking_config", DEFAULT_CONFIG_PATH.resolve()))
    )
    output_root = project_path(str(workflow.get("output_root", "results")))
    return ProjectContext(source, root, values, bundle, ranking_config, output_root)


def _yaml_text(values: dict[str, Any]) -> str:
    return yaml.safe_dump(values, sort_keys=False, allow_unicode=True)


def initialize_project(
    out_dir: str | Path,
    *,
    project_id: str,
    display_name: str,
    reference_accession: str = "CHANGE_ME",
    force: bool = False,
) -> Path:
    """Create a self-contained, editable new-virus project skeleton."""
    if not project_id or any(character.isspace() for character in project_id):
        raise ValueError("project_id must be non-empty and contain no whitespace")
    root = Path(out_dir).resolve()
    project_file = root / "project.yaml"
    if project_file.exists() and not force:
        raise FileExistsError(f"Project already exists: {project_file}; pass --force to replace")
    for directory in (
        root / "data",
        root / "evidence",
        root / "external" / "host",
        root / "profiles",
        root / "results",
    ):
        directory.mkdir(parents=True, exist_ok=True)

    project = {
        "schema_version": PROJECT_SCHEMA_VERSION,
        "project_type": PROJECT_TYPE,
        "id": project_id,
        "display_name": display_name,
        "project_root": ".",
        "profiles": {
            "virus": "profiles/virus.yaml",
            "host": "profiles/host.yaml",
            "nuclease": "profiles/nuclease.yaml",
        },
        "workflow": {
            "ranking_config": "profiles/ranking.yaml",
            "output_root": "results",
            "minimum_strain_coverage": 0.95,
            "balanced_top_per_gene": 50,
            "global_top": 500,
            "maximum_host_candidates": None,
        },
    }
    virus = {
        "schema_version": "1.0",
        "profile_type": "virus",
        "id": project_id,
        "display_name": display_name,
        "reference_accession": reference_accession,
        "reference_fasta": "data/reference.fasta",
        "reference_genbank": None,
        "annotation_gff": "data/reference.gff3",
        "strain_alignment": "data/strains.aligned.fasta",
        "evidence_table": "evidence/gene_evidence.tsv",
        "domain_table": None,
        "disorder_table": None,
        "external_validation_table": None,
        "circular_genome": False,
        "notes": "Replace CHANGE_ME values and add reference-matched research inputs.",
    }
    host = {
        "schema_version": "1.0",
        "profile_type": "host",
        "id": "host",
        "display_name": "Host assembly — edit this profile",
        "assembly_name": "CHANGE_ME",
        "assembly_accession": "CHANGE_ME",
        "fasta_root": "external/host",
        "annotation_gff": None,
    }
    nuclease = {
        "schema_version": "1.0",
        "profile_type": "nuclease",
        "id": "spcas9",
        "display_name": "Streptococcus pyogenes Cas9",
        "name": "SpCas9",
        "protospacer_length": 20,
        "pam_pattern": "NGG",
        "pam_orientation": "3prime",
        "cut_offset": 3,
        "mismatch_search_threshold": 3,
        "tested": True,
        "notes": "Bundled scanner validation currently covers this 3-prime-PAM profile.",
    }
    project_file.write_text(_yaml_text(project), encoding="utf-8")
    (root / "profiles" / "virus.yaml").write_text(_yaml_text(virus), encoding="utf-8")
    (root / "profiles" / "host.yaml").write_text(_yaml_text(host), encoding="utf-8")
    (root / "profiles" / "nuclease.yaml").write_text(_yaml_text(nuclease), encoding="utf-8")
    shutil.copyfile(DEFAULT_CONFIG_PATH, root / "profiles" / "ranking.yaml")
    pd.DataFrame(columns=SOURCE_LINKED_EVIDENCE_COLUMNS).to_csv(
        root / "evidence" / "gene_evidence.tsv", sep="\t", index=False
    )
    (root / "data" / "README.md").write_text(
        "# Required inputs\n\n"
        "- `reference.fasta`: the frozen reference record.\n"
        "- `reference.gff3`: annotation using the same reference identifier and coordinates.\n"
        "- `strains.aligned.fasta`: an already aligned multi-FASTA containing the reference.\n",
        encoding="utf-8",
    )
    (root / "README.md").write_text(
        f"# {display_name}\n\n"
        "1. Add the three files documented in `data/README.md`.\n"
        "2. Edit the profiles and source-linked evidence table.\n"
        "3. Run `vst project validate --project project.yaml`.\n"
        "4. Run `vst project run --project project.yaml`.\n"
        "5. If the host stage is external, run Cas-OFFinder and then use "
        "`vst project resume --project project.yaml`.\n\n"
        "Outputs are computational hypotheses, not validated interventions.\n",
        encoding="utf-8",
    )
    return project_file


def _add_check(rows: list[dict[str, str]], component: str, status: str, detail: str) -> None:
    rows.append({"component": component, "status": status, "detail": detail})


def validate_project(
    project: ProjectContext | str | Path,
    *,
    require_host_reference: bool = False,
) -> pd.DataFrame:
    context = load_project(project) if not isinstance(project, ProjectContext) else project
    rows = validate_profile_bundle(
        context.profiles,
        require_large_host_reference=require_host_reference,
        require_virus_inputs=True,
    ).to_dict("records")
    _add_check(rows, "project specification", "pass", str(context.source))
    _add_check(
        rows,
        "ranking configuration",
        "pass" if context.ranking_config.is_file() else "fail",
        str(context.ranking_config),
    )
    alignment_path = context.profiles.resolve(context.profiles.virus.get("strain_alignment"))
    reference_path = context.profiles.resolve(context.profiles.virus.get("reference_fasta"))
    gff_path = context.profiles.resolve(context.profiles.virus.get("annotation_gff"))
    reference_id = str(context.profiles.virus.get("reference_accession", ""))
    if alignment_path and alignment_path.is_file():
        try:
            records = read_fasta(alignment_path)
            require_aligned(records)
            status = "pass" if reference_id in records else "fail"
            detail = f"records={len(records)}; reference={reference_id!r}"
            _add_check(rows, "aligned viral panel", status, detail)
            if reference_path and reference_path.is_file() and reference_id in records:
                reference_records = read_fasta(reference_path)
                reference_sequence = next(
                    (
                        sequence
                        for accession, sequence in reference_records.items()
                        if accession == reference_id or accession.startswith(reference_id)
                    ),
                    None,
                )
                matches = reference_sequence is not None and (
                    records[reference_id].replace("-", "").upper()
                    == reference_sequence.replace("-", "").upper()
                )
                _add_check(
                    rows,
                    "reference/alignment identity",
                    "pass" if matches else "fail",
                    "ungapped aligned reference matches reference FASTA"
                    if matches
                    else "reference record missing or nucleotide content differs",
                )
        except (ValueError, KeyError) as error:
            _add_check(rows, "aligned viral panel", "fail", str(error))
    if gff_path and gff_path.is_file():
        features = read_gff3(gff_path)
        seqids = sorted(features.get("seqid", pd.Series(dtype=str)).astype(str).unique())
        _add_check(
            rows,
            "annotation/reference identity",
            "pass" if reference_id in seqids else "fail",
            f"reference={reference_id!r}; GFF seqids={seqids[:5]}",
        )
    return pd.DataFrame(rows, columns=["component", "status", "detail"])


def _settings(context: ProjectContext) -> dict[str, Any]:
    settings = copy.deepcopy(load_config(context.ranking_config))
    editor = context.profiles.editor
    settings["editor"] = {
        "name": editor.name,
        "protospacer_length": editor.protospacer_length,
        "pam_pattern": editor.pam_pattern,
        "pam_orientation": editor.pam_orientation,
        "cut_offset": editor.cut_offset,
        "mismatch_search_threshold": editor.mismatch_search_threshold,
        "tested": editor.tested,
        "notes": editor.notes,
    }
    settings["off_target"]["human_assembly"] = context.profiles.host.get(
        "assembly_name", "unknown host assembly"
    )
    settings["off_target"]["human_assembly_accession"] = context.profiles.host.get(
        "assembly_accession", ""
    )
    return settings


def _signature(paths: list[Path], parameters: dict[str, Any]) -> str:
    digest = hashlib.sha256(json.dumps(parameters, sort_keys=True, default=str).encode())
    for path in sorted(paths, key=lambda item: str(item)):
        digest.update(str(path.resolve()).encode())
        digest.update(sha256_file(path).encode() if path.is_file() else b"MISSING")
    return digest.hexdigest()


class _State:
    def __init__(self, context: ProjectContext, *, restart: bool = False):
        self.context = context
        self.path = context.output_root / "workflow_state.json"
        context.output_root.mkdir(parents=True, exist_ok=True)
        if self.path.is_file() and not restart:
            self.values = json.loads(self.path.read_text(encoding="utf-8"))
        else:
            self.values = {
                "schema_version": "1.0",
                "project_id": context.values["id"],
                "project_file": str(context.source),
                "created_utc": _utc_now(),
                "stages": {},
            }
        self.write()

    def write(self) -> None:
        self.values["updated_utc"] = _utc_now()
        self.path.write_text(json.dumps(self.values, indent=2, sort_keys=True), encoding="utf-8")

    def reusable(self, stage: str, signature: str, outputs: list[Path]) -> bool:
        record = self.values["stages"].get(stage, {})
        return (
            record.get("status") == "completed"
            and record.get("signature") == signature
            and all(path.is_file() for path in outputs)
        )

    def set(
        self,
        stage: str,
        status: str,
        *,
        signature: str = "",
        outputs: list[Path] | None = None,
        message: str = "",
    ) -> None:
        self.values["stages"][stage] = {
            "status": status,
            "signature": signature,
            "outputs": [str(path.resolve()) for path in outputs or []],
            "message": message,
            "updated_utc": _utc_now(),
        }
        self.write()


def _run_cached_stage(
    state: _State,
    stage: str,
    signature: str,
    outputs: list[Path],
    action: Callable[[], None],
) -> str:
    if state.reusable(stage, signature, outputs):
        return "cached"
    state.set(stage, "running", signature=signature, outputs=outputs)
    try:
        action()
    except Exception as error:
        state.set(stage, "failed", signature=signature, outputs=outputs, message=str(error))
        raise
    state.set(stage, "completed", signature=signature, outputs=outputs)
    return "completed"


def _cas_offinder_executable() -> str | None:
    configured = os.environ.get("CAS_OFFINDER_PATH")
    if configured and Path(configured).is_file() and os.access(configured, os.X_OK):
        return configured
    return shutil.which("cas-offinder")


def _host_reference_available(path: Path | None) -> bool:
    if path is None:
        return False
    if path.is_file():
        return path.stat().st_size > 0
    if path.is_dir():
        return any(
            candidate.is_file() and candidate.stat().st_size > 0
            for candidate in path.rglob("*")
            if candidate.suffix.lower() in {".fa", ".fna", ".fasta", ".fas"}
        )
    return False


def _run_cas_offinder(executable: str, input_path: Path, output_path: Path) -> None:
    completed = subprocess.run(
        [executable, str(input_path), "C", str(output_path)],
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(
            f"Cas-OFFinder failed ({completed.returncode}): "
            f"{(completed.stderr or completed.stdout).strip()}"
        )
    if not output_path.is_file():
        raise RuntimeError("Cas-OFFinder returned success without creating its output")


def _read_optional_csv(path: Path) -> pd.DataFrame:
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def run_project(
    project: str | Path,
    *,
    run_external: bool = False,
    restart: bool = False,
    stop_after: str | None = None,
) -> dict[str, Any]:
    """Run or resume the generic sequence-level project workflow."""
    if stop_after and stop_after not in STAGE_ORDER:
        raise ValueError("stop_after must be one of: " + ", ".join(STAGE_ORDER))
    context = load_project(project)
    state = _State(context, restart=restart)
    checks = validate_project(context, require_host_reference=False)
    validation_path = context.output_root / "profile_validation.csv"
    checks.to_csv(validation_path, index=False)
    if checks["status"].eq("fail").any():
        state.set(
            "validate",
            "failed",
            outputs=[validation_path],
            message="Project validation failed; inspect profile_validation.csv",
        )
        raise ValueError("Project validation failed; inspect " + str(validation_path))
    state.set("validate", "completed", outputs=[validation_path])
    if stop_after == "validate":
        return project_status(context)

    workflow = context.values.get("workflow") or {}
    settings = _settings(context)
    alignment_path = context.profiles.resolve(context.profiles.virus["strain_alignment"])
    gff_path = context.profiles.resolve(context.profiles.virus["annotation_gff"])
    reference_id = str(context.profiles.virus["reference_accession"])
    assert alignment_path is not None and gff_path is not None
    discovery_dir = context.output_root / "discovery"
    discovery_dir.mkdir(parents=True, exist_ok=True)
    discovery_outputs = [
        discovery_dir / "candidates_ranked_pre_host.csv",
        discovery_dir / "candidates_rejected_pre_host.csv",
        discovery_dir / "candidate_feature_map.csv",
        discovery_dir / "discovery_panel.csv",
        discovery_dir / "selection_audit.csv",
    ]
    discovery_signature = _signature(
        [alignment_path, gff_path, context.ranking_config, context.profiles.source_paths[2]],
        {
            "minimum_strain_coverage": workflow.get("minimum_strain_coverage", 0.95),
            "balanced_top_per_gene": workflow.get("balanced_top_per_gene", 50),
            "global_top": workflow.get("global_top", 500),
        },
    )

    def discover() -> None:
        records = read_fasta(alignment_path)
        features = read_gff3(gff_path)
        raw = scan_editor_candidates(
            records,
            reference_id,
            context.profiles.editor,
            float(workflow.get("minimum_strain_coverage", 0.95)),
        )
        primary = annotate_candidates(raw, features, reference_id)
        ranked = rank_pre_human_candidates(primary, settings, evidence_path=None)
        retained = ranked[ranked["rejection_reasons"].fillna("").eq("")].copy()
        rejected = ranked[ranked["rejection_reasons"].fillna("").ne("")].copy()
        feature_map = build_candidate_feature_map(
            ranked, features, config=settings, annotation_source=gff_path
        )
        selection = select_balanced_discovery_panel(
            retained,
            feature_map,
            features,
            top_per_gene=int(workflow.get("balanced_top_per_gene", 50)),
            global_top=int(workflow.get("global_top", 500)),
        )
        retained.to_csv(discovery_outputs[0], index=False)
        rejected.to_csv(discovery_outputs[1], index=False)
        feature_map.to_csv(discovery_outputs[2], index=False)
        selection.panel.to_csv(discovery_outputs[3], index=False)
        selection.audit.to_csv(discovery_outputs[4], index=False)

    _run_cached_stage(state, "discover", discovery_signature, discovery_outputs, discover)
    if stop_after == "discover":
        return project_status(context)

    panel = pd.read_csv(discovery_outputs[3])
    host_dir = context.output_root / "host_screen"
    host_dir.mkdir(parents=True, exist_ok=True)
    host_input = host_dir / "cas_offinder_input.txt"
    host_manifest = host_dir / "selected_candidates.csv"
    host_output = host_dir / "cas_offinder_output.tsv"
    post_host = host_dir / "candidates_ranked_post_host.csv"
    predicted_hits_path = host_dir / "predicted_host_hits.csv"
    host_root = context.profiles.resolve(context.profiles.host.get("fasta_root"))
    host_outputs = [host_input, host_manifest]
    host_signature_paths = [discovery_outputs[3], context.profiles.source_paths[1]]
    if host_root and host_root.is_file():
        host_signature_paths.append(host_root)
    host_signature = _signature(
        host_signature_paths,
        {
            "assembly": context.profiles.host.get("assembly_accession", ""),
            "maximum_candidates": workflow.get("maximum_host_candidates"),
            "mismatch_threshold": context.profiles.editor.mismatch_search_threshold,
        },
    )
    if not _host_reference_available(host_root):
        state.set(
            "host_screen",
            "external_required",
            signature=host_signature,
            outputs=[],
            message="Host reference is missing; no host-screen result was inferred.",
        )
    else:
        if not state.reusable("host_screen_prepare", host_signature, host_outputs):
            maximum = workflow.get("maximum_host_candidates")
            maximum = len(panel) if maximum in (None, "") else int(maximum)
            build_cas_offinder_input(
                panel,
                host_root,
                host_input,
                host_manifest,
                maximum_candidates=maximum,
                stratify_by_gene=True,
                config=settings,
            )
            state.set(
                "host_screen_prepare",
                "completed",
                signature=host_signature,
                outputs=host_outputs,
            )
        executable = _cas_offinder_executable()
        if run_external and not host_output.is_file():
            if not executable:
                state.set(
                    "host_screen",
                    "external_required",
                    signature=host_signature,
                    outputs=host_outputs,
                    message=(
                        "Cas-OFFinder is unavailable. Set CAS_OFFINDER_PATH, run it against "
                        "cas_offinder_input.txt, and resume."
                    ),
                )
            else:
                state.set("host_screen", "running", signature=host_signature, outputs=[host_output])
                _run_cas_offinder(executable, host_input, host_output)
        if host_output.is_file():
            selected = pd.read_csv(host_manifest)
            hits = read_cas_offinder_output(host_output)
            summarized = summarize_cas_offinder_hits(
                panel,
                hits,
                max_mismatches=context.profiles.editor.mismatch_search_threshold,
                selected_manifest=host_manifest,
                config=settings,
            )
            predicted_hits = summarized.attrs.get("predicted_human_hits", pd.DataFrame())
            ranked_post = rank_post_human_candidates(summarized, settings)
            # Only selected candidates have a completed host screen. Preserve that boundary.
            selected_ids = set(selected["candidate_id"].astype(str))
            ranked_post = ranked_post[
                ranked_post["candidate_id"].astype(str).isin(selected_ids)
            ].copy()
            ranked_post.to_csv(post_host, index=False)
            predicted_hits.to_csv(predicted_hits_path, index=False)
            state.set(
                "host_screen",
                "completed",
                signature=host_signature,
                outputs=[host_input, host_manifest, host_output, post_host, predicted_hits_path],
                message=(
                    "Completed model-bounded host search; zero predicted hits is not proof "
                    "of safety."
                ),
            )
        elif not run_external or not executable:
            state.set(
                "host_screen",
                "external_required",
                signature=host_signature,
                outputs=host_outputs,
                message=(
                    "Input prepared. Run Cas-OFFinder to create cas_offinder_output.tsv, "
                    "then execute `vst project resume`."
                ),
            )
    if stop_after == "host_screen":
        return project_status(context)

    pair_dir = context.output_root / "pairs"
    pair_dir.mkdir(parents=True, exist_ok=True)
    pair_outputs = [
        pair_dir / "pair_hypotheses_same_gene.csv",
        pair_dir / "pair_hypotheses_multi_target.csv",
    ]
    pair_signature = _signature(
        [discovery_outputs[3], gff_path, alignment_path, context.ranking_config], {}
    )

    def pairs() -> None:
        features = read_gff3(gff_path)
        records = read_fasta(alignment_path)
        common = {
            "features": features,
            "aligned_records": records,
            "reference_id": reference_id,
            "config": settings,
        }
        same = rank_candidate_pairs(panel, same_feature_only=True, **common)
        all_pairs = rank_candidate_pairs(panel, same_feature_only=False, **common)
        if "hypothesis_type" in all_pairs:
            multi = all_pairs[all_pairs["hypothesis_type"].eq("multi_target_hypothesis")].copy()
        else:
            multi = all_pairs.copy()
        same.to_csv(pair_outputs[0], index=False)
        multi.to_csv(pair_outputs[1], index=False)

    _run_cached_stage(state, "pairs", pair_signature, pair_outputs, pairs)
    if stop_after == "pairs":
        return project_status(context)

    report_dir = context.output_root / "report"
    report_dir.mkdir(parents=True, exist_ok=True)
    report_outputs = [
        report_dir / "report.html",
        report_dir / "methods.md",
        report_dir / "limitations.md",
        context.output_root / "run_manifest.json",
    ]
    report_inputs = [
        post_host if post_host.is_file() else discovery_outputs[3],
        *pair_outputs,
        context.source,
        *context.profiles.source_paths,
    ]
    report_signature = _signature(report_inputs, {"project_id": context.values["id"]})

    def report() -> None:
        candidates = pd.read_csv(post_host if post_host.is_file() else discovery_outputs[3])
        same = _read_optional_csv(pair_outputs[0])
        multi = _read_optional_csv(pair_outputs[1])
        hits = _read_optional_csv(predicted_hits_path)
        write_methods_and_limitations(report_dir)
        write_html_report(
            candidates,
            report_outputs[0],
            title=f"{context.values.get('display_name', context.values['id'])} research report",
            rejected=pd.read_csv(discovery_outputs[1]),
            pairs=pd.concat([same, multi], ignore_index=True),
            predicted_hits=hits,
            output_links=[path.name for path in report_outputs[:3]],
        )
        manifest = {
            "schema_version": "1.0",
            "project_id": context.values["id"],
            "created_utc": _utc_now(),
            "project_file": str(context.source),
            "project_file_sha256": sha256_file(context.source),
            "profile_sources": [
                {"path": str(path), "sha256": sha256_file(path)}
                for path in context.profiles.source_paths
            ],
            "input_sources": [
                {"path": str(path), "sha256": sha256_file(path)}
                for path in (alignment_path, gff_path, context.ranking_config)
            ],
            "stage_status": {
                stage: record.get("status")
                for stage, record in state.values.get("stages", {}).items()
            },
            "interpretation": (
                "Computational target-prioritization hypotheses only; no wet-lab protocol, "
                "editing, safety, efficacy, delivery, or therapeutic claim."
            ),
        }
        report_outputs[3].write_text(
            json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
        )

    _run_cached_stage(state, "report", report_signature, report_outputs, report)
    return project_status(context)


def project_status(project: ProjectContext | str | Path) -> dict[str, Any]:
    context = load_project(project) if not isinstance(project, ProjectContext) else project
    state_path = context.output_root / "workflow_state.json"
    values = (
        json.loads(state_path.read_text(encoding="utf-8"))
        if state_path.is_file()
        else {"stages": {}}
    )
    rows = []
    for stage in STAGE_ORDER:
        record = values.get("stages", {}).get(stage, {})
        rows.append(
            {
                "stage": stage,
                "status": record.get("status", "pending"),
                "message": record.get("message", ""),
            }
        )
    return {
        "project_id": context.values["id"],
        "project_file": str(context.source),
        "output_root": str(context.output_root),
        "state_file": str(state_path),
        "stages": rows,
        "report": str(context.output_root / "report" / "report.html"),
    }
