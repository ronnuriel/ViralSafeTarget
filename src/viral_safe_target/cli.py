"""Command-line interface for the research pipeline."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import webbrowser
from pathlib import Path

import pandas as pd

from . import __version__
from .annotations import annotate_candidates, read_gff3
from .benchmarking import run_benchmark, write_benchmark_outputs
from .config import DEFAULT_CONFIG_PATH, get_editor, load_config
from .crispr import scan_editor_candidates
from .disruption import rank_candidate_pairs
from .experimental import import_crispresso2_results
from .integrations import (
    CasOffinderAdapter,
    CrispritzAdapter,
    ExternalImportAdapter,
    MafftAdapter,
    load_external_results,
)
from .io_utils import read_fasta
from .offtarget import (
    build_cas_offinder_input,
    read_cas_offinder_output,
    screen_against_small_fasta,
    summarize_cas_offinder_hits,
    write_offtarget_metadata,
)
from .project_workflow import STAGE_ORDER
from .provenance import sha256_file, write_run_manifest
from .reporting import write_html_report, write_methods_and_limitations
from .scoring import rank_post_human_candidates, rank_pre_human_candidates
from .sdk import load_run


def _config(args: argparse.Namespace) -> dict:
    return load_config(getattr(args, "config", None))


def _split_and_write(ranked: pd.DataFrame, out_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    rejected = ranked[ranked["rejection_reasons"].fillna("").ne("")].copy()
    retained = ranked[ranked["rejection_reasons"].fillna("").eq("")].copy()
    retained.to_csv(out_dir / "candidates_ranked_pre_human.csv", index=False)
    rejected.to_csv(out_dir / "candidates_rejected_pre_human.csv", index=False)
    return retained, rejected


def _scan(args: argparse.Namespace) -> None:
    settings = _config(args)
    editor = get_editor(settings)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    records = read_fasta(args.virus_alignment)
    candidates = scan_editor_candidates(records, args.reference_id, editor, args.min_coverage)
    if args.gff:
        candidates = annotate_candidates(candidates, read_gff3(args.gff), args.reference_id)
    ranked = rank_pre_human_candidates(candidates, settings, args.gene_evidence)
    retained, rejected = _split_and_write(ranked, out_dir)
    result = retained
    predicted_hits = pd.DataFrame()
    if args.small_host_fasta:
        result = screen_against_small_fasta(
            retained, read_fasta(args.small_host_fasta), args.max_mismatches
        )
        result = rank_post_human_candidates(result, settings)
        result.to_csv(out_dir / "candidates_ranked_post_human.csv", index=False)
    result.to_csv(out_dir / "candidates.csv", index=False)
    links = [
        "candidates_ranked_pre_human.csv",
        "candidates_rejected_pre_human.csv",
        "candidates.csv",
        "run_manifest.json",
        "methods.md",
        "limitations.md",
    ]
    write_methods_and_limitations(out_dir)
    write_html_report(
        result,
        out_dir / "report.html",
        title="ViralSafeTarget scan",
        rejected=rejected,
        predicted_hits=predicted_hits,
        output_links=links,
    )
    inputs = [args.virus_alignment, *([args.gff] if args.gff else [])]
    if args.small_host_fasta:
        inputs.append(args.small_host_fasta)
    write_run_manifest(
        out_dir / "run_manifest.json",
        inputs,
        {"reference_id": args.reference_id, "minimum_coverage": args.min_coverage},
        config_path=settings["_config_path"],
        editor_profile=settings["editor"],
        accepted_accessions=records.keys(),
        human_assembly_identifier=(
            "synthetic_or_small_host_fasta" if args.small_host_fasta else None
        ),
        random_seed=int(settings["random_seed"]),
        output_paths=[out_dir / link for link in links if (out_dir / link).exists()],
    )
    print(f"Wrote {len(retained)} retained and {len(rejected)} rejected candidates to {out_dir}")


def _rank(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ranked = rank_pre_human_candidates(
        pd.read_csv(args.candidates), _config(args), args.gene_evidence
    )
    retained, rejected = _split_and_write(ranked, out_dir)
    print(f"Wrote {len(retained)} retained and {len(rejected)} rejected candidates")


def _build_offtarget(args: argparse.Namespace) -> None:
    settings = _config(args)
    maximum = args.max_candidates or int(settings["off_target"]["maximum_candidates"])
    selected = build_cas_offinder_input(
        pd.read_csv(args.candidates),
        args.human_fasta_directory,
        args.output,
        args.manifest,
        maximum_candidates=maximum,
        genes=args.genes,
        stratify_by_gene=not args.no_stratify,
        config=settings,
    )
    metadata = Path(args.manifest).with_suffix(".metadata.json")
    write_offtarget_metadata(metadata, settings)
    print(f"Wrote {len(selected)} queries to {args.output}; manifest: {args.manifest}")


def _summarize_offtargets(args: argparse.Namespace) -> None:
    settings = _config(args)
    candidates = pd.read_csv(args.candidates)
    hits = read_cas_offinder_output(args.cas_output)
    summarized = summarize_cas_offinder_hits(
        candidates,
        hits,
        max_mismatches=int(settings["editor"]["mismatch_search_threshold"]),
        selected_manifest=args.manifest,
        human_gff=args.human_gff,
        config=settings,
    )
    predicted_hits = summarized.attrs.get("predicted_human_hits", pd.DataFrame()).copy()
    ranked = rank_post_human_candidates(summarized, settings)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(out_dir / "candidates_ranked_post_human.csv", index=False)
    predicted_hits.to_csv(out_dir / "predicted_human_hits.csv", index=False)
    print(f"Wrote post-human ranking for {len(ranked)} candidates to {out_dir}")


def _simulate_pairs(args: argparse.Namespace) -> None:
    candidates = pd.read_csv(args.candidates)
    features = read_gff3(args.gff) if args.gff else None
    alignment = read_fasta(args.virus_alignment) if args.virus_alignment else None
    common = dict(
        features=features,
        aligned_records=alignment,
        reference_id=args.reference_id,
        min_distance_bp=args.min_distance,
        max_distance_bp=args.max_distance,
        max_candidates=args.max_candidates,
        config=_config(args),
        genes=args.genes,
        feature_types=args.feature_types,
        minimum_conservation=args.minimum_conservation,
        maximum_viral_occurrence_count=args.maximum_viral_occurrence_count,
        maximum_candidates_per_gene=args.maximum_candidates_per_gene,
        stratify_by_gene=not args.no_stratify,
    )
    if args.out_dir:
        out_dir = Path(args.out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        same = rank_candidate_pairs(candidates, same_feature_only=True, **common)
        all_pairs = rank_candidate_pairs(candidates, same_feature_only=False, **common)
        multi = all_pairs[all_pairs["hypothesis_type"] == "multi_target_hypothesis"].copy()
        same.to_csv(out_dir / "pair_hypotheses_same_gene.csv", index=False)
        multi.to_csv(out_dir / "pair_hypotheses_multi_target.csv", index=False)
        print(f"Wrote {len(same)} same-gene and {len(multi)} multi-target hypotheses")
    else:
        pairs = rank_candidate_pairs(
            candidates, same_feature_only=not args.allow_cross_feature, **common
        )
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        pairs.to_csv(output, index=False)
        print(f"Wrote {len(pairs)} pair hypotheses to {output}")


def _benchmark(args: argparse.Namespace) -> None:
    detail, summary = run_benchmark(pd.read_csv(args.candidates), pd.read_csv(args.known_targets))
    write_benchmark_outputs(detail, summary, args.out_dir)
    print(json.dumps(summary, indent=2, sort_keys=True))


def _report(args: argparse.Namespace) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(args.candidates)
    rejected = pd.read_csv(args.rejected) if args.rejected else pd.DataFrame()
    pairs = pd.read_csv(args.pairs) if args.pairs else pd.DataFrame()
    if args.multi_pairs:
        pairs = pd.concat([pairs, pd.read_csv(args.multi_pairs)], ignore_index=True)
    hits = pd.read_csv(args.predicted_hits) if args.predicted_hits else pd.DataFrame()
    write_methods_and_limitations(out_dir)
    links = [path.name for path in out_dir.iterdir() if path.is_file()]
    write_html_report(
        candidates,
        out_dir / "report.html",
        title=args.title,
        rejected=rejected,
        pairs=pairs,
        predicted_hits=hits,
        output_links=sorted(set(links + ["report.html"])),
    )
    print(f"Wrote {out_dir / 'report.html'}")


def _tool_version(command: str, arguments: list[str]) -> str:
    executable = shutil.which(command)
    if not executable:
        return "not found"
    try:
        completed = subprocess.run(
            [executable, *arguments], capture_output=True, text=True, timeout=10, check=False
        )
        text = (completed.stdout or completed.stderr).strip().splitlines()
        return f"{executable} — {text[0] if text else 'available'}"
    except (OSError, subprocess.TimeoutExpired) as error:
        return f"{executable} — version check failed: {error}"


def _doctor(args: argparse.Namespace) -> None:
    settings = _config(args)
    root = Path.cwd()
    disk = shutil.disk_usage(root)
    checks = {
        "Python": f"{sys.executable} — {sys.version.split()[0]}",
        "package import": "available",
        "MAFFT": _tool_version("mafft", ["--version"]),
        "NCBI datasets": _tool_version("datasets", ["version"]),
        "Cas-OFFinder": _tool_version("cas-offinder", ["--help"]),
        "OpenCL": _tool_version("clinfo", ["--version"]),
        "input directory": str((root / "data").resolve()),
        "output directory": str((root / "reports").resolve()),
        "disk free": f"{disk.free / (1024**3):.2f} GiB",
        "active configuration": settings["_config_path"],
        "editor": settings["editor"]["name"],
    }
    print("ViralSafeTarget doctor")
    for label, value in checks.items():
        print(f"- {label}: {value}")


def _tools_doctor(args: argparse.Namespace) -> None:
    del args
    adapters = [MafftAdapter(), CasOffinderAdapter(), CrispritzAdapter()]
    adapters.extend(ExternalImportAdapter(name) for name in ("crispor", "chopchop", "guidescan2"))
    print("ViralSafeTarget external-tool adapters")
    for adapter in adapters:
        status = adapter.detect()
        label = "available" if status.available else "pending"
        details = status.version or status.message
        print(f"- {status.name}: {label} — {details}")


def _crispritz_build(args: argparse.Namespace) -> None:
    run = load_run(args.run_dir)
    candidates = run.candidates
    if "human_total_predicted_hits" in candidates:
        candidates = candidates[candidates["human_total_predicted_hits"].eq(0)].copy()
    manifest = CrispritzAdapter(args.docker_image).build_input(
        candidates,
        {
            "reference_genome": args.reference_genome,
            "genome_or_assembly": args.assembly,
            "pam_file": args.pam_file,
            "mismatches": args.mismatches,
            "dna_bulges": args.dna_bulges,
            "rna_bulges": args.rna_bulges,
            "variant_aware": args.variant_aware,
            "variant_file": args.variant_file,
        },
        args.out_dir,
    )
    print(f"Wrote CRISPRitz input bundle for {len(candidates)} candidates: {manifest}")


def _crispritz_run(args: argparse.Namespace) -> None:
    manifest = Path(args.input_dir) / "crispritz_manifest.json"
    execution = CrispritzAdapter(args.docker_image).run(
        manifest, args.out_dir, dry_run=args.dry_run
    )
    if args.dry_run:
        print("Dry-run command:", " ".join(execution.command))
    else:
        print(f"CRISPRitz completed with return code {execution.returncode}")


def _crispritz_import(args: argparse.Namespace) -> None:
    adapter = CrispritzAdapter(args.docker_image)
    manifest_path = Path(args.manifest)
    settings = json.loads(manifest_path.read_text(encoding="utf-8"))
    guides = manifest_path.parent / settings["guides_file"]
    candidates = pd.read_csv(guides, sep="\t")
    parsed = adapter.parse(args.results, guides)
    normalized = adapter.normalize(
        parsed,
        candidates=candidates,
        source_file=args.results,
        assembly=settings.get("genome_or_assembly", ""),
        command="vst tools crispritz import",
    )
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    parsed.to_csv(output / "crispritz_parsed_raw.csv", index=False)
    normalized.to_csv(output / "tool_results_long.csv", index=False)
    print(f"Imported {len(parsed)} CRISPRitz hits for {len(candidates)} candidates")


def _external_import(args: argparse.Namespace) -> None:
    candidates = pd.read_csv(args.candidates)
    imported = load_external_results(args.tool, args.input, args.mapping, candidates)
    output = Path(args.out_dir)
    output.mkdir(parents=True, exist_ok=True)
    imported.results.dataframe.to_csv(output / "tool_results_long.csv", index=False)
    imported.raw_rows.to_csv(output / "raw_external_rows.csv", index=False)
    unmatched = pd.concat([imported.unmatched_rows, imported.ambiguous_rows], ignore_index=True)
    unmatched.to_csv(output / "unmatched_external_rows.csv", index=False)
    print(
        f"Imported {len(imported.results)} normalized metric rows; "
        f"{len(imported.unmatched_rows)} unmatched and {len(imported.ambiguous_rows)} ambiguous"
    )


def _crispresso_import(args: argparse.Namespace) -> None:
    imported = import_crispresso2_results(args.input, args.candidate_map)
    measurements, manifest = imported.write(args.out_dir)
    print(f"Wrote measured metrics to {measurements}; provenance: {manifest}")


def _discover_genome_wide(args: argparse.Namespace) -> None:
    from .discovery_workflow import run_genome_wide_discovery

    result = run_genome_wide_discovery(
        virus=args.virus,
        run_dir=args.run_dir,
        config_path=args.config,
        out_dir=args.out_dir,
        top_per_gene=args.top_per_gene,
        global_top=args.global_top,
        batch_size=args.batch_size,
        run_cas_offinder=args.run_cas_offinder,
        analysis_only=args.analysis_only,
        exhaustive=args.exhaustive,
        confirm_exhaustive=args.confirm_exhaustive,
    )
    print(json.dumps(result["answers"], indent=2, sort_keys=True))
    if args.open_report:
        webbrowser.open((result["output_dir"] / "report.html").as_uri())


def _analyze_gene_function(args: argparse.Namespace) -> None:
    from .gene_function_workflow import run_gene_function_analysis

    result = run_gene_function_analysis(
        genome_wide_dir=args.genome_wide_dir,
        hsv2_genbank=args.hsv2_genbank,
        hsv1_genbank=args.hsv1_genbank,
        virus_alignment=args.virus_alignment,
        domain_table=args.domain_table,
        disorder_table=args.disorder_table,
        evidence_table=args.evidence_table,
        out_dir=args.out_dir,
        top_per_gene=args.top_per_gene,
    )
    print(json.dumps({key: str(value) for key, value in result.items()}, indent=2))
    if args.open_report:
        webbrowser.open((result["output_dir"] / "gene_evidence_report.html").as_uri())


def _analyze_population(args: argparse.Namespace) -> None:
    from .population_reporting import build_population_comparison, write_population_report
    from .population_validation import (
        exact_guide_presence_by_accession,
        map_population_to_reference,
        summarize_locus_aware_population_validation,
    )

    candidates = pd.read_csv(args.candidates)
    records = read_fasta(args.population_fasta)
    editor = get_editor(load_config(args.config))
    presence = exact_guide_presence_by_accession(candidates, records, editor)
    alignments = map_population_to_reference(
        records,
        args.reference_fasta,
        minimum_mapq=args.minimum_mapq,
        minimum_identity=args.minimum_identity,
    )
    validation = summarize_locus_aware_population_validation(
        candidates, presence, alignments, list(records)
    )
    comparison, genes = build_population_comparison(candidates, validation)
    output = Path(args.out_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    alignments.to_csv(output / "population_reference_alignments.csv", index=False)
    validation.to_csv(output / "candidate_locus_population_validation.csv", index=False)
    provenance = {
        "schema_version": "1.0",
        "population_record_count": len(records),
        "candidate_count": len(candidates),
        "reference_alignment_count": len(alignments),
        "records_with_reference_alignment": int(alignments["accession"].nunique())
        if not alignments.empty
        else 0,
        "minimum_mapq": args.minimum_mapq,
        "minimum_identity": args.minimum_identity,
        "population_fasta": str(Path(args.population_fasta).resolve()),
        "population_fasta_sha256": sha256_file(args.population_fasta),
        "reference_fasta": str(Path(args.reference_fasta).resolve()),
        "reference_fasta_sha256": sha256_file(args.reference_fasta),
        "candidate_source": str(Path(args.candidates).resolve()),
        "candidate_source_sha256": sha256_file(args.candidates),
        "score_integration": "none; population evidence is reported separately",
        "interpretation": (
            "Held-out population-genomic validation only; no editing efficacy, safety, "
            "delivery, or therapeutic claim."
        ),
    }
    write_population_report(comparison, genes, output, provenance)
    print(json.dumps(provenance, indent=2, sort_keys=True))
    if args.open_report:
        webbrowser.open((output / "population_validation_report.html").as_uri())


def _analyze_virtual_knockout(args: argparse.Namespace) -> None:
    from .virtual_analysis_workflow import run_virtual_knockout_analysis

    result = run_virtual_knockout_analysis(args.project)
    print(json.dumps({key: str(value) for key, value in result.items()}, indent=2))


def _analyze_escape(args: argparse.Namespace) -> None:
    from .virtual_analysis_workflow import run_escape_analysis

    result = run_escape_analysis(args.project)
    print(json.dumps({key: str(value) for key, value in result.items()}, indent=2))


def _analyze_multiplex(args: argparse.Namespace) -> None:
    from .virtual_analysis_workflow import run_full_virtual_analysis

    result = run_full_virtual_analysis(args.project)
    print(json.dumps(result, indent=2, default=str))
    if args.open_report:
        webbrowser.open(Path(result["report"]).as_uri())


def _tool_benchmark(args: argparse.Namespace) -> None:
    from .tool_benchmark import run_tool_benchmark

    result = run_tool_benchmark(args.config)
    print(json.dumps(result, indent=2, default=str))
    if args.open_report:
        webbrowser.open(Path(result["report"]).resolve().as_uri())


def _profiles_validate(args: argparse.Namespace) -> None:
    from .profiles import load_profile_bundle, validate_profile_bundle

    bundle = load_profile_bundle(
        args.virus_profile,
        args.host_profile,
        args.nuclease_profile,
        project_root=args.project_root,
    )
    checks = validate_profile_bundle(
        bundle,
        require_large_host_reference=args.require_host_reference,
        require_virus_inputs=args.require_virus_inputs,
    )
    print(checks.to_string(index=False))
    if checks["status"].eq("fail").any():
        raise SystemExit(2)


def _showcase_build(args: argparse.Namespace) -> None:
    from .showcase_workflow import run_showcase

    result = run_showcase(
        virus_profile=args.virus_profile,
        host_profile=args.host_profile,
        nuclease_profile=args.nuclease_profile,
        out_dir=args.out_dir,
        project_root=args.project_root,
        per_gene=args.per_gene,
    )
    print(json.dumps({key: str(value) for key, value in result.items()}, indent=2))
    if args.open_report:
        webbrowser.open((result["output_dir"] / "FINAL_REPORT.html").as_uri())


def _project_init(args: argparse.Namespace) -> None:
    from .project_workflow import initialize_project

    project_file = initialize_project(
        args.out_dir,
        project_id=args.id,
        display_name=args.display_name or args.id,
        reference_accession=args.reference_accession,
        force=args.force,
    )
    print(f"Created research project: {project_file}")
    print(f"Next: add inputs, then run `vst project validate --project {project_file}`")


def _project_validate(args: argparse.Namespace) -> None:
    from .project_workflow import validate_project

    checks = validate_project(args.project, require_host_reference=args.require_host_reference)
    print(checks.to_string(index=False))
    if checks["status"].eq("fail").any():
        raise SystemExit(2)


def _project_run(args: argparse.Namespace) -> None:
    from .project_workflow import run_project

    result = run_project(
        args.project,
        run_external=args.run_external,
        restart=args.restart,
        stop_after=args.stop_after,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    if args.open_report and Path(result["report"]).is_file():
        webbrowser.open(Path(result["report"]).resolve().as_uri())


def _project_status(args: argparse.Namespace) -> None:
    from .project_workflow import project_status

    print(json.dumps(project_status(args.project), indent=2, sort_keys=True))


def _prompt(label: str, default: str = "") -> str:
    suffix = f" [{default}]" if default else ""
    value = input(f"{label}{suffix}: ").strip()
    return value or default


def _simple_init(args: argparse.Namespace) -> None:
    from .researcher import create_project

    project_name = args.project_name
    virus_name = args.virus_name
    tax_id = args.tax_id
    reference_accession = args.reference_accession
    reference_fasta = args.reference_fasta
    annotation_gff = args.annotation_gff
    strains_fasta = args.strains_fasta
    host_fasta = args.host_fasta
    sequence_only = args.sequence_only
    if args.interactive:
        virus_name = virus_name or _prompt("Virus scientific name", project_name)
        tax_id = tax_id or _prompt("NCBI tax ID (optional)") or None
        if not reference_accession and not reference_fasta:
            source = _prompt("Reference accession, or enter 'local'", "local")
            if source.lower() == "local":
                reference_fasta = _prompt("Local reference FASTA")
            else:
                reference_accession = source
        if reference_fasta and not annotation_gff:
            annotation_gff = _prompt("Local GFF3 (blank for sequence-only)") or None
            sequence_only = not annotation_gff
        if not strains_fasta:
            strains_fasta = _prompt("Strain FASTA (optional; reference-only if blank)") or None
        if not host_fasta:
            host_fasta = _prompt("Local host FASTA (optional)") or None
    virus_name = virus_name or project_name
    output = Path(args.output_dir or project_name)
    project_file = create_project(
        output,
        project_name=project_name,
        virus_name=virus_name,
        tax_id=tax_id,
        reference_accession=reference_accession,
        reference_fasta=reference_fasta,
        annotation_gff=annotation_gff,
        strains_fasta=strains_fasta,
        strains_aligned=not args.strains_unaligned,
        host_profile=args.host_profile,
        host_fasta=host_fasta,
        nuclease_profile=args.nuclease_profile,
        evidence_enabled=args.enable_evidence,
        run_external=args.run_external,
        sequence_only=sequence_only,
        force=args.force,
    )
    print(f"Created research project: {project_file}")
    print(f"Next: vst plan {project_file}")


def _quickstart(args: argparse.Namespace) -> None:
    from .researcher import create_demo_project

    project = create_demo_project(args.out, force=args.force)
    print(f"Created synthetic demo project: {project}")
    print(f"Next: vst run {project}")


def _plan(args: argparse.Namespace) -> None:
    from .researcher import plan_project

    result = plan_project(args.project)
    if args.json:
        print(json.dumps(result, indent=2, sort_keys=True))
        return
    print(f"Project: {result['project']}")
    print(f"Candidate-count estimate: {result['candidate_count_estimate']}")
    for row in result["stages"]:
        estimate = row["estimated_seconds"]
        rendered = "unavailable" if estimate is None else f"{estimate:.2f} s"
        cached = " (cached)" if row["cached"] else ""
        print(f"- {row['stage']}: {rendered}; {row['confidence']}{cached} — {row['reason']}")
    total = result["total_estimated_seconds"]
    print(f"Total estimate: {'unavailable' if total is None else f'{total:.2f} s'}")
    print("Runtime estimates are best effort, hardware-dependent, and not guarantees.")


def _simple_run(args: argparse.Namespace) -> None:
    from .project_workflow import run_project

    result = run_project(
        args.project,
        run_external=args.run_external,
        restart=args.restart,
        stop_after=args.stop_after,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


def _simple_status(args: argparse.Namespace) -> None:
    from .project_workflow import project_status

    print(json.dumps(project_status(args.project), indent=2, sort_keys=True))


def _open_results(args: argparse.Namespace) -> None:
    from .researcher import open_results

    print(open_results(args.path, no_browser=args.no_browser))


def _export_results(args: argparse.Namespace) -> None:
    from .researcher import export_project

    print(
        export_project(
            args.project,
            output=args.output,
            include_large_raw=args.include_large_raw,
        )
    )


def _researcher_doctor(args: argparse.Namespace) -> None:
    from .researcher import doctor_report

    report = doctor_report(args.project)
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return
    print(f"ViralSafeTarget {report['viral_safe_target_version']}")
    print(f"Python: {report['python']['version']} ({report['python']['executable']})")
    print(f"Memory: {report['memory_gib']} GiB; disk free: {report['disk_free_gib']} GiB")
    for name, value in report["tools"].items():
        print(f"- {name}: {value['version']}")
    if "project" in report:
        print(f"Project sequence stages: {report['project']['can_run_sequence_stages']}")
        print(f"Project external host stage: {report['project']['can_run_external_host_stage']}")


def _tools_setup(args: argparse.Namespace) -> None:
    del args
    print("External tools are explicit optional dependencies.")
    print("MAFFT: https://mafft.cbrc.jp/alignment/software/")
    print("Cas-OFFinder: https://github.com/snugel/cas-offinder")
    print("CRISPRitz: https://github.com/pinellolab/CRISPRitz")
    print("After installation, run: vst tools status")


def _reproduce_hsv2(args: argparse.Namespace) -> None:
    from .reproduction import reproduce_hsv2

    result = reproduce_hsv2(
        args.project_root,
        execute=args.execute,
        skip_population=args.skip_population,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    if args.open_report and result.get("final_report"):
        report = Path(result["final_report"])
        if report.is_file():
            webbrowser.open(report.resolve().as_uri())


def _evidence_discover(args: argparse.Namespace) -> None:
    from .evidence_agent import discover_evidence
    from .project_workflow import load_project

    context = load_project(args.project)
    gff_path = context.profiles.resolve(context.profiles.virus.get("annotation_gff"))
    if gff_path is None or not gff_path.is_file():
        raise FileNotFoundError("The project virus annotation GFF is missing")
    out_dir = Path(args.out_dir).resolve() if args.out_dir else context.output_root / "evidence"
    result = discover_evidence(
        gff_path=gff_path,
        virus_profile=context.profiles.virus,
        out_dir=out_dir,
        sources=args.sources,
        maximum_results_per_query=args.max_results_per_query,
        email=args.email,
        api_key=args.api_key,
        offline=args.offline,
        genes=args.genes,
    )
    print(json.dumps({key: str(value) for key, value in result.items()}, indent=2))
    if args.open_report:
        webbrowser.open((Path(result["output_dir"]) / "evidence_review_report.html").as_uri())


def _evidence_apply(args: argparse.Namespace) -> None:
    from .evidence_agent import apply_reviewed_evidence
    from .project_workflow import load_project

    context = load_project(args.project)
    review_queue = (
        Path(args.review_queue).resolve()
        if args.review_queue
        else context.output_root / "evidence" / "review_queue.tsv"
    )
    configured = context.profiles.resolve(context.profiles.virus.get("evidence_table"))
    output_table = Path(args.output).resolve() if args.output else configured
    if output_table is None:
        raise ValueError("No evidence output is configured; pass --output")
    result = apply_reviewed_evidence(
        review_queue,
        output_table,
        append=args.append,
        reference_accession=str(context.profiles.virus.get("reference_accession", "")),
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    print("Only approved rows were exported. Re-run the project to refresh evidence-aware outputs.")


def _add_config(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="vst",
        description="Virus-first computational target discovery for research use",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    simple_init = subparsers.add_parser("init", help="create a new-virus research project")
    simple_init.add_argument("project_name")
    simple_init.add_argument("--interactive", action="store_true")
    simple_init.add_argument("--virus-name")
    simple_init.add_argument("--tax-id")
    simple_init.add_argument("--reference-accession")
    simple_init.add_argument("--reference-fasta")
    simple_init.add_argument("--annotation-gff")
    simple_init.add_argument("--strains-fasta")
    simple_init.add_argument("--strains-unaligned", action="store_true")
    simple_init.add_argument("--host-profile", default="human_grch38")
    simple_init.add_argument("--host-fasta")
    simple_init.add_argument("--nuclease-profile", default="spcas9")
    simple_init.add_argument("--output-dir")
    simple_init.add_argument("--enable-evidence", action="store_true")
    simple_init.add_argument("--run-external", action="store_true")
    simple_init.add_argument("--sequence-only", action="store_true")
    simple_init.add_argument("--force", action="store_true")
    simple_init.set_defaults(func=_simple_init)

    quickstart = subparsers.add_parser("quickstart", help="create the bundled synthetic demo")
    quickstart.add_argument("--out", required=True)
    quickstart.add_argument("--force", action="store_true")
    quickstart.set_defaults(func=_quickstart)

    plan = subparsers.add_parser("plan", help="validate inputs and estimate a project run")
    plan.add_argument("project")
    plan.add_argument("--json", action="store_true")
    plan.set_defaults(func=_plan)

    for command, handler, help_text in (
        ("run", _simple_run, "run a project with resumable stage caching"),
        ("resume", _simple_run, "resume a project from valid cached stages"),
    ):
        simple_run = subparsers.add_parser(command, help=help_text)
        simple_run.add_argument("project")
        simple_run.add_argument("--run-external", action="store_true")
        simple_run.add_argument("--restart", action="store_true")
        simple_run.add_argument("--stop-after", choices=[*STAGE_ORDER, "bundle"])
        simple_run.set_defaults(func=handler)

    simple_status = subparsers.add_parser("status", help="show project stage status")
    simple_status.add_argument("project")
    simple_status.set_defaults(func=_simple_status)

    open_command = subparsers.add_parser("open", help="open a project's START_HERE report")
    open_command.add_argument("path")
    open_command.add_argument("--no-browser", action="store_true")
    open_command.set_defaults(func=_open_results)

    export_command = subparsers.add_parser("export", help="create a portable result ZIP")
    export_command.add_argument("project")
    export_command.add_argument("--output")
    export_command.add_argument("--include-large-raw", action="store_true")
    export_command.set_defaults(func=_export_results)

    project = subparsers.add_parser(
        "project", help="single-entry workflows for reproducible virus projects"
    )
    project_commands = project.add_subparsers(dest="project_command", required=True)
    project_init = project_commands.add_parser(
        "init", help="create a self-contained new-virus research project"
    )
    project_init.add_argument("--id", required=True)
    project_init.add_argument("--display-name")
    project_init.add_argument("--reference-accession", default="CHANGE_ME")
    project_init.add_argument("--out-dir", required=True)
    project_init.add_argument("--force", action="store_true")
    project_init.set_defaults(func=_project_init)

    project_validate = project_commands.add_parser(
        "validate", help="validate project inputs, profiles, and coordinate identity"
    )
    project_validate.add_argument("--project", required=True)
    project_validate.add_argument("--require-host-reference", action="store_true")
    project_validate.set_defaults(func=_project_validate)

    for command, help_text in (
        ("run", "run the project and cache each completed stage"),
        ("resume", "resume a project without repeating valid completed stages"),
    ):
        project_run = project_commands.add_parser(command, help=help_text)
        project_run.add_argument("--project", required=True)
        project_run.add_argument(
            "--run-external",
            action="store_true",
            help="run Cas-OFFinder when available; otherwise preserve external_required status",
        )
        project_run.add_argument(
            "--restart",
            action="store_true",
            help="discard workflow state and recompute stages",
        )
        project_run.add_argument(
            "--stop-after",
            choices=[
                "validate",
                "discover",
                "host_screen",
                "pairs",
                "virtual_knockout",
                "escape",
                "multiplex",
                "report",
                "bundle",
            ],
        )
        project_run.add_argument("--open-report", action="store_true")
        project_run.set_defaults(func=_project_run)

    project_status = project_commands.add_parser(
        "status", help="show completed, pending, failed, and external stages"
    )
    project_status.add_argument("--project", required=True)
    project_status.set_defaults(func=_project_status)

    reproduce = subparsers.add_parser(
        "reproduce", help="plan or execute a versioned computational case-study reproduction"
    )
    reproduce_commands = reproduce.add_subparsers(dest="reproduce_command", required=True)
    reproduce_hsv2 = reproduce_commands.add_parser(
        "hsv2", help="plan or run the complete public-data HSV-2 case study"
    )
    reproduce_hsv2.add_argument("--project-root", default=".")
    reproduce_hsv2.add_argument(
        "--execute",
        action="store_true",
        help="perform downloads and long external computations; without this flag print a plan",
    )
    reproduce_hsv2.add_argument("--skip-population", action="store_true")
    reproduce_hsv2.add_argument("--open-report", action="store_true")
    reproduce_hsv2.set_defaults(func=_reproduce_hsv2)

    evidence = subparsers.add_parser(
        "evidence", help="discover source-linked gene evidence with mandatory human review"
    )
    evidence_commands = evidence.add_subparsers(dest="evidence_command", required=True)
    evidence_discover = evidence_commands.add_parser(
        "discover", help="query official sources and create a review-pending evidence queue"
    )
    evidence_discover.add_argument("--project", required=True)
    evidence_discover.add_argument("--out-dir")
    evidence_discover.add_argument(
        "--sources",
        nargs="+",
        choices=["pubmed", "europepmc", "uniprot", "ncbi_refseq"],
        default=["pubmed", "europepmc", "uniprot", "ncbi_refseq"],
    )
    evidence_discover.add_argument("--max-results-per-query", type=int, default=5)
    evidence_discover.add_argument(
        "--genes", nargs="*", help="optional annotated gene subset; default: all genes"
    )
    evidence_discover.add_argument("--email", default="")
    evidence_discover.add_argument("--api-key", default="")
    evidence_discover.add_argument(
        "--offline", action="store_true", help="use only previously cached API responses"
    )
    evidence_discover.add_argument("--open-report", action="store_true")
    evidence_discover.set_defaults(func=_evidence_discover)

    evidence_apply = evidence_commands.add_parser(
        "apply", help="export only researcher-approved proposals to gene_evidence.tsv"
    )
    evidence_apply.add_argument("--project", required=True)
    evidence_apply.add_argument("--review-queue")
    evidence_apply.add_argument("--output")
    evidence_apply.add_argument("--append", action="store_true")
    evidence_apply.set_defaults(func=_evidence_apply)

    discover = subparsers.add_parser(
        "discover", help="advanced case-study discovery commands; prefer `vst project`"
    )
    discover_commands = discover.add_subparsers(dest="discover_command", required=True)
    genome_wide = discover_commands.add_parser(
        "genome-wide", help="run or resume balanced HSV-2 genome-wide discovery"
    )
    genome_wide.add_argument("--virus", default="hsv2")
    genome_wide.add_argument("--run-dir")
    genome_wide.add_argument("--config", default="configs/hsv2_genome_wide.yaml")
    genome_wide.add_argument("--out-dir", default="reports/hsv2_genome_wide")
    genome_wide.add_argument("--top-per-gene", type=int)
    genome_wide.add_argument("--global-top", type=int)
    genome_wide.add_argument("--batch-size", type=int)
    genome_wide.add_argument("--run-cas-offinder", action="store_true")
    genome_wide.add_argument("--analysis-only", action="store_true")
    genome_wide.add_argument("--exhaustive", action="store_true")
    genome_wide.add_argument("--confirm-exhaustive", action="store_true")
    genome_wide.add_argument("--open-report", action="store_true")
    genome_wide.set_defaults(func=_discover_genome_wide)

    analyze = subparsers.add_parser("analyze", help="post-discovery biological analyses")
    analyze_commands = analyze.add_subparsers(dest="analyze_command", required=True)
    gene_function = analyze_commands.add_parser(
        "gene-function",
        help="HSV-2 case-study protein adapter; generic projects do not infer this stage",
    )
    gene_function.add_argument("--genome-wide-dir", default="reports/hsv2_genome_wide")
    gene_function.add_argument("--hsv2-genbank", default="data/raw/hsv2_reference.gb")
    gene_function.add_argument("--hsv1-genbank", default="data/raw/hsv1_reference.gb")
    gene_function.add_argument("--virus-alignment", default="data/processed/hsv2_aligned_25.fasta")
    gene_function.add_argument("--domain-table", default="data/curated/hsv2_target_domains.tsv")
    gene_function.add_argument("--disorder-table", default="data/curated/hsv2_target_disorder.tsv")
    gene_function.add_argument(
        "--evidence-table", default="data/curated/hsv_gene_function_evidence.tsv"
    )
    gene_function.add_argument("--out-dir", default="reports/hsv2_gene_function")
    gene_function.add_argument("--top-per-gene", type=int, default=10)
    gene_function.add_argument("--open-report", action="store_true")
    gene_function.set_defaults(func=_analyze_gene_function)

    population = analyze_commands.add_parser(
        "population", help="run locus-aware held-out viral population validation"
    )
    population.add_argument("--population-fasta", required=True)
    population.add_argument("--reference-fasta", required=True)
    population.add_argument("--candidates", required=True)
    population.add_argument("--out-dir", required=True)
    population.add_argument("--minimum-mapq", type=int, default=20)
    population.add_argument("--minimum-identity", type=float, default=0.9)
    population.add_argument("--config", default="configs/hsv2_pilot.yaml")
    population.add_argument("--open-report", action="store_true")
    population.set_defaults(func=_analyze_population)

    virtual_knockout = analyze_commands.add_parser(
        "virtual-knockout",
        help="enumerate bounded coding-sequence hypotheses for a generic virus project",
    )
    virtual_knockout.add_argument("--project", required=True)
    virtual_knockout.set_defaults(func=_analyze_virtual_knockout)

    escape = analyze_commands.add_parser(
        "escape",
        help="measure observed target support and exact-target sequence counterfactuals",
    )
    escape.add_argument("--project", required=True)
    escape.set_defaults(func=_analyze_escape)

    multiplex = analyze_commands.add_parser(
        "multiplex",
        help="compare configured panels using a sequence-level escape barrier",
    )
    multiplex.add_argument("--project", required=True)
    multiplex.add_argument("--open-report", action="store_true")
    multiplex.set_defaults(func=_analyze_multiplex)

    profiles = subparsers.add_parser("profiles", help="generic research-profile operations")
    profile_commands = profiles.add_subparsers(dest="profiles_command", required=True)
    profile_validate = profile_commands.add_parser(
        "validate", help="validate virus, host, and nuclease profiles"
    )
    profile_validate.add_argument("--virus-profile", required=True)
    profile_validate.add_argument("--host-profile", required=True)
    profile_validate.add_argument("--nuclease-profile", required=True)
    profile_validate.add_argument("--project-root", default=".")
    profile_validate.add_argument("--require-host-reference", action="store_true")
    profile_validate.add_argument("--require-virus-inputs", action="store_true")
    profile_validate.set_defaults(func=_profiles_validate)

    showcase = subparsers.add_parser(
        "showcase", help="build a presentation-ready standardized case study"
    )
    showcase_commands = showcase.add_subparsers(dest="showcase_command", required=True)
    showcase_build = showcase_commands.add_parser(
        "build", help="build multi-objective panels, figures, and reports"
    )
    showcase_build.add_argument("--virus-profile", required=True)
    showcase_build.add_argument("--host-profile", required=True)
    showcase_build.add_argument("--nuclease-profile", required=True)
    showcase_build.add_argument("--project-root", default=".")
    showcase_build.add_argument("--out-dir", default="reports/hsv2_showcase")
    showcase_build.add_argument("--per-gene", type=int, default=4)
    showcase_build.add_argument("--open-report", action="store_true")
    showcase_build.set_defaults(func=_showcase_build)

    scan = subparsers.add_parser("scan", help="scan and rank an aligned viral genome collection")
    scan.add_argument("--virus-alignment", required=True)
    scan.add_argument("--reference-id", required=True)
    scan.add_argument("--gff")
    scan.add_argument("--gene-evidence")
    scan.add_argument("--small-host-fasta", help="teaching-scale FASTA only; guarded at 5 Mb")
    scan.add_argument("--out-dir", required=True)
    scan.add_argument("--min-coverage", type=float, default=0.0)
    scan.add_argument("--max-mismatches", type=int, default=3)
    _add_config(scan)
    scan.set_defaults(func=_scan)

    rank = subparsers.add_parser("rank", help="rank an existing candidate CSV")
    rank.add_argument("--candidates", required=True)
    rank.add_argument("--gene-evidence")
    rank.add_argument("--out-dir", required=True)
    _add_config(rank)
    rank.set_defaults(func=_rank)

    build = subparsers.add_parser("build-offtarget-input", help="build Cas-OFFinder input")
    build.add_argument("--candidates", required=True)
    build.add_argument("--human-fasta-directory", required=True)
    build.add_argument("--output", required=True)
    build.add_argument("--manifest", required=True)
    build.add_argument("--max-candidates", type=int)
    build.add_argument("--genes", nargs="*")
    build.add_argument("--no-stratify", action="store_true")
    _add_config(build)
    build.set_defaults(func=_build_offtarget)

    summarize = subparsers.add_parser("summarize-offtargets", help="summarize Cas-OFFinder hits")
    summarize.add_argument("--candidates", required=True)
    summarize.add_argument("--cas-output", required=True)
    summarize.add_argument("--manifest")
    summarize.add_argument("--human-gff")
    summarize.add_argument("--out-dir", required=True)
    _add_config(summarize)
    summarize.set_defaults(func=_summarize_offtargets)

    pairs = subparsers.add_parser("simulate-pairs", help="rank explicit pair hypotheses")
    pairs.add_argument("--candidates", required=True)
    pairs.add_argument("--gff")
    pairs.add_argument("--virus-alignment")
    pairs.add_argument("--reference-id")
    output_group = pairs.add_mutually_exclusive_group(required=True)
    output_group.add_argument("--output")
    output_group.add_argument("--out-dir")
    pairs.add_argument("--min-distance", type=int)
    pairs.add_argument("--max-distance", type=int)
    pairs.add_argument("--max-candidates", type=int)
    pairs.add_argument("--maximum-candidates-per-gene", type=int)
    pairs.add_argument("--maximum-viral-occurrence-count", type=int)
    pairs.add_argument("--minimum-conservation", type=float)
    pairs.add_argument("--genes", nargs="*")
    pairs.add_argument("--feature-types", nargs="*")
    pairs.add_argument("--allow-cross-feature", action="store_true")
    pairs.add_argument("--no-stratify", action="store_true")
    _add_config(pairs)
    pairs.set_defaults(func=_simulate_pairs)

    benchmark = subparsers.add_parser("benchmark", help="measure known-target recovery")
    benchmark.add_argument("--candidates", required=True)
    benchmark.add_argument("--known-targets", required=True)
    benchmark.add_argument("--out-dir", required=True)
    benchmark.set_defaults(func=_benchmark)

    report = subparsers.add_parser("report", help="generate researcher-facing report outputs")
    report.add_argument("--candidates", required=True)
    report.add_argument("--rejected")
    report.add_argument("--pairs")
    report.add_argument("--multi-pairs")
    report.add_argument("--predicted-hits")
    report.add_argument("--out-dir", required=True)
    report.add_argument("--title", default="ViralSafeTarget research report")
    report.set_defaults(func=_report)

    tools = subparsers.add_parser("tools", help="external-tool adapters and imports")
    tool_commands = tools.add_subparsers(dest="tools_command", required=True)
    tools_status = tool_commands.add_parser("status", help="show external-tool availability")
    tools_status.add_argument("--project")
    tools_status.add_argument("--json", action="store_true")
    tools_status.set_defaults(func=_researcher_doctor)
    tools_setup = tool_commands.add_parser(
        "setup", help="show supported external-tool installation resources"
    )
    tools_setup.set_defaults(func=_tools_setup)
    tools_doctor = tool_commands.add_parser("doctor", help="detect supported external tools")
    tools_doctor.set_defaults(func=_tools_doctor)

    tools_benchmark = tool_commands.add_parser(
        "benchmark",
        help="run a frozen-panel multi-tool benchmark with explicit missingness",
    )
    tools_benchmark.add_argument("--config", required=True)
    tools_benchmark.add_argument("--open-report", action="store_true")
    tools_benchmark.set_defaults(func=_tool_benchmark)

    crispritz = tool_commands.add_parser("crispritz", help="CRISPRitz integration")
    crispritz_commands = crispritz.add_subparsers(dest="crispritz_command", required=True)
    crispritz_build = crispritz_commands.add_parser(
        "build-input", help="build a resumable CRISPRitz input bundle"
    )
    crispritz_build.add_argument("--run-dir", required=True)
    crispritz_build.add_argument("--out-dir", required=True)
    crispritz_build.add_argument("--reference-genome")
    crispritz_build.add_argument("--assembly", default="GRCh38.p14")
    crispritz_build.add_argument("--pam-file")
    crispritz_build.add_argument("--mismatches", type=int, default=3)
    crispritz_build.add_argument("--dna-bulges", type=int, default=0)
    crispritz_build.add_argument("--rna-bulges", type=int, default=0)
    crispritz_build.add_argument("--variant-aware", action="store_true")
    crispritz_build.add_argument("--variant-file")
    crispritz_build.add_argument("--docker-image", default="pinellolab/crispritz:latest")
    crispritz_build.set_defaults(func=_crispritz_build)

    crispritz_run = crispritz_commands.add_parser("run", help="run native or Docker CRISPRitz")
    crispritz_run.add_argument("--input-dir", required=True)
    crispritz_run.add_argument("--out-dir", required=True)
    crispritz_run.add_argument("--docker-image", default="pinellolab/crispritz:latest")
    crispritz_run.add_argument("--dry-run", action="store_true")
    crispritz_run.set_defaults(func=_crispritz_run)

    crispritz_import = crispritz_commands.add_parser(
        "import", help="import an existing CRISPRitz output"
    )
    crispritz_import.add_argument("--results", required=True)
    crispritz_import.add_argument("--manifest", required=True)
    crispritz_import.add_argument("--out-dir", required=True)
    crispritz_import.add_argument("--docker-image", default="pinellolab/crispritz:latest")
    crispritz_import.set_defaults(func=_crispritz_import)

    external_import = tool_commands.add_parser(
        "import", help="import CRISPOR, CHOPCHOP, or GuideScan2 exports"
    )
    external_import.add_argument(
        "--tool", required=True, choices=["crispor", "chopchop", "guidescan2"]
    )
    external_import.add_argument("--input", required=True)
    external_import.add_argument("--mapping", required=True)
    external_import.add_argument("--candidates", required=True)
    external_import.add_argument("--out-dir", required=True)
    external_import.set_defaults(func=_external_import)

    experimental = subparsers.add_parser(
        "experimental", help="measured-result imports kept separate from predictions"
    )
    experimental_commands = experimental.add_subparsers(dest="experimental_command", required=True)
    crispresso = experimental_commands.add_parser(
        "import-crispresso2", help="import an existing CRISPResso2 result directory"
    )
    crispresso.add_argument("--input", required=True)
    crispresso.add_argument("--candidate-map", required=True)
    crispresso.add_argument("--out-dir", required=True)
    crispresso.set_defaults(func=_crispresso_import)

    doctor = subparsers.add_parser("doctor", help="inspect the local research environment")
    doctor.add_argument("--project")
    doctor.add_argument("--json", action="store_true")
    doctor.set_defaults(func=_researcher_doctor)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
