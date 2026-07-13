"""Command-line interface for the research pipeline."""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

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
from .provenance import write_run_manifest
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


def _add_config(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--config", default=str(DEFAULT_CONFIG_PATH))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="vst",
        description="Virus-first computational target discovery for research use",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

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
    tools_doctor = tool_commands.add_parser("doctor", help="detect supported external tools")
    tools_doctor.set_defaults(func=_tools_doctor)

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
    _add_config(doctor)
    doctor.set_defaults(func=_doctor)
    return parser


def main() -> None:
    args = build_parser().parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
