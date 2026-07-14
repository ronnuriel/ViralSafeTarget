#!/usr/bin/env python3
"""Write final provenance after the HSV-2 pilot outputs are complete."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from viral_safe_target import load_config, write_run_manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--reports", type=Path, required=True)
    parser.add_argument("--alignment", type=Path, required=True)
    parser.add_argument("--human-fasta", type=Path, required=True)
    parser.add_argument("--human-gff", type=Path)
    parser.add_argument("--qc-report", type=Path, required=True)
    args = parser.parse_args()

    config = load_config(args.config)
    qc = pd.read_csv(args.qc_report)
    accepted = qc.loc[qc["decision"] == "accepted", "accession"].astype(str).tolist()
    rejected = qc.loc[
        qc["decision"] == "rejected", ["accession", "rejection_reason"]
    ].to_dict("records")
    post_path = args.reports / "candidates_ranked_post_human.csv"
    post = pd.read_csv(post_path) if post_path.is_file() else pd.DataFrame()
    hit_path = args.reports / "predicted_human_hits.csv"
    hits = pd.read_csv(hit_path) if hit_path.is_file() else pd.DataFrame()
    inputs = [
        Path("reports/real_hsv2/candidates_ranked_pre_human.csv"),
        args.alignment,
        Path("data/processed/hsv2_reference.gff3"),
        args.human_fasta,
        args.reports / "cas_offinder_input.txt",
    ]
    cas_output = args.reports / "cas_offinder_output.tsv"
    if cas_output.is_file():
        inputs.append(cas_output)
    if args.human_gff and args.human_gff.is_file():
        inputs.append(args.human_gff)
    outputs = [
        path
        for path in args.reports.iterdir()
        if path.is_file() and path.name != "run_manifest.json"
    ]
    write_run_manifest(
        args.reports / "run_manifest.json",
        inputs,
        {
            "pilot_genes": ["UL19", "UL30"],
            "maximum_candidates_per_gene": 100,
            "selected_candidate_count": 200,
            "post_human_candidate_count": len(post),
            "predicted_human_hit_count": len(hits),
            "maximum_mismatches": int(config["editor"]["mismatch_search_threshold"]),
            "pair_ranking_stage": "post_human" if not post.empty else "pre_human",
        },
        config_path=args.config,
        editor_profile=config["editor"],
        accepted_accessions=accepted,
        rejected_accessions=rejected,
        human_assembly_identifier=(
            f"{config['off_target']['human_assembly_accession']} / "
            f"{config['off_target']['human_assembly']}"
        ),
        command_line=["bash", "scripts/run_hsv2_pilot.sh"],
        random_seed=int(config["random_seed"]),
        output_paths=outputs,
    )
    print(f"Wrote {args.reports / 'run_manifest.json'}")


if __name__ == "__main__":
    main()
