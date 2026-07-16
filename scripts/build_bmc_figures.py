#!/usr/bin/env python3
"""Rebuild BMC submission figures from committed machine-readable sources."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
FINAL = ROOT / "paper/bmc_bioinformatics/final"
FIGURES = FINAL / "figures"
ADDITIONAL = FINAL / "additional_files"

BLUE = "#1769AA"
NAVY = "#16324F"
TEAL = "#168C8C"
ORANGE = "#E07A2D"
LIGHT = "#EAF2F8"
GRAY = "#6B7280"
RED = "#B42318"

EXPECTED_COUNTS = {
    "initial_candidates": 28578,
    "eligible_candidates": 23108,
    "unique_guides": 21654,
    "completed_batches": 109,
    "total_batches": 109,
    "human_matches": 440341,
    "zero_hit_rows": 2668,
    "deep_guides": 257,
    "cds_mappings": 271,
    "mapped_guides": 250,
    "indel_hypotheses": 5691,
    "snv_counterfactuals": 17733,
}

RAW_SOURCE_SHA256 = {
    "reports/hsv2_genome_wide_exhaustive/provenance.json": (
        "dfe88939a5bccfc7b9d1b524f16c04a4f35bcbd5108f900aa99640d1d1e75cb7"
    ),
    "reports/hsv2_genome_wide_exhaustive/genome_wide_candidates_post_human.csv": (
        "165c9706a3c8bc7a64b0cbb6821b0c637e7b2bd9546ac530afbdf37c70aefeac"
    ),
    "reports/hsv2_genome_wide_exhaustive/genome_wide_human_hits.csv": (
        "4fadded91d229e21026bb675b3f228ddfb2fffcc1482a0cb9ddc8547a73918b8"
    ),
}


def _style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 11,
            "axes.labelsize": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _save(fig: plt.Figure, number: int, *, supplementary: bool = False) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    ADDITIONAL.mkdir(parents=True, exist_ok=True)
    if supplementary:
        pdf = ADDITIONAL / "Additional_file_1_Supplementary_Figure_S1.pdf"
        png = FIGURES / "Supplementary_Figure_S1.png"
    else:
        pdf = FIGURES / f"Figure_{number}.pdf"
        png = FIGURES / f"Figure_{number}.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_sources() -> dict[str, float | int]:
    raw_paths = {relative: ROOT / relative for relative in RAW_SOURCE_SHA256}
    raw_available = {relative: path.exists() for relative, path in raw_paths.items()}
    if any(raw_available.values()) and not all(raw_available.values()):
        missing = [relative for relative, available in raw_available.items() if not available]
        raise FileNotFoundError(f"Partial raw publication source set; missing: {missing}")

    if all(raw_available.values()):
        for relative, expected_hash in RAW_SOURCE_SHA256.items():
            observed_hash = _sha256(raw_paths[relative])
            if observed_hash != expected_hash:
                raise AssertionError(
                    f"Raw publication source checksum changed: {relative}: {observed_hash}"
                )
        provenance = json.loads(raw_paths[next(iter(RAW_SOURCE_SHA256))].read_text())
        post = pd.read_csv(
            raw_paths["reports/hsv2_genome_wide_exhaustive/genome_wide_candidates_post_human.csv"],
            usecols=["guide_sequence", "human_total_predicted_hits"],
        )
        hits = pd.read_csv(
            raw_paths["reports/hsv2_genome_wide_exhaustive/genome_wide_human_hits.csv"],
            usecols=["candidate_id"],
        )
        primary = {
            "initial_candidates": int(provenance["initial_candidate_count"]),
            "eligible_candidates": len(post),
            "unique_guides": int(post["guide_sequence"].nunique()),
            "completed_batches": int(provenance["completed_batches"]),
            "total_batches": int(provenance["total_batches"]),
            "human_matches": len(hits),
            "zero_hit_rows": int(pd.to_numeric(post["human_total_predicted_hits"]).eq(0).sum()),
        }
    else:
        snapshot_path = FINAL / "verified_statistics.json"
        if not snapshot_path.exists():
            raise FileNotFoundError(
                "Raw sources are not present and the committed verified snapshot is missing"
            )
        snapshot = json.loads(snapshot_path.read_text())
        primary = {
            key: int(snapshot[key])
            for key in (
                "initial_candidates",
                "eligible_candidates",
                "unique_guides",
                "completed_batches",
                "total_batches",
                "human_matches",
                "zero_hit_rows",
            )
        }

    virtual = json.loads(
        (ROOT / "reports/hsv2_virtual_knockout_escape/run_manifest.json").read_text()
    )
    benchmark = pd.read_csv(ROOT / "reports/hsv2_tool_benchmark/rank_agreement.csv")
    overlap = pd.read_csv(ROOT / "reports/hsv2_tool_benchmark/top_k_overlap.csv")

    observed = {
        **primary,
        "deep_guides": int(virtual["outputs"]["virtual_knockout"]["candidate_count"]),
        "cds_mappings": int(virtual["outputs"]["virtual_knockout"]["mapping_row_count"]),
        "mapped_guides": int(virtual["outputs"]["virtual_knockout"]["mapped_candidate_count"]),
        "indel_hypotheses": int(virtual["outputs"]["virtual_knockout"]["hypothesis_count"]),
        "snv_counterfactuals": int(virtual["outputs"]["escape"]["counterfactual_count"]),
    }
    if observed != EXPECTED_COUNTS:
        raise AssertionError(
            f"Publication source mismatch: observed={observed}, expected={EXPECTED_COUNTS}"
        )

    row = benchmark[(benchmark.tool_a == "cas-offinder") & (benchmark.tool_b == "crispritz")].iloc[
        0
    ]
    if not np.isclose(row.spearman_rank_correlation, 0.880413, atol=5e-7):
        raise AssertionError("Cas-OFFinder/CRISPRitz correlation changed")
    prepost = benchmark[
        (benchmark.tool_a == "viral_safe_target_post_human")
        & (benchmark.tool_b == "viral_safe_target_pre_human")
    ].iloc[0]
    if not np.isclose(prepost.spearman_rank_correlation, 0.961975, atol=5e-7):
        raise AssertionError("ViralSafeTarget pre/post correlation changed")
    expected_overlap = {10: 9, 25: 24, 50: 49}
    subset = overlap[(overlap.tool_a == "cas-offinder") & (overlap.tool_b == "crispritz")]
    actual_overlap = dict(
        zip(subset.top_k.astype(int), subset.overlap_count.astype(int), strict=True)
    )
    if actual_overlap != expected_overlap:
        raise AssertionError(f"Top-k overlap changed: {actual_overlap}")

    summary = {**observed, "cas_crispritz_spearman": float(row.spearman_rank_correlation)}
    summary["vst_pre_post_spearman"] = float(prepost.spearman_rank_correlation)
    (FINAL / "verified_statistics.json").parent.mkdir(parents=True, exist_ok=True)
    (FINAL / "verified_statistics.json").write_text(json.dumps(summary, indent=2) + "\n")
    (FINAL / "raw_source_checksums.json").write_text(json.dumps(RAW_SOURCE_SHA256, indent=2) + "\n")
    return summary


def figure_1() -> None:
    fig, ax = plt.subplots(figsize=(12, 5.8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6)
    ax.axis("off")
    ax.text(
        6,
        5.72,
        "ViralSafeTarget: auditable virus-first decision support",
        ha="center",
        fontsize=15,
        weight="bold",
        color=NAVY,
    )
    columns = [
        (
            0.3,
            "INPUTS",
            [
                "Reference FASTA",
                "Viral strains",
                "GFF3 annotation",
                "Host assembly",
                "Evidence sources",
            ],
            BLUE,
        ),
        (
            3.1,
            "SEQUENCE",
            [
                "Candidate discovery",
                "Cross-strain support",
                "One-to-many mapping",
                "Host-search status",
            ],
            TEAL,
        ),
        (
            6.0,
            "INTERPRETATION",
            [
                "Guide rank",
                "Gene targetability",
                "Coding hypotheses",
                "Escape counterfactuals",
                "Evidence review",
            ],
            ORANGE,
        ),
        (
            9.1,
            "OUTPUTS",
            [
                "Research shortlist",
                "Multiplex panels",
                "Per-guide explanation",
                "START_HERE report",
                "Portable export",
            ],
            NAVY,
        ),
    ]
    widths = [2.2, 2.25, 2.45, 2.35]
    for (x, title, items, color), width in zip(columns, widths, strict=True):
        rect = FancyBboxPatch(
            (x, 1.1),
            width,
            3.95,
            boxstyle="round,pad=0.05,rounding_size=0.08",
            fc="white",
            ec=color,
            lw=1.8,
        )
        ax.add_patch(rect)
        ax.add_patch(
            FancyBboxPatch(
                (x, 4.55),
                width,
                0.5,
                boxstyle="round,pad=0.02,rounding_size=0.06",
                fc=color,
                ec=color,
            )
        )
        ax.text(x + width / 2, 4.8, title, ha="center", va="center", color="white", weight="bold")
        for i, item in enumerate(items):
            y = 4.18 - i * 0.62
            ax.text(x + 0.18, y, "●", color=color, va="center", fontsize=7)
            ax.text(x + 0.38, y, item, color=NAVY, va="center", fontsize=9)
    for x1, x2 in [(2.5, 3.1), (5.35, 6.0), (8.45, 9.1)]:
        ax.add_patch(
            FancyArrowPatch(
                (x1, 3.05), (x2, 3.05), arrowstyle="-|>", mutation_scale=14, lw=1.5, color=GRAY
            )
        )
    ax.add_patch(
        FancyBboxPatch((0.9, 0.25), 10.2, 0.55, boxstyle="round,pad=0.03", fc=LIGHT, ec="none")
    )
    ax.text(
        6,
        0.52,
        "Guide quality  ≠  gene targetability  ≠  biological evidence  ≠  escape robustness",
        ha="center",
        va="center",
        weight="bold",
        color=NAVY,
        fontsize=11,
    )
    _save(fig, 1)


def figure_2(stats: dict[str, float | int]) -> None:
    labels = ["Initial\ncandidates", "Eligible\nrows", "Unique\nguides", "Zero-hit\nrows*"]
    values = [
        stats["initial_candidates"],
        stats["eligible_candidates"],
        stats["unique_guides"],
        stats["zero_hit_rows"],
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    bars = ax.bar(labels, values, color=[NAVY, BLUE, TEAL, ORANGE])
    ax.set_ylabel("Count")
    ax.set_title("Exhaustive HSV-2 computational funnel")
    ax.grid(axis="y", alpha=0.2)
    for bar, value in zip(bars, values, strict=True):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            value + max(values) * 0.025,
            f"{int(value):,}",
            ha="center",
            weight="bold",
        )
    ax.text(
        0.01,
        -0.18,
        "*No predicted GRCh38.p14 hit under the configured SpCas9 search; not evidence of safety.",
        transform=ax.transAxes,
        fontsize=8,
        color=GRAY,
    )
    _save(fig, 2)


def figure_3() -> None:
    genes = pd.read_csv(ROOT / "reports/hsv2_genome_wide_exhaustive/gene_rankings.csv")
    genes = genes.sort_values("targetability_rank").head(15).copy()
    fig, ax = plt.subplots(figsize=(8.0, 5.5))
    colors = [ORANGE if g in {"UL3", "UL36", "UL30", "UL52"} else BLUE for g in genes.gene_name]
    ax.barh(genes.gene_name[::-1], genes.targetability_score[::-1], color=colors[::-1])
    ax.set_xlabel("Gene targetability score")
    ax.set_title("Leading HSV-2 gene-targetability portfolios")
    ax.set_xlim(0, 1)
    for y, (_, row) in enumerate(genes.iloc[::-1].iterrows()):
        ax.text(
            row.targetability_score + 0.012,
            y,
            f"rank {int(row.targetability_rank)}",
            va="center",
            fontsize=8,
        )
    ax.text(
        0.02,
        -0.15,
        "Orange: discussed focus genes. Top individual guide: UL36; "
        "UL36 gene portfolio rank: 21 (outside plotted top 15).",
        transform=ax.transAxes,
        color=GRAY,
        fontsize=8,
    )
    _save(fig, 3)


def figure_4() -> None:
    gene = pd.read_csv(ROOT / "reports/hsv2_virtual_knockout_escape/gene_virtual_knockout.csv")
    focus = ["UL3", "UL10", "UL52", "UL30", "UL18", "UL36"]
    d = gene[gene.gene_name.isin(focus)].set_index("gene_name").reindex(focus)
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    x = np.arange(len(d))
    ax.bar(x - 0.18, d.mapped_guide_count, 0.36, label="Mapped guides", color=BLUE)
    ax.bar(
        x + 0.18,
        d.guides_cutting_annotated_domain,
        0.36,
        label="Domain-overlapping guides",
        color=ORANGE,
    )
    ax.set_xticks(x, d.index)
    ax.set_ylabel("Guide mappings")
    ax.set_title("Coding and domain context in focus genes")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)
    _save(fig, 4)


def figure_5() -> None:
    d = pd.read_csv(ROOT / "reports/hsv2_virtual_knockout_escape/strategy_comparison.csv")
    fig, ax = plt.subplots(figsize=(8.3, 4.7))
    labels = [s.replace("-", "\n") for s in d.strategy]
    bars = ax.bar(labels, d.sequence_escape_barrier, color=[NAVY, BLUE, TEAL, ORANGE])
    ax.set_ylim(0, max(d.sequence_escape_barrier) + 1)
    ax.set_ylabel("Minimum distinct substitutions")
    ax.set_title("Configured exact-target multiplex barriers")
    for bar, row in zip(bars, d.itertuples(), strict=True):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            row.sequence_escape_barrier + 0.08,
            f"barrier {row.sequence_escape_barrier}\nheld-out min {row.heldout_support_min:.3f}",
            ha="center",
            fontsize=8,
        )
    ax.text(
        0.01,
        -0.23,
        "Sequence-level exact-target set cover only; not an evolutionary escape probability.",
        transform=ax.transAxes,
        fontsize=8,
        color=GRAY,
    )
    _save(fig, 5)


def figure_6() -> None:
    agreement = pd.read_csv(ROOT / "reports/hsv2_tool_benchmark/rank_agreement.csv")
    overlap = pd.read_csv(ROOT / "reports/hsv2_tool_benchmark/top_k_overlap.csv")
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.5))
    display_names = {
        "cas-offinder": "Cas-OFFinder",
        "crispritz": "CRISPRitz",
        "viral_safe_target_post_human": "VST post-host",
        "viral_safe_target_pre_human": "VST pre-host",
    }
    labels = [
        f"{display_names[a]}\nvs\n{display_names[b]}"
        for a, b in zip(agreement.tool_a, agreement.tool_b, strict=True)
    ]
    axes[0].bar(
        range(len(agreement)),
        agreement.spearman_rank_correlation,
        color=[ORANGE] + [BLUE] * (len(agreement) - 2) + [TEAL],
    )
    axes[0].set_xticks(range(len(agreement)), labels, rotation=35, ha="right", fontsize=7)
    axes[0].set_ylim(0, 1)
    axes[0].set_ylabel("Spearman correlation")
    axes[0].set_title("Pairwise rank agreement")
    sub = overlap[(overlap.tool_a == "cas-offinder") & (overlap.tool_b == "crispritz")]
    axes[1].bar(sub.top_k.astype(str), sub.overlap_count, color=TEAL)
    axes[1].plot(sub.top_k.astype(str), sub.top_k, color=NAVY, marker="o", label="Maximum possible")
    axes[1].set_ylabel("Shared guides")
    axes[1].set_xlabel("Top-k")
    axes[1].set_title("Cas-OFFinder–CRISPRitz overlap")
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    _save(fig, 6)


def supplementary_1() -> None:
    d = pd.read_csv(ROOT / "reports/hsv2_tool_benchmark/ablation_summary.csv")
    d = d[d.variant != "all_components"].sort_values("median_absolute_rank_shift")
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    labels = [x.replace("without_", "− ").replace("_", " ") for x in d.variant]
    ax.barh(labels, d.median_absolute_rank_shift, color=[BLUE] * (len(d) - 1) + [ORANGE])
    ax.set_xlabel("Median absolute rank shift")
    ax.set_title("Leave-one-component-out sensitivity")
    ax.grid(axis="x", alpha=0.2)
    _save(fig, 1, supplementary=True)


def main() -> None:
    _style()
    stats = validate_sources()
    figure_1()
    figure_2(stats)
    figure_3()
    figure_4()
    figure_5()
    figure_6()
    supplementary_1()
    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
