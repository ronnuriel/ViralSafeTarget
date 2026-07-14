"""Researcher-facing report for gene function and predicted disruption."""

# HTML prose and column lists are intentionally kept as readable complete strings.
# ruff: noqa: E501

from __future__ import annotations

import html
from pathlib import Path

import pandas as pd


def _table(frame: pd.DataFrame, columns: list[str], rows: int = 50) -> str:
    available = [column for column in columns if column in frame]
    if frame.empty or not available:
        return '<p class="unknown">No result is available for this table.</p>'
    return (
        frame[available]
        .head(rows)
        .to_html(index=False, border=0, classes="data", render_links=True)
    )


def write_gene_function_report(
    output_path: str | Path,
    *,
    gene_scores: pd.DataFrame,
    evidence: pd.DataFrame,
    mapping: pd.DataFrame,
    outcomes: pd.DataFrame,
    domains: pd.DataFrame,
    disorder: pd.DataFrame,
    domain_overlap: pd.DataFrame,
    evolution: pd.DataFrame,
) -> None:
    paired = outcomes[outcomes["event_class"].eq("paired_guide_theoretical_deletion")]
    single = outcomes[
        outcomes["event_class"].eq("single_guide_indel") & outcomes["indel_size_bp"].ne(0)
    ]
    direct_hsv2 = evidence[evidence["virus_type"].eq("HSV-2")]
    hsv1 = evidence[evidence["virus_type"].eq("HSV-1")]
    domain_inside = domain_overlap[domain_overlap["relation_to_cut"].eq("cut_inside_domain")]
    sections = [
        (
            "Scope and answer",
            f"<p>This report maps {len(mapping):,} top coordinate-level candidates across nine HSV-2 genes to coding and protein coordinates, evaluates sequence conservation and imported annotations, and simulates {len(single):,} size-defined single-guide outcomes plus {len(paired):,} theoretical paired deletions.</p>"
            "<p><strong>These are model-bounded computational hypotheses.</strong> They do not establish editing, repair frequencies, gene essentiality, viral inhibition, safety, delivery, or clinical efficacy.</p>",
        ),
        (
            "Four scores kept separate",
            _table(
                gene_scores,
                [
                    "gene_name",
                    "top_candidate_count",
                    "sequence_targetability_score",
                    "evidence_based_essentiality_score",
                    "evidence_based_essentiality_scope",
                    "hsv2_evidence_based_essentiality_score",
                    "hsv2_essentiality_status",
                    "hsv1_ortholog_essentiality_score",
                    "hsv1_essentiality_status",
                    "predicted_protein_disruption_score",
                    "evidence_coverage_score",
                ],
            )
            + "<p>Scores are not combined into a single therapeutic rank. Missing evidence stays missing; it is not converted to zero. HSV-1 ortholog evidence is displayed separately and never relabeled as direct HSV-2 evidence.</p>",
        ),
        (
            "Candidate-to-protein mapping",
            _table(
                mapping.sort_values(["gene_name", "within_gene_rank"]),
                [
                    "gene_name",
                    "within_gene_rank",
                    "candidate_id",
                    "post_human_rank",
                    "cut_position",
                    "cds_strand",
                    "cds_position_start_1based",
                    "cds_position_end_1based",
                    "cut_cds_offset_0based",
                    "cut_phase",
                    "reference_codon",
                    "reference_amino_acid",
                    "amino_acid_position_1based",
                    "protein_length_aa",
                    "relative_protein_position",
                    "predicted_region_class",
                    "cut_domain_accessions",
                ],
                rows=100,
            ),
        ),
        (
            "Cross-strain and ortholog conservation",
            _table(
                evolution,
                [
                    "gene_name",
                    "hsv2_strain_count",
                    "mean_cds_nucleotide_entropy_bits",
                    "mean_amino_acid_conservation",
                    "invariant_amino_acid_fraction",
                    "median_pairwise_dN",
                    "median_pairwise_dS",
                    "median_pairwise_dN_dS",
                    "dN_dS_estimable_comparison_count",
                    "hsv1_ortholog_accession",
                    "hsv1_ortholog_protein_identity",
                    "hsv1_ortholog_reference_coverage",
                ],
            )
            + "<p>dN/dS is a pairwise NG86 summary only where unambiguous codons and a finite dS permit estimation. It is descriptive for this small convenience panel and is not a population-selection estimate.</p>",
        ),
        (
            "Direct HSV-2 evidence",
            _table(
                direct_hsv2,
                [
                    "gene_name",
                    "evidence_category",
                    "essentiality_call",
                    "essentiality_score",
                    "evidence_strength",
                    "experimental_system",
                    "finding",
                    "source_identifier",
                    "source_url",
                ],
            )
            + "<p>A knockdown phenotype is not treated as proof of null-mutant essentiality.</p>",
        ),
        (
            "HSV-1 ortholog evidence",
            _table(
                hsv1,
                [
                    "gene_name",
                    "evidence_category",
                    "essentiality_call",
                    "essentiality_score",
                    "evidence_strength",
                    "experimental_system",
                    "finding",
                    "source_identifier",
                    "source_url",
                ],
                rows=100,
            )
            + "<p>These studies concern HSV-1 and are retained as ortholog evidence, not direct HSV-2 essentiality assignments.</p>",
        ),
        (
            "Imported InterPro protein annotations",
            _table(
                domains,
                [
                    "gene_name",
                    "uniprot_accession",
                    "interpro_accession",
                    "entry_type",
                    "domain_name",
                    "protein_start_1based",
                    "protein_end_1based",
                    "source_url",
                ],
                rows=100,
            ),
        ),
        (
            "MobiDB-lite disordered regions",
            _table(
                disorder,
                [
                    "gene_name",
                    "uniprot_accession",
                    "region_type",
                    "protein_start_1based",
                    "protein_end_1based",
                    "method",
                    "source_url",
                ],
                rows=100,
            )
            + "<p>For proteins with a returned MobiDB-lite region set, the complement is labeled predicted structured complement. Absence of a usable region set remains unknown and is not called structured.</p>",
        ),
        (
            "Cuts inside imported domains",
            _table(
                domain_inside.sort_values(["gene_name", "candidate_id", "interpro_accession"]),
                [
                    "gene_name",
                    "candidate_id",
                    "amino_acid_position_1based",
                    "interpro_accession",
                    "interpro_entry_type",
                    "domain_name",
                    "domain_start_1based",
                    "domain_end_1based",
                    "relation_to_cut",
                ],
                rows=100,
            ),
        ),
        (
            "Small-indel model",
            "<p>Each top candidate is evaluated for sizes -10 through +10 bp. Negative values delete bases immediately downstream of the cut in coding orientation. Positive values insert an unspecified N placeholder. Frame status follows length modulo three. A downstream stop is reported only when it is predictable outside the ambiguous inserted segment.</p>"
            + _table(
                single.sort_values(["gene_name", "candidate_id", "indel_size_bp"]),
                [
                    "gene_name",
                    "candidate_id",
                    "indel_size_bp",
                    "frameshift",
                    "premature_stop_position_aa",
                    "premature_stop_status",
                    "retained_protein_fraction",
                    "affected_domain_accessions",
                ],
                rows=100,
            ),
        ),
        (
            "Paired-guide theoretical deletions",
            _table(
                paired.sort_values(["gene_name", "deleted_coding_bp"], ascending=[True, False]),
                [
                    "gene_name",
                    "candidate_a",
                    "candidate_b",
                    "deleted_coding_bp",
                    "theoretical_deletion_start_cds_1based",
                    "theoretical_deletion_end_cds_1based",
                    "frameshift",
                    "premature_stop_position_aa",
                    "retained_protein_fraction",
                    "affected_domain_accessions",
                ],
                rows=100,
            )
            + "<p>These rows describe cut-to-cut sequence consequences only. They do not predict that both cuts occur or that the deletion is produced.</p>",
        ),
        (
            "Limitations",
            "<ul><li>The 14-genome HSV-2 alignment is a convenience set, not a representative clinical population.</li><li>Candidate entropy and amino-acid conservation are alignment-bound.</li><li>InterPro and MobiDB-lite annotations are imported predictions/curation and may be incomplete.</li><li>Size-only insertion simulations cannot know inserted bases; sequence-dependent stops near the insertion are therefore unknown.</li><li>Essentiality calls require cited evidence and remain context-specific. Unknown is not evidence of nonessentiality.</li><li>No wet-lab protocol, clinical recommendation, safety conclusion, or claim of cure is provided.</li></ul>",
        ),
    ]
    body = "".join(
        f"<section><h2>{html.escape(title)}</h2>{content}</section>" for title, content in sections
    )
    document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>HSV-2 gene function and predicted disruption</title>
<style>
body{{font:16px/1.5 system-ui,sans-serif;max-width:1440px;margin:auto;padding:2rem;color:#15202b}}
h1,h2{{color:#183b56}} section{{margin:2.2rem 0}} table{{border-collapse:collapse;width:100%;font-size:.82rem;display:block;overflow:auto}}
th,td{{padding:.45rem;border:1px solid #d9e2ec;text-align:left;vertical-align:top;white-space:nowrap}} th{{background:#eef4f8}}
.warning{{background:#fff4e5;padding:1rem;border-left:5px solid #d97706}} .unknown{{color:#8a4b08}}
</style></head><body>
<h1>HSV-2 gene function and predicted disruption</h1>
<div class="warning"><strong>Research-use computational report.</strong> No outcome shown here is evidence of safety, efficacy, or a therapeutic intervention.</div>
{body}
</body></html>"""
    Path(output_path).write_text(document, encoding="utf-8")
