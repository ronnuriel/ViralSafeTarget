"""ViralSafeTarget computational research package."""

__version__ = "0.5.0"

from .annotations import annotate_candidates, read_gff3
from .config import EditorProfile, load_config
from .consensus import (
    ComparisonResult,
    build_consensus,
    candidate_metrics_as_tool_results,
    compare_tools,
)
from .conservation import conservation_profile, find_conserved_runs
from .crispr import (
    reverse_complement,
    scan_editor_candidates,
    scan_spcas9_candidates,
    stable_candidate_id,
)
from .disruption import (
    exact_pair_coverage,
    hypothetical_indel_consequences,
    rank_candidate_pairs,
    simulate_candidate_pair,
    spcas9_cut_after_1based,
)
from .integrations import ToolAdapter, load_external_results
from .io_utils import read_fasta, write_fasta
from .offtarget import (
    build_cas_offinder_input,
    read_cas_offinder_output,
    screen_against_small_fasta,
    summarize_cas_offinder_hits,
    write_cas_offinder_input,
)
from .provenance import sha256_file, write_run_manifest
from .reporting import write_html_report
from .scorers import CandidateScorer, ExampleRuleScorer
from .scoring import rank_candidates, rank_post_human_candidates, rank_pre_human_candidates
from .sdk import ResearchRun, load_run
from .tables import CandidateTable, ToolResultTable

__all__ = [
    "__version__",
    "read_fasta",
    "write_fasta",
    "conservation_profile",
    "find_conserved_runs",
    "scan_spcas9_candidates",
    "scan_editor_candidates",
    "stable_candidate_id",
    "EditorProfile",
    "load_config",
    "reverse_complement",
    "read_gff3",
    "annotate_candidates",
    "screen_against_small_fasta",
    "write_cas_offinder_input",
    "build_cas_offinder_input",
    "read_cas_offinder_output",
    "summarize_cas_offinder_hits",
    "rank_candidates",
    "rank_pre_human_candidates",
    "rank_post_human_candidates",
    "write_html_report",
    "spcas9_cut_after_1based",
    "hypothetical_indel_consequences",
    "simulate_candidate_pair",
    "rank_candidate_pairs",
    "exact_pair_coverage",
    "sha256_file",
    "write_run_manifest",
    "ResearchRun",
    "load_run",
    "CandidateTable",
    "ToolResultTable",
    "ToolAdapter",
    "CandidateScorer",
    "ExampleRuleScorer",
    "ComparisonResult",
    "compare_tools",
    "build_consensus",
    "candidate_metrics_as_tool_results",
    "load_external_results",
]
