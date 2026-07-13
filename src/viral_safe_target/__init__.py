"""ViralSafeTarget computational research package."""

__version__ = "0.2.0"

from .annotations import annotate_candidates, read_gff3
from .conservation import conservation_profile, find_conserved_runs
from .crispr import reverse_complement, scan_spcas9_candidates
from .disruption import (
    exact_pair_coverage,
    hypothetical_indel_consequences,
    rank_candidate_pairs,
    simulate_candidate_pair,
    spcas9_cut_after_1based,
)
from .io_utils import read_fasta, write_fasta
from .offtarget import (
    read_cas_offinder_output,
    screen_against_small_fasta,
    summarize_cas_offinder_hits,
    write_cas_offinder_input,
)
from .provenance import sha256_file, write_run_manifest
from .reporting import write_html_report
from .scoring import rank_candidates

__all__ = [
    "__version__",
    "read_fasta",
    "write_fasta",
    "conservation_profile",
    "find_conserved_runs",
    "scan_spcas9_candidates",
    "reverse_complement",
    "read_gff3",
    "annotate_candidates",
    "screen_against_small_fasta",
    "write_cas_offinder_input",
    "read_cas_offinder_output",
    "summarize_cas_offinder_hits",
    "rank_candidates",
    "write_html_report",
    "spcas9_cut_after_1based",
    "hypothetical_indel_consequences",
    "simulate_candidate_pair",
    "rank_candidate_pairs",
    "exact_pair_coverage",
    "sha256_file",
    "write_run_manifest",
]
