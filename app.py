"""Local Streamlit interface for ViralSafeTarget.

Run with:
    streamlit run app.py
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import streamlit as st

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from viral_safe_target import (  # noqa: E402
    annotate_candidates,
    rank_candidate_pairs,
    rank_candidates,
    read_fasta,
    read_gff3,
    scan_spcas9_candidates,
    screen_against_small_fasta,
)

st.set_page_config(page_title="ViralSafeTarget", layout="wide")
st.title("ViralSafeTarget")
st.caption("Conserved viral target prioritization — computational research only")
st.warning(
    "This interface does not prove that a virus is disabled, does not establish safety, "
    "and is not a protocol for laboratory or clinical use."
)

with st.expander("What files do I upload?", expanded=True):
    st.markdown(
        """
- **Required:** an already aligned multi-FASTA containing several viral genomes.
- **Optional:** a GFF3 annotation matching the selected reference genome.
- **Optional demo only:** a small host FASTA (up to 5 Mb). For a human genome use
  Cas-OFFinder/CRISPRitz outside this app.
        """
    )

alignment_file = st.file_uploader("Aligned viral genomes (FASTA)", type=["fa", "fasta", "fna"])
gff_file = st.file_uploader("Reference annotation (GFF3, optional)", type=["gff", "gff3"])
host_file = st.file_uploader(
    "Small host FASTA for a smoke test (optional, not GRCh38)", type=["fa", "fasta", "fna"]
)
min_coverage = st.slider("Minimum exact 23-nt site coverage", 0.0, 1.0, 0.95, 0.01)

if alignment_file:
    with tempfile.TemporaryDirectory() as temp_directory:
        temp = Path(temp_directory)
        alignment_path = temp / alignment_file.name
        alignment_path.write_bytes(alignment_file.getvalue())
        records = read_fasta(alignment_path)
        reference_id = st.selectbox("Reference record", list(records.keys()))

        if st.button("Scan candidates", type="primary"):
            candidates = scan_spcas9_candidates(records, reference_id, min_coverage)
            features = None
            if gff_file:
                gff_path = temp / gff_file.name
                gff_path.write_bytes(gff_file.getvalue())
                features = read_gff3(gff_path)
                candidates = annotate_candidates(candidates, features, seqid=reference_id)
            if host_file:
                host_path = temp / host_file.name
                host_path.write_bytes(host_file.getvalue())
                host = read_fasta(host_path)
                candidates = screen_against_small_fasta(candidates, host)
            candidates = rank_candidates(candidates)

            st.session_state["candidates"] = candidates
            st.session_state["features"] = features
            st.session_state["records"] = records
            st.session_state["reference_id"] = reference_id

if "candidates" in st.session_state:
    candidates = st.session_state["candidates"]
    st.subheader(f"Candidate sites: {len(candidates):,}")
    st.dataframe(candidates, use_container_width=True, height=420)
    st.download_button(
        "Download candidates.csv",
        data=candidates.to_csv(index=False).encode("utf-8"),
        file_name="candidates.csv",
        mime="text/csv",
    )

    st.subheader("Idealized two-cut sequence simulation")
    st.caption(
        "This only describes the deletion coordinates that would result if both canonical cuts "
        "occurred. It is not a prediction that editing or viral inactivation will occur."
    )
    max_candidates = st.number_input(
        "Maximum ranked, gene-stratified candidates considered for pairing",
        min_value=10,
        max_value=500,
        value=100,
        step=10,
    )
    if st.button("Simulate candidate pairs"):
        pairs = rank_candidate_pairs(
            candidates,
            features=st.session_state.get("features"),
            aligned_records=st.session_state.get("records"),
            reference_id=st.session_state.get("reference_id"),
            max_candidates=int(max_candidates),
        )
        st.session_state["pairs"] = pairs

if "pairs" in st.session_state:
    pairs = st.session_state["pairs"]
    st.dataframe(pairs, use_container_width=True, height=420)
    st.download_button(
        "Download simulated_pairs.csv",
        data=pairs.to_csv(index=False).encode("utf-8"),
        file_name="simulated_pairs.csv",
        mime="text/csv",
    )
