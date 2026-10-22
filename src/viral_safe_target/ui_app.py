"""Researcher-facing Streamlit interface.

The interface orchestrates the public project workflow. It does not add a
scientific score and does not convert missing evidence into favorable values.
"""

from __future__ import annotations

import json
import tempfile
import time
from html import escape
from pathlib import Path
from typing import Any

import pandas as pd


def _streamlit():
    try:
        import streamlit as st
    except ImportError as error:  # pragma: no cover - exercised through the launcher
        raise RuntimeError(
            'Install the optional interface with: pip install "viral-safe-target[ui]"'
        ) from error
    return st


def _inject_styles(st: Any) -> None:
    st.markdown(
        """
<style>
    :root {
        --vst-navy: #102a43;
        --vst-blue: #176b87;
        --vst-teal: #20a39e;
        --vst-mist: #eef6f8;
        --vst-muted: #627d98;
        --vst-line: #d9e2ec;
    }
    .stApp { background: linear-gradient(180deg, #f7fbfc 0, #ffffff 320px); }
    .block-container { max-width: 1280px; padding-top: 1.7rem; padding-bottom: 4rem; }
    h1, h2, h3 { color: var(--vst-navy); letter-spacing: -0.02em; }
    .vst-hero {
        background: linear-gradient(128deg, #102a43 0%, #176b87 68%, #20a39e 100%);
        color: white; padding: 2.2rem 2.4rem; border-radius: 22px;
        box-shadow: 0 14px 34px rgba(16,42,67,.16); margin-bottom: 1.4rem;
    }
    .vst-hero h1 { color: white; margin: 0 0 .45rem; font-size: 2.45rem; }
    .vst-hero p { margin: 0; max-width: 850px; opacity: .92; font-size: 1.08rem; }
    .vst-eyebrow { text-transform: uppercase; letter-spacing: .12em; font-size: .76rem;
        font-weight: 700; opacity: .78; margin-bottom: .55rem; }
    .vst-boundary {
        background: #fff8e8; border-left: 5px solid #f0a202; color: #594a1f;
        padding: .85rem 1rem; border-radius: 8px; margin: .8rem 0 1.35rem;
    }
    .vst-grid { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr));
        gap: .75rem; margin: .8rem 0 1.2rem; }
    .vst-stage { background: white; border: 1px solid var(--vst-line); border-radius: 14px;
        padding: 1rem; min-height: 130px; box-shadow: 0 5px 16px rgba(16,42,67,.055); }
    .vst-stage .num { display: inline-flex; width: 28px; height: 28px; align-items: center;
        justify-content: center; border-radius: 50%; background: var(--vst-mist);
        color: var(--vst-blue); font-weight: 700; margin-bottom: .65rem; }
    .vst-stage strong { display: block; color: var(--vst-navy); margin-bottom: .3rem; }
    .vst-stage span { color: var(--vst-muted); font-size: .9rem; line-height: 1.35; }
    .vst-note { color: var(--vst-muted); font-size: .92rem; }
    .vst-status { display:inline-block; border-radius:999px; padding:.2rem .58rem;
        font-weight:650; font-size:.78rem; }
    .vst-completed { color:#116149; background:#e6f6f0; }
    .vst-running { color:#0b5793; background:#e6f2ff; }
    .vst-pending { color:#765c10; background:#fff6d9; }
    .vst-required { color:#873c32; background:#fdece9; }
    .vst-muted-status { color:#52606d; background:#edf2f7; }
    .vst-stage.running { border-color:#2b7bbb; box-shadow:0 0 0 2px rgba(43,123,187,.12); }
    .vst-stage.completed { border-color:#73c6ad; }
    .vst-stage.external_required, .vst-stage.failed { border-color:#ef9a8d; }
    .vst-empty { text-align:center; padding:2rem; border:1px dashed var(--vst-line);
        border-radius:14px; color:var(--vst-muted); background:rgba(255,255,255,.7); }
    [data-testid="stMetric"] { background:white; border:1px solid var(--vst-line);
        border-radius:14px; padding:.8rem 1rem; }
    [data-testid="stSidebar"] { background:#f1f7f8; border-right:1px solid var(--vst-line); }
    .stButton > button, .stDownloadButton > button { border-radius:10px; font-weight:650; }
    .stButton button[kind="primary"], .stFormSubmitButton button[kind="primary"] {
        background: var(--vst-blue); border-color: var(--vst-blue); color: white;
    }
    @media (max-width: 900px) { .vst-grid { grid-template-columns: repeat(2, minmax(0,1fr)); } }
    @media (max-width: 560px) { .vst-grid { grid-template-columns: 1fr; }
        .vst-hero { padding:1.45rem; } .vst-hero h1 { font-size:1.9rem; } }
</style>
        """,
        unsafe_allow_html=True,
    )


def _hero(st: Any) -> None:
    st.markdown(
        """
<section class="vst-hero">
  <div class="vst-eyebrow">Virus-first research workflow</div>
  <h1>ViralSafeTarget</h1>
  <p>Move from viral genomes to an auditable research shortlist while keeping
  conservation, host-risk, gene context, disruption, escape, and evidence as
  separate, inspectable decision layers.</p>
</section>
<div class="vst-boundary"><b>Research boundary.</b> Outputs are computational
hypotheses for independent evaluation. They do not establish editing, safety,
viral inhibition, delivery, treatment, clearance, or cure.</div>
        """,
        unsafe_allow_html=True,
    )


def _pipeline(st: Any) -> None:
    stages = [
        ("1", "Inputs & QC", "Reference, strain panel, annotation, host, and editor profile."),
        ("2", "Guide discovery", "Enumerate editor-compatible sites across the viral genome."),
        ("3", "Conservation", "Measure exact guide and PAM support across viral strains."),
        ("4", "Annotation", "Map sites to genes, CDS, proteins, and documented domains."),
        ("5", "Host-risk", "Preserve explicit off-target status and external-tool boundaries."),
        ("6", "Disruption & escape", "Model bounded sequence consequences and target loss."),
        ("7", "Evidence review", "Separate direct-virus, ortholog, pending, and unknown evidence."),
        ("8", "Shortlist & export", "Compare multiplex panels and package auditable outputs."),
    ]
    cards = "".join(
        f'<div class="vst-stage"><div class="num">{number}</div><strong>{title}</strong>'
        f"<span>{detail}</span></div>"
        for number, title, detail in stages
    )
    st.markdown(f'<div class="vst-grid">{cards}</div>', unsafe_allow_html=True)


def _format_seconds(value: float | None) -> str:
    if value is None:
        return "Unavailable"
    if value < 60:
        return f"{value:.0f} sec"
    if value < 3600:
        return f"{value / 60:.1f} min"
    return f"{value / 3600:.1f} hr"


def _resolve_project(value: str) -> Path | None:
    if not value.strip():
        return None
    path = Path(value).expanduser().resolve()
    return path if path.is_file() else None


def _save_upload(upload: Any, destination: Path) -> Path | None:
    if upload is None:
        return None
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(upload.getvalue())
    return destination


def _status_badge(status: str) -> str:
    normalized = status.lower().replace("_", " ")
    css = (
        "vst-completed"
        if status == "completed"
        else "vst-running"
        if status == "running"
        else "vst-required"
        if status in {"external_required", "failed"}
        else "vst-muted-status"
        if status in {"disabled", "unavailable"}
        else "vst-pending"
    )
    return f'<span class="vst-status {css}">{normalized}</span>'


STAGE_PRESENTATION = {
    "validate": ("Inputs & QC", "Validate project files, profiles, and declared boundaries."),
    "discover": (
        "Guide discovery & conservation",
        "Enumerate editor-compatible viral sites and measure strain support.",
    ),
    "host_screen": (
        "Host-risk search",
        "Run or await the declared Cas-OFFinder/host-reference search.",
    ),
    "pairs": ("Guide pairs", "Build same-gene and multi-target pair hypotheses."),
    "virtual_knockout": (
        "Predicted disruption",
        "Map bounded indel hypotheses to CDS, protein, and domain context.",
    ),
    "escape": (
        "Exact-target escape",
        "Enumerate sequence changes that remove an exact guide/PAM target.",
    ),
    "multiplex": (
        "Multiplex comparison",
        "Compare configured multi-guide strategies without claiming efficacy.",
    ),
    "report": (
        "Report & evidence assembly",
        "Assemble results and preserve evidence-review status and missingness.",
    ),
    "bundle": ("Shortlist & export", "Create START_HERE, tables, explanations, and ZIP."),
}


def _workflow_cards(
    statuses: dict[str, str],
    messages: dict[str, str] | None = None,
    estimates: dict[str, float | None] | None = None,
) -> str:
    messages = messages or {}
    estimates = estimates or {}
    cards = []
    for stage, (title, description) in STAGE_PRESENTATION.items():
        status = statuses.get(stage, "pending")
        runtime = estimates.get(stage)
        estimate = (
            f"Estimated: {_format_seconds(runtime)}" if runtime is not None else "ETA unavailable"
        )
        detail = messages.get(stage) or description
        cards.append(
            f"<div class='vst-stage {escape(status)}'>"
            f"<strong>{escape(title)}</strong>{_status_badge(status)}<br><br>"
            f"<span>{escape(detail)}</span><br><span><b>{escape(estimate)}</b></span></div>"
        )
    return "<div class='vst-grid'>" + "".join(cards) + "</div>"


def _remaining_estimate(
    statuses: dict[str, str], estimates: dict[str, float | None]
) -> float | None:
    terminal = {"completed", "external_required", "failed", "unavailable", "disabled"}
    remaining = [
        estimates.get(stage)
        for stage in STAGE_PRESENTATION
        if statuses.get(stage, "pending") not in terminal
    ]
    if not remaining or any(value is None for value in remaining):
        return 0.0 if not remaining else None
    return float(sum(value for value in remaining if value is not None))


def _render_overview(st: Any) -> None:
    _hero(st)
    st.subheader("From sequence data to a reviewable research decision")
    st.caption(
        "Every layer remains visible. Missing external searches or biological evidence are "
        "reported as missing, never converted into a favorable score."
    )
    _pipeline(st)
    left, right = st.columns(2)
    with left:
        st.markdown("### Start a new virus")
        st.markdown(
            "Create a self-contained project from local FASTA/GFF3 files or an NCBI "
            "accession, inspect the plan, then run or resume cached stages."
        )
        st.code(
            "vst init my-virus --interactive\n"
            "vst plan my-virus/project.yaml\n"
            "vst run my-virus/project.yaml",
            language="bash",
        )
    with right:
        st.markdown("### What you receive")
        st.markdown(
            "- One `START_HERE.html` report\n"
            "- Ranked guide and gene tables\n"
            "- A research shortlist and multiplex comparison\n"
            "- Human-review evidence queues\n"
            "- Timings, provenance, manifests, and a portable ZIP"
        )


def _render_create_project(st: Any) -> None:
    from viral_safe_target.researcher import create_demo_project, create_project

    st.header("Create a project")
    st.write(
        "Use the bundled demonstration for a two-minute orientation, or create a project "
        "from your own reference, annotation, and strain panel."
    )
    demo_col, note_col = st.columns([1, 2])
    with demo_col:
        demo_dir = st.text_input("Demo project directory", "demo-project")
        if st.button("Create bundled demo", type="primary", width="stretch"):
            try:
                project = create_demo_project(demo_dir, force=False)
                st.session_state["project_path"] = str(project)
                st.success(f"Created {project}")
            except Exception as error:  # noqa: BLE001 - actionable UI failure
                st.error(str(error))
    with note_col:
        st.info(
            "The demo is synthetic and small. It verifies installation and output structure; "
            "it is not biological evidence."
        )

    st.divider()
    st.subheader("New-virus project")
    source_mode = st.radio(
        "Reference source",
        ["Local files", "NCBI accession"],
        horizontal=True,
        help="Accession mode downloads public reference and annotation records from NCBI.",
    )
    with st.form("new_project"):
        a, b = st.columns(2)
        project_name = a.text_input("Project name", "my-virus")
        virus_name = b.text_input("Virus scientific name", "")
        output_dir = a.text_input("Project directory", "my-virus")
        tax_id = b.text_input("NCBI tax ID (optional)", "")
        reference_accession = ""
        reference_file = annotation_file = strains_file = host_file = None
        reference_local = annotation_local = strains_local = host_local = ""
        sequence_only = False
        if source_mode == "NCBI accession":
            reference_accession = st.text_input("NCBI reference accession")
            st.caption("Network access is used only when you submit this form.")
        else:
            reference_file = st.file_uploader(
                "Reference genome FASTA", type=["fa", "fasta", "fna"], key="reference"
            )
            annotation_file = st.file_uploader(
                "Matching annotation GFF3 (recommended)", type=["gff", "gff3"]
            )
            sequence_only = st.checkbox(
                "Continue in sequence-only mode without annotation",
                help="Gene, protein, domain, and evidence stages will be unavailable.",
            )
        strains_file = st.file_uploader(
            "Viral strain panel FASTA (optional)", type=["fa", "fasta", "fna"]
        )
        strains_aligned = st.checkbox("Strain panel is already aligned", value=True)
        host_file = st.file_uploader(
            "Local host FASTA (optional)", type=["fa", "fasta", "fna"], key="host"
        )
        with st.expander("Use existing local file paths"):
            st.caption(
                "Recommended for large assemblies that should not be copied through the browser. "
                "A path overrides the matching upload."
            )
            reference_local = st.text_input(
                "Reference FASTA path", disabled=source_mode != "Local files"
            )
            annotation_local = st.text_input(
                "Annotation GFF3 path", disabled=source_mode != "Local files"
            )
            strains_local = st.text_input("Strain panel FASTA path")
            host_local = st.text_input("Host FASTA path")
        evidence_enabled = st.checkbox("Prepare literature evidence discovery", value=False)
        submitted = st.form_submit_button("Create research project", type="primary")

    if not submitted:
        return
    try:
        with tempfile.TemporaryDirectory(prefix="vst-ui-") as temporary:
            staging = Path(temporary)
            reference_path = _save_upload(reference_file, staging / "reference.fasta")
            annotation_path = _save_upload(annotation_file, staging / "reference.gff3")
            strains_path = _save_upload(strains_file, staging / "strains.fasta")
            host_path = _save_upload(host_file, staging / "host.fasta")
            reference_source = (
                Path(reference_local).expanduser() if reference_local else reference_path
            )
            annotation_source = (
                Path(annotation_local).expanduser() if annotation_local else annotation_path
            )
            strains_source = Path(strains_local).expanduser() if strains_local else strains_path
            host_source = Path(host_local).expanduser() if host_local else host_path
            project = create_project(
                output_dir,
                project_name=project_name,
                virus_name=virus_name or project_name,
                tax_id=tax_id or None,
                reference_accession=reference_accession or None,
                reference_fasta=reference_source,
                annotation_gff=annotation_source,
                strains_fasta=strains_source,
                strains_aligned=strains_aligned,
                host_fasta=host_source,
                evidence_enabled=evidence_enabled,
                sequence_only=sequence_only,
            )
        st.session_state["project_path"] = str(project)
        st.success(f"Project created: {project}")
        st.info("Next: open Project workspace, review the plan, then run the workflow.")
    except Exception as error:  # noqa: BLE001 - actionable UI error
        st.error(f"Project creation stopped: {error}")


def _render_project_workspace(st: Any) -> None:
    from viral_safe_target.project_workflow import project_status, run_project
    from viral_safe_target.researcher import export_project, plan_project

    st.header("Project workspace")
    st.caption("Plan first, run only the stages you can support, and resume from valid caches.")
    default = st.session_state.get("project_path", "")
    project_value = st.text_input("Path to project.yaml", value=default)
    project = _resolve_project(project_value)
    if project is None:
        st.markdown(
            '<div class="vst-empty">Choose a valid <code>project.yaml</code>, or create a '
            "project from the Create project page.</div>",
            unsafe_allow_html=True,
        )
        return
    st.session_state["project_path"] = str(project)
    plan_col, run_col, resume_col, export_col = st.columns(4)
    show_plan = plan_col.button("Build plan", width="stretch")
    run = run_col.button("Run", type="primary", width="stretch")
    resume = resume_col.button("Resume", width="stretch")
    export = export_col.button("Export ZIP", width="stretch")
    external = st.toggle(
        "Run configured external host-search tools",
        value=False,
        help="Leave off when Cas-OFFinder/host references are unavailable. Missing stays explicit.",
    )
    current_status = project_status(project)
    live_statuses = {row["stage"]: row["status"] for row in current_status["stages"]}
    live_messages = {row["stage"]: row["message"] for row in current_status["stages"]}
    try:
        if show_plan or run or resume:
            st.session_state["project_plan"] = plan_project(project)
        plan = st.session_state.get("project_plan")
        estimates = (
            {row["stage"]: row.get("estimated_seconds") for row in plan["stages"]}
            if plan and plan.get("project") == str(project)
            else {}
        )
        st.subheader("Live workflow")
        current_stage = st.empty()
        progress = st.progress(
            sum(
                value
                in {"completed", "external_required", "failed", "unavailable", "disabled"}
                for value in live_statuses.values()
            )
            / len(STAGE_PRESENTATION),
            text="Ready to run or resume",
        )
        eta = st.empty()
        workflow_view = st.empty()
        workflow_view.markdown(
            _workflow_cards(live_statuses, live_messages, estimates), unsafe_allow_html=True
        )
        remaining = _remaining_estimate(live_statuses, estimates)
        eta.caption(
            "Best-effort remaining estimate: "
            + _format_seconds(remaining)
            + ". Hardware, external tools, and cache state can change this estimate."
        )
        if run or resume:
            started = time.perf_counter()

            def update_progress(event: dict[str, Any]) -> None:
                live_statuses.update(event["stages"])
                if event.get("message"):
                    live_messages[event["stage"]] = str(event["message"])
                title = STAGE_PRESENTATION.get(
                    event["stage"], (event["stage"].replace("_", " ").title(), "")
                )[0]
                current_stage.info(
                    f"Current stage: {title} — {event['status'].replace('_', ' ')}"
                )
                workflow_view.markdown(
                    _workflow_cards(live_statuses, live_messages, estimates),
                    unsafe_allow_html=True,
                )
                progress.progress(
                    min(1.0, float(event["progress_fraction"])),
                    text=(
                        f"Resolved {event['resolved_count']} of {event['total_count']} stages; "
                        f"elapsed {_format_seconds(time.perf_counter() - started)}"
                    ),
                )
                live_remaining = _remaining_estimate(live_statuses, estimates)
                eta.caption(
                    "Best-effort remaining estimate: "
                    + _format_seconds(live_remaining)
                    + ". A missing external tool resolves as external_required, not success."
                )

            with st.status("Running the project workflow…", expanded=True) as status_box:
                result = run_project(
                    project,
                    run_external=external,
                    restart=False,
                    progress_callback=update_progress,
                )
                status_box.write("Completed available stages and refreshed the result bundle.")
                status_box.update(label="Workflow finished", state="complete")
            st.session_state["project_result"] = result
            current_stage.success(
                "Workflow finished. Review completed and external-required stages below."
            )
        if export:
            bundle = export_project(project)
            st.success(f"Portable bundle created: {bundle}")
    except Exception as error:  # noqa: BLE001 - UI boundary
        st.error(str(error))

    plan = st.session_state.get("project_plan")
    if plan and plan.get("project") == str(project):
        c1, c2, c3 = st.columns(3)
        c1.metric("Estimated candidates", plan.get("candidate_count_estimate") or "Unavailable")
        disk = plan.get("estimated_disk_bytes")
        c2.metric("Estimated disk", f"{disk / 1024**2:.1f} MiB" if disk else "Unavailable")
        c3.metric("Estimated runtime", _format_seconds(plan.get("total_estimated_seconds")))
        stages = pd.DataFrame(plan["stages"])
        stages["estimated runtime"] = stages["estimated_seconds"].map(_format_seconds)
        st.dataframe(
            stages[["stage", "cached", "estimated runtime", "confidence", "reason"]],
            width="stretch",
            hide_index=True,
        )
        st.caption(plan["warning"])

    status = project_status(project)
    st.subheader("Recorded stage status")
    st.markdown(
        _workflow_cards(
            {row["stage"]: row["status"] for row in status["stages"]},
            {row["stage"]: row["message"] for row in status["stages"]},
            {
                row["stage"]: row.get("estimated_seconds")
                for row in (plan or {}).get("stages", [])
            },
        ),
        unsafe_allow_html=True,
    )
    st.markdown("### Evidence Agent (optional, separately reviewed)")
    evidence_queue = Path(status["output_root"]) / "evidence_review_queue.tsv"
    if evidence_queue.is_file() and evidence_queue.stat().st_size:
        st.success(
            "An evidence review queue exists. Agent proposals remain pending until a named "
            "human reviewer records a decision."
        )
    else:
        st.info(
            "Evidence discovery has not produced a review queue for this project. It is a "
            "separate, opt-in workflow using PubMed, Europe PMC, and UniProt; absence is not "
            "negative evidence."
        )
    report = Path(status["report"])
    if report.is_file():
        st.success(f"Researcher report is ready: {report}")


def _load_table(root: Path, name: str, *, separator: str = ",") -> pd.DataFrame:
    path = root / name
    return pd.read_csv(path, sep=separator) if path.is_file() else pd.DataFrame()


def _render_results(st: Any) -> None:
    st.header("Results explorer")
    st.caption("Review completed outputs without hiding unavailable stages or evidence gaps.")
    default_project = _resolve_project(st.session_state.get("project_path", ""))
    default_results = ""
    if default_project:
        try:
            from viral_safe_target.project_workflow import load_project

            default_results = str(load_project(default_project).output_root)
        except Exception:  # noqa: BLE001 - empty default is sufficient
            pass
    results_value = st.text_input("Results directory", value=default_results)
    root = Path(results_value).expanduser().resolve() if results_value.strip() else None
    if root is None or not root.is_dir():
        st.markdown(
            '<div class="vst-empty">Select a project results directory containing '
            "<code>summary.json</code> and the exported tables.</div>",
            unsafe_allow_html=True,
        )
        return
    summary_path = root / "summary.json"
    summary = json.loads(summary_path.read_text()) if summary_path.is_file() else {}
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Candidates", summary.get("candidate_count", "Unknown"))
    m2.metric("Shortlist", summary.get("shortlist_count", "Unknown"))
    m3.metric("Genes", summary.get("gene_count", "Unknown"))
    m4.metric("Host-risk", summary.get("host_risk_status", "Unknown"))

    guides = _load_table(root, "top_guides.csv")
    genes = _load_table(root, "top_genes.csv")
    shortlist = _load_table(root, "research_shortlist.csv")
    multiplex = _load_table(root, "multiplex_panels.csv")
    interventions = _load_table(root, "intervention_landscape.tsv", separator="\t")
    evidence = _load_table(root, "evidence_review_queue.tsv", separator="\t")
    tabs = st.tabs(
        ["Research shortlist", "Top genes", "Multiplex", "Evidence", "Files"]
    )
    with tabs[0]:
        if shortlist.empty:
            st.info("No shortlist is available yet.")
        else:
            st.dataframe(shortlist, width="stretch", hide_index=True, height=460)
            st.download_button(
                "Download research_shortlist.csv",
                shortlist.to_csv(index=False),
                "research_shortlist.csv",
                "text/csv",
            )
    with tabs[1]:
        data = genes if not genes.empty else guides
        if data.empty:
            st.info("No gene or guide ranking is available yet.")
        else:
            st.dataframe(data, width="stretch", hide_index=True, height=460)
    with tabs[2]:
        if multiplex.empty:
            st.info("Multiplex analysis is unavailable or pending.")
        else:
            st.dataframe(multiplex, width="stretch", hide_index=True)
    with tabs[3]:
        if not interventions.empty:
            st.markdown("#### Target-to-intervention landscape")
            st.dataframe(interventions, width="stretch", hide_index=True)
        if evidence.empty:
            st.info("No evidence proposals are awaiting human review.")
        else:
            st.markdown("#### Human review queue")
            st.dataframe(evidence, width="stretch", hide_index=True)
        st.caption("Pending or unknown evidence never increases a target ranking.")
    with tabs[4]:
        important = [
            "START_HERE.html",
            "SUMMARY.md",
            "summary.json",
            "top_guides.csv",
            "top_genes.csv",
            "research_shortlist.csv",
            "multiplex_panels.csv",
            "stage_timings.json",
            "run_manifest.json",
            "export.zip",
        ]
        rows = [
            {"file": name, "status": "available" if (root / name).is_file() else "missing"}
            for name in important
        ]
        st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)
        report = root / "START_HERE.html"
        if report.is_file():
            st.link_button("Open START_HERE report", report.as_uri())
        bundle = root / "export.zip"
        if bundle.is_file():
            st.download_button(
                "Download portable result bundle",
                bundle.read_bytes(),
                file_name=bundle.name,
                mime="application/zip",
            )


def _render_tool_hub(st: Any) -> None:
    from viral_safe_target.integrations.external_import import load_external_results
    from viral_safe_target.researcher import doctor_report, resource_path

    st.header("Tool & evidence hub")
    st.caption(
        "ViralSafeTarget orchestrates documented tools and exports. It does not silently "
        "replace their algorithms or scrape their public web interfaces."
    )
    project = _resolve_project(st.session_state.get("project_path", ""))
    detected = doctor_report(project).get("tools", {})
    rows = [
        {
            "system": "MAFFT",
            "role": "Viral strain alignment",
            "connection": "Native executable",
            "status": "available" if detected.get("mafft", {}).get("available") else "optional",
        },
        {
            "system": "Cas-OFFinder",
            "role": "Host-genome similarity search",
            "connection": "Native executable + parser",
            "status": "available"
            if detected.get("cas_offinder", {}).get("available")
            else "external required",
        },
        {
            "system": "CRISPRitz",
            "role": "Mismatch, bulge, and variant-aware host search",
            "connection": "Native/container execution + import",
            "status": "available"
            if detected.get("crispritz", {}).get("available")
            else "external required",
        },
        {
            "system": "CRISPOR",
            "role": "Independent guide scores and off-target metrics",
            "connection": "Researcher-supplied CSV/TSV/JSON export",
            "status": "import ready",
        },
        {
            "system": "CHOPCHOP",
            "role": "Independent guide ranking and design context",
            "connection": "Researcher-supplied CSV/TSV/JSON export",
            "status": "import ready",
        },
        {
            "system": "GuideScan2",
            "role": "Independent specificity and guide-database results",
            "connection": "Researcher-supplied CSV/TSV/JSON export",
            "status": "import ready",
        },
        {
            "system": "CRISPResso2",
            "role": "Measured editing outcomes from an existing result directory",
            "connection": "Measured-result import; kept separate from predictions",
            "status": "import ready",
        },
    ]
    st.dataframe(pd.DataFrame(rows), width="stretch", hide_index=True)

    st.markdown("### Import an independent design-tool export")
    st.write(
        "Match a documented export back to stable ViralSafeTarget candidate IDs. Raw values, "
        "unmatched rows, ambiguity, hashes, and missingness remain visible."
    )
    tool = st.selectbox("External tool", ["crispor", "chopchop", "guidescan2"])
    candidates_path = st.text_input(
        "Candidate table path",
        help="For example, results/discovery/discovery_panel.csv",
    )
    exported = st.file_uploader("CSV, TSV, or JSON export", type=["csv", "tsv", "txt", "json"])
    if st.button("Validate and normalize export", width="stretch"):
        candidate_file = Path(candidates_path).expanduser().resolve()
        if not candidate_file.is_file() or exported is None:
            st.error("Provide an existing candidate table and an external-tool export.")
        else:
            with tempfile.TemporaryDirectory(prefix="vst-import-") as temporary:
                suffix = Path(exported.name).suffix or ".csv"
                source = Path(temporary) / f"external{suffix}"
                source.write_bytes(exported.getvalue())
                result = load_external_results(
                    tool,
                    source,
                    resource_path("import_maps", f"{tool}.yaml"),
                    pd.read_csv(candidate_file),
                )
                normalized = result.results.dataframe
                c1, c2, c3 = st.columns(3)
                c1.metric("Normalized metrics", len(normalized))
                c2.metric("Unmatched rows", len(result.unmatched_rows))
                c3.metric("Ambiguous rows", len(result.ambiguous_rows))
                st.dataframe(normalized.head(100), width="stretch", hide_index=True)
                st.download_button(
                    "Download normalized long-format results",
                    normalized.to_csv(index=False).encode("utf-8"),
                    file_name=f"{tool}_normalized_results.csv",
                    mime="text/csv",
                )
                if not result.unmatched_rows.empty:
                    st.download_button(
                        "Download unmatched rows",
                        result.unmatched_rows.to_csv(index=False).encode("utf-8"),
                        file_name=f"{tool}_unmatched_rows.csv",
                        mime="text/csv",
                    )

    st.markdown("### Evidence Agent")
    evidence_rows = [
        ("PubMed / NCBI E-utilities", "Article discovery and PMID-linked metadata", "connected"),
        ("Europe PMC REST API", "Life-science publications, abstracts, and OA links", "connected"),
        (
            "UniProt REST API",
            "Protein function, domains, organism, and literature links",
            "connected",
        ),
        ("Human reviewer", "Approve, reject, or qualify every biological claim", "required"),
    ]
    st.dataframe(
        pd.DataFrame(evidence_rows, columns=["source", "role", "status"]),
        width="stretch",
        hide_index=True,
    )
    st.warning(
        "The agent discovers and structures candidate evidence. It cannot approve a claim, "
        "turn ortholog evidence into direct-virus evidence, or treat no search result as "
        "evidence of absence."
    )
    st.markdown(
        "See `docs/reference/TOOL_ECOSYSTEM_AND_INTEROPERABILITY.md` for official sources, "
        "exchange formats, current support, and the audited roadmap."
    )


def _render_environment(st: Any) -> None:
    from viral_safe_target.researcher import doctor_report

    st.header("Environment")
    st.caption("The local machine is inspected; no retired server or remote compute is assumed.")
    project = _resolve_project(st.session_state.get("project_path", ""))
    report = doctor_report(project)
    c1, c2, c3 = st.columns(3)
    c1.metric("ViralSafeTarget", report["viral_safe_target_version"])
    c2.metric("Python", report["python"]["version"])
    c3.metric("Disk free", f"{report['disk_free_gib']} GiB")
    tools = []
    for name, status in report["tools"].items():
        tools.append(
            {
                "tool": name.replace("_", " ").title(),
                "available": status["available"],
                "version": status["version"],
            }
        )
    st.dataframe(pd.DataFrame(tools), width="stretch", hide_index=True)
    st.info(report["scientific_boundary"])
    st.code('pip install "viral-safe-target[ui]"\nvst ui', language="bash")


def render_app() -> None:
    """Render the complete researcher interface."""
    st = _streamlit()
    st.set_page_config(
        page_title="ViralSafeTarget Research Workspace",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="expanded",
    )
    _inject_styles(st)
    with st.sidebar:
        st.markdown("## ViralSafeTarget")
        st.caption("Auditable viral target research")
        page = st.radio(
            "Workspace",
            [
                "Overview",
                "Create project",
                "Project workspace",
                "Results explorer",
                "Tool & evidence hub",
                "Environment",
            ],
            label_visibility="collapsed",
        )
        st.divider()
        st.caption(
            "Guide quality, host-risk, gene targetability, disruption, evidence, and escape "
            "remain separate decision axes."
        )
    if page == "Overview":
        _render_overview(st)
    elif page == "Create project":
        _render_create_project(st)
    elif page == "Project workspace":
        _render_project_workspace(st)
    elif page == "Results explorer":
        _render_results(st)
    elif page == "Tool & evidence hub":
        _render_tool_hub(st)
    else:
        _render_environment(st)


if __name__ == "__main__":  # Streamlit executes the module as a script.
    render_app()
