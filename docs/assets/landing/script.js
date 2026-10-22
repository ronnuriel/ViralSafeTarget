const workflowDetails = [
  {
    number: "01",
    title: "Build the viral population model",
    copy: "Validate the reference, strain panel, annotation, host profile, and editor. Record accessions, versions, checksums, and rejected genomes.",
    items: ["FASTA and GFF3 validation", "Reference-coordinate control", "Input provenance and QC"],
    nodes: ["INPUT", "QC"],
  },
  {
    number: "02",
    title: "Enumerate every compatible guide site",
    copy: "Scan the complete viral reference for editor-compatible guide and PAM combinations instead of starting from a preselected gene.",
    items: ["Genome-wide enumeration", "Deterministic guide IDs", "Explicit rejection reasons"],
    nodes: ["GENOME", "GUIDES"],
  },
  {
    number: "03",
    title: "Measure population support and biological context",
    copy: "Quantify exact guide/PAM conservation across strains, then preserve one-to-many mappings into genes, coding sequences, proteins, and domains.",
    items: ["Multi-strain conservation", "Overlapping annotations", "Protein/domain context"],
    nodes: ["STRAINS", "CONTEXT"],
  },
  {
    number: "04",
    title: "Screen the declared host model",
    copy: "Connect Cas-OFFinder or CRISPRitz outputs to stable candidate identities. Missing tools remain external_required—never zero risk.",
    items: ["Mismatch-bounded search", "Bulge/variant extensions", "Raw hit provenance"],
    nodes: ["GUIDES", "HOST"],
  },
  {
    number: "05",
    title: "Model bounded disruption and exact-target loss",
    copy: "Generate clearly labelled sequence hypotheses for coding disruption and enumerate single-nucleotide counterfactuals that remove exact guide/PAM targets.",
    items: ["Size-defined indels", "CDS/protein mapping", "Multiplex escape barrier"],
    nodes: ["CUT", "OUTCOME"],
  },
  {
    number: "06",
    title: "Let agents propose and humans decide",
    copy: "Retrieve structured records from PubMed, Europe PMC, and UniProt. Separate direct-virus, ortholog, pending, and unknown evidence before human review.",
    items: ["Documented APIs", "Source-linked proposals", "Named reviewer required"],
    nodes: ["SOURCES", "REVIEW"],
  },
  {
    number: "07",
    title: "Package an auditable research shortlist",
    copy: "Deliver guide and gene tables, multiplex comparisons, explanation pages, logs, timings, manifests, and one portable collaborator bundle.",
    items: ["START_HERE report", "Machine-readable tables", "Portable export ZIP"],
    nodes: ["DECISION", "EXPORT"],
  },
];

const header = document.querySelector("[data-header]");
const menuButton = document.querySelector(".menu-button");
const navigation = document.querySelector(".site-nav");
const workflowSteps = [...document.querySelectorAll(".workflow-step")];
const detailNumber = document.querySelector("[data-detail-number]");
const detailTitle = document.querySelector("[data-detail-title]");
const detailCopy = document.querySelector("[data-detail-copy]");
const detailList = document.querySelector("[data-detail-list]");
const flowNodes = [...document.querySelectorAll(".flow-node")];
const reducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;

function updateHeader() {
  header?.classList.toggle("scrolled", window.scrollY > 24);
}

function activateStep(index) {
  const detail = workflowDetails[index];
  if (!detail) return;
  workflowSteps.forEach((step, stepIndex) => {
    const active = stepIndex === index;
    step.classList.toggle("active", active);
    step.setAttribute("aria-selected", String(active));
  });
  detailNumber.textContent = detail.number;
  detailTitle.textContent = detail.title;
  detailCopy.textContent = detail.copy;
  detailList.replaceChildren(...detail.items.map((item) => {
    const listItem = document.createElement("li");
    listItem.textContent = item;
    return listItem;
  }));
  flowNodes.forEach((node, nodeIndex) => { node.textContent = detail.nodes[nodeIndex]; });
}

workflowSteps.forEach((step, index) => {
  step.addEventListener("click", () => activateStep(index));
});

menuButton?.addEventListener("click", () => {
  const open = navigation.classList.toggle("open");
  menuButton.setAttribute("aria-expanded", String(open));
});

navigation?.querySelectorAll("a").forEach((link) => {
  link.addEventListener("click", () => {
    navigation.classList.remove("open");
    menuButton?.setAttribute("aria-expanded", "false");
  });
});

const revealObserver = new IntersectionObserver((entries) => {
  entries.forEach((entry) => {
    if (entry.isIntersecting) {
      entry.target.classList.add("visible");
      revealObserver.unobserve(entry.target);
    }
  });
}, { threshold: 0.12 });

document.querySelectorAll(".reveal").forEach((element) => revealObserver.observe(element));
window.addEventListener("scroll", updateHeader, { passive: true });
updateHeader();

if (!reducedMotion) {
  let workflowIndex = 0;
  window.setInterval(() => {
    const workflowSection = document.querySelector("#workflow");
    const rect = workflowSection?.getBoundingClientRect();
    if (rect && rect.top < window.innerHeight * .45 && rect.bottom > window.innerHeight * .4) {
      workflowIndex = (workflowIndex + 1) % workflowDetails.length;
      activateStep(workflowIndex);
    }
  }, 5000);
}
