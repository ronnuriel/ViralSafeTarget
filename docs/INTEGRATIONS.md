# External-tool integrations

ViralSafeTarget v0.4 integrates tools through adapters; it does not reimplement their search,
alignment, scoring, or sequencing-analysis algorithms.

| Tool | Support | Role |
|---|---|---|
| MAFFT | executable adapter | Alignment command, version, input and output provenance |
| Cas-OFFinder | executable and native-output parser | Reference-genome predicted-hit search |
| CRISPRitz | executable, Docker-command, input-only and import adapter | Mismatches, DNA/RNA bulges, variants and population metadata |
| CRISPOR | documented export import | Known exported metrics only |
| CHOPCHOP | documented export import | Known exported metrics only |
| GuideScan2 | documented export import | Known exported metrics only |
| CRISPResso2 | measured-result directory import | Experimental sequencing metrics, kept separate |

Run `vst tools doctor` for availability. An unavailable executable is a pending stage, never a
zero-risk result. Adapters preserve source files, hashes, versions, command/import method and errors.
Do not scrape public web interfaces or rely on undocumented endpoints.

To add a tool, implement `ToolAdapter`, emit the normalized long schema, document metric direction,
use synthetic fixtures, and make missing output explicit.
