from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def write_html_report(
    candidates: pd.DataFrame,
    output_path: str | Path,
    title: str = "ViralSafeTarget report",
) -> Path:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    preferred_columns = [
        "candidate_id",
        "gene_name",
        "product",
        "strand",
        "guide_sequence",
        "pam",
        "virus_site_coverage",
        "host_exact_matches",
        "host_min_mismatches",
        "demo_score",
        "pre_human_score",
        "decision",
    ]
    display_columns = [column for column in preferred_columns if column in candidates.columns]
    if display_columns:
        table = candidates[display_columns].to_html(
            index=False,
            escape=True,
            float_format=lambda value: f"{value:.3f}",
        )
    else:
        table = "<p>No candidates were produced.</p>"
    timestamp = datetime.now(timezone.utc).isoformat()
    notice = (
        "Computational research only: scores do not establish editing, viral inactivation, "
        "safety, or clinical efficacy."
    )
    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
body {{ font-family: Arial, sans-serif; max-width: 1400px; margin: 30px auto; line-height: 1.5; }}
table {{ border-collapse: collapse; width: 100%; display: block; overflow-x: auto; }}
th, td {{ border: 1px solid #ddd; padding: 7px; text-align: left; white-space: nowrap; }}
th {{ background: #f2f2f2; }}
.notice {{ padding: 12px; border: 1px solid #a76; background: #fff8e8; }}
</style>
</head>
<body>
<h1>{title}</h1>
<p>Generated: {timestamp}</p>
<div class="notice"><strong>Scope:</strong> {notice}</div>
<h2>Candidates</h2>
{table}
</body>
</html>"""
    output_path.write_text(html, encoding="utf-8")
    return output_path
