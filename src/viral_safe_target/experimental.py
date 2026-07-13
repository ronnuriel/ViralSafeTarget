"""Import measured CRISPResso2 outputs separately from predicted tool metrics."""

from __future__ import annotations

import json
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from .provenance import sha256_file

MEASURED_COLUMNS = [
    "candidate_id",
    "guide_sequence",
    "amplicon_name",
    "aligned_read_count",
    "modified_read_percentage",
    "insertion_percentage",
    "deletion_percentage",
    "substitution_percentage",
    "frameshift_percentage",
    "in_frame_percentage",
    "quantification_window",
    "measurement_type",
    "source_file",
    "source_file_sha256",
    "timestamp",
]


@dataclass(frozen=True)
class ExperimentalImportResult:
    measurements: pd.DataFrame
    metadata: dict[str, Any]

    def write(self, output_directory: str | Path) -> tuple[Path, Path]:
        output = Path(output_directory)
        output.mkdir(parents=True, exist_ok=True)
        measurements_path = output / "crispresso2_measured_metrics.csv"
        metadata_path = output / "crispresso2_import_manifest.json"
        self.measurements.to_csv(measurements_path, index=False)
        metadata_path.write_text(
            json.dumps(self.metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        return measurements_path, metadata_path


def _find_column(frame: pd.DataFrame, names: tuple[str, ...]) -> pd.Series:
    lowered = {str(column).lower(): column for column in frame.columns}
    for name in names:
        if name.lower() in lowered:
            return frame[lowered[name.lower()]]
    return pd.Series(pd.NA, index=frame.index)


def import_crispresso2_results(
    input_directory: str | Path, candidate_map: str | Path
) -> ExperimentalImportResult:
    """Read an existing CRISPResso2 directory; no experimental procedure is provided."""
    directory = Path(input_directory)
    mapping_path = Path(candidate_map)
    if not directory.is_dir():
        raise FileNotFoundError(f"CRISPResso2 output directory is missing: {directory}")
    candidates = pd.read_csv(mapping_path)
    quantification_files = sorted(directory.rglob("*quantification*.txt")) + sorted(
        directory.rglob("*Quantification*.txt")
    )
    if not quantification_files:
        quantification_files = sorted(directory.rglob("*.tsv"))
    rows = []
    imported_files = []
    timestamp = datetime.now(timezone.utc).isoformat()
    for path in dict.fromkeys(quantification_files):
        frame = pd.read_csv(path, sep=None, engine="python", comment="#")
        amplicon = _find_column(frame, ("Amplicon", "Amplicon_Name", "Name"))
        for index, _record in frame.iterrows():
            amplicon_name = (
                str(amplicon.loc[index]) if pd.notna(amplicon.loc[index]) else path.parent.name
            )
            matches = candidates
            if "amplicon_name" in candidates:
                matches = candidates[candidates["amplicon_name"].astype(str).eq(amplicon_name)]
            if len(matches) != 1:
                continue
            candidate = matches.iloc[0]
            rows.append(
                {
                    "candidate_id": candidate["candidate_id"],
                    "guide_sequence": candidate.get("guide_sequence", ""),
                    "amplicon_name": amplicon_name,
                    "aligned_read_count": _find_column(
                        frame, ("Reads_aligned", "Aligned_Reads")
                    ).loc[index],
                    "modified_read_percentage": _find_column(
                        frame, ("Modified%", "Modified_percentage")
                    ).loc[index],
                    "insertion_percentage": _find_column(
                        frame, ("Insertions%", "Insertion_percentage")
                    ).loc[index],
                    "deletion_percentage": _find_column(
                        frame, ("Deletions%", "Deletion_percentage")
                    ).loc[index],
                    "substitution_percentage": _find_column(
                        frame, ("Substitutions%", "Substitution_percentage")
                    ).loc[index],
                    "frameshift_percentage": _find_column(
                        frame, ("Frameshift%", "Frameshift_percentage")
                    ).loc[index],
                    "in_frame_percentage": _find_column(
                        frame, ("In-frame%", "In_frame_percentage")
                    ).loc[index],
                    "quantification_window": _find_column(
                        frame, ("Quantification_window", "QuantificationWindow")
                    ).loc[index],
                    "measurement_type": "measured experimental metric",
                    "source_file": str(path.resolve()),
                    "source_file_sha256": sha256_file(path),
                    "timestamp": timestamp,
                }
            )
        imported_files.append(
            {"path": str(path.resolve()), "sha256": sha256_file(path), "rows": len(frame)}
        )
    measurements = pd.DataFrame(rows, columns=MEASURED_COLUMNS)
    metadata = {
        "adapter": "CRISPResso2 import 0.4",
        "input_directory": str(directory.resolve()),
        "candidate_map": str(mapping_path.resolve()),
        "candidate_map_sha256": sha256_file(mapping_path),
        "imported_files": imported_files,
        "measurement_count": len(measurements),
        "separation_notice": "Measured experimental metrics are not prediction scores.",
        "timestamp": timestamp,
    }
    return ExperimentalImportResult(measurements, metadata)
