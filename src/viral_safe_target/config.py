"""Versioned configuration loading and validation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

_REPOSITORY_CONFIG = Path(__file__).resolve().parents[2] / "configs" / "research_v0.3.yaml"
_PACKAGED_CONFIG = Path(__file__).resolve().parent / "data" / "research_v0.3.yaml"
DEFAULT_CONFIG_PATH = _REPOSITORY_CONFIG if _REPOSITORY_CONFIG.is_file() else _PACKAGED_CONFIG


@dataclass(frozen=True)
class EditorProfile:
    """Sequence-search properties for one editor profile.

    Only the bundled SpCas9 profile is currently covered by scanner tests.
    """

    name: str
    protospacer_length: int
    pam_pattern: str
    pam_orientation: str
    cut_offset: int
    mismatch_search_threshold: int
    notes: str = ""
    tested: bool = False

    @classmethod
    def from_mapping(cls, values: dict[str, Any]) -> EditorProfile:
        required = {
            "name",
            "protospacer_length",
            "pam_pattern",
            "pam_orientation",
            "cut_offset",
            "mismatch_search_threshold",
        }
        missing = sorted(required - values.keys())
        if missing:
            raise ValueError(f"Editor profile is missing fields: {', '.join(missing)}")
        profile = cls(
            name=str(values["name"]),
            protospacer_length=int(values["protospacer_length"]),
            pam_pattern=str(values["pam_pattern"]).upper(),
            pam_orientation=str(values["pam_orientation"]),
            cut_offset=int(values["cut_offset"]),
            mismatch_search_threshold=int(values["mismatch_search_threshold"]),
            notes=str(values.get("notes", "")),
            tested=bool(values.get("tested", False)),
        )
        profile.validate()
        return profile

    def validate(self) -> None:
        if self.protospacer_length <= 0:
            raise ValueError("protospacer_length must be positive")
        if not self.pam_pattern or set(self.pam_pattern) - set("ACGTNRYSWKMBDHV"):
            raise ValueError(f"Invalid IUPAC PAM pattern: {self.pam_pattern!r}")
        if self.pam_orientation not in {"3prime", "5prime"}:
            raise ValueError("pam_orientation must be '3prime' or '5prime'")
        if not 0 <= self.cut_offset <= self.protospacer_length:
            raise ValueError("cut_offset must fall within the protospacer")
        if self.mismatch_search_threshold < 0:
            raise ValueError("mismatch_search_threshold cannot be negative")


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """Load and minimally validate the versioned research configuration."""
    config_path = Path(path) if path else DEFAULT_CONFIG_PATH
    with config_path.open(encoding="utf-8") as handle:
        values = yaml.safe_load(handle) or {}
    if values.get("schema_version") != "0.3":
        raise ValueError("Configuration schema_version must be '0.3'")
    for section in ("editor", "ranking", "pair_selection", "pair_scoring", "off_target"):
        if section not in values:
            raise ValueError(f"Configuration is missing section: {section}")
    EditorProfile.from_mapping(values["editor"])
    ranking = values["ranking"]
    weights = ranking.get("weights", {})
    required_weights = {
        "conservation",
        "viral_uniqueness",
        "gc",
        "sequence_complexity",
        "annotation",
        "gene_evidence",
    }
    if required_weights - weights.keys():
        raise ValueError("Ranking weights are incomplete")
    if sum(float(value) for value in weights.values()) <= 0:
        raise ValueError("At least one ranking weight must be positive")
    values["_config_path"] = str(config_path.resolve())
    return values


def get_editor(config: dict[str, Any]) -> EditorProfile:
    return EditorProfile.from_mapping(config["editor"])
