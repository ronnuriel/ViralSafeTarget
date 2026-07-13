"""Checksum-based stage caching helpers."""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path

from .provenance import sha256_file


def stage_signature(input_paths: Iterable[str | Path], parameters: dict[str, object]) -> dict:
    return {
        "inputs": [
            {"path": str(Path(path).resolve()), "sha256": sha256_file(path)} for path in input_paths
        ],
        "parameters": parameters,
    }


def stage_is_current(
    stamp_path: str | Path,
    output_paths: Iterable[str | Path],
    input_paths: Iterable[str | Path],
    parameters: dict[str, object],
) -> bool:
    stamp = Path(stamp_path)
    outputs = [Path(path) for path in output_paths]
    if (
        not stamp.is_file()
        or not outputs
        or any(not path.is_file() or path.stat().st_size == 0 for path in outputs)
    ):
        return False
    try:
        stored = json.loads(stamp.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return stored == stage_signature(input_paths, parameters)


def write_stage_stamp(
    stamp_path: str | Path,
    input_paths: Iterable[str | Path],
    parameters: dict[str, object],
) -> Path:
    stamp = Path(stamp_path)
    stamp.parent.mkdir(parents=True, exist_ok=True)
    stamp.write_text(
        json.dumps(stage_signature(input_paths, parameters), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return stamp
