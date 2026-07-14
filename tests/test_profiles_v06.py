from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from viral_safe_target.profiles import load_profile_bundle, validate_profile_bundle

ROOT = Path(__file__).resolve().parents[1]


def _bundle():
    return load_profile_bundle(
        ROOT / "configs/viruses/hsv2.yaml",
        ROOT / "configs/hosts/human_grch38.yaml",
        ROOT / "configs/nucleases/spcas9.yaml",
        project_root=ROOT,
    )


def test_hsv2_profiles_are_configuration_driven_and_source_linked() -> None:
    bundle = _bundle()
    checks = validate_profile_bundle(bundle)
    assert bundle.virus["id"] == "hsv2"
    assert bundle.host["assembly_accession"] == "GCF_000001405.40"
    assert bundle.editor.pam_pattern == "NGG"
    assert not checks["status"].eq("fail").any()
    evidence_check = checks.loc[checks["component"].eq("gene evidence schema")].iloc[0]
    assert "uncited_rows=0" in evidence_check["detail"]


def test_large_host_reference_can_remain_explicitly_external() -> None:
    checks = validate_profile_bundle(_bundle(), require_large_host_reference=False)
    host_status = checks.loc[checks["component"].eq("host reference"), "status"].iloc[0]
    assert host_status in {"pass", "external_pending"}


def test_missing_virus_inputs_are_pending_unless_strictly_required() -> None:
    bundle = _bundle()
    bundle.virus["reference_fasta"] = "data/raw/definitely_missing_reference.fasta"
    relaxed = validate_profile_bundle(bundle)
    strict = validate_profile_bundle(bundle, require_virus_inputs=True)
    component = "virus path: reference_fasta"
    assert relaxed.loc[relaxed["component"].eq(component), "status"].iloc[0] == "input_pending"
    assert strict.loc[strict["component"].eq(component), "status"].iloc[0] == "fail"


def test_profile_type_mismatch_is_rejected(tmp_path: Path) -> None:
    invalid = tmp_path / "not-a-virus.yaml"
    invalid.write_text(
        yaml.safe_dump(
            {
                "schema_version": "1.0",
                "profile_type": "host",
                "id": "wrong_type",
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="not a 'virus' profile"):
        load_profile_bundle(
            invalid,
            ROOT / "configs/hosts/human_grch38.yaml",
            ROOT / "configs/nucleases/spcas9.yaml",
            project_root=ROOT,
        )
