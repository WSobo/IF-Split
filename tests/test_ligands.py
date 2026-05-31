"""Phase 4 (early) tests: ligand classification + purification-artifact curation."""

from __future__ import annotations

from pathlib import Path

from ifsplit.config import Config, load_config
from ifsplit.ligands import (
    classify_components,
    elements_in_formula,
    has_histag,
    is_metal_ion,
    is_purification_artifact,
    longest_residue_run,
)
from ifsplit.schema import CandidateRecord, NonpolymerComp

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config" / "default.yaml"


def _cfg() -> Config:
    return load_config(DEFAULT_CONFIG)


def test_elements_in_formula():
    assert elements_in_formula("C34 H32 Fe N4 O4") == {"C", "H", "FE", "N", "O"}
    assert elements_in_formula("Zn") == {"ZN"}
    assert elements_in_formula("Ni 2+") == {"NI"}
    assert elements_in_formula(None) == set()


def test_is_metal_ion_distinguishes_ion_from_cofactor():
    zn = NonpolymerComp(comp_id="ZN", formula="Zn")
    hem = NonpolymerComp(comp_id="HEM", formula="C34 H32 Fe N4 O4")
    assert is_metal_ion(zn) is True
    assert is_metal_ion(hem) is False  # contains Fe but is an organic cofactor


def test_longest_residue_run_and_histag():
    assert longest_residue_run("AAAHHHHHHAAA", "H") == 6
    assert longest_residue_run("HAHAHA", "H") == 1
    assert has_histag("GSGSGHHHHHHGS", 6) is True
    assert has_histag("GSGSGHHHHHGS", 6) is False  # only 5


def test_zinc_finger_not_flagged_as_artifact(sample_entries):
    # 1A1F has Zn (a real, non-purification metal) -> never an artifact.
    rec = CandidateRecord.from_data_api(sample_entries["1A1F"])
    assert (
        is_purification_artifact(rec, purification_metals={"NI", "CO"}, histag_min_run=6) is False
    )


def test_histag_nickel_flagged_as_artifact(artifact_entry):
    rec = CandidateRecord.from_data_api(artifact_entry)
    assert is_purification_artifact(rec, purification_metals={"NI", "CO"}, histag_min_run=6) is True


def test_classify_drops_artifact_metal_by_default(artifact_entry):
    rec = CandidateRecord.from_data_api(artifact_entry)
    result = classify_components(rec, _cfg())
    assert result["purification_artifact"] is True
    # NI was the only metal and gets dropped -> no metal class.
    assert "metal" not in result["classes"]
    assert result["metals"] == []


def test_classify_keeps_artifact_metal_when_disabled(artifact_entry):
    rec = CandidateRecord.from_data_api(artifact_entry)
    cfg = _cfg().model_copy(update={"exclude_purification_artifacts": False})
    result = classify_components(rec, cfg)
    assert result["purification_artifact"] is True  # still flagged
    assert result["metals"] == ["NI"]  # but retained
    assert "metal" in result["classes"]


def test_classify_4hhb_small_molecule_and_no_metal(sample_entries):
    # HEM is small-molecule (organic cofactor); PO4 is a blacklisted additive.
    rec = CandidateRecord.from_data_api(sample_entries["4HHB"])
    result = classify_components(rec, _cfg())
    assert result["small_molecules"] == ["HEM"]
    assert result["metals"] == []
    assert result["has_nucleotide"] is False


def test_classify_1a1f_metal_and_nucleotide(sample_entries):
    rec = CandidateRecord.from_data_api(sample_entries["1A1F"])
    result = classify_components(rec, _cfg())
    assert "metal" in result["classes"]  # real Zn
    assert "nucleotide" in result["classes"]  # DNA chains
    assert result["metals"] == ["ZN"]
