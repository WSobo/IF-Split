"""Phase 4 tests: ligand tiering + classification + purification-artifact curation."""

from __future__ import annotations

from pathlib import Path

from ifsplit.config import Config, load_config
from ifsplit.ligands import (
    TIER_AMBIGUOUS,
    TIER_ARTIFACT,
    TIER_FUNCTIONAL,
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


def _record(comps, *, bound=None, affinity=None, seq="ACDEFGHIKLMNPQRSTVWY", ptype="Protein"):
    """Build a CandidateRecord from a crafted Data-API-shaped dict."""
    entry = {
        "rcsb_id": "TEST",
        "exptl": [{"method": "X-RAY DIFFRACTION"}],
        "rcsb_entry_info": {
            "resolution_combined": [2.0],
            "deposited_polymer_monomer_count": len(seq),
            "nonpolymer_bound_components": bound or [],
        },
        "rcsb_accession_info": {"initial_release_date": "2020-01-01T00:00:00Z"},
        "rcsb_binding_affinity": [{"comp_id": c} for c in (affinity or [])],
        "polymer_entities": [
            {
                "rcsb_id": "TEST_1",
                "entity_poly": {
                    "rcsb_entity_polymer_type": ptype,
                    "pdbx_seq_one_letter_code_can": seq,
                },
            }
        ],
        "nonpolymer_entities": [
            {"nonpolymer_comp": {"chem_comp": {"id": cid, "formula": f}}} for cid, f in comps
        ],
        "assemblies": [
            {"rcsb_id": "TEST-1", "rcsb_assembly_info": {"polymer_monomer_count": len(seq)}}
        ],
    }
    return CandidateRecord.from_data_api(entry)


# ----------------------------- low-level helpers --------------------------- #
def test_elements_in_formula():
    assert elements_in_formula("C34 H32 Fe N4 O4") == {"C", "H", "FE", "N", "O"}
    assert elements_in_formula("Zn") == {"ZN"}
    assert elements_in_formula("Ni 2+") == {"NI"}
    assert elements_in_formula(None) == set()


def test_is_metal_ion_distinguishes_ion_from_cofactor():
    assert is_metal_ion(NonpolymerComp(comp_id="ZN", formula="Zn")) is True
    assert is_metal_ion(NonpolymerComp(comp_id="HEM", formula="C34 H32 Fe N4 O4")) is False


def test_longest_residue_run_and_histag():
    assert longest_residue_run("AAAHHHHHHAAA", "H") == 6
    assert longest_residue_run("HAHAHA", "H") == 1
    assert has_histag("GSGSGHHHHHHGS", 6) is True
    assert has_histag("GSGSGHHHHHGS", 6) is False  # only 5


# ------------------------------- tiering ----------------------------------- #
def test_bound_ligand_is_functional():
    rec = _record([("STI", "C29 H31 N7 O")], bound=["STI"])
    res = classify_components(rec, _cfg())
    assert res["tiers"]["STI"]["tier"] == TIER_FUNCTIONAL
    assert res["small_molecules"] == ["STI"]
    assert "small_molecule" in res["classes"]


def test_unbound_ligand_is_ambiguous_not_functional():
    rec = _record([("STI", "C29 H31 N7 O")], bound=[])  # present but not contacting
    res = classify_components(rec, _cfg())
    assert res["tiers"]["STI"]["tier"] == TIER_AMBIGUOUS
    assert res["small_molecules"] == []  # not labelled
    assert "small_molecule" in res["ambiguous_classes"]
    assert "small_molecule" not in res["classes"]


def test_affinity_forces_functional_even_if_unbound_list_empty():
    rec = _record([("STI", "C29 H31 N7 O")], bound=[], affinity=["STI"])
    res = classify_components(rec, _cfg())
    assert res["tiers"]["STI"]["tier"] == TIER_FUNCTIONAL
    assert res["small_molecules"] == ["STI"]


def test_additive_is_artifact():
    rec = _record([("GOL", "C3 H8 O3")], bound=["GOL"])  # glycerol, even if "bound"
    res = classify_components(rec, _cfg())
    assert res["tiers"]["GOL"]["tier"] == TIER_ARTIFACT
    assert res["small_molecules"] == []


def test_counterion_metal_is_artifact():
    rec = _record([("NA", "Na")], bound=["NA"])
    res = classify_components(rec, _cfg())
    assert res["tiers"]["NA"]["tier"] == TIER_ARTIFACT
    assert res["metals"] == []


def test_unbound_metal_is_ambiguous():
    rec = _record([("MG", "Mg")], bound=[])
    res = classify_components(rec, _cfg())
    assert res["tiers"]["MG"]["tier"] == TIER_AMBIGUOUS
    assert res["metals"] == []
    assert "metal" in res["ambiguous_classes"]


def test_bound_metal_is_functional():
    rec = _record([("MG", "Mg")], bound=["MG"])
    res = classify_components(rec, _cfg())
    assert res["tiers"]["MG"]["tier"] == TIER_FUNCTIONAL
    assert res["metals"] == ["MG"]
    assert "metal" in res["classes"]


# --------------------- purification-artifact curation ---------------------- #
def test_zinc_finger_not_flagged_as_artifact(sample_entries):
    rec = CandidateRecord.from_data_api(sample_entries["1A1F"])
    assert (
        is_purification_artifact(rec, purification_metals={"NI", "CO"}, histag_min_run=6) is False
    )


def test_histag_nickel_flagged_as_artifact(artifact_entry):
    rec = CandidateRecord.from_data_api(artifact_entry)
    assert is_purification_artifact(rec, purification_metals={"NI", "CO"}, histag_min_run=6) is True


def test_classify_drops_artifact_metal_by_default(artifact_entry):
    rec = CandidateRecord.from_data_api(artifact_entry)
    res = classify_components(rec, _cfg())
    assert res["purification_artifact"] is True
    assert res["tiers"]["NI"]["tier"] == TIER_ARTIFACT
    assert res["tiers"]["NI"]["reason"] == "histag_metal"
    assert "metal" not in res["classes"]
    assert res["metals"] == []


def test_classify_keeps_artifact_metal_when_disabled(artifact_entry):
    rec = CandidateRecord.from_data_api(artifact_entry)
    cfg = _cfg().model_copy(update={"exclude_purification_artifacts": False})
    res = classify_components(rec, cfg)
    assert res["purification_artifact"] is True  # still flagged
    # With exclusion off and Ni bound by the His-tag, it counts as a metal again.
    assert res["metals"] == ["NI"]
    assert "metal" in res["classes"]


# --------------------------- sample-entry checks --------------------------- #
def test_classify_4hhb_small_molecule_and_no_metal(sample_entries):
    rec = CandidateRecord.from_data_api(sample_entries["4HHB"])
    res = classify_components(rec, _cfg())
    assert res["small_molecules"] == ["HEM"]  # bound cofactor
    assert res["tiers"]["PO4"]["tier"] == TIER_ARTIFACT  # buffer
    assert res["metals"] == []
    assert res["has_nucleic_acid"] is False


def test_classify_1a1f_metal_and_nucleic_acid(sample_entries):
    rec = CandidateRecord.from_data_api(sample_entries["1A1F"])
    res = classify_components(rec, _cfg())
    assert "metal" in res["classes"]  # real bound Zn
    assert "nucleic_acid" in res["classes"]  # DNA chains, verified protein/NA interface
    assert res["metals"] == ["ZN"]


# -------------------- nucleic_acid holo gate (interface) ------------------- #
def _protein_na_record(prot_na_interfaces: int | None) -> CandidateRecord:
    """A protein + DNA entry; optional protein<->NA interface count on assembly 1.

    ``None`` omits the field entirely (RCSB had no interface data); an int sets
    ``num_prot_na_interface_entities`` (0 = no protein/NA contact, >0 = holo).
    """
    info = {"polymer_monomer_count": 50}
    if prot_na_interfaces is not None:
        info["num_prot_na_interface_entities"] = prot_na_interfaces
    asm = {"rcsb_id": "TNA-1", "rcsb_assembly_info": info}
    return CandidateRecord.from_data_api(
        {
            "rcsb_id": "TNA",
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {
                "resolution_combined": [2.0],
                "deposited_polymer_monomer_count": 50,
            },
            "rcsb_accession_info": {"initial_release_date": "2020-01-01T00:00:00Z"},
            "polymer_entities": [
                {
                    "rcsb_id": "TNA_1",
                    "entity_poly": {
                        "rcsb_entity_polymer_type": "Protein",
                        "pdbx_seq_one_letter_code_can": "ACDEFGHIKLMNPQRSTVWY",
                    },
                },
                {
                    "rcsb_id": "TNA_2",
                    "entity_poly": {
                        "rcsb_entity_polymer_type": "DNA",
                        "pdbx_seq_one_letter_code_can": "ATGCATGC",
                    },
                },
            ],
            "nonpolymer_entities": [],
            "assemblies": [asm],
        }
    )


def test_nucleic_acid_functional_with_protein_na_interface():
    res = classify_components(_protein_na_record(2), _cfg())  # 2 protein/NA interfaces
    assert res["has_nucleic_acid"] is True
    assert "nucleic_acid" in res["classes"]
    assert res["tiers"]["nucleic_acid"]["tier"] == TIER_FUNCTIONAL


def test_nucleic_acid_ambiguous_without_interface():
    # DNA present but zero protein/NA interfaces (co-deposited, not holo).
    res = classify_components(_protein_na_record(0), _cfg())
    assert res["has_nucleic_acid"] is True
    assert "nucleic_acid" not in res["classes"]
    assert "nucleic_acid" in res["ambiguous_classes"]
    assert res["tiers"]["nucleic_acid"]["tier"] == TIER_AMBIGUOUS


def test_nucleic_acid_ambiguous_when_no_interface_data():
    res = classify_components(_protein_na_record(None), _cfg())
    assert "nucleic_acid" not in res["classes"]
    assert "nucleic_acid" in res["ambiguous_classes"]
