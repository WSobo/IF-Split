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


def _record(
    comps,
    *,
    bound=None,
    affinity=None,
    investigated=None,
    annotations=None,
    comp_types=None,
    seq="ACDEFGHIKLMNPQRSTVWY",
    ptype="Protein",
):
    """Build a CandidateRecord from a crafted Data-API-shaped dict.

    ``investigated`` = comp ids RCSB flagged subject-of-investigation (SOI='Y').
    ``annotations`` = GO/InterPro/Pfam annotation names on the protein entity
    (e.g. ["nickel cation binding"]).
    """
    investigated = set(investigated or [])
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
                "rcsb_polymer_entity_annotation": [
                    {"type": "GO", "name": n} for n in (annotations or [])
                ],
            }
        ],
        "nonpolymer_entities": [
            {
                "nonpolymer_comp": {
                    "chem_comp": {"id": cid, "formula": f, "type": (comp_types or {}).get(cid)}
                },
                "nonpolymer_entity_instances": [
                    {
                        "rcsb_nonpolymer_instance_validation_score": [
                            {"is_subject_of_investigation": "Y" if cid in investigated else "N"}
                        ]
                    }
                ],
            }
            for cid, f in comps
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


def test_terminal_histag_catches_partial_tag():
    # A C-terminal HHH (a 6xHis with 3 residues unmodeled/trimmed) is NOT caught
    # by the full-run rule, but IS caught by the terminal rule.
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGEFHHH"  # terminal run of 3
    assert has_histag(seq, 6) is False
    assert has_histag(seq, 6, terminal_min_run=3) is True
    # An *internal* HHH (a real metal motif, far from either end) is NOT caught.
    internal = "M" + "A" * 30 + "HHH" + "A" * 30 + "K"
    assert has_histag(internal, 6, terminal_min_run=3) is False


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


def test_investigated_cofactor_is_functional_even_if_unbound():
    # FAD-like case: non-covalently bound (absent from bound_components) but RCSB
    # flags it subject-of-investigation -> functional, not ambiguous.
    rec = _record([("FAD", "C27 H33 N9 O15 P2")], bound=[], investigated=["FAD"])
    res = classify_components(rec, _cfg())
    assert res["tiers"]["FAD"]["tier"] == TIER_FUNCTIONAL
    assert res["tiers"]["FAD"]["reason"] == "ligand_investigated"
    assert "small_molecule" in res["classes"]


def test_uninvestigated_unbound_ligand_stays_ambiguous():
    rec = _record([("STI", "C29 H31 N7 O")], bound=[], investigated=[])
    res = classify_components(rec, _cfg())
    assert res["tiers"]["STI"]["tier"] == TIER_AMBIGUOUS


def test_investigated_buffer_still_artifact():
    # A blacklisted additive is an artifact even if (oddly) flagged investigated:
    # the additive gate precedes the SOI signal.
    rec = _record([("GOL", "C3 H8 O3")], investigated=["GOL"])
    res = classify_components(rec, _cfg())
    assert res["tiers"]["GOL"]["tier"] == TIER_ARTIFACT


def test_investigated_nickel_corroborates_not_artifact():
    # A lone Ni flagged subject-of-investigation is corroborated -> functional,
    # not demoted to ambiguous.
    rec = _record([("NI", "Ni")], bound=["NI"], investigated=["NI"], seq="ACDEFGHIKLMNPQRSTVWY")
    res = classify_components(rec, _cfg())
    assert res["tiers"]["NI"]["tier"] == TIER_FUNCTIONAL
    assert res["metals"] == ["NI"]


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


def test_bound_halide_anion_is_counterion_not_ligand():
    # Regression: a bound halide (F/I) is a counterion, not a functional small
    # molecule. Previously the counterion check sat inside the metal branch, which
    # is False for anions, so a bound F leaked through as ligand_bound.
    for anion in ("F", "I", "IOD"):
        rec = _record([(anion, anion)], bound=[anion])
        res = classify_components(rec, _cfg())
        assert res["tiers"][anion]["tier"] == TIER_ARTIFACT
        assert res["tiers"][anion]["reason"] == "counterion"
        assert res["small_molecules"] == []


def test_bound_glycan_is_glycan_not_small_molecule():
    # A bound N-acetylglucosamine with no SOI/affinity is decorative glycosylation,
    # not a ligand pocket -> reported glycan, not a small-molecule target.
    rec = _record(
        [("NAG", "C8 H15 N O6")],
        bound=["NAG"],
        comp_types={"NAG": "D-saccharide, beta linking"},
    )
    res = classify_components(rec, _cfg())
    assert res["tiers"]["NAG"]["tier"] == TIER_AMBIGUOUS
    assert res["tiers"]["NAG"]["reason"] == "glycan"
    assert res["small_molecules"] == []
    assert "small_molecule" in res["ambiguous_classes"]


def test_investigated_glycan_is_still_glycan_soi_does_not_rescue():
    # RCSB's SOI flag is noisy for carbohydrates (it flags glycosylation + detergents),
    # so an SOI-only glycan is NOT rescued -> stays glycan, not a small-molecule target.
    rec = _record(
        [("MAN", "C6 H12 O6")],
        bound=["MAN"],
        investigated=["MAN"],
        comp_types={"MAN": "D-saccharide, alpha linking"},
    )
    res = classify_components(rec, _cfg())
    assert res["tiers"]["MAN"]["tier"] == TIER_AMBIGUOUS
    assert res["tiers"]["MAN"]["reason"] == "glycan"
    assert res["small_molecules"] == []


def test_glycan_with_affinity_is_functional():
    rec = _record(
        [("GAL", "C6 H12 O6")],
        bound=["GAL"],
        affinity=["GAL"],
        comp_types={"GAL": "D-saccharide, beta linking"},
    )
    res = classify_components(rec, _cfg())
    assert res["tiers"]["GAL"]["tier"] == TIER_FUNCTIONAL
    assert res["small_molecules"] == ["GAL"]


def test_non_saccharide_bound_ligand_still_functional():
    # A bound non-sugar ligand is unaffected by the glycan gate.
    rec = _record([("STI", "C29 H31 N7 O")], bound=["STI"], comp_types={"STI": "non-polymer"})
    res = classify_components(rec, _cfg())
    assert res["tiers"]["STI"]["tier"] == TIER_FUNCTIONAL
    assert res["tiers"]["STI"]["reason"] == "ligand_bound"


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


def test_lone_bound_nickel_no_tag_is_ambiguous_not_functional():
    # Ni is the only metal, bound, but no His-tag and no affinity: most such lone
    # Ni have their tag trimmed from the sequence -> demote to ambiguous.
    rec = _record([("NI", "Ni")], bound=["NI"], seq="ACDEFGHIKLMNPQRSTVWY")
    res = classify_components(rec, _cfg())
    assert res["tiers"]["NI"]["tier"] == TIER_AMBIGUOUS
    assert res["tiers"]["NI"]["reason"] == "purification_metal_uncorroborated"
    assert res["metals"] == []
    assert "metal" in res["ambiguous_classes"]


def test_lone_nickel_with_affinity_stays_functional():
    # A measured binding affinity corroborates a real Ni site -> functional.
    rec = _record([("NI", "Ni")], bound=["NI"], affinity=["NI"], seq="ACDEFGHIKLMNPQRSTVWY")
    res = classify_components(rec, _cfg())
    assert res["tiers"]["NI"]["tier"] == TIER_FUNCTIONAL
    assert res["metals"] == ["NI"]


def test_lone_nickel_with_matching_metal_annotation_is_rescued():
    # The protein is annotated (GO) as a nickel-binding enzyme (urease-like): a lone
    # Ni with no tag/affinity is rescued to functional, not demoted.
    rec = _record([("NI", "Ni")], bound=["NI"], annotations=["nickel cation binding"])
    res = classify_components(rec, _cfg())
    assert res["tiers"]["NI"]["tier"] == TIER_FUNCTIONAL
    assert res["tiers"]["NI"]["reason"] == "metal_annotated"
    assert res["metals"] == ["NI"]


def test_lone_nickel_nonnative_metal_annotation_is_ambiguous_nonnative():
    # A real metalloprotein whose annotated metal is Zn (Ni is a substitute): stays
    # ambiguous, but with the distinguishing metal_site_nonnative reason.
    rec = _record([("NI", "Ni")], bound=["NI"], annotations=["zinc ion binding"])
    res = classify_components(rec, _cfg())
    assert res["tiers"]["NI"]["tier"] == TIER_AMBIGUOUS
    assert res["tiers"]["NI"]["reason"] == "metal_site_nonnative"
    assert res["metals"] == []
    assert "metal" in res["ambiguous_classes"]


def test_lone_nickel_generic_metal_annotation_is_nonnative_not_uncorroborated():
    # A generic "transition metal ion binding" term marks a metalloprotein without
    # naming Ni -> metal_site_nonnative, distinct from the no-annotation case.
    rec = _record([("NI", "Ni")], bound=["NI"], annotations=["transition metal ion binding"])
    res = classify_components(rec, _cfg())
    assert res["tiers"]["NI"]["reason"] == "metal_site_nonnative"


def test_metal_annotation_overrides_histag_demotion():
    # A His-tagged construct of a real nickel enzyme: positive annotation beats the
    # purification heuristic (comp-level tiering can't split tag-Ni from catalytic-Ni,
    # so real metal biology present -> functional).
    tagged = "MGHHHHHHSSG" + "ACDEFGIKLMNPQRSTVWY" * 2
    rec = _record([("NI", "Ni")], bound=["NI"], annotations=["nickel cation binding"], seq=tagged)
    res = classify_components(rec, _cfg())
    assert res["purification_artifact"] is True  # a His-tag IS detected
    assert res["tiers"]["NI"]["tier"] == TIER_FUNCTIONAL  # but annotation overrides it
    assert res["tiers"]["NI"]["reason"] == "metal_annotated"


def test_nickel_with_real_metal_stays_functional():
    # Ni alongside a real biological metal (Zn) is not a lone-Ni case -> both
    # bound metals stay functional (entry_metals is not a subset of {NI, CO}).
    rec = _record([("NI", "Ni"), ("ZN", "Zn")], bound=["NI", "ZN"], seq="ACDEFGHIKLMNPQRSTVWY")
    res = classify_components(rec, _cfg())
    assert res["tiers"]["NI"]["tier"] == TIER_FUNCTIONAL
    assert res["tiers"]["ZN"]["tier"] == TIER_FUNCTIONAL
    assert set(res["metals"]) == {"NI", "ZN"}


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
