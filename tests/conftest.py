"""Shared test fixtures: sample Data-API entry objects + a fake RCSB client.

Includes cluster-membership, bound-component / binding-affinity curation signals,
and an extended-PDB-ID entry so tests cover both the legacy 4-char and the
``pdb_xxxxxxxx`` identifier forms.
"""

from __future__ import annotations

import pytest


def _cluster_membership(c30: int) -> list[dict]:
    # Mirror the live shape: one record per identity level.
    return [
        {"cluster_id": c30 + 4, "identity": 100},
        {"cluster_id": c30 + 3, "identity": 95},
        {"cluster_id": c30 + 2, "identity": 90},
        {"cluster_id": c30 + 1, "identity": 70},
        {"cluster_id": c30 + 5, "identity": 50},
        {"cluster_id": c30, "identity": 30},
    ]


def _raw_4hhb() -> dict:
    return {
        "rcsb_id": "4HHB",
        "exptl": [{"method": "X-RAY DIFFRACTION"}],
        "rcsb_entry_info": {
            "resolution_combined": [1.74],
            "deposited_polymer_monomer_count": 574,
            # HEM contacts the protein; the PO4 buffer does not (so PO4 -> artifact).
            "nonpolymer_bound_components": ["HEM"],
        },
        "rcsb_accession_info": {"initial_release_date": "1984-07-17T00:00:00Z"},
        "polymer_entities": [
            {
                "rcsb_id": "4HHB_2",
                "entity_poly": {
                    "rcsb_entity_polymer_type": "Protein",
                    # deliberately wrapped to test whitespace stripping
                    "pdbx_seq_one_letter_code_can": "VHLTPEEKSAVTALWGK\nVNVDEVGGEALGR",
                },
                "rcsb_cluster_membership": _cluster_membership(48),
            },
            {
                "rcsb_id": "4HHB_1",
                "entity_poly": {
                    "rcsb_entity_polymer_type": "Protein",
                    "pdbx_seq_one_letter_code_can": "VLSPADKTNVKAAWGKVGAHAGEYGAEALE",
                },
                "rcsb_cluster_membership": _cluster_membership(49),
            },
        ],
        "nonpolymer_entities": [
            {
                "nonpolymer_comp": {
                    "chem_comp": {
                        "id": "PO4",
                        "name": "PHOSPHATE ION",
                        "formula": "O4 P",
                        "type": "non-polymer",
                    }
                }
            },
            {
                "nonpolymer_comp": {
                    "chem_comp": {
                        "id": "HEM",
                        "name": "PROTOPORPHYRIN IX CONTAINING FE",
                        "formula": "C34 H32 Fe N4 O4",
                        "type": "non-polymer",
                    }
                }
            },
        ],
        "assemblies": [
            {
                "rcsb_id": "4HHB-1",
                "rcsb_assembly_info": {"polymer_monomer_count": 574},
                "interfaces": [{"rcsb_interface_info": {"polymer_composition": "Protein (only)"}}],
            }
        ],
        # Real 4HHB validation values: a 1984 X-ray entry with a dreadful
        # clashscore (142) and no recomputed diffraction summary.
        "pdbx_vrpt_summary_geometry": [
            {
                "clashscore": 142.32,
                "percent_ramachandran_outliers": 1.24,
                "percent_rotamer_outliers": 9.52,
            }
        ],
        "pdbx_vrpt_summary_diffraction": None,
        "pdbx_vrpt_summary_em": None,
    }


def _raw_1a1f() -> dict:
    return {
        "rcsb_id": "1A1F",
        "exptl": [{"method": "X-RAY DIFFRACTION"}],
        "rcsb_entry_info": {
            "resolution_combined": [2.1],
            "deposited_polymer_monomer_count": 112,
            # The zinc-finger Zn is a real, bound, structural metal.
            "nonpolymer_bound_components": ["ZN"],
        },
        "rcsb_accession_info": {"initial_release_date": "1998-06-10T00:00:00Z"},
        "polymer_entities": [
            {
                "rcsb_id": "1A1F_3",
                "entity_poly": {
                    "rcsb_entity_polymer_type": "Protein",
                    "pdbx_seq_one_letter_code_can": "MERPYACPVESCDRRFSDSSNLTRHIRIHT",
                },
                "rcsb_cluster_membership": _cluster_membership(1000),
            },
            {
                "rcsb_id": "1A1F_1",
                "entity_poly": {
                    "rcsb_entity_polymer_type": "DNA",
                    "pdbx_seq_one_letter_code_can": "AGCGTGGGACC",
                },
                # DNA entities have no cluster membership.
            },
        ],
        "nonpolymer_entities": [
            {
                "nonpolymer_comp": {
                    "chem_comp": {
                        "id": "ZN",
                        "name": "ZINC ION",
                        "formula": "Zn",
                        "type": "non-polymer",
                    }
                }
            }
        ],
        "assemblies": [
            {
                "rcsb_id": "1A1F-1",
                "rcsb_assembly_info": {"polymer_monomer_count": 112},
                # The zinc-finger protein interfaces the DNA -> verified holo NA.
                "interfaces": [{"rcsb_interface_info": {"polymer_composition": "Protein/NA"}}],
            }
        ],
        # A well-refined X-ray entry: good geometry + a diffraction summary.
        "pdbx_vrpt_summary_geometry": [
            {
                "clashscore": 4.5,
                "percent_ramachandran_outliers": 0.0,
                "percent_rotamer_outliers": 1.0,
            }
        ],
        "pdbx_vrpt_summary_diffraction": [{"DCC_Rfree": 0.21, "percent_RSRZ_outliers": 2.5}],
        "pdbx_vrpt_summary_em": None,
    }


def _raw_extended() -> dict:
    """A His-tagged, Ni-only entry under an extended PDB ID (purification artifact)."""
    return {
        "rcsb_id": "pdb_00009xyz",
        "exptl": [{"method": "X-RAY DIFFRACTION"}],
        "rcsb_entry_info": {
            "resolution_combined": [2.0],
            "deposited_polymer_monomer_count": 250,
            # The His-tag does coordinate the Ni, so it shows as bound -- yet the
            # purification-artifact rule must still demote it (precedence test).
            "nonpolymer_bound_components": ["NI"],
        },
        "rcsb_accession_info": {"initial_release_date": "2025-01-15T00:00:00Z"},
        "polymer_entities": [
            {
                "rcsb_id": "pdb_00009xyz_1",
                "entity_poly": {
                    "rcsb_entity_polymer_type": "Protein",
                    # C-terminal hexa-His tag
                    "pdbx_seq_one_letter_code_can": "MKTAYIAKQRQISFVKSHFSRQLEERLGHHHHHH",
                },
                "rcsb_cluster_membership": _cluster_membership(7777),
            }
        ],
        "nonpolymer_entities": [
            {
                "nonpolymer_comp": {
                    "chem_comp": {
                        "id": "NI",
                        "name": "NICKEL (II) ION",
                        "formula": "Ni",
                        "type": "non-polymer",
                    }
                }
            }
        ],
        "assemblies": [
            {"rcsb_id": "pdb_00009xyz-1", "rcsb_assembly_info": {"polymer_monomer_count": 250}}
        ],
        "pdbx_vrpt_summary_geometry": [
            {
                "clashscore": 8.0,
                "percent_ramachandran_outliers": 0.5,
                "percent_rotamer_outliers": 2.0,
            }
        ],
        "pdbx_vrpt_summary_diffraction": [{"DCC_Rfree": 0.25, "percent_RSRZ_outliers": 4.0}],
        "pdbx_vrpt_summary_em": None,
    }


@pytest.fixture
def sample_entries() -> dict[str, dict]:
    return {"4HHB": _raw_4hhb(), "1A1F": _raw_1a1f()}


@pytest.fixture
def artifact_entry() -> dict:
    return _raw_extended()


class FakeRcsbClient:
    """In-memory stand-in for RcsbClient (no network)."""

    def __init__(self, entries: dict[str, dict]) -> None:
        self._entries = entries

    def search_entry_ids(self, cfg, limit=None, *, progress=None):
        ids = sorted(self._entries)
        return ids[:limit] if limit is not None else ids

    def count_entries(self, cfg):
        return len(self._entries)

    def fetch_entries(self, ids):
        for i in ids:
            yield self._entries[i]

    def close(self):
        pass


@pytest.fixture
def fake_client(sample_entries):
    return FakeRcsbClient(sample_entries)
