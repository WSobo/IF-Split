"""Phase 2 tests: candidate parsing + canonical serialization (offline)."""

from __future__ import annotations

from copy import deepcopy

from ifsplit.schema import (
    CandidateRecord,
    PolymerEntity,
    canonical_jsonl_bytes,
    metal_symbols_in_annotation,
    read_candidates_jsonl,
    sha256_hex,
    structural_families_from_instances,
)


def _pe(polymer_type: str) -> PolymerEntity:
    return PolymerEntity(entity_id="X_1", polymer_type=polymer_type, seq_len=4, seq="ACGU")


def test_is_nucleic_covers_dna_rna_and_hybrid():
    assert _pe("DNA").is_nucleic is True
    assert _pe("RNA").is_nucleic is True
    # rcsb_entity_polymer_type "NA-hybrid" (a single mixed DNA/RNA strand) contains
    # neither "DNA" nor "RNA" — it must still count as nucleic, or a protein/hybrid
    # complex silently loses its nucleic_acid class.
    assert _pe("NA-hybrid").is_nucleic is True
    assert _pe("Protein").is_nucleic is False
    assert _pe("Other").is_nucleic is False
    assert _pe("Protein").is_protein is True
    assert _pe("NA-hybrid").is_protein is False


def test_structural_families_from_instances():
    instances = [
        {
            "rcsb_polymer_instance_annotation": [
                {"type": "CATH", "annotation_id": "1.10.490.10", "name": "Globins"},
                {"type": "ECOD", "annotation_id": "e101mA1", "name": "Bcl-2"},
                {"type": "SCOP2", "annotation_id": "8035604", "name": "Globin-like"},
                {"type": "Pfam", "annotation_id": "PF00042", "name": "Globin"},  # ignored
            ]
        }
    ]
    fams = structural_families_from_instances(instances)
    assert fams["cath"] == ["1.10.490.10"]  # CATH keyed on the superfamily code
    assert fams["ecod"] == ["Bcl-2"]  # ECOD keyed on the family name
    assert fams["scop2"] == ["Globin-like"]  # SCOP2 keyed on the family name
    assert "pfam" not in fams  # non-structural annotations are not grouped on
    assert structural_families_from_instances(None) == {}
    assert structural_families_from_instances([]) == {}


def test_structural_families_roundtrip_through_data_api(sample_entries):
    # A parsed record carries per-method structural families when instances have them.
    entry = deepcopy(sample_entries["4HHB"])
    for pe in entry["polymer_entities"]:
        pe["polymer_entity_instances"] = [
            {"rcsb_polymer_instance_annotation": [{"type": "CATH", "annotation_id": "1.10.490.10"}]}
        ]
    rec = CandidateRecord.from_data_api(entry)
    assert rec.polymer_entities[0].structural_families.get("cath") == ["1.10.490.10"]
    # Deterministic serialization includes the new field.
    assert "structural_families" in rec.to_canonical_json()


def test_from_data_api_core_fields(sample_entries):
    rec = CandidateRecord.from_data_api(sample_entries["4HHB"])
    assert rec.entry_id == "4HHB"
    assert rec.methods == ["X-RAY DIFFRACTION"]
    assert rec.resolution_A == 1.74
    assert rec.release_date == "1984-07-17"
    assert rec.deposited_residues == 574
    assert rec.assemblies == {"4HHB-1": 574}


def test_polymer_entities_sorted_and_whitespace_stripped(sample_entries):
    rec = CandidateRecord.from_data_api(sample_entries["4HHB"])
    ids = [e.entity_id for e in rec.polymer_entities]
    assert ids == ["4HHB_1", "4HHB_2"]  # sorted regardless of input order
    e2 = next(e for e in rec.polymer_entities if e.entity_id == "4HHB_2")
    assert "\n" not in e2.seq
    assert e2.seq == "VHLTPEEKSAVTALWGKVNVDEVGGEALGR"
    assert e2.seq_len == len(e2.seq)


def test_cluster_membership_extracted(sample_entries):
    rec = CandidateRecord.from_data_api(sample_entries["4HHB"])
    e1 = next(e for e in rec.polymer_entities if e.entity_id == "4HHB_1")
    # _cluster_membership(49) -> identity 30 maps to cluster_id 49.
    assert e1.cluster_ids[30] == 49
    assert e1.cluster_ids[100] == 53
    assert set(e1.cluster_ids) == {30, 50, 70, 90, 95, 100}


def test_nucleic_entity_has_no_cluster_ids(sample_entries):
    rec = CandidateRecord.from_data_api(sample_entries["1A1F"])
    dna = next(e for e in rec.polymer_entities if e.is_nucleic)
    assert dna.cluster_ids == {}


def test_metal_symbols_in_annotation():
    assert metal_symbols_in_annotation("nickel cation binding") == {"NI"}
    assert metal_symbols_in_annotation("Urease nickel binding site") == {"NI"}
    assert metal_symbols_in_annotation("cobalt ion binding") == {"CO"}
    assert metal_symbols_in_annotation("zinc ion binding") == {"ZN"}
    # Heavy / lanthanide metals now covered so a native site is rescued, not demoted.
    assert metal_symbols_in_annotation("mercuric reductase activity") == {"HG"}
    assert metal_symbols_in_annotation("lanthanum-dependent methanol dehydrogenase") == {"LA"}
    assert metal_symbols_in_annotation("cerium ion binding") == {"CE"}
    # Whole-word match: "ytterbium"/"terbium" must NOT also fire the substrings
    # "terbium"/"erbium" (they are substrings of one another).
    assert metal_symbols_in_annotation("ytterbium binding protein") == {"YB"}
    assert metal_symbols_in_annotation("terbium luminescence site") == {"TB"}
    # Generic sentinel still uses a substring test (so "metallopeptidase" counts).
    assert metal_symbols_in_annotation("metallopeptidase activity") == {"METAL"}
    # A generic metal term with no named element -> the METAL sentinel.
    assert metal_symbols_in_annotation("transition metal ion binding") == {"METAL"}
    # No metal mentioned at all.
    assert metal_symbols_in_annotation("serine-type endopeptidase activity") == set()
    assert metal_symbols_in_annotation(None) == set()


def test_metal_annotations_extracted_from_data_api():
    entry = {
        "rcsb_id": "TEST",
        "polymer_entities": [
            {
                "rcsb_id": "TEST_1",
                "entity_poly": {
                    "rcsb_entity_polymer_type": "Protein",
                    "pdbx_seq_one_letter_code_can": "ACDE",
                },
                "rcsb_polymer_entity_annotation": [
                    {"type": "GO", "name": "nickel cation binding"},
                    {"type": "InterPro", "name": "Urease, alpha subunit"},
                ],
            }
        ],
    }
    rec = CandidateRecord.from_data_api(entry)
    assert rec.polymer_entities[0].metal_annotations == ["NI"]
    # A protein with no metal annotation gets an empty list (backward compatible).
    entry["polymer_entities"][0]["rcsb_polymer_entity_annotation"] = []
    assert CandidateRecord.from_data_api(entry).polymer_entities[0].metal_annotations == []


def test_nonpolymer_comps_sorted(sample_entries):
    rec = CandidateRecord.from_data_api(sample_entries["4HHB"])
    assert [c.comp_id for c in rec.nonpolymer_comps] == ["HEM", "PO4"]


def test_nucleic_acid_typed_in_metadata(sample_entries):
    # The "DNA-as-ATOM-records" gotcha dissolves: it's a typed entity here.
    rec = CandidateRecord.from_data_api(sample_entries["1A1F"])
    types = {e.polymer_type for e in rec.polymer_entities}
    assert "DNA" in types and "Protein" in types


def test_extended_pdb_id_stored_verbatim(artifact_entry):
    # Extended ids (pdb_xxxxxxxx) are kept exactly as RCSB returns them —
    # not sliced, not case-folded.
    rec = CandidateRecord.from_data_api(artifact_entry)
    assert rec.entry_id == "pdb_00009xyz"
    assert rec.polymer_entities[0].entity_id == "pdb_00009xyz_1"
    assert "pdb_00009xyz-1" in rec.assemblies


def test_quality_metrics_parsed(sample_entries):
    hhb = CandidateRecord.from_data_api(sample_entries["4HHB"])
    assert hhb.quality.clashscore == 142.32
    assert hhb.quality.ramachandran_outlier_pct == 1.24
    assert hhb.quality.rfree is None  # 4HHB has no recomputed diffraction summary
    assert hhb.quality.has_report is True

    a1f = CandidateRecord.from_data_api(sample_entries["1A1F"])
    assert a1f.quality.rfree == 0.21
    assert a1f.quality.rsrz_outlier_pct == 2.5


def test_protein_na_interface_count_parsed(sample_entries):
    a1f = CandidateRecord.from_data_api(sample_entries["1A1F"])
    assert a1f.protein_na_interface_count == 1  # zinc-finger / DNA
    hhb = CandidateRecord.from_data_api(sample_entries["4HHB"])
    assert hhb.protein_na_interface_count == 0  # protein-only


def test_quality_metrics_in_canonical_bytes(sample_entries):
    # Metrics are serialized into candidates.jsonl so consumers can post-filter
    # (full read-back is covered by test_jsonl_roundtrip).
    rec = CandidateRecord.from_data_api(sample_entries["1A1F"])
    data = canonical_jsonl_bytes([rec])
    assert b'"clashscore":4.5' in data
    assert b'"rfree":0.21' in data


def test_canonical_jsonl_is_order_independent(sample_entries):
    a = CandidateRecord.from_data_api(sample_entries["4HHB"])
    b = CandidateRecord.from_data_api(sample_entries["1A1F"])
    assert canonical_jsonl_bytes([a, b]) == canonical_jsonl_bytes([b, a])


def test_canonical_jsonl_sorted_by_entry_id(sample_entries):
    a = CandidateRecord.from_data_api(sample_entries["4HHB"])
    b = CandidateRecord.from_data_api(sample_entries["1A1F"])
    data = canonical_jsonl_bytes([a, b]).decode()
    first_line = data.splitlines()[0]
    assert '"entry_id":"1A1F"' in first_line  # 1A1F sorts before 4HHB


def test_jsonl_roundtrip(tmp_path, sample_entries):
    recs = [CandidateRecord.from_data_api(e) for e in sample_entries.values()]
    data = canonical_jsonl_bytes(recs)
    path = tmp_path / "candidates.jsonl"
    path.write_bytes(data)
    back = read_candidates_jsonl(path)
    assert canonical_jsonl_bytes(back) == data
    assert sha256_hex(canonical_jsonl_bytes(back)) == sha256_hex(data)
    # cluster_ids survive the roundtrip (json stringifies int keys -> restored).
    e1 = next(e for e in back[0].polymer_entities if e.cluster_ids)
    assert all(isinstance(k, int) for k in e1.cluster_ids)
