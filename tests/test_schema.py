"""Phase 2 tests: candidate parsing + canonical serialization (offline)."""

from __future__ import annotations

from ifsplit.schema import (
    CandidateRecord,
    canonical_jsonl_bytes,
    read_candidates_jsonl,
    sha256_hex,
)


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
