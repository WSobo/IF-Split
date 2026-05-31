"""Phases 3/5/6/7 tests: filter, cluster, split, manifest, determinism (offline)."""

from __future__ import annotations

from pathlib import Path

from ifsplit.cluster import build_clusters
from ifsplit.config import Config, load_config
from ifsplit.dataset import load_dataset
from ifsplit.ligands import classify_components
from ifsplit.manifest import build_manifest, write_manifest
from ifsplit.parse import (
    DROP_NO_PROTEIN,
    DROP_TOO_LARGE,
    drop_summary,
    filter_candidates,
)
from ifsplit.schema import CandidateRecord
from ifsplit.split import assert_no_cluster_leakage, assign_splits, bucket, split_for_key

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config" / "default.yaml"


def _cfg(**over) -> Config:
    cfg = load_config(DEFAULT_CONFIG)
    return cfg.model_copy(update=over) if over else cfg


def _records(sample_entries, artifact_entry) -> list[CandidateRecord]:
    recs = [CandidateRecord.from_data_api(e) for e in sample_entries.values()]
    recs.append(CandidateRecord.from_data_api(artifact_entry))
    return recs


# ----------------------------- Stage 3: filter ----------------------------- #
def test_filter_keeps_protein_entries(sample_entries):
    recs = [CandidateRecord.from_data_api(e) for e in sample_entries.values()]
    kept, drops = filter_candidates(recs, _cfg())
    assert {r.entry_id for r in kept} == {"1A1F", "4HHB"}
    assert drops == []


def test_filter_drops_no_protein(sample_entries):
    # Strip 1A1F down to its DNA entities only.
    rec = CandidateRecord.from_data_api(sample_entries["1A1F"])
    rec = rec.model_copy(
        update={"polymer_entities": [e for e in rec.polymer_entities if e.is_nucleic]}
    )
    kept, drops = filter_candidates([rec], _cfg())
    assert kept == []
    assert drops[0]["reason"] == DROP_NO_PROTEIN


def test_filter_drops_too_large(sample_entries):
    kept, drops = filter_candidates(
        [CandidateRecord.from_data_api(sample_entries["4HHB"])],
        _cfg(max_total_residues=100),  # 4HHB assembly has 574
    )
    assert kept == []
    assert drops[0]["reason"] == DROP_TOO_LARGE
    assert drops[0]["residues"] == 574


def test_drop_summary_counts():
    drops = [
        {"entry_id": "A", "reason": DROP_NO_PROTEIN},
        {"entry_id": "B", "reason": DROP_NO_PROTEIN},
        {"entry_id": "C", "reason": DROP_TOO_LARGE},
    ]
    assert drop_summary(drops) == {DROP_NO_PROTEIN: 2, DROP_TOO_LARGE: 1}


# ---------------------------- Stage 5: cluster ----------------------------- #
def test_cluster_groups_by_membership(sample_entries, artifact_entry):
    recs = _records(sample_entries, artifact_entry)
    kept, _ = filter_candidates(recs, _cfg())
    cr = build_clusters(kept, _cfg())
    # 4HHB(two protein clusters)->longest chain; 1A1F; artifact: 3 distinct entries.
    assert set(cr.entry_to_cluster) == {"4HHB", "1A1F", "pdb_00009xyz"}
    # canonical keys are entity ids (smallest member), not raw integers.
    assert all(":" not in k or k.startswith("singleton:") for k in cr.cluster_members)


def test_cluster_multichain_detected(sample_entries):
    # 4HHB has two different protein clusters (alpha/beta) -> multichain.
    kept, _ = filter_candidates([CandidateRecord.from_data_api(sample_entries["4HHB"])], _cfg())
    cr = build_clusters(kept, _cfg())
    assert "4HHB" in cr.multichain_entries


# ------------------------------ Stage 6: split ----------------------------- #
def test_bucket_is_deterministic_and_unit_range():
    b = bucket("4HHB_1", "snapsplit-v1")
    assert b == bucket("4HHB_1", "snapsplit-v1")
    assert 0.0 <= b < 1.0
    assert bucket("4HHB_1", "other-salt") != b


def test_split_assignment_deterministic(sample_entries, artifact_entry):
    recs = _records(sample_entries, artifact_entry)
    kept, _ = filter_candidates(recs, _cfg())
    cr = build_clusters(kept, _cfg())
    a = assign_splits(cr, _cfg())
    b = assign_splits(cr, _cfg())
    assert a.cluster_split == b.cluster_split
    assert a.entry_split == b.entry_split


def test_no_cluster_leakage_invariant(sample_entries, artifact_entry):
    recs = _records(sample_entries, artifact_entry)
    kept, _ = filter_candidates(recs, _cfg())
    cr = build_clusters(kept, _cfg())
    res = assign_splits(cr, _cfg())
    assert_no_cluster_leakage(res, cr)  # raises on leakage


def test_registry_pins_assignment(sample_entries, artifact_entry):
    recs = _records(sample_entries, artifact_entry)
    kept, _ = filter_candidates(recs, _cfg())
    cr = build_clusters(kept, _cfg())
    # Force every cluster to "test" via a registry; hash is overridden.
    reg = {k: "test" for k in cr.cluster_members}
    res = assign_splits(cr, _cfg(), registry=reg)
    assert set(res.cluster_split.values()) == {"test"}


def test_fractions_roughly_respected_on_many_keys():
    # Synthetic: hash 3000 distinct keys, check broad proportions hold.
    cfg = _cfg()
    buckets = [split_for_key(f"K{i}", cfg) for i in range(3000)]
    train = buckets.count("train") / len(buckets)
    assert 0.74 < train < 0.86  # ~0.80 with sampling slack


# --------------------- Stages 6/7: manifest + loader ----------------------- #
def _full_manifest(sample_entries, artifact_entry, cfg):
    recs = _records(sample_entries, artifact_entry)
    kept, drops = filter_candidates(recs, cfg)
    class_map = {r.entry_id: classify_components(r, cfg) for r in kept}
    cr = build_clusters(kept, cfg)
    sp = assign_splits(cr, cfg)
    return build_manifest(
        cfg,
        candidates_sha256="deadbeef",
        n_candidates=len(recs),
        drops=drops,
        drop_counts=drop_summary(drops),
        clusters=cr,
        splits=sp,
        class_map=class_map,
    )


def test_manifest_is_deterministic(sample_entries, artifact_entry):
    cfg = _cfg()
    import json

    m1 = json.dumps(_full_manifest(sample_entries, artifact_entry, cfg), sort_keys=True)
    m2 = json.dumps(_full_manifest(sample_entries, artifact_entry, cfg), sort_keys=True)
    assert m1 == m2  # no wall-clock fields -> byte-identical


def test_manifest_has_all_entries(sample_entries, artifact_entry):
    m = _full_manifest(sample_entries, artifact_entry, _cfg())
    entries = m["splits"]["entries"]
    total = sum(len(v) for v in entries.values())
    assert total == 3  # 4HHB, 1A1F, pdb_00009xyz


def test_loader_roundtrip(tmp_path, sample_entries, artifact_entry):
    m = _full_manifest(sample_entries, artifact_entry, _cfg())
    path = write_manifest(m, tmp_path)
    ds = load_dataset(path)
    total = len(ds.train) + len(ds.val) + len(ds.test)
    assert total == 3
    assert ds.config_hash == _cfg().config_hash()


# ----------------------------- Phase 7: growth ----------------------------- #
def test_existing_cluster_does_not_move_when_dataset_grows(sample_entries, artifact_entry):
    cfg = _cfg()
    # Snapshot A: just the two sample entries.
    recs_a = [CandidateRecord.from_data_api(e) for e in sample_entries.values()]
    kept_a, _ = filter_candidates(recs_a, cfg)
    cr_a = build_clusters(kept_a, cfg)
    sp_a = assign_splits(cr_a, cfg)

    # Snapshot B: adds a third entry (growth).
    recs_b = _records(sample_entries, artifact_entry)
    kept_b, _ = filter_candidates(recs_b, cfg)
    cr_b = build_clusters(kept_b, cfg)
    sp_b = assign_splits(cr_b, cfg, registry=sp_a.cluster_split)

    # Every cluster present in A keeps its split in B.
    for key, split in sp_a.cluster_split.items():
        assert sp_b.cluster_split[key] == split
