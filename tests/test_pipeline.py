"""Phases 3/5/6/7 tests: filter, cluster, split, manifest, determinism (offline)."""

from __future__ import annotations

from pathlib import Path

from ifsplit.cluster import build_clusters
from ifsplit.config import Config, load_config
from ifsplit.dataset import load_dataset
from ifsplit.ligands import classify_components
from ifsplit.manifest import build_manifest, write_manifest
from ifsplit.parse import (
    DROP_CLASHSCORE,
    DROP_NO_PROTEIN,
    DROP_NO_VALIDATION,
    DROP_RFREE,
    DROP_TOO_LARGE,
    drop_summary,
    filter_candidates,
)
from ifsplit.schema import CandidateRecord, PolymerEntity
from ifsplit.split import assign_splits, bucket, check_no_leakage, split_for_key

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


def test_filter_drops_high_clashscore(sample_entries):
    # 4HHB's real clashscore is 142; a 40 cap drops it but keeps 1A1F (4.5).
    recs = [CandidateRecord.from_data_api(e) for e in sample_entries.values()]
    kept, drops = filter_candidates(recs, _cfg(max_clashscore=40.0))
    assert {r.entry_id for r in kept} == {"1A1F"}
    assert drops[0]["reason"] == DROP_CLASHSCORE
    assert drops[0]["value"] == 142.32


def test_filter_keeps_when_metric_absent(sample_entries):
    # 4HHB has no diffraction summary -> rfree is None -> an rfree cap can't drop it.
    rec = CandidateRecord.from_data_api(sample_entries["4HHB"])
    assert rec.quality.rfree is None
    kept, drops = filter_candidates([rec], _cfg(max_rfree=0.25))
    assert [r.entry_id for r in kept] == ["4HHB"]
    assert drops == []


def test_filter_drops_high_rfree(sample_entries):
    # 1A1F has DCC_Rfree 0.21; a 0.20 cap drops it.
    rec = CandidateRecord.from_data_api(sample_entries["1A1F"])
    kept, drops = filter_candidates([rec], _cfg(max_rfree=0.20))
    assert kept == []
    assert drops[0]["reason"] == DROP_RFREE


def test_require_validation_report_drops_reportless_entry():
    # _protein_record() builds a record with no validation summary at all.
    rec = _protein_record("AAA1", [10])
    assert rec.quality.has_report is False
    kept, drops = filter_candidates([rec], _cfg(require_validation_report=True))
    assert kept == []
    assert drops[0]["reason"] == DROP_NO_VALIDATION


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
    # 4HHB (two protein clusters merge into one component); 1A1F; artifact: 3 entries.
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
    check_no_leakage(res, cr)  # raises on leakage


def test_registry_pins_assignment(sample_entries, artifact_entry):
    recs = _records(sample_entries, artifact_entry)
    kept, _ = filter_candidates(recs, _cfg())
    cr = build_clusters(kept, _cfg())
    # Force every cluster to "test" via a registry; hash is overridden.
    reg = {k: "test" for k in cr.cluster_members}
    res = assign_splits(cr, _cfg(), registry=reg)
    assert set(res.cluster_split.values()) == {"test"}


def test_test_minimums_recruit_components_no_leakage(sample_entries, artifact_entry):
    # 1A1F carries a functional metal (bound Zn). With the pure hash it may not be
    # in test; a metal floor of 1 must pull its whole component into test.
    recs = _records(sample_entries, artifact_entry)
    kept, _ = filter_candidates(recs, _cfg())
    class_map = {r.entry_id: classify_components(r, _cfg()) for r in kept}
    entry_classes = {eid: info["classes"] for eid, info in class_map.items()}
    cr = build_clusters(kept, _cfg())
    cfg_min = _cfg(test_min_per_class={"metal": 1})
    res = assign_splits(cr, cfg_min, entry_classes=entry_classes)
    # The floor is met and the structural no-leakage invariant still holds.
    metal_in_test = sum(
        1 for e, s in res.entry_split.items() if s == "test" and "metal" in entry_classes.get(e, [])
    )
    assert metal_in_test >= 1
    assert res.minimum_shortfalls == {}
    check_no_leakage(res, cr)


def test_test_minimums_report_shortfall_when_supply_short(sample_entries, artifact_entry):
    recs = _records(sample_entries, artifact_entry)
    kept, _ = filter_candidates(recs, _cfg())
    class_map = {r.entry_id: classify_components(r, _cfg()) for r in kept}
    entry_classes = {eid: info["classes"] for eid, info in class_map.items()}
    cr = build_clusters(kept, _cfg())
    # Demand far more metal entries than exist -> shortfall reported, not forced.
    cfg_min = _cfg(test_min_per_class={"metal": 999})
    res = assign_splits(cr, cfg_min, entry_classes=entry_classes)
    assert res.minimum_shortfalls.get("metal", 0) > 0
    check_no_leakage(res, cr)


def test_minimums_off_by_default_matches_pure_hash(sample_entries, artifact_entry):
    recs = _records(sample_entries, artifact_entry)
    kept, _ = filter_candidates(recs, _cfg())
    cr = build_clusters(kept, _cfg())
    base = assign_splits(cr, _cfg())  # no entry_classes, default empty minimums
    assert base.minimum_shortfalls == {}
    # Providing classes but no minimums must not change the assignment.
    class_map = {r.entry_id: classify_components(r, _cfg()) for r in kept}
    ec = {eid: info["classes"] for eid, info in class_map.items()}
    same = assign_splits(cr, _cfg(), entry_classes=ec)
    assert same.entry_split == base.entry_split


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
    # The manifest holds counts only (per-entry lists live in train/val/test.json).
    total = sum(m["splits"]["entry_counts"].values())
    assert total == 3  # 4HHB, 1A1F, pdb_00009xyz


def test_manifest_is_lightweight(sample_entries, artifact_entry):
    import json

    m = _full_manifest(sample_entries, artifact_entry, _cfg())
    # No per-entry arrays in the manifest itself — only counts + a files index.
    assert "entries" not in m["splits"]
    assert "entry_clusters" not in m["splits"]
    assert "classes" not in m["ligands"]
    assert "tiers" not in m["ligands"]
    assert set(m["files"]["splits"]) == {"train", "val", "test"}
    # Sanity: the whole manifest is tiny (well under 10 KB for 3 entries).
    assert len(json.dumps(m)) < 10_000


def test_loader_roundtrip(tmp_path, sample_entries, artifact_entry):
    from ifsplit.manifest import write_classes, write_clusters, write_split_files

    cfg = _cfg()
    recs = _records(sample_entries, artifact_entry)
    kept, drops = filter_candidates(recs, cfg)
    class_map = {r.entry_id: classify_components(r, cfg) for r in kept}
    cr = build_clusters(kept, cfg)
    sp = assign_splits(cr, cfg)
    m = build_manifest(
        cfg,
        candidates_sha256="deadbeef",
        n_candidates=len(recs),
        drops=drops,
        drop_counts=drop_summary(drops),
        clusters=cr,
        splits=sp,
        class_map=class_map,
    )
    write_split_files(sp, class_map, tmp_path)
    write_clusters(cr.entry_to_cluster, tmp_path)
    write_classes(class_map, tmp_path)
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


# --------------- union-find: structural leakage prevention ----------------- #
def _protein_record(entry_id: str, cluster30_ids: list[int]) -> CandidateRecord:
    """A record whose protein chains sit in the given raw clusters (id at 30%)."""
    pes = [
        PolymerEntity(
            entity_id=f"{entry_id}_{i + 1}",
            polymer_type="Protein",
            seq_len=100,
            seq="A" * 100,
            cluster_ids={30: cid},
        )
        for i, cid in enumerate(cluster30_ids)
    ]
    return CandidateRecord(
        entry_id=entry_id,
        methods=["X-RAY DIFFRACTION"],
        resolution_A=2.0,
        release_date="2020-01-01",
        deposited_residues=100,
        assemblies={f"{entry_id}-1": 100},
        polymer_entities=pes,
        nonpolymer_comps=[],
        bound_components=[],
        affinity_comp_ids=[],
    )


def test_union_find_merges_bridged_clusters_no_leakage():
    cfg = _cfg()
    # X bridges raw clusters 1 and 2 via two chains; Y is in 1, Z is in 2.
    # Without union-find, clusters 1 and 2 could hash to different splits and Y/Z
    # would leak X's sequences across splits. With it, {1,2} is one component.
    recs = [
        _protein_record("X1AA", [1, 2]),
        _protein_record("Y2BB", [1]),
        _protein_record("Z3CC", [2]),
    ]
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    assert cr.n_raw_clusters == 2
    assert cr.n_clusters == 1  # the two raw clusters merged into one component
    res = assign_splits(cr, cfg)
    check_no_leakage(res, cr)  # would raise if 1 and 2 split apart
    assert len(set(res.entry_split.values())) == 1  # all three co-assigned


def test_independent_clusters_can_differ_and_check_passes():
    cfg = _cfg()
    # Two unrelated single-chain entries: separate components, no shared sequence.
    recs = [_protein_record("AAA1", [10]), _protein_record("BBB2", [20])]
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    assert cr.n_clusters == 2
    res = assign_splits(cr, cfg)
    check_no_leakage(res, cr)  # passes regardless of which splits they land in


def _fold_record(entry_id: str, cluster30: int, families: dict[str, list[str]]) -> CandidateRecord:
    """Single-chain record in raw cluster ``cluster30`` carrying structural ``families``."""
    pe = PolymerEntity(
        entity_id=f"{entry_id}_1",
        polymer_type="Protein",
        seq_len=100,
        seq="A" * 100,
        cluster_ids={30: cluster30},
        structural_families=families,
    )
    return CandidateRecord(
        entry_id=entry_id,
        methods=["X-RAY DIFFRACTION"],
        resolution_A=2.0,
        release_date="2020-01-01",
        deposited_residues=100,
        assemblies={f"{entry_id}-1": 100},
        polymer_entities=[pe],
        nonpolymer_comps=[],
        bound_components=[],
        affinity_comp_ids=[],
    )


def test_structural_clustering_merges_same_fold():
    # Two entries in DIFFERENT sequence clusters (10, 20) but the same CATH
    # superfamily. Sequence-only leaves them separable (a fold-leakage risk);
    # cath clustering folds them into one leakage-safe component.
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.490.10"]}),
        _fold_record("BBB2", 20, {"cath": ["1.10.490.10"]}),
    ]
    cfg = _cfg(structural_clustering="cath")
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    assert cr.n_raw_clusters == 2
    assert cr.n_seq_only_components == 2  # sequence edges alone: two components
    assert cr.n_clusters == 1  # ... folded into one by the shared superfamily
    assert cr.structural_method == "cath"
    assert cr.n_structural_families == 1

    # Off -> prior behavior, two separate components.
    cr_off = build_clusters(kept, _cfg(structural_clustering="off"))
    assert cr_off.n_clusters == 2
    assert cr_off.structural_method == "off"


def test_structural_method_is_selectable():
    # Same ECOD family name but different CATH codes: 'cath' keeps them apart,
    # 'ecod' merges them. The classification method is the config's to choose.
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.1.1"], "ecod": ["Globin-like"]}),
        _fold_record("BBB2", 20, {"cath": ["2.20.2.2"], "ecod": ["Globin-like"]}),
    ]
    kept, _ = filter_candidates(recs, _cfg())
    assert build_clusters(kept, _cfg(structural_clustering="cath")).n_clusters == 2
    assert build_clusters(kept, _cfg(structural_clustering="ecod")).n_clusters == 1


def test_structural_clustering_keeps_split_leakage_safe():
    # A fold shared across two entries must land in ONE split, never straddle.
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.490.10"]}),
        _fold_record("BBB2", 20, {"cath": ["1.10.490.10"]}),
    ]
    cfg = _cfg(structural_clustering="cath")
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    res = assign_splits(cr, cfg)
    check_no_leakage(res, cr)
    assert len(set(res.entry_split.values())) == 1  # same fold -> same split
