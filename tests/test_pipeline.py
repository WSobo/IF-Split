"""Phases 3/5/6/7 tests: filter, cluster, split, manifest, determinism (offline)."""

from __future__ import annotations

from pathlib import Path

import pytest

from ifsplit.cluster import build_clusters
from ifsplit.config import Config, load_config
from ifsplit.dataset import load_dataset
from ifsplit.ligands import classify_components
from ifsplit.manifest import build_manifest, summarize_manifest, write_manifest
from ifsplit.parse import (
    DROP_CLASHSCORE,
    DROP_EM_INCLUSION,
    DROP_NO_PROTEIN,
    DROP_NO_SEQUENCE,
    DROP_NO_VALIDATION,
    DROP_RESOLUTION,
    DROP_RFREE,
    DROP_SEQUENCE_TOO_SHORT,
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


def _seq_record(entry_id: str, seq: str) -> CandidateRecord:
    """A single-protein-chain record carrying the given canonical sequence."""
    return CandidateRecord(
        entry_id=entry_id,
        methods=["X-RAY DIFFRACTION"],
        resolution_A=2.0,
        release_date="2020-01-01",
        deposited_residues=len(seq),
        assemblies={f"{entry_id}-1": len(seq)},
        polymer_entities=[
            PolymerEntity(
                entity_id=f"{entry_id}_1",
                polymer_type="Protein",
                seq_len=len(seq),
                seq=seq,
                cluster_ids={30: 1},
            )
        ],
        nonpolymer_comps=[],
        bound_components=[],
        affinity_comp_ids=[],
    )


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


def test_size_cap_boundary_keeps_exactly_max_residues():
    # max_total_residues is the max KEPT (drop if > it), so an entry with exactly
    # that many residues is kept — LigandMPNN's "< 6000" is keep <= 5999, not <= 5998.
    rec = _seq_record("BND1", "A" * 100)  # assembly-1 count == 100
    kept, _ = filter_candidates([rec], _cfg(max_total_residues=100))
    assert [r.entry_id for r in kept] == ["BND1"]
    # One residue over the cap is dropped.
    kept2, drops = filter_candidates([rec], _cfg(max_total_residues=99))
    assert kept2 == []
    assert drops[0]["reason"] == DROP_TOO_LARGE


def test_filter_drops_poly_unk_sequence():
    # Every protein chain is all-X (poly-UNK): no known residue identities, so no
    # learnable inverse-folding label -> always dropped (even at the default min=0).
    kept, drops = filter_candidates([_seq_record("UNK1", "X" * 80)], _cfg())
    assert kept == []
    assert drops[0]["reason"] == DROP_NO_SEQUENCE


def test_filter_keeps_partially_modeled_sequence():
    # A chain with even a few modeled residues is usable at the default (min=0).
    kept, _ = filter_candidates([_seq_record("OK1", "X" * 70 + "ACDEFGHIKL")], _cfg())
    assert [r.entry_id for r in kept] == ["OK1"]


def test_min_modeled_residues_drops_short_chain():
    rec = _seq_record("SHRT", "ACDEFGHIKLMNPQR")  # 15 modeled residues
    kept, drops = filter_candidates([rec], _cfg(min_modeled_residues=20))
    assert kept == []
    assert drops[0]["reason"] == DROP_SEQUENCE_TOO_SHORT
    assert drops[0]["modeled"] == 15
    # Off by default (min=0): the same record is kept.
    kept2, _ = filter_candidates([rec], _cfg())
    assert [r.entry_id for r in kept2] == ["SHRT"]


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


# ------------------------ Stage 3: resolution re-filter -------------------- #
def test_filter_resolution_refilter_is_auditable():
    # Stage 3 re-derives the resolution cut (Search applied it too) so it is auditable
    # from candidates.jsonl. An entry over the cap is dropped; at the cap it is kept.
    over = _seq_record("RES1", "A" * 50).model_copy(update={"resolution_A": 3.8})
    kept, drops = filter_candidates([over], _cfg())  # default cap 3.5
    assert kept == []
    assert drops[0]["reason"] == DROP_RESOLUTION
    assert drops[0]["resolution"] == 3.8
    at_cap = _seq_record("RES2", "A" * 50).model_copy(update={"resolution_A": 3.5})
    assert [r.entry_id for r in filter_candidates([at_cap], _cfg())[0]] == ["RES2"]


def test_filter_resolution_missing_is_kept():
    rec = _seq_record("RESN", "A" * 50).model_copy(update={"resolution_A": None})
    assert [r.entry_id for r in filter_candidates([rec], _cfg())[0]] == ["RESN"]


def test_per_method_resolution_cap():
    cfg = _cfg(resolution_max_A_by_method={"ELECTRON MICROSCOPY": 3.0})
    # A 3.2 A cryo-EM entry is dropped by the tighter EM cap...
    em = _seq_record("EM1", "A" * 50).model_copy(
        update={"resolution_A": 3.2, "methods": ["ELECTRON MICROSCOPY"]}
    )
    kept, drops = filter_candidates([em], cfg)
    assert kept == []
    assert drops[0]["reason"] == DROP_RESOLUTION
    # ...but a 3.2 A X-ray entry passes (its cap is still the global 3.5).
    xr = _seq_record("XR1", "A" * 50).model_copy(
        update={"resolution_A": 3.2, "methods": ["X-RAY DIFFRACTION"]}
    )
    assert [r.entry_id for r in filter_candidates([xr], cfg)[0]] == ["XR1"]


def test_search_resolution_cap_is_loosest():
    # The Search query must pull a superset: the loosest cap across enabled methods.
    assert (
        _cfg(resolution_max_A_by_method={"ELECTRON MICROSCOPY": 3.0}).search_resolution_cap() == 3.5
    )
    assert (
        _cfg(
            resolution_max_A_by_method={"X-RAY DIFFRACTION": 4.0, "ELECTRON MICROSCOPY": 3.0}
        ).search_resolution_cap()
        == 4.0
    )


# ---------------------- Stage 3: cryo-EM map-fit floor --------------------- #
def _with_em_inclusion(entry_id: str, value: float | None) -> CandidateRecord:
    rec = _seq_record(entry_id, "A" * 50)
    return rec.model_copy(
        update={"quality": rec.quality.model_copy(update={"em_backbone_inclusion": value})}
    )


def test_em_backbone_inclusion_floor_drops_low_fit():
    cfg = _cfg(min_em_backbone_inclusion=0.7)
    kept, drops = filter_candidates([_with_em_inclusion("EMLO", 0.6)], cfg)
    assert kept == []
    assert drops[0]["reason"] == DROP_EM_INCLUSION
    assert drops[0]["value"] == 0.6
    # At/above the floor is kept.
    assert [r.entry_id for r in filter_candidates([_with_em_inclusion("EMHI", 0.85)], cfg)[0]] == [
        "EMHI"
    ]


def test_em_floor_ignores_entries_without_the_metric():
    # X-ray entries have no em_backbone_inclusion -> the floor never drops them.
    cfg = _cfg(min_em_backbone_inclusion=0.7)
    assert [r.entry_id for r in filter_candidates([_with_em_inclusion("XNOEM", None)], cfg)[0]] == [
        "XNOEM"
    ]


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


# ------------------------- resplit (offline, no RCSB) ---------------------- #
def test_resplit_reproduces_build_offline(tmp_path, fake_client):
    import argparse

    from ifsplit.cli import _run_pipeline, cmd_resplit
    from ifsplit.enumerate import enumerate_candidates
    from ifsplit.manifest import read_lock, read_manifest

    cfg = _cfg()
    # Reference: enumerate to `ref` (writes candidates.jsonl) + run Stages 3-7 there.
    ref = tmp_path / "ref"
    records, cand_path, sha = enumerate_candidates(cfg, ref, client=fake_client)
    _run_pipeline(cfg, records, sha, ref, limit=None, registry_path=None)
    man_ref = read_manifest(ref / "manifest.json")

    # Resplit re-derives from the SAME candidates.jsonl offline (no client).
    out = tmp_path / "out"
    args = argparse.Namespace(
        config=str(DEFAULT_CONFIG), candidates=str(cand_path), out=str(out), registry=None
    )
    assert cmd_resplit(args) == 0
    man_out = read_manifest(out / "manifest.json")
    # Same snapshot + config -> identical split output; the offline sha matches the
    # sha enumerate computed (the file bytes hash identically).
    assert man_out["splits"]["entry_counts"] == man_ref["splits"]["entry_counts"]
    out_lock = read_lock(out / "dataset.lock")
    assert out_lock["candidates"]["sha256"] == sha
    # The lock records how it was produced: resplit vs a live build.
    assert out_lock["source"] == "resplit"
    assert read_lock(ref / "dataset.lock")["source"] == "build"


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


def _clusters_from(members):
    """Build a ClusterResult from a {component_key: [entry_ids]} map (for split tests)."""
    from ifsplit.cluster import ClusterResult

    e2c, eraw = {}, {}
    for key, ents in members.items():
        for e in ents:
            e2c[e] = key
            eraw[e] = [key]
    return ClusterResult(
        30, dict(sorted(e2c.items())), dict(sorted(members.items())), dict(sorted(eraw.items()))
    )


def _members(prefix, n, start=0):
    """n components with tail sizes 1..20, so the balanced val/test fill is exercised."""
    return {
        f"{prefix}{i:05d}": [f"{prefix}{i:05d}_{j}" for j in range(1 + i % 20)]
        for i in range(start, start + n)
    }


def test_balanced_growth_stability_needs_registry():
    # A balanced split's val/test fill boundaries scale with total entries, so WITHOUT a
    # registry a grown snapshot moves prior components across splits; the registry pins them.
    cfg = _cfg(split_strategy="balanced")
    a_mem = _members("A", 1000)
    clusters_a = _clusters_from(a_mem)
    clusters_b = _clusters_from({**a_mem, **_members("B", 1000)})  # A + 1000 new components

    sp_a = assign_splits(clusters_a, cfg)
    sp_b_noreg = assign_splits(clusters_b, cfg)  # registry defaults to {}
    moved = sum(1 for k in a_mem if sp_a.cluster_split[k] != sp_b_noreg.cluster_split[k])
    assert moved > 0, "balanced must NOT be assumed growth-stable without a registry"

    sp_b_reg = assign_splits(clusters_b, cfg, registry=sp_a.cluster_split)
    moved_reg = sum(1 for k in a_mem if sp_a.cluster_split[k] != sp_b_reg.cluster_split[k])
    assert moved_reg == 0, "the registry must restore growth stability for balanced"


def test_balanced_rebuild_auto_pins_registry(tmp_path, sample_entries, artifact_entry, capsys):
    # The fix: a balanced rebuild into the same --out auto-adopts the prior registry when
    # the config matches, so the lock records it and the manifest reports growth_stable.
    from ifsplit.cli import _run_pipeline
    from ifsplit.manifest import read_lock, read_manifest

    cfg = _cfg(split_strategy="balanced")
    recs = _records(sample_entries, artifact_entry)
    out = tmp_path / "d"

    _run_pipeline(cfg, recs, "sha", out, limit=None, registry_path=None)  # first build
    assert read_lock(out / "dataset.lock")["split"]["registry_sha256"] is None  # registry-free

    capsys.readouterr()
    _run_pipeline(cfg, recs, "sha", out, limit=None, registry_path=None)  # in-place rebuild
    assert "pinning" in capsys.readouterr().out
    assert read_lock(out / "dataset.lock")["split"]["registry_sha256"] is not None
    assert read_manifest(out / "manifest.json")["splits"]["growth_stable"] is True

    _run_pipeline(cfg, recs, "sha", out, limit=None, registry_path=None, fresh=True)  # opt out
    assert read_lock(out / "dataset.lock")["split"]["registry_sha256"] is None


def test_hash_rebuild_stays_registry_free(tmp_path, sample_entries, artifact_entry):
    # No regression: hash is input-independent, so it is never auto-pinned and stays
    # registry-free (verify can still certify the split output).
    from ifsplit.cli import _run_pipeline
    from ifsplit.manifest import read_lock, read_manifest

    cfg = _cfg()  # hash (default)
    recs = _records(sample_entries, artifact_entry)
    out = tmp_path / "d"
    _run_pipeline(cfg, recs, "sha", out, limit=None, registry_path=None)
    _run_pipeline(cfg, recs, "sha", out, limit=None, registry_path=None)  # rebuild
    assert read_lock(out / "dataset.lock")["split"]["registry_sha256"] is None
    assert read_manifest(out / "manifest.json")["splits"]["growth_stable"] is True


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


def test_balanced_strategy_caps_dominant_fold_and_balances_entries():
    # 400 entries share one dominant cluster (a mega-fold); 200 are singletons.
    # Per-component hashing would let the 400-entry fold balloon a split; balanced
    # caps it to train and fills val/test to their ENTRY targets from the tail.
    cfg = _cfg(split_strategy="balanced")
    recs = [_protein_record(f"D{i:04d}", [1]) for i in range(400)]
    recs += [_protein_record(f"S{i:04d}", [1000 + i]) for i in range(200)]
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    res = assign_splits(cr, cfg)
    check_no_leakage(res, cr)
    assert res.strategy == "balanced"
    assert res.capped_folds == 1  # the 400-entry fold -> train
    c = res.counts
    assert (c["train"], c["val"], c["test"]) == (480, 60, 60)  # 80/10/10 by entries
    assert not res.balance_gaps


def test_balanced_strategy_reports_thin_tail_gap():
    # Almost everything in one fold: the tail can't fill val+test to 10% each.
    # The gap is reported, not forced (never breaks leakage safety to hit a target).
    cfg = _cfg(split_strategy="balanced")
    recs = [_protein_record(f"D{i:04d}", [1]) for i in range(580)]
    recs += [_protein_record(f"S{i:04d}", [1000 + i]) for i in range(20)]
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    res = assign_splits(cr, cfg)
    check_no_leakage(res, cr)
    assert res.balance_gaps  # tail too thin -> reported
    assert "val" in res.balance_gaps


# ---------- negative leakage: the guard MUST fire on a corrupted partition ------- #
# These prove check_no_leakage is a real invariant, not a happy-path pass — it
# actually raises when a sequence cluster or a fold (super)family spans two splits.


def test_check_no_leakage_fires_on_shared_sequence_cluster():
    # X bridges raw clusters 1 & 2; Y is in 1. A partition that puts X and Y in
    # different splits leaks X's sequence via raw cluster 1 — the guard must catch it.
    cfg = _cfg()
    recs = [
        _protein_record("X1AA", [1, 2]),
        _protein_record("Y2BB", [1]),
        _protein_record("Z3CC", [2]),
    ]
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    res = assign_splits(cr, cfg)
    check_no_leakage(res, cr)  # the real (valid) partition passes
    # Corrupt it: move Y to a different split from X (they share raw cluster 1).
    res.entry_split = dict(res.entry_split)
    res.entry_split["Y2BB"] = "test" if res.entry_split["X1AA"] != "test" else "train"
    with pytest.raises(AssertionError, match="raw cluster"):
        check_no_leakage(res, cr)


def test_check_no_leakage_fires_on_fold_span_when_structural_on():
    # Two entries with DIFFERENT sequence clusters but the SAME CATH superfamily are
    # one component under structural clustering. Splitting them apart shares no
    # sequence cluster (sequence check passes) — only the fold-level guard catches
    # it, proving structural_clustering's guarantee is actually enforced.
    cfg = _cfg(structural_clustering="cath")
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.490.10"]}),
        _fold_record("BBB2", 20, {"cath": ["1.10.490.10"]}),
    ]
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    assert cr.n_clusters == 1  # merged into one component by shared fold
    res = assign_splits(cr, cfg)
    check_no_leakage(res, cr)
    res.entry_split = {"AAA1": "train", "BBB2": "test"}  # force same fold across splits
    with pytest.raises(AssertionError, match="fold leakage"):
        check_no_leakage(res, cr)


def test_fold_leakage_guard_is_a_noop_when_structural_off():
    # Same two same-fold entries, structural OFF: distinct sequence clusters, so
    # they may legitimately land in different splits. entry_families is empty, so
    # the fold guard must NOT fire (no over-merging beyond what the user asked for).
    cfg = _cfg(structural_clustering="off")
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.490.10"]}),
        _fold_record("BBB2", 20, {"cath": ["1.10.490.10"]}),
    ]
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    assert cr.entry_families == {}  # off -> no fold edges recorded
    res = assign_splits(cr, cfg)
    res.entry_split = {"AAA1": "train", "BBB2": "test"}
    check_no_leakage(res, cr)  # must NOT raise — distinct folds may differ when off


def test_single_chain_only_filter(sample_entries):
    # 1A1F is a protein+DNA complex (multiple polymer entities); 4HHB has two protein
    # entities (alpha/beta). single_chain_only keeps only single-protein-entity records.
    recs = [CandidateRecord.from_data_api(e) for e in sample_entries.values()]
    kept, drops = filter_candidates(recs, _cfg(single_chain_only=True))
    assert all(len(r.polymer_entities) == 1 for r in kept)
    reasons = {d["entry_id"]: d["reason"] for d in drops}
    assert reasons.get("1A1F") == "not_single_chain"  # protein+DNA -> dropped
    # A genuine single-chain record passes.
    kept2, _ = filter_candidates([_seq_record("ONE1", "A" * 80)], _cfg(single_chain_only=True))
    assert [r.entry_id for r in kept2] == ["ONE1"]
    # Off by default: the complex is kept.
    kept3, _ = filter_candidates(recs, _cfg())
    assert "1A1F" in {r.entry_id for r in kept3}


def test_manifest_tier_reason_histogram(sample_entries, artifact_entry):
    # The tier-reason histogram summarizes every curation call (a "tier:reason" key
    # per component) so the distribution is auditable without the per-component file.
    m = _full_manifest(sample_entries, artifact_entry, _cfg())
    trc = m["ligands"]["tier_reason_counts"]
    assert isinstance(trc, dict) and trc and all(isinstance(v, int) for v in trc.values())
    assert all(":" in k for k in trc)  # keys are "tier:reason"


def test_manifest_fold_coverage_counts_distinct_folds():
    # per_split_fold_coverage counts the distinct structural families held in each
    # split, plus the unclassified count — the residual-leakage ceiling: entries no
    # fold taxonomy classifies are held out by sequence only, not by fold.
    cfg = _cfg(structural_clustering="cath")
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.1.1"]}),
        _fold_record("BBB2", 20, {"cath": ["2.20.2.2"]}),
        _fold_record("CCC3", 30, {}),  # no fold classification -> unclassified
    ]
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
    cov = m["splits"]["per_split_fold_coverage"]
    assert set(cov) == {"train", "val", "test"}
    for c in cov.values():
        assert set(c) == {
            "total_entries",
            "classified_entries",
            "unclassified_entries",
            "n_distinct_folds",
        }
        assert c["unclassified_entries"] == c["total_entries"] - c["classified_entries"]
    assert sum(c["n_distinct_folds"] for c in cov.values()) == 2  # two distinct folds
    assert sum(c["classified_entries"] for c in cov.values()) == 2
    assert sum(c["total_entries"] for c in cov.values()) == 3
    assert sum(c["unclassified_entries"] for c in cov.values()) == 1  # CCC3, unclassified


def test_summarize_manifest_reports_residual_ceiling(tmp_path, capsys):
    # `stats` surfaces the unclassified fraction per split (the residual-leakage
    # ceiling) whenever fold-aware clustering is on.
    cfg = _cfg(structural_clustering="cath")
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.1.1"]}),
        _fold_record("CCC3", 30, {}),  # unclassified
    ]
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
    path = write_manifest(m, tmp_path)
    assert summarize_manifest(path) == 0
    assert "unclassified" in capsys.readouterr().out
