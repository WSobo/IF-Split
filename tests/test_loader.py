"""Stage 8 tests: manifest loader + cluster-balanced sampling (offline)."""

from __future__ import annotations

from pathlib import Path

from ifsplit.cluster import build_clusters
from ifsplit.config import load_config
from ifsplit.dataset import load_dataset
from ifsplit.ligands import classify_components
from ifsplit.manifest import (
    build_manifest,
    build_targets,
    build_tiers_doc,
    write_classes,
    write_clusters,
    write_manifest,
    write_split_files,
    write_targets,
    write_tiers,
)
from ifsplit.parse import drop_summary, filter_candidates
from ifsplit.schema import CandidateRecord
from ifsplit.split import assign_splits

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config" / "default.yaml"


def _build(tmp_path, sample_entries, artifact_entry):
    cfg = load_config(DEFAULT_CONFIG)
    recs = [CandidateRecord.from_data_api(e) for e in sample_entries.values()]
    recs.append(CandidateRecord.from_data_api(artifact_entry))
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
    # Write the full file set the loader reads (split lists + supporting maps).
    write_split_files(sp, class_map, tmp_path)
    write_clusters(cr.entry_to_cluster, tmp_path)
    write_classes(class_map, tmp_path)
    write_tiers(build_tiers_doc(class_map), tmp_path)
    write_targets(build_targets(class_map, sp, cr), tmp_path)
    return write_manifest(m, tmp_path)


def _all_targets(ds, **kw):
    """Conditioning targets across all splits."""
    return [t for s in ("train", "val", "test") for t in ds.split(s).conditioning_targets(**kw)]


def test_loader_views_total(tmp_path, sample_entries, artifact_entry):
    ds = load_dataset(_build(tmp_path, sample_entries, artifact_entry))
    assert len(ds.train) + len(ds.val) + len(ds.test) == 3
    assert ds.config_hash == load_config(DEFAULT_CONFIG).config_hash()


def test_entry_clusters_present(tmp_path, sample_entries, artifact_entry):
    ds = load_dataset(_build(tmp_path, sample_entries, artifact_entry))
    # Every entry across all splits maps to some cluster key.
    for name in ("train", "val", "test"):
        view = ds.split(name)
        for e in view.entry_ids:
            assert view.entry_clusters[e]  # non-empty key


def test_sample_by_cluster_one_per_cluster(tmp_path, sample_entries, artifact_entry):
    ds = load_dataset(_build(tmp_path, sample_entries, artifact_entry))
    for name in ("train", "val", "test"):
        view = ds.split(name)
        sample = view.sample_by_cluster(seed=0)
        assert len(sample) == len(view.cluster_groups())
        assert len(set(sample)) == len(sample)  # no duplicates
        assert set(sample) <= set(view.entry_ids)


def test_sample_by_cluster_is_deterministic(tmp_path, sample_entries, artifact_entry):
    ds = load_dataset(_build(tmp_path, sample_entries, artifact_entry))
    a = ds.train.sample_by_cluster(seed=7)
    b = ds.train.sample_by_cluster(seed=7)
    assert a == b


def test_split_files_are_plain_id_lists(tmp_path, sample_entries, artifact_entry):
    import json

    mpath = _build(tmp_path, sample_entries, artifact_entry)
    d = mpath.parent
    # train/val/test.json each parse as a flat list of ids; together they total 3.
    total = 0
    for name in ("train", "val", "test"):
        ids = json.loads((d / f"{name}.json").read_text())
        assert isinstance(ids, list)
        total += len(ids)
    assert total == 3
    # 1A1F (zinc-finger/DNA + bound Zn) is a metal-class test/train member; wherever
    # it landed, that split's per-class file must list it if it's in test.


def test_manifest_lean_tiers_in_sidecar(tmp_path, sample_entries, artifact_entry):
    from ifsplit.manifest import TIERS_FILENAME, read_manifest, read_tiers

    mpath = _build(tmp_path, sample_entries, artifact_entry)
    m = read_manifest(mpath)
    # The manifest stays lean: no per-entry arrays, only a files index + counts.
    assert "tiers" not in m["ligands"]
    assert "classes" not in m["ligands"]
    assert "entries" not in m["splits"]
    assert m["files"]["ligand_tiers"] == TIERS_FILENAME
    assert "per_split_ambiguous_counts" in m["splits"]

    # The audit detail lives in the sidecar next to the manifest.
    tiers = read_tiers(mpath.parent / TIERS_FILENAME)
    tiers_4hhb = tiers.get("4HHB", {})
    assert tiers_4hhb.get("PO4", {}).get("tier") == "artifact"
    assert tiers_4hhb.get("HEM", {}).get("tier") == "functional"


def test_conditioning_targets_are_functional_only(tmp_path, sample_entries, artifact_entry):
    ds = load_dataset(_build(tmp_path, sample_entries, artifact_entry))
    targets = _all_targets(ds)
    assert all(t.tier == "functional" for t in targets)
    # 4HHB -> HEM (small_molecule); 1A1F -> ZN (metal) + nucleic_acid.
    tset = {(t.entry_id, t.ligand_class, t.comp_id) for t in targets}
    assert ("4HHB", "small_molecule", "HEM") in tset
    assert ("1A1F", "metal", "ZN") in tset
    assert ("1A1F", "nucleic_acid", None) in tset
    # The His-tag/Ni artifact entry contributes NO conditioning target.
    assert not any(t.entry_id == "pdb_00009xyz" for t in targets)


def test_backbones_include_every_structure(tmp_path, sample_entries, artifact_entry):
    ds = load_dataset(_build(tmp_path, sample_entries, artifact_entry))
    backbones = {e for s in ("train", "val", "test") for e in ds.split(s).backbones}
    assert backbones == {"4HHB", "1A1F", "pdb_00009xyz"}  # every kept structure
    conditioned = {t.entry_id for t in _all_targets(ds)}
    # The artifact-only entry is a backbone but not a conditioning target.
    assert "pdb_00009xyz" in backbones
    assert "pdb_00009xyz" not in conditioned
    assert conditioned < backbones


def test_targets_by_entry_groups_1a1f(tmp_path, sample_entries, artifact_entry):
    ds = load_dataset(_build(tmp_path, sample_entries, artifact_entry))
    for s in ("train", "val", "test"):
        grouped = ds.split(s).targets_by_entry()
        if "1A1F" in grouped:
            classes = {t.ligand_class for t in grouped["1A1F"]}
            assert classes == {"metal", "nucleic_acid"}  # both of 1A1F's targets, grouped
            return
    raise AssertionError("1A1F not found in any split")


def test_class_filter_on_conditioning_targets(tmp_path, sample_entries, artifact_entry):
    ds = load_dataset(_build(tmp_path, sample_entries, artifact_entry))
    metals = _all_targets(ds, classes=["metal"])
    assert metals and all(t.ligand_class == "metal" for t in metals)
    assert all(t.entry_id == "1A1F" for t in metals)  # only 1A1F has a functional metal


def test_manifest_training_summary(tmp_path, sample_entries, artifact_entry):
    from ifsplit.manifest import read_manifest

    m = read_manifest(_build(tmp_path, sample_entries, artifact_entry))
    tr = m["training"]
    assert tr["n_backbones"] == 3
    assert tr["n_conditioning_targets"] == 3  # HEM + ZN + nucleic_acid
    assert m["files"]["targets"] == "targets.jsonl"


def test_per_class_test_files_written(tmp_path, sample_entries, artifact_entry):
    # Force everything into test so the per-class test files are populated.
    import json

    from ifsplit.cluster import build_clusters as _bc
    from ifsplit.manifest import TEST_SUBDIR, write_split_files

    cfg = load_config(DEFAULT_CONFIG)
    recs = [CandidateRecord.from_data_api(e) for e in sample_entries.values()]
    recs.append(CandidateRecord.from_data_api(artifact_entry))
    kept, _ = filter_candidates(recs, cfg)
    class_map = {r.entry_id: classify_components(r, cfg) for r in kept}
    cr = _bc(kept, cfg)
    reg = {k: "test" for k in cr.cluster_members}
    sp = assign_splits(cr, cfg, registry=reg)
    paths = write_split_files(sp, class_map, tmp_path)

    # 1A1F has a functional metal (bound Zn) -> metal_test.json lists it.
    metal_file = tmp_path / TEST_SUBDIR / "metal_test.json"
    assert metal_file.exists()
    assert "1A1F" in json.loads(metal_file.read_text())
    assert any(k.startswith("test:") for k in paths)
