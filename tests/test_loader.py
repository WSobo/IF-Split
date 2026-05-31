"""Stage 8 tests: manifest loader + cluster-balanced sampling (offline)."""

from __future__ import annotations

from pathlib import Path

from ifsplit.cluster import build_clusters
from ifsplit.config import load_config
from ifsplit.dataset import load_dataset
from ifsplit.ligands import classify_components
from ifsplit.manifest import build_manifest, write_manifest
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
    return write_manifest(m, tmp_path)


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


def test_manifest_has_tiers_and_ambiguous_counts(tmp_path, sample_entries, artifact_entry):
    from ifsplit.manifest import read_manifest

    m = read_manifest(_build(tmp_path, sample_entries, artifact_entry))
    assert "tiers" in m["ligands"]
    assert "per_split_ambiguous_counts" in m["splits"]
    # 4HHB's PO4 must be tiered as an artifact somewhere in the ligand tiers.
    tiers_4hhb = m["ligands"]["tiers"].get("4HHB", {})
    assert tiers_4hhb.get("PO4", {}).get("tier") == "artifact"
    assert tiers_4hhb.get("HEM", {}).get("tier") == "functional"
