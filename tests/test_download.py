"""Stage 2 tests: sharding, paths, fetch orchestration, index (offline).

A FakeFetcher stands in for the network so the hydration packaging (layout,
index, dataset card, resume) is tested without downloading anything.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from ifsplit.cluster import build_clusters
from ifsplit.config import load_config
from ifsplit.download import (
    FetchResult,
    core_id,
    filename_for,
    rel_path_for,
    shard_for,
    url_for,
)
from ifsplit.hydrate import hydrate, select_targets
from ifsplit.ligands import classify_components
from ifsplit.manifest import build_manifest, write_manifest
from ifsplit.parse import drop_summary, filter_candidates
from ifsplit.schema import CandidateRecord
from ifsplit.split import assign_splits

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config" / "default.yaml"


# ------------------------------ id / path logic ---------------------------- #
def test_core_id_legacy_and_extended():
    assert core_id("4HHB") == "4hhb"
    assert core_id("pdb_00004hhb") == "4hhb"
    assert core_id("pdb_00009xyz") == "9xyz"


def test_shard_is_pdb_divided_scheme():
    assert shard_for("4HHB") == "hh"  # middle two chars
    assert shard_for("1ABC") == "ab"
    assert shard_for("pdb_00004hhb") == "hh"  # extended maps the same


def test_filename_and_url_assembly_vs_au():
    assert filename_for("4HHB", assembly=True) == "4hhb-assembly1.cif.gz"
    assert filename_for("4HHB", assembly=False) == "4hhb.cif.gz"
    assert url_for("4HHB", assembly=True).endswith("/4hhb-assembly1.cif.gz")


def test_rel_path_is_split_partitioned_and_sharded():
    p = rel_path_for("4HHB", "train", assembly=True)
    assert p == Path("structures/train/hh/4hhb-assembly1.cif.gz")


# ------------------------------ fake fetcher ------------------------------- #
class FakeFetcher:
    """In-memory stand-in: writes tiny valid .cif.gz files, no network."""

    def __init__(self, assembly=True):
        self.assembly = assembly

    def estimate_bytes(self, entry_ids, sample=12):
        return 1234 * len(entry_ids)

    def fetch(self, targets, root, *, progress=None):
        res = FetchResult()
        for eid, split in targets:
            rel = rel_path_for(eid, split, assembly=self.assembly)
            dest = Path(root) / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            payload = gzip.compress(f"data_{core_id(eid)}\n".encode())
            if dest.exists():
                res.skipped.append(eid)
                status = "skipped"
            else:
                dest.write_bytes(payload)
                res.fetched.append(eid)
                status = "fetched"
            import hashlib

            res.index_rows.append(
                {
                    "entry_id": eid,
                    "split": split,
                    "path": str(rel),
                    "sha256": hashlib.sha256(dest.read_bytes()).hexdigest(),
                    "status": status,
                }
            )
        res.index_rows.sort(key=lambda r: (r["split"], r["entry_id"]))
        return res

    def close(self):
        pass


@pytest.fixture
def manifest_path(tmp_path, sample_entries, artifact_entry):
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
    return write_manifest(m, tmp_path / "src")


def test_select_targets_orders_and_scopes(manifest_path):
    from ifsplit.manifest import read_manifest

    m = read_manifest(manifest_path)
    all_targets = select_targets(m, ["train", "val", "test"])
    assert len(all_targets) == 3  # 4HHB, 1A1F, pdb_00009xyz
    # Subsetting to one split returns only that split's entries.
    test_only = select_targets(m, ["test"])
    assert all(s == "test" for _, s in test_only)


def test_hydrate_builds_ml_tree(tmp_path, manifest_path):
    root = tmp_path / "ds"
    summary = hydrate(manifest_path, root, splits=["train", "val", "test"], fetcher=FakeFetcher())
    assert summary["fetched"] == 3
    # Files land split-partitioned + sharded.
    cifs = list(root.glob("structures/*/*/*.cif.gz"))
    assert len(cifs) == 3
    for c in cifs:
        assert gzip.decompress(c.read_bytes())  # valid gzip
    # Index + card + manifest copy all present.
    assert (root / "index.jsonl").exists()
    assert (root / "manifest.json").exists()
    assert (root / "DATASET_CARD.md").exists()
    card = (root / "DATASET_CARD.md").read_text()
    assert "structures/<split>/<shard>" in card


def test_hydrate_is_resumable(tmp_path, manifest_path):
    # Fetch all splits so the test doesn't depend on which split the fixture
    # entries hash into (with the default salt they all land in one split).
    root = tmp_path / "ds"
    splits = ["train", "val", "test"]
    a = hydrate(manifest_path, root, splits=splits, fetcher=FakeFetcher())
    b = hydrate(manifest_path, root, splits=splits, fetcher=FakeFetcher())
    assert a["fetched"] >= 1
    assert b["fetched"] == 0  # second run skips everything
    assert b["skipped"] == a["fetched"]


def test_index_jsonl_rows_have_integrity_fields(tmp_path, manifest_path):
    root = tmp_path / "ds"
    hydrate(manifest_path, root, splits=["train", "val", "test"], fetcher=FakeFetcher())
    rows = [json.loads(line) for line in (root / "index.jsonl").read_text().splitlines()]
    assert rows, "index.jsonl is empty"
    for r in rows:
        assert set(r) >= {"entry_id", "split", "path", "sha256", "cluster", "ligand_classes"}
        assert len(r["sha256"]) == 64
