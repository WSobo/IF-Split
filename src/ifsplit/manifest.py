"""Stage 7 - Snapshot lock, manifest, split registry, and verify/stats commands.

Artifacts written to the output dir:

- ``dataset.lock``   - reproduction anchor: embedded config + canonical
  ``candidates.jsonl`` hash + entry-id list. ``verify`` re-enumerates from it and
  reports drift (added/removed entries, hash match), warning not failing.
- ``manifest.json``  - small (~KB) provenance record: config, drop log, per-split
  + per-class counts, cluster/component stats, and a ``files`` index pointing at
  the data files below. No per-entry arrays, so it stays tiny at any scale. Built
  as a pure function of its inputs (no wall-clock fields) -> byte-identical across
  runs of the same config.

The split itself is plain lists of PDB ids, each its own file:
- ``train.json`` / ``val.json`` / ``test.json`` - the entry ids in each split
  (one id per line; grepable and trivially loadable).
- ``test/<class>_test.json`` - the test ids carrying each functional ligand class
  (``metal`` / ``small_molecule`` / ``nucleic_acid``), for per-class evaluation.

Supporting maps (only needed for sampling / curation, not to read the split):
- ``clusters.json`` - entry_id -> component key (for cluster-balanced sampling).
- ``ligands.classes.json`` - entry_id -> functional class labels.
- ``ligands.tiers.json`` - per-component curation *audit trail* (tier + reason);
  bulky (~24 MB at full-PDB scale), read only by ``fetch`` and curation audits.
- ``splits.registry.json`` - canonical_key -> split, so a later, larger snapshot
  reuses prior assignments instead of re-hashing (growth stability).
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Any

from . import __version__
from .config import Config

LOCK_SCHEMA = "if-split/lock@1"
MANIFEST_SCHEMA = "if-split/manifest@2"
REGISTRY_SCHEMA = "if-split/registry@1"
TIERS_SCHEMA = "if-split/tiers@1"

# Data files written next to manifest.json (referenced by manifest["files"]).
TIERS_FILENAME = "ligands.tiers.json"
CLASSES_FILENAME = "ligands.classes.json"
CLUSTERS_FILENAME = "clusters.json"
SPLIT_FILES = {"train": "train.json", "val": "val.json", "test": "test.json"}
TEST_SUBDIR = "test"  # per-class test-id lists live here: test/<class>_test.json


# --------------------------------------------------------------------------- #
# Lock (Stage 1 reproduction anchor)
# --------------------------------------------------------------------------- #
def build_lock(
    cfg: Config,
    *,
    entry_ids: list[str],
    candidates_sha256: str,
    limit: int | None,
) -> dict[str, Any]:
    """Assemble the lock document (pure; does not touch disk)."""
    return {
        "lock_schema": LOCK_SCHEMA,
        "dataset_version": cfg.dataset_version,
        "if_split_version": __version__,
        "config_hash": cfg.config_hash(),
        "config": cfg.canonical_dict(),
        "selection": {"limit": limit},
        "candidates": {
            "count": len(entry_ids),
            "sha256": candidates_sha256,
            "entry_ids": sorted(entry_ids),
        },
    }


def _write_json(obj: dict[str, Any], path: Path, *, compact: bool = False) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    if compact:
        text = json.dumps(obj, sort_keys=True, separators=(",", ":"))
    else:
        text = json.dumps(obj, indent=2, sort_keys=True)
    path.write_text(text + "\n", encoding="utf-8")
    return path


def write_lock(lock: dict[str, Any], out_dir: str | Path) -> Path:
    return _write_json(lock, Path(out_dir) / "dataset.lock")


def read_lock(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Lock file not found: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------- #
# Split registry (growth stability)
# --------------------------------------------------------------------------- #
def read_registry(path: str | Path) -> dict[str, str]:
    """Load a canonical_key -> split registry, or {} if absent."""
    path = Path(path)
    if not path.exists():
        return {}
    doc = json.loads(path.read_text(encoding="utf-8"))
    return dict(doc.get("assignments", {}))


def write_registry(cluster_split: dict[str, str], out_dir: str | Path) -> Path:
    doc = {
        "registry_schema": REGISTRY_SCHEMA,
        "assignments": dict(sorted(cluster_split.items())),
    }
    return _write_json(doc, Path(out_dir) / "splits.registry.json")


# --------------------------------------------------------------------------- #
# Manifest (human-facing, deterministic)
# --------------------------------------------------------------------------- #
def build_manifest(
    cfg: Config,
    *,
    candidates_sha256: str,
    n_candidates: int,
    drops: list[dict],
    drop_counts: dict[str, int],
    clusters,
    splits,
    class_map: dict[str, dict],
) -> dict[str, Any]:
    """Assemble manifest.json as a pure function of the build outputs."""
    from .split import SPLITS

    # Per-split entry lists.
    per_split: dict[str, list[str]] = {s: [] for s in SPLITS}
    for entry, s in splits.entry_split.items():
        per_split[s].append(entry)
    for s in per_split:
        per_split[s].sort()

    # Per-split, per-class counts at the functional tier (the test-quality view),
    # plus the ambiguous counts so under-/over-confidence is visible.
    def _class_counts(entries, key):
        counts: dict[str, int] = {}
        for entry in entries:
            for cls in class_map.get(entry, {}).get(key, []):
                counts[cls] = counts.get(cls, 0) + 1
        return dict(sorted(counts.items()))

    per_split_class_counts = {s: _class_counts(per_split[s], "classes") for s in SPLITS}
    per_split_ambiguous_counts = {
        s: _class_counts(per_split[s], "ambiguous_classes") for s in SPLITS
    }
    n_artifacts = sum(1 for info in class_map.values() if info.get("purification_artifact"))

    # Per-class test-id files that will be written (only classes that occur).
    test_class_files = {
        cls: f"{TEST_SUBDIR}/{cls}_test.json" for cls in per_split_class_counts["test"]
    }

    # The manifest is small provenance only: NO per-entry arrays live here. The
    # split membership and supporting maps are separate files (see "files").
    return {
        "manifest_schema": MANIFEST_SCHEMA,
        "dataset_version": cfg.dataset_version,
        "if_split_version": __version__,
        "config_hash": cfg.config_hash(),
        "config": cfg.canonical_dict(),
        "candidates": {"count": n_candidates, "sha256": candidates_sha256},
        "filter": {
            "kept": len(splits.entry_split),
            "dropped": len(drops),
            "drop_counts": dict(sorted(drop_counts.items())),
        },
        "clustering": {
            "backend": cfg.clustering_backend,
            "identity": clusters.identity,
            "n_clusters": clusters.n_clusters,
            "n_raw_clusters": clusters.n_raw_clusters,
            "multichain_entries": len(clusters.multichain_entries),
            "unclustered_entries": len(clusters.unclustered_entries),
        },
        "splits": {
            "entry_counts": dict(sorted(splits.counts.items())),
            "cluster_counts": dict(sorted(splits.cluster_counts.items())),
            "per_split_class_counts": per_split_class_counts,
            "per_split_ambiguous_counts": per_split_ambiguous_counts,
            "test_minimum_shortfalls": dict(sorted(splits.minimum_shortfalls.items())),
        },
        "ligands": {"n_purification_artifacts": n_artifacts},
        # Pointers to the data files written alongside this manifest.
        "files": {
            "splits": dict(SPLIT_FILES),
            "test_by_class": dict(sorted(test_class_files.items())),
            "clusters": CLUSTERS_FILENAME,
            "ligand_classes": CLASSES_FILENAME,
            "ligand_tiers": TIERS_FILENAME,
        },
    }


# --------------------------------------------------------------------------- #
# Split data files (the actual lists of PDB ids + supporting maps)
# --------------------------------------------------------------------------- #
def _write_id_list(ids: list[str], path: Path) -> Path:
    """Write a JSON array of ids, one per line (compact yet grepable)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    body = ",\n".join(json.dumps(i) for i in ids)
    path.write_text(f"[\n{body}\n]\n" if ids else "[]\n", encoding="utf-8")
    return path


def write_split_files(splits, class_map: dict[str, dict], out_dir: str | Path) -> dict[str, Path]:
    """Write train/val/test id lists + per-class test lists + supporting maps.

    Returns a name->path map of everything written. Pure function of the inputs:
    ids are sorted, so output is byte-stable.
    """
    out = Path(out_dir)
    per_split: dict[str, list[str]] = {s: [] for s in SPLIT_FILES}
    for entry, s in splits.entry_split.items():
        per_split[s].append(entry)
    for s in per_split:
        per_split[s].sort()

    written: dict[str, Path] = {}
    for s, fname in SPLIT_FILES.items():
        written[s] = _write_id_list(per_split[s], out / fname)

    # Per-class test-id lists: test entries carrying each functional class.
    test_ids = per_split["test"]
    class_to_ids: dict[str, list[str]] = {}
    for eid in test_ids:
        for cls in class_map.get(eid, {}).get("classes", []):
            class_to_ids.setdefault(cls, []).append(eid)
    for cls, ids in class_to_ids.items():
        written[f"test:{cls}"] = _write_id_list(sorted(ids), out / TEST_SUBDIR / f"{cls}_test.json")

    return written


def write_clusters(entry_to_cluster: dict[str, str], out_dir: str | Path) -> Path:
    """entry_id -> component key, for cluster-balanced sampling."""
    doc = {
        "clusters_schema": "if-split/clusters@1",
        "entry_clusters": dict(sorted(entry_to_cluster.items())),
    }
    return _write_json(doc, Path(out_dir) / CLUSTERS_FILENAME, compact=True)


def read_clusters(path: str | Path) -> dict[str, str]:
    path = Path(path)
    if not path.exists():
        return {}
    return dict(json.loads(path.read_text(encoding="utf-8")).get("entry_clusters", {}))


def write_classes(class_map: dict[str, dict], out_dir: str | Path) -> Path:
    """entry_id -> functional class labels."""
    classes = {eid: info["classes"] for eid, info in sorted(class_map.items())}
    doc = {"classes_schema": "if-split/classes@1", "classes": classes}
    return _write_json(doc, Path(out_dir) / CLASSES_FILENAME, compact=True)


def read_classes(path: str | Path) -> dict[str, list[str]]:
    path = Path(path)
    if not path.exists():
        return {}
    return dict(json.loads(path.read_text(encoding="utf-8")).get("classes", {}))


def read_id_list(path: str | Path) -> list[str]:
    """Read a split id-list file (train.json etc.) into a list of ids."""
    path = Path(path)
    if not path.exists():
        return []
    return list(json.loads(path.read_text(encoding="utf-8")))


def write_manifest(manifest: dict[str, Any], out_dir: str | Path) -> Path:
    # Pretty-print: the manifest is now small (KB) provenance, meant to be read.
    return _write_json(manifest, Path(out_dir) / "manifest.json")


def read_manifest(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


# --------------------------------------------------------------------------- #
# Ligand-tier audit sidecar (bulky; off the load path)
# --------------------------------------------------------------------------- #
def build_tiers_doc(class_map: dict[str, dict]) -> dict[str, Any]:
    """Per-entry, per-component tier + reason — the curation audit trail.

    Pure function of ``class_map`` (the same input as the manifest), so it stays
    deterministic and byte-stable. Lives in its own file because it is large and
    read by nobody on the load path.
    """
    tiers = {eid: info.get("tiers", {}) for eid, info in sorted(class_map.items())}
    return {"tiers_schema": TIERS_SCHEMA, "tiers": tiers}


def write_tiers(doc: dict[str, Any], out_dir: str | Path) -> Path:
    return _write_json(doc, Path(out_dir) / TIERS_FILENAME, compact=True)


def read_tiers(path: str | Path) -> dict[str, dict]:
    """Load the tier map from a sidecar file, or {} if absent."""
    path = Path(path)
    if not path.exists():
        return {}
    doc = json.loads(path.read_text(encoding="utf-8"))
    return dict(doc.get("tiers", {}))


# --------------------------------------------------------------------------- #
# verify / stats commands
# --------------------------------------------------------------------------- #
def verify_lock(lock_path: str | Path, *, client=None) -> int:
    """Re-enumerate from a lock's embedded config and report drift.

    Returns a process exit code: 0 = reproduced exactly, 1 = drift detected.
    ``client`` is injectable for offline testing; production passes None.
    """
    from .enumerate import enumerate_candidates

    lock = read_lock(lock_path)
    if lock.get("lock_schema") != LOCK_SCHEMA:
        print(f"warning: unexpected lock_schema {lock.get('lock_schema')!r}")

    cfg = Config.model_validate(lock["config"])
    limit = (lock.get("selection") or {}).get("limit")
    locked = lock["candidates"]
    locked_ids = set(locked["entry_ids"])
    locked_sha = locked["sha256"]

    print(f"verifying {lock['dataset_version']} (config {cfg.config_hash()})")
    print(f"  locked: {locked['count']} entries, candidates sha256={locked_sha[:12]}...")

    with tempfile.TemporaryDirectory() as tmp:
        records, _, sha = enumerate_candidates(
            cfg, tmp, limit=limit, client=client, progress=lambda m: print(f"  {m}")
        )

    now_ids = {r.entry_id for r in records}
    added = sorted(now_ids - locked_ids)
    removed = sorted(locked_ids - now_ids)  # obsoleted / withdrawn

    if sha == locked_sha and not added and not removed:
        print(f"OK: reproduced exactly ({len(records)} entries, hashes match).")
        return 0

    print("DRIFT detected:")
    if sha != locked_sha:
        print(f"  candidates sha256 differs: now {sha[:12]}... vs locked {locked_sha[:12]}...")
    if removed:
        print(f"  {len(removed)} entries no longer present (obsoleted/withdrawn):")
        print(f"    {', '.join(removed[:20])}{' ...' if len(removed) > 20 else ''}")
    if added:
        print(f"  {len(added)} new entries match the snapshot filters:")
        print(f"    {', '.join(added[:20])}{' ...' if len(added) > 20 else ''}")
    if not added and not removed:
        print("  entry set unchanged, but per-entry metadata changed (see hash).")
    return 1


def summarize_manifest(manifest_path: str | Path) -> int:
    """`stats` command: print split sizes and per-class (functional) test counts."""
    m = read_manifest(manifest_path)
    print(f"{m['dataset_version']}  (config {m['config_hash']})")
    flt = m["filter"]
    print(
        f"  candidates: {m['candidates']['count']}  kept: {flt['kept']}  dropped: {flt['dropped']}"
    )
    for reason, n in flt["drop_counts"].items():
        print(f"    - {reason}: {n}")
    cl = m["clustering"]
    print(
        f"  clustering: {cl['backend']} @ {cl['identity']}%  "
        f"components={cl['n_clusters']} (from {cl.get('n_raw_clusters', '?')} raw)  "
        f"multichain={cl['multichain_entries']}"
    )
    sp = m["splits"]
    print("  splits (entries / components):")
    for s in ("train", "val", "test"):
        ec = sp["entry_counts"].get(s, 0)
        cc = sp["cluster_counts"].get(s, 0)
        print(f"    {s:5s}: {ec:>7} entries  {cc:>7} components")
    print("  test set by ligand class (functional tier):")
    for cls, n in sp["per_split_class_counts"].get("test", {}).items():
        print(f"    {cls}: {n}")
    amb = sp.get("per_split_ambiguous_counts", {}).get("test", {})
    if amb:
        print("  test set ambiguous (reported, not labelled):")
        for cls, n in amb.items():
            print(f"    {cls}: {n}")
    shortfalls = sp.get("test_minimum_shortfalls", {})
    if shortfalls:
        print("  test minimum shortfalls (floor exceeded available supply):")
        for cls, n in shortfalls.items():
            print(f"    {cls}: short by {n}")
    lig = m["ligands"]
    n_arts = lig.get("n_purification_artifacts", len(lig.get("purification_artifacts", [])))
    print(f"  His-tag/Ni purification artifacts flagged: {n_arts}")
    files = m.get("files", {})
    if files:
        sf = files.get("splits", {})
        tbc = files.get("test_by_class", {})
        print(f"  split files: {', '.join(sf.values())}")
        if tbc:
            print(f"  per-class test files: {', '.join(tbc.values())}")
    return 0
