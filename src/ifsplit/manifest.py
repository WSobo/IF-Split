"""Stage 7 - Snapshot lock, manifest, split registry, and verify/stats commands.

Artifacts written to the output dir:

- ``dataset.lock``   - reproduction anchor: embedded config + canonical
  ``candidates.jsonl`` hash + entry-id list. ``verify`` re-enumerates from it and
  reports drift (added/removed entries, hash match), warning not failing.
- ``manifest.json``  - human-facing run record: config, drop log, per-split entry
  lists, ligand-class tags, per-class test counts, cluster/leakage stats. Built
  as a pure function of its inputs (no wall-clock fields) so two runs of the same
  config produce byte-identical manifests (Phase 7).
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
MANIFEST_SCHEMA = "if-split/manifest@1"
REGISTRY_SCHEMA = "if-split/registry@1"


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


def _write_json(obj: dict[str, Any], path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n", encoding="utf-8")
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

    # Ligand-class tags per entry + per-class counts within each split.
    ligand_classes = {eid: info["classes"] for eid, info in sorted(class_map.items())}
    purification_artifacts = sorted(
        eid for eid, info in class_map.items() if info.get("purification_artifact")
    )

    per_split_class_counts: dict[str, dict[str, int]] = {}
    for s in SPLITS:
        counts: dict[str, int] = {}
        for entry in per_split[s]:
            for cls in class_map.get(entry, {}).get("classes", []):
                counts[cls] = counts.get(cls, 0) + 1
        per_split_class_counts[s] = dict(sorted(counts.items()))

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
            "multichain_entries": len(clusters.multichain_entries),
            "unclustered_entries": len(clusters.unclustered_entries),
        },
        "splits": {
            "entry_counts": dict(sorted(splits.counts.items())),
            "cluster_counts": dict(sorted(splits.cluster_counts.items())),
            "leakage_entries": len(splits.leakage_entries),
            "per_split_class_counts": per_split_class_counts,
            "entries": per_split,
        },
        "ligands": {
            "classes": ligand_classes,
            "purification_artifacts": purification_artifacts,
        },
    }


def write_manifest(manifest: dict[str, Any], out_dir: str | Path) -> Path:
    return _write_json(manifest, Path(out_dir) / "manifest.json")


def read_manifest(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Manifest not found: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


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
    """`stats` command: print split sizes and per-class test counts."""
    m = read_manifest(manifest_path)
    print(f"{m['dataset_version']}  (config {m['config_hash']})")
    flt = m["filter"]
    print(
        f"  candidates: {m['candidates']['count']}  kept: {flt['kept']}  dropped: {flt['dropped']}"
    )
    if flt["drop_counts"]:
        for reason, n in flt["drop_counts"].items():
            print(f"    - {reason}: {n}")
    cl = m["clustering"]
    print(
        f"  clustering: {cl['backend']} @ {cl['identity']}%  "
        f"clusters={cl['n_clusters']}  multichain={cl['multichain_entries']}"
    )
    sp = m["splits"]
    print("  splits (entries / clusters):")
    for s in ("train", "val", "test"):
        ec = sp["entry_counts"].get(s, 0)
        cc = sp["cluster_counts"].get(s, 0)
        print(f"    {s:5s}: {ec:>7} entries  {cc:>7} clusters")
    print(f"  cross-split secondary-chain overlap: {sp['leakage_entries']} entries")
    print("  test set by ligand class:")
    for cls, n in sp["per_split_class_counts"].get("test", {}).items():
        print(f"    {cls}: {n}")
    arts = m["ligands"]["purification_artifacts"]
    print(f"  His-tag/Ni purification artifacts flagged: {len(arts)}")
    return 0
