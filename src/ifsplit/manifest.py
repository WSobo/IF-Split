"""Stage 7 - Snapshot lock, manifest, split registry, and verify/stats commands.

Artifacts written to the output dir:

- ``dataset.lock``   - reproduction anchor: embedded config + canonical
  ``candidates.jsonl`` hash + entry-id list, plus a ``split`` block hashing the
  entry->split partition (``@2`` locks). ``verify`` re-enumerates from it and reports
  candidate drift (added/removed entries, hash), and — when the candidate set
  reproduced byte-for-byte — recomputes the split and certifies the split *output*
  matches (registry-free builds), so a curation/split-logic change is caught even
  when the inputs are identical.
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

LOCK_SCHEMA = "if-split/lock@2"  # @2 adds the optional `split` block (split_sha256)
KNOWN_LOCK_SCHEMAS = frozenset({"if-split/lock@1", "if-split/lock@2"})
MANIFEST_SCHEMA = "if-split/manifest@2"
REGISTRY_SCHEMA = "if-split/registry@1"
TIERS_SCHEMA = "if-split/tiers@1"

# Data files written next to manifest.json (referenced by manifest["files"]).
TIERS_FILENAME = "ligands.tiers.json"
CLASSES_FILENAME = "ligands.classes.json"
CLUSTERS_FILENAME = "clusters.json"
TARGETS_FILENAME = "targets.jsonl"
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
    split_sha256: str | None = None,
    registry_sha256: str | None = None,
    split_strategy: str | None = None,
    source: str = "build",
) -> dict[str, Any]:
    """Assemble the lock document (pure; does not touch disk).

    When ``split_sha256`` is given, a ``split`` block is added so ``verify`` can
    certify the split *output* reproduced, not just the Stage-1 candidate inputs
    (see :func:`ifsplit.split.split_fingerprint`). Omitted for candidate-only locks
    (older ``@1`` locks and callers that build before the split exists).

    ``source`` records how the candidate set was produced: ``"build"`` (a live Stage-1
    enumeration — verifiable online) or ``"resplit"`` (re-derived from a cached
    ``candidates.jsonl`` whose config may differ in Stage-1 filters — so online
    re-enumeration would misreport drift; ``verify`` steers such locks to offline
    ``--candidates`` verification).
    """
    lock: dict[str, Any] = {
        "lock_schema": LOCK_SCHEMA,
        "dataset_version": cfg.dataset_version,
        "if_split_version": __version__,
        "config_hash": cfg.config_hash(),
        "config": cfg.canonical_dict(),
        "source": source,
        "selection": {"limit": limit},
        "candidates": {
            "count": len(entry_ids),
            "sha256": candidates_sha256,
            "entry_ids": sorted(entry_ids),
        },
    }
    if split_sha256 is not None:
        lock["split"] = {
            "sha256": split_sha256,
            "registry_sha256": registry_sha256,  # None for a registry-free build
            "strategy": split_strategy,
        }
    return lock


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

    # Training corpora: every kept structure is a backbone; the functional-tier
    # ligands are the conditioning targets (see targets.jsonl / build_targets).
    targets = build_targets(class_map, splits, clusters)
    functional_targets = [t for t in targets if t["tier"] == "functional"]

    def _target_counts(split):
        counts: dict[str, int] = {}
        for t in functional_targets:
            if t["split"] == split:
                counts[t["class"]] = counts.get(t["class"], 0) + 1
        return dict(sorted(counts.items()))

    training = {
        "n_backbones": len(splits.entry_split),
        "n_conditioning_targets": len(functional_targets),
        "targets_per_split_class": {s: _target_counts(s) for s in SPLITS},
        "n_optional_nonnative_targets": sum(1 for t in targets if t["tier"] == "ambiguous"),
    }

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
            "structural_method": clusters.structural_method,
            "n_seq_only_components": clusters.n_seq_only_components,
            "structural_bridging_families": clusters.n_structural_families,
        },
        "splits": {
            "strategy": splits.strategy,
            "capped_folds": splits.capped_folds,
            "balance_gaps": dict(sorted(splits.balance_gaps.items())),
            "entry_counts": dict(sorted(splits.counts.items())),
            "cluster_counts": dict(sorted(splits.cluster_counts.items())),
            "per_split_class_counts": per_split_class_counts,
            "per_split_ambiguous_counts": per_split_ambiguous_counts,
            "test_minimum_shortfalls": dict(sorted(splits.minimum_shortfalls.items())),
        },
        "ligands": {"n_purification_artifacts": n_artifacts},
        "training": training,
        # Pointers to the data files written alongside this manifest.
        "files": {
            "splits": dict(SPLIT_FILES),
            "test_by_class": dict(sorted(test_class_files.items())),
            "clusters": CLUSTERS_FILENAME,
            "ligand_classes": CLASSES_FILENAME,
            "ligand_tiers": TIERS_FILENAME,
            "targets": TARGETS_FILENAME,
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


# --------------------------------------------------------------------------- #
# Conditioning-target corpus (the ligand-conditioned training view)
# --------------------------------------------------------------------------- #
def build_targets(class_map: dict[str, dict], splits, clusters) -> list[dict]:
    """One row per (structure, conditioning target) for ligand-conditioned training.

    A *target* is a ligand a ligand-conditioned inverse-folding model should condition
    on. ``functional``-tier ligands (metal / small_molecule / nucleic_acid) are targets
    by default. A ``metal_site_nonnative`` site (a real metal pocket whose native metal
    isn't Ni/Co) and a ``glycan`` (a carbohydrate: glycosylation or a lectin ligand) are
    emitted as *opt-in* targets at tier ``ambiguous`` so a consumer can add them.
    Artifact / uncorroborated ligands are never
    targets -- their structures are still usable *backbones*, just with nothing to
    condition on. Finest grain (one row per target) so a consumer can condition on all
    of a structure's targets at once (group by entry_id) or one at a time. Pure and
    deterministic (sorted).
    """
    entry_split = splits.entry_split
    entry_cluster = clusters.entry_to_cluster

    def row(entry, cls, comp_id, tier, reason):
        return {
            "entry_id": entry,
            "split": entry_split[entry],
            "cluster": entry_cluster.get(entry, entry),
            "class": cls,
            "comp_id": comp_id,
            "tier": tier,
            "reason": reason,
        }

    rows: list[dict] = []
    for entry, info in class_map.items():
        if entry not in entry_split:
            continue
        tiers = info.get("tiers", {})
        for comp in info.get("metals", []):
            rows.append(
                row(entry, "metal", comp, "functional", tiers.get(comp, {}).get("reason", ""))
            )
        for comp in info.get("small_molecules", []):
            rows.append(
                row(
                    entry,
                    "small_molecule",
                    comp,
                    "functional",
                    tiers.get(comp, {}).get("reason", ""),
                )
            )
        if "nucleic_acid" in info.get("classes", []):
            reason = tiers.get("nucleic_acid", {}).get("reason", "protein_na_interface")
            rows.append(row(entry, "nucleic_acid", None, "functional", reason))
        # Opt-in ambiguous targets, recoverable if a consumer wants them: a real metal
        # site whose native metal isn't Ni/Co, or a carbohydrate (glycosylation vs a
        # genuine lectin/glycosidase ligand).
        for comp, t in tiers.items():
            reason = t.get("reason")
            if reason == "metal_site_nonnative":
                rows.append(row(entry, "metal", comp, "ambiguous", reason))
            elif reason == "glycan":
                rows.append(row(entry, "small_molecule", comp, "ambiguous", reason))

    rows.sort(key=lambda r: (r["entry_id"], r["class"], r["comp_id"] or "", r["tier"]))
    return rows


def write_targets(targets: list[dict], out_dir: str | Path) -> Path:
    """Write targets.jsonl (one compact JSON object per line; already sorted)."""
    path = Path(out_dir) / TARGETS_FILENAME
    path.parent.mkdir(parents=True, exist_ok=True)
    body = "".join(json.dumps(t, sort_keys=True, separators=(",", ":")) + "\n" for t in targets)
    path.write_text(body, encoding="utf-8")
    return path


def read_targets(path: str | Path) -> list[dict]:
    """Read a targets.jsonl conditioning corpus into a list of target rows."""
    path = Path(path)
    if not path.exists():
        return []
    out: list[dict] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out


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
def verify_lock(lock_path: str | Path, *, client=None, candidates_path=None) -> int:
    """Re-derive from a lock's embedded config and report drift.

    Sources the candidate set one of two ways and compares it to the lock: by
    re-enumerating Stage 1 from the live PDB (default), or — when ``candidates_path``
    is given — OFFLINE by hashing a local ``candidates.jsonl`` (a distributed dataset
    can then be integrity-checked with no network). If the candidate set reproduces
    byte-for-byte AND the lock records a ``split`` hash (``@2`` locks), it also
    recomputes Stages 3-6 and compares the split *output* — so a curation/split-logic
    change is caught even when the inputs are identical. Returns a process exit code:
    0 = reproduced exactly (split certified when possible), 1 = drift (candidate set
    differs, or the split output changed). ``client`` is injectable for offline testing.
    """
    lock = read_lock(lock_path)
    if lock.get("lock_schema") not in KNOWN_LOCK_SCHEMAS:
        print(f"warning: unexpected lock_schema {lock.get('lock_schema')!r}")

    cfg = Config.model_validate(lock["config"])
    limit = (lock.get("selection") or {}).get("limit")
    locked = lock["candidates"]
    locked_ids = set(locked["entry_ids"])
    locked_sha = locked["sha256"]
    locked_version = lock.get("if_split_version")
    version_match = locked_version == __version__
    locked_split = lock.get("split")  # None for @1 (candidate-only) locks

    print(f"verifying {lock['dataset_version']} (config {cfg.config_hash()})")
    print(f"  locked: {locked['count']} entries, candidates sha256={locked_sha[:12]}...")
    print(f"  if-split version: locked {locked_version}, running {__version__}")

    # A resplit lock pins a cached candidate set whose Stage-1 filters may differ from
    # this lock's config, so a live re-enumeration would misreport drift. Steer it to
    # offline verification (the candidate set it was built from) instead of misleading.
    if lock.get("source") == "resplit" and candidates_path is None:
        print(
            "  this lock was produced by `resplit` from a cached snapshot; its config's "
            "Stage-1 filters may not reproduce that snapshot, so an online re-enumerate would "
            "misreport drift. Verify it OFFLINE against the candidates.jsonl it was built from:\n"
            "    if-split verify <lock> --candidates <candidates.jsonl>"
        )
        return 2

    if candidates_path is not None:
        from .schema import read_candidates_jsonl, sha256_hex

        print(f"  offline: hashing local candidates {candidates_path}")
        sha = sha256_hex(Path(candidates_path).read_bytes())
        try:
            records = read_candidates_jsonl(candidates_path)
        except (ValueError, OSError) as exc:
            # A corrupt / truncated / non-canonical candidates.jsonl is exactly what an
            # offline integrity check should CATCH -> report as an integrity failure,
            # not let it surface as an unrelated "invalid config" error.
            print(f"  candidates file is corrupt or unparseable ({exc.__class__.__name__}): {exc}")
            print("INTEGRITY CHECK FAILED: candidates.jsonl could not be parsed.")
            return 1
    else:
        from .enumerate import enumerate_candidates

        with tempfile.TemporaryDirectory() as tmp:
            records, _, sha = enumerate_candidates(
                cfg, tmp, limit=limit, client=client, progress=lambda m: print(f"  {m}")
            )

    now_ids = {r.entry_id for r in records}
    added = sorted(now_ids - locked_ids)
    removed = sorted(locked_ids - now_ids)  # obsoleted / withdrawn
    candidates_match = sha == locked_sha and not added and not removed

    # Candidate drift: report and fail. The split is NOT compared here — a
    # grown/shrunk snapshot legitimately changes the split, so a split-hash check
    # only makes sense once the candidate set reproduced byte-for-byte.
    if not candidates_match:
        print("DRIFT detected:")
        if not version_match:
            print(
                f"  if-split version differs: locked {locked_version}, running {__version__} — "
                "curation/split logic may differ; install the locked version to reproduce exactly."
            )
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

    # Candidate set reproduced byte-for-byte. If the lock records a split hash,
    # certify the split OUTPUT directly by recomputing Stages 3-6 and comparing.
    if locked_split and locked_split.get("sha256"):
        return _verify_split(cfg, records, locked_split, version_match, locked_version)

    # Legacy @1 lock (no split hash): fall back to the version-string caveat.
    if version_match:
        print(f"OK: reproduced exactly ({len(records)} entries, hashes + version match).")
        return 0
    print(f"OK: candidate set reproduced exactly ({len(records)} entries, hashes match).")
    print(
        f"  WARNING: if-split version differs (locked {locked_version}, running "
        f"{__version__}). This lock predates split_sha256, so the split output is not "
        "certified; curation/split logic is version-specific and the rebuilt split may "
        "differ. Rebuild the lock with this version to certify the split, or install the "
        "locked version."
    )
    return 0


def _verify_split(cfg, records, locked_split, version_match, locked_version) -> int:
    """Recompute the split from reproduced candidates and compare its fingerprint.

    Called only when the candidate set already reproduced byte-for-byte, so any
    difference here is a genuine split-output change (curation/split logic drift),
    not snapshot growth. Recomputes registry-free — the build's conditions for a
    registry-free lock; a lock whose split used a registry is reported as not
    certifiable rather than compared against a registry-blind rebuild.
    """
    from .cluster import build_clusters
    from .ligands import classify_components
    from .parse import filter_candidates
    from .split import assign_splits, split_fingerprint

    if locked_split.get("registry_sha256") is not None:
        print(f"OK: candidate set reproduced exactly ({len(records)} entries).")
        print(
            "  NOTE: the build used a split registry, so the split output cannot be "
            "certified without it (candidate set verified). Re-run against a "
            "registry-free lock to certify the split."
        )
        return 0

    kept, _ = filter_candidates(records, cfg)
    class_map = {r.entry_id: classify_components(r, cfg) for r in kept}
    clusters = build_clusters(kept, cfg)
    splits = assign_splits(
        clusters, cfg, entry_classes={e: i["classes"] for e, i in class_map.items()}
    )
    now_split = split_fingerprint(splits.entry_split)
    locked_hash = locked_split["sha256"]

    if now_split == locked_hash:
        print(f"OK: reproduced exactly ({len(records)} entries; candidates + split verified).")
        if not version_match:
            print(
                f"  (if-split version differs: locked {locked_version}, running {__version__}; "
                "the split output is nonetheless byte-identical.)"
            )
        return 0

    print("DRIFT detected:")
    print(
        f"  split output differs: candidates are identical but the split hash changed "
        f"(now {now_split[:12]}... vs locked {locked_hash[:12]}...)."
    )
    if not version_match:
        print(
            f"  if-split version differs (locked {locked_version}, running {__version__}); "
            "the curation/split logic changed the assignment. Install the locked version "
            "to reproduce the split exactly."
        )
    else:
        print("  same tool version — an uncommitted or local change altered the split logic.")
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
    smethod = cl.get("structural_method", "off")
    if smethod != "off":
        seq_only = cl.get("n_seq_only_components", cl["n_clusters"])
        merged = seq_only - cl["n_clusters"]
        print(
            f"    fold-leakage control: {smethod}  "
            f"seq-only components={seq_only} -> {cl['n_clusters']} "
            f"({merged} folded by {cl.get('structural_bridging_families', 0)} shared superfamilies)"
        )
    sp = m["splits"]
    strat = sp.get("strategy", "hash")
    if strat != "hash":
        extra = f", {sp.get('capped_folds', 0)} dominant folds -> train"
        gaps = sp.get("balance_gaps") or {}
        if gaps:
            extra += f"; TAIL TOO THIN, val/test short by {gaps}"
        print(f"  split strategy: {strat}{extra}")
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
    tr = m.get("training")
    if tr:
        print("  training corpora:")
        print(f"    backbones (all structures):        {tr['n_backbones']}")
        print(f"    conditioning targets (functional): {tr['n_conditioning_targets']}")
        by = tr.get("targets_per_split_class", {})
        for s in ("train", "val", "test"):
            per = by.get(s, {})
            if per:
                print(f"      {s}: {', '.join(f'{c}={n}' for c, n in per.items())}")
        n_opt = tr.get("n_optional_nonnative_targets", 0)
        if n_opt:
            print(f"    opt-in non-native metal sites:     {n_opt}")
    files = m.get("files", {})
    if files:
        sf = files.get("splits", {})
        tbc = files.get("test_by_class", {})
        print(f"  split files: {', '.join(sf.values())}")
        if tbc:
            print(f"  per-class test files: {', '.join(tbc.values())}")
    return 0
