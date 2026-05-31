"""Stage 7 - Snapshot lock + manifest, and the verify/stats commands.

Phase 2 implements the snapshot **lock** and `verify`:

- ``dataset.lock`` (JSON) records the two things needed to reproduce the
  candidate set: the embedded config (so verify is self-contained) and the
  canonical ``candidates.jsonl`` hash + entry-id list. Later phases extend it
  with the locked cluster file and the split assignment.
- ``verify`` re-enumerates from the embedded config and reports drift
  (added/removed entries, hash match), warning rather than failing so
  reproductions are honest about what changed in the live PDB.

The human-facing ``manifest.json`` and ``stats`` land in Phase 6.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from . import __version__
from .config import Config

LOCK_SCHEMA = "if-split/lock@1"


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


def write_lock(lock: dict[str, Any], out_dir: str | Path) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "dataset.lock"
    # Pretty-printed but stable: the candidates.sha256 is the integrity anchor,
    # not the lock file's own bytes.
    path.write_text(json.dumps(lock, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return path


def read_lock(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Lock file not found: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def verify_lock(lock_path: str | Path, *, client=None) -> int:
    """Re-enumerate from a lock's embedded config and report drift.

    Returns a process exit code: 0 = reproduced exactly, 1 = drift detected.
    ``client`` is injectable for offline testing; production passes None.
    """
    # Lazy import avoids an enumerate <-> manifest cycle.
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

    import tempfile

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


def summarize_manifest(manifest_path, *args, **kwargs):
    raise NotImplementedError("Stage 7 stats lands in Phase 6.")
