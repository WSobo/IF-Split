"""Stage 6 - Deterministic cluster -> split assignment (the reproducibility core).

Each cluster's split is decided by ``blake2b(salt + ':' + canonical_key)`` mapped
onto the cumulative ``split_fractions``. Same salt + same key -> same split,
forever, independent of how many other clusters exist - so a larger snapshot only
*adds* clusters and never moves existing ones.

An optional ``registry`` (canonical_key -> split) pins prior assignments: if a
key is already in the registry its recorded split wins over the hash. Phase 7
persists this so growth is provably stable even in the edge case where a
smaller-id member later changes a cluster's canonical key.

Leakage: by construction each cluster maps to exactly one split (hard invariant,
asserted). Because an entry is assigned via its longest chain, a *secondary*
chain may belong to a cluster placed in another split; that residual
cross-split sequence overlap is measured and reported (not silently ignored).
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field

from .cluster import ClusterResult
from .config import Config

SPLITS = ("train", "val", "test")


def bucket(key: str, salt: str) -> float:
    """Uniform float in [0, 1) from a stable hash of ``salt:key``."""
    digest = hashlib.blake2b(f"{salt}:{key}".encode(), digest_size=8).digest()
    return int.from_bytes(digest, "big") / 2**64


def split_for_key(key: str, cfg: Config) -> str:
    """Map a cluster key to a split via the cumulative fractions."""
    b = bucket(key, cfg.split_salt)
    sf = cfg.split_fractions
    if b < sf.train:
        return "train"
    if b < sf.train + sf.val:
        return "val"
    return "test"


@dataclass
class SplitResult:
    cluster_split: dict[str, str]  # canonical key -> split
    entry_split: dict[str, str]  # entry_id -> split
    counts: dict[str, int]  # split -> entry count
    cluster_counts: dict[str, int]  # split -> cluster count
    leakage_entries: list[str] = field(default_factory=list)  # secondary-chain overlap


def assign_splits(
    clusters: ClusterResult, cfg: Config, registry: dict[str, str] | None = None
) -> SplitResult:
    """Assign every cluster (and thus every entry) to a split."""
    registry = registry or {}

    cluster_split: dict[str, str] = {}
    for key in clusters.cluster_members:
        cluster_split[key] = registry.get(key) or split_for_key(key, cfg)

    # Hard invariant: a cluster key resolves to exactly one split. (True by
    # construction - dict is single-valued - but assert the mapping is total.)
    missing = [k for k in clusters.cluster_members if k not in cluster_split]
    if missing:
        raise AssertionError(f"clusters without a split: {missing[:5]}")

    entry_split: dict[str, str] = {}
    counts = {s: 0 for s in SPLITS}
    for entry, key in clusters.entry_to_cluster.items():
        s = cluster_split[key]
        entry_split[entry] = s
        counts[s] += 1

    cluster_counts = {s: 0 for s in SPLITS}
    for s in cluster_split.values():
        cluster_counts[s] += 1

    # Residual leakage audit: entries whose secondary-chain clusters land in a
    # split other than the entry's own.
    leakage: list[str] = []
    for entry, keys in clusters.entry_all_clusters.items():
        own = entry_split[entry]
        if any(cluster_split.get(k, own) != own for k in keys):
            leakage.append(entry)

    return SplitResult(
        cluster_split=dict(sorted(cluster_split.items())),
        entry_split=dict(sorted(entry_split.items())),
        counts=counts,
        cluster_counts=cluster_counts,
        leakage_entries=sorted(leakage),
    )


def assert_no_cluster_leakage(result: SplitResult, clusters: ClusterResult) -> None:
    """Fail loudly if any cluster maps to more than one split.

    Verified by checking every entry's split equals its assigned cluster's split.
    """
    for entry, key in clusters.entry_to_cluster.items():
        expected = result.cluster_split[key]
        actual = result.entry_split[entry]
        if expected != actual:
            raise AssertionError(
                f"leakage: entry {entry} in {actual} but its cluster {key} is {expected}"
            )
