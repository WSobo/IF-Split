"""Stage 6 - Deterministic component -> split assignment (the reproducibility core).

Each *component* (a leakage-safe group of sequence clusters joined by shared
multi-chain entries; see cluster.py) is assigned to a split by
``blake2b(salt + ':' + component_key)`` mapped onto the cumulative
``split_fractions``. Same salt + same key -> same split, forever, independent of
how many other components exist - so a larger snapshot only *adds* components and
never moves existing ones.

An optional ``registry`` (component_key -> split) pins prior assignments: if a
key is already in the registry its recorded split wins over the hash, so growth
is stable even if a component's canonical key shifts (e.g. a smaller-id member
joins later).

**No-leakage is structural, not heuristic.** Because every entity an entry
touches lives in the same component (union-find merged them), and a component
maps to exactly one split, two splits can never share a sequence cluster. The
``check_no_leakage`` invariant re-derives this from the cluster membership and
fails loudly if it is ever violated - a real guard, not a tautology.
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
    """Map a component key to a split via the cumulative fractions."""
    b = bucket(key, cfg.split_salt)
    sf = cfg.split_fractions
    if b < sf.train:
        return "train"
    if b < sf.train + sf.val:
        return "val"
    return "test"


@dataclass
class SplitResult:
    cluster_split: dict[str, str]  # component key -> split
    entry_split: dict[str, str]  # entry_id -> split
    counts: dict[str, int]  # split -> entry count
    cluster_counts: dict[str, int]  # split -> component count
    # Per-class test floors that could not be fully met (class -> shortfall). Empty
    # when no minimums were requested or all were satisfied. Reported, never forced.
    minimum_shortfalls: dict[str, int] = field(default_factory=dict)


def _enforce_test_minimums(
    cluster_split: dict[str, str],
    clusters: ClusterResult,
    cfg: Config,
    entry_classes: dict[str, list[str]],
    registry: dict[str, str],
) -> tuple[dict[str, str], dict[str, int]]:
    """Recruit whole components into test until per-class floors are met.

    Leakage-safe (moves components, never entries), deterministic (hash-ordered),
    and growth-stable (never overrides a registry-pinned component). Returns the
    updated ``cluster_split`` and a ``{class: shortfall}`` map for any floor that
    could not be fully satisfied from the available supply.
    """
    minimums = {c: n for c, n in cfg.test_min_per_class.items() if n > 0}
    if not minimums:
        return cluster_split, {}

    cluster_split = dict(cluster_split)
    # Per-component count of entries carrying each class (so test totals update O(1)).
    comp_class_counts: dict[str, dict[str, int]] = {}
    for key, entries in clusters.cluster_members.items():
        counts: dict[str, int] = {}
        for e in entries:
            for cls in entry_classes.get(e, []):
                counts[cls] = counts.get(cls, 0) + 1
        comp_class_counts[key] = counts

    test_totals: dict[str, int] = {}
    for key, split in cluster_split.items():
        if split == "test":
            for cls, n in comp_class_counts[key].items():
                test_totals[cls] = test_totals.get(cls, 0) + n

    shortfalls: dict[str, int] = {}
    for cls in sorted(minimums):
        need = minimums[cls]
        eligible = [
            key
            for key in clusters.cluster_members
            if cluster_split[key] != "test"
            and comp_class_counts[key].get(cls, 0) > 0
            and key not in registry  # respect pinned assignments (growth stability)
        ]
        eligible.sort(key=lambda k: (bucket(k, cfg.split_salt), k))
        for key in eligible:
            if test_totals.get(cls, 0) >= need:
                break
            cluster_split[key] = "test"
            for c, n in comp_class_counts[key].items():
                test_totals[c] = test_totals.get(c, 0) + n
        deficit = need - test_totals.get(cls, 0)
        if deficit > 0:
            shortfalls[cls] = deficit
    return cluster_split, shortfalls


def assign_splits(
    clusters: ClusterResult,
    cfg: Config,
    registry: dict[str, str] | None = None,
    entry_classes: dict[str, list[str]] | None = None,
) -> SplitResult:
    """Assign every component (and thus every entry) to a split.

    With ``cfg.test_min_per_class`` set and ``entry_classes`` provided, a
    deterministic top-up recruits whole components into test to meet per-class
    floors (see :func:`_enforce_test_minimums`).
    """
    registry = registry or {}

    cluster_split: dict[str, str] = {
        key: registry.get(key, split_for_key(key, cfg)) for key in clusters.cluster_members
    }

    shortfalls: dict[str, int] = {}
    if cfg.test_min_per_class and entry_classes is not None:
        cluster_split, shortfalls = _enforce_test_minimums(
            cluster_split, clusters, cfg, entry_classes, registry
        )

    counts = {s: 0 for s in SPLITS}
    entry_split: dict[str, str] = {}
    for entry, key in clusters.entry_to_cluster.items():
        s = cluster_split[key]
        entry_split[entry] = s
        counts[s] += 1

    cluster_counts = {s: 0 for s in SPLITS}
    for s in cluster_split.values():
        cluster_counts[s] += 1

    return SplitResult(
        cluster_split=dict(sorted(cluster_split.items())),
        entry_split=dict(sorted(entry_split.items())),
        counts=counts,
        cluster_counts=cluster_counts,
        minimum_shortfalls=dict(sorted(shortfalls.items())),
    )


def check_no_leakage(result: SplitResult, clusters: ClusterResult) -> None:
    """Verify no sequence cluster spans two splits. Raises on violation.

    Genuine check (not a tautology): for every entry, every *raw* sequence
    cluster it touches must resolve to the entry's own split. Union-find
    guarantees this, so a failure means a real bug upstream.
    """
    raw_to_split: dict[str, str] = {}
    for entry, raw_keys in clusters.entry_raw_clusters.items():
        split = result.entry_split[entry]
        for rk in raw_keys:
            prior = raw_to_split.setdefault(rk, split)
            if prior != split:
                raise AssertionError(
                    f"leakage: raw cluster {rk} appears in both {prior!r} and {split!r}"
                )
