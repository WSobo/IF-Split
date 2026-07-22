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

The ``balanced`` strategy (``split_strategy="balanced"``) exists because
per-component hashing balances *components*, not *entries*: with heavy-tailed
component sizes (a dominant fold under structural clustering, or the antibody
mega-cluster even in sequence-only mode) one component balloons a split. Balanced
caps dominant folds to train and fills val/test to their *entry* targets from the
tail of smaller folds in hash order — restoring ~80/10/10 by entries with diverse,
fold-honest val/test sets. It stays leakage-safe (whole components) and
growth-stable via the registry (an in-place balanced rebuild auto-adopts
``<out>/splits.registry.json`` when the config matches; ``--fresh`` opts out), and
reports a gap if the tail was too thin.

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

# A component holding more than this fraction of all entries is a "dominant fold"
# (antibodies, TIM barrels): in the "balanced" strategy it is sent to train rather
# than allowed to balloon a val/test split. Small enough that no single component
# can overshoot a 10% val/test target by much; large enough to leave a rich tail.
BALANCE_MAX_COMPONENT_FRAC = 0.002


def bucket(key: str, salt: str) -> float:
    """Uniform float in [0, 1) from a stable hash of ``salt:key``."""
    digest = hashlib.blake2b(f"{salt}:{key}".encode(), digest_size=8).digest()
    return int.from_bytes(digest, "big") / 2**64


def split_fingerprint(entry_split: dict[str, str]) -> str:
    """Content hash of the ``entry_id -> split`` partition — the split *output*.

    ``sha256`` over sorted ``"<entry_id>\\t<split>\\n"`` lines. This is the actual
    deliverable (train/val/test membership), and it is pure ASCII: no float
    formatting, JSON key ordering, or dict-insertion order to drift across
    platforms — so the same partition always hashes identically. ``verify``
    compares this to prove the split reproduced, not just its Stage-1 inputs.
    Deliberately NOT the manifest (which stores only counts — two different
    partitions can share counts).
    """
    body = "".join(f"{eid}\t{s}\n" for eid, s in sorted(entry_split.items()))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


def registry_fingerprint(registry: dict[str, str]) -> str | None:
    """Hash of the split registry used at build time, or ``None`` if none was used.

    The split is a function of (candidates, config, code, **registry**); the lock
    pins the first three. Recording this lets ``verify`` distinguish a registry-free
    build — which it can reproduce and certify — from a ``--registry`` build, which
    it cannot without the same registry, so it reports *"not certified"* instead of
    a false split drift.
    """
    if not registry:
        return None
    body = "".join(f"{k}\t{v}\n" for k, v in sorted(registry.items()))
    return hashlib.sha256(body.encode("utf-8")).hexdigest()


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
    # "balanced" strategy diagnostics. strategy echoes cfg.split_strategy;
    # capped_folds is how many dominant folds were pinned to train; balance_gaps is
    # {split: entries below target} when the fold tail was too thin to fill val/test
    # (a signal the structural method is too aggressive — reported, never forced).
    strategy: str = "hash"
    capped_folds: int = 0
    balance_gaps: dict[str, int] = field(default_factory=dict)


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


def _balanced_assign(
    clusters: ClusterResult, cfg: Config, registry: dict[str, str]
) -> tuple[dict[str, str], int, dict[str, int]]:
    """Cap dominant folds to train; fill val/test to ENTRY targets from the tail.

    Leakage-safe (whole components move, never entries), deterministic (hash-ordered
    fill), and growth-stable via ``registry`` (pinned components keep their split;
    only new components fill the remaining budget). Returns
    ``(cluster_split, n_capped, gaps)`` where ``gaps`` records any val/test entry
    target the fold tail was too thin to reach.
    """
    sizes = {k: len(v) for k, v in clusters.cluster_members.items()}
    n_entries = sum(sizes.values())
    sf = cfg.split_fractions
    targets = {"val": sf.val * n_entries, "test": sf.test * n_entries}
    cap = BALANCE_MAX_COMPONENT_FRAC * n_entries

    cluster_split: dict[str, str] = {}
    totals = {"val": 0, "test": 0}
    n_capped = 0
    for key in clusters.cluster_members:
        if key in registry:  # pinned by a prior build (growth stability)
            s = registry[key]
            cluster_split[key] = s
            if s in totals:
                totals[s] += sizes[key]
        elif sizes[key] > cap:
            cluster_split[key] = "train"  # dominant fold -> train by design
            n_capped += 1

    eligible = sorted(
        (k for k in clusters.cluster_members if k not in cluster_split),
        key=lambda k: (bucket(k, cfg.split_salt), k),
    )
    for key in eligible:
        if totals["test"] < targets["test"]:
            cluster_split[key] = "test"
            totals["test"] += sizes[key]
        elif totals["val"] < targets["val"]:
            cluster_split[key] = "val"
            totals["val"] += sizes[key]
        else:
            cluster_split[key] = "train"

    gaps = {s: int(targets[s] - totals[s]) for s in ("val", "test") if totals[s] + 0.5 < targets[s]}
    return cluster_split, n_capped, gaps


def assign_splits(
    clusters: ClusterResult,
    cfg: Config,
    registry: dict[str, str] | None = None,
    entry_classes: dict[str, list[str]] | None = None,
) -> SplitResult:
    """Assign every component (and thus every entry) to a split.

    ``cfg.split_strategy`` selects "hash" (per-component hash onto the fractions) or
    "balanced" (cap dominant folds to train, fill val/test to entry targets from the
    fold tail; see :func:`_balanced_assign`). With ``cfg.test_min_per_class`` set and
    ``entry_classes`` provided, a deterministic top-up then recruits whole components
    into test to meet per-class floors (see :func:`_enforce_test_minimums`).
    """
    registry = registry or {}

    n_capped = 0
    balance_gaps: dict[str, int] = {}
    if cfg.split_strategy == "balanced":
        cluster_split, n_capped, balance_gaps = _balanced_assign(clusters, cfg, registry)
    else:
        cluster_split = {
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
        strategy=cfg.split_strategy,
        capped_folds=n_capped,
        balance_gaps=dict(sorted(balance_gaps.items())),
    )


def check_no_leakage(result: SplitResult, clusters: ClusterResult) -> None:
    """Verify no sequence cluster — and, with structural clustering on, no fold
    (super)family — spans two splits. Raises ``AssertionError`` on violation.

    Genuine check (not a tautology): for every entry, every *raw* sequence cluster
    it touches must resolve to the entry's own split, and (when
    ``structural_clustering`` is on) so must every structural (super)family it
    carries. Union-find guarantees both, so a failure means a real bug upstream.
    The fold check is a no-op for sequence-only builds (``entry_families`` empty),
    where two distinct folds may legitimately land in different splits.
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

    # Fold-level: no homologous (super)family may straddle two splits either — the
    # guarantee that makes structural_clustering meaningful. Empty (no-op) when off.
    fam_to_split: dict[str, str] = {}
    for entry, fams in clusters.entry_families.items():
        split = result.entry_split[entry]
        for fam in fams:
            prior = fam_to_split.setdefault(fam, split)
            if prior != split:
                raise AssertionError(
                    f"fold leakage: family {fam!r} appears in both {prior!r} and {split!r}"
                )
