"""Stage 5 - Cluster protein entities (default: RCSB precomputed membership).

The ``precomputed`` backend reads each protein entity's cluster id at the
configured identity level straight from ``PolymerEntity.cluster_ids`` (captured
in Stage 1 from the Data API ``rcsb_cluster_membership`` field) - no file
download, no mmseqs2 binary.

A *cluster* is the set of protein entities that share a cluster id at the level.
Its **canonical key** is the lexicographically smallest entity id among the
snapshot's members of that cluster. Keying the split hash on a member id (rather
than RCSB's volatile integer cluster id) keeps assignments stable as the dataset
grows (see PLAN.md §6; Phase 7 adds a registry for the residual edge case where
a smaller-id member joins later).

An entry may touch several clusters via different protein chains. We assign the
entry to the cluster of its **longest** protein chain (recording the others), so
each entry maps to exactly one cluster and the leakage invariant is well-defined.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .config import Config
from .schema import CandidateRecord

SINGLETON_PREFIX = "singleton:"


@dataclass
class ClusterResult:
    """Outcome of Stage 5 clustering at a single identity level."""

    identity: int
    entry_to_cluster: dict[str, str]  # entry_id -> canonical cluster key
    cluster_members: dict[str, list[str]]  # canonical key -> sorted entry_ids
    entry_all_clusters: dict[str, list[str]]  # entry_id -> all touched cluster keys
    multichain_entries: list[str] = field(default_factory=list)
    unclustered_entries: list[str] = field(default_factory=list)

    @property
    def n_clusters(self) -> int:
        return len(self.cluster_members)


def build_clusters(records: list[CandidateRecord], cfg: Config) -> ClusterResult:
    """Cluster the (already-filtered) records at ``cfg.identity_level``."""
    if cfg.clustering_backend != "precomputed":
        raise NotImplementedError(
            f"clustering_backend {cfg.clustering_backend!r} not implemented "
            "(only 'precomputed' is available)."
        )
    level = cfg.identity_level

    # 1. Group entities by RCSB cluster id -> their entity ids, to derive a
    #    stable canonical key (the smallest member id) per cluster.
    cluster_entities: dict[int, set[str]] = {}
    for r in records:
        for e in r.polymer_entities:
            if e.is_protein and level in e.cluster_ids:
                cluster_entities.setdefault(e.cluster_ids[level], set()).add(e.entity_id)
    canon: dict[int, str] = {cid: min(ents) for cid, ents in cluster_entities.items()}

    entry_to_cluster: dict[str, str] = {}
    entry_all_clusters: dict[str, list[str]] = {}
    members: dict[str, set[str]] = {}
    multichain: list[str] = []
    unclustered: list[str] = []

    for r in records:
        proteins = [e for e in r.polymer_entities if e.is_protein]
        if not proteins:
            continue  # defensive; Stage 3 already drops these

        keys = {canon[e.cluster_ids[level]] for e in proteins if level in e.cluster_ids}

        if not keys:
            # No protein chain has a cluster at this level (e.g. <10-aa peptides
            # that RCSB excludes). Treat as its own singleton cluster.
            key = SINGLETON_PREFIX + min(e.entity_id for e in proteins)
            entry_to_cluster[r.entry_id] = key
            entry_all_clusters[r.entry_id] = [key]
            members.setdefault(key, set()).add(r.entry_id)
            unclustered.append(r.entry_id)
            continue

        if len(keys) > 1:
            multichain.append(r.entry_id)

        # Assign to the cluster of the longest protein chain (tie-break on id).
        longest = max(proteins, key=lambda e: (e.seq_len, e.entity_id))
        key = canon[longest.cluster_ids[level]] if level in longest.cluster_ids else sorted(keys)[0]

        entry_to_cluster[r.entry_id] = key
        entry_all_clusters[r.entry_id] = sorted(keys)
        members.setdefault(key, set()).add(r.entry_id)

    return ClusterResult(
        identity=level,
        entry_to_cluster=entry_to_cluster,
        cluster_members={k: sorted(v) for k, v in sorted(members.items())},
        entry_all_clusters=dict(sorted(entry_all_clusters.items())),
        multichain_entries=sorted(multichain),
        unclustered_entries=sorted(unclustered),
    )
