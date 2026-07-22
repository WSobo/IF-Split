"""Stage 5 - Cluster protein entities into leakage-safe groups.

Clustering reuses RCSB's published polymer-entity clusters: each protein entity's
RCSB cluster id at the configured identity level is read from
``PolymerEntity.cluster_ids`` (captured in Stage 1 from the Data API
``rcsb_cluster_membership`` field) - no file download, no external binary. These
are the same mmseqs2-computed 30% clusters ProteinMPNN/LigandMPNN used, but locked
via the snapshot so the split stays byte-for-byte reproducible.

A *raw cluster* is the set of protein entities sharing an RCSB cluster id. But an
entry with several protein chains can touch several raw clusters, so raw clusters
alone are NOT a leakage-safe split unit: if entry X has chain a (raw cluster A)
and chain b (raw cluster B), then A and B must land in the same split or X's b
sequence leaks across splits. We therefore merge raw clusters joined by a shared
entry into **components** (connected components, union-find). The component is the
unit Stage 6 assigns to a split, which makes cross-split sequence overlap
impossible by construction - no heuristic, no after-the-fact audit.

A component's canonical key is the lexicographically smallest entity id across all
its members. (Equivalently the smallest raw-cluster key, since each raw key is
itself a min-entity-id - so min-of-mins = global min.) Keying the split hash on a
stable member id, not RCSB's volatile integer cluster id, keeps assignments
stable as the dataset grows (PLAN.md §6). Sub-10-aa peptides that RCSB does not
cluster become their own singleton components.

**Fold-level leakage control (opt-in via ``structural_clustering``).** Sequence
clustering alone misses structural redundancy: two chains below the identity
threshold can be the same fold, which an inverse-folding model (structure ->
sequence) would leak across splits. When enabled, raw clusters whose entities
share a structural (super)family - from RCSB's precomputed CATH/ECOD/SCOP2
annotations, metadata only - are union-merged too, so the same fold cannot
straddle train/test. It is purely additive (only merges, never splits) and
degrades gracefully: a chain with no classification simply adds no structural
edge.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .config import Config
from .schema import CandidateRecord

SINGLETON_PREFIX = "singleton:"


@dataclass
class ClusterResult:
    """Stage 5 output: raw sequence clusters merged into leakage-safe components."""

    identity: int
    entry_to_cluster: dict[str, str]  # entry_id -> component key (the split unit)
    cluster_members: dict[str, list[str]]  # component key -> sorted entry_ids
    entry_raw_clusters: dict[str, list[str]]  # entry_id -> raw cluster keys it touches
    multichain_entries: list[str] = field(default_factory=list)
    unclustered_entries: list[str] = field(default_factory=list)
    n_raw_clusters: int = 0
    # Fold-level leakage control (Stage 5). "off" when disabled; otherwise the
    # classification method used. n_seq_only_components is the component count from
    # sequence edges alone, so (n_seq_only_components - n_clusters) is how many
    # components structural merging folded together. n_structural_families counts
    # the distinct (super)families that actually bridged >=2 sequence clusters.
    structural_method: str = "off"
    n_seq_only_components: int = 0
    n_structural_families: int = 0
    # Per-entry structural (super)family keys under the active method (empty when
    # structural_method == "off"). Lets check_no_leakage assert the fold-level
    # guarantee — no homologous (super)family straddles two splits — directly,
    # matching the fold-leakage claim (not just sequence-cluster leakage).
    entry_families: dict[str, list[str]] = field(default_factory=dict)

    @property
    def n_clusters(self) -> int:
        """Number of components (the split units)."""
        return len(self.cluster_members)


def build_clusters(records: list[CandidateRecord], cfg: Config) -> ClusterResult:
    """Cluster filtered records at ``cfg.identity_level``, merged into components."""
    level = cfg.identity_level

    # 1. Raw clusters: RCSB cluster id -> member entity ids -> canonical raw key
    #    (the smallest member entity id).
    raw_entities: dict[int, set[str]] = {}
    for r in records:
        for e in r.polymer_entities:
            if e.is_protein and level in e.cluster_ids:
                raw_entities.setdefault(e.cluster_ids[level], set()).add(e.entity_id)
    raw_key = {cid: min(ents) for cid, ents in raw_entities.items()}

    # 2. Each entry -> the raw cluster keys it touches (a singleton key if no
    #    protein chain is clustered at this level). Also record each protein
    #    entity's raw key, so structural families (step 3b) can be attached to the
    #    same union-find nodes the sequence clusters use.
    method = cfg.structural_clustering
    entry_raw: dict[str, list[str]] = {}
    multichain: list[str] = []
    unclustered: list[str] = []
    all_keys: set[str] = set(raw_key.values())
    family_raw: dict[str, set[str]] = {}  # structural family -> raw keys sharing it
    entry_families: dict[str, list[str]] = {}  # entry -> fold family keys (method on)
    for r in records:
        proteins = [e for e in r.polymer_entities if e.is_protein]
        if not proteins:
            continue  # defensive; Stage 3 already drops no-protein entries
        keys = sorted({raw_key[e.cluster_ids[level]] for e in proteins if level in e.cluster_ids})
        if not keys:
            singleton = SINGLETON_PREFIX + min(e.entity_id for e in proteins)
            keys = [singleton]
            all_keys.add(singleton)
            unclustered.append(r.entry_id)
        elif len(keys) > 1:
            multichain.append(r.entry_id)
        entry_raw[r.entry_id] = keys
        if method != "off":
            efams: set[str] = set()
            for e in proteins:
                fams = e.structural_families.get(method)
                if not fams:
                    continue
                efams.update(fams)
                rk = raw_key[e.cluster_ids[level]] if level in e.cluster_ids else keys[0]
                for fam in fams:
                    family_raw.setdefault(fam, set()).add(rk)
            if efams:
                entry_families[r.entry_id] = sorted(efams)

    # 3. Union-find: merge raw clusters joined by a shared entry into components.
    #    The smaller key is always made the root, so a component's root is its
    #    global-minimum key (order-independent -> deterministic).
    parent = {k: k for k in all_keys}

    def find(x: str) -> str:
        root = x
        while parent[root] != root:
            root = parent[root]
        while parent[x] != root:  # path compression
            parent[x], x = root, parent[x]
        return root

    def union(a: str, b: str) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            lo, hi = (ra, rb) if ra < rb else (rb, ra)
            parent[hi] = lo

    for keys in entry_raw.values():
        for k in keys[1:]:
            union(keys[0], k)

    # 3b. Structural edges: union raw clusters that share a fold (homologous
    #     superfamily). Sequence-only components are counted first so the manifest
    #     can report how many the structural pass folded together. Purely additive.
    n_seq_only = len({find(k) for k in all_keys})
    n_bridging_families = 0
    for fam in sorted(family_raw):
        rks = sorted(family_raw[fam])
        if len({find(k) for k in rks}) > 1:
            n_bridging_families += 1
        for k in rks[1:]:
            union(rks[0], k)

    # 4. Materialize components: component key -> entries; entry -> component.
    entry_to_cluster: dict[str, str] = {}
    members: dict[str, set[str]] = {}
    for entry, keys in entry_raw.items():
        comp = find(keys[0])
        entry_to_cluster[entry] = comp
        members.setdefault(comp, set()).add(entry)

    return ClusterResult(
        identity=level,
        entry_to_cluster=dict(sorted(entry_to_cluster.items())),
        cluster_members={k: sorted(v) for k, v in sorted(members.items())},
        entry_raw_clusters=dict(sorted(entry_raw.items())),
        multichain_entries=sorted(multichain),
        unclustered_entries=sorted(unclustered),
        n_raw_clusters=len(raw_key),
        structural_method=method,
        n_seq_only_components=n_seq_only,
        n_structural_families=n_bridging_families,
        entry_families=dict(sorted(entry_families.items())),
    )
