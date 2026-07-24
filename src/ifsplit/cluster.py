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

A component's canonical key is the lexicographically smallest raw-cluster key across
its members. (Each RCSB raw key is itself a min-entity-id; an unclustered chain's raw
key is ``singleton:<hash of its sequence>``.) Keying the split hash on a stable,
content-derived key rather than RCSB's volatile integer cluster id keeps assignments
stable as the dataset grows (PLAN.md §6). A chain RCSB does not cluster (a short peptide)
is keyed by a hash of its sequence: a fully-unclustered entry's chains always key and merge
(bounded — the component holds only entries that are entirely that sequence), while an
unclustered chain inside an otherwise-clustered entry adds a merge edge only when it is long
enough to be a real chain (``MIN_UNCLUSTERED_MERGE_LEN``), so a promiscuous short peptide
cannot fan out into a spurious mega-component. An identical unclustered sequence that keys
thus cannot straddle two splits (check_no_leakage, comparing raw keys, sees it).

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

import hashlib
from dataclasses import dataclass, field

from .config import Config
from .schema import CandidateRecord

SINGLETON_PREFIX = "singleton:"

# Minimum sequence length for an UNCLUSTERED chain *inside an otherwise-clustered entry*
# to add a merge edge (union its entry's component with every other entry carrying that
# exact sequence). Below this, the chain is a short peptide/tag RCSB does not cluster;
# letting it union would fan out — one promiscuous peptide would merge every host protein
# it appears in into a spurious mega-component (catastrophic under `hash`, where that
# component lands in a salt-chosen split). Gating on LENGTH — intrinsic to the sequence and
# therefore growth-stable — rather than on how many components a sequence touches (a
# snapshot-dependent count that would make the split input-dependent) is what keeps the
# rule deterministic. A *fully*-unclustered entry is exempt: its sequence hash is its only
# identity, and its component can only hold entries that are entirely that sequence (no
# fan-out). PROVISIONAL (2026-07-24): a targeted RCSB probe found 9-mer peptides unclustered
# and the smallest clustered protein (crambin) at 46 aa. Confirm the exact knee against the
# full snapshot with scripts/measure_unclustered_fanout.py before treating this as final.
MIN_UNCLUSTERED_MERGE_LEN = 40


def _seq_singleton_key(seq: str) -> str:
    """Component key for an unclustered protein chain, derived from its SEQUENCE.

    Keying on a sequence hash (not the entity id) makes two entries carrying an
    identical unclustered sequence collapse into one component, so an identical
    sequence cannot straddle two splits — and ``check_no_leakage`` (which compares
    raw keys) can actually see it. It is also more growth-stable than an entity-id
    key: the same sequence always yields the same key regardless of when or under
    what id it was deposited.
    """
    return SINGLETON_PREFIX + hashlib.blake2b(seq.encode("utf-8"), digest_size=16).hexdigest()


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
    # Fold-benchmark labels (cfg.fold_benchmark_method): per-entry (super)family
    # membership DECOUPLED from union-find — never merges components, so it can label
    # folds even on a fold-leaky split. Empty when fold_benchmark_method == "off".
    entry_fold_labels: dict[str, list[str]] = field(default_factory=dict)

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
    bench_method = cfg.fold_benchmark_method
    entry_raw: dict[str, list[str]] = {}
    multichain: list[str] = []
    unclustered: list[str] = []
    all_keys: set[str] = set(raw_key.values())
    family_raw: dict[str, set[str]] = {}  # structural family -> raw keys sharing it
    entry_families: dict[str, list[str]] = {}  # entry -> fold family keys (method on)
    entry_fold_labels: dict[str, list[str]] = {}  # entry -> fold labels (benchmark; no merge)
    for r in records:
        proteins = [e for e in r.polymer_entities if e.is_protein]
        if not proteins:
            continue  # defensive; Stage 3 already drops no-protein entries
        clustered = {raw_key[e.cluster_ids[level]] for e in proteins if level in e.cluster_ids}
        uncl = [e for e in proteins if level not in e.cluster_ids]
        if clustered:
            # Already identified by its clustered chain(s); an unclustered chain adds a
            # MERGE EDGE only if it is long enough to be a real chain, never a short peptide
            # (which would fan out — see MIN_UNCLUSTERED_MERGE_LEN).
            singletons = {
                _seq_singleton_key(e.seq) for e in uncl if len(e.seq) >= MIN_UNCLUSTERED_MERGE_LEN
            }
        else:
            # Fully unclustered: the sequence hash is the chain's ONLY identity, so it must
            # always key (gating here would leave a short-peptide-only entry with no key,
            # reopening the straddle bug). Merging is bounded (the component holds only
            # entries that are entirely this sequence), so there is nothing to gate.
            singletons = {_seq_singleton_key(e.seq) for e in uncl}
        keys = sorted(clustered | singletons)
        key_set = set(keys)
        all_keys.update(keys)
        if not clustered:
            unclustered.append(r.entry_id)  # every protein chain is unclustered
        if len(keys) > 1:
            multichain.append(r.entry_id)
        entry_raw[r.entry_id] = keys
        if method != "off":
            efams: set[str] = set()
            for e in proteins:
                fams = e.structural_families.get(method)
                if not fams:
                    continue
                efams.update(fams)
                if level in e.cluster_ids:
                    rk = raw_key[e.cluster_ids[level]]
                else:
                    sk = _seq_singleton_key(e.seq)
                    rk = (
                        sk if sk in key_set else keys[0]
                    )  # short uncl chain -> entry's own component
                for fam in fams:
                    family_raw.setdefault(fam, set()).add(rk)
            if efams:
                entry_families[r.entry_id] = sorted(efams)
        if bench_method != "off":
            bfams: set[str] = set()
            for e in proteins:
                fams = e.structural_families.get(bench_method)
                if fams:
                    bfams.update(fams)
            if bfams:
                entry_fold_labels[r.entry_id] = sorted(bfams)

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
        entry_fold_labels=dict(sorted(entry_fold_labels.items())),
    )
