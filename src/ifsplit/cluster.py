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
unclustered chain inside an otherwise-clustered entry adds a merge edge only when it carries
enough MODELED sequence to be a real chain (``MIN_UNCLUSTERED_MERGE_MODELED``), so an
unmodeled (poly-'X') or low-complexity fragment cannot fan out into a spurious
mega-component. An identical unclustered sequence that keys thus cannot straddle two splits
(check_no_leakage, comparing raw keys, sees it).

**Exact-sequence identity is keyed independently of RCSB's clustering.** RCSB's cluster
file is *not* identity-complete: byte-identical sequences can be assigned DIFFERENT 30%
cluster ids (measured on the 2026-07-22 snapshot: 69 protein sequences across 497 entries,
including a 621-residue chain that landed in test *and* val, and a 532-residue chain in test
*and* train). A cluster id alone is therefore not a sound identity key, so every protein
chain carrying real modeled sequence also keys by its sequence hash. This is safe by
construction — exact identity is a strict subset of 30% identity, so the edge can only merge
what a correct 30% clustering would already have merged (measured: 19,593 -> 19,395
components, largest 43.2% -> 44.1%). Chains below ``MIN_UNCLUSTERED_MERGE_MODELED`` modeled
residues are excluded, so leakage of *identical* protein chains is eliminated above that
bound and explicitly not guaranteed below it.

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
from .schema import STRUCTURAL_METHODS, CandidateRecord

SINGLETON_PREFIX = "singleton:"

# An UNCLUSTERED chain *inside an otherwise-clustered entry* adds a merge edge (unions its
# entry's component with every other entry carrying that exact sequence) only if it has at
# least this many MODELED (non-'X') residues. Below it the chain is unmodeled (poly-'X') or
# a short/low-complexity fragment — RCSB does not cluster these, and they are exactly what
# fans out: a shared such chain would union every host protein it appears in into a spurious
# mega-component (catastrophic under `hash`, where it lands in a salt-chosen split). Gating
# on modeled sequence CONTENT is intrinsic → growth-stable (unlike a snapshot-dependent
# occurrence count, which would make the split input-dependent). A *fully*-unclustered entry
# is exempt: its component is bounded to entries that are entirely that sequence, and all-'X'
# fully-unclustered entries are already dropped in Stage 3.
#
# MEASURED (2026-07-22 full snapshot, scripts/measure_unclustered_fanout.py): among
# unclustered partial-entry chains the max merge fan-out collapses from 429 (all-'X') / 177
# (short real) to 2 at >= 12 modeled residues — a clean knee. RAW LENGTH does NOT separate
# them (a 72-'X' chain still bridges 283 clusters); the driver is unmodeled/low-complexity
# sequence, not length. See PLAN.md §5.
MIN_UNCLUSTERED_MERGE_MODELED = 12


def _modeled_len(seq: str) -> int:
    """Count of MODELED (non-'X') residues — the usable sequence content of a chain."""
    return sum(1 for c in seq if c != "X")


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
            # Already identified by its clustered chain(s); an unclustered chain adds a MERGE
            # EDGE only if it carries enough MODELED sequence to be a real chain — never an
            # unmodeled (poly-'X') or short/low-complexity fragment, which would fan out (see
            # MIN_UNCLUSTERED_MERGE_MODELED).
            singletons = {
                _seq_singleton_key(e.seq)
                for e in uncl
                if _modeled_len(e.seq) >= MIN_UNCLUSTERED_MERGE_MODELED
            }
        else:
            # Fully unclustered: the sequence hash is the chain's ONLY identity, so a real
            # chain must always key (gating here would leave a short-peptide-only entry with
            # no key, reopening the straddle bug) — merging is bounded (the component holds
            # only entries that are entirely this sequence). But skip poly-'X' (0-modeled)
            # chains: they carry no sequence and would merge unrelated entries through a
            # shared unmodeled trace. A kept fully-unclustered entry always has a non-poly-'X'
            # chain (Stage 3 drops all-'X' entries); the fallback is purely defensive.
            singletons = {_seq_singleton_key(e.seq) for e in uncl if _modeled_len(e.seq) > 0}
            if not singletons:
                singletons = {_seq_singleton_key(e.seq) for e in uncl}
        # RCSB's cluster file is not identity-complete (see the module docstring): identical
        # sequences can carry different 30% cluster ids, and a cluster id alone would then let
        # the SAME protein straddle two splits. Key every clustered chain by its sequence too,
        # so exact identity always merges regardless of what the cluster file says. Purely
        # additive, and bounded by the same modeled-length gate that guards fan-out.
        identity = {
            _seq_singleton_key(e.seq)
            for e in proteins
            if level in e.cluster_ids and _modeled_len(e.seq) >= MIN_UNCLUSTERED_MERGE_MODELED
        }
        # NB `multichain` counts entries that *bridge* clusters, so it is computed BEFORE the
        # identity keys are folded in — otherwise every single-chain entry would look bridging.
        bridging = clustered | singletons
        keys = sorted(bridging | identity)
        key_set = set(keys)
        all_keys.update(keys)
        if not clustered:
            unclustered.append(r.entry_id)  # every protein chain is unclustered
        if len(bridging) > 1:
            multichain.append(r.entry_id)
        entry_raw[r.entry_id] = keys
        if method != "off":
            # "union" merges on ANY of the three authorities (namespaced so a CATH code
            # can't collide with an ECOD/SCOP2 name) — the strictest metadata fold control.
            fam_methods = STRUCTURAL_METHODS if method == "union" else (method,)
            efams: set[str] = set()
            for e in proteins:
                if level in e.cluster_ids:
                    rk = raw_key[e.cluster_ids[level]]
                else:
                    sk = _seq_singleton_key(e.seq)
                    rk = sk if sk in key_set else keys[0]  # short uncl chain -> entry's component
                for fm in fam_methods:
                    for fam in e.structural_families.get(fm, []):
                        fam_key = f"{fm}:{fam}" if method == "union" else fam
                        efams.add(fam_key)
                        family_raw.setdefault(fam_key, set()).add(rk)
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
