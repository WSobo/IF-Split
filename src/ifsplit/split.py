"""Stage 6 - Deterministic cluster -> split assignment (the reproducibility core).

Each cluster is assigned to a split by ``blake2b(canonical_cluster_id + salt)``
mapped onto the cumulative split fractions. The cluster id must be a stable,
input-independent representative (e.g. the lexicographically smallest member),
NOT whatever mmseqs2 happens to pick as the representative on a given run --
otherwise clusters can silently move splits as the dataset grows. The test set
is then stratified by ligand class. Invariant: no cluster spans two splits.
Lands in Phase 5.
"""

from __future__ import annotations

from .config import Config


def assign_splits(cfg: Config, *args, **kwargs):
    raise NotImplementedError("Stage 6 (split assignment) lands in Phase 5.")
