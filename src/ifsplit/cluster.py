"""Stage 5 - Sequence clustering with mmseqs2.

Pools unique protein chain sequences and runs ``mmseqs easy-cluster`` at
``--min-seq-id identity_threshold`` with logged coverage flags. Inputs are
sorted and the mmseqs2 version is pinned/logged for determinism. Lands in
Phase 5.
"""

from __future__ import annotations

from .config import Config


def cluster_sequences(cfg: Config, *args, **kwargs):
    raise NotImplementedError("Stage 5 (mmseqs2 clustering) lands in Phase 5.")
