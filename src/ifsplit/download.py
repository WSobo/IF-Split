"""Stage 2 - Download mmCIF (biological assembly 1 when configured).

Fetches from RCSB file-download services, caches by ID, stores a SHA-256 per
file, and is resumable from the cache. Lands in Phase 2.
"""

from __future__ import annotations

from .config import Config


def download_entries(cfg: Config, *args, **kwargs):
    raise NotImplementedError("Stage 2 (download) lands in Phase 2.")
