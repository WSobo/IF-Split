"""Opt-in network test: a real RCSB round-trip at tiny scale.

Skipped unless IFSPLIT_NETWORK_TESTS=1, so the default suite stays offline/fast.
Run with:  IFSPLIT_NETWORK_TESTS=1 uv run pytest -q tests/test_integration.py
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from ifsplit.config import load_config
from ifsplit.enumerate import enumerate_candidates
from ifsplit.manifest import build_lock, verify_lock, write_lock

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config" / "default.yaml"

pytestmark = pytest.mark.skipif(
    os.environ.get("IFSPLIT_NETWORK_TESTS") != "1",
    reason="set IFSPLIT_NETWORK_TESTS=1 to run network tests",
)


def test_small_build_and_verify_roundtrip(tmp_path):
    cfg = load_config(DEFAULT_CONFIG)
    records, path, sha = enumerate_candidates(cfg, tmp_path, limit=5)
    assert len(records) == 5
    assert path.exists()

    # Every record should carry real metadata pulled without coordinates.
    for r in records:
        assert r.release_date <= cfg.snapshot_date.isoformat()
        assert r.polymer_entities  # at least one polymer entity

    lock = build_lock(cfg, entry_ids=[r.entry_id for r in records], candidates_sha256=sha, limit=5)
    lock_path = write_lock(lock, tmp_path)
    # A back-to-back verify should reproduce exactly (no drift in seconds).
    assert verify_lock(lock_path) == 0
