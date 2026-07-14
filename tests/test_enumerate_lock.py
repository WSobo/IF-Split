"""Phase 2 tests: Stage 1 enumerate determinism + snapshot lock/verify (offline).

These use the FakeRcsbClient fixture, so no network is touched. A real
network round-trip is exercised by the (opt-in) test in test_integration.py.
"""

from __future__ import annotations

from pathlib import Path

from ifsplit.config import load_config
from ifsplit.enumerate import enumerate_candidates
from ifsplit.manifest import build_lock, read_lock, verify_lock, write_lock

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config" / "default.yaml"


def _cfg():
    return load_config(DEFAULT_CONFIG)


def test_enumerate_writes_candidates(tmp_path, fake_client):
    records, path, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    assert path.exists()
    assert {r.entry_id for r in records} == {"4HHB", "1A1F"}
    assert len(sha) == 64


def test_enumerate_is_deterministic(tmp_path, fake_client):
    _, _, sha1 = enumerate_candidates(_cfg(), tmp_path / "a", client=fake_client)
    _, _, sha2 = enumerate_candidates(_cfg(), tmp_path / "b", client=fake_client)
    assert sha1 == sha2
    # And byte-identical files.
    a = (tmp_path / "a" / "candidates.jsonl").read_bytes()
    b = (tmp_path / "b" / "candidates.jsonl").read_bytes()
    assert a == b


def test_limit_is_reproducible_first_n(tmp_path, fake_client):
    records, _, _ = enumerate_candidates(_cfg(), tmp_path, limit=1, client=fake_client)
    # Sorted ascending -> "1A1F" is the first entry.
    assert [r.entry_id for r in records] == ["1A1F"]


def test_lock_roundtrip(tmp_path, fake_client):
    records, _, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    lock = build_lock(
        _cfg(),
        entry_ids=[r.entry_id for r in records],
        candidates_sha256=sha,
        limit=None,
    )
    lock_path = write_lock(lock, tmp_path)
    back = read_lock(lock_path)
    assert back["candidates"]["sha256"] == sha
    assert back["candidates"]["entry_ids"] == ["1A1F", "4HHB"]
    assert back["config_hash"] == _cfg().config_hash()


def test_verify_no_drift_returns_zero(tmp_path, fake_client):
    records, _, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    lock = build_lock(
        _cfg(),
        entry_ids=[r.entry_id for r in records],
        candidates_sha256=sha,
        limit=None,
    )
    lock_path = write_lock(lock, tmp_path)
    assert verify_lock(lock_path, client=fake_client) == 0


def test_verify_warns_on_version_mismatch(tmp_path, fake_client, capsys):
    # Same candidates, but the lock was written by a different tool version. The
    # lock pins the candidate set (Stage 1), not the split labels, so a version
    # bump is a WARNING (the split may differ), not data drift: verify still
    # succeeds because the thing the lock actually pins reproduced exactly.
    records, _, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    lock = build_lock(
        _cfg(),
        entry_ids=[r.entry_id for r in records],
        candidates_sha256=sha,
        limit=None,
    )
    lock["if_split_version"] = "0.0.1-ancient"  # simulate an older build
    lock_path = write_lock(lock, tmp_path)

    assert verify_lock(lock_path, client=fake_client) == 0
    out = capsys.readouterr().out
    assert "WARNING" in out
    assert "version differs" in out
    # The candidate set still reproduced exactly — only the version drifted.
    assert "candidate set reproduced exactly" in out
    assert "DRIFT detected" not in out


def test_verify_detects_removed_entry(tmp_path, fake_client, sample_entries):
    records, _, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    lock = build_lock(
        _cfg(),
        entry_ids=[r.entry_id for r in records],
        candidates_sha256=sha,
        limit=None,
    )
    lock_path = write_lock(lock, tmp_path)

    # Simulate obsolescence: an entry withdrawn since the lock was written.
    from conftest import FakeRcsbClient

    shrunk = FakeRcsbClient({"1A1F": sample_entries["1A1F"]})
    assert verify_lock(lock_path, client=shrunk) == 1


def test_verify_detects_added_entry(tmp_path, fake_client, sample_entries):
    # Lock with only one entry, then verify against a client that has two.
    one = enumerate_candidates(_cfg(), tmp_path, limit=1, client=fake_client)
    records, _, sha = one
    lock = build_lock(
        _cfg(),
        entry_ids=[r.entry_id for r in records],
        candidates_sha256=sha,
        limit=None,  # no limit on re-verify -> the second entry shows as "added"
    )
    lock_path = write_lock(lock, tmp_path)
    assert verify_lock(lock_path, client=fake_client) == 1
