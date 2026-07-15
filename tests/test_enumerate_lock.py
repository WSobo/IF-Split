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


def test_enumerate_warns_on_resolutionless_method(tmp_path, fake_client):
    # Adding a method RCSB gives no resolution (NMR) under the resolution predicate
    # silently yields no such entries -> warn loudly instead of failing silently.
    import pytest

    cfg = _cfg().model_copy(update={"experimental_methods": ["X-RAY DIFFRACTION", "SOLUTION NMR"]})
    with pytest.warns(UserWarning, match="no resolution"):
        enumerate_candidates(cfg, tmp_path, client=fake_client)


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


# --------------------------------------------------------------------------- #
# Split-output certification (@2 locks, issue #5)
# --------------------------------------------------------------------------- #
def _lock_with_split(cfg, records, sha, *, registry_sha256=None):
    """Build a lock that pins the split output, by running Stages 3-6 offline."""
    from ifsplit.cluster import build_clusters
    from ifsplit.ligands import classify_components
    from ifsplit.parse import filter_candidates
    from ifsplit.split import assign_splits, split_fingerprint

    kept, _ = filter_candidates(records, cfg)
    class_map = {r.entry_id: classify_components(r, cfg) for r in kept}
    clusters = build_clusters(kept, cfg)
    splits = assign_splits(
        clusters, cfg, entry_classes={e: i["classes"] for e, i in class_map.items()}
    )
    return build_lock(
        cfg,
        entry_ids=[r.entry_id for r in records],
        candidates_sha256=sha,
        limit=None,
        split_sha256=split_fingerprint(splits.entry_split),
        registry_sha256=registry_sha256,
        split_strategy=cfg.split_strategy,
    )


def test_verify_certifies_reproduced_split(tmp_path, fake_client, capsys):
    # Candidates AND split reproduce -> verify certifies the split output.
    records, _, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    lock_path = write_lock(_lock_with_split(_cfg(), records, sha), tmp_path)
    assert verify_lock(lock_path, client=fake_client) == 0
    assert "split verified" in capsys.readouterr().out


# --------------- offline verify + resplit (no RCSB) ------------------------ #
def test_verify_offline_against_local_candidates(tmp_path, fake_client, capsys):
    # A distributed candidates.jsonl + lock can be verified with NO network.
    records, cand_path, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    lock_path = write_lock(_lock_with_split(_cfg(), records, sha), tmp_path)
    assert verify_lock(lock_path, candidates_path=cand_path) == 0
    out = capsys.readouterr().out
    assert "offline" in out
    assert "split verified" in out


def test_verify_offline_detects_wrong_candidates(tmp_path, fake_client):
    # A local candidates.jsonl that isn't the locked one -> drift (offline).
    records, _, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    lock_path = write_lock(_lock_with_split(_cfg(), records, sha), tmp_path)
    _, shrunk_path, _ = enumerate_candidates(_cfg(), tmp_path / "one", limit=1, client=fake_client)
    assert verify_lock(lock_path, candidates_path=shrunk_path) == 1


def test_resplit_warns_when_config_widens(tmp_path, fake_client, capsys):
    # resplit re-derives a FIXED snapshot; a config that would enumerate more must warn.
    from datetime import date

    from ifsplit.cli import _warn_if_config_would_reenumerate

    records, cand_path, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    write_lock(_lock_with_split(_cfg(), records, sha), tmp_path)  # lock beside candidates
    wider = _cfg().model_copy(update={"snapshot_date": date(2099, 1, 1)})
    _warn_if_config_would_reenumerate(wider, cand_path, sha)
    assert "WARNING" in capsys.readouterr().out


def test_resplit_notes_fixed_snapshot_without_lock(tmp_path, fake_client, capsys):
    from ifsplit.cli import _warn_if_config_would_reenumerate

    _records, cand_path, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    _warn_if_config_would_reenumerate(_cfg(), cand_path, sha)  # no lock beside candidates
    assert "note:" in capsys.readouterr().out


def test_resplit_guard_survives_malformed_lock(tmp_path, fake_client, capsys):
    # A dataset.lock of valid JSON but the wrong shape (a bare list) must degrade to
    # the generic caveat, never crash the purely-informational guard.
    from ifsplit.cli import _warn_if_config_would_reenumerate

    _records, cand_path, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    (tmp_path / "dataset.lock").write_text("[]", encoding="utf-8")
    _warn_if_config_would_reenumerate(_cfg(), cand_path, sha)  # must not raise
    assert "note:" in capsys.readouterr().out


def test_verify_resplit_lock_online_is_refused(tmp_path, fake_client, capsys):
    # A resplit-sourced lock can't be verified online (its config may not reproduce the
    # cached snapshot) -> verify steers to offline and returns 2, never a false drift.
    records, _, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    lock = _lock_with_split(_cfg(), records, sha)
    lock["source"] = "resplit"
    lock_path = write_lock(lock, tmp_path)
    assert verify_lock(lock_path, client=fake_client) == 2
    assert "--candidates" in capsys.readouterr().out


def test_verify_resplit_lock_offline_ok(tmp_path, fake_client):
    # ...but WITH --candidates the resplit lock verifies fine offline.
    records, cand_path, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    lock = _lock_with_split(_cfg(), records, sha)
    lock["source"] = "resplit"
    lock_path = write_lock(lock, tmp_path)
    assert verify_lock(lock_path, candidates_path=cand_path) == 0


def test_verify_offline_corrupt_candidates_is_integrity_failure(tmp_path, fake_client, capsys):
    # Offline verify is pitched as an integrity check: a corrupt candidates.jsonl must
    # be reported as an integrity FAILURE (exit 1), not an unrelated "invalid config".
    records, cand_path, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    lock_path = write_lock(_lock_with_split(_cfg(), records, sha), tmp_path)
    Path(cand_path).write_text('{"not a valid": "candidate record"}\n', encoding="utf-8")
    assert verify_lock(lock_path, candidates_path=cand_path) == 1
    assert "INTEGRITY CHECK FAILED" in capsys.readouterr().out


def test_verify_detects_split_output_drift(tmp_path, fake_client, capsys):
    # Same candidates, but the locked split hash doesn't match what the code now
    # produces (simulates a curation/split-logic change). Must be hard DRIFT.
    records, _, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    lock = _lock_with_split(_cfg(), records, sha)
    lock["split"]["sha256"] = "0" * 64  # a partition the current code won't produce
    lock_path = write_lock(lock, tmp_path)
    assert verify_lock(lock_path, client=fake_client) == 1
    assert "split output differs" in capsys.readouterr().out


def test_verify_split_certification_overrides_version_warning(tmp_path, fake_client, capsys):
    # A version bump no longer just warns: if the split is proven byte-identical,
    # verify certifies it (no scary "may differ" WARNING).
    records, _, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    lock = _lock_with_split(_cfg(), records, sha)
    lock["if_split_version"] = "0.0.1-ancient"
    lock_path = write_lock(lock, tmp_path)
    assert verify_lock(lock_path, client=fake_client) == 0
    out = capsys.readouterr().out
    assert "byte-identical" in out
    assert "WARNING" not in out


def test_verify_skips_split_when_registry_used(tmp_path, fake_client, capsys):
    # A --registry build can't be certified registry-blind -> report, don't false-drift.
    records, _, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    lock = _lock_with_split(_cfg(), records, sha, registry_sha256="deadbeef")
    lock_path = write_lock(lock, tmp_path)
    assert verify_lock(lock_path, client=fake_client) == 0
    assert "registry" in capsys.readouterr().out.lower()


def test_verify_candidate_drift_skips_split_check(tmp_path, fake_client, sample_entries, capsys):
    # When candidates drift (a grown/shrunk snapshot legitimately changes the
    # split), verify reports candidate drift and does NOT report split drift.
    records, _, sha = enumerate_candidates(_cfg(), tmp_path, client=fake_client)
    lock_path = write_lock(_lock_with_split(_cfg(), records, sha), tmp_path)
    from conftest import FakeRcsbClient

    shrunk = FakeRcsbClient({"1A1F": sample_entries["1A1F"]})
    assert verify_lock(lock_path, client=shrunk) == 1
    out = capsys.readouterr().out
    assert "DRIFT detected" in out
    assert "split output differs" not in out
