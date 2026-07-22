"""CLI exit-code + error-path tests (offline; no network).

Exercise the documented exit codes and the actionable error messages `main()`
emits for bad input — the surface a researcher hits first when a file is missing,
malformed, old-schema, or an argument is out of range.
"""

from __future__ import annotations

from pathlib import Path

from ifsplit.cli import EXIT_BAD_INPUT, main

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config" / "default.yaml"


def test_build_missing_config_is_bad_input(tmp_path, capsys):
    rc = main(["build", "--config", str(tmp_path / "nope.yaml"), "--out", str(tmp_path / "o")])
    assert rc == EXIT_BAD_INPUT
    assert "not found" in capsys.readouterr().err.lower()


def test_build_invalid_config_is_bad_input(tmp_path, capsys):
    # 0.40 -> 40% is not an RCSB precomputed cluster level -> ValidationError (exit 2),
    # caught before any network access.
    bad = tmp_path / "bad.yaml"
    bad.write_text(
        DEFAULT_CONFIG.read_text().replace("identity_threshold: 0.30", "identity_threshold: 0.40")
    )
    rc = main(["build", "--config", str(bad), "--out", str(tmp_path / "o")])
    assert rc == EXIT_BAD_INPUT
    assert "invalid config" in capsys.readouterr().err.lower()


def test_stats_missing_manifest_is_bad_input(tmp_path):
    assert main(["stats", str(tmp_path / "missing.json")]) == EXIT_BAD_INPUT


def test_stats_malformed_json_is_bad_input(tmp_path, capsys):
    m = tmp_path / "manifest.json"
    m.write_text("{ this is not valid json")
    rc = main(["stats", str(m)])
    assert rc == EXIT_BAD_INPUT
    assert "malformed json" in capsys.readouterr().err.lower()


def test_stats_old_schema_manifest_is_bad_input(tmp_path):
    # Valid JSON, but missing every expected field: a clean bad-input exit, not a crash.
    m = tmp_path / "manifest.json"
    m.write_text("{}")
    assert main(["stats", str(m)]) == EXIT_BAD_INPUT


def test_verify_missing_lock_is_bad_input(tmp_path):
    assert main(["verify", str(tmp_path / "missing.lock")]) == EXIT_BAD_INPUT


def test_fetch_zero_workers_is_bad_input(tmp_path, capsys):
    m = tmp_path / "manifest.json"
    m.write_text("{}")
    rc = main(["fetch", str(m), "--all", "--workers", "0", "--out", str(tmp_path / "s")])
    assert rc == EXIT_BAD_INPUT
    assert "workers" in capsys.readouterr().err.lower()


def test_fetch_requires_a_scope(tmp_path, capsys):
    m = tmp_path / "manifest.json"
    m.write_text("{}")
    rc = main(["fetch", str(m), "--out", str(tmp_path / "s")])
    assert rc == EXIT_BAD_INPUT
    assert "scope" in capsys.readouterr().err.lower()


def test_fetch_rejects_both_scopes(tmp_path):
    m = tmp_path / "manifest.json"
    m.write_text("{}")
    rc = main(["fetch", str(m), "--all", "--split", "test", "--out", str(tmp_path / "s")])
    assert rc == EXIT_BAD_INPUT
