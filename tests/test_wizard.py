"""Tests for ``if-split init`` — the config-scaffolding wizard (offline, no network).

Covers: the embedded recipes stay byte-identical to ``config/*.yaml`` (the drift
guard the wheel depends on), the comment-preserving override transform, the
interactive prompt flow, non-interactive verbatim emission, and the guardrails
(overwrite protection, invalid-input re-prompt, always-valid output).
"""

from __future__ import annotations

from datetime import date
from pathlib import Path

import yaml

from ifsplit._recipes import RECIPES
from ifsplit.cli import EXIT_BAD_INPUT, EXIT_OK, build_parser, main
from ifsplit.config import Config, load_config
from ifsplit.wizard import apply_overrides, run_init

ROOT = Path(__file__).resolve().parents[1]


def _scripted(answers):
    """A fake ``input`` that returns each answer in turn, then raises EOF (keep default)."""
    it = iter(answers)
    calls: list[str] = []

    def fake(prompt: str) -> str:
        calls.append(prompt)
        try:
            return next(it)
        except StopIteration:
            raise EOFError from None

    fake.calls = calls
    return fake


# --------------------------------------------------------------------------- #
# The embedded recipes must match the source-of-truth config files.            #
# --------------------------------------------------------------------------- #
def test_embedded_recipes_match_config_files():
    for name, fname in [("default", "default.yaml"), ("fold-aware", "fold-aware.yaml")]:
        disk = (ROOT / "config" / fname).read_text(encoding="utf-8")
        assert RECIPES[name] == disk, f"{name} recipe drifted from config/{fname}"


def test_embedded_recipes_are_valid_configs():
    for text in RECIPES.values():
        Config.model_validate(yaml.safe_load(text))


# --------------------------------------------------------------------------- #
# apply_overrides: rewrite the value, keep the comment; append when commented.  #
# --------------------------------------------------------------------------- #
def test_apply_overrides_rewrites_value_and_keeps_comment():
    out = apply_overrides(RECIPES["default"], {"structural_clustering": '"scop2"'})
    line = next(ln for ln in out.splitlines() if ln.startswith("structural_clustering:"))
    assert line.startswith('structural_clustering: "scop2"')
    assert "# fold-leakage control" in line  # inline comment preserved
    assert Config.model_validate(yaml.safe_load(out)).structural_clustering == "scop2"


def test_apply_overrides_appends_key_absent_as_active_line():
    # fold-aware.yaml leaves fold_benchmark_method commented out (defaults to "off"),
    # so an override has no active line to rewrite and must be appended.
    out = apply_overrides(RECIPES["fold-aware"], {"fold_benchmark_method": '"scop2"'})
    cfg = Config.model_validate(yaml.safe_load(out))
    assert cfg.fold_benchmark_method == "scop2"
    assert cfg.structural_clustering == "scop2"  # untouched
    active = [ln for ln in out.splitlines() if ln.startswith("fold_benchmark_method:")]
    assert active == ['fold_benchmark_method: "scop2"']


def test_apply_overrides_noop_is_identity():
    assert apply_overrides(RECIPES["default"], {}) == RECIPES["default"]


# --------------------------------------------------------------------------- #
# Interactive flow (scripted input; interactive forced on).                     #
# --------------------------------------------------------------------------- #
def test_interactive_full_override(tmp_path, capsys):
    out = tmp_path / "cfg.yaml"
    args = build_parser().parse_args(["init", "--out", str(out)])
    # benchmark on an INDEPENDENT authority (ecod) vs the scop2 merge — a meaningful,
    # non-tautological pairing (the scop2==scop2 warning path has its own test).
    fake = _scripted(
        ["fold-aware", "2026-07-24", "2.8", "0.7 0.15 0.15", "scop2", "ecod", "my-salt"]
    )
    rc = run_init(args, input_fn=fake, interactive=True)
    assert rc == EXIT_OK
    assert len(fake.calls) == 7  # recipe + six knobs

    cfg = load_config(out)
    assert cfg.snapshot_date == date(2026, 7, 24)
    assert cfg.resolution_max_A == 2.8
    sf = cfg.split_fractions
    assert (sf.train, sf.val, sf.test) == (0.7, 0.15, 0.15)
    assert cfg.structural_clustering == "scop2"
    assert cfg.fold_benchmark_method == "ecod"
    assert cfg.split_salt == "my-salt"
    # Base recipe's teaching comments survive.
    assert 'IF-Split "fold-aware" config' in out.read_text()
    assert "if-split build --config" in capsys.readouterr().out


def test_interactive_blank_keeps_recipe_defaults(tmp_path):
    out = tmp_path / "cfg.yaml"
    args = build_parser().parse_args(["init", "--out", str(out)])
    # Choose default, change only the split fractions, keep everything else.
    fake = _scripted(["default", "", "", "0.9 0.05 0.05"])
    assert run_init(args, input_fn=fake, interactive=True) == EXIT_OK

    text = out.read_text()
    cfg = load_config(out)
    sf = cfg.split_fractions
    assert (sf.train, sf.val, sf.test) == (0.9, 0.05, 0.05)
    assert cfg.structural_clustering == "all"  # untouched default
    assert cfg.split_salt == "snapsplit-v1"  # untouched default
    # The unchanged snapshot_date line still carries its explanatory comment.
    snap = next(ln for ln in text.splitlines() if ln.startswith("snapshot_date:"))
    assert "reproducibility anchor" in snap


def test_wizard_warns_on_tautological_benchmark(tmp_path, capsys):
    out = tmp_path / "cfg.yaml"
    args = build_parser().parse_args(["init", "--out", str(out)])
    # fold-aware recipe keeps structural_clustering=scop2; set fold_benchmark_method=scop2 too
    # (blanks keep snapshot/resolution/fractions/structural) — that pairing is tautological.
    fake = _scripted(["fold-aware", "", "", "", "", "scop2"])
    assert run_init(args, input_fn=fake, interactive=True) == EXIT_OK
    assert "tautological" in capsys.readouterr().out


def test_interactive_reprompts_on_invalid_then_accepts(tmp_path, capsys):
    out = tmp_path / "cfg.yaml"
    args = build_parser().parse_args(["init", "--out", str(out)])
    fake = _scripted(["default", "not-a-date", "2026-01-01"])  # bad date, then valid
    assert run_init(args, input_fn=fake, interactive=True) == EXIT_OK
    assert load_config(out).snapshot_date == date(2026, 1, 1)
    assert "  ! " in capsys.readouterr().out  # the re-prompt error line


# --------------------------------------------------------------------------- #
# Non-interactive + CLI wiring + guardrails.                                    #
# --------------------------------------------------------------------------- #
def test_non_interactive_writes_recipe_verbatim(tmp_path, capsys):
    for recipe in ("default", "fold-aware"):
        out = tmp_path / f"{recipe}.yaml"
        rc = main(["init", "--out", str(out), "--recipe", recipe, "--non-interactive"])
        assert rc == EXIT_OK
        assert out.read_text() == RECIPES[recipe]  # byte-for-byte the base recipe
        load_config(out)  # loads + validates
    assert "if-split build --config" in capsys.readouterr().out


def test_non_interactive_default_recipe_and_out(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rc = main(["init", "--non-interactive"])
    assert rc == EXIT_OK
    written = tmp_path / "if-split.yaml"
    assert written.exists() and written.read_text() == RECIPES["default"]


def test_refuses_to_overwrite_without_force(tmp_path, capsys):
    out = tmp_path / "cfg.yaml"
    out.write_text("do not clobber")
    rc = main(["init", "--out", str(out), "--recipe", "default", "--non-interactive"])
    assert rc == EXIT_BAD_INPUT
    assert "already exists" in capsys.readouterr().err
    assert out.read_text() == "do not clobber"  # untouched


def test_force_overwrites(tmp_path):
    out = tmp_path / "cfg.yaml"
    out.write_text("stale")
    rc = main(["init", "--out", str(out), "--recipe", "fold-aware", "--non-interactive", "--force"])
    assert rc == EXIT_OK
    assert out.read_text() == RECIPES["fold-aware"]
    assert load_config(out).split_strategy == "balanced"
