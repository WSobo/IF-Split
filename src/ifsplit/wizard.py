"""``if-split init`` — interactively scaffold a build ``config.yaml`` from a recipe.

The wizard picks one of the two base recipes (``default`` or ``fold-aware``),
prompts for the handful of highest-signal knobs (each defaulting to the recipe's
value — press Enter to keep it), and writes a config the user can then ``build``.
It NEVER runs a build and downloads nothing; its only side effect is writing the
config file. Every other knob keeps the recipe's inline comment for hand-editing.

The base recipes are embedded (``_recipes.py``) because the wheel omits ``config/``,
so ``init`` works from an installed package with no repo checkout. Overrides are
applied by surgically rewriting only the touched lines, so the emitted config keeps
all of the recipe's explanatory comments; the result is validated as a
:class:`~ifsplit.config.Config` before it is written, so ``init`` never emits an
invalid file.
"""

from __future__ import annotations

import re
import sys
from collections.abc import Callable
from datetime import date
from pathlib import Path

import yaml
from pydantic import ValidationError

from ._recipes import RECIPES
from .config import Config, SplitFractions, is_tautological_benchmark

DEFAULT_INIT_OUT = "if-split.yaml"
_FOLD_LEVELS = ("off", "cath", "ecod", "scop2")
# A top-level YAML key at column 0 (never a comment or an indented/flow-mapping key).
_KEY_RE = re.compile(r"^(?P<key>[A-Za-z0-9_]+):")


# --------------------------------------------------------------------------- #
# Pure text transform: apply value overrides to a recipe, comments intact.     #
# --------------------------------------------------------------------------- #
def apply_overrides(text: str, overrides: dict[str, str]) -> str:
    """Return ``text`` with each ``key -> rendered-YAML-value`` override applied.

    A key present as an active top-level line has only its value rewritten (its
    inline comment is preserved); a key with no active line (e.g. a recipe leaves
    it commented out to fall back to the default) is appended at the end. Values
    are pre-rendered YAML fragments (e.g. ``'"scop2"'`` or ``'3.0'``).
    """
    if not overrides:
        return text
    remaining = dict(overrides)
    out: list[str] = []
    for line in text.splitlines(keepends=True):
        m = _KEY_RE.match(line)
        if m and m.group("key") in remaining:
            key = m.group("key")
            out.append(_rewrite_line(line, key, remaining.pop(key)))
        else:
            out.append(line)
    if remaining:
        if out and not out[-1].endswith("\n"):
            out[-1] += "\n"
        for key, value in remaining.items():
            out.append(f"{key}: {value}\n")
    return "".join(out)


def _rewrite_line(line: str, key: str, newval: str) -> str:
    """Rewrite one ``key: value  # comment`` line's value, keeping the comment."""
    suffix = "\n" if line.endswith("\n") else ""
    core = line[: -len(suffix)] if suffix else line
    _, _, rest = core.partition(":")  # rest = ' oldvalue   # comment'
    if "#" in rest:
        valpart, comment = rest.split("#", 1)
        gap = valpart[len(valpart.rstrip()) :]  # whitespace between value and comment
        return f"{key}: {newval}{gap}#{comment}{suffix}"
    return f"{key}: {newval}{suffix}"


def _validate(text: str) -> Config:
    """Parse + validate the config text; raises on anything :class:`Config` rejects."""
    raw = yaml.safe_load(text)
    if not isinstance(raw, dict):
        raise ValueError("config is not a YAML mapping")
    return Config.model_validate(raw)


# --------------------------------------------------------------------------- #
# Prompt helpers (parsers return the rendered YAML value to substitute).        #
# --------------------------------------------------------------------------- #
def _ask(input_fn: Callable[[str], str], label: str, default_display: object, parse) -> str | None:
    """Prompt until ``parse`` accepts; blank input or EOF returns ``None`` (keep default)."""
    prompt = f"{label} [{default_display}]: "
    while True:
        try:
            raw = input_fn(prompt).strip()
        except EOFError:
            return None
        if not raw:
            return None
        try:
            return parse(raw)
        except (ValueError, ValidationError) as exc:
            print(f"  ! {exc}")


def _parse_date(raw: str) -> str:
    date.fromisoformat(raw)  # raises ValueError on a bad date
    return f'"{raw}"'


def _parse_pos_float(raw: str) -> str:
    if float(raw) <= 0:
        raise ValueError(f"must be > 0, got {raw}")
    return raw


def _parse_fractions(raw: str) -> str:
    parts = [p for p in re.split(r"[,\s]+", raw) if p]
    if len(parts) != 3:
        raise ValueError("give three numbers: train val test (e.g. 0.8 0.1 0.1)")
    try:
        t, v, te = (float(p) for p in parts)
        SplitFractions(train=t, val=v, test=te)  # reuse the real range + sum-to-1 validator
    except (ValueError, ValidationError):
        raise ValueError("each fraction must be in (0, 1) and the three must sum to 1.0") from None
    return f"{{train: {t:g}, val: {v:g}, test: {te:g}}}"


def _parse_salt(raw: str) -> str:
    if '"' in raw or "#" in raw:
        raise ValueError('salt may not contain a quote (") or hash (#)')
    return f'"{raw}"'


def _enum_parser(allowed: tuple[str, ...]):
    def parse(raw: str) -> str:
        low = raw.strip().lower()
        if low not in allowed:
            raise ValueError(f"choose one of: {', '.join(allowed)}")
        return f'"{low}"'

    return parse


def _ask_recipe(input_fn: Callable[[str], str]) -> str:
    print("Choose a base recipe:")
    print("  default     - hash split, sequence-only clustering (standard reproducible split)")
    print("  fold-aware  - scop2 fold-merge + balanced (fold-honest val/test)")

    def parse(raw: str) -> str:
        low = raw.strip().lower()
        if low not in RECIPES:
            raise ValueError(f"choose one of: {', '.join(RECIPES)}")
        return low

    chosen = _ask(input_fn, "recipe", "default", parse)
    return chosen if chosen is not None else "default"


def _collect_overrides(recipe_text: str, input_fn: Callable[[str], str]) -> dict[str, str]:
    """Prompt for the highest-signal knobs; return only the ones the user changed."""
    d = yaml.safe_load(recipe_text)
    sf = d["split_fractions"]
    frac_display = f"{sf['train']:g} {sf['val']:g} {sf['test']:g}"
    fold_bench_default = d.get("fold_benchmark_method", "off")

    questions = [
        (
            "snapshot_date",
            "snapshot_date (release_date <= this; YYYY-MM-DD)",
            d["snapshot_date"],
            _parse_date,
        ),
        (
            "resolution_max_A",
            "resolution_max_A (Angstrom cap)",
            d["resolution_max_A"],
            _parse_pos_float,
        ),
        ("split_fractions", "split_fractions (train val test)", frac_display, _parse_fractions),
        (
            "structural_clustering",
            "structural_clustering (fold-merge: off|cath|ecod|scop2)",
            d["structural_clustering"],
            _enum_parser(_FOLD_LEVELS),
        ),
        (
            "fold_benchmark_method",
            "fold_benchmark_method (novel-fold export: off|cath|ecod|scop2)",
            fold_bench_default,
            _enum_parser(_FOLD_LEVELS),
        ),
        ("split_salt", "split_salt (bump to reshuffle)", d["split_salt"], _parse_salt),
    ]
    overrides: dict[str, str] = {}
    for key, label, default_display, parse in questions:
        value = _ask(input_fn, label, default_display, parse)
        if value is not None:
            overrides[key] = value
    return overrides


# --------------------------------------------------------------------------- #
# CLI entry point.                                                             #
# --------------------------------------------------------------------------- #
def run_init(
    args, *, input_fn: Callable[[str], str] | None = None, interactive: bool | None = None
):
    """Scaffold a config from a recipe. Returns a process exit code (0 or bad-input)."""
    from .cli import EXIT_BAD_INPUT, EXIT_OK

    if input_fn is None:
        input_fn = input
    if interactive is None:
        interactive = sys.stdin.isatty() and not getattr(args, "non_interactive", False)

    out_path = Path(getattr(args, "out", None) or DEFAULT_INIT_OUT)
    if out_path.exists() and not getattr(args, "force", False):
        print(
            f"error: {out_path} already exists; pass --force to overwrite or --out for a new path",
            file=sys.stderr,
        )
        return EXIT_BAD_INPUT

    recipe = getattr(args, "recipe", None)
    if interactive and recipe is None:
        print("if-split init - scaffold a build config (press Enter to keep each [default]).")
        recipe = _ask_recipe(input_fn)
    if recipe is None:
        recipe = "default"
    if recipe not in RECIPES:
        print(f"error: unknown recipe {recipe!r} (choose: {', '.join(RECIPES)})", file=sys.stderr)
        return EXIT_BAD_INPUT

    text = RECIPES[recipe]
    overrides = _collect_overrides(text, input_fn) if interactive else {}
    text = apply_overrides(text, overrides)

    # Belt-and-suspenders: never write a config that would fail to load.
    try:
        cfg = _validate(text)
    except (ValueError, ValidationError) as exc:
        print(f"error: the scaffolded config is invalid ({exc})", file=sys.stderr)
        return EXIT_BAD_INPUT

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")

    if is_tautological_benchmark(cfg):
        print(
            f"  ! note: fold_benchmark_method == structural_clustering "
            f"({cfg.fold_benchmark_method}) — the novel-fold benchmark will read a tautological "
            f"~100%. Use a different authority for one, or structural_clustering: off, to make it "
            f"meaningful. (Editing the file to change it is fine; the split is unaffected.)"
        )

    print()
    print(f"Wrote {out_path}  (recipe: {recipe}, config_hash {cfg.config_hash()})")
    print("Build the split with:")
    print(f"  if-split build --config {out_path} --out data/out")
    print(f"Preview the match count first:  if-split build --config {out_path} --count")
    print(
        "The file keeps an inline comment for every other knob - edit and re-run. "
        "`if-split spec` stamps shareable provenance."
    )
    return EXIT_OK
