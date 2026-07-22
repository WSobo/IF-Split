# IF-Split — agent guide

Reproducible, date-pinned, ligand-aware train/val/test splitter for the PDB. It
borrows LigandMPNN's *split logic* but generates it on demand from today's PDB.
Read [PLAN.md](PLAN.md) for the design rationale and [README.md](README.md) for
usage. This file is the orientation for working in the repo.

## The one load-bearing idea

**The split is computed from metadata + sequences only — `build` never downloads
structure coordinates.** Everything needed (resolution, method, release date,
residue counts, per-entity sequences, ligand chem-comp + bound-component signals,
RCSB sequence-cluster membership, and CATH/ECOD/SCOP2 structural classifications)
comes from the RCSB Search + Data APIs. Coordinates
(mmCIF) are large and only needed downstream, so `fetch` (Stage 2) is optional.
Keep it that way: do not add coordinate access to the build path.

## Environment (WSL over a Windows UNC mount)

- Code lives in WSL at `~/projects/IF-Split`, opened over `\\wsl.localhost\ubuntu\...`.
- **Run everything through WSL**, not the Bash tool's git-bash (that's Windows
  Python and lacks the deps): `wsl -d ubuntu bash -lc '...'`.
- The deps live in a **uv** venv; `uv` is at `$HOME/.local/bin`, so prefix:
  `export PATH="$HOME/.local/bin:$PATH"`.
- Use **single quotes** for `bash -lc '...'`. Avoid backticks and `python3 -c "..."`
  inside it — nested quotes/backticks get shell-evaluated and corrupt the command
  (this has mangled commit messages and grep filters). For commit messages, write
  the text to a file and use `git commit -F file`.
- **Do not put a command that may exit nonzero in a parallel tool batch** — one
  failure cancels every sibling call in that batch. Run risky/dependent commands
  one at a time.

## Commands

```bash
wsl -d ubuntu bash -lc 'cd ~/projects/IF-Split && export PATH="$HOME/.local/bin:$PATH" && uv run pytest -q'
uv run ruff check .      # lint (must pass)
uv run ruff format .     # format
uv run if-split build --limit 50 --out /tmp/ifs   # dev build (small, live RCSB)
uv run if-split build --config config/fold-aware.yaml --out /tmp/mc  # fold-honest split (scop2 + balanced)
uv run if-split resplit --candidates data/out/candidates.jsonl --config X.yaml --out /tmp/rs  # re-derive Stages 3-7 offline (no RCSB)
```

- `uv sync` sets up the env; `uv sync --extra mlops` adds pyarrow for `fetch`'s
  parquet index. Dev tools (ruff, pytest) are a PEP 735 dependency-group.
- The offline test suite needs no network. One opt-in live test runs only with
  `IFSPLIT_NETWORK_TESTS=1`.

## Architecture (src/ifsplit/, one module per stage)

`enumerate.py`+`rcsb.py` (Stage 1, Search+Data API → candidates.jsonl) →
`parse.py` (3, metadata filters) → `ligands.py` (4, confidence tiering) →
`cluster.py` (5, union-find components: sequence + optional fold-level structural
clustering) → `split.py` (6, split assignment: `hash` | `balanced`) →
`manifest.py` (7, lock + manifest + registry, verify/stats) → `dataset.py` (8,
loader). `download.py`+`hydrate.py` are the optional Stage 2 `fetch`. `cli.py`'s
`resplit` re-runs Stages 3-7 from a cached `candidates.jsonl` (no Stage 1) via the
shared `_run_pipeline`; `verify --candidates` does the same for offline checking.

Invariants that must not regress:
- **Determinism:** same config → byte-identical `manifest.json` (no wall-clock
  fields). `test_manifest_is_deterministic` guards this.
- **No cross-split leakage:** sequence clusters joined by a shared multi-chain
  entry (and, with `structural_clustering` on, by a shared fold superfamily) are
  union-find–merged into one component; a component maps to exactly one split, so
  overlap is impossible by construction. This holds for both split strategies
  (`hash` and `balanced`, which only chooses *which* split a whole component lands
  in). `check_no_leakage` is a real invariant (not a tautology) — keep it that way.
- **Growth stability:** `hash` maps a component to `hash(salt + canonical_key)`
  into cumulative fractions, keyed on the global-min member id (not RCSB's volatile
  integer id) — input-independent and registry-free, so a larger snapshot only
  *adds* components and never moves existing ones. `balanced` differs: its val/test
  fill boundaries scale with the snapshot's total entries, so prior components move
  unless pinned. `splits.registry.json` pins them; an in-place `balanced` rebuild
  auto-adopts `<out>/splits.registry.json` when the prior `dataset.lock`
  `config_hash` matches (`--fresh` opts out), and the manifest records
  `splits.growth_stable`.
- **Annotate, never destroy:** ligand quality is a per-component *tier*
  (functional / ambiguous / artifact) in the manifest; structures are never
  dropped for ligand quality. Class labels derive from the functional tier.
- **PDB-ID compatibility:** store entry/entity ids verbatim from `rcsb_id` (legacy
  `4HHB` and extended `pdb_xxxxxxxx`); never slice/length-validate/case-fold them.

## Conventions

- Python ≥ 3.11, uv + ruff (line length 100). Keep ruff clean and tests green
  before committing.
- Don't commit on `main` without the user asking; end commit messages with the
  `Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>` trailer.
- GitHub remote: `github.com/WSobo/IF-Split` (public).
