# IF-Split

A reproducible, date-pinned, ligand-aware train/val/test splitter for the PDB.

IF-Split reproduces the *methodology* of the LigandMPNN data split (Dauparas et
al., *Nature Methods* 2025) — 30% sequence-identity clustering with
ligand-categorized test sets — but generates it on demand from a current PDB
snapshot, with a manifest + lock file so a collaborator can reproduce the exact
dataset later. See [PLAN.md](PLAN.md) for the full spec.

## Two reproducibility guarantees

1. **Snapshot by release date, not query time.** Entries are selected by
   `release_date <= snapshot_date`, so re-running with the same `snapshot_date`
   yields the same candidate set regardless of when you run it (modulo obsoleted
   entries, which are tracked explicitly).
2. **Deterministic cluster → split assignment.** A cluster's split is decided by
   `hash(canonical_cluster_id + salt) mod N`, so existing clusters never move
   when the dataset grows — the property that prevents train/test leakage on
   regeneration.

## Status

**Phases 1–2 complete.** Config layer + Stage 1 (enumerate via RCSB Search +
Data API → `candidates.jsonl`) + the snapshot lock and `verify` loop work
end-to-end against live RCSB — **no structure coordinates are downloaded**
(see [PLAN.md](PLAN.md) §1.5). Stages 3–8 (filter → ligands → cluster → split →
manifest → loader) are stubbed; see [PLAN.md](PLAN.md) §8 for the build order.

## Install

Requires Python ≥ 3.11 and [`uv`](https://docs.astral.sh/uv/). `build` itself
needs only network access to RCSB (no external binaries). The optional `mmseqs2`
clustering backend and the optional coordinate/featurization path (`gemmi`) are
Linux-native, so build/run under Linux/WSL.

```bash
uv sync          # create .venv from uv.lock and install deps + dev tools (ruff, pytest)
```

`uv sync` installs the `dev` dependency group (ruff, pytest) by default. The
lockfile `uv.lock` is committed for reproducible environments.

## Usage

```bash
# Enumerate the snapshot from RCSB (metadata only) -> candidates.jsonl + dataset.lock.
uv run if-split build --config config/default.yaml --out data/out
# Dev: cap to the first N candidates by sorted entry id (reproducible).
uv run if-split build --limit 50 --out /tmp/ifs

# Reproduce-check: re-enumerate from a lock and report drift vs the live PDB.
uv run if-split verify data/out/dataset.lock
```

`build` writes two artifacts: `candidates.jsonl` (the snapshot definition — one
canonical JSON record per entry) and `dataset.lock` (embedded config + the
candidates' SHA-256 + entry list). `verify` re-runs Stage 1 and reports added /
removed entries and whether the hash still matches — warning rather than failing
so reproductions are honest about what changed upstream. `stats <manifest.json>`
lands with the manifest in a later phase.

## Develop

```bash
uv run pytest              # tests
uv run ruff check .        # lint
uv run ruff format .       # format
```

## Layout

```
config/default.yaml   # single source of truth for a run (hashed into the manifest)
src/ifsplit/          # config.py + one module per pipeline stage (1–8)
data/cache/           # downloaded mmCIF (gitignored)
data/out/             # generated manifests + lock files
tests/
```
