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

**All phases complete.** The full pipeline runs end-to-end against live RCSB and
is deterministic — **no structure coordinates are downloaded** (see
[PLAN.md](PLAN.md) §1.5):

Stage 1 enumerate (Search + Data API → `candidates.jsonl`) → Stage 3 filter
(residue cap, no-protein, drop log) → Stage 4 ligand **confidence tiering**
(`functional`/`ambiguous`/`artifact` from bound-components + binding-affinity
signals, incl. His-tag/Ni purification-artifact curation — structures are
annotated, never dropped) → Stage 5 cluster (RCSB precomputed per-entity
membership) → Stage 6 deterministic hash split (no-cluster-leakage invariant
asserted, growth-stable via a split registry) → Stage 7 `manifest.json` +
`dataset.lock` + `splits.registry.json` → Stage 8 loader with
**cluster-balanced sampling**.

Verified: two `build` runs produce byte-identical manifests; `verify`
round-trips the lock; live tiering correctly calls e.g. 101M `{HEM: functional,
SO4: artifact}`. Optional remaining work: on-demand coordinate fetch, the
`mmseqs2` clustering backend, and the opt-in `--enforce-minimums` test
stratification top-up (all outside the core split path).

### Quality model: annotate, don't destroy

IF-Split is a *tool*, not a single frozen dataset, so it won't make an
irreversible quality call for you. Every ligand gets a **confidence tier** + a
reason (e.g. `metal_bound`, `ligand_affinity`, `histag_metal`, `additive`,
`counterion`, `metal_unbound`) recorded in the manifest. Class labels derive from
the `functional` tier by default; `ambiguous` is reported but unlabelled;
`artifact` is excluded from labels — **but the structure always stays in its
split**. A consumer wanting "pristine metal sites only" vs "maximum scale, I'll
filter myself" changes a threshold, not the build. The same per-component tier is
what a downstream featurizer reads to decide what counts as real ligand context.

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
# Build the full split from RCSB (metadata only).
uv run if-split build --config config/default.yaml --out data/out
# Dev: cap to the first N candidates by sorted entry id (reproducible).
uv run if-split build --limit 50 --out /tmp/ifs

# Reproduce-check: re-enumerate from a lock and report drift vs the live PDB.
uv run if-split verify data/out/dataset.lock
# Summarize a built manifest (split sizes, per-class test counts, curation).
uv run if-split stats data/out/manifest.json

# Growth-stable regeneration: pin prior cluster->split assignments.
uv run if-split build --registry data/out/splits.registry.json --out data/out2
```

`build` writes: `candidates.jsonl` (snapshot definition — one canonical JSON
record per entry), `dataset.lock` (embedded config + candidates SHA-256 + entry
list), `manifest.json` (per-split entry lists, ligand-class tags, per-class test
counts, drop log, cluster/leakage stats), and `splits.registry.json`
(canonical-key → split, for growth-stable regeneration). `verify` re-runs Stage 1
and reports added/removed entries + hash match — warning rather than failing so
reproductions are honest about upstream changes.

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
