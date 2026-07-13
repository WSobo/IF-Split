# IF-Split

[![CI](https://github.com/WSobo/IF-Split/actions/workflows/ci.yml/badge.svg)](https://github.com/WSobo/IF-Split/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/if-split.svg)](https://pypi.org/project/if-split/)
[![Python](https://img.shields.io/pypi/pyversions/if-split.svg)](https://pypi.org/project/if-split/)

**A reproducible, date-pinned, ligand-aware train/val/test splitter for the PDB.**

IF-Split borrows the *split logic* of LigandMPNN (Dauparas et al., *Nature
Methods* 2025) — cluster proteins at 30% sequence identity, partition so no
cluster spans two splits, categorize the test set by ligand class — but instead
of inheriting a frozen 2022 snapshot it generates the split **on demand from
today's PDB**, and emits a lock file so a collaborator can reproduce the exact
dataset later. See [PLAN.md](PLAN.md) for the full design spec.

It is built entirely on RCSB **metadata** (the Search + Data APIs): **no
structure coordinates are downloaded** to build a split — only tiny per-entry
records and sequences. Coordinates are an optional, downstream concern.

---

## Why it's different

| | |
|---|---|
| **Fresh** | Builds from the current PDB, not a years-old frozen copy. |
| **Reproducible** | A `dataset.lock` pins the snapshot; `verify` re-derives it and reports any drift. |
| **Cheap** | Metadata-only — a split is megabytes of JSON, not a terabyte of mmCIF. |
| **Honest about quality** | Every ligand is tiered (`functional` / `ambiguous` / `artifact`) with a reason; nothing is silently dropped. |

### Two reproducibility guarantees

1. **Snapshot by release date, not query time.** Entries are selected by
   `release_date <= snapshot_date`. Re-running with the same `snapshot_date`
   yields the same candidate set no matter *when* you run it (obsoleted entries
   are tracked, not silently dropped).
2. **Deterministic cluster → split assignment.** A cluster's split is decided by
   hashing a stable cluster key into the cumulative split fractions — independent
   of how many other clusters exist. Existing clusters never move when the PDB
   grows, which is what prevents train/test leakage on regeneration. A
   `splits.registry.json` pins prior assignments to make this exact even across
   re-clustering.

---

## Install

Requires Python ≥ 3.11. `build` needs only network access to RCSB — no external
binaries.

```bash
pip install if-split          # from PyPI
```

Or for development, with [`uv`](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/WSobo/IF-Split && cd IF-Split
uv sync          # creates .venv from uv.lock, installs deps + dev tools (ruff, pytest)
```

`uv.lock` is committed, so dev environments are reproducible. (The optional
coordinate/featurization path via `gemmi` is Linux-native, so run under
Linux/WSL if you use `fetch`.)

## Quickstart

```bash
# Build the full split from today's PDB (metadata only).
uv run if-split build --config config/default.yaml --out data/out

# Dev: cap to the first N candidates (by sorted entry id — still reproducible).
uv run if-split build --limit 50 --out /tmp/ifs

# Summarize a build: split sizes, per-class test counts, curation tiers.
uv run if-split stats data/out/manifest.json

# Reproduce-check: re-derive from a lock and report drift vs the live PDB.
uv run if-split verify data/out/dataset.lock

# Growth-stable regeneration: pin prior cluster→split assignments.
uv run if-split build --registry data/out/splits.registry.json --out data/out2

# OPTIONAL: download the actual structures for a built split (see below).
uv run if-split fetch data/out/manifest.json --split test --out data/structures

# Emit a portable, shareable split spec (see "Sharing a split spec" below).
uv run if-split spec data/out/manifest.json --name my-split --out my-split.ifsplit.yaml
```

### Outputs (`--out` directory)

| File | Purpose |
|---|---|
| `candidates.jsonl` | The snapshot definition — one canonical JSON record per entry. Hashed into the lock. |
| `dataset.lock` | Reproduction anchor: embedded config + candidates SHA-256 + entry list. |
| `manifest.json` | Human-facing run record: per-split entry lists, ligand classes + tiers, per-class (and ambiguous) counts, drop log, cluster/leakage stats, entry→cluster map. |
| `splits.registry.json` | `cluster key → split`, for growth-stable regeneration. |

## Downloading structures (`fetch`)

`build` produces a tiny, coordinate-free split. When you actually want the mmCIF
files — to featurize or train — `fetch` hydrates a *built manifest* into a
clean, ML-ready tree. It is **opt-in and downstream**: nothing about a split
requires coordinates.

```bash
# Scope is explicit by design (no accidental terabyte): choose splits or --all.
uv run if-split fetch data/out/manifest.json --split test                 # just test
uv run if-split fetch data/out/manifest.json --split train --split val    # repeatable
uv run if-split fetch data/out/manifest.json --all --yes --workers 16     # everything
uv run if-split fetch data/out/manifest.json --all --asymmetric-unit      # AU not assembly 1
```

`fetch` prints an estimated download size first and refuses pulls over ~1000
structures without `--yes`. It is **resumable** (existing, valid files are
skipped) and parallel (`--workers`).

### Layout — browsable *and* scalable

Files are split-partitioned (so you can `ls` a split) and sharded by the PDB
"divided" scheme — the middle two characters of the entry id — so no single
directory holds an unwieldy number of files:

```
data/structures/
  structures/
    train/  hh/4hhb-assembly1.cif.gz   01/101m-assembly1.cif.gz   02/102l-… 102m-…
    val/    …
    test/   0a/10ad-assembly1.cif.gz
  index.jsonl            # one row per structure (zero-dep, greppable)
  index.parquet          # same, columnar (written if pyarrow is installed)
  manifest.json          # copy of the source split manifest
  DATASET_CARD.md        # provenance + how-to-load
```

The **index** is the ML entry point — one row per structure with `entry_id`,
`split`, `path`, **`sha256`** (integrity + dedupe), `cluster` (for
cluster-balanced batches), and `ligand_classes` / `ligand_tiers`:

```python
import pandas as pd
df = pd.read_parquet("data/structures/index.parquet")   # or read_json(..., lines=True)
train = df[df.split == "train"]
metal_train = train[train.ligand_classes.str.contains("metal")]
# de-redundified epoch: one structure per sequence cluster
epoch = train.sort_values("entry_id").groupby("cluster").head(1)
```

The columnar `index.parquet` needs `pyarrow`: `uv sync --extra mlops` (the
zero-dependency `index.jsonl` is always written regardless).

## How it works

A `build` runs eight stages; none touch coordinates.

| Stage | Module | What it does |
|---|---|---|
| 1 — enumerate | `enumerate.py`, `rcsb.py` | RCSB Search → entry IDs; Data API (GraphQL, batched) → sequences, ligands, residue counts, cluster membership → `candidates.jsonl`. |
| 3 — filter | `parse.py` | Drop no-protein / no-sequence / oversized entries (assembly-1 residue count vs `max_total_residues`), plus optional wwPDB validation-report quality caps (clashscore, R-free, Ramachandran/rotamer/RSRZ) — all from metadata. Every drop is logged with its reason. |
| 4 — ligands | `ligands.py` | Tier each non-protein component `functional`/`ambiguous`/`artifact`; derive class labels (metal / small-molecule / nucleic-acid). `nucleic_acid` = a protein↔DNA/RNA *complex* (verified assembly interface), **not** a bound mononucleotide. **Annotate, never drop.** |
| 5 — cluster | `cluster.py` | Group protein entities by RCSB precomputed cluster id at `identity_threshold`; canonical key = smallest member id. |
| 6 — split | `split.py` | Deterministic hash → train/val/test; assert no cluster spans two splits; audit residual secondary-chain overlap. |
| 7 — manifest | `manifest.py` | Emit lock + manifest + registry (all deterministic, no wall-clock fields). |
| 8 — loader | `dataset.py` | Read a manifest into train/val/test views with cluster-balanced sampling. |
| 2 — fetch *(opt-in)* | `download.py`, `hydrate.py` | Download mmCIF for a built manifest into a sharded, indexed, ML-ready tree. |

> Stage 2 (mmCIF coordinate download) is **optional and downstream** — only
> needed to extract ligand context or feed a model, never to build a split. See
> [Downloading structures](#downloading-structures-fetch) for the `fetch` command.

### Structure quality (validation report)

For the highest-quality backbones, `build` can filter on the **wwPDB validation
report** — fetched as metadata, so the no-download invariant still holds. The
metrics come straight from the deposited report:

| Cap | Metric | Applies to |
|---|---|---|
| `max_clashscore` | all-atom clashscore | X-ray + cryo-EM |
| `max_ramachandran_outlier_pct` | % backbone Ramachandran outliers | X-ray + cryo-EM |
| `max_rotamer_outlier_pct` | % sidechain rotamer outliers | X-ray + cryo-EM |
| `max_rfree` | R-free (DCC) | X-ray |
| `max_rsrz_outlier_pct` | % real-space-R Z-score outliers | X-ray |

Two rules keep it honest: a cap fires **only when the metric is present**, so a
cryo-EM entry is never dropped for a missing R-free; and every cap is **off by
default**, so the snapshot is unchanged until you opt in. `require_validation_report`
drops entries with no report at all. Each drop is logged with its reason and
value (e.g. `clashscore_too_high`) and is summarised by `if-split stats`.

> Strict starting point: `max_clashscore: 40`, `max_rfree: 0.30`,
> `max_ramachandran_outlier_pct: 1.0`. Some classic low-quality depositions drop
> out — e.g. the 1984 entry `4HHB` has a clashscore of 142.

### Ligand quality: annotate, don't destroy

IF-Split is a *tool*, not one frozen dataset, so it won't make an irreversible
quality call for you. Every non-protein component is tiered, with a
machine-readable reason, from RCSB metadata signals:

| Tier | Meaning | Example reasons |
|---|---|---|
| `functional` | Real ligand/site → gets a class label | `metal_bound`/`ligand_bound` (contacts protein), `*_affinity` (measured), `*_investigated` (RCSB SOI), `metal_annotated` (protein annotated to bind this metal) |
| `ambiguous` | Present but uncorroborated → reported, **not** labelled | `metal_unbound`, `ligand_unbound`, `metal_site_nonnative`, `purification_metal_uncorroborated` |
| `artifact` | Buffer / counterion / purification tag → excluded from labels | `additive`, `counterion`, `histag_metal` |

**Holo gating (metadata-only).** Presence isn't enough. A small molecule or metal
is `functional` only if RCSB reports it *contacting* the protein (`bound_components`)
or it has a measured binding affinity; an unbound one is `ambiguous`. A DNA/RNA
chain is `functional` `nucleic_acid` only when the biological assembly has a verified
protein↔nucleic-acid interface (`num_prot_na_interface_entities > 0`) — a
co-deposited but non-contacting oligo is reported `ambiguous`, never silently
labelled. (Interfaces are RCSB-computed metadata, available for X-ray *and* cryo-EM,
so no coordinates are downloaded.)

> The `nucleic_acid` class is the protein–nucleic-acid **complex** category (DNA/RNA
> polymer chains), matching LigandMPNN's "nucleotide" split. Bound *mononucleotide*
> ligands (ATP, GTP, NAD, SAM, …) are not this class — they fall under
> `small_molecule`.

The His-tag/Ni curation catches a known blemish in the LigandMPNN metal set:
structures whose only "metal site" is a poly-His tag chelating Ni/Co from
affinity purification. A poly-His run anywhere — or a short run at a chain
terminus (`histag_terminal_min_run`, catching 6×His tags left partial by
unmodeled or trimmed residues) — flags the entry's Ni/Co as an `artifact`.

But an audit (reproducible via [`scripts/audit_nico_histag.py`](scripts/audit_nico_histag.py))
showed a subtler issue: **~82% of lone Ni/Co entries carry no detectable His-tag
in the deposited sequence** — IMAC tags are frequently absent from the SEQRES
record, not just unmodeled, so a sequence scan can't recover them. So even with no
detectable tag, a *lone* Ni/Co (the entry's only metal) with no corroboration is
demoted from `functional` to `ambiguous` — reported, not labelled.

To avoid over-firing on genuine bare-Ni/Co enzymes (urease, cobalt methionine
aminopeptidase, nitrile hydratase, …), a lone Ni/Co is **rescued to `functional`
(`metal_annotated`)** when the protein's RCSB GO/InterPro/Pfam annotation says it
binds that metal. A protein that binds a *different* native metal (Ni/Co as an
isomorphous substitute — e.g. Co in a Mg enzyme) is reported `metal_site_nonnative`
so a consumer can choose to keep it; one with no metal annotation at all stays
`purification_metal_uncorroborated`. All from RCSB's own metadata (no extra
UniProt call). Real metals (Zn, Mg, Fe, …), and Ni/Co with affinity/SOI or beside
a genuine metal, are untouched. Rerun [`scripts/eval_metal_tiering.py`](scripts/eval_metal_tiering.py)
to measure the tier distribution over the whole lone-Ni/Co set.

Crucially, **the structure always stays in its split** — a protein with a junk
ion is still a good backbone; we just don't label the junk. A consumer wanting
"pristine metal sites only" vs "maximum scale, I'll filter myself" changes a
threshold, not the build. The same per-component tier is what a downstream
featurizer reads to decide what counts as real ligand context.

> **Per-instance is a featurizer concern.** These tiers are per *component* — they
> establish whether a structure contains a real Ni/Co site, not *which* of several
> same-element ions is it. A deposition can hold both a catalytic Ni and a surface
> crystallization Ni under one `NI` id, and no metadata separates them (adventitious
> Ni binds surface His/Asp with the same geometry as a catalytic site). Deciding
> which individual ion to featurize is left to the coordinate-level featurizer.

### Test-set representation

The split is a pure deterministic hash, so the test set's ligand mix is reported
but not forced by default: `manifest.json` carries per-split, per-class
`functional` counts plus `ambiguous` counts, so under-representation is visible.
An opt-in `--enforce-minimums` top-up (recruit `functional`-only ligand clusters
into test in deterministic order) is scoped for a future release.

### Using a split (loader)

```python
from ifsplit.dataset import load_dataset

ds = load_dataset("data/out/manifest.json")
print(len(ds.train), len(ds.val), len(ds.test))

# Ligand-class views.
metal_test = ds.test.with_class("metal")

# Cluster-balanced sampling: one representative per sequence cluster per epoch,
# so over-represented folds (lysozyme, common kinases) don't dominate.
for epoch in range(3):
    batch_ids = ds.train.sample_by_cluster(seed=epoch)
```

## Sharing a split spec

The config **is** the shareable recipe. Everything that affects the split lives in
one small YAML file with a content hash, so you can hand someone that file and they
reproduce your methodology exactly — like `params.yaml` in DVC. `if-split spec`
emits a portable, self-identifying version from any build or config:

```bash
# Extract a stand-alone spec from a finished build (config is embedded in the manifest):
uv run if-split spec data/out/manifest.json --name "my-split" --author "you" \
    --out my-split.ifsplit.yaml

# Anyone reproduces your split from just that file:
uv run if-split build --config my-split.ifsplit.yaml --out their/out
```

The emitted file carries a `spec:` header that announces what it is and pins the
expected hash:

```yaml
spec:
  ifsplit_spec: ifsplit/config@1          # schema id — the file says what it is
  name: my-split
  author: you
  created_with: if-split 0.1.0
  expected_config_hash: 3b63318286fd2ac4994f34d10936be05
snapshot_date: '2026-05-31'
resolution_max_A: 3.5
# ... all output-affecting settings ...
```

On load, if `expected_config_hash` no longer matches the settings (someone edited
them after stamping), IF-Split warns. The `spec:` metadata is **excluded from the
hash**, so name/author/description never change the split identity — two specs that
differ only in their labels produce byte-identical outputs.

| Artifact | Question it answers | Size |
|---|---|--:|
| `*.ifsplit.yaml` (or `config.yaml`) | *"How did you make this split?"* — the recipe | ~KB |
| `manifest.json` | *"What's in it?"* — counts, provenance, file index | ~KB |
| `dataset.lock` | *"Reproduce the exact bytes"* — pins entry set + candidates SHA | ~MB |

## Configuration

Everything that affects the output lives in one YAML file
([`config/default.yaml`](config/default.yaml)); its canonical hash is embedded
in every manifest, so two builds with the same hash used identical settings. It
doubles as a shareable **split spec** — see [Sharing a split spec](#sharing-a-split-spec).

| Key | Default | Meaning |
|---|---|---|
| `snapshot_date` | `2026-05-30` | `release_date <= this` — the reproducibility anchor. |
| `experimental_methods` | X-ray, EM | Allowed `exptl.method` values. |
| `resolution_max_A` | `3.5` | Resolution cutoff. |
| `max_total_residues` | `5999` | Size cap (LigandMPNN used `< 6000`). |
| `excluded_het` | waters + common ions | Extra components forced to `artifact`. |
| `use_biological_assembly` | `true` | Count residues from assembly 1, not the deposited asymmetric unit. |
| `purification_metals` | `[NI, CO]` | Metals treated as IMAC tags; `[]` disables the heuristic. |
| `histag_min_run` | `6` | His-run length (anywhere) that marks a purification tag. |
| `histag_terminal_min_run` | `3` | Shorter His-run at a chain terminus that also counts as a tag (partial/unmodeled 6×His). |
| `exclude_purification_artifacts` | `true` | Demote His-tag metals to `artifact`; lone uncorroborated Ni/Co → `ambiguous`. |
| `identity_threshold` | `0.30` | Clustering cutoff (RCSB levels: 30/50/70/90/95/100). |
| `clustering_backend` | `precomputed` | `precomputed` (RCSB clusters) or `mmseqs2` (run your own). |
| `split_fractions` | 0.80 / 0.10 / 0.10 | train / val / test. |
| `split_salt` | `snapsplit-v1` | Bump to intentionally reshuffle the split. |
| `max_clashscore`, `max_rfree`, `max_ramachandran_outlier_pct`, `max_rotamer_outlier_pct`, `max_rsrz_outlier_pct`, `require_validation_report` | off | Optional validation-report quality caps — see [Structure quality](#structure-quality-validation-report). |
| `ligand_context_radius_A`, `max_ligand_atoms` | `8.0`, `25` | Featurization only (not part of the split). |

## Develop

```bash
uv run pytest              # tests (offline; 1 opt-in network test, see below)
uv run ruff check .        # lint
uv run ruff format .       # format

# Run the opt-in live RCSB round-trip test.
IFSPLIT_NETWORK_TESTS=1 uv run pytest tests/test_integration.py
```

## Layout

```
config/default.yaml      # single source of truth for a run (hashed into the manifest)
src/ifsplit/             # config.py + one module per pipeline stage
  enumerate.py rcsb.py   #   Stage 1: RCSB Search + Data API
  parse.py               #   Stage 3: metadata filters
  ligands.py             #   Stage 4: ligand tiering + classification
  cluster.py split.py    #   Stages 5-6: clustering + deterministic split
  manifest.py            #   Stage 7: lock + manifest + registry, verify/stats
  dataset.py             #   Stage 8: loader + cluster-balanced sampling
  download.py            #   Stage 2: optional mmCIF fetch (featurization only)
data/cache/              # downloaded mmCIF, if ever used (gitignored)
data/out/                # generated manifests + lock files
tests/
```

## License

MIT — see [LICENSE](LICENSE).
