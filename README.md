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
| **Reproducible** | A `dataset.lock` pins the snapshot **and the split output**; `verify` re-derives both and certifies they reproduced byte-for-byte (or reports exactly what drifted). |
| **Cheap** | Metadata-only — a split is megabytes of JSON, not a terabyte of mmCIF. |
| **Honest about quality** | Every ligand is tiered (`functional` / `ambiguous` / `artifact`) with a reason; nothing is silently dropped. |
| **Fold-aware** | *Reduces and measures* structural leakage, not just sequence: for chains the chosen taxonomy classifies, same-family chains can't straddle train/test. Measured effect (controlled, ECOD-scored): +50% novel-fold test entries vs a sequence-only split, but the test set is still 98.6% ECOD-fold-*seen*. Unclassified chains are held out by sequence only — a disclosed residual-leakage *ceiling*, not elimination. |

### Two reproducibility guarantees

1. **Snapshot by release date, not query time.** Entries are selected by
   `release_date <= snapshot_date`, so the same `snapshot_date` yields the same
   *candidate set* (obsoleted entries are tracked, not silently dropped). The exact
   *split*, though, reproduces from the `dataset.lock` + `candidates.jsonl`, **not**
   from `snapshot_date` alone: RCSB recomputes its sequence clusters and CATH/ECOD/SCOP2
   annotations over time (no public history), so a *fresh* re-enumeration months later
   can differ — share the lock + candidates to reproduce a split byte-for-byte.
2. **Deterministic cluster → split assignment.** Under the `hash` strategy a
   cluster's split is decided by hashing a stable cluster key into the cumulative
   split fractions — independent of how many other clusters exist, so a cluster whose
   key is unchanged never moves as the PDB grows (and `verify` can certify a `hash`
   build). The one exception is a *merge*: a later bridging multi-chain entry (or a
   shared fold under `structural_clustering`) can unite two prior clusters, and the
   absorbed one's entries then follow the survivor. With a `splits.registry.json` this is
   both **prevented and reported**: the registry pins prior assignments matched on *any*
   key a component covers (so a held-out cluster stays held out across a merge; conflicts
   resolve `test > val > train`), and any pin it must override is **counted** in
   `splits.pinned_reassignments` (and, in entries, `splits.pinned_entries_reassigned`).
   One thing outranks a pin: under `maximal` the holdout ceiling. The dominant component
   absorbs held-out components as the snapshot grows and so can inherit their held-out pin;
   honoring that would put ~98% of the PDB in the holdout, so a component that does not fit
   the budget is capped to train whatever it was pinned to. Read
   `pinned_entries_reassigned`, not `pinned_reassignments`: one override on the giant moves
   most of the corpus while the component count reads `1`. Without a registry — the `hash`
   reassignment still happens and is **not** counted (`pinned_reassignments` needs a
   registry). An in-place rebuild instead reports it at the **entry level**: how many prior
   entries changed split, and how many were absorbed **into train** — the direction
   registry-free `hash` merges are biased toward (the survivor's bucket, and train owns 80%
   of it), i.e. held-out data eroding into train. Aggregate fractions can't detect this
   (they're conserved by construction). Pass `--registry` (or use `balanced`, which
   auto-adopts one in place when the config matches; `--fresh` opts out) to pin merges.

### Fold-level leakage control

Sequence clustering alone is not enough for inverse folding. A model learns
**structure → sequence**, so two chains below the 30% identity threshold that
nonetheless share a **fold** (TIM barrels, Rossmann folds, globins…) leak
structural information across the split — the model has effectively seen the test
backbone during training. `structural_clustering` **narrows** this: protein entities
sharing a structural **(super)family** *in the chosen authority* are union-merged into
the same component in addition to shared sequence clusters, so a family that authority
names cannot straddle train/test. It does not eliminate structural leakage — scored with
ECOD as an independent authority, a SCOP2 fold-aware test set is still 98.6%
ECOD-fold-*seen*. It reduces and, more usefully, **measures** the leak.

It uses RCSB's precomputed **CATH / ECOD / SCOP2** classifications — still metadata
only, no coordinates — selectable per build (`off | cath | ecod | scop2 | union`). It is
**purely additive**: it can only merge components, never split them, and a chain
with no classification simply contributes no structural edge. `if-split stats`
reports how many components the structural pass folded together, so the effect is
always measurable (`scripts/eval_structural_clustering.py` compares the methods).

Coverage is partial by nature: CATH ≈ 41%, ECOD ≈ 77%, SCOP2 ≈ 51% of protein
entities are classified (83.9% for their union; measured on the full 2026-07-22
snapshot — one entity per distinct sequence, *not* per chain instance), and
coverage falls sharply for recent depositions — 97% of entries
released through 2018 are classified, but only 60% of 2025 and 42% of 2026 ones.
The rest fall back to sequence-only. Chains with no classification are held out by sequence
only, so their fold-level hold-out is not guaranteed — the unclassified fraction
is a residual-leakage *ceiling* (an upper bound), which `if-split stats` reports
per split. It is a bound, not a measured number: IF-Split never measures fold
leakage from coordinates.

**Why it needs a balance-aware split.** On its own, fold-merging collapses the
dominant superfamilies (antibodies, TIM barrels) into mega-components that land
wholesale in one split, skewing the *entry* balance to ~95/3/2 (the *component*
split stays ~80/10/10). `split_strategy: "balanced"` fixes this: it **caps the
dominant folds to train and fills val/test to their entry targets from the tail of
smaller folds** — restoring ~80/10/10 by entries. Be precise about what lands there:
in the shipped fold-aware run only 12.2% of test and 15.1% of val entries carry a SCOP2
label at all (473 and 519 distinct families, ~992 held entirely out of train), so val/test
are a minority of genuinely held-out folds plus a majority of *unclassified* chains held
out by sequence alone. It stays leakage-safe (whole components) and
growth-stable via the registry (an in-place rebuild auto-adopts
`<out>/splits.registry.json` when its config matches; `--fresh` opts out — see
[Quickstart](#quickstart)), and reports a gap if a method's tail is too thin.
`balanced` also fixes the plain sequence-only skew (88/6/6 → 80/10/10) from the
antibody mega-cluster.

**Two honest caveats.** (1) *No fold authority is ground truth, and a "fold-disjoint"
split is authority-relative.* The PDB's fold-similarity graph **percolates**: under any
authority a single connected fold-component (built by multi-domain bridging) swallows most
of the PDB — measured on the 2026-07-22 snapshot, the largest component holds **79%** under
SCOP2, **87%** under CATH, **97% under ECOD** (their *union* percolates further still — see
`structural_clustering: "union"`). A fold-disjoint 80/10/10 exists only by capping that
giant component to train and holding out the small residual tail, and the tail size *is*
the authority's percolation: ECOD (fullest coverage) percolates most, leaving too thin a
tail to fill even one 10% holdout. `config/fold-aware.yaml` ships `scop2` because it
percolates least and *can* fill 80/10/10 — i.e. because it merges least, not because it is
the strongest criterion. **Mind the denominator:** SCOP2's ~62% entry coverage is
concentrated in *train* (73.7% of train entries classified); the held-out sets are the
opposite — only **12.2% of test** entries (2,649/21,683) and 15.1% of val carry a SCOP2
assignment, because classified entries are exactly the ones that merged into the capped
giant. So the guarantee covers ~12% of test and is silent on the other ~88%. Scored with
ECOD as an independent authority, a SCOP2 split is **98.6% ECOD-fold-*seen* in test**
(98.4% in val) — barely better than the sequence-only splits this critiques. Report a split
as *SCOP2*-fold-disjoint **with its per-split coverage**, never as fold-disjoint; and score
it with an authority you did *not* split on (`fold_benchmark_method`), since scoring with
the merge criterion is circular and will certify a leaky set as clean. (2) *`balanced`
shifts the val/test distribution by construction:* dominant folds are capped to train, so
val/test hold small, rare folds *and* — mostly — chains the authority never classified. That
is still a harder test than a sequence split, but
recovery on a `balanced` test set is **not comparable** to published LigandMPNN/ProteinMPNN
numbers or across strategies — report it as its own measurement.

### The fold-aware split

```bash
uv run if-split build --config config/fold-aware.yaml --out data/mc
```

`config/fold-aware.yaml` = default **+ `structural_clustering: scop2` + `split_strategy: balanced`**: a fold-aware ~80/10/10 split holding 992 distinct SCOP2 families entirely out of train (473 test + 519 val). It is the strongest fold control the tool can produce from metadata — *not* a fold-clean split: scored with ECOD it is still 98.6% fold-seen in test.

| recipe | strategy | structural | entry balance | val/test |
|---|---|---|--:|---|
| `default.yaml` | maximal | all | 98.6 / 0.7 / 0.7 | the whole leakage-safe tail; 888 certified-novel |
| `certified.yaml` | maximal | all | 98.6 / 0.7 / 0.7 | the published split (`config_hash f7d4203586df3dc7b10d2948e76d20d8`) |
| `fold-aware.yaml` | balanced | scop2 | 80 / 10 / 10 | 992 held-out SCOP2 families; still 98.6% ECOD-fold-seen |

A fixed 80/10/10 is available (`split_strategy: balanced`) but is not the default, because
under strong merging the leakage-safe tail is smaller than a 10% target and a proportional
assignment then fills test first and empties val. `maximal` reads `split_fractions` as a
**ceiling** and holds out whatever the tail allows.

### Novel-fold benchmark (fold-seen vs novel-fold)

Set `fold_benchmark_method: cath|ecod|scop2` (default `off`) to export a
**fold-seen vs novel-fold** partition of the test set — a turnkey way to
re-benchmark an *existing* checkpoint on folds it never saw in training. Fold
**labels** are decoupled from the fold **merging** used by `structural_clustering`:
they never feed union-find, so the split output and `check_no_leakage` are
byte-identical whether it's on or off, and labels attach even to a fold-*leaky*
split (e.g. the sequence-only split a public checkpoint was trained on).

A test entry is **novel-fold** iff it is fold-classified and *every* one of its
fold (super)families is absent from train. Unclassified test entries are neither
novel nor seen — the unclassified fraction is a residual-leakage *ceiling*, not a
measured number (IF-Split never measures fold leakage from coordinates).

It writes three sidecars into `--out`: `folds.json` (per-entry fold labels +
`novel_fold` flag, all splits), `fold_groups.json` (per-superfamily TEST groups
for per-family reweighting), and `novel_fold_test.json` (just the novel-fold test
ids).

```python
from ifsplit.dataset import load_dataset

ds = load_dataset("data/out/manifest.json")
ds.test.novel_fold_entries()   # test ids whose every fold is unseen in train
ds.test.is_novel_fold("1abc")  # bool
ds.test.folds_of("1abc")       # ["3.40.50.300", ...]
ds.fold_groups()               # {superfamily -> {novel, test_entries}} for reweighting
```

`if-split stats` prints the novel-fold count. The export is opt-in metadata only:
enabling it never downloads coordinates and never changes the split.

## Known limits

Things worth knowing before you build a split on this, stated here rather than left to be
discovered:

- **Identity thresholds are RCSB's, not arbitrary.** `identity_threshold` accepts only the
  precomputed levels 30/50/70/90/95/100%. There is no 25% or 40%. This follows from the
  no-mmseqs2 decision (no local binary, no version drift); the validator fails loudly rather
  than silently producing singletons. If you need to sweep identity thresholds, this is the
  wrong tool.
- **RCSB's 30% cluster file is not identity-complete.** Identical sequences can carry
  different cluster ids, so every chain with enough modeled residues also keys by a sequence
  hash. Below that gate, identical short peptides are deliberately *not* guaranteed to share
  a split.
- **Reproducible is not the same as stable under growth.** A pinned snapshot rebuilds
  byte-identically forever (`dataset.lock` + `candidates.jsonl` + `verify`). Regenerating on
  a *later* snapshot does not: under `maximal`, replaying a 2023-cutoff build on the
  2026-07-22 snapshot left only 43.7% of test entries in test (52.4% moved to val, 3.9% into
  train). Pin and cite a snapshot; treat a rebuild as a new benchmark.
- **The reproducibility anchor is the candidates file, not the date.** RCSB recomputes
  clusters and CATH/ECOD/SCOP2/Pfam/InterPro assignments over time with no public history,
  so `snapshot_date` alone does not pin a split. Share the lock plus `candidates.jsonl`.
- **A `maximal` holdout is the tail, not a miniature of the PDB.** On 2026-07-22 the
  held-out set skews short (median longest chain 158 residues against 295 in train; 27.1%
  have no chain over 50 residues) because the authorities classify short chains least often,
  so those are the components the merge does not absorb. Raise `min_modeled_residues` if that
  matters for your task, at the cost of the archived config hash.
- **A fold guarantee is only as wide as the configured authority's coverage.** A `scop2`
  build carries a SCOP2 label on ~12% of its test entries and leaves the rest constrained by
  sequence alone. Quote the denominator with the guarantee, always.
- **PDB only.** No AlphaFold DB or other predicted structures; the union-find has no fold
  assignment for them. If your training corpus includes predicted models, this covers part
  of your leakage problem, not all of it.
- **It stops at split lists and labels.** No coordinate parsing, no tensorization. `fetch`
  is an optional downstream convenience; the featurizer is yours.

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
# Scaffold a config interactively — pick a recipe (default | fold-aware) and tweak the
# key knobs. It writes the config (comments intact) and prints the build command; it
# NEVER runs a build. Works from an installed wheel too. --non-interactive for scripts.
uv run if-split init --out my-split.yaml

# Build the full split from today's PDB (metadata only).
uv run if-split build --config config/default.yaml --out data/out

# Dev: cap to the first N candidates (by sorted entry id — still reproducible).
uv run if-split build --limit 50 --out /tmp/ifs

# Summarize a build: split sizes + entry-fraction skew vs target, per-class test
# counts, curation tiers, growth-stability status, fold coverage, and (if exported)
# the novel-fold count.
uv run if-split stats data/out/manifest.json

# Reproduce-check: re-derive from a lock, report candidate drift vs the live PDB,
# and (when candidates reproduce) certify the split output matches its hash.
uv run if-split verify data/out/dataset.lock

# Offline verify: integrity-check a distributed candidates.jsonl + lock, no network.
uv run if-split verify data/out/dataset.lock --candidates data/out/candidates.jsonl

# Re-derive the split from a CACHED candidates.jsonl — no RCSB. Ablate curation,
# clustering, split strategy, or TIGHTEN a filter (e.g. resolution) on a fixed
# snapshot in ~seconds instead of re-enumerating the whole PDB.
uv run if-split resplit --candidates data/out/candidates.jsonl \
    --config config/fold-aware.yaml --out data/mc

# Growth-stable regeneration (balanced strategy): rebuilding IN PLACE auto-adopts the
# prior <out>/splits.registry.json when its config matches, preserving prior
# component→split assignments as the PDB grows. --fresh starts a new lineage;
# --registry <path> pins from a different prior build. (`hash` is growth-stable on
# its own and ignores the registry, so `verify` still certifies it.)
uv run if-split build --config config/fold-aware.yaml --out data/out          # re-run: auto-pins
uv run if-split build --config config/fold-aware.yaml --out data/out --fresh   # opt out

# OPTIONAL: download the actual structures for a built split (see below).
uv run if-split fetch data/out/manifest.json --split test --out data/structures

# Emit a portable, shareable split spec (see "Sharing a split spec" below).
uv run if-split spec data/out/manifest.json --name my-split --out my-split.ifsplit.yaml
```

### Outputs (`--out` directory)

| File | Purpose |
|---|---|
| `train.json` / `val.json` / `test.json` | **The split** — plain JSON arrays of PDB entry ids (one per line, grepable and trivially loadable). |
| `test/<class>_test.json` | Test ids carrying each functional ligand class (`metal` / `small_molecule` / `nucleic_acid`), for per-class evaluation. |
| `manifest.json` | Human-facing run record: config, drop log, per-split + per-class (and ambiguous) counts, cluster/leakage stats, per-split fold coverage, growth-stability flag, and a `files` index. |
| `dataset.lock` | Reproduction anchor: embedded config + candidates SHA-256 + entry list + a `split` hash of the entry→split partition (so `verify` certifies the split output, not just the inputs). |
| `candidates.jsonl` | The snapshot definition — one canonical JSON record per entry. Hashed into the lock. |
| `clusters.json` | `entry_id → component key`, for cluster-balanced sampling. |
| `ligands.classes.json` | `entry_id → functional ligand class labels`. |
| `ligands.tiers.json` | Per-component ligand curation audit trail (tier + reason). |
| `targets.jsonl` | Conditioning-target corpus: one row per (structure, functional ligand) for ligand-conditioned training. |
| `splits.registry.json` | `component key → split`, for growth-stable regeneration. |
| `folds.json`, `fold_groups.json`, `novel_fold_test.json` | **Novel-fold benchmark** export (only with `fold_benchmark_method` set): per-entry fold labels + novel-fold flag, per-superfamily test groups (for reweighting), and the novel-fold test subset. |

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
| 3 — filter | `parse.py` | Drop no-protein / no-usable-sequence (empty **or** all-`X` poly-UNK) / too-short (opt-in `min_modeled_residues`) / oversized entries (assembly-1 residue count `> max_total_residues`) / over-resolution (re-derived here so it is auditable; per-method caps via `resolution_max_A_by_method`), plus optional wwPDB validation-report quality caps (clashscore, R-free, Ramachandran/rotamer/RSRZ, cryo-EM map-fit floor) — all from metadata. Every drop is logged with its reason. |
| 4 — ligands | `ligands.py` | Tier each non-protein component `functional`/`ambiguous`/`artifact`; derive class labels (metal / small-molecule / nucleic-acid). `nucleic_acid` = a protein↔DNA/RNA *complex* (verified assembly interface), **not** a bound mononucleotide. **Annotate, never drop.** |
| 5 — cluster | `cluster.py` | Group protein entities by RCSB precomputed cluster id at `identity_threshold`; canonical key = smallest member id. Optionally union same-fold entities (CATH/ECOD/SCOP2) for structural-leakage control. |
| 6 — split | `split.py` | Assign components → train/val/test (`maximal` (default) for the data-sized holdout, `hash` for per-component hashing, `balanced` for entry-balanced fold-tail val/test); assert no cluster spans two splits; audit residual secondary-chain overlap. |
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
| `min_em_backbone_inclusion` | backbone atom-in-density (a **floor** — higher is better) | cryo-EM |

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
| `ambiguous` | Present but uncorroborated → reported, **not** labelled | `metal_unbound`, `ligand_unbound`, `metal_site_nonnative`, `glycan`, `purification_metal_uncorroborated` |
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

**Glycans aren't ligand pockets.** A carbohydrate (RCSB CCD type `*saccharide*` —
NAG/BMA/MAN/…, and sugar-detergents like LMT) is overwhelmingly decorative
glycosylation or a purification detergent, not a site an inverse-folding model
conditions on. So a carbohydrate is tiered `glycan` (reported, not a small-molecule
target) unless it has a *measured binding affinity* — RCSB's `is_subject_of_investigation`
(SOI) flag is too prevalent on sugars to rescue here (see below). A genuine
lectin/glycosidase ligand is recoverable as an opt-in target (`include_ambiguous=True`).
Real cofactors (ATP, NAD, HEM, FAD) and structural lipids (cardiolipin, phosphatidyl-*)
are `non-polymer`, not saccharides, so they're untouched.

**Curation beats a comp-id rule — where being wrong is expensive.** Since v0.6.2 the SOI
flag outranks the *additive blacklist* (order: affinity → glycan → SOI → blacklist). Every
comp id is an additive *somewhere*, so a blanket blacklist convicts exactly the
depositions where the additive is the point — lauric acid and dodecyl sulfate bound in the
β-lactoglobulin calyx, since a lipocalin's function **is** binding fatty acids and
detergents. The flag deliberately does **not** outrank the glycan gate, because the two
gates differ on both axes that matter (3,000-entry stride sample, 2026-08 snapshot):

| gate | SOI prevalence | a wrong call is |
|---|---|---|
| additive blacklist | 4.1% (87/2,114) | `artifact` — **never** reaches `targets.jsonl` |
| glycan | 85.7% (209/244; 106 are NAG) | `ambiguous` — emitted, one `include_ambiguous=True` away |

So SOI is used where it carries signal and where being wrong is unrecoverable, and kept
out where it carries almost none and being wrong costs a flag rather than the data.

### Test-set representation

The split is deterministic (a per-component hash, or the `balanced` fold-tail
fill), so the test set's ligand mix is reported but not forced by default:
`manifest.json` carries per-split, per-class `functional` counts plus `ambiguous`
counts, so under-representation is visible.
An opt-in `test_min_per_class` floor (e.g. `{metal: 500}`) tops up
under-represented classes by recruiting WHOLE `functional`-tier sequence
components into test in deterministic hash order — never individual entries (no
leakage) and skipping registry-pinned components (growth stays stable). A floor
beyond the available supply is met as far as possible and the shortfall is
reported in the manifest (and by `if-split stats`), not forced.

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

### Training strategies for inverse folding

An inverse-folding model consumes two different things, with opposite scale/quality
tradeoffs. IF-Split emits **both from the same leakage-safe split**, so you pick a
strategy without re-deriving the split:

| Corpus | What | Use it for |
|---|---|---|
| **Backbones** — every kept structure | `ds.train.backbones` | ProteinMPNN-style, ligand-agnostic. The scale lever — a structure with only junk ions or no ligand is still a good backbone. |
| **Conditioning targets** — the `functional`-tier ligands | `ds.train.conditioning_targets()` | LigandMPNN-style. One row per `(structure, ligand)`; junk is never a target. The quality lever. |

```python
ds = load_dataset("data/out/manifest.json")

# 1. Backbone-only training (max data): every structure.
backbones = ds.train.backbones

# 2. Ligand-conditioned training: condition on the real ligands only.
targets = ds.train.conditioning_targets()                 # metal / small_molecule / nucleic_acid
metal_targets = ds.train.conditioning_targets(classes=["metal"])

# 3. Condition on ALL of a structure's ligands at once (group by entry), or one at a time:
for entry_id, ligs in ds.train.targets_by_entry().items():
    ctx = [(t.ligand_class, t.comp_id) for t in ligs]     # e.g. [("small_molecule","HEM")]

# 4. Only structures that actually carry a conditioning target:
conditioned = ds.train.conditioned_entry_ids()            # subset of backbones

# 5. Opt in to ambiguous targets — non-native metal pockets (Ni/Co substituting the
#    native metal) and glycans (glycosylation / lectin ligands) — off by default:
any_site = ds.train.conditioning_targets(include_ambiguous=True)
```

The full corpus is also written to `targets.jsonl` (one row per target: entry, split,
cluster, class, comp_id, tier, reason) and mirrored into `index.parquet`'s
`conditioning_targets` column after `fetch`.

#### From a target to its pocket (your featurizer owns this)

IF-Split stops at the **labels** — it never parses coordinates. The last mile is
yours: join a target's `comp_id` to a fetched structure with your own parser and
pull the ligand atoms + pocket. This is deliberate — featurization is
model-specific (ligand atoms? SMILES? a pocket mask? all-atom context?), and every
inverse-folding model already has its own pipeline. `comp_id` is the join key:

```python
from pathlib import Path

import gemmi  # or biotite / Biopython — same pattern
from ifsplit.dataset import load_dataset
from ifsplit.download import rel_path_for

ds = load_dataset("data/out/manifest.json")
for entry_id, targets in ds.test.targets_by_entry().items():
    path = Path("data/structures") / rel_path_for(entry_id, "test", assembly=True)
    model = gemmi.read_structure(str(path))[0]
    for t in targets:
        # A structure may hold several copies of one ligand (plus adventitious ions).
        # Surface ALL copies; pick the instance to condition on per your model.
        copies = [r for ch in model for r in ch if r.name == t.comp_id]
        # ... extract residues within config.ligand_context_radius_A of each copy ...
```

When an entry has several functional ligands (e.g. a cofactor *and* an inhibitor,
or a catalytic metal alongside an adventitious one), that shows up as **multiple
target rows** — condition on all (group by `entry_id`) or one per example, your
call. A runnable, copy-and-adapt version with pocket extraction is
[`scripts/consume_split.py`](scripts/consume_split.py).

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
  created_with: if-split 0.5.0
  expected_config_hash: 3b63318286fd2ac4994f34d10936be05
snapshot_date: '2026-07-14'
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
| `dataset.lock` | *"Reproduce the exact bytes"* — pins entry set + candidates SHA + split-output hash | ~MB |

## Configuration

Everything that affects the output lives in one YAML file
([`config/default.yaml`](config/default.yaml)); its canonical hash is embedded
in every manifest, so two builds with the same hash used identical settings. It
doubles as a shareable **split spec** — see [Sharing a split spec](#sharing-a-split-spec).

| Key | Default | Meaning |
|---|---|---|
| `snapshot_date` | `2026-05-30` | `release_date <= this` — the reproducibility anchor. |
| `experimental_methods` | X-ray, EM | Allowed `exptl.method` values. |
| `resolution_max_A` | `3.5` | Resolution cutoff (re-derived in Stage 3, so it is auditable from `candidates.jsonl`). |
| `resolution_max_A_by_method` | `{}` | Optional per-method resolution overrides, e.g. `{ELECTRON MICROSCOPY: 3.0}` — cryo-EM 3.5 Å ≠ X-ray 3.5 Å. Empty = one cap for all. |
| `max_total_residues` | `5999` | Max residues **kept** (drop if `> this`) — LigandMPNN kept `< 6000`, i.e. `<= 5999`. |
| `min_modeled_residues` | `0` | Opt-in floor on modeled (non-`X`) residues in a protein chain. `0` = off; only the always-on empty/all-`X` (poly-UNK) drop applies. `~20` also drops tiny/mostly-unknown chains. |
| `single_chain_only` | `false` | Opt-in: keep only single-protein-entity structures (no complex, no nucleic-acid partner) — a metadata proxy for single-chain design targets (ProteinMPNN's single-chain CATH setup). |
| `excluded_het` | waters + common ions | Extra components forced to `artifact`. |
| `use_biological_assembly` | `true` | Count residues from assembly 1, not the deposited asymmetric unit. |
| `purification_metals` | `[NI, CO]` | Metals treated as IMAC tags; `[]` disables the heuristic. |
| `histag_min_run` | `6` | His-run length (anywhere) that marks a purification tag. |
| `histag_terminal_min_run` | `3` | Shorter His-run at a chain terminus that also counts as a tag (partial/unmodeled 6×His). |
| `exclude_purification_artifacts` | `true` | Demote His-tag metals to `artifact`; lone uncorroborated Ni/Co → `ambiguous`. |
| `identity_threshold` | `0.30` | Clustering cutoff (RCSB levels: 30/50/70/90/95/100). |
| `clustering_backend` | `precomputed` | Reuse RCSB's published 30% clusters (the only backend; locked via the snapshot). |
| `structural_clustering` | `all` | Fold-level leakage control: `off` \| `cath` \| `ecod` \| `scop2` \| `pfam` \| `interpro` \| `union` (the curated three) \| `all` (all five). Union-merges same-(super)family protein chains so a family *that authority classifies* can't straddle train/test (it reduces and measures fold leakage; it does not eliminate it) — see [Fold-level leakage control](#fold-level-leakage-control). Additive: only merges, never splits. |
| `split_fractions` | 0.80 / 0.10 / 0.10 | train / val / test. |
| `split_strategy` | `maximal` | `maximal` (size the holdout from the data: cap any component too large for the budget to train, hold out the whole remaining tail, `split_fractions` acting as a **ceiling** not a target), `hash` (balance components; input-independent and registry-free, so `verify` can certify it), or `balanced` (balance *entries*: cap dominant folds to train, fill val/test from the fold tail). |
| `fold_benchmark_method` | `off` | Novel-fold benchmark export: `off` \| `cath` \| `ecod` \| `scop2`. When set, emits `folds.json` / `fold_groups.json` / `novel_fold_test.json`; fold *labels* are decoupled from fold *merging*, so it never changes the split or `check_no_leakage`. Omitted from the config hash when `off`. |
| `test_min_per_class` | `{}` | Optional per-class floor of test entries carrying each functional ligand class; recruits WHOLE sequence components into test in deterministic order (never individual entries → no leakage), skipping registry-pinned components. Empty = off. |
| `split_salt` | `snapsplit-v1` | Bump to intentionally reshuffle the split. |
| `max_clashscore`, `max_rfree`, `max_ramachandran_outlier_pct`, `max_rotamer_outlier_pct`, `max_rsrz_outlier_pct`, `min_em_backbone_inclusion`, `require_validation_report` | off | Optional validation-report quality caps — see [Structure quality](#structure-quality-validation-report). |
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

## Contributing

Bug reports, feature requests, and pull requests are welcome — see
[CONTRIBUTING.md](CONTRIBUTING.md) for the dev setup (`uv sync`, `pytest`, `ruff`)
and the two load-bearing invariants (determinism; no cross-split leakage). All
participation is under our [Code of Conduct](CODE_OF_CONDUCT.md).

## Citation

If you use IF-Split, please cite it — see [CITATION.cff](CITATION.cff).

## Changelog

Release history is in [CHANGELOG.md](CHANGELOG.md). The latest tagged release is **0.6.0**; **0.6.1** is prepared but not yet tagged. Previously **0.5.0**
(the novel-fold benchmark: a fold-seen vs novel-fold export for re-benchmarking existing
checkpoints, a growth-stability fix for the balanced split, and an entry-skew stats view).

## License

MIT — see [LICENSE](LICENSE).
