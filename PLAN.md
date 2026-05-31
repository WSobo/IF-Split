# IF-Split — Build Plan

A reproducible, date-pinned, ligand-aware train/val/test splitter for the PDB.
Produces LigandMPNN-style splits (30% sequence-identity clustering,
ligand-categorized test sets) but generated on demand from a current PDB
snapshot, with a manifest/lock file that lets a collaborator reproduce the exact
dataset later. This document is the spec. Build it in the phase order at the
bottom.

## 1. Goal and the one hard constraint

Replicate the methodology of the LigandMPNN data split (Dauparas et al., *Nature
Methods* 2025) without inheriting its frozen 2022 snapshot. The output is:

1. A set of cleaned, parsed protein(+ligand) structures drawn from the PDB.
2. A train/val/test partition with no sequence-cluster leakage across splits.
3. Test sets categorized by ligand class: small-molecule, metal, nucleotide.
4. A manifest + lock file that makes the dataset reproducible byte-for-byte.

The non-negotiable design constraint: **fresh + reproducible at the same time.**
Two mechanisms deliver this, and they must both be in from the start:

- **Snapshot by release-date cutoff, not query time.** Entries are selected by
  `release_date <= snapshot_date`. Re-running with the same `snapshot_date`
  yields the same candidate set regardless of the calendar date the script runs
  (modulo obsoleted entries, which we track explicitly — see §6). Do not define
  the snapshot as "whatever the API returns today."
- **Deterministic cluster→split assignment by hash.** A cluster's split is
  decided by `hash(cluster_representative_id + salt) mod N`. New clusters added
  in a later, larger snapshot get assigned by the same function; existing
  clusters never move. This gives stable splits under dataset growth and is the
  property that prevents train/test contamination when you regenerate.

  > **Implementation note (Stage 6):** the "existing clusters never move"
  > guarantee holds only if `cluster_representative_id` is *input-independent*.
  > mmseqs2 picks a representative from whatever is in the input set, so as the
  > snapshot grows a cluster's representative can change and its hash with it.
  > Hash a canonical member (e.g. the lexicographically smallest member id over
  > the cluster's *full* membership in the locked cluster file, not just the
  > date-surviving members) and persist a cluster→split registry so growth only
  > ever *adds* clusters.

## 1.5 Architecture: metadata-first, no bulk structure downloads

The split is computed entirely from **metadata + sequences**, which are tiny and
come from RCSB APIs. 3D coordinates (mmCIF) are large and are needed only for
two things, both strictly *downstream of the split*: extracting ligand context
atoms, and feeding a model. So **`build` never downloads structures** — coordinate
download is an optional, on-demand featurization step (Stage 2, demoted).

What each concern actually consumes:

| Concern                                  | Needs                          | Source (no coordinates) |
|------------------------------------------|--------------------------------|-------------------------|
| entry selection (method/res/date)        | metadata                       | RCSB Search + Data API  |
| residue-count / chain filters            | metadata                       | Data API                |
| ligand classification (SM/metal/nucleo)  | chem-comp ids + formulas, polymer types | Data API       |
| sequence clustering                       | clusters (or sequences)        | precomputed cluster file (or Data API seqs) |
| cluster → split assignment                | cluster ids                    | derived                 |
| **ligand context atoms** (featurization) | **coordinates**                | mmCIF (on demand, opt.) |
| **model input** (featurization)          | **coordinates**                | mmCIF (on demand, opt.) |

This makes the snapshot rest on **two pinned anchors**, both fetched fresh at
build time and both lockable:

1. `snapshot_date` — selects entries by `release_date <= date`.
2. the **sequence cluster file** — RCSB recomputes `clusters-by-entity-30.txt`
   weekly and keeps no history, so we store it (≈17 MB) + its SHA-256 in the
   lock. Reproduction uses the *locked* file, never a re-download.

Same `snapshot_date` + same locked cluster file → identical split, forever.

**Clustering backend** is configurable: `precomputed` (default — reuse RCSB's
published polymer-entity clusters at the configured identity; no external binary,
instant) or `mmseqs2` (run our own over the snapshot's sequences for full
control). Both feed the same Stage 6 hash assignment.

> **Design stance:** LigandMPNN is the reference for the *split logic*
> (cluster-by-identity, no cross-split leakage, ligand-categorized test sets),
> **not** a byte-for-byte fidelity target. IF-Split builds on *today's RCSB*. So
> the default `precomputed` backend reuses RCSB's current clustering (computed
> with **DIAMOND** — `--cluster`, member-coverage 0.8, BLOSUM62, DNA/RNA and
> <10-aa peptides filtered, recomputed weekly), which is perfectly fine for the
> logic we're after and is made reproducible by locking the cluster file. The
> optional `mmseqs2` backend exists for anyone who wants to run their own
> clustering on the snapshot. Record the backend + parameters in the manifest so
> the choice is always explicit.

## 2. Repository layout

```
IF-Split/
  README.md
  PLAN.md                      # this file
  pyproject.toml               # uv / pip; pin deps
  config/
    default.yaml               # the single source of truth for a run
  src/ifsplit/
    __init__.py
    config.py                  # load + validate + hash the config
    enumerate.py               # Stage 1: RCSB Search + Data API -> candidates.jsonl
    download.py                # Stage 2: OPTIONAL on-demand mmCIF fetch (featurization)
    parse.py                   # Stage 3: metadata filters + drop log
    ligands.py                 # Stage 4: classify non-protein entities from metadata
    cluster.py                 # Stage 5: precomputed RCSB clusters (default) | mmseqs2
    split.py                   # Stage 6: deterministic hash assignment + stratification
    manifest.py                # Stage 7: emit manifest + lock file
    dataset.py                 # Stage 8: loader / torch Dataset consuming a manifest
    cli.py                     # `if-split build`, `if-split verify`, `if-split stats`
  data/
    cache/                     # downloaded mmCIF (gitignored)
    out/                       # generated manifests + lock files
  tests/
```

> Note: the importable package is `ifsplit` (Python identifiers can't contain
> `-`); the distribution and CLI are `if-split`.

### PDB identifier compatibility (legacy 4-char + extended)

wwPDB is migrating from 4-character PDB IDs to the **extended** form `pdb_` +
8 chars (legacy `4HHB` → `pdb_00004hhb`, case-insensitive). 4-char IDs are
expected to be exhausted around 2027–2028; entries issued after the switch get
only extended IDs. CCD chemical-component codes are likewise extending beyond
3 chars. Rules IF-Split follows so both forms work:

- **Store identifiers verbatim** from the Data API `rcsb_id` — never slice to 4
  chars, length-validate, or upper-case entry/entity ids. (Verified live: the
  API returns legacy `4HHB`/`4HHB_1` today; the schema also accepts
  `pdb_00009xyz`/`pdb_00009xyz_1`.)
- **Chemical-component (CCD) codes** stay upper-cased (CCD codes are
  case-insensitive uppercase), but are not length-restricted.
- **Caveat (verified 2026-05):** the Data API does *not yet* resolve an extended
  id passed as a query argument (`pdb_00004hhb` → empty result); it accepts
  legacy ids and returns legacy ids. So we query with whatever the Search API
  hands us and store what comes back — no client-side id reformatting either
  direction. When RCSB flips inputs/outputs to extended ids, the verbatim policy
  means no code change is needed.

## 3. Configuration (config/default.yaml)

Every parameter that affects the output lives here. The config is hashed and the
hash is embedded in the manifest, so two manifests with the same config-hash are
guaranteed to have used identical settings.

```yaml
snapshot_date: "2026-05-30"        # release_date <= this. The reproducibility anchor.
experimental_methods: ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"]
resolution_max_A: 3.5
max_total_residues: 5999           # LigandMPNN used "< 6000"
excluded_het: ["HOH", "NA", "CL", "K", "BR"]   # waters + common crystallization ions
use_biological_assembly: true      # biounits, as in LigandMPNN (assembly 1)
# purification-artifact curation (Stage 4): His-tag + Ni/Co only -> not a metal site
purification_metals: ["NI", "CO"]
histag_min_run: 6
exclude_purification_artifacts: true
identity_threshold: 0.30           # clustering cutoff (Data API levels: 30/50/70/90/95/100)
clustering_backend: "precomputed"  # "precomputed" (reuse RCSB) | "mmseqs2" (run our own)
split_fractions: {train: 0.80, val: 0.10, test: 0.10}
split_salt: "snapsplit-v1"         # bump to reshuffle intentionally
seed: 0
# ligand context (downstream-optional; see §4 on separation of concerns)
ligand_context_radius_A: 8.0
max_ligand_atoms: 25
```

## 4. Pipeline stages

Keep a clean separation between what defines the split (entry selection,
clustering, partition) and featurization (context radius, atom caps, cropping).
The split should be model-agnostic; featurization is a consumer of the split.
Emit both, but let downstream users ignore the featurization layer and
re-featurize from the cleaned structures.

**Stage 1 — Enumerate candidates (`enumerate.py`)**

- Query the RCSB Search API (v2) for entries matching: method in
  `experimental_methods`, `resolution <= resolution_max_A`, and
  `release_date <= snapshot_date`. (Verify the current Search API schema/field
  names against `search.rcsb.org` docs before hardcoding — endpoints are stable
  but field paths drift.)
- Then enrich each hit via the **Data API** (GraphQL, batched) — the same call
  proven in the probe: per-entity polymer type + canonical sequence, non-polymer
  `chem_comp` ids/formulas/types, residue counts, assembly info. No coordinates.
- Record per entry: id, release date, resolution, method, polymer entities
  (id/type/sequence), non-polymer components, residue counts. Write
  `candidates.jsonl`. This file is the snapshot definition and carries
  everything Stages 3–6 need.

**Stage 2 — Download (`download.py`) — OPTIONAL, featurization-only**

> Not part of `build`. The split is computed without coordinates (§1.5). This
> stage exists so a downstream consumer can fetch coordinates on demand for the
> entries in a split — e.g. to extract ligand context (Stage 4) or feed a model.

- Fetch mmCIF (not legacy PDB — large/modern entries have no legacy PDB file).
  If `use_biological_assembly`, fetch assembly 1 (`-assembly1.cif`). Use RCSB
  file download services. Cache by ID; store a SHA-256 of each fetched file.
- Respect rate limits; make downloads resumable from the cache. Driven by an
  explicit `if-split fetch` (or the loader), never implicitly by `build`.

**Stage 3 — Filter (`parse.py`)**

- Operate on the metadata in `candidates.jsonl` (no coordinate parsing). Apply
  filters: drop entries with `total_residues >= max_total_residues` (use the
  assembly residue count when `use_biological_assembly`), drop entries whose only
  non-protein components are in `excluded_het`, drop entries with no protein
  polymer entity. Record drop reasons + counts.
- Sequences come from the Data API canonical one-letter code, which already maps
  modified residues to canonical parents (e.g. MSE→MET). Coordinate-level
  re-parsing with `gemmi` belongs to the optional featurization path, not here.

**Stage 4 — Ligand classification & context (`ligands.py`)**

- Classify from metadata (chem-comp ids/formulas/types + polymer types), no
  coordinates:
  - **nucleotide**: entities whose `rcsb_entity_polymer_type` is DNA/RNA. In the
    metadata path these are cleanly typed (the probe showed 1A1F's two chains as
    `DNA`), so the ATOM-vs-HETATM gotcha that put nucleotides out of scope for
    the UMA-Inverse parser *dissolves* — we never touch ATOM records here.
  - **metal**: non-polymer comps whose formula is composed *only* of metal
    element(s) (derived from the `chem_comp` formula). A cofactor like HEM
    (`C34 H32 Fe N4 O4`) contains Fe but also C/H/N/O → it is a small-molecule,
    not a metal ion.
  - **small-molecule**: remaining non-polymer comps after removing waters, the
    `excluded_het` set, and a curated crystallization-additive ignore-list
    (glycerol, PEG, sulfate, EDO, etc.). Distinguishing biologically relevant
    ligands from buffer junk is the genuinely hard curation problem — start with
    a published additive blacklist and make it config-extensible.
- **Confidence tiering — annotate, never destroy (implemented).** Curation runs
  *before* classification: every non-protein component is tagged into one of
  three tiers with a machine-readable reason, instead of silently dropping
  anything.
  - `functional` — real ligand/site: appears in
    `rcsb_entry_info.nonpolymer_bound_components` (actually contacts the protein)
    **or** has a measured `rcsb_binding_affinity`.
  - `ambiguous` — present but uncorroborated (e.g. an unbound metal/ligand);
    reported per-class but **not** given a class label.
  - `artifact` — additive/buffer blacklist, a monatomic counterion (Na⁺/Cl⁻/…),
    or a His-tag/Ni|Co purification metal (reason `histag_metal`).
  Class labels (`metal`/`small_molecule`) derive from the `functional` tier only;
  `nucleotide` = DNA/RNA chains. **No structure is ever dropped for ligand
  quality** — a protein with a junk ion is still a good training backbone; we
  just don't label the junk. Live-verified: 101M → `{HEM: functional, SO4:
  artifact}`; 102L → `{BME: artifact, CL: artifact}`.
- **Purification-artifact curation (the LigandMPNN metal blemish).** LigandMPNN's
  metal test set included structures where the only "metal site" was a poly-His
  purification tag chelating Ni/Co — an artifact of IMAC purification, not
  biology. The tier rule demotes Ni/Co to `artifact` only when it's the entry's
  *sole* metal **and** a protein chain carries a His-run ≥ `histag_min_run`
  (default 6); a genuine catalytic Zn alongside a His-tag is *not* demoted.
  Always recorded in the manifest. Toggle via `exclude_purification_artifacts`.
- Tag each structure with its functional-tier classes, ambiguous classes
  (reported), and all per-component tiers + reasons. The same per-component tier
  is what a downstream featurizer reads to decide real ligand context — the lever
  that improves *training* quality, not just test reporting.
- (Featurization-optional) extract ligand heavy atoms within
  `ligand_context_radius_A` of any protein atom, capped at `max_ligand_atoms`.

> **Test-set stratification (see §6).** Default is a **report-only floor**: the
> pure hash split is untouched; the manifest carries per-split, per-class
> *functional* counts plus *ambiguous* counts so under-representation is visible.
> An opt-in `--enforce-minimums N` top-up (recruit `functional`-only ligand
> clusters into test in deterministic hash order, registry-pinned, shortfall
> logged) is scoped but deferred — tiering first guarantees a quota can only ever
> recruit *functional* ligands, never junk.

**Stage 5 — Sequence clustering (`cluster.py`) — two backends**

Default backend — **`precomputed`** (reuse RCSB clusters):

- **Per-entity cluster ids come from the Data API itself**, not a separate file.
  The `rcsb_cluster_membership { cluster_id identity }` field on each
  `polymer_entity` returns its cluster id at every identity level
  (30/50/70/90/95/100) — verified live (4HHB_1 → cluster 90 at 30%). Stage 1
  already captures this into `PolymerEntity.cluster_ids`, so the default backend
  needs **no extra download at all**: clustering data rides along with the
  snapshot metadata and is locked by the same `candidates.jsonl` hash.
- This supersedes the earlier "download + lock the 21 MB cluster file" plan: the
  cluster file is redundant when the Data API hands us membership per entity. (We
  keep the bulk file only as a fallback if the field is ever unavailable.) Note
  the available levels are 30/50/70/90/95/100 (no 40) — `identity_threshold`
  must map to one of these for the precomputed backend.
- A cluster's split-hash key is the lexicographically smallest entity id among
  the snapshot members of that cluster id (input-independent → see §6 caveat: to
  be fully growth-stable we persist the cluster→split registry, since the Data
  API only reports membership among *current* entities).

Optional backend — **`mmseqs2`** (run our own):

- Pool the snapshot's protein entity sequences. Run mmseqs2 `easy-cluster` at
  `--min-seq-id identity_threshold` with coverage flags matching LigandMPNN's
  intent (record exact flags). Pin and log the mmseqs2 version. Sort inputs for
  determinism.

Both backends output cluster id → member entities, and entry → cluster(s). An
entry may touch multiple clusters via different chains; assign the entry to the
cluster of its longest protein chain (record all) so split assignment is
unambiguous.

**Stage 6 — Split assignment (`split.py`) — the reproducibility core**

- Assign each cluster (not each entry) to a split by deterministic hash:
  `bucket = int(blake2b(cluster_repr_id + split_salt)) ...` mapped to cumulative
  `split_fractions`. Same salt + same cluster IDs → same assignment, forever.
  (See the Stage 6 note in §1 on making `cluster_repr_id` input-independent.)
- Stratify the test set by ligand class so SM/metal/nucleotide are all
  represented (LigandMPNN's test sets are deliberately ligand-containing).
  Implement as: within the test-bucketed clusters, label structures by ligand
  class and report per-class counts; optionally enforce minimum per-class counts
  by pulling additional ligand-bearing clusters into test via the same hash
  ordering.
- Assert invariant: no cluster appears in more than one split. Fail loudly
  otherwise.

**Stage 7 — Manifest & lock file (`manifest.py`)**

Emit two artifacts in `data/out/`:

- `manifest.json` — human-facing: snapshot_date, config hash, tool versions
  (mmseqs2, gemmi), per-split entry lists, per-structure metadata, ligand-class
  tags, per-class test counts, drop log.
- `dataset.lock` — reproduction-facing: the two snapshot anchors. (1) every
  entry ID + obsolescence status (the candidate set is reproduced from the Data
  API by id + `release_date <= snapshot_date`); (2) the **cluster file**: its
  SHA-256, RCSB `Last-Modified`, and a stored copy (≈17 MB, the only sizeable
  artifact). `if-split verify dataset.lock` re-fetches metadata + cluster file
  and confirms the cluster-file hash matches and the entry set is unchanged.
  (mmCIF SHA-256s belong to the *optional* featurization fetch, not this lock.)
- Version the dataset as `IF-Split-<snapshot_date>` (e.g. `IF-Split-2026.05.30`).

**Stage 8 — Loader (`dataset.py`)**

- A thin `Dataset` that reads a manifest and exposes train/val/test views (entry
  ids, ligand classes, entry→cluster map). Featurization stays pluggable; the
  loader carries no coordinates.
- **Cluster-balanced sampling (implemented).** The PDB is heavily redundant
  (thousands of near-identical lysozyme/kinase co-crystals); sampling entries
  uniformly drowns the model in over-represented folds.
  `SplitView.sample_by_cluster(seed)` draws one representative per sequence
  cluster per epoch — deterministic given the seed (stable hash, no global RNG),
  so an epoch is reproducible and varying the seed rotates which member is drawn.
  Bigger *training-quality* lever than perfecting ligand tiers, and free because
  the clusters already exist.

## 5. CLI surface (`cli.py`)

- `if-split build --config config/default.yaml` → runs Stages 1–7, writes
  manifest+lock.
- `if-split verify data/out/dataset.lock` → re-downloads by ID, checks hashes,
  reports drift.
- `if-split stats data/out/manifest.json` → split sizes, per-class test counts,
  identity audit (sanity-check cross-split max identity is below threshold on a
  sample).

## 6. Gotchas to handle (these will bite otherwise)

- **Obsoleted/superseded entries.** A snapshot reproduced a year later may find
  some IDs withdrawn. Record obsolescence in the lock file; `verify` should
  warn, not silently drop, so reproductions are honest about what changed.
- **Cluster-file drift.** RCSB recomputes the cluster files weekly and keeps no
  history, so they are *not* reproducible by URL. This is why we lock a stored
  copy + hash (§1.5); never reproduce a split by re-downloading the live file.
- **DNA/RNA as ATOM records, not HETATM.** This bit coordinate-parsing pipelines
  (the reason the UMA-Inverse paper marked the nucleotide split out-of-scope).
  In the metadata path it's a non-issue — nucleic acids are typed entities
  (Stage 4) — but it *returns* if you take the optional coordinate/featurization
  path, so handle ATOM-record nucleic acids explicitly there.
- **Crystallization additives vs real ligands** — the curation judgment call.
  Ship a default blacklist; make it overridable.
- **Biological assembly vs asymmetric unit** — pick one (config), apply
  consistently; assemblies can duplicate chains (affects length and clustering).
- **Covalent / multi-residue ligands**, metals with coordinating waters, NMR
  entries (excluded by method filter), and entries with missing resolution.
- **His-tag/Ni purification artifacts** masquerading as metal binders — handled
  by Stage 4 curation above. This is a concrete defect inherited examples would
  carry over from LigandMPNN; we detect it from sequence + comp metadata.
- **Determinism of mmseqs2** — clustering can be sensitive to input order and
  version; sort inputs, pin the version, and log flags.

Do not add training-time coordinate noise here (LigandMPNN's 0.1 Å Gaussian
noise is a model-side regularizer applied at load time, not a property of the
split).

## 7. Tech stack

- Python ≥ 3.11, **`uv`** (env + lockfile), **`ruff`** (lint + format).
- `httpx` (RCSB Search + Data API + cluster file), `pyyaml`, `pydantic` (config
  validation). `gemmi` (mmCIF parsing + assemblies) is needed only for the
  *optional* featurization path (Stage 2/4 context, Stage 8 loader). `mmseqs2`
  (external binary, version-pinned) is needed **only** for the optional
  `mmseqs2` clustering backend — the default `precomputed` backend has no
  external-binary dependency. `torch` optional for the loader. `pytest` for tests.

## 8. Build order (phases)

1. **Skeleton + config:** layout, `config.py` with validation + config hashing,
   CLI stub. ← *done*
2. **Enumerate + verify loop:** Stage 1 (Search + Data API → `candidates.jsonl`)
   and the `verify` command first, so snapshot reproducibility is provable before
   any science. Determinism is the riskiest part — nail it early. (Stage 2
   download is optional/featurization and can come last.) ← *done*
3. **Filter:** Stage 3 filters over the metadata in `candidates.jsonl` (residue
   counts, excluded HET, no-protein) with a drop log and counts. (mmCIF parsing
   moves to the optional featurization path.) ← *done*
4. **Ligands:** Stage 4 classification from metadata (chem-comp + polymer types)
   + additive blacklist + nucleotide-as-polymer handling + His-tag/Ni
   purification-artifact curation. ← *done*
5. **Cluster + split:** Stage 5 `precomputed` backend (per-entity Data API
   cluster membership) and Stage 6 hash assignment, with the leakage invariant
   as a hard test. ← *done*
6. **Manifest + stats + loader:** Stages 7 (manifest + split registry) and 8
   (loader). Optional Stage 2 coordinate fetch still pending. ← *done (fetch pending)*
7. **Validation:** regenerate twice from the same config, assert identical
   manifests (✓ `test_manifest_is_deterministic` + live `IDENTICAL_MANIFESTS`);
   regenerate with growth, assert existing clusters didn't move splits (✓
   `test_existing_cluster_does_not_move_when_dataset_grows`, via the registry).
   ← *done*

**Status: all phases complete.** `build` is deterministic given a config (proven
live: two `--limit 60` builds produced byte-identical manifests), `verify`
round-trips the lock, and the no-cluster-leakage invariant is asserted every
build. Remaining optional work: Stage 2 on-demand coordinate fetch + the
`mmseqs2` clustering backend (both explicitly out of the core split path).

Each phase ends with a runnable CLI command and a test. The dataset is "done"
when `build` is deterministic given a config, and `verify` round-trips the lock
file.
