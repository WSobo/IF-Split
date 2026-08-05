# Changelog

All notable changes to IF-Split are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The **split is always computed from metadata + sequences only** — `build` never
downloads structure coordinates. That invariant holds across every release below.

## [0.6.3] — 2026-08-04

Fixes a growth-path bug in the `maximal` strategy and makes `config/certified.yaml`
reproduce the split the preprint reports. The archived 2026-07-22 split is unchanged and
still verifies byte-identically against its lock, so no figure in the preprint moves.

### Fixed

- **A registry pin can no longer hold out the dominant component under `maximal`.** The
  holdout ceiling now beats a registry pin: a component that does not fit the budget is
  capped to train whatever it was pinned to. Previously `_maximal_assign` honored pins
  unconditionally before capping, so as the snapshot grew and the giant absorbed
  previously-held-out components it inherited their `test` pin by `test > val > train`
  precedence and was placed in the holdout. Replaying a 2023-cutoff registry on the
  2026-07-22 snapshot produced **train=1,164 / val=862 / test=214,796**, and reported it as
  `growth_stable: true` with `pinned_reassignments: 1`. `check_no_leakage` passed
  throughout, correctly: the split was component-consistent, just inverted. The same
  rebuild now yields the expected 213,890 / 1,467 / 1,465.

  Only the `maximal` path with a registry is affected; `hash` and `balanced` are unchanged,
  and the archived 2026-07-22 split still verifies byte-identically against its lock.

- **`growth_stable` is derived from the outcome, not from the presence of a registry.** It
  was `strategy not in (balanced, maximal) or bool(registry)`, which asserted the property
  a registry is *meant* to provide rather than the one it delivered. A new
  `splits.pinned_entries_reassigned` counts overridden pins in **entries** (one override on
  the dominant component moves most of the corpus while `pinned_reassignments` reads 1),
  `growth_stable` keys off it, and Stage 6 prints a warning naming the entry count.

### Changed

- **`config/certified.yaml` now reproduces the published split.** It pinned
  `snapshot_date: 2026-05-30`, `min_modeled_residues: 20`, cryo-EM 3.0 Å and
  `min_em_backbone_inclusion: 0.7`, none of which match the preprint's numbers or the Zenodo
  deposit, while asserting those numbers in its own header comment. It now builds
  `config_hash f7d4203586df3dc7b10d2948e76d20d8` (train 213,890 / val 1,467 / test 1,465).
  The three quality knobs move to a commented block that states they break hash
  reproduction, and the header records the holdout's measured composition.

### Documentation

- **The README described a default the tool no longer ships.** It called `hash` the default
  strategy and `off` the default `structural_clustering`; `config/default.yaml` has shipped
  `maximal` and `all` since v0.6.0. Fixed in the reproducibility section, the recipe table,
  the stage table and the config reference, and `maximal` is now documented alongside the
  other two strategies rather than omitted.
- **New "Known limits" section**, stating up front what users otherwise hit by surprise:
  identity thresholds are restricted to RCSB's precomputed levels (no 25% or 40%), RCSB's
  cluster file is not identity-complete, reproducible-from-a-lock is not the same as stable
  under growth (a `maximal` rebuild retains 43.7% of test), the reproducibility anchor is
  `candidates.jsonl` rather than `snapshot_date`, a `maximal` holdout is length-skewed, a
  fold guarantee is only as wide as its authority's coverage, PDB only, and the tool stops
  at split lists.
- Stale figures corrected: `PLAN.md`'s co-occurrence giant 43.2% → **44.1%**, its
  frontier series → **98.9% (≤2019) / 85.7% (2024) / 60.3% (2025) / 41.8% (2026)** with the
  population named, and 63.1% → **63.3%**; the audit README's "95 more distinct folds" →
  **97 (573 → 670)**, which is what the paper and a re-derivation both give.
- `CITATION.cff` was left at 0.6.2 by the version bump; now 0.6.3 with its release date.
- `schema.py` now carries an NB documenting the open ECOD empty-F-group bug (entries whose
  ECOD annotation has an empty `name`, e.g. `3OLT`, are counted unclassified). Measured
  upper bound on the 2026-07-22 snapshot: 28,239 of 437,572 protein entities (6.5%) lack
  ECOD but carry CATH/SCOP2, so a full recovery would lift ECOD coverage 77.4% → at most
  83.9%. Left unfixed on purpose: the true figure needs a live RCSB sample, and the fix
  moves the published ceiling table.

## [0.6.2] — 2026-08-03

Class labels change (a rescued additive/glycan now sets `small_molecule`), so this needs a
new version for the same reason 0.6.1 did: same `config_hash`, different code, different
output. Split *assignment* is untouched — tiers only ever set labels, never drop entries —
so no partition moves and no leakage invariant is affected.

### Changed

- **RCSB's curated `is_subject_of_investigation` (SOI) flag now outranks the additive
  blacklist** in Stage 4 (`ligands.py`). It does **not** outrank the glycan gate; the new
  order is affinity → glycan → SOI → blacklist. Previously both comp-id rules fired before
  SOI, on the argument that the flag is noisy. Checking the tier's own output against
  deposited coordinates showed that convicts the depositions where the "additive" *is* the
  ligand: lauric acid and dodecyl sulfate in the β-lactoglobulin calyx (`3UEU`, `4GNY` — a
  lipocalin's function is binding fatty acids and detergents). Both were among six
  LigandMPNN small-molecule test entries the tier wrongly flagged as contaminated.

  The two gates are treated differently because their error costs and their noise rates run
  in opposite directions, measured on a 3,000-entry stride sample:

  | gate | SOI prevalence | a wrong call is |
  |---|---|---|
  | additive blacklist | 4.1% (87/2,114) | `artifact` — **never** emitted to `targets.jsonl` |
  | glycan | 85.7% (209/244, 106 of them NAG) | `ambiguous` — emitted, one `include_ambiguous=True` away |

  So SOI buys the most where it is most informative, and is kept out of the gate where it
  is nearly uninformative and where being wrong costs little. `5YFS`/`5YFT`
  (ribose-1,5-bisphosphate, a phosphorylated sugar substrate with a saccharide CCD type)
  therefore stay tiered `glycan` — reported and recoverable, not discarded.

### Added

- **`scripts/verify_flagged_ligands.py`** — re-checks every ligand-tiering flag against the
  deposited mmCIF (contacts ≤2.8 Å for metals, ≤4.0 Å for organics; coordinating residues
  inside a terminal poly-His run are marked `[tag]`). A *validation* script only: `build`
  still never downloads coordinates. This is what caught the ordering bug above, and what
  cleared `3I9Z` — its Cu(I) sits 2.25/2.29 Å from the CxxC thiolates of a copper chaperone,
  and was called `metal_unbound` only because that deposition carries no connectivity
  records for the bond-based signal to read. Two of the six small-molecule false positives
  (`2B4L`, `6I67`) share that blind spot and are **not** fixed by this release: they carry
  no SOI flag at all, so no metadata signal reaches them. Net: of the six, two are fixed in
  code, two are `glycan`-tiered by design (recoverable), and two are a disclosed limit of
  the metadata-only approach.

## [0.6.1] — 2026-08-03 (tagged `v0.6.1`, merge b8f6646)

Split output changes (the unclustered-chain keying below), so this **must** carry a new
version even pre-release: two builds with the same `config_hash` but different code would
otherwise share a version string and silently produce different splits.

### Added

- **Pfam/InterPro as merge authorities** (`structural_clustering: pfam | interpro | all`,
  and `fold_benchmark_method` likewise). CATH/ECOD/SCOP2 are curated retrospectively, so
  their coverage lags deposition badly — 0.8% of pre-2020 entries are unclassified versus
  36.5% of 2025 and **55.3% of 2026** releases — which means a holdout built from them alone
  is "fold-disjoint" largely because the databases have not caught up. Measured: **75.9% of
  such held-out entries share a Pfam/InterPro family with train**, leakage entirely invisible
  to the structural authorities, and **no cutoff avoids it** (60–87% in every chain-length and
  release-year stratum, *rising* with length). Pfam/InterPro are HMMs over sequence, so they
  carry no lag and see homology below 30% identity: InterPro alone covers **94.9%** of entries
  against 92.3% for all three structural authorities combined. Costs one extra field
  (`annotation_id`) on a Data API request already made — still metadata-only, no coordinates.
  `scripts/fetch_domain_annotations.py --apply` back-fills an existing `candidates.jsonl`, so
  adopting this needs no Stage-1 re-enumeration.
- **`split_strategy: "maximal"`** — size the holdout from the data instead of demanding a
  ratio. Leakage-safety forces the giant component into one split; putting it in train
  maximizes train and leaves the tail free, so the largest leakage-safe holdout *is* the tail.
  `split_fractions` becomes a **ceiling**, never a target, and val/test fill toward whichever
  is smaller — fixing the `balanced` failure where a thin tail empties val entirely (measured:
  `all` + `balanced` gives val=0, test=2,932). A component is capped to train exactly when it
  cannot fit the holdout budget, which is self-scaling and needs no fixed dominance threshold.
- **`config/certified.yaml`** — `all` + `maximal`, the most trustworthy holdout the metadata
  supports. On 2026-07-22: **train 213,890 (98.6%) / val 1,467 / test 1,465**, with 888
  held-out entries certified novel under all five authorities and **0% domain-family leakage
  into train** (down from 75.9%). The holdout is ~1.35% because that is how much genuinely
  novel structure the PDB contains, not a tuning choice.
- **`structural_clustering: "union"`** — merge on **any** of CATH/ECOD/SCOP2 (namespaced),
  the strictest fold control the metadata can express and the highest coverage. Still
  metadata-only: all three authorities are already captured per entity, so this is a Stage-5
  recombination — no new fetch, no coordinates. It is a **measured-ceiling diagnostic**, not
  a production config: it merges the most and so percolates the most (see the fold-graph
  percolation note in the README/PLAN), leaving too thin a tail to fill val/test.

### Fixed

- **Exact-sequence identity now merges, regardless of RCSB's cluster ids** (real leak, found
  by measurement). A *clustered* chain was identified **solely** by its RCSB 30% cluster id,
  and RCSB's cluster file turns out not to be identity-complete: byte-identical sequences can
  carry **different** cluster ids. The same protein could therefore straddle two splits.
  Measured on the 2026-07-22 snapshot: **74 protein sequences across 497 entries** straddled,
  **38 of them contaminating train**, including a **621-residue** chain in test *and* val and
  a 532-residue chain in test *and* train. Every protein chain with at least
  `MIN_UNCLUSTERED_MERGE_MODELED` modeled residues now also keys by its sequence hash, so
  identical chains always co-key. Safe by construction — exact identity is a strict subset of
  30% identity, so the edge can only merge what a correct 30% clustering would already have
  merged; measured effect on the sequence-only build is 19,593 → 19,395 components and a
  43.2% → 44.1% largest component. **After the fix the straddle count is 0.** This is the
  project's own recurring failure mode (asserting an invariant at the level the code makes
  true by construction — cluster ids — rather than the level users depend on: sequences), so
  `check_no_leakage`'s guarantee is now stated over *sequences* above the modeled gate, and
  explicitly not guaranteed below it.
- **Unclustered-chain keying, gated by modeled content.** An unclustered protein chain is
  keyed by its sequence hash so two entries sharing an identical such chain co-key and cannot
  straddle splits. A *fully*-unclustered entry always keys + merges (bounded — its component
  holds only entries that are entirely that sequence). An unclustered chain *inside an
  otherwise-clustered entry* adds a merge edge only when it carries at least
  `MIN_UNCLUSTERED_MERGE_MODELED` (**12**) modeled (non-'X') residues, so an unmodeled or
  low-complexity fragment cannot fan out into a spurious mega-component (catastrophic under
  `hash`, where it lands in a salt-chosen split). The gate is on modeled sequence **content**
  — intrinsic and growth-stable — never on a snapshot-dependent occurrence count. **Measured
  on the full 2026-07-22 snapshot** (`scripts/measure_unclustered_fanout.py`): the fan-out is
  driven by unmodeled poly-'X' / low-complexity sequence, **not length** (a 72-'X' chain
  bridges 283 clusters), and collapses from 429 to a max of 2 at ≥ 12 modeled residues — a
  clean knee, and a finding about the PDB's unclustered tail.
- **Entry-level rebuild diff (the faithful growth signal).** An in-place `build`/`resplit`
  now reports how many prior entries **changed split** vs the build already in `--out`, and
  how many were **absorbed into train** — the direction registry-free `hash` merges are
  biased toward (the survivor's bucket, and train owns 80% of it), i.e. held-out data eroding
  into train. Aggregate fractions can't detect this (they are conserved by construction — a
  simulation churned ~10% of entries while the 80/10/10 barely moved), so the report is at the
  entry level. Diagnostic only: it reads the prior output, never the assignment, so
  `verify`-from-config-alone and the deterministic manifest are untouched.
- **Honest growth reporting.** `splits.pinned_reassignments` counts merge-overridden pins
  only when a registry is in use; the docs no longer claim the registry-free `hash` path's
  merge migration is "reported" via `pinned_reassignments` or aggregate drift (it isn't — see
  the rebuild diff above). `stats` also warns when a `balanced` split's realized entry
  fractions drift from target (a coarse, balanced-only signal for the `test > val > train`
  ratchet — not a substitute for the entry-level diff). The `examples/IF-Split-2026.07.14`
  README no longer claims a fresh `build` reproduces the split byte-for-byte (that needs the
  locked `candidates.jsonl` + `dataset.lock`) and notes the example is the fold-leaky default.

## [0.6.0] — 2026-07-24

A correctness-hardening release: five confirmed silent-failure fixes (each found by
probing the real split functions, each now guarded by a red-without-the-fix test),
plus the `if-split init` config wizard. The default `hash` split output is unchanged
except the singleton-keying fix below, which only affects unclustered short peptides.

### Added

- **`if-split init` — a config wizard.** TTY-aware; scaffolds a `config.yaml` from the
  `default` or `fold-aware` recipe, prompts for the highest-signal knobs (snapshot date,
  resolution, split fractions, `structural_clustering`, `fold_benchmark_method`, salt),
  preserves every recipe comment, validates before writing, and **never runs a build**.
  The two recipes are embedded, so it works from an installed wheel (which omits
  `config/`); `--recipe` / `--non-interactive` / `--force` make it scriptable.

### Fixed

- **Growth stability across a merge (the old claim was false).** A later snapshot's
  bridging multi-chain entry can union two prior components into one; the absorbed
  component's entries used to silently follow the survivor's split (a held-out test entry
  could become train), and a registry pin under the vanished key was ignored. The registry
  now matches a pin on **any** key a component covers, so a held-out component stays held
  out across a merge (conflicts resolve `test > val > train`); any unavoidable reassignment
  is **counted** in `splits.pinned_reassignments` and surfaced by `stats`. Affects both
  `hash` and `balanced`; the docs no longer claim `hash` "never moves existing" components.
- **Identical unclustered sequences could straddle splits.** A fully-unclustered protein
  chain was keyed on its entity id, so two entries with the same peptide sequence got
  different singleton components and could land in different splits while `check_no_leakage`
  stayed blind. Fully-unclustered singletons are now keyed on a **hash of the sequence**, so
  two such entries with an identical sequence share one component (a genuine, if small,
  default-config leak — closed for the fully-unclustered case). *Note:* this changes the
  singleton key format, so a `splits.registry.json` from ≤ v0.5.0 does not carry over pins
  for unclustered-peptide components on the first in-place `balanced` rebuild — they re-key
  and may reshuffle (use `--fresh` for a clean lineage). Clustered components are unaffected.
- **`test_min_per_class` could silently blow up a `balanced` split.** The per-class top-up
  could recruit a dominant fold (capped to train) into test to meet a small floor, pushing
  test far past its target with no shortfall reported. The top-up now recruits
  **smallest-sufficient-component first** and never pulls an above-cap fold under
  `balanced`; an unmeetable floor is reported as a shortfall.
- **Novel-fold benchmark could be tautological.** With `fold_benchmark_method ==
  structural_clustering`, the merge already holds every fold out of test, so the novel-fold
  fraction is ~100% by construction. `load_config` now **warns**, the manifest records
  `tautological_with_merge`, `stats` flags it, the wizard notes it, and
  `config/fold-aware.yaml` suggests an independent authority (`ecod`) instead.
- **Packaging.** Dropped the stale `mmseqs2` keyword and dependency comment (the backend
  was removed in v0.4.0).

### Changed / documented

- **ECOD/SCOP2 fold keys are free-text names** (their `annotation_id` is per-domain), so a
  *fresh* re-enumeration could merge a fold differently if RCSB renames a superfamily.
  Documented as a fresh-rebuild caveat (a locked build reproduces exactly via
  `candidates.jsonl`; CATH is stable); the stable-lineage-id fix is scoped for a follow-up.
- **Honesty pass on the claims.** The README feature table now qualifies fold hold-out by
  coverage (SCOP2 covers 47.7% of chains / 61.7% of entries corpus-wide, but only ~12% of
  *test* entries in a fold-aware run — mind the denominator); reproducibility guarantee #1 states the exact split reproduces
  from the lock + `candidates.jsonl`, not `snapshot_date` alone; the `balanced` covariate
  shift (val/test hold small, rare folds plus a majority of fold-unclassified chains — not
  comparable to published numbers) and
  the `scop2`-vs-`ecod` trade-off are stated plainly.

## [0.5.0] — 2026-07-22

Toward "The Novel-Fold Benchmark". No change to the default (`hash`) split output.

### Added

- **Novel-fold benchmark export** (opt-in `fold_benchmark_method: cath|ecod|scop2`).
  Emits the fold-seen vs novel-fold TEST partition as turnkey lists + labels —
  `novel_fold_test.json` (the novel-fold test subset), `fold_groups.json` (per-superfamily
  test groups, for per-family reweighting), and `folds.json` (per-entry fold labels +
  novel-fold flag) — so a model developer can score native recovery on the novel-fold
  subset and per-superfamily-reweighted on an existing checkpoint. Fold *labels* are
  decoupled from fold *merging*, so they attach even to a fold-leaky split (the split a
  checkpoint was trained on) and never change the split or `check_no_leakage`. `stats` and
  the loader (`SplitView.novel_fold_entries()`, `IFSplitDataset.fold_groups()`) expose it.
- **`stats` entry-skew view**: each split prints its entry fraction against the configured
  target (e.g. `train: 95.0% / target 80.0%`), so the entry-balance skew the `balanced`
  strategy corrects is visible. The README Outputs table now lists every build output, and
  the hydrated `DATASET_CARD.md` integrity snippet is filled in.

### Fixed

- **Growth-stability for the `balanced` strategy.** A `balanced` split's val/test
  fill boundaries scale with the snapshot's total entries, so a growing snapshot
  could move a few percent of prior components across train/val/test (including
  train→val/test contamination) unless a registry pinned them — and the CLI never
  self-pinned. An in-place rebuild now auto-adopts `<out>/splits.registry.json` when
  the prior build used the same config (its `dataset.lock` `config_hash` matches);
  `--fresh` opts out. `hash` is unchanged (already input-independent and registry-free,
  so `verify` can still certify it). The manifest records `splits.growth_stable` and
  `stats` prints it.

## [0.4.0] — 2026-07-22 (hardening)

Reliability, correctness-guard, and publication-readiness pass. No change to the
default split output.

### Added

- **Fold-level leakage guard.** `check_no_leakage` now also asserts that no
  structural (super)family *the configured authority classifies* straddles two splits (not
  just sequence clusters) when `structural_clustering` is on. It is blind to unclassified
  chains, which are the majority of a fold-aware val/test. Backed by
  new *negative* tests that construct leaky partitions and prove the guard fires.
- **`single_chain_only`** filter (opt-in): keep only single-protein-entity
  structures — a metadata proxy for the single-chain CATH setup.
- **`build --count`**: preview how many entries the snapshot matches (one fast
  Search API call) before committing to a full build.
- **Manifest observability**: a ligand tier-reason histogram and per-split fold
  coverage — distinct held-out folds *and* the unclassified fraction per split (the
  **residual-leakage ceiling**, set by the *configured* authority: entries it does not classify are
  held out by sequence only, so fold-level hold-out is not guaranteed for them).
  `stats` prints it whenever fold-aware clustering is on.
- **CLI test suite** (`tests/test_cli.py`) covering exit codes and error paths.

### Changed

- **Removed the `mmseqs2` clustering backend.** RCSB's precomputed clusters (the
  same 30% clustering ProteinMPNN/LigandMPNN used, locked via the snapshot) are the
  sole backend. `clustering_backend: mmseqs2` was an unimplemented stub that crashed
  mid-build; it is now rejected at config validation.
- **Robust CLI error handling**: malformed JSON, old-schema files, bad values, and
  network failures now produce actionable one-line messages with documented exit
  codes (2 bad input, 3 not implemented, 4 network, 130 interrupted) instead of a
  traceback. `fetch --workers` is validated `>= 1`.

### Fixed

- **Atomic writes** (temp file + rename) for the manifest, lock, and split lists —
  a crash mid-write can no longer leave a partial file that crashes every reader.
- **Stale per-class test files**: rebuilding into a used `--out` now clears the
  managed `test/` subtree, so a `test/<class>_test.json` can no longer linger with
  an entry that has since moved to train (which read as leakage).
- **Loader fails loudly** on a missing split file instead of silently returning an
  empty (wrong) partition.
- `count_entries` no longer crashes on a zero-match (HTTP 204) Search response.

## [0.3.0] — 2026-07-14

A large release: fold-aware splitting, split-output certification, a two-corpus
training model, a metadata-only curation overhaul, and offline re-derivability.

### Added

- **Fold-level structural leakage control** (opt-in `structural_clustering`:
  `off` | `cath` | `ecod` | `scop2`). Same-fold protein chains are union-merged into
  one leakage-safe component in addition to shared sequence clusters, so a family the
  configured authority names cannot straddle train/test (families it does not classify are
  unconstrained) — using RCSB's precomputed CATH/ECOD/SCOP2 classifications
  (metadata only, no coordinates).
- **Balance-aware split strategy** (`split_strategy: balanced`). Caps dominant folds
  to train and fills val/test to their *entry* targets from the fold tail, restoring
  ~80/10/10 by entries, holding 992 distinct SCOP2 families out of train. `config/fold-aware.yaml`
  ships the fold-aware recipe (`scop2` + `balanced`) — fold *measurement* over the classified
  fraction, not fold-clean: its test set is still 98.6% ECOD-fold-seen.
- **Split-output certification.** The `@2` `dataset.lock` records `split_sha256` (a
  hash of the entry→split partition); `verify` re-derives Stages 3–6 and certifies the
  split *output* reproduced, not just the Stage-1 candidate set.
- **Two training corpora from one split**: all kept structures as design *backbones*,
  plus a functional-ligand *conditioning-target* corpus (`targets.jsonl`, one row per
  ligand keyed to entry + split + class + tier). `SplitView` exposes both views.
- **Offline `resplit`** (`if-split resplit --candidates candidates.jsonl --config X`):
  re-derives Stages 3–7 from a cached snapshot with no RCSB — ablate curation /
  clustering / split settings, or tighten a filter, in seconds instead of
  re-enumerating the PDB. The lock records `source` (`build` | `resplit`).
- **Offline `verify`** (`verify LOCK --candidates candidates.jsonl`): integrity-check
  a distributed dataset with no network; a corrupt candidates file is reported as an
  integrity failure. A `resplit` lock is steered to offline verification.
- **Per-method resolution caps** (`resolution_max_A_by_method`) and a **cryo-EM
  map-fit floor** (`min_em_backbone_inclusion`, wiring in the previously-unused
  `em_backbone_inclusion` metric). Resolution is now re-derived in Stage 3, so the cut
  is auditable from `candidates.jsonl` and tightenable offline.
- **Opt-in sequence-usability floor** (`min_modeled_residues`) and an always-on drop of
  empty / all-`X` (poly-UNK) protein chains, which carry no learnable label.
- RCSB **metal-binding annotations** (GO/InterPro/Pfam) captured to rescue native
  metalloenzymes; `if-split spec` to emit a portable, self-identifying split spec.

### Changed / curation

- **Metal tiering**: heavy-atom / lanthanide **phasing derivatives** (Hg/Au/Pt/Pb/Tl/…)
  demoted to `ambiguous` (reported, recoverable) rather than counted as functional
  metal sites; inorganic **Fe-S / metal-oxo / FeMo clusters** (SF4/FES, the OEC) now
  classed `metal`; native Ni/Co (and heavy/lanthanide) sites rescued via annotation,
  affinity, or subject-of-investigation. The lone-Ni/Co His-tag figure was corrected
  (~96% → ~82%).
- **Glycans** (RCSB CCD `type` = *saccharide*) with no measured affinity are tiered
  `glycan` (decorative / detergent), recoverable via an opt-in tier — not counted as
  small-molecule conditioning targets.
- **Small molecules**: a measured binding affinity now overrides the additive
  blacklist, so a blacklisted comp that is the real measured ligand stays functional.
- **Nucleic acids**: `is_nucleic` now recognizes the `NA-hybrid` polymer type; the
  ligand class was renamed `nucleotide` → `nucleic_acid`.
- The size cap keeps `< 6000` residues correctly (`> max_total_residues`, not `>=`).
- Adding a resolution-less method (NMR/SAXS) now warns instead of silently returning
  zero entries.

### Fixed

- `verify` warns (rather than fails) on a version-only lock mismatch.
- `fetch` reads split id-lists from the manifest directory, not the current directory.
- `identity_threshold` is validated against RCSB's precomputed cluster levels
  (30/50/70/90/95/100) so an unsupported level can't silently disable clustering.
- A bound halide is tiered a counterion, not a functional small molecule.

## [0.2.0] — 2026

- Recover non-covalently bound cofactors (FAD/NAD/FMN/NADP, inhibitors) via RCSB's
  `is_subject_of_investigation` flag.
- Harden Ni/Co metal curation against His-tags absent from the deposited sequence.
- Shareable split spec (`if-split spec`) and a self-identifying config header.
- Rename the ligand class `nucleotide` → `nucleic_acid`; PyPI/CI badges + install docs.

## [0.1.0] — 2026

- Initial release: a reproducible, date-pinned, ligand-aware train/val/test splitter
  for the PDB. Enumerate → filter → tier ligands → cluster (union-find, leakage-safe)
  → deterministic split → manifest + lock, all from RCSB Search + Data API metadata
  (no coordinates). Optional `fetch` downloads structures for a built split.
