# Changelog

All notable changes to IF-Split are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The **split is always computed from metadata + sequences only** — `build` never
downloads structure coordinates. That invariant holds across every release below.

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
  coverage (SCOP2 ≈ 52%); reproducibility guarantee #1 states the exact split reproduces
  from the lock + `candidates.jsonl`, not `snapshot_date` alone; the `balanced` covariate
  shift (val/test hold only small, rare folds — not comparable to published numbers) and
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
  structural (super)family straddles two splits (not just sequence clusters) when
  `structural_clustering` is on — matching the fold-leakage guarantee. Backed by
  new *negative* tests that construct leaky partitions and prove the guard fires.
- **`single_chain_only`** filter (opt-in): keep only single-protein-entity
  structures — a metadata proxy for the single-chain CATH setup.
- **`build --count`**: preview how many entries the snapshot matches (one fast
  Search API call) before committing to a full build.
- **Manifest observability**: a ligand tier-reason histogram and per-split fold
  coverage — distinct held-out folds *and* the unclassified fraction per split (the
  **residual-leakage ceiling**: entries no CATH/ECOD/SCOP2 taxonomy classifies are
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

A large release: fold-honest splitting, split-output certification, a two-corpus
training model, a metadata-only curation overhaul, and offline re-derivability.

### Added

- **Fold-level structural leakage control** (opt-in `structural_clustering`:
  `off` | `cath` | `ecod` | `scop2`). Same-fold protein chains are union-merged into
  one leakage-safe component in addition to shared sequence clusters, so a fold cannot
  straddle train/test — using RCSB's precomputed CATH/ECOD/SCOP2 classifications
  (metadata only, no coordinates).
- **Balance-aware split strategy** (`split_strategy: balanced`). Caps dominant folds
  to train and fills val/test to their *entry* targets from the fold tail, restoring
  ~80/10/10 by entries with thousands of held-out folds. `config/fold-aware.yaml`
  ships the fold-honest recipe (`scop2` + `balanced`).
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
