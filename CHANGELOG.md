# Changelog

All notable changes to IF-Split are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The **split is always computed from metadata + sequences only** — `build` never
downloads structure coordinates. That invariant holds across every release below.

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
  ~80/10/10 by entries with thousands of held-out folds. `config/masterclass.yaml`
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
