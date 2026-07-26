---
name: if-split-development
description: >
  Correctness discipline for working on IF-Split (the reproducible, metadata-only PDB
  train/val/test splitter). Use whenever a change touches split assignment, clustering,
  components, the growth registry, ligand tiering, fold leakage, the novel-fold
  benchmark, growth stability, the manifest/lock, or any README/PLAN/paper claim about
  those — and whenever adding a config knob or reviewing a change or claim in this repo.
  Triggers: "if-split", "ifsplit", "split.py", "cluster.py", "check_no_leakage",
  "structural_clustering", "fold leakage", "growth stability", "splits.registry.json",
  "novel fold", "balanced strategy", "test_min_per_class", "the splitter". Read CLAUDE.md
  first for environment/commands/layout — this file is the part that is NOT written down
  there: how not to be wrong.
---

# Working on IF-Split

`CLAUDE.md` covers environment, commands, and module layout. This covers correctness —
the discipline distilled from six confirmed silent-failure bugs (all fixed in v0.6.0),
every one of which shipped green.

## The one idea everything follows from

Every correctness bug in this codebase was **silent**: tests green, ruff clean, manifest
reporting no shortfall and no gap — split wrong. They share one root cause:

> **The invariant was asserted at the level the implementation makes true by
> construction, instead of the level the user depends on.**

Union-find guarantees a component maps to one split, so an assertion *about components*
can never fail — that is a restatement, not a guarantee. Before writing any check or
test, decide which column it lives in:

| Implementation's level (trivially true) | User's level (what actually breaks) |
|---|---|
| a component key keeps its split | an **entry id** keeps its split across snapshots |
| a raw cluster key is in one split | an identical **sequence** is in one split |
| fold labels used for merging are disjoint | an **independent** authority says test is novel |
| the candidate-set hash matches | the **entry→split partition** matches |
| a per-class floor was requested | the floor was met **and** the balance survived |
| an aggregate fraction / a counter moved | the **entries that changed** are reported |

**This applies to observability, not just invariants — and it is the trap this repo keeps
falling into.** The reporting layer has failed the same way three times: `cluster_split`
instead of `entry_split`, then `pinned_reassignments` (a registry-only counter) instead of
entries-moved, then aggregate entry-fractions instead of entries-moved. Each claim got
narrower and more careful and still landed on a quantity the implementation can compute
cheaply rather than the one a reader needs. Aggregates are the worst offenders: an 80/10/10
is *conserved by construction*, so ~10% of entries can migrate (held-out → train, under
`hash`, the unsafe direction) while the aggregate barely moves. **Specify a diagnostic by
what would go wrong, not by what is convenient to count** — for growth that means the
entry-level rebuild diff (`_report_rebuild_migration`), not a fraction or a counter.

If a check only ever reads the left column, it is decoration. Move it right.

## Hard constraints (settled — do not relitigate in a feature PR)

- **Metadata only in the build path.** No coordinate download in Stages 1, 3–7. `fetch`
  is optional and downstream. Every feature must be expressible in RCSB metadata.
- **No mmseqs2, no foldseek.** RCSB precomputed clustering is the only backend, on
  purpose (no version drift, no build dependency). Do not reintroduce a local binary.
- **Annotate, never destroy.** Ligand quality is a tier + machine-readable reason;
  structures are never dropped for ligand quality, only relabeled.
- **Determinism.** Same config + same `candidates.jsonl` → byte-identical `manifest.json`.
  No wall-clock, no set-iteration order, no float formatting in anything hashed.
- **PDB ids verbatim.** Never slice, length-validate, or case-fold `rcsb_id`.

## The failure taxonomy (verified; #1-#5 fixed in v0.6.0, **#6 still live**)

Check every change against these six patterns. Severities are the *verified* ones —
several were narrower than a first read suggests.

1. **Component merge under growth.** A later snapshot's bridging multi-chain entry unites
   two prior components; the absorbed one's entries follow the survivor's key. This is a
   **cross-snapshot migration, NOT within-build leakage** — `check_no_leakage` still holds
   in any single build. The registry pin under the vanished key was silently ignored.
   *Fix:* pin on **any** key a component covers, resolve conflicts `test > val > train`,
   count overrides in `splits.pinned_reassignments`. Affects hash and balanced.
2. **A benchmark scored with the labels it was split on.** `fold_benchmark_method ==
   structural_clustering` ⇒ ~100% novel by construction — it corrupts a *reported
   statistic*, not the split (labels are decoupled from union-find). (The audit's "~0%
   when off" is data-dependent; a synthetic gave 40% — the ~0% is the real-PDB LigandMPNN
   context.) *Fix:* warn when the two methods match; the shipped config suggests an
   independent authority. Generalize: any metric computed from the signal used to build
   the split measures the config, not the data.
3. **Keys derived from identifiers, not content.** Unclustered singletons keyed on
   entity id let an identical sequence straddle splits while `check_no_leakage` (which
   compares raw keys) stayed blind. Real, default-config leak — but bounded to RCSB's
   unclustered tail (short peptides); normal-length identical sequences co-cluster at 30%.
   *Fix:* key singletons on a sequence hash. Any key standing for "the same biological
   thing" must derive from the thing, not an accession.
4. **Knobs that quietly cancel.** `test_min_per_class` top-up could recruit a capped
   mega-fold into test, blowing `balanced` to ~35% test with `shortfalls={}` and
   `balance_gaps={}`. Needs BOTH balanced AND a floor set (no shipped recipe does). *Fix:*
   recruit smallest-sufficient-first, exclude above-cap folds under balanced, report the
   shortfall. **Every new knob needs a test against each existing knob it can interact
   with** — the single-knob test never catches this.
5. **Coverage ceilings quoted as guarantees.** On the 2026-07-22 snapshot SCOP2 classifies
   47.7% of chains (CATH 38.1, ECOD 71.8, union 77.8); by *entries* 61.7/65.3/87.1/92.3.
   And the corpus rate is the WRONG denominator for a hold-out claim: in the fold-aware run
   only **12.2% of TEST entries** carry a SCOP2 label (train is 73.7%), because classified
   entries are exactly the ones that merged into the capped giant. **Never state a fold
   guarantee without the denominator that matches the set you are describing.**
6. **Keys on mutable strings.** ECOD/SCOP2 fold-merge keys on free-text names (their
   `annotation_id` is per-domain), so a *fresh* rebuild can merge differently if RCSB
   renames a superfamily. A locked build reproduces exactly via `candidates.jsonl`; CATH
   is stable. The stable-lineage-id fix is deferred (choosing the level that preserves
   merge granularity needs full-PDB validation) and documented, not hidden.

## Testing discipline

> **A test that still passes when you delete the feature it tests is broken.**

Before adding a test, delete or stub the code path and confirm it goes red. If it stays
green, it asserts the left column. Specifics for this repo:

- **Growth tests must include a *bridging* entry** (`_protein_record("C", [1, 2])`), not
  only disjoint additions, and **assert on `entry_split`, not `cluster_split`** (a merge
  makes the old component key vanish, so a component-level assertion is skipped, not
  failed). See `test_growth_bridging_merge_is_honest_and_registry_stable`.
- **Interaction matrix.** A test touching `split_strategy` should also run with
  `test_min_per_class` set, and vice versa.
- **Reported diagnostics are part of the contract.** `balance_gaps`, `minimum_shortfalls`,
  `pinned_reassignments`, `tautological_with_merge` — a *wrong silence* is a bug of equal
  severity to a wrong split. Assert the diagnostic fires.
- Keep the suite offline; the one live RCSB test stays behind `IFSPLIT_NETWORK_TESTS=1`.

## Claims discipline (README / PLAN / paper are part of the artifact)

Sort every claim into a bucket and use the matching verb:

| Bucket | Verb | Example |
|---|---|---|
| Guaranteed by construction | "cannot" / "never" | a sequence cluster cannot straddle two splits |
| Bounded by coverage | "bounded by" / "ceiling of" — **with the denominator that matches the set** | in a fold-aware build, 87.8% of *test* entries carry no SCOP2 label, so fold hold-out is guaranteed for ~12% of test |
| Measured | "measured at X on the \<snapshot\>, by \<authority\>, on \<which set\>" | LigandMPNN's *published* test set is 0.2% novel-fold; IF-Split's own fold-aware test set is 1.4% ECOD-novel (229/16,569) on 2026-07-22 |

- **Reproducibility = lock + `candidates.jsonl`, not config/`snapshot_date` alone.** RCSB
  recomputes clusters and CATH/ECOD/SCOP2 annotations over time with no public history.
- **Do not pick a fold method by the answer it gives.** SCOP2 is chosen over ECOD because
  ECOD's merging starves val — report ECOD's numbers alongside; "ECOD can't fill a
  fold-disjoint 80/10/10" is a finding, not a reason to switch authorities silently.
- **`balanced` introduces covariate shift by construction** (val/test hold only small,
  rare folds), so its recovery numbers are not comparable across strategies or to
  published LigandMPNN/ProteinMPNN figures. Say so where balanced is described.

## Review recipe

Read the diff, then read `split.py` and `cluster.py` as they will exist after it — then
**stop reading and start probing.** Construct ~10 synthetic components (or a handful of
`CandidateRecord`s) and loop over `split_salt` against the *real* functions
(`build_clusters`, `assign_splits`, `check_no_leakage`, `build_fold_benchmark`). Every one
of the six bugs fell out of a probe like that; none was visible by reading. For growth,
always add a bridging entry and diff `entry_split` across the two snapshots.

## Adding a config knob

1. Field + validator in `config.py`. Decide explicitly whether it enters `config_hash`;
   omit only when the "off" value must stay hash-compatible with older builds, and say why.
2. Implement the output-affecting logic in the stage module, not `cli.py`.
3. Test the knob alone. Then test it against every existing knob it can interact with.
4. Confirm its manifest diagnostic actually fires when it misbehaves.
5. Update the README config table, PLAN rationale, and CHANGELOG.
6. `uv run ruff check .` and `uv run pytest -q` both clean.
