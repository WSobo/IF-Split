# Auditing the LigandMPNN / ProteinMPNN train–test split

**What the split holds out, and what it lets through.** LigandMPNN and ProteinMPNN
(Dauparas et al., *Science* 2022; *Nature Methods* 2025) partition the PDB by
**30 % sequence-identity clusters only**, assign whole clusters to train/val/test,
and evaluate native-sequence recovery on the held-out test clusters. That stops a
model from having literally seen a >30 %-identical chain. It does **nothing** about
the three things a structure→sequence (inverse-folding) model actually keys on:

- **structural (fold) homology below 30 % identity** — the "twilight zone" where
  sequence similarity is undetectable but the backbone is conserved;
- **shared ligand context** — the same cofactor/metal pocket in train and test;
- **multi-chain / complex bridges** — the same assembly's context on both sides.

This directory reports an audit of the **actual, published** LigandMPNN and
ProteinMPNN splits — not a re-derivation of their logic — using an **independent
structural authority** (RCSB's precomputed CATH / ECOD / SCOP2 classifications and
its own 30 % sequence clusters) applied to the partition they shipped.

> The split criterion, verbatim (LigandMPNN, *Nat. Methods* 2025, Methods →
> "Training data"): *"Protein sequences were clustered at 30 % sequence identity
> cutoff using mmseqs2. We held out a nonoverlapping subset of proteins that have
> small-molecule contexts (a total of 317), nucleotide contexts (74) and metal
> contexts (83)."* The words **CATH, fold, topology, homolog, leak** appear
> nowhere in either paper.

---

## Key findings

Applied to the models' **real** split lists (`train.json` / `valid.json` /
`test_*.json` for LigandMPNN; `list.csv` + `{valid,test}_clusters.txt` for
ProteinMPNN), annotated with RCSB metadata for every entry:

We audit **both** held-out sets, because they fail differently (see below):

| | LMPNN val (7 495) | LMPNN test (469) | PMPNN val (8 015) | PMPNN test (7 444) |
|---|--:|--:|--:|--:|
| Genuinely **novel-fold** held-out structures | **57 (0.8 %)** | **1 (0.2 %)** | **32 (0.4 %)** | **38 (0.5 %)** |
| Fold-seen among fold-classified | 99.2 % | 99.8 % | 99.6 % | 99.5 % |
| Median training structures sharing the fold (SCOP2) | 189 | 665 | 288 | 311 |
| Shares a **leakage-safe fold-component** with train | 67.1 % | 90.6 % | 79.7 % | 72.3 % |
| Co-clusters with train (**independent 30 %** clustering) | 20.9 % | 46.9 % | 34.6 % | 34.2 % |

LigandMPNN test **(fold × ligand-class)** contexts also in training: **91.7 %**.
Metal-test contamination (IF-Split tiering, then checked against the deposited
coordinates): 6 flagged, **5 / 83 (6.0 %) confirmed** non-functional.

**The one-sentence version:** *the sets on which these models are selected and
reported hold out essentially zero folds — the "novel-fold" fraction is 0.2–0.8 %
— so the numbers measure fold and ligand-context **recall**, not generalization to
unseen structure.* The result reproduces on two independently constructed splits
differing 16× in size, so it is a property of the **sequence-only-clustering
methodology**, not of one dataset. The single genuinely novel-fold structure in
LigandMPNN's entire *test* set is **1I3J**, an I-TevI homing-endonuclease
zinc-finger–DNA complex.

### Validation leakage is the more consequential failure

Test leakage inflates a *reported* number after training. **Validation** loss is
computed *during* training and drives early stopping, checkpoint selection and
hyperparameter tuning — so a fold-leaked validation set (0.4–0.8 % novel-fold here)
means every one of those decisions selects the model that best reproduces
training-fold statistics. The model that gets **shipped** is optimized for fold
recall, not generalization.

### The recommendation (a call to the field)

The fix needs no new architecture, only a better-controlled split and an honest metric:
**retrain with a fold-held-out validation set**, **re-benchmark with
per-superfamily-reweighted** recovery reported separately on the novel-fold subset,
**score the split with an authority it was not built from**, and **publish the split**
as a reproducible spec. Note that a fully fold-clean split is not available: see the
metadata ceiling below. We deliberately do *not* perform the retraining/re-benchmarking — it belongs to
the model developers, whose pipelines and published numbers a corrected split would
revise — but we release the fold-seen vs. novel-fold partition so Δrecovery can be
measured directly on existing checkpoints.

---

## Why this is non-circular

The obvious wrong way to "audit" a 30 %-mmseqs2 split is to re-run 30 % mmseqs2 and
check it against itself — that finds nothing by construction. Instead we hold the
partition **fixed** (exactly the entries each model shipped to each split) and score
it with a **different, stronger, external** criterion:

- **fold identity** from RCSB's CATH, ECOD and SCOP2 — three *independent*
  structural taxonomies built by structure comparison, which by design group
  <30 %-identity domains that share a fold (precisely the homology a sequence
  clusterer cannot see). All three agree: 95–99 % of test folds are in train.
- **a second 30 % clustering** (RCSB's DIAMOND-based clusters) — a sanity check on
  the split's *own* stated guarantee.

No coordinates are downloaded; every signal is RCSB metadata, pulled with
IF-Split's own `RcsbClient`, and the leakage measurement reuses IF-Split's
`cluster.build_clusters` union-find and `ligands.classify_components` tiering.

---

## Results in detail

Full numbers in [`lmpnn_audit_summary.json`](lmpnn_audit_summary.json) and
[`pmpnn_audit_summary.json`](pmpnn_audit_summary.json).

### A. Fold leakage — three independent taxonomies agree

LigandMPNN test, fraction fold-seen in train among **classified** entries:
CATH 96.6 %, ECOD 99.3 %, SCOP2 97.2 %. Combined across methods: **467 / 469
fold-seen, 1 novel, 1 unclassified.** By ligand class: metal 82/83 seen (0 novel),
nucleotide 68/69 (1 novel), small-molecule 317/317 (0 novel). ProteinMPNN: CATH
95.2 %, ECOD 99.1 %, SCOP2 96.0 %; **38 / 7 444 novel-fold**.

### B. Memorization pressure — folds seen hundreds of times

Not merely *present* in training — abundant. The **median** LigandMPNN test
protein's SCOP2 superfamily appears in **665** training structures (mean 903, p90
2 157); CATH median 613, ECOD median 222. A model can reach high native recovery on
such a test protein by reproducing the amino-acid propensities of a fold it has
been shown hundreds to thousands of times.

### C. The "30 % cluster" boundary is clusterer-dependent

Under RCSB's own independent 30 % clustering, **46.9 %** of LigandMPNN test entries
(34.2 % for ProteinMPNN) co-cluster with a training entry. This does *not* mean
their mmseqs2 split was internally broken; it means "clustered at 30 % identity" is
not a robust separation boundary — a second reputable 30 % clusterer disagrees
about roughly a third to a half of the held-out set. *(Caveat: RCSB's DIAMOND
clustering, coverage params, and 2026 membership differ from the Dec-2022 mmseqs2
run; read this as "does the split survive an independent 30 % clustering," not as a
bug in their pipeline.)*

### D. Ligand-context leakage — the LigandMPNN-specific channel

LigandMPNN conditions on the ligand, so the relevant question is whether the
*conditioning context* was seen. **91.7 %** of test **(SCOP2 fold × ligand-class)**
contexts already occur in training. At the level of the exact ligand molecule,
44.9 % of test functional-ligand comp-ids also appear as functional ligands in
training — i.e. even when the *molecule* is novel, the *fold + ligand-class* pocket
almost always is not. (ProteinMPNN, whose test is not ligand-curated, shows 84.7 %
comp-id and 74.8 % fold×class overlap.)

### E. Complex / fold-bridge leakage

Rebuilding leakage-safe components with IF-Split's union-find (a component = raw
30 % clusters merged by shared multi-chain entries, optionally by shared fold):
**47.1 %** of LigandMPNN test entries already share a *sequence+complex* component
with a training entry, rising to **90.6 %** once same-fold entities are merged.
ProteinMPNN, split at the **chain** level, additionally places **14 individual PDB
entries** on both sides of its own split: **6** have one chain in train and another
in test (5a20, 5a21, 5im6, 6gk2, 6td6, 6vfi) and **8** span train and validation
(4ejx, 4v6u, 5grs, 5li2, 5oid, 6b5b, 6eny, 7abi). The same deposition sits on both
sides in each case. LigandMPNN's entry-level scheme avoids this specific bug, though
not the duplicate-id one below. Reproduce with
`uv run python scripts/count_published_splits.py`.

### E2. Duplicate ids in LigandMPNN's own lists

The same script checks the release files for internal consistency, and they fail it.
The five lists hold **157,492** ids but only **157,485** distinct ones. Five of the
seven repeats sit in two ligand-class *test* files at once (1qum, 1u3e, 2nq9, 6wdz,
7kii), which is defensible — a structure can genuinely carry two classes, and
IF-Split independently assigns 1qum both `metal` and `nucleic_acid`. It matters only
if the per-class counts are read as a disjoint partition. The other two are not
defensible: **`2zio` and `3olt` appear
verbatim in both `train.json` and `test_nucleotide.json`**. Two of the 74
nucleotide-test structures (2.7 %) are therefore also training structures — direct
entry-level contamination of the set the 50.5 % nucleotide recovery is computed on.

They are also the **only two of the 74 that contain no nucleic acid at all**. The
other 72 each have a DNA, RNA or hybrid chain (median 13 nt, shortest 4); these two
are protein-only:

| entry | what it is | ligands |
|---|---|---|
| `2zio` | pyrrolysyl-tRNA synthetase catalytic domain, no tRNA in the crystal | AYB (Lys-AMP analog), 2PN |
| `3olt` | R513H murine COX-2 with arachidonic acid | ACD, COH, NAG, BOG, EDO |

`2zio` is at least nucleotide-adjacent — AYB carries an adenosine, so a substructure
rule would catch it — but there is no nucleic acid to condition on. `3olt` has nothing
nucleotide about it whatsoever. Since the conditioning signal scored for this set *is*
the nucleic acid, neither can contribute a nucleotide context, and the effective set
is 72.

Three unrelated checks pick out the same pair: the duplicate-id scan, the polymer
composition, and the ligand tiering (§F).

**The other two sets pass the same check**, which is what keeps this a defect in one
list rather than a claim about all three. Read "pass" narrowly: this asks only
whether the class is *present*. The metal set passes it entry for entry and is
still the most contaminated of the three once you ask whether the metal is a real
site (5/83 adventitious, §F).

| set | must contain | checked | result |
|---|---|--:|---|
| `test_nucleotide` | a DNA/RNA/hybrid chain | 74 | **2 fail** (`2zio`, `3olt`) |
| `test_metal` | a metal-bearing component | 83 | clean |
| `test_small_molecule` | any non-polymer component | 317 | clean |

Metals count whether they are bare ions or carried by a cofactor, using IF-Split's
own `METAL_ELEMENTS`. Note `3olt` *would* pass the metal test (its Co-protoporphyrin
carries cobalt) — it is simply filed under the wrong class. Reproduce with
`uv run python scripts/count_published_splits.py --check-composition`.

Note the ids are **lowercase** in the released files; a case-sensitive comparison
against upper-cased RCSB ids will not surface this.

**Why IF-Split can't make this particular mistake.** The root cause looks like the
word: "nucleotide" reads equally as "nucleic acid" and as "ligand containing a
nucleotide", and `2zio` satisfies the second reading. IF-Split's `nucleic_acid`
class is keyed on whether the entry has a DNA/RNA/hybrid **polymer entity**
(`ligands.py`, `has_nucleic_acid = any(e.is_nucleic for e in polymer_entities)`), so
a bound mononucleotide goes to `small_molecule` instead; a second gate requires RCSB
to report a protein↔NA assembly interface, so a nucleic acid that never contacts the
protein is tiered *ambiguous* rather than labeled. Fed `2zio` and `3olt`, both come
back `small_molecule`, and a real protein/DNA complex (`1qum`) comes back
`['metal', 'nucleic_acid']`. That is a claim about one class definition, not about
the tool in general.

### E3. The test set is internally redundant

Separately from leakage across the split, the 469 distinct test entries are not 469
independent measurements: they reduce to **297 distinct protein sequence sets**, so
**36.7%** share a sequence set with another test entry and the largest identical group holds
21 entries. Averaging native recovery per entry therefore over-weights whichever scaffolds
were deposited most often.

This is a property of the PDB rather than a defect peculiar to this split — IF-Split's own
maximal holdout is 23.3% redundant on the same measure, and the training set 44.4% — but it
means the effective size of the benchmark is smaller than its entry count suggests. Reproduce
with `scripts/measure_split_redundancy.py` (see Reproduce, step 4).

### F. Test-set contamination (metal)

IF-Split's metadata ligand-tiering flags **6 of the 83** LigandMPNN metal-test
structures as *not* a functional metal site. The tiering is metadata-only, so every
flag is a proxy; `scripts/verify_flagged_ligands.py` re-checks each one against the
deposited coordinates (the only place in this audit that reads a structure). **Five
of the six hold.**

| entry | tiering said | coordinates say | verdict |
|---|---|---|---|
| `2CFV` | His-tag Ni | Ni ligated by **tag** His, 1.97 / 2.19 Å | confirmed |
| `2NZ6` | His-tag Ni | Ni ligated by **six tag** His, 1.91–2.39 Å | confirmed |
| `3HG9` | His-tag Ni | only contact ≤2.8 Å is *native* His31 (2.74 Å); tag 3.7 Å away | adventitious (PilM is an actin-fold pilus protein), but **not** tag-ligated |
| `1T31` | uncorroborated Co | His/Glu at 2.07 / 2.08 Å — on a serine protease | confirmed |
| `4X68` | uncorroborated Ni | Asp/His at 2.32–2.76 Å — on a serine β-lactamase | confirmed |
| `3I9Z` | unbound Cu | Cu(I) **2.25 / 2.29 Å from the CxxC thiolates** | **withdrawn** |

`3I9Z` is *Bacillus subtilis* CopZ, a copper chaperone holding its cargo. RCSB
records no bound component for it only because the deposition carries no
connectivity records at all, so the bond-based signal has nothing to read.

The measured distances behind both tables are archived in
[`flagged_ligand_verification.json`](flagged_ligand_verification.json) — every
component of all 14 entries, not just the flagged ones, since a Stage 4 gate-order
change moves components in and out of the flagged set while the distances never move.
Its `_provenance` block stamps the `ifsplit` version that produced the `tier`/`reason`
fields (the only part of the file that is not a property of the PDB). Regenerate with
`uv run python scripts/verify_flagged_ligands.py --all-comps --json <path>`.

**The six small-molecule flags are all false positives, and we withdraw them.**
Every one is a deposition whose title names the flagged component as its ligand:
glycine betaine in OpuAC (`2B4L`, 58 contacts ≤4 Å), lauric acid and dodecyl
sulfate in the β-lactoglobulin calyx (`3UEU`, `4GNY`), ribose-1,5-bisphosphate in
its isomerase's active site (`5YFS`/`5YFT`, >320 contacts), tetrahydronaphthalenol
in the ERRγ pocket (`6I67`). Two were called unbound because RCSB records no *bond*
for a non-covalent binder; two because a fatty acid and a detergent are on the
additive blacklist that a lipocalin exists to bind; two because a phosphorylated
sugar carries a carbohydrate CCD type. Stage 4 ranked those comp-id rules above RCSB's
curated ligand-of-interest flag, and on this set that trade is wrong 6/6 — every comp
id is an additive or a sugar *somewhere*, so a comp-id rule convicts precisely the
depositions where the comp is the point. **v0.6.2 puts that flag above the additive
blacklist** (fixing `3UEU`/`4GNY`) but deliberately not above the glycan gate, where
the flag covers 85.7% of sugars and a wrong call is only `ambiguous` — so `5YFS`/`5YFT`
stay tiered `glycan`, reported and one `include_ambiguous=True` away rather than
discarded. `2B4L` and `6I67` carry no curated flag at all, so no metadata signal
reaches them (they are the bond-based blind spot that also produced the `3I9Z` error).
The small-molecule test set is clean: **317 / 317**. This was a precision limit of our
tier's flagged tail, not a defect in LigandMPNN's list — the tier still called 311
of 317 functional, and it never removes a structure from a split.

> This is *complementary* to — not the same as — the crystallization-additive
> contamination reported for the metal set in the UMA-Inverse work (e.g. **1F35 /
> 1JOB**, Zn from 200 mM zinc-acetate). IF-Split's metadata tiers those Zn ions as
> functional (they *are* bound); catching them needs the deposition's
> crystallization conditions. The two methods flag *different* contaminated
> entries, so the true metal-test contamination is at least the union of both.

---

## Honest scope — what this does and does not show

This audit shows the **benchmark cannot detect novel-fold generalization** and that
aggregate/near-ligand recovery is dominated by folds and ligand-contexts seen many
times in training. It does **not** show that the models "only memorize":

- **The models demonstrably generalize.** ProteinMPNN sequences *de novo* backbones
  (RFdiffusion binders, hallucinated oligomers) that are in no training set, with
  experimental validation. Pure memorization cannot explain that.
- **High recovery is largely physics.** Buried-core residues are near-determined by
  local packing (ProteinMPNN core recovery 90–95 %, surface ~35 %); Rosetta, an
  energy function with no "memory," reaches 32.9 %. Elevated aggregate recovery is
  expected of any competent method.
- **Fold-aware splits move the number only modestly — but nobody has measured it
  for these models.** ESM-IF's *own* leakage audit found that even a CATH-*topology*
  holdout still had **54 % of test structures matching a training structure at
  TM-score > 0.5**, and TM-filtering dropped recovery only ~2 pp (42.2→40.4,
  51.6→49.5). A **sequence-only** split (MPNN) leaks strictly more than a topology
  holdout, so the inflation is plausibly larger — but its magnitude for
  LigandMPNN/ProteinMPNN has never been reported. **That Δrecovery is the experiment
  this audit motivates** (run the model on the fold-seen vs. novel-fold partition
  this audit produces).
- **Native recovery is arguably the wrong metric** anyway; the field increasingly
  prefers self-consistency / designability. "High recovery = memorization" attacks a
  number practitioners already discount.

The defensible, primary-sourced claim is therefore: **the reported recovery
overstates generalization to novel folds and is dominated by redundant, abundant
families; the test set holds out ~0 % of folds and is not a generalization
benchmark.** The design tool remains useful; the *benchmark* is the problem.

The PDB's fold skew makes an aggregate number especially misleading: ~half of
non-redundant CATH domains fall into four "superfolds" (Rossmann, αβ-plait, TIM
barrel, immunoglobulin) and the top-20 fold groups cover 46 % of them (Cuff et al.
2009). A **per-superfamily-reweighted** recovery is the honest metric — and the one
nobody publishes.

---

## The fix (why IF-Split exists)

IF-Split's `structural_clustering: scop2` + `split_strategy: balanced`
("fold-aware") holds **992 distinct SCOP2 families** entirely out of training (473 test
+ 519 val) by union-merging same-(super)family entities before assigning splits, while
staying leakage-safe (whole components), growth-stable, and metadata-only.

**It narrows the channel this audit exposes; it does not close it.** Measured on one
snapshot with fold-merging as the only variable, scored by ECOD (an authority neither
split was built from): a sequence-only split is 99.18% ECOD-fold-seen in test, the
fold-aware split 98.62% — 0.56pp, though it does isolate 50% more novel-fold test
entries (153 -> 229) across 97 more distinct folds (573 -> 670). Only 12.2% of its test entries carry
a SCOP2 label at all. The same fold classification used to *audit* the leak is used to
*constrain and measure* it — elimination is not available from metadata (see PLAN.md,
"the PDB collapses into one giant component").

---

## Reproduce

```bash
# 1. Download LigandMPNN's published split + fetch RCSB metadata for every entry
#    (metadata only; ~157k entries; resumable).
uv run python scripts/fetch_external_split.py --out /tmp/lmpnn_audit --name lmpnn

# 2. Measure the leakage.
uv run python scripts/audit_ligandmpnn_split.py \
    /tmp/lmpnn_audit/lmpnn_candidates.jsonl \
    /tmp/lmpnn_audit/lmpnn_splits.json \
    /tmp/lmpnn_audit/lmpnn_audit_summary.json

# 3. Recount both published splits from their release files (no RCSB needed).
uv run python scripts/count_published_splits.py

# 4. Internal redundancy (E3): entries vs distinct proteins vs components, for the
#    published test set and for IF-Split's own holdout, on the same measure.
uv run python scripts/measure_split_redundancy.py \
    --candidates data/run-2026.07.22/candidates-annotated.jsonl \
    --split-dir data/rs-recommended \
    --ids "lmpnn_test=/tmp/lmpnn_audit/test_small_molecule.json,\
/tmp/lmpnn_audit/test_nucleotide.json,/tmp/lmpnn_audit/test_metal.json"
```

Step 3 is the check on everything this audit says about the split *files* rather
than about the structures in them: the train/valid/cluster totals (which neither
paper states), the duplicate ids, and the entries that appear on both sides of a
split. ProteinMPNN's chain-level split (`list.csv` + `{valid,test}_clusters.txt`,
inside the 47 MB `pdb_2021aug02_sample.tar.gz` on
`files.ipd.uw.edu/pub/training_sets/`) is rolled up to entry level there the same
way it is here: an entry is held out if any of its chains is, test over validation.

## Sources

- Dauparas et al., "Robust deep learning–based protein sequence design using
  ProteinMPNN," *Science* 378:49–56 (2022). Test recovery **52.4 %** (402 monomers).
- Dauparas et al., "Atomic context-conditioned protein sequence design using
  LigandMPNN," *Nature Methods* 22:717–723 (2025). Near-ligand recovery: small
  molecule **63.3 %**, nucleotide **50.5 %**, metal **77.5 %** (test 317/74/83).
- Hsu et al., "Learning inverse folding from millions of predicted structures"
  (ESM-IF), ICML 2022 — CATH **topology-level** split; Appendix B leakage audit
  (54 % of test at TM > 0.5 to train).
- Rost, "Twilight zone of protein sequence alignments," *Protein Eng.* 12:85–94
  (1999). Cuff et al., "The CATH classification revisited," *Structure* 17:1051–1062
  (2009) — PDB superfold over-representation.
