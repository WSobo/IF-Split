# IF-Split pinned PDB snapshot, 2026-07-22

DOI: [10.5281/zenodo.21782554](https://doi.org/10.5281/zenodo.21782554)

> **Depositing?** Paste `DESCRIPTION.txt` into Zenodo's description box. It is the
> same content as this file with no Markdown, no tables and no column alignment, so
> it reads correctly whether the box renders in a proportional or monospace font.
> This `README.md` is for repository readers.

The metadata snapshot and split behind *Fold Leakage in Inverse-Folding Benchmarks,
and the Limit on Fixing It: Auditing the ProteinMPNN and LigandMPNN Splits*.

Everything here is **metadata only** — sequences, cluster memberships, fold and domain
family assignments, ligand signals. No structure coordinates. That is the point of the
method: the split is computed without downloading a single mmCIF.

Code: https://github.com/WSobo/IF-Split (MIT), tag `v0.6.1`.

## Contents

| file | size | what |
|---|--:|---|
| `candidates-annotated.jsonl` | 471 MB | the pinned snapshot: 225,653 entries, one JSON record each, **including Pfam/InterPro domain families** |
| `dataset.lock` | 3.2 MB | pins the candidate set (sha256), the full config, and the split output hash |
| `manifest.json` | 5.4 KB | the reported statistics for this split |
| `train.json` / `val.json` / `test.json` | 1.7 MB | the split, as entry-id lists |
| `test/{metal,nucleic_acid,small_molecule}_test.json` | | the ligand-class-stratified test sets |
| `config.yaml` | 890 B | the exact config that produced this split (`config_hash f7d4203586df3dc7b10d2948e76d20d8`) |
| `SHA256SUMS` | | checksums for all of the above |

## Reproduce the split from the snapshot

```bash
pip install if-split==0.6.1     # or: git clone && uv sync, at tag v0.6.1
if-split resplit --candidates candidates-annotated.jsonl --config config.yaml --out ./out
if-split verify dataset.lock --candidates candidates-annotated.jsonl
```

The last command re-derives the split offline and checks it against the lock. Expect:

```
OK: reproduced exactly (225653 entries; candidates + split verified).
```

It exits non-zero and prints `DRIFT detected` if either the snapshot or the derived
split differs, so it is safe to use as a CI gate.

## The split

On this snapshot, with all five merge authorities (CATH, ECOD, SCOP2, Pfam, InterPro)
and the data-derived `maximal` strategy:

| | entries | leakage-safe components |
|---|--:|--:|
| train | 213,890 | 1 |
| val | 1,467 | 813 |
| test | 1,465 | 864 |

Training is a **single component**. That is not a configuration choice, it is the
finding: merging on shared families collapses the PDB into one component holding 98.6%
of entries, so at most 1.35% can be held out at the fold level, and this split takes
essentially all of it. Of the 2,932 held-out entries, 888 are *certified* novel — they
carry a family under at least one authority and share none with training. The remaining
2,044 carry no annotation under any authority and are fold-disjoint but unverifiable;
they are reported separately rather than counted as novel.

## Caveats worth reading before you use this

- **`candidates.jsonl` vs `candidates-annotated.jsonl`.** Only the annotated file is
  here, deliberately. An unannotated snapshot carries no Pfam/InterPro families and
  cannot reproduce anything above; it also fails `verify` against this lock.
- **Snapshot-dependent, in a known direction.** CATH, ECOD and SCOP2 are curated
  retrospectively, so their coverage grows over time. As they annotate the recent PDB,
  more chains acquire merge edges and the fold-disjoint residual *shrinks*. The 1.35%
  here is an upper bound on what a later snapshot will yield.
- **Domain-family disjointness is not fold disjointness.** Sharing an InterPro family
  is a claim about a domain, not a fold. The holdout is disjoint under the criteria
  named; residual leakage below those criteria is bounded and disclosed, not zero.

## Citation

Cite the preprint for the analysis and this DOI for the snapshot. The lock embeds the
full config, so a reader who has only `dataset.lock` can recover the exact settings
without this README.
