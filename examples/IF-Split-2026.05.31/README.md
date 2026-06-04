# Example split — IF-Split-2026.05.31

A full-PDB split at a **today's-date cutoff** (snapshot `2026-05-31`), produced
entirely from RCSB metadata — **no structure coordinates downloaded**. The split
itself is committed here: plain lists of PDB ids, KB-to-MB in size. Reproduce it
byte-for-byte with:

```bash
uv run if-split build --config examples/IF-Split-2026.05.31/config.yaml --out data/out
```

## The split (committed)

| File | Size | What it is |
|---|--:|---|
| [`train.json`](train.json) | 1.5 MB | training-set PDB ids (one per line) |
| [`val.json`](val.json) | 108 KB | validation-set PDB ids |
| [`test.json`](test.json) | 100 KB | test-set PDB ids (all of them) |
| [`test/metal_test.json`](test/metal_test.json) | 36 KB | test ids with a functional **metal** site |
| [`test/small_molecule_test.json`](test/small_molecule_test.json) | 28 KB | test ids with a functional **small molecule** |
| [`test/nucleic_acid_test.json`](test/nucleic_acid_test.json) | 8 KB | test ids that are protein↔**nucleic-acid** complexes (DNA/RNA chains) |
| [`manifest.json`](manifest.json) | 4 KB | provenance: config, counts, clustering stats, file index |
| [`config.yaml`](config.yaml) | — | the exact config used (= `config/default.yaml`, cutoff pinned) |
| [`STATS.txt`](STATS.txt) | — | `if-split stats` output |

Each `*.json` split file is a flat JSON array of PDB ids — `grep`-friendly and
loadable in one line:

```python
import json
train = json.load(open("examples/IF-Split-2026.05.31/train.json"))
metal_test = json.load(open("examples/IF-Split-2026.05.31/test/metal_test.json"))
```

## Not committed (bulky, regenerable)

These are produced by the same `build` but kept out of git; distribute via a
GitHub Release / Zenodo if you want them downloadable:

| Artifact | Size | What it is |
|---|--:|---|
| `candidates.jsonl` | ~335 MB | snapshot definition: per-entity sequences + curation signals (the hash anchor) |
| `ligands.tiers.json` | ~24 MB | per-component curation audit trail (tier + reason) |
| `ligands.classes.json` | ~3.7 MB | entry → functional class labels |
| `clusters.json` | ~3.4 MB | entry → sequence-cluster component (for `sample_by_cluster`) |
| `dataset.lock` | ~3 MB | reproduction anchor: embedded config + all entry ids + candidates SHA-256 |

`build` regenerates everything from `config.yaml` (same config hash → identical
output).

## Headline numbers

- **223,422** candidates (X-ray + EM, ≤ 3.5 Å, released ≤ 2026-05-31)
- **214,794 kept** — dropped 3,078 no-protein, 5,550 over-size (≥ 6000 residues)
- **34,222** raw RCSB sequence clusters @ 30% identity → **19,587 leakage-safe
  components** after union-find merged **38,834** multi-chain bridging entries

| Split | Entries | Components | Component % |
|---|--:|--:|--:|
| train | 188,672 | 15,613 | 79.7% |
| val   |  13,726 |  1,993 | 10.2% |
| test  |  12,396 |  1,981 | 10.1% |

**The split is balanced on sequence *components*, not entry counts** — that's why
train holds ~88% of entries (redundant families like lysozyme carry many entries
per component). Splitting on components is what prevents cross-split leakage; use
`SplitView.sample_by_cluster()` to draw a de-redundified, one-per-cluster epoch.

## Curation highlights (holo-gated, annotate-never-destroy)

- **Test set, functional tier:** metal 3,984 · small-molecule 3,530 · nucleic-acid 586
- **Test set, ambiguous (reported, not labelled):** small-molecule 3,596 · metal 187 · nucleic-acid 1
  - The small-molecule ambiguous count ≈ the functional one: roughly half of
    bound-looking small molecules aren't corroborated by contact or a measured
    affinity, so they're flagged rather than silently labelled.
  - The metal ambiguous count includes lone, uncorroborated Ni/Co (likely IMAC
    artifacts whose His-tag is absent from the deposited sequence) — demoted from
    functional, not dropped.
- **415** His-tag/Ni(Co) purification artifacts flagged and demoted from the metal
  class — the LigandMPNN metal-set blemish, caught automatically (full His run or
  a partial terminal tag).

Every structure stays in its split regardless of ligand quality; only the labels
and confidence tiers change.
