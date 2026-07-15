# Example split — IF-Split-2026.07.14

A full-PDB split at a **today's-date cutoff** (snapshot `2026-07-14`), produced
entirely from RCSB metadata — **no structure coordinates downloaded**. The split
itself is committed here: plain lists of PDB ids, KB-to-MB in size. Reproduce it
byte-for-byte with:

```bash
uv run if-split build --config examples/IF-Split-2026.07.14/config.yaml --out data/out
```

## The split (committed)

| File | Size | What it is |
|---|--:|---|
| [`train.json`](train.json) | 1.5 MB | training-set PDB ids (one per line) |
| [`val.json`](val.json) | 108 KB | validation-set PDB ids |
| [`test.json`](test.json) | 98 KB | test-set PDB ids (all of them) |
| [`test/metal_test.json`](test/metal_test.json) | 32 KB | test ids with a functional **metal** site |
| [`test/small_molecule_test.json`](test/small_molecule_test.json) | 50 KB | test ids with a functional **small molecule** |
| [`test/nucleic_acid_test.json`](test/nucleic_acid_test.json) | 4.6 KB | test ids that are protein↔**nucleic-acid** complexes (DNA/RNA chains) |
| [`manifest.json`](manifest.json) | 4 KB | provenance: config, counts, clustering stats, file index |
| [`config.yaml`](config.yaml) | — | the exact config used (= `config/default.yaml`, cutoff pinned) |
| [`STATS.txt`](STATS.txt) | — | `if-split stats` output |

Each `*.json` split file is a flat JSON array of PDB ids — `grep`-friendly and
loadable in one line:

```python
import json
train = json.load(open("examples/IF-Split-2026.07.14/train.json"))
metal_test = json.load(open("examples/IF-Split-2026.07.14/test/metal_test.json"))
```

## Not committed (bulky, regenerable)

These are produced by the same `build` but kept out of git; distribute via a
GitHub Release / Zenodo if you want them downloadable:

| Artifact | Size | What it is |
|---|--:|---|
| `candidates.jsonl` | ~392 MB | snapshot definition: per-entity sequences + curation signals (the hash anchor) |
| `targets.jsonl` | ~34 MB | conditioning-target corpus: one row per functional (or opt-in ambiguous) ligand, with its split, class, tier, and entry — the *what to design for* half of the two-corpus output |
| `ligands.tiers.json` | ~24 MB | per-component curation audit trail (tier + reason) |
| `ligands.classes.json` | ~4.5 MB | entry → functional class labels |
| `clusters.json` | ~3.3 MB | entry → sequence-cluster component (for `sample_by_cluster`) |
| `dataset.lock` | ~3 MB | reproduction anchor: embedded config + all entry ids + candidates SHA-256 + split hash |

`build` regenerates everything from `config.yaml` (same config hash → identical
output); `verify … --candidates candidates.jsonl` re-derives and certifies it offline.

## Headline numbers

- **225,106** candidates (X-ray + EM, ≤ 3.5 Å, released ≤ 2026-07-14)
- **216,323 kept** — dropped 3,096 no-protein, 15 poly-UNK (no usable sequence),
  5,672 over-size (> 5999 residues)
- **34,419** raw RCSB sequence clusters @ 30% identity → **19,661 leakage-safe
  components** after union-find merged **39,359** multi-chain bridging entries

| Split | Entries | Components | Component % |
|---|--:|--:|--:|
| train | 189,933 | 15,676 | 79.7% |
| val   |  13,831 |  1,993 | 10.1% |
| test  |  12,559 |  1,992 | 10.1% |

**The split is balanced on sequence *components*, not entry counts** — that's why
train holds ~88% of entries (redundant families like lysozyme carry many entries
per component). Splitting on components is what prevents cross-split leakage; use
`SplitView.sample_by_cluster()` to draw a de-redundified, one-per-cluster epoch.

## Curation highlights (holo-gated, annotate-never-destroy)

- **Test set, functional tier:** metal 4,135 · small-molecule 6,356 · nucleic-acid 593
- **Test set, ambiguous (reported, not labelled):** small-molecule 735 · metal 186 · nucleic-acid 2
  - Functional small molecules far exceed ambiguous (6,356 vs 735): the
    `is_subject_of_investigation` gate recovers non-covalently bound cofactors
    (FAD/NAD/FMN/NADP) and inhibitors that the bond-based contact field misses.
  - **Glycans** (RCSB CCD `type` = *saccharide*: NAG/BMA/MAN/…, and sugar-detergents
    like LMT/BOG) with no measured affinity are tiered `glycan` — decorative
    N-glycosylation and cryo additives, not a ligand pocket — so they sit in the
    ambiguous small-molecule count rather than inflating the functional one.
  - The metal ambiguous count includes lone, uncorroborated Ni/Co (likely IMAC
    artifacts whose His-tag is absent from the deposited sequence) and heavy-atom
    phasing derivatives (below) — demoted from functional, not dropped. RCSB
    GO/InterPro metal annotations rescue native metalloenzymes back to functional.
- **415** His-tag/Ni(Co) purification artifacts flagged and demoted from the metal
  class — the LigandMPNN metal-set blemish, caught automatically (full His run or
  a partial terminal tag).
- **Heavy-atom phasing derivatives** (Hg/Au/Pt/Pb/Tl/lanthanide MAD soaks) are
  demoted to `ambiguous` (reported, recoverable) unless a measured affinity, curated
  SOI, or a matching metal annotation vouches for them — so they leave the functional
  metal set (into the ambiguous metal count, 186) without being destroyed. Conversely
  **inorganic Fe-S, metal-oxo, and FeMo clusters** (SF4/FES, the Mn₄CaO₅ OEC) are
  classed `metal` rather than small-molecule. A blacklisted additive that is the
  *measured* ligand (a malonate/citrate with a Ki) stays functional.

Every structure stays in its split regardless of ligand quality; only the labels
and confidence tiers change.

## Two training corpora (backbones + conditioning targets)

The same build emits **two** corpora keyed off one split, so an inverse-folding
model can train on scale *and* on quality:

- **Backbones — all 216,323 kept structures.** Every entry is a valid design
  target regardless of ligand quality; this is the *sequence-design* corpus.
- **Conditioning targets — 224,071 functional ligands** (`targets.jsonl`, one row
  per ligand, keyed to entry + split + class + tier). This is the *what to design
  for* corpus: when a structure holds several ligands, you pick the right one at
  training/inference time instead of conditioning on all of them.
  - per split (functional): train metal 67,143 · sm 118,808 · na 10,117 ·
    val metal 5,078 · sm 7,952 · na 556 · test metal 4,650 · sm 9,174 · na 593
  - **18,504 opt-in** ambiguous targets (non-native metal sites + glycans) are
    available via `include_ambiguous=True` for a lectin/glycosidase or
    metal-site consumer, and excluded by default.

`SplitView` in `ifsplit.dataset` exposes both views (`backbones`,
`conditioning_targets(...)`); `manifest.json` carries a `training` summary.
