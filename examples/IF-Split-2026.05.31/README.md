# Example split — IF-Split-2026.05.31

A full-PDB build at a **today's-date cutoff** (snapshot `2026-05-31`), produced
entirely from RCSB metadata — **no structure coordinates downloaded**. This is a
worked example of what `if-split build` emits; reproduce it with:

```bash
uv run if-split build --config examples/IF-Split-2026.05.31/config.yaml --out data/out
```

## What's here

| File | What it is |
|---|---|
| [`config.yaml`](config.yaml) | The exact config used (= `config/default.yaml` with `snapshot_date: 2026-05-31`). |
| [`STATS.txt`](STATS.txt) | `if-split stats` output for the build. |
| [`manifest.summary.json`](manifest.summary.json) | The manifest's **aggregate** sections only (config, filter, clustering, split counts). |

The **full artifacts are not in git** — they're large and regenerable:

| Artifact | Size | |
|---|--:|---|
| `candidates.jsonl` | ~689 MB | snapshot definition: per-entity sequences + curation signals |
| `dataset.lock` | ~18 MB | reproduction anchor: embedded config + all entry IDs + candidates SHA-256 |
| `manifest.json` | ~23 MB | full record incl. per-split entry lists, entry→component map, per-entry ligand tiers |
| `splits.registry.json` | ~4.7 MB | canonical_key → split, for growth-stable regeneration |

To distribute these as a citable, downloadable split, attach them to a GitHub
Release (or Zenodo) rather than committing them. `build` reproduces them exactly
from `config.yaml` (same config hash → byte-identical manifest).

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

- **Test set, functional tier:** metal 4,099 · small-molecule 3,530 · nucleotide 586
- **Test set, ambiguous (reported, not labelled):** small-molecule 3,596 · metal 73 · nucleotide 1
  - The small-molecule ambiguous count ≈ the functional one: roughly half of
    bound-looking small molecules aren't corroborated by contact or a measured
    affinity, so they're flagged rather than silently labelled.
- **404** His-tag/Ni(Co) purification artifacts flagged and demoted from the metal
  class — the LigandMPNN metal-set blemish, caught automatically.

Every structure stays in its split regardless of ligand quality; only the labels
and confidence tiers change.
