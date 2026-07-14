"""Measure the fold-level leakage control (Stage 5 structural clustering).

For each classification method (cath / ecod / scop2) over a built snapshot's
``candidates.jsonl``, reports:

  - coverage: fraction of protein entities RCSB actually classifies (the ceiling
    on how much fold-leakage the method can catch),
  - merge effect: how many sequence-only components the method folds together via
    shared (super)families, and how many distinct families did the bridging.

The merge counts need the whole snapshot (fold-sharing is cross-cluster), so this
reads a candidates.jsonl rather than a stride sample. Run after any structural
change to compare methods on real data:

  uv run python scripts/eval_structural_clustering.py path/to/candidates.jsonl
"""

from __future__ import annotations

import sys

from ifsplit.cluster import build_clusters
from ifsplit.config import load_config
from ifsplit.parse import filter_candidates
from ifsplit.schema import STRUCTURAL_METHODS, read_candidates_jsonl


def main() -> None:
    path = sys.argv[1] if len(sys.argv) > 1 else "data/out/candidates.jsonl"
    cfg = load_config("config/default.yaml")
    records = read_candidates_jsonl(path)
    kept, _ = filter_candidates(records, cfg)
    proteins = [e for r in kept for e in r.polymer_entities if e.is_protein]
    n_prot = len(proteins)
    print(f"{path}: {len(kept)} kept entries, {n_prot} protein entities\n")

    baseline = build_clusters(kept, cfg.model_copy(update={"structural_clustering": "off"}))
    print(f"sequence-only components: {baseline.n_clusters}\n")
    print(f"{'method':7s} {'coverage':>18s} {'components':>12s} {'folded':>8s} {'families':>9s}")
    for m in STRUCTURAL_METHODS:
        covered = sum(1 for e in proteins if e.structural_families.get(m))
        cr = build_clusters(kept, cfg.model_copy(update={"structural_clustering": m}))
        folded = baseline.n_clusters - cr.n_clusters
        pct = 100 * covered / n_prot if n_prot else 0
        print(
            f"{m:7s} {covered:>8d} ({pct:4.1f}%) {cr.n_clusters:>12d} "
            f"{folded:>8d} {cr.n_structural_families:>9d}"
        )


if __name__ == "__main__":
    main()
