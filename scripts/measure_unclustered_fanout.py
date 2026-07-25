"""Measure the unclustered-chain merge fan-out from a candidates.jsonl — no build.

Derives the MIN_UNCLUSTERED_MERGE_MODELED knee (cluster.py) from the data. For every
UNCLUSTERED protein chain that appears inside an otherwise-clustered entry (a partial
entry), it counts how many distinct sequence clusters that chain's sequence would bridge
if allowed to union — the merge fan-out, i.e. the size of the spurious mega-component a
promiscuous chain would seed. Cross-tabulated by MODELED (non-'X') residue count, because
the full-PDB finding is that the driver is unmodeled/low-complexity sequence, not raw
length: poly-'X' chains fan out at all lengths, so the gate must be on modeled content.

    uv run python scripts/measure_unclustered_fanout.py path/to/candidates.jsonl

Set MIN_UNCLUSTERED_MERGE_MODELED to the modeled count above which max fan-out collapses to
~1, and record the number + snapshot date in PLAN.md.
"""

from __future__ import annotations

import sys
from pathlib import Path

from ifsplit.schema import read_candidates_jsonl

LEVEL = 30  # 30% identity — IF-Split's clustering level


def _modeled(seq: str) -> int:
    return sum(1 for c in seq if c != "X")


def main(path: str) -> int:
    records = read_candidates_jsonl(Path(path))
    # unclustered sequence -> set of distinct clustered cluster-ids it co-occurs with across
    # partial entries. |set| is the fan-out (mega-component seed size).
    fanout: dict[str, set[int]] = {}
    for r in records:
        proteins = [e for e in r.polymer_entities if e.is_protein]
        clustered = {e.cluster_ids[LEVEL] for e in proteins if LEVEL in e.cluster_ids}
        uncl = [e for e in proteins if LEVEL not in e.cluster_ids]
        if not clustered or not uncl:
            continue  # only partial entries fan out
        for e in uncl:
            fanout.setdefault(e.seq, set()).update(clustered)

    if not fanout:
        print("No partial entries with unclustered chains found. Nothing to gate.")
        return 0

    polyx = [s for s in fanout if _modeled(s) == 0]
    ranked = sorted(fanout.items(), key=lambda kv: len(kv[1]), reverse=True)
    print(f"distinct unclustered sequences in partial entries: {len(fanout)}")
    print(
        f"  all-'X' (0 modeled residues): {len(polyx)}  "
        f"max fan-out {max((len(fanout[s]) for s in polyx), default=0)}"
    )

    print("\ntop merge fan-out (distinct clusters a shared unclustered chain would bridge):")
    print("  fan-out  len  modeled  seq[:20]")
    for seq, clusters in ranked[:20]:
        print(f"  {len(clusters):>7}  {len(seq):>3}  {_modeled(seq):>7}  {seq[:20]}")

    print("\nmax fan-out among sequences with MODELED (non-'X') residues >= threshold (the knee):")
    print("  modeled>=  max_fanout  n_seqs")
    for thr in (1, 5, 8, 10, 11, 12, 15, 20, 30, 50):
        elig = [len(c) for s, c in fanout.items() if _modeled(s) >= thr]
        print(f"  {thr:>9}  {max(elig) if elig else 0:>10}  {len(elig):>6}")
    print(
        "\nSet MIN_UNCLUSTERED_MERGE_MODELED to the smallest modeled count at which max_fanout"
        "\ncollapses to ~1-2; record it + the snapshot date in PLAN.md."
    )
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1]))
