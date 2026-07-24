"""Measure the unclustered-chain merge fan-out from a candidates.jsonl — no build.

Derives the MIN_UNCLUSTERED_MERGE_LEN knee (cluster.py) from the data instead of
intuition. For every UNCLUSTERED protein chain that appears inside an otherwise-clustered
entry (a partial entry), it counts how many distinct sequence clusters that chain's
sequence would bridge if it were allowed to union — the merge fan-out, i.e. the size of
the spurious mega-component a promiscuous peptide would seed. Cross-tabulated by sequence
length so you can see the bimodal split (short peptides fan out; long chains don't) and
pick the length threshold at the knee.

    uv run python scripts/measure_unclustered_fanout.py path/to/candidates.jsonl

Run this during a full-PDB validation, then set MIN_UNCLUSTERED_MERGE_LEN to the length
above which max fan-out collapses to ~1, and record the number + date in PLAN.md.
"""

from __future__ import annotations

import sys
from pathlib import Path

from ifsplit.schema import read_candidates_jsonl

LEVEL = 30  # 30% identity — IF-Split's clustering level


def main(path: str) -> int:
    records = read_candidates_jsonl(Path(path))
    # unclustered sequence -> (length, set of distinct clustered cluster-ids it co-occurs
    # with across partial entries). |set| is the fan-out (mega-component seed size).
    fanout: dict[str, set[int]] = {}
    length: dict[str, int] = {}
    n_partial = 0
    for r in records:
        proteins = [e for e in r.polymer_entities if e.is_protein]
        clustered_ids = {e.cluster_ids[LEVEL] for e in proteins if LEVEL in e.cluster_ids}
        uncl = [e for e in proteins if LEVEL not in e.cluster_ids]
        if not clustered_ids or not uncl:
            continue  # only partial entries (some clustered + some unclustered) fan out
        n_partial += 1
        for e in uncl:
            fanout.setdefault(e.seq, set()).update(clustered_ids)
            length[e.seq] = len(e.seq)

    if not fanout:
        print("No partial entries with unclustered chains found. Nothing to gate.")
        return 0

    ranked = sorted(fanout.items(), key=lambda kv: len(kv[1]), reverse=True)
    print(f"partial entries with >=1 unclustered chain: {n_partial}")
    print(f"distinct unclustered sequences in partial entries: {len(fanout)}\n")

    print("top merge fan-out (distinct sequence clusters a shared unclustered chain bridges):")
    print("  fan-out  length  seq[:24]")
    for seq, clusters in ranked[:20]:
        print(f"  {len(clusters):>7}  {length[seq]:>6}  {seq[:24]}")

    print("\nmax fan-out among sequences with length >= threshold (find the knee):")
    print("  threshold  max_fanout  n_seqs>=thr")
    for thr in (10, 12, 15, 20, 25, 30, 35, 40, 45, 50, 60):
        elig = [len(c) for s, c in fanout.items() if length[s] >= thr]
        mx = max(elig) if elig else 0
        print(f"  {thr:>9}  {mx:>10}  {len(elig):>11}")
    print(
        "\nPick the smallest threshold at which max_fanout collapses to ~1 (no mega-component"
        "\nseeds survive the gate) and set MIN_UNCLUSTERED_MERGE_LEN there; record it in PLAN.md."
    )
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1]))
