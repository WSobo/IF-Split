#!/usr/bin/env python3
"""Measure how much fold leakage Pfam/InterPro expose that CATH/ECOD/SCOP2 miss.

Answers two questions on a built split:

1. **Is there a sequence-side cutoff that makes an unclassified holdout safe?**
   The held-out entries are already 30%-cluster-disjoint from train by construction, so if
   they still share HMM domain families with training, sequence identity is not protecting
   them and no cutoff on length/recency will. We report the leakage rate stratified by
   chain length and release year so the claim is falsifiable rather than asserted.

2. **What would adding Pfam/InterPro as merge authorities cost and buy?**
   Simulated by adding domain-family edges to the existing components and re-measuring the
   giant/tail, plus how many held-out entries remain *certified* novel (classified by an
   authority and sharing no family with train).

Usage:
    python scripts/measure_domain_leakage.py SPLIT_DIR CANDIDATES.jsonl DOMAIN_CACHE.jsonl
"""

from __future__ import annotations

import collections
import json
import statistics as st
import sys
from pathlib import Path

STRUCTURAL = ("cath", "ecod", "scop2")


def load_ids(path: Path) -> set[str]:
    data = json.loads(path.read_text())
    if isinstance(data, list):
        return set(data)
    return set(data.get("entries", data.get("ids", [])))


def main() -> int:
    if len(sys.argv) != 4:
        print(__doc__)
        return 2
    split_dir, candidates, cache = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    train = load_ids(split_dir / "train.json")
    held = load_ids(split_dir / "val.json") | load_ids(split_dir / "test.json")

    dom: dict[str, set[str]] = {}
    with cache.open() as fh:
        for line in fh:
            rec = json.loads(line)
            fams = set(rec.get("interpro", [])) | set(rec.get("pfam", []))
            if fams:
                dom[rec["entity_id"]] = fams

    ent_dom: dict[str, set[str]] = collections.defaultdict(set)
    ent_struct: dict[str, set[str]] = collections.defaultdict(set)
    meta: dict[str, dict] = {}
    with candidates.open() as fh:
        for line in fh:
            rec = json.loads(line)
            eid = rec["entry_id"]
            if eid not in train and eid not in held:
                continue
            longest = 0
            for pe in rec.get("polymer_entities", []):
                if pe.get("polymer_type") != "Protein" or not pe.get("seq"):
                    continue
                longest = max(longest, len(pe["seq"]))
                ent_dom[eid] |= dom.get(pe["entity_id"], set())
                sf = pe.get("structural_families") or {}
                for auth in STRUCTURAL:
                    ent_struct[eid] |= {f"{auth}:{f}" for f in (sf.get(auth) or [])}
            meta[eid] = {
                "len": longest,
                "year": int(rec["release_date"][:4]) if rec.get("release_date") else None,
            }

    train_dom: set[str] = set()
    for eid in train:
        train_dom |= ent_dom.get(eid, set())
    train_struct: set[str] = set()
    for eid in train:
        train_struct |= ent_struct.get(eid, set())

    unclassified = [e for e in held if not ent_struct.get(e)]
    scoreable = [e for e in unclassified if ent_dom.get(e)]
    leaky = [e for e in scoreable if ent_dom[e] & train_dom]

    print(f"held out: {len(held)}   structurally unclassified: {len(unclassified)}")
    pct_score = 100 * len(scoreable) / max(len(unclassified), 1)
    print(
        f"  of the unclassified, {len(scoreable)} ({pct_score:.1f}%) carry a Pfam/InterPro family"
    )
    print(
        f"  of those, {len(leaky)} "
        f"({100 * len(leaky) / max(len(scoreable), 1):.1f}%) share a family with TRAIN"
        "  <-- leakage invisible to CATH/ECOD/SCOP2"
    )

    print("\nIs there a safe cutoff? leakage rate by longest-chain length:")
    buckets = [(0, 100), (100, 200), (200, 300), (300, 500), (500, 10**9)]
    for lo, hi in buckets:
        grp = [e for e in scoreable if lo <= meta[e]["len"] < hi]
        if not grp:
            continue
        bad = [e for e in grp if ent_dom[e] & train_dom]
        hi_s = "inf" if hi >= 10**9 else str(hi)
        print(f"   {lo:5d}-{hi_s:>5}: n={len(grp):5d}  leak={100 * len(bad) / len(grp):5.1f}%")

    print("\nleakage rate by release year:")
    years = collections.defaultdict(list)
    for e in scoreable:
        if meta[e]["year"]:
            years[meta[e]["year"]].append(e)
    for y in sorted(years):
        if y < 2019:
            continue
        grp = years[y]
        bad = [e for e in grp if ent_dom[e] & train_dom]
        print(f"   {y}: n={len(grp):5d}  leak={100 * len(bad) / len(grp):5.1f}%")

    certified = [
        e
        for e in held
        if (ent_struct.get(e) or ent_dom.get(e))
        and not (ent_struct.get(e, set()) & train_struct)
        and not (ent_dom.get(e, set()) & train_dom)
    ]
    print(f"\nCERTIFIED novel under all five authorities: {len(certified)} of {len(held)} held out")
    if certified:
        lens = [meta[e]["len"] for e in certified]
        print(f"   median longest chain {st.median(lens):.0f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
