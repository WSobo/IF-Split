#!/usr/bin/env python3
"""Measure a split's INTERNAL redundancy: entries vs distinct proteins vs components.

A test set is a set of entries but not a set of independent measurements. The PDB deposits
the same protein many times (a ligand-soaking series is one scaffold with n bound ligands),
so a per-entry mean over-weights whatever was crystallised most often. This reports the
three counts that bracket the honest denominator:

    entries  >=  distinct protein sequence sets  >=  leakage-safe components

Works on an IF-Split output directory and on an EXTERNAL split given as id lists, so the
same measurement can be applied to a published benchmark. Metadata only; no coordinates.

Usage:
    # IF-Split build (reads train/val/test.json next to the manifest)
    python scripts/measure_split_redundancy.py --split-dir data/rs-recommended \\
        --candidates data/run-2026.07.22/candidates-annotated.jsonl

    # External split, e.g. LigandMPNN's three ligand-class test files
    python scripts/measure_split_redundancy.py --ids lmpnn_test=a.json,b.json,c.json \\
        --candidates data/run-2026.07.22/candidates-annotated.jsonl
"""

from __future__ import annotations

import argparse
import collections
import json
from pathlib import Path


def _load_ids(path: Path) -> list[str]:
    data = json.loads(path.read_text())
    if isinstance(data, list):
        return [str(x) for x in data]
    return [str(x) for x in data.get("entries", data.get("ids", []))]


def _sequence_key(rec: dict) -> tuple[str, ...]:
    """Sorted tuple of an entry's modeled protein sequences (its identity as a backbone)."""
    return tuple(
        sorted(
            e["seq"]
            for e in rec.get("polymer_entities", [])
            if e.get("polymer_type") == "Protein"
            and e.get("seq")
            and any(c != "X" for c in e["seq"])
        )
    )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--candidates", type=Path, required=True)
    ap.add_argument("--split-dir", type=Path, help="IF-Split output dir (train/val/test.json)")
    ap.add_argument(
        "--ids",
        action="append",
        default=[],
        help="NAME=file1.json[,file2.json...] — an external split's entry-id lists",
    )
    args = ap.parse_args()

    sets: dict[str, set[str]] = {}
    if args.split_dir:
        for name in ("train", "val", "test"):
            p = args.split_dir / f"{name}.json"
            if p.exists():
                sets[name] = {i.upper() for i in _load_ids(p)}
    for spec in args.ids:
        name, _, files = spec.partition("=")
        ids: set[str] = set()
        for f in files.split(","):
            ids |= {i.upper() for i in _load_ids(Path(f))}
        sets[name] = ids
    if not sets:
        ap.error("give --split-dir and/or --ids")

    wanted = set().union(*sets.values())
    seq_of: dict[str, tuple[str, ...]] = {}
    comp_of: dict[str, str] = {}
    with args.candidates.open() as f:
        for line in f:
            rec = json.loads(line)
            eid = rec["entry_id"].upper()
            if eid in wanted:
                seq_of[eid] = _sequence_key(rec)
                comp_of[eid] = rec.get("component") or ""

    # Component ids are not carried in candidates.jsonl; read them from the build if present.
    if args.split_dir and (args.split_dir / "clusters.json").exists():
        doc = json.loads((args.split_dir / "clusters.json").read_text())
        ec = {k.upper(): v for k, v in doc.get("entry_clusters", {}).items()}
        comp_of = {e: ec.get(e, "") for e in wanted}

    print(f"{'split':14s} {'entries':>8s} {'distinct':>9s} {'redundant':>10s} {'components':>11s}")
    for name, ids in sets.items():
        resolved = [e for e in ids if e in seq_of]
        groups = collections.Counter(seq_of[e] for e in resolved)
        comps = {comp_of.get(e, "") for e in resolved if comp_of.get(e)}
        n, d = len(resolved), len(groups)
        miss = len(ids) - n
        print(
            f"{name:14s} {n:8d} {d:9d} {(n - d) / n:10.1%} "
            f"{len(comps) if comps else '-':>11}" + (f"   ({miss} unresolvable)" if miss else "")
        )
        big = groups.most_common(1)
        if big and big[0][1] > 1:
            rep = sorted(e for e in resolved if seq_of[e] == big[0][0])
            print(f"{'':14s} largest identical group: {big[0][1]} entries, e.g. {rep[0]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
