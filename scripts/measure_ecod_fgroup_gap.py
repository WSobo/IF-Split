#!/usr/bin/env python3
"""Measure the ECOD empty-F-group gap (the open bug in ``schema.py``).

``structural_families_from_instances`` keys ECOD on the annotation's ``name``, which is
its F-group. Some entries carry an ECOD annotation whose ``name`` is empty (3OLT returns
``name: None`` with lineage ``"F:"``) while still carrying a real ``annotation_id``. Those
chains are silently counted ECOD-*unclassified*, so ECOD coverage is understated by an
unknown amount and a few "novel-fold" calls may be spurious.

The magnitude cannot be read off a built snapshot, because the snapshot was parsed with the
buggy rule and the empty-name annotations are already gone. So this samples the *suspect*
population live: protein entities that CATH or SCOP2 classifies but ECOD apparently does
not. Any ECOD annotation found there with an id but no name is a chain the parser dropped.

Network: read-only RCSB Data API, one batched GraphQL query per 50 entries.

Usage:
    python scripts/measure_ecod_fgroup_gap.py CANDIDATES.jsonl [--sample 400] [--seed 0]
"""

from __future__ import annotations

import argparse
import json
import random
from pathlib import Path

from ifsplit.rcsb import RcsbClient


def suspect_entries(path: Path) -> tuple[list[str], set[str], int, int]:
    """Entry ids holding >=1 *suspect* protein entity, plus the suspect entity ids.

    A suspect entity is one the build recorded as ECOD-unclassified while CATH or SCOP2
    did classify it: ECOD covers most of what they cover, so its absence there is the
    signature of a dropped annotation rather than of genuine non-coverage. Keeping the
    entity ids matters because the projection denominator is the suspect *entities*, not
    every entity that happens to share an entry with one.

    Also returns (protein entities seen, entities ECOD classified).
    """
    suspects: list[str] = []
    suspect_entities: set[str] = set()
    n_prot = n_ecod = 0
    with path.open() as f:
        for line in f:
            r = json.loads(line)
            hit = False
            for e in r.get("polymer_entities", []):
                if e.get("polymer_type") != "Protein":
                    continue
                n_prot += 1
                fams = e.get("structural_families") or {}
                has_ecod = bool(fams.get("ecod"))
                n_ecod += has_ecod
                if not has_ecod and (fams.get("cath") or fams.get("scop2")):
                    hit = True
                    suspect_entities.add(e["entity_id"])
            if hit:
                suspects.append(r["entry_id"])
    return suspects, suspect_entities, n_prot, n_ecod


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("candidates", type=Path)
    ap.add_argument("--sample", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    suspects, suspect_ents, n_prot, n_ecod = suspect_entries(args.candidates)
    print(f"protein entities              : {n_prot:,}")
    print(f"  ECOD-classified (as built)  : {n_ecod:,} ({n_ecod / n_prot:.1%})")
    print(f"  suspect (CATH/SCOP2, no ECOD): {len(suspect_ents):,}")
    print(f"entries holding a suspect chain: {len(suspects):,}")

    rng = random.Random(args.seed)
    sample = sorted(rng.sample(suspects, min(args.sample, len(suspects))))
    print(f"sampling {len(sample)} of them live...\n")

    # Denominator is the SUSPECT entities only. Other entities in the same entry are
    # not the population being projected and would dilute the rate.
    checked = recoverable = truly_absent = 0
    examples: list[tuple[str, str]] = []
    client = RcsbClient()
    for entry in client.fetch_entries(sample):
        for pe in entry.get("polymer_entities") or []:
            if pe.get("rcsb_id") not in suspect_ents:
                continue
            checked += 1
            anns = [
                a
                for inst in (pe.get("polymer_entity_instances") or [])
                for a in (inst.get("rcsb_polymer_instance_annotation") or [])
                if a.get("type") == "ECOD"
            ]
            named = [a for a in anns if a.get("name")]
            ided = [a for a in anns if a.get("annotation_id")]
            if not named and ided:
                recoverable += 1
                if len(examples) < 5:
                    examples.append((pe["rcsb_id"], ided[0]["annotation_id"]))
            elif not anns:
                truly_absent += 1

    print(f"suspect entities re-checked live          : {checked}")
    print(f"  ECOD id present but EMPTY F-group name  : {recoverable}  <- parser drops these")
    print(f"  no ECOD annotation at all (genuine gap) : {truly_absent}")
    if examples:
        print("\nexamples of dropped (entity, annotation_id):")
        for eid, aid in examples:
            print(f"  {eid}  {aid}")

    if checked:
        rate = recoverable / checked
        recovered = rate * len(suspect_ents)
        print(f"\nrecoverable share of the suspect population: {rate:.1%}")
        print(f"projected entities recovered              : {recovered:,.0f}")
        print(f"ECOD entity coverage {n_ecod / n_prot:.1%} -> {(n_ecod + recovered) / n_prot:.1%}")
        print(
            "\nNB this projects over the SUSPECT population only (entities CATH/SCOP2\n"
            "classify but ECOD apparently does not). Entities no authority classifies are\n"
            "not sampled and are far less likely to hide an ECOD annotation, so treat the\n"
            "figure as the recoverable share of the suspect set, not of all ECOD-unclassified\n"
            "entities."
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
