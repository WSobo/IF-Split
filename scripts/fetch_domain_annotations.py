#!/usr/bin/env python3
"""Fetch Pfam/InterPro domain annotations for every protein entity in a candidates.jsonl.

Why this exists
---------------
CATH/ECOD/SCOP2 are curated *retrospectively* by the classification groups, so their
coverage lags deposition badly: measured on the 2026-07-22 snapshot, 0.8% of pre-2020
entries are unclassified versus 36.5% of 2025 and 55.3% of 2026 releases. A "fold-disjoint"
holdout built only from those three authorities is therefore disjoint largely because the
databases have not reached the newest entries yet -- not because the folds are novel.

Pfam/InterPro are HMM models applied to *sequence*, so they carry no curation lag and they
detect remote homology well below the 30% identity a sequence clusterer can see. They cover
61.3% of the entries CATH/ECOD/SCOP2 cannot classify, and 88.3% of those turn out to share a
family with training -- i.e. they expose leakage the structural authorities cannot.

This script is the *prototype* path: it annotates an existing candidates.jsonl in place-ish
(writing a sidecar cache) so the merge effect can be sized without re-enumerating the whole
snapshot from Stage 1. It is metadata-only and downloads no coordinates.

Usage
-----
    python scripts/fetch_domain_annotations.py CANDIDATES.jsonl OUT_CACHE.jsonl
    python scripts/fetch_domain_annotations.py CANDIDATES.jsonl CACHE.jsonl --apply OUT.jsonl

Resumable: re-running skips entities already present in OUT_CACHE.jsonl. With ``--apply`` the
cache is folded into a copy of CANDIDATES as each entity's ``domain_families``, which is the
migration path for an existing snapshot -- adopting ``structural_clustering: pfam|interpro|all``
then needs no Stage-1 re-enumeration.
"""

from __future__ import annotations

import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

GRAPHQL = "https://data.rcsb.org/graphql"
BATCH = 1000
SLEEP = 0.15  # be a polite client; RCSB is a shared public resource
WANTED = ("Pfam", "InterPro")

QUERY = (
    "query($ids:[String!]!){polymer_entities(entity_ids:$ids){"
    "rcsb_id rcsb_polymer_entity_annotation{type annotation_id}}}"
)


def protein_entity_ids(candidates: Path) -> list[str]:
    ids: list[str] = []
    with candidates.open() as fh:
        for line in fh:
            rec = json.loads(line)
            for pe in rec.get("polymer_entities", []):
                if pe.get("polymer_type") == "Protein" and pe.get("seq"):
                    ids.append(pe["entity_id"])
    return ids


def fetch(batch: list[str], attempts: int = 4) -> list[dict]:
    payload = json.dumps({"query": QUERY, "variables": {"ids": batch}}).encode()
    for attempt in range(attempts):
        req = urllib.request.Request(
            GRAPHQL, data=payload, headers={"Content-Type": "application/json"}
        )
        try:
            with urllib.request.urlopen(req, timeout=180) as resp:
                data = json.load(resp)
            return [e for e in (data.get("data", {}).get("polymer_entities") or []) if e]
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            if attempt == attempts - 1:
                print(f"  !! giving up on a batch of {len(batch)}: {exc}", flush=True)
                return []
            time.sleep(2 * (attempt + 1))
    return []


def apply_cache(candidates: Path, cache: Path, dest: Path) -> None:
    """Write a copy of ``candidates`` with each entity's ``domain_families`` filled in."""
    dom: dict[str, dict[str, list[str]]] = {}
    with cache.open() as fh:
        for line in fh:
            rec = json.loads(line)
            fams = {k: rec[k] for k in ("pfam", "interpro") if rec.get(k)}
            if fams:
                dom[rec["entity_id"]] = fams
    n_ent = 0
    with candidates.open() as src, dest.open("w") as out:
        for line in src:
            rec = json.loads(line)
            for pe in rec.get("polymer_entities", []):
                fams = dom.get(pe.get("entity_id"))
                if fams:
                    pe["domain_families"] = fams
                    n_ent += 1
            out.write(json.dumps(rec) + "\n")
    print(f"applied domain families to {n_ent} entities -> {dest}", flush=True)


def main() -> int:
    if len(sys.argv) == 5 and sys.argv[3] == "--apply":
        apply_cache(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[4]))
        return 0
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    candidates, out = Path(sys.argv[1]), Path(sys.argv[2])
    done: set[str] = set()
    if out.exists():
        with out.open() as fh:
            for line in fh:
                try:
                    done.add(json.loads(line)["entity_id"])
                except (json.JSONDecodeError, KeyError):
                    continue
    ids = [i for i in protein_entity_ids(candidates) if i not in done]
    print(f"protein entities: {len(ids)} to fetch ({len(done)} already cached)", flush=True)

    started = time.time()
    with out.open("a") as fh:
        for start in range(0, len(ids), BATCH):
            batch = ids[start : start + BATCH]
            for ent in fetch(batch):
                ann = ent.get("rcsb_polymer_entity_annotation") or []
                rec = {"entity_id": ent["rcsb_id"]}
                for kind in WANTED:
                    hits = sorted({a["annotation_id"] for a in ann if a.get("type") == kind})
                    if hits:
                        rec[kind.lower()] = hits
                fh.write(json.dumps(rec) + "\n")
            fh.flush()
            seen = start + len(batch)
            rate = seen / max(time.time() - started, 1e-6)
            eta = (len(ids) - seen) / rate if rate else 0
            print(f"  {seen}/{len(ids)}  ({rate:.0f} ent/s, ETA {eta / 60:.1f} min)", flush=True)
            time.sleep(SLEEP)
    print("done", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
