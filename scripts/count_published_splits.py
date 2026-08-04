"""Recount the published ProteinMPNN / LigandMPNN splits from their release files.

Both papers state their test-set sizes and snapshot dates but not their train,
validation or cluster totals, so every such figure quoted in the audit is counted
here rather than cited. This script downloads the release files and prints those
counts, so the numbers can be checked without trusting the write-up.

It also reports the integrity checks the audit relies on: how many entry ids the
five LigandMPNN lists contain versus how many are distinct, and which ids appear
in more than one list. Two of them (``2zio``, ``3olt``) are in both ``train.json``
and ``test_nucleotide.json``, i.e. the same entry is trained on and tested on.

``--check-composition`` asks RCSB whether each ligand-class test entry contains
anything of the class it is filed under. Two do not: ``2zio`` (pyrrolysyl-tRNA
synthetase with a Lys-AMP analog, no tRNA in the crystal) and ``3olt`` (COX-2 with
arachidonic acid, nothing nucleotide about it) are in ``test_nucleotide.json``
with no nucleic acid between them, while the other 72 all have one. They are the
same two entries that also appear in ``train.json``, which is why they are worth
naming: three unrelated checks single out one pair.

The metal and small-molecule sets pass the same check, so this is a defect in one
list rather than a fault in how all three were assembled. Note the check is weak by
design -- it asks only whether the class is present, not whether it is functional,
which is the ligand tiering's job and a question 12 further entries fail.

Sources (all public, no credentials):
  LigandMPNN  github.com/dauparas/LigandMPNN/training/*.json
  ProteinMPNN files.ipd.uw.edu/pub/training_sets/pdb_2021aug02_sample.tar.gz
                (~47 MB; only list.csv + {valid,test}_clusters.txt are read, and
                 they are the full split, not a sample)
  RCSB        data.rcsb.org/graphql (only with --check-composition; metadata only)

Usage:
  uv run python scripts/count_published_splits.py
  uv run python scripts/count_published_splits.py --lmpnn-dir path/to/jsons \
      --pmpnn-tar path/to/pdb_2021aug02_sample.tar.gz   # offline, from local copies
  uv run python scripts/count_published_splits.py --skip-pmpnn --check-composition
  uv run python scripts/count_published_splits.py --json counts.json          # machine-readable
"""

from __future__ import annotations

import argparse
import csv
import io
import itertools
import json
import re
import tarfile
from collections import Counter
from pathlib import Path
from typing import Any

import httpx

from ifsplit.ligands import METAL_ELEMENTS

LMPNN_RAW = "https://raw.githubusercontent.com/dauparas/LigandMPNN/main/training"
LMPNN_API = "https://api.github.com/repos/dauparas/LigandMPNN/contents/training"
# Fallback only. The directory is listed at run time so that a list added upstream
# is scanned too, rather than silently skipped by a hardcoded set of five names.
LMPNN_FILES_FALLBACK = [
    "train.json",
    "valid.json",
    "test_small_molecule.json",
    "test_nucleotide.json",
    "test_metal.json",
]
PMPNN_TAR = "https://files.ipd.uw.edu/pub/training_sets/pdb_2021aug02_sample.tar.gz"
# Which lists are training rather than held out. Used only to label an overlap: a
# train/held-out pair is a leak, a held-out/held-out pair is double-counting.
TRAIN_LISTS = {"train.json", "valid.json"}


def _list_lmpnn_files(client: httpx.Client) -> list[str]:
    try:
        r = client.get(LMPNN_API)
        r.raise_for_status()
        names = sorted(e["name"] for e in r.json() if e["type"] == "file")
    except (httpx.HTTPError, KeyError, ValueError):
        return LMPNN_FILES_FALLBACK
    found = [n for n in names if n.endswith(".json")]
    return found or LMPNN_FILES_FALLBACK


def _load_lmpnn(local: Path | None) -> dict[str, list[str]]:
    if local is not None:
        found = sorted(p.name for p in local.glob("*.json"))
        if not found:
            raise SystemExit(f"no .json files in {local}")
        return {f: json.loads((local / f).read_text()) for f in found}
    out: dict[str, list[str]] = {}
    with httpx.Client(timeout=120, follow_redirects=True) as client:
        for fname in _list_lmpnn_files(client):
            r = client.get(f"{LMPNN_RAW}/{fname}")
            r.raise_for_status()
            payload = r.json()
            if isinstance(payload, list):
                out[fname] = payload
    return out


def count_lmpnn(lists: dict[str, list[str]]) -> dict[str, Any]:
    sizes = {f: len(v) for f, v in lists.items()}
    all_ids = [i for v in lists.values() for i in v]
    counts = Counter(all_ids)
    as_sets = {f: set(v) for f, v in lists.items()}

    # Every pair, not just the ones we expected to collide.
    pairwise = []
    for a, b in itertools.combinations(sorted(as_sets), 2):
        shared = sorted(as_sets[a] & as_sets[b])
        leak = (a in TRAIN_LISTS) != (b in TRAIN_LISTS)
        pairwise.append({"a": a, "b": b, "n": len(shared), "shared": shared, "leak": leak})

    overlaps = []
    for entry_id in sorted(k for k, n in counts.items() if n > 1):
        where = sorted(f for f, v in as_sets.items() if entry_id in v)
        in_train = [f for f in where if f in TRAIN_LISTS]
        in_heldout = [f for f in where if f not in TRAIN_LISTS]
        overlaps.append(
            {
                "entry_id": entry_id,
                "lists": where,
                # trained on AND tested on; the held-out/held-out case is harmless
                "train_test_straddle": bool(in_train) and bool(in_heldout),
            }
        )

    return {
        "sizes": sizes,
        "total_ids": len(all_ids),
        "distinct_ids": len(counts),
        "pairwise": pairwise,
        "overlaps": overlaps,
        "straddles": [o["entry_id"] for o in overlaps if o["train_test_straddle"]],
        # ids are lowercase in the released files; the audit upper-cases before
        # querying RCSB, which is why a case-sensitive check would miss these
        "ids_are_lowercase": all(i == i.lower() for i in all_ids),
    }


NUCLEIC_TYPES = {"DNA", "RNA", "NA-hybrid"}
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
# What each ligand-class test set has to contain for its own conditioning signal to
# exist at all. The requirement is deliberately weak: not "is this a good example of
# the class" (that is the tiering's job) but "is anything of the class present".
CLASS_REQUIREMENT = {
    "test_nucleotide.json": "a DNA, RNA or hybrid chain",
    "test_metal.json": "a metal-bearing component",
    "test_small_molecule.json": "any non-polymer component",
}


def _elements(formula: str) -> set[str]:
    return {m.upper() for m in re.findall(r"[A-Z][a-z]?", formula or "")}


def _satisfies(list_name: str, entry: dict[str, Any]) -> bool:
    polymers = {
        (pe.get("entity_poly") or {}).get("rcsb_entity_polymer_type")
        for pe in (entry.get("polymer_entities") or [])
    }
    comps = [ne["nonpolymer_comp"]["chem_comp"] for ne in (entry.get("nonpolymer_entities") or [])]
    if list_name == "test_nucleotide.json":
        return bool(polymers & NUCLEIC_TYPES)
    if list_name == "test_metal.json":
        # A metal counts whether it is a bare ion or carried by a cofactor such as
        # heme, so match on the formula rather than on a list of ion comp ids.
        return any(_elements(c.get("formula") or "") & METAL_ELEMENTS for c in comps)
    return bool(comps)


def check_class_composition(lists: dict[str, list[str]]) -> dict[str, Any]:
    """Does each ligand-class test entry contain anything of the class it is filed under?

    LigandMPNN is scored on each set for residues near that class of context, so an
    entry with none of it contributes nothing to its own metric no matter what else it
    holds. This is a composition check, not a quality one: an entry can pass here and
    still be a purification artifact, which is what the ligand tiering catches.
    """
    out: dict[str, Any] = {}
    for list_name, requirement in CLASS_REQUIREMENT.items():
        ids = sorted({i.upper() for i in lists.get(list_name, [])})
        if not ids:
            continue
        query = (
            f"{{ entries(entry_ids: {json.dumps(ids)}) {{ rcsb_id struct {{ title }} "
            "polymer_entities { entity_poly { rcsb_entity_polymer_type } } "
            "nonpolymer_entities { nonpolymer_comp { chem_comp { id name formula } } } } }"
        )
        with httpx.Client(timeout=180, follow_redirects=True) as client:
            r = client.post(RCSB_GRAPHQL, json={"query": query})
            r.raise_for_status()
            entries = r.json()["data"]["entries"]

        failing = [
            {
                "entry_id": e["rcsb_id"],
                "title": (e.get("struct") or {}).get("title", ""),
                "polymer_types": sorted(
                    t
                    for t in {
                        (pe.get("entity_poly") or {}).get("rcsb_entity_polymer_type")
                        for pe in (e.get("polymer_entities") or [])
                    }
                    if t
                ),
                "ligands": [
                    (c["id"], c["name"])
                    for c in (
                        n["nonpolymer_comp"]["chem_comp"]
                        for n in (e.get("nonpolymer_entities") or [])
                    )
                ],
            }
            for e in entries
            if not _satisfies(list_name, e)
        ]
        out[list_name] = {
            "requirement": requirement,
            "listed": len(ids),
            "resolved": len(entries),
            "unresolved": sorted(set(ids) - {e["rcsb_id"] for e in entries}),
            "failing": failing,
        }
    return out


def _load_pmpnn(tar_path: Path | None) -> dict[str, bytes]:
    if tar_path is not None:
        blob = tar_path.read_bytes()
    else:
        with httpx.Client(timeout=600, follow_redirects=True) as client:
            r = client.get(PMPNN_TAR)
            r.raise_for_status()
            blob = r.content
    wanted = ("list.csv", "valid_clusters.txt", "test_clusters.txt")
    out: dict[str, bytes] = {}
    with tarfile.open(fileobj=io.BytesIO(blob), mode="r:gz") as tf:
        for member in tf:
            name = Path(member.name).name
            if name in wanted and member.isfile():
                fh = tf.extractfile(member)
                if fh is not None:
                    out[name] = fh.read()
    missing = [w for w in wanted if w not in out]
    if missing:
        raise SystemExit(f"archive is missing {missing}; expected the pdb_2021aug02 layout")
    return out


def count_pmpnn(tar_path: Path | None) -> dict[str, Any]:
    blobs = _load_pmpnn(tar_path)
    rows = list(csv.DictReader(io.StringIO(blobs["list.csv"].decode())))
    valid = {ln.strip() for ln in blobs["valid_clusters.txt"].decode().splitlines() if ln.strip()}
    test = {ln.strip() for ln in blobs["test_clusters.txt"].decode().splitlines() if ln.strip()}

    # Roll the chain-level split up to entries the way the audit does: an entry is
    # held out if ANY of its chains is, test winning over validation.
    entry_split: dict[str, str] = {}
    for r in rows:
        entry_id = r["CHAINID"].split("_")[0]
        cluster = r["CLUSTER"]
        split = "test" if cluster in test else "valid" if cluster in valid else "train"
        prior = entry_split.get(entry_id)
        if prior is None or split == "test" or (split == "valid" and prior == "train"):
            entry_split[entry_id] = split

    # An entry with chains on both sides of the split is the same deposition in train
    # and in a held-out set at once — possible here only because the split is per chain.
    per_entry_splits: dict[str, set[str]] = {}
    for r in rows:
        entry_id = r["CHAINID"].split("_")[0]
        cluster = r["CLUSTER"]
        split = "test" if cluster in test else "valid" if cluster in valid else "train"
        per_entry_splits.setdefault(entry_id, set()).add(split)
    straddling = sorted(e for e, s in per_entry_splits.items() if {"train", "test"} <= s)
    straddling_valid = sorted(e for e, s in per_entry_splits.items() if {"train", "valid"} <= s)

    return {
        "chains": len(rows),
        "distinct_chain_ids": len({r["CHAINID"] for r in rows}),
        "clusters": len({r["CLUSTER"] for r in rows}),
        "entries": len(entry_split),
        "valid_clusters": len(valid),
        "test_clusters": len(test),
        "entry_level": dict(Counter(entry_split.values())),
        "train_test_straddling_entries": straddling,
        "train_valid_straddling_entries": straddling_valid,
    }


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--lmpnn-dir", type=Path, default=None, help="local dir of split JSONs")
    ap.add_argument("--pmpnn-tar", type=Path, default=None, help="local pdb_2021aug02 tarball")
    ap.add_argument("--skip-pmpnn", action="store_true", help="skip the 47 MB download")
    ap.add_argument(
        "--check-composition",
        action="store_true",
        help="ask RCSB whether each ligand-class test entry contains anything of its class",
    )
    ap.add_argument("--json", type=Path, default=None, help="also write the counts here")
    args = ap.parse_args()

    report: dict[str, Any] = {}

    lists = _load_lmpnn(args.lmpnn_dir)
    lm = count_lmpnn(lists)
    report["ligandmpnn"] = lm
    print(f"LigandMPNN ({LMPNN_RAW.split('//')[1]}) -- {len(lm['sizes'])} lists")
    for fname in sorted(lm["sizes"]):
        print(f"  {fname:<28} {lm['sizes'][fname]:>7,}")
    print(f"  {'total ids':<28} {lm['total_ids']:>7,}")
    print(f"  {'distinct ids':<28} {lm['distinct_ids']:>7,}")

    print("\n  every pair of lists:")
    for p in lm["pairwise"]:
        if p["n"] == 0:
            verdict = "clean"
        elif p["leak"]:
            verdict = "TRAINED ON AND TESTED ON: " + ", ".join(p["shared"])
        else:
            verdict = "double-counted: " + ", ".join(p["shared"])
        print(f"    {p['a']:<26} n {p['b']:<26} {p['n']:>3}  {verdict}")

    straddles = lm["straddles"]
    print(
        f"\n  {len(straddles)} id(s) in both a training and a held-out list"
        + (f": {', '.join(straddles)}" if straddles else "")
    )
    if not lm["ids_are_lowercase"]:
        print("  note: ids are not uniformly lowercase in this copy")

    if args.check_composition:
        comp = check_class_composition(lists)
        report["class_composition"] = comp
        print("\n  does each ligand-class test entry contain anything of its class?")
        for list_name, res in comp.items():
            failing = res["failing"]
            verdict = "clean" if not failing else f"{len(failing)} FAIL"
            print(
                f"    {list_name:<26} needs {res['requirement']:<27}"
                f" {res['resolved']:>3} checked  {verdict}"
            )
            for m in failing:
                print(f"      {m['entry_id']}  polymers={m['polymer_types']}  {m['title'][:62]}")
                for cid, name in m["ligands"]:
                    print(f"          {cid}: {name[:66]}")
            if res["unresolved"]:
                print(f"      unresolved by RCSB: {', '.join(res['unresolved'])}")
        print(
            "    NB 'clean' here means present, not functional. It is the weaker question:\n"
            "    the metal set passes this entry for entry and is still the most contaminated\n"
            "    of the three (6 of 83 adventitious). Run scripts/audit_ligandmpnn_split.py\n"
            "    for the tiering that asks whether a site is real."
        )

    if not args.skip_pmpnn:
        pm = count_pmpnn(args.pmpnn_tar)
        report["proteinmpnn"] = pm
        print("\nProteinMPNN (files.ipd.uw.edu/pub/training_sets, snapshot 2021-08-02)")
        print(f"  {'chains (list.csv rows)':<28} {pm['chains']:>7,}")
        print(f"  {'sequence clusters':<28} {pm['clusters']:>7,}")
        print(f"  {'distinct entries':<28} {pm['entries']:>7,}")
        n_v, n_t = pm["valid_clusters"], pm["test_clusters"]
        print(f"  {'valid / test clusters':<28} {n_v:>7,} / {n_t:,}")
        el = pm["entry_level"]
        print(
            f"  rolled up to entries:        train {el.get('train', 0):,} / "
            f"valid {el.get('valid', 0):,} / test {el.get('test', 0):,}"
        )
        tt = pm["train_test_straddling_entries"]
        tv = pm["train_valid_straddling_entries"]
        print(f"  entries with chains in train AND a held-out set: {len(tt) + len(tv)}")
        print(f"    train + test  ({len(tt)}): {', '.join(tt) if tt else '-'}")
        print(f"    train + valid ({len(tv)}): {', '.join(tv) if tv else '-'}")

    if args.json:
        args.json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
