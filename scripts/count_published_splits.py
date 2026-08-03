"""Recount the published ProteinMPNN / LigandMPNN splits from their release files.

Both papers state their test-set sizes and snapshot dates but not their train,
validation or cluster totals, so every such figure quoted in the audit is counted
here rather than cited. This script downloads the release files and prints those
counts, so the numbers can be checked without trusting the write-up.

It also reports the integrity checks the audit relies on: how many entry ids the
five LigandMPNN lists contain versus how many are distinct, and which ids appear
in more than one list. Two of them (``2zio``, ``3olt``) are in both ``train.json``
and ``test_nucleotide.json``, i.e. the same entry is trained on and tested on.

Sources (both public, no credentials):
  LigandMPNN  github.com/dauparas/LigandMPNN/training/*.json
  ProteinMPNN files.ipd.uw.edu/pub/training_sets/pdb_2021aug02_sample.tar.gz
                (~47 MB; only list.csv + {valid,test}_clusters.txt are read, and
                 they are the full split, not a sample)

Usage:
  uv run python scripts/count_published_splits.py
  uv run python scripts/count_published_splits.py --lmpnn-dir path/to/jsons \
      --pmpnn-tar path/to/pdb_2021aug02_sample.tar.gz   # offline, from local copies
  uv run python scripts/count_published_splits.py --skip-pmpnn                # JSONs only
  uv run python scripts/count_published_splits.py --json counts.json          # machine-readable
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import tarfile
from collections import Counter
from pathlib import Path
from typing import Any

import httpx

LMPNN_RAW = "https://raw.githubusercontent.com/dauparas/LigandMPNN/main/training"
LMPNN_FILES = [
    "train.json",
    "valid.json",
    "test_small_molecule.json",
    "test_nucleotide.json",
    "test_metal.json",
]
PMPNN_TAR = "https://files.ipd.uw.edu/pub/training_sets/pdb_2021aug02_sample.tar.gz"
# Held-out lists take precedence when an id is listed twice, which is what makes a
# train/test double-listing a leak rather than a bookkeeping quirk.
TRAIN_LISTS = {"train.json", "valid.json"}


def _load_lmpnn(local: Path | None) -> dict[str, list[str]]:
    if local is not None:
        return {f: json.loads((local / f).read_text()) for f in LMPNN_FILES}
    out: dict[str, list[str]] = {}
    with httpx.Client(timeout=120, follow_redirects=True) as client:
        for fname in LMPNN_FILES:
            r = client.get(f"{LMPNN_RAW}/{fname}")
            r.raise_for_status()
            out[fname] = r.json()
    return out


def count_lmpnn(local: Path | None) -> dict[str, Any]:
    lists = _load_lmpnn(local)
    sizes = {f: len(v) for f, v in lists.items()}
    all_ids = [i for v in lists.values() for i in v]
    counts = Counter(all_ids)

    overlaps = []
    for entry_id in sorted(k for k, n in counts.items() if n > 1):
        where = sorted(f for f, v in lists.items() if entry_id in set(v))
        held_out = [f for f in where if f not in TRAIN_LISTS]
        overlaps.append(
            {
                "entry_id": entry_id,
                "lists": where,
                # a train/test straddle is a leak; test/test is double-counting
                "train_test_straddle": bool(held_out) and len(held_out) < len(where),
            }
        )

    return {
        "sizes": sizes,
        "total_ids": len(all_ids),
        "distinct_ids": len(counts),
        "overlaps": overlaps,
        "straddles": [o["entry_id"] for o in overlaps if o["train_test_straddle"]],
        # ids are lowercase in the released files; the audit upper-cases before
        # querying RCSB, which is why a case-sensitive check would miss these
        "ids_are_lowercase": all(i == i.lower() for i in all_ids),
    }


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
    ap.add_argument("--lmpnn-dir", type=Path, default=None, help="local dir with the 5 JSONs")
    ap.add_argument("--pmpnn-tar", type=Path, default=None, help="local pdb_2021aug02 tarball")
    ap.add_argument("--skip-pmpnn", action="store_true", help="skip the 47 MB download")
    ap.add_argument("--json", type=Path, default=None, help="also write the counts here")
    args = ap.parse_args()

    report: dict[str, Any] = {}

    lm = count_lmpnn(args.lmpnn_dir)
    report["ligandmpnn"] = lm
    print("LigandMPNN (github.com/dauparas/LigandMPNN/training)")
    for fname in LMPNN_FILES:
        print(f"  {fname:<28} {lm['sizes'][fname]:>7,}")
    print(f"  {'total ids':<28} {lm['total_ids']:>7,}")
    print(f"  {'distinct ids':<28} {lm['distinct_ids']:>7,}")
    if lm["overlaps"]:
        print(f"  ids in more than one list:   {len(lm['overlaps'])}")
        for o in lm["overlaps"]:
            flag = "  <-- TRAIN/TEST STRADDLE" if o["train_test_straddle"] else ""
            print(f"    {o['entry_id']}  {', '.join(o['lists'])}{flag}")

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
