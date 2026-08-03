"""Fetch RCSB metadata for an EXTERNAL split's entry ids (for auditing it).

Given the published train/val/test entry-id lists of another dataset (LigandMPNN's
train.json / valid.json / test_*.json by default), download them, then pull each
entry's RCSB metadata — folds (CATH/ECOD/SCOP2), 30% sequence clusters, ligand
signals — into a canonical ``candidates.jsonl`` using IF-Split's own RcsbClient.
Metadata only; no coordinates. Resumable + threaded.

The output feeds ``scripts/audit_ligandmpnn_split.py``, which measures the fold-,
ligand-context- and complex-bridge leakage the external split admits.

Outputs (in --out):
  <name>_splits.json      entry_id(UPPER) -> {"split": ..., "test_class": ...}
  <name>_candidates.jsonl canonical CandidateRecord per resolvable entry
  <name>_missing.txt      ids RCSB did not return (obsoleted/superseded)

Usage:
  # default: audit LigandMPNN's published split (downloads the 5 JSONs)
  uv run python scripts/fetch_external_split.py --out /tmp/lmpnn_audit
  # or point at a local dir already holding the 5 JSON files
  uv run python scripts/fetch_external_split.py --local path/to/jsons --out /tmp/lmpnn_audit
"""

from __future__ import annotations

import argparse
import contextlib
import json
import queue
import threading
from pathlib import Path

import httpx

from ifsplit.rcsb import DATA_BATCH_SIZE, RcsbClient
from ifsplit.schema import CandidateRecord

# LigandMPNN's published split (github.com/dauparas/LigandMPNN/training). Entry-level,
# ligand-class-stratified test set — the set the paper's Fig. 2a recovery is computed on.
LMPNN_RAW = "https://raw.githubusercontent.com/dauparas/LigandMPNN/main/training"
SPLIT_FILES = {
    # local_name: (remote_file, split_label, test_class)
    "train": ("train.json", "train", None),
    "valid": ("valid.json", "valid", None),
    "test_sm": ("test_small_molecule.json", "test", "small_molecule"),
    "test_nuc": ("test_nucleotide.json", "test", "nucleotide"),
    "test_metal": ("test_metal.json", "test", "metal"),
}


def load_or_download(local: Path | None) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    with httpx.Client(timeout=60) as client:
        for _key, (fname, _split, _klass) in SPLIT_FILES.items():
            if local is not None:
                with open(local / fname) as lf:
                    out[fname] = json.load(lf)
            else:
                r = client.get(f"{LMPNN_RAW}/{fname}")
                r.raise_for_status()
                out[fname] = r.json()
    return out


def build_splits(files: dict[str, list[str]]) -> dict[str, dict]:
    """entry_id(UPPER) -> {split, test_class}. test files override train/valid."""
    labels: dict[str, dict] = {}
    # apply in an order where later (more specific) entries override earlier ones
    order = ["train", "valid", "test_sm", "test_nuc", "test_metal"]
    for key in order:
        fname, split, klass = SPLIT_FILES[key]
        for raw in files[fname]:
            labels[raw.upper()] = {"split": split, "test_class": klass}
    return labels


def fetch(all_ids: list[str], cand_path: Path, n_workers: int) -> list[str]:
    done: set[str] = set()
    if cand_path.exists():
        with open(cand_path) as rf:
            for line in rf:
                line = line.strip()
                if line:
                    with contextlib.suppress(Exception):
                        done.add(json.loads(line)["entry_id"].upper())
    todo = [i for i in all_ids if i not in done]
    print(f"todo {len(todo)} / {len(all_ids)} (resume: {len(done)})", flush=True)

    batches: queue.Queue = queue.Queue()
    for i in range(0, len(todo), DATA_BATCH_SIZE):
        batches.put(todo[i : i + DATA_BATCH_SIZE])
    lock = threading.Lock()
    got_all: set[str] = set()
    counter = {"n": len(done)}

    with open(cand_path, "a", encoding="utf-8") as fh:

        def worker() -> None:
            client = RcsbClient()
            try:
                while True:
                    try:
                        batch = batches.get_nowait()
                    except queue.Empty:
                        return
                    try:
                        for entry in client.fetch_entries(batch):
                            rec = CandidateRecord.from_data_api(entry)
                            with lock:
                                fh.write(rec.to_canonical_json() + "\n")
                                got_all.add(rec.entry_id.upper())
                    except Exception as exc:
                        with lock:
                            print(f"  batch error: {exc}", flush=True)
                    with lock:
                        counter["n"] += len(batch)
                        if counter["n"] % 8000 < DATA_BATCH_SIZE:
                            print(f"  {counter['n']}/{len(all_ids)}", flush=True)
                        fh.flush()
                    batches.task_done()
            finally:
                client.close()

        threads = [threading.Thread(target=worker, daemon=True) for _ in range(n_workers)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()
    return sorted(set(all_ids) - done - got_all)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--local", type=Path, default=None, help="dir with the 5 split JSONs")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--name", default="lmpnn")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    files = load_or_download(args.local)
    labels = build_splits(files)
    with open(args.out / f"{args.name}_splits.json", "w") as sf:
        json.dump(labels, sf, sort_keys=True)
    all_ids = sorted(labels)
    print(f"{args.name}: {len(all_ids)} unique entry ids", flush=True)

    missing = fetch(all_ids, args.out / f"{args.name}_candidates.jsonl", args.workers)
    (args.out / f"{args.name}_missing.txt").write_text("\n".join(missing) + "\n")
    print(f"DONE. unresolved (obsoleted): {len(missing)}", flush=True)


if __name__ == "__main__":
    main()
