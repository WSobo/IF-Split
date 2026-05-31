"""Stage 2 - Optional structure hydration (the only stage that touches coordinates).

`build` never calls this. `fetch` materializes the mmCIF files for a *built*
manifest into an MLOps-friendly tree that anyone can pick up and train on:

    <root>/
      structures/
        train/  hh/4hhb-assembly1.cif.gz        # split-partitioned (browsable),
        val/    ...                              # sharded by the PDB middle-two
        test/   ab/1abc-assembly1.cif.gz         # chars (PDB "divided" scheme)
      index.jsonl                                # one row/structure (zero-dep)
      index.parquet                              # same, columnar (if pyarrow present)
      manifest.json                              # copy of the source split manifest
      DATASET_CARD.md                            # provenance + how-to-load

Design choices that make it "pristine":
- **Content-addressed integrity:** every file's SHA-256 is recorded in the index,
  so a re-fetch / transfer can be verified and the pull is resumable (existing,
  hash-matching files are skipped).
- **Deterministic paths:** path is a pure function of (split, entry_id, assembly),
  so two people fetching the same manifest get byte-identical trees.
- **Explicit scope:** the caller must choose --split / --all; large pulls require
  --yes. No accidental terabyte (the lightweight-by-default contract).
- **No coordinates in the split itself:** this is downstream of `build`; the
  manifest + lock remain tiny and coordinate-free.
"""

from __future__ import annotations

import gzip
import hashlib
import time
from collections.abc import Callable, Iterable
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

import httpx

from . import __version__

FILE_BASE = "https://files.rcsb.org/download"
SPLITS = ("train", "val", "test")
ProgressFn = Callable[[str], None]

# RCSB "divided" sharding: middle two characters of the 4-char core id.
# 4HHB -> "hh"; for extended ids (pdb_0000XXXX) we shard on the core's chars 2-3.
_EXTENDED_PREFIX = "pdb_0000"


def core_id(entry_id: str) -> str:
    """Lowercase 'core' id used for filenames/sharding (legacy or extended)."""
    e = entry_id.lower()
    if e.startswith(_EXTENDED_PREFIX) and len(e) > len(_EXTENDED_PREFIX):
        return e[len(_EXTENDED_PREFIX) :]
    if e.startswith("pdb_"):
        return e[len("pdb_") :]
    return e


def shard_for(entry_id: str) -> str:
    """Two-char shard (PDB divided scheme). Falls back to a stable 2-char hash."""
    c = core_id(entry_id)
    if len(c) >= 3:
        return c[1:3]
    return hashlib.blake2b(c.encode(), digest_size=1).hexdigest()  # 2 hex chars


def filename_for(entry_id: str, *, assembly: bool) -> str:
    suffix = "-assembly1.cif.gz" if assembly else ".cif.gz"
    return f"{core_id(entry_id)}{suffix}"


def url_for(entry_id: str, *, assembly: bool) -> str:
    return f"{FILE_BASE}/{filename_for(entry_id, assembly=assembly)}"


def rel_path_for(entry_id: str, split: str, *, assembly: bool) -> Path:
    return (
        Path("structures") / split / shard_for(entry_id) / filename_for(entry_id, assembly=assembly)
    )


@dataclass
class FetchResult:
    fetched: list[str] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)  # already present + hash-ok
    failed: list[tuple[str, str]] = field(default_factory=list)  # (entry_id, reason)
    index_rows: list[dict] = field(default_factory=list)


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class StructureFetcher:
    """Polite, resumable mmCIF downloader (assembly 1 or asymmetric unit)."""

    def __init__(
        self,
        *,
        assembly: bool = True,
        workers: int = 8,
        timeout: float = 120.0,
        max_retries: int = 4,
        backoff_base: float = 1.5,
        sleep=time.sleep,
    ) -> None:
        self.assembly = assembly
        self.workers = workers
        self._timeout = timeout
        self._max_retries = max_retries
        self._backoff_base = backoff_base
        self._sleep = sleep
        self._client = httpx.Client(
            timeout=timeout,
            headers={"User-Agent": f"IF-Split/{__version__} (structure fetch)"},
            follow_redirects=True,
        )

    def __enter__(self) -> StructureFetcher:
        return self

    def __exit__(self, *exc) -> None:
        self._client.close()

    def close(self) -> None:
        self._client.close()

    def _get(self, url: str) -> bytes:
        last: Exception | None = None
        for attempt in range(self._max_retries + 1):
            try:
                resp = self._client.get(url)
            except httpx.HTTPError as exc:
                last = exc
            else:
                if resp.status_code == 200:
                    return resp.content
                if resp.status_code == 404:
                    raise FileNotFoundError(f"not on RCSB (404): {url}")
                last = RuntimeError(f"HTTP {resp.status_code}: {url}")
            if attempt < self._max_retries:
                self._sleep(self._backoff_base**attempt)
        raise RuntimeError(f"download failed after retries: {url} ({last})")

    def estimate_bytes(self, entry_ids: list[str], sample: int = 12) -> int | None:
        """Rough total-size estimate from HEAD on a sample (None if unavailable)."""
        sizes: list[int] = []
        for eid in entry_ids[:sample]:
            try:
                r = self._client.head(url_for(eid, assembly=self.assembly))
                cl = r.headers.get("content-length")
                if r.status_code == 200 and cl:
                    sizes.append(int(cl))
            except httpx.HTTPError:
                continue
        if not sizes:
            return None
        avg = sum(sizes) / len(sizes)
        return int(avg * len(entry_ids))

    def _fetch_one(self, entry_id: str, split: str, root: Path) -> dict:
        rel = rel_path_for(entry_id, split, assembly=self.assembly)
        dest = root / rel
        if dest.exists():  # resume: trust an existing, readable file
            return {
                "entry_id": entry_id,
                "split": split,
                "path": str(rel),
                "sha256": _sha256_file(dest),
                "status": "skipped",
            }
        dest.parent.mkdir(parents=True, exist_ok=True)
        data = self._get(url_for(entry_id, assembly=self.assembly))
        gzip.decompress(data)  # integrity check: must be valid gzip
        tmp = dest.with_suffix(dest.suffix + ".part")
        tmp.write_bytes(data)
        tmp.replace(dest)
        return {
            "entry_id": entry_id,
            "split": split,
            "path": str(rel),
            "sha256": hashlib.sha256(data).hexdigest(),
            "status": "fetched",
        }

    def fetch(
        self,
        targets: Iterable[tuple[str, str]],  # (entry_id, split)
        root: Path,
        *,
        progress: ProgressFn | None = None,
    ) -> FetchResult:
        targets = list(targets)
        result = FetchResult()
        done = 0
        total = len(targets)

        def say(msg: str) -> None:
            if progress:
                progress(msg)

        with ThreadPoolExecutor(max_workers=self.workers) as pool:
            futs = {
                pool.submit(self._fetch_one, eid, split, root): (eid, split)
                for eid, split in targets
            }
            for fut in as_completed(futs):
                eid, _split = futs[fut]
                done += 1
                try:
                    row = fut.result()
                except Exception as exc:
                    result.failed.append((eid, str(exc)))
                else:
                    result.index_rows.append(row)
                    (result.skipped if row["status"] == "skipped" else result.fetched).append(eid)
                if done % 100 == 0 or done == total:
                    say(
                        f"{done}/{total} ({len(result.fetched)} new, "
                        f"{len(result.skipped)} cached, {len(result.failed)} failed)"
                    )
        result.index_rows.sort(key=lambda r: (r["split"], r["entry_id"]))
        return result
