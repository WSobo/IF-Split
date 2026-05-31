"""Stage 1 - Enumerate candidates from RCSB (Search v2 + Data API).

Selects entries by ``release_date <= snapshot_date`` plus the method/resolution
filters, enriches each via the Data API (sequences, ligand comps, residue
counts, assemblies), and writes the byte-stable ``candidates.jsonl`` -- the
snapshot definition. No coordinates are downloaded (PLAN.md §1.5).
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from .config import Config
from .rcsb import RcsbClient
from .schema import CandidateRecord, canonical_jsonl_bytes, sha256_hex

ProgressFn = Callable[[str], None]


def enumerate_candidates(
    cfg: Config,
    out_dir: str | Path,
    *,
    limit: int | None = None,
    client: RcsbClient | None = None,
    progress: ProgressFn | None = None,
) -> tuple[list[CandidateRecord], Path, str]:
    """Run Stage 1.

    Returns ``(records, candidates_path, sha256)``. ``candidates.jsonl`` is
    written to ``out_dir`` in canonical (byte-stable) form.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    def say(msg: str) -> None:
        if progress:
            progress(msg)

    owns_client = client is None
    client = client or RcsbClient()
    try:
        ids = client.search_entry_ids(cfg, limit=limit)
        say(f"search: {len(ids)} entries match the snapshot")

        records: list[CandidateRecord] = []
        for raw in client.fetch_entries(ids):
            records.append(CandidateRecord.from_data_api(raw))
            if len(records) % 1000 == 0:
                say(f"enriched: {len(records)}/{len(ids)}")
        say(f"enriched: {len(records)}/{len(ids)} (done)")
    finally:
        if owns_client:
            client.close()

    data = canonical_jsonl_bytes(records)
    sha = sha256_hex(data)
    candidates_path = out_dir / "candidates.jsonl"
    candidates_path.write_bytes(data)
    say(f"wrote {candidates_path} ({len(records)} records, sha256={sha[:12]}...)")

    # Return in canonical (entry_id-sorted) order so callers see what was written.
    records.sort(key=lambda r: r.entry_id)
    return records, candidates_path, sha
