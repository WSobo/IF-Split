"""Stage 1 - Enumerate candidates from RCSB (Search v2 + Data API).

Selects entries by ``release_date <= snapshot_date`` plus the method/resolution
filters, enriches each via the Data API (sequences, ligand comps, residue
counts, assemblies), and writes the byte-stable ``candidates.jsonl`` -- the
snapshot definition. No coordinates are downloaded (PLAN.md §1.5).
"""

from __future__ import annotations

import sys
import time
import warnings
from collections.abc import Callable
from pathlib import Path

from .config import Config
from .rcsb import RcsbClient
from .schema import CandidateRecord, canonical_jsonl_bytes, sha256_hex

ProgressFn = Callable[[str], None]

# How often (in enriched records) to emit a progress line.
_REPORT_EVERY = 1000

# Methods RCSB assigns no resolution to. The Search query always carries a resolution
# predicate, so an entry of one of these methods can never match — requesting one alone
# silently yields an empty snapshot. We WARN (excluding NMR is a defensible inverse-
# folding choice, per PLAN.md) rather than fail, so the footgun is at least visible.
RESOLUTIONLESS_METHODS: frozenset[str] = frozenset(
    {
        "SOLUTION NMR", "SOLID-STATE NMR", "SOLUTION SCATTERING", "EPR",
        "FLUORESCENCE TRANSFER", "INFRARED SPECTROSCOPY", "THEORETICAL MODEL",
    }
)  # fmt: skip


def fmt_duration(seconds: float) -> str:
    """Human-friendly duration: ``45s`` / ``3m29s`` / ``1h02m``."""
    s = int(max(0, seconds))
    if s < 60:
        return f"{s}s"
    m, s = divmod(s, 60)
    if m < 60:
        return f"{m}m{s:02d}s"
    h, m = divmod(m, 60)
    return f"{h}h{m:02d}m"


def progress_line(label: str, n: int, total: int, t0: float) -> str:
    """A ``label: n/total (pct)  rate/s  eta ...`` line from a monotonic start."""
    elapsed = time.monotonic() - t0
    pct = (100 * n / total) if total else 0.0
    msg = f"{label}: {n}/{total} ({pct:.0f}%)"
    if elapsed > 0 and n > 0:
        rate = n / elapsed
        eta = (total - n) / rate if rate > 0 else 0.0
        msg += f"  {rate:.0f}/s  eta {fmt_duration(eta)}"
    return msg


def make_console_progress(stream=None) -> ProgressFn:
    """A timestamped, line-flushed progress printer for long CLI runs.

    The flush is the important part: when stdout is redirected to a file, Python
    block-buffers it, so unflushed progress lines stay invisible until the process
    exits. This forces each line out immediately.
    """
    out = stream or sys.stdout

    def say(msg: str) -> None:
        print(f"  [{time.strftime('%H:%M:%S')}] {msg}", file=out, flush=True)

    return say


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

    resless = sorted(set(cfg.experimental_methods) & RESOLUTIONLESS_METHODS)
    if resless:
        warnings.warn(
            f"experimental_methods includes {resless}, which RCSB assigns no resolution; "
            f"the resolution filter (<= {cfg.search_resolution_cap()} A) excludes them, so "
            f"they contribute no entries. Remove them, or expect them to be absent.",
            stacklevel=2,
        )
        say(f"WARNING: {resless} have no RCSB resolution and will contribute no entries")

    owns_client = client is None
    client = client or RcsbClient()
    try:
        ids = client.search_entry_ids(cfg, limit=limit, progress=progress)
        say(f"search: {len(ids)} entries match the snapshot")

        records: list[CandidateRecord] = []
        total = len(ids)
        t0 = time.monotonic()
        for raw in client.fetch_entries(ids):
            records.append(CandidateRecord.from_data_api(raw))
            if len(records) % _REPORT_EVERY == 0:
                say(progress_line("enriched", len(records), total, t0))
        say(progress_line("enriched", len(records), total, t0) + " (done)")
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
