"""Stage 3 - Filter candidates on metadata (no coordinates).

Operates on the records in ``candidates.jsonl``. Drops entries that:

- have no protein polymer entity (``no_protein_entity``),
- have a protein entity but no usable sequence (``no_protein_sequence``),
- exceed the residue cap (``too_large``): ``total_residues >= max_total_residues``.

When ``use_biological_assembly`` the residue count is taken from assembly 1
(``<entry>-1``); otherwise the deposited polymer monomer count is used. Every
drop is recorded with its reason so the build is auditable.
"""

from __future__ import annotations

from .config import Config
from .schema import CandidateRecord

DROP_NO_PROTEIN = "no_protein_entity"
DROP_NO_SEQUENCE = "no_protein_sequence"
DROP_TOO_LARGE = "too_large"


def assembly1_residue_count(record: CandidateRecord) -> int | None:
    """Residue count of biological assembly 1 (id ending ``-1``), else smallest."""
    if not record.assemblies:
        return None
    for aid in sorted(record.assemblies):
        if aid.endswith("-1"):
            return record.assemblies[aid]
    return record.assemblies[sorted(record.assemblies)[0]]


def total_residues(record: CandidateRecord, cfg: Config) -> int | None:
    """Residue count used for the size filter, per the assembly config."""
    if cfg.use_biological_assembly:
        count = assembly1_residue_count(record)
        if count is not None:
            return count
    return record.deposited_residues


def filter_candidates(
    records: list[CandidateRecord], cfg: Config
) -> tuple[list[CandidateRecord], list[dict]]:
    """Return ``(kept, drops)`` where drops is a list of ``{entry_id, reason, ...}``."""
    kept: list[CandidateRecord] = []
    drops: list[dict] = []
    for r in records:
        proteins = [e for e in r.polymer_entities if e.is_protein]
        if not proteins:
            drops.append({"entry_id": r.entry_id, "reason": DROP_NO_PROTEIN})
            continue
        if not any(e.seq for e in proteins):
            drops.append({"entry_id": r.entry_id, "reason": DROP_NO_SEQUENCE})
            continue
        tr = total_residues(r, cfg)
        if tr is not None and tr >= cfg.max_total_residues:
            drops.append({"entry_id": r.entry_id, "reason": DROP_TOO_LARGE, "residues": tr})
            continue
        kept.append(r)
    kept.sort(key=lambda r: r.entry_id)
    drops.sort(key=lambda d: d["entry_id"])
    return kept, drops


def drop_summary(drops: list[dict]) -> dict[str, int]:
    """Count drops by reason (deterministic order via sorted keys downstream)."""
    out: dict[str, int] = {}
    for d in drops:
        out[d["reason"]] = out.get(d["reason"], 0) + 1
    return out
