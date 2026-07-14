"""Stage 3 - Filter candidates on metadata (no coordinates).

Operates on the records in ``candidates.jsonl``. Drops entries that:

- have no protein polymer entity (``no_protein_entity``),
- have a protein entity but no usable sequence (``no_protein_sequence``): the
  canonical sequence is empty, or every protein chain is all-``X`` (poly-UNK, no
  known residue identities), so there is no learnable inverse-folding label,
- (opt-in ``min_modeled_residues``) have no protein chain with at least that many
  modeled (non-``X``) residues (``sequence_too_short``): tiny peptide fragments and
  mostly-unknown chains,
- exceed the residue cap (``too_large``): ``total_residues > max_total_residues``,
- exceed the (per-method) resolution cap (``resolution_too_low``): re-derived here
  from ``resolution_A`` so the cut is auditable from ``candidates.jsonl`` and can be
  tightened offline; ``resolution_max_A_by_method`` lets X-ray and cryo-EM differ,
- violate an (optional) wwPDB validation-report quality cap — clashscore, R-free,
  Ramachandran/rotamer/RSRZ outliers, EM map-fit floor — or lack a report when
  one is required.

When ``use_biological_assembly`` the residue count is taken from assembly 1
(``<entry>-1``); otherwise the deposited polymer monomer count is used. Quality
metrics come from the metadata API (no coordinates); a cap fires only when both
the cap and the metric are present, so a missing metric never drops an entry.
Every drop is recorded with its reason so the build is auditable.
"""

from __future__ import annotations

from .config import Config
from .schema import CandidateRecord

DROP_NO_PROTEIN = "no_protein_entity"
DROP_NO_SEQUENCE = "no_protein_sequence"
DROP_SEQUENCE_TOO_SHORT = "sequence_too_short"
DROP_TOO_LARGE = "too_large"
DROP_RESOLUTION = "resolution_too_low"
DROP_CLASHSCORE = "clashscore_too_high"
DROP_RFREE = "rfree_too_high"
DROP_RAMACHANDRAN = "ramachandran_outliers_too_high"
DROP_ROTAMER = "rotamer_outliers_too_high"
DROP_RSRZ = "rsrz_outliers_too_high"
DROP_EM_INCLUSION = "em_backbone_inclusion_too_low"
DROP_NO_VALIDATION = "no_validation_report"


def effective_resolution_cap(record: CandidateRecord, cfg: Config) -> float:
    """The resolution cap that applies to ``record``, honoring per-method overrides.

    An entry is kept if it passes the *loosest* cap among its enabled methods (so a
    dual X-ray/EM deposition qualifies if it is a good structure by either route).
    With no ``resolution_max_A_by_method`` overrides this is just ``resolution_max_A``.
    """
    caps = [cfg.method_resolution_cap(m) for m in record.methods if m in cfg.experimental_methods]
    return max(caps) if caps else cfg.resolution_max_A


def modeled_residue_count(seq: str) -> int:
    """Number of residues with a known identity (non-``X``) in a canonical sequence.

    ``pdbx_seq_one_letter_code_can`` maps modified residues to their standard parent
    and leaves only genuine unknowns as ``X``, so an all-``X`` chain (poly-UNK) has
    zero modeled residues and no learnable inverse-folding label.
    """
    return sum(1 for ch in seq if ch not in ("X", "x"))


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


def quality_drop(record: CandidateRecord, cfg: Config) -> tuple[str, float] | None:
    """First validation-report cap this record violates, else ``None``.

    A cap fires only when both the cap and the metric are present and the metric
    exceeds the cap; a missing metric never drops the entry. With
    ``require_validation_report`` an entry that has no report at all is dropped.
    Returns ``(reason, value)`` so the drop log records the offending number.
    """
    q = record.quality
    if cfg.require_validation_report and not q.has_report:
        return (DROP_NO_VALIDATION, 0.0)
    checks = (
        (cfg.max_clashscore, q.clashscore, DROP_CLASHSCORE),
        (cfg.max_rfree, q.rfree, DROP_RFREE),
        (cfg.max_ramachandran_outlier_pct, q.ramachandran_outlier_pct, DROP_RAMACHANDRAN),
        (cfg.max_rotamer_outlier_pct, q.rotamer_outlier_pct, DROP_ROTAMER),
        (cfg.max_rsrz_outlier_pct, q.rsrz_outlier_pct, DROP_RSRZ),
    )
    for cap, value, reason in checks:
        if cap is not None and value is not None and value > cap:
            return (reason, value)
    # EM map-model agreement is a FLOOR (higher is better), so it drops BELOW the cap.
    if (
        cfg.min_em_backbone_inclusion is not None
        and q.em_backbone_inclusion is not None
        and q.em_backbone_inclusion < cfg.min_em_backbone_inclusion
    ):
        return (DROP_EM_INCLUSION, q.em_backbone_inclusion)
    return None


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
        # Usable sequence = at least one protein chain with >=1 modeled residue. An
        # empty string OR an all-X (poly-UNK) chain has none, so the label is
        # unlearnable — always dropped (categorically, like an absent sequence).
        best_modeled = max((modeled_residue_count(e.seq) for e in proteins), default=0)
        if best_modeled == 0:
            drops.append({"entry_id": r.entry_id, "reason": DROP_NO_SEQUENCE})
            continue
        # Opt-in floor on modeled length: drop tiny fragments / mostly-unknown chains.
        if cfg.min_modeled_residues > 0 and best_modeled < cfg.min_modeled_residues:
            drops.append(
                {"entry_id": r.entry_id, "reason": DROP_SEQUENCE_TOO_SHORT, "modeled": best_modeled}
            )
            continue
        tr = total_residues(r, cfg)
        if tr is not None and tr > cfg.max_total_residues:
            drops.append({"entry_id": r.entry_id, "reason": DROP_TOO_LARGE, "residues": tr})
            continue
        # Re-derive the resolution cut here (Search applied it too) so it is auditable
        # from candidates.jsonl and can be tightened offline; per-method caps apply.
        if r.resolution_A is not None and r.resolution_A > effective_resolution_cap(r, cfg):
            drops.append(
                {"entry_id": r.entry_id, "reason": DROP_RESOLUTION, "resolution": r.resolution_A}
            )
            continue
        qd = quality_drop(r, cfg)
        if qd is not None:
            reason, value = qd
            drops.append({"entry_id": r.entry_id, "reason": reason, "value": value})
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
