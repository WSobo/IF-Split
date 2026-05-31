"""Stage 3 - Parse mmCIF (gemmi) into structured records and apply filters.

Extracts protein/nucleic chains and non-polymer entities, normalizes modified
residues to canonical parents (e.g. MSE->MET), drops oversized/empty entries,
and records drop reasons. Lands in Phase 3.
"""

from __future__ import annotations

from .config import Config


def parse_and_filter(cfg: Config, *args, **kwargs):
    raise NotImplementedError("Stage 3 (parse + filter) lands in Phase 3.")
