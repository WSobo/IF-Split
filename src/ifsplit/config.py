"""Load, validate, and hash an IF-Split run configuration.

The config is the single source of truth for a build. Its canonical hash is
embedded in the manifest so that two manifests sharing a config hash are
guaranteed to have used identical, output-affecting settings.

The hash is computed over the *validated, normalized* settings (not the raw YAML
text), so comments and formatting do not change it — only values do.
"""

from __future__ import annotations

import hashlib
import json
from datetime import date
from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class SplitFractions(BaseModel):
    """Train/val/test partition fractions; must sum to 1.0."""

    model_config = ConfigDict(extra="forbid")

    train: float = Field(gt=0, lt=1)
    val: float = Field(gt=0, lt=1)
    test: float = Field(gt=0, lt=1)

    @model_validator(mode="after")
    def _sum_to_one(self) -> SplitFractions:
        total = self.train + self.val + self.test
        if abs(total - 1.0) > 1e-9:
            raise ValueError(f"split_fractions must sum to 1.0, got {total}")
        return self


class Config(BaseModel):
    """A fully-validated IF-Split build configuration."""

    model_config = ConfigDict(extra="forbid")

    # --- snapshot definition (reproducibility anchor) ---
    snapshot_date: date
    experimental_methods: list[str] = Field(min_length=1)
    resolution_max_A: float = Field(gt=0)
    max_total_residues: int = Field(gt=0)
    excluded_het: list[str] = Field(default_factory=list)
    use_biological_assembly: bool = True

    # --- curation: purification-artifact detection (Stage 4) ---
    # A poly-His tag coordinating Ni/Co is a purification artifact, not a
    # biological metal site (a known blemish in the LigandMPNN metal set). An
    # entry whose *only* metal is a purification metal AND that carries a His-tag
    # is flagged so it can be excluded from the metal class. Empty
    # purification_metals disables the heuristic.
    purification_metals: list[str] = Field(default_factory=lambda: ["NI", "CO"])
    histag_min_run: int = Field(default=6, gt=0)
    exclude_purification_artifacts: bool = True

    # --- clustering + split ---
    identity_threshold: float = Field(gt=0, le=1)
    # "precomputed": reuse RCSB's entity clusters (default, no external binary).
    # "mmseqs2": run our own over the snapshot's sequences.
    clustering_backend: Literal["precomputed", "mmseqs2"] = "precomputed"
    split_fractions: SplitFractions
    split_salt: str = Field(min_length=1)
    seed: int = Field(ge=0)

    # --- test-set minimums (opt-in stratification top-up) ---
    # Floor on the number of test entries carrying each functional ligand class,
    # e.g. {"metal": 500, "nucleotide": 200}. Empty (default) = pure deterministic
    # hash, no top-up. When set, after the base assignment any class below its floor
    # recruits *whole components* (never individual entries, so no leakage) into
    # test in deterministic hash order, skipping components already pinned by a
    # registry (so growth stays stable). A floor larger than the available supply
    # is satisfied as far as possible and the shortfall is reported, not forced.
    test_min_per_class: dict[str, int] = Field(default_factory=dict)

    # --- quality filters (Stage 3): wwPDB validation-report metrics ---
    # Metadata only (no coordinates). Each cap is optional (None disables it). An
    # entry is dropped only when the metric is present AND violates the cap;
    # entries whose report lacks a metric are kept (never penalized for an absent
    # metric). Geometry caps (clashscore/Ramachandran/rotamer) apply to X-ray and
    # cryo-EM; diffraction caps (R-free/RSRZ) naturally no-op on EM entries.
    max_clashscore: float | None = Field(default=None, gt=0)
    max_rfree: float | None = Field(default=None, gt=0)
    max_ramachandran_outlier_pct: float | None = Field(default=None, ge=0)
    max_rotamer_outlier_pct: float | None = Field(default=None, ge=0)
    max_rsrz_outlier_pct: float | None = Field(default=None, ge=0)
    require_validation_report: bool = False

    # --- featurization (downstream-optional; not part of the split definition) ---
    ligand_context_radius_A: float = Field(gt=0)
    max_ligand_atoms: int = Field(gt=0)

    @field_validator("experimental_methods")
    @classmethod
    def _normalize_methods(cls, v: list[str]) -> list[str]:
        return [m.strip().upper() for m in v]

    @field_validator("excluded_het", "purification_metals")
    @classmethod
    def _normalize_codes(cls, v: list[str]) -> list[str]:
        return [h.strip().upper() for h in v]

    @field_validator("test_min_per_class")
    @classmethod
    def _check_minimums(cls, v: dict[str, int]) -> dict[str, int]:
        for k, n in v.items():
            if n < 0:
                raise ValueError(f"test_min_per_class[{k!r}] must be >= 0, got {n}")
        return v

    @property
    def dataset_version(self) -> str:
        """Versioned dataset name, e.g. 'IF-Split-2026.05.30'."""
        return f"IF-Split-{self.snapshot_date:%Y.%m.%d}"

    @property
    def identity_level(self) -> int:
        """``identity_threshold`` as an integer percent (e.g. 0.30 -> 30)."""
        return round(self.identity_threshold * 100)

    def canonical_dict(self) -> dict[str, Any]:
        """JSON-mode dump (dates -> ISO strings) used for hashing and manifests."""
        return self.model_dump(mode="json")

    def config_hash(self) -> str:
        """Deterministic, formatting-independent hash of the settings."""
        canonical = json.dumps(self.canonical_dict(), sort_keys=True, separators=(",", ":"))
        return hashlib.blake2b(canonical.encode("utf-8"), digest_size=16).hexdigest()


def load_config(path: str | Path) -> Config:
    """Load and validate a YAML config file into a :class:`Config`."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config not found: {path}")
    with path.open("r", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    if not isinstance(raw, dict):
        raise ValueError(f"Config must be a YAML mapping, got {type(raw).__name__}: {path}")
    return Config.model_validate(raw)
