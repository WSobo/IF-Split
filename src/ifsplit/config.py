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

SPEC_SCHEMA = "ifsplit/config@1"

# RCSB publishes precomputed polymer-entity clusters at only these identity
# percentages, so the "precomputed" backend can look up a cluster id only at one
# of them. Any other level would match nothing and silently make every entity a
# singleton (no clustering, hence cross-split leakage) — so we reject it up front.
RCSB_IDENTITY_LEVELS: frozenset[int] = frozenset({30, 50, 70, 90, 95, 100})


class SpecMeta(BaseModel):
    """Self-identifying header + human metadata for a shareable split spec.

    These fields are *descriptive only* and are deliberately EXCLUDED from
    ``config_hash`` — two identical splits with different names/authors must still
    share a hash. ``expected_config_hash`` lets a shared spec self-verify: on load,
    if it is set and does not match the computed hash, the user is warned.
    """

    model_config = ConfigDict(extra="forbid")

    ifsplit_spec: str = SPEC_SCHEMA  # schema id, so the file announces what it is
    name: str | None = None
    description: str | None = None
    author: str | None = None
    created_with: str | None = None  # e.g. "if-split 0.1.0"
    expected_config_hash: str | None = None


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
    # Opt-in sequence-usability floor (Stage 3). Keep an entry only if some protein
    # chain has at least this many modeled (non-'X') residues. 0 (default) = off, so
    # only the always-on empty/all-'X' (poly-UNK) drop applies. A modest value (e.g.
    # 20) removes tiny peptide fragments and mostly-unknown chains whose sequence is a
    # poor inverse-folding label; it only ever drops entries whose *every* protein
    # chain is that short (any() over chains), so multi-chain complexes are unaffected.
    min_modeled_residues: int = Field(default=0, ge=0)
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
    # A shorter His run counts as a tag *only if it sits at a chain terminus*:
    # a 6xHis tag with a few residues unmodeled or trimmed from the deposited
    # sequence still leaves a short terminal His run. Internal His clusters (real
    # metalloprotein motifs) are not flagged. 0 disables the terminal rule.
    histag_terminal_min_run: int = Field(default=3, ge=0)
    exclude_purification_artifacts: bool = True

    # --- clustering + split ---
    identity_threshold: float = Field(gt=0, le=1)
    # "precomputed": reuse RCSB's entity clusters (default, no external binary).
    # "mmseqs2": run our own over the snapshot's sequences.
    clustering_backend: Literal["precomputed", "mmseqs2"] = "precomputed"
    # Fold-level leakage control (Stage 5). Sequence clustering alone misses
    # structural redundancy: two chains under the identity threshold can still be
    # the same fold, which an inverse-folding model (structure -> sequence) would
    # leak across splits. When set, protein entities sharing a structural
    # (super)family are union-merged into the same component in ADDITION to shared
    # sequence clusters — so the same fold cannot straddle train/test. Metadata
    # only (RCSB's precomputed classifications; no coordinates). "cath" keys on the
    # homologous-superfamily code (e.g. 1.10.490.10); "ecod"/"scop2" key on the
    # (super)family name. "off" = sequence-only (prior behavior). Purely additive:
    # it can only merge components, never split them, and chains lacking the chosen
    # classification simply add no structural edge. Off by default: fold-merging
    # the dominant superfamilies (antibodies, TIM barrels) collapses them into
    # mega-components that land wholesale in one split, skewing the ENTRY-level
    # train/val/test balance (~95/3/2 at superfamily grain) even though the
    # COMPONENT-level split stays ~80/10/10. Off by default; pair it with
    # split_strategy="balanced" (below) to restore entry balance — that is the
    # "masterclass" recipe (structural_clustering="scop2" + balanced).
    structural_clustering: Literal["off", "cath", "ecod", "scop2"] = "off"
    split_fractions: SplitFractions
    # Component -> split assignment strategy.
    #   "hash": each component hashed onto the cumulative fractions (balances
    #     COMPONENTS). Simple and registry-free-stable, but heavy-tailed component
    #     sizes skew ENTRY balance (one dominant fold => that split balloons).
    #   "balanced": cap the dominant folds to train and fill val/test to their
    #     ENTRY targets from the tail of smaller folds (hash-ordered). Restores
    #     ~80/10/10 by entries and yields diverse, fold-honest val/test sets. Best
    #     paired with structural_clustering (esp. "scop2"); also fixes the plain
    #     sequence-only skew from the antibody mega-cluster. Growth stability comes
    #     from splits.registry.json pinning prior assignments (like test minimums).
    split_strategy: Literal["hash", "balanced"] = "hash"
    split_salt: str = Field(min_length=1)
    seed: int = Field(ge=0)

    # --- test-set minimums (opt-in stratification top-up) ---
    # Floor on the number of test entries carrying each functional ligand class,
    # e.g. {"metal": 500, "nucleic_acid": 200}. Empty (default) = pure deterministic
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

    # --- shareable-spec metadata (descriptive only; EXCLUDED from config_hash) ---
    # A self-identifying header + optional human metadata, so a config.yaml doubles
    # as a portable "split spec" you can hand someone to reproduce your methodology.
    spec: SpecMeta | None = None

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

    @model_validator(mode="after")
    def _identity_level_supported(self) -> Config:
        """Reject identity thresholds the precomputed backend can't actually use.

        The precomputed backend looks up ``identity_level`` in each entity's RCSB
        ``rcsb_cluster_membership``, which exists only at RCSB_IDENTITY_LEVELS. An
        unsupported level (e.g. 0.40 -> 40%) matches nothing, so every entity would
        fall into its own singleton cluster: no clustering, and cross-split sequence
        leakage that ``check_no_leakage`` cannot detect. Fail loudly instead. The
        mmseqs2 backend clusters at an arbitrary threshold, so it is exempt.
        """
        if (
            self.clustering_backend == "precomputed"
            and self.identity_level not in RCSB_IDENTITY_LEVELS
        ):
            levels = ", ".join(str(x) for x in sorted(RCSB_IDENTITY_LEVELS))
            raise ValueError(
                f"identity_threshold={self.identity_threshold} -> {self.identity_level}% is not an "
                f"RCSB precomputed cluster level. For clustering_backend='precomputed' use one of "
                f"{levels}% (i.e. 0.30/0.50/0.70/0.90/0.95/1.00); otherwise no clustering happens."
            )
        return self

    @property
    def dataset_version(self) -> str:
        """Versioned dataset name, e.g. 'IF-Split-2026.05.30'."""
        return f"IF-Split-{self.snapshot_date:%Y.%m.%d}"

    @property
    def identity_level(self) -> int:
        """``identity_threshold`` as an integer percent (e.g. 0.30 -> 30)."""
        return round(self.identity_threshold * 100)

    def canonical_dict(self) -> dict[str, Any]:
        """JSON-mode dump of the output-affecting settings (for hashing/manifests).

        The ``spec`` metadata is descriptive, not output-affecting, so it is
        excluded here — two splits identical but for name/author hash the same.
        """
        d = self.model_dump(mode="json")
        d.pop("spec", None)
        return d

    def config_hash(self) -> str:
        """Deterministic, formatting-independent hash of the settings."""
        canonical = json.dumps(self.canonical_dict(), sort_keys=True, separators=(",", ":"))
        return hashlib.blake2b(canonical.encode("utf-8"), digest_size=16).hexdigest()

    def to_spec_dict(self, *, stamp_hash: bool = True) -> dict[str, Any]:
        """A portable split-spec mapping: a self-identifying header + all settings.

        Suitable for ``yaml.safe_dump`` and re-loading via :func:`load_config`. The
        ``spec`` header carries the schema id and (optionally) the expected
        config-hash so the shared file self-verifies on reload.
        """
        from . import __version__

        meta = (self.spec or SpecMeta()).model_dump(mode="json", exclude_none=True)
        meta["ifsplit_spec"] = SPEC_SCHEMA
        meta.setdefault("created_with", f"if-split {__version__}")
        if stamp_hash:
            meta["expected_config_hash"] = self.config_hash()
        body = self.canonical_dict()  # settings only (no spec)
        return {"spec": meta, **body}


def load_config(path: str | Path) -> Config:
    """Load and validate a YAML config file into a :class:`Config`."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config not found: {path}")
    with path.open("r", encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    if not isinstance(raw, dict):
        raise ValueError(f"Config must be a YAML mapping, got {type(raw).__name__}: {path}")
    cfg = Config.model_validate(raw)
    # If a shared spec stamped an expected hash, verify we reproduce it. A mismatch
    # means the settings were edited after stamping (or a tool-version difference) —
    # warn, don't fail, so the file stays usable.
    expected = cfg.spec.expected_config_hash if cfg.spec else None
    if expected and expected != cfg.config_hash():
        import warnings

        warnings.warn(
            f"spec.expected_config_hash {expected} != computed {cfg.config_hash()}: "
            f"the settings in {path} differ from when the spec was stamped.",
            stacklevel=2,
        )
    return cfg
