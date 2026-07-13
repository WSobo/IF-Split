"""Stage 8 - Thin loader over a manifest, with cluster-balanced sampling.

Reads ``manifest.json`` and exposes train/val/test views. Coordinate/featurization
loading is intentionally *not* here - it's the optional downstream concern
(PLAN.md §1.5); a model repo plugs its own featurizer onto these entry lists.

Two training views over the *same* leakage-safe split, for the two things an
inverse-folding model actually consumes:

- **Backbones** (``entry_ids`` / ``backbones``): every kept structure, ligand-
  agnostic - the ProteinMPNN-style scale lever. A structure with only junk ions or
  no ligand is still a good backbone.
- **Conditioning targets** (``conditioning_targets``): the ``functional``-tier
  ligands to condition on (LigandMPNN-style), one per (structure, ligand). Junk is
  never a target. ``metal_site_nonnative`` pockets are opt-in via ``include_nonnative``.

The PDB is heavily redundant (thousands of near-identical lysozyme / kinase
co-crystals). ``sample_by_cluster`` draws one entry per sequence cluster per epoch
so over-represented folds don't dominate - deterministic given a seed (no global
RNG), so an epoch is reproducible.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

from .manifest import (
    TARGETS_FILENAME,
    read_classes,
    read_clusters,
    read_id_list,
    read_manifest,
    read_targets,
)

SPLITS = ("train", "val", "test")


@dataclass(frozen=True)
class ConditioningTarget:
    """A ligand for a ligand-conditioned model to condition on."""

    entry_id: str
    split: str
    cluster: str
    ligand_class: str  # metal | small_molecule | nucleic_acid
    comp_id: str | None  # CCD id (None for a nucleic_acid complex)
    tier: str  # functional | ambiguous (ambiguous only for opt-in non-native sites)
    reason: str


def _stable_rank(key: str, seed: int) -> int:
    """Deterministic pseudo-random rank for ``key`` under ``seed`` (no global RNG)."""
    digest = hashlib.blake2b(f"{seed}:{key}".encode(), digest_size=8).digest()
    return int.from_bytes(digest, "big")


@dataclass
class SplitView:
    name: str
    entry_ids: list[str]
    ligand_classes: dict[str, list[str]]  # entry_id -> classes
    entry_clusters: dict[str, str]  # entry_id -> cluster key
    targets: list[ConditioningTarget] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.entry_ids)

    def with_class(self, cls: str) -> list[str]:
        """Entry ids in this split tagged with ligand class ``cls``."""
        return [e for e in self.entry_ids if cls in self.ligand_classes.get(e, [])]

    # ----------------------------- training views -------------------------- #
    @property
    def backbones(self) -> list[str]:
        """All entry ids in this split (backbone-only training corpus)."""
        return list(self.entry_ids)

    def conditioning_targets(
        self, classes: list[str] | None = None, *, include_nonnative: bool = False
    ) -> list[ConditioningTarget]:
        """Ligand-conditioning targets in this split (LigandMPNN-style).

        ``functional``-tier only by default; ``include_nonnative`` also returns the
        opt-in ``metal_site_nonnative`` pockets (real sites whose native metal isn't
        Ni/Co). ``classes`` filters to given ligand classes (metal / small_molecule /
        nucleic_acid). One entry per (structure, ligand) - group by ``entry_id`` to
        condition on all of a structure's ligands at once.
        """
        out = []
        for t in self.targets:
            if t.tier != "functional" and not include_nonnative:
                continue
            if classes is not None and t.ligand_class not in classes:
                continue
            out.append(t)
        return out

    def targets_by_entry(
        self, classes: list[str] | None = None, *, include_nonnative: bool = False
    ) -> dict[str, list[ConditioningTarget]]:
        """Conditioning targets grouped by structure (condition-on-all view)."""
        groups: dict[str, list[ConditioningTarget]] = {}
        for t in self.conditioning_targets(classes, include_nonnative=include_nonnative):
            groups.setdefault(t.entry_id, []).append(t)
        return dict(sorted(groups.items()))

    def conditioned_entry_ids(
        self, classes: list[str] | None = None, *, include_nonnative: bool = False
    ) -> list[str]:
        """Sorted entry ids that carry >=1 conditioning target (ligand-conditioned set)."""
        return sorted(self.targets_by_entry(classes, include_nonnative=include_nonnative))

    def cluster_groups(self) -> dict[str, list[str]]:
        """Map cluster key -> sorted entry ids within this split."""
        groups: dict[str, list[str]] = {}
        for e in self.entry_ids:
            key = self.entry_clusters.get(e, e)
            groups.setdefault(key, []).append(e)
        return {k: sorted(v) for k, v in sorted(groups.items())}

    def sample_by_cluster(self, seed: int = 0) -> list[str]:
        """One entry per cluster, chosen deterministically by ``seed``.

        De-redundifies the split: each sequence cluster contributes exactly one
        representative, so over-represented folds don't dominate. Vary ``seed``
        across epochs to rotate which member of each cluster is drawn. Returns a
        deterministically ordered list (sorted by the same stable rank).
        """
        chosen: list[tuple[int, str]] = []
        for key, members in self.cluster_groups().items():
            rep = min(members, key=lambda e: (_stable_rank(e, seed), e))
            chosen.append((_stable_rank(key, seed), rep))
        return [e for _, e in sorted(chosen)]


class IFSplitDataset:
    """Read-only view over a built manifest's train/val/test partition."""

    def __init__(self, manifest_path: str | Path) -> None:
        self._m = read_manifest(manifest_path)
        self._dir = Path(manifest_path).parent
        self.dataset_version: str = self._m["dataset_version"]
        self.config_hash: str = self._m["config_hash"]
        files = self._m.get("files", {})
        self._split_files: dict[str, str] = files.get("splits", {})
        # Supporting maps live in sidecar files referenced by the manifest.
        self._classes = read_classes(
            self._dir / files.get("ligand_classes", "ligands.classes.json")
        )
        self._entry_clusters = read_clusters(self._dir / files.get("clusters", "clusters.json"))
        # Conditioning-target corpus, grouped by split (absent -> empty, older builds).
        self._targets_by_split: dict[str, list[ConditioningTarget]] = {s: [] for s in SPLITS}
        for row in read_targets(self._dir / files.get("targets", TARGETS_FILENAME)):
            split = row.get("split")
            if split in self._targets_by_split:
                self._targets_by_split[split].append(
                    ConditioningTarget(
                        entry_id=row["entry_id"],
                        split=split,
                        cluster=row.get("cluster", row["entry_id"]),
                        ligand_class=row["class"],
                        comp_id=row.get("comp_id"),
                        tier=row.get("tier", "functional"),
                        reason=row.get("reason", ""),
                    )
                )

    def split(self, name: str) -> SplitView:
        if name not in SPLITS:
            raise KeyError(f"unknown split {name!r}; expected one of {SPLITS}")
        fname = self._split_files.get(name, f"{name}.json")
        ids = read_id_list(self._dir / fname)
        return SplitView(
            name=name,
            entry_ids=ids,
            ligand_classes={e: self._classes.get(e, []) for e in ids},
            entry_clusters={e: self._entry_clusters.get(e, e) for e in ids},
            targets=self._targets_by_split.get(name, []),
        )

    @property
    def train(self) -> SplitView:
        return self.split("train")

    @property
    def val(self) -> SplitView:
        return self.split("val")

    @property
    def test(self) -> SplitView:
        return self.split("test")


def load_dataset(manifest_path: str | Path) -> IFSplitDataset:
    return IFSplitDataset(manifest_path)
