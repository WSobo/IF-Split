"""Stage 8 - Thin loader over a manifest, with cluster-balanced sampling.

Reads ``manifest.json`` and exposes train/val/test views (entry ids, ligand
classes, and the entry->cluster map). Coordinate/featurization loading is
intentionally *not* here - it's the optional downstream concern (PLAN.md §1.5);
a model repo plugs its own featurizer onto these entry lists.

The PDB is heavily redundant (thousands of near-identical lysozyme / kinase
co-crystals). Training by sampling entries uniformly drowns the model in
over-represented folds. ``sample_by_cluster`` draws one entry per sequence
cluster per epoch, which is the bigger training-quality lever than ligand tiering
- and it is free here because the clusters are already computed. Sampling is
deterministic given a seed (no global RNG), so an epoch is reproducible.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

from .manifest import read_classes, read_clusters, read_id_list, read_manifest

SPLITS = ("train", "val", "test")


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

    def __len__(self) -> int:
        return len(self.entry_ids)

    def with_class(self, cls: str) -> list[str]:
        """Entry ids in this split tagged with ligand class ``cls``."""
        return [e for e in self.entry_ids if cls in self.ligand_classes.get(e, [])]

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
