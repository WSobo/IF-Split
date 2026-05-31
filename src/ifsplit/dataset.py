"""Stage 8 - Thin loader over a manifest.

Reads ``manifest.json`` and exposes train/val/test views (entry ids + ligand
classes). Coordinate/featurization loading is intentionally *not* here - it's the
optional downstream concern (PLAN.md §1.5); a model repo plugs its own
featurizer onto these entry lists.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .manifest import read_manifest

SPLITS = ("train", "val", "test")


@dataclass
class SplitView:
    name: str
    entry_ids: list[str]
    ligand_classes: dict[str, list[str]]  # entry_id -> classes

    def __len__(self) -> int:
        return len(self.entry_ids)

    def with_class(self, cls: str) -> list[str]:
        """Entry ids in this split tagged with ligand class ``cls``."""
        return [e for e in self.entry_ids if cls in self.ligand_classes.get(e, [])]


class IFSplitDataset:
    """Read-only view over a built manifest's train/val/test partition."""

    def __init__(self, manifest_path: str | Path) -> None:
        self._m = read_manifest(manifest_path)
        self.dataset_version: str = self._m["dataset_version"]
        self.config_hash: str = self._m["config_hash"]
        self._classes: dict[str, list[str]] = self._m["ligands"]["classes"]
        self._entries: dict[str, list[str]] = self._m["splits"]["entries"]

    def split(self, name: str) -> SplitView:
        if name not in SPLITS:
            raise KeyError(f"unknown split {name!r}; expected one of {SPLITS}")
        ids = self._entries.get(name, [])
        return SplitView(
            name=name,
            entry_ids=ids,
            ligand_classes={e: self._classes.get(e, []) for e in ids},
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
