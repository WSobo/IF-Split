"""Stage 8 - Thin loader / torch Dataset over a manifest.

Reads a manifest, lazily loads parsed structures + (optional) ligand context,
and exposes train/val/test views. Featurization is pluggable so a model repo can
supply its own. Lands in Phase 6.
"""

from __future__ import annotations


def load_dataset(manifest_path, *args, **kwargs):
    raise NotImplementedError("Stage 8 (loader) lands in Phase 6.")
