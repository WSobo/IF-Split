"""Phase 1 tests: config loading, validation, and deterministic hashing."""

from __future__ import annotations

from pathlib import Path

import pytest
from pydantic import ValidationError

from ifsplit.config import Config, load_config

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config" / "default.yaml"


def _good_dict() -> dict:
    return {
        "snapshot_date": "2026-05-30",
        "experimental_methods": ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"],
        "resolution_max_A": 3.5,
        "max_total_residues": 5999,
        "excluded_het": ["HOH", "NA", "CL", "K", "BR"],
        "use_biological_assembly": True,
        "identity_threshold": 0.30,
        "clustering_backend": "precomputed",
        "split_fractions": {"train": 0.80, "val": 0.10, "test": 0.10},
        "split_salt": "snapsplit-v1",
        "seed": 0,
        "ligand_context_radius_A": 8.0,
        "max_ligand_atoms": 25,
    }


def test_default_config_loads():
    cfg = load_config(DEFAULT_CONFIG)
    assert isinstance(cfg, Config)
    assert cfg.identity_threshold == 0.30
    assert cfg.use_biological_assembly is True


def test_dataset_version():
    cfg = load_config(DEFAULT_CONFIG)
    assert cfg.dataset_version == "IF-Split-2026.05.30"


def test_methods_are_normalized_uppercase():
    d = _good_dict()
    d["experimental_methods"] = ["x-ray diffraction"]
    cfg = Config.model_validate(d)
    assert cfg.experimental_methods == ["X-RAY DIFFRACTION"]


def test_config_hash_is_deterministic():
    cfg = load_config(DEFAULT_CONFIG)
    assert cfg.config_hash() == cfg.config_hash()
    assert len(cfg.config_hash()) == 32  # blake2b digest_size=16 -> 32 hex chars


def test_config_hash_is_formatting_independent():
    # Same values, different declaration order -> identical hash.
    a = Config.model_validate(_good_dict())
    d = _good_dict()
    d = {k: d[k] for k in reversed(list(d))}
    b = Config.model_validate(d)
    assert a.config_hash() == b.config_hash()


def test_config_hash_changes_with_salt():
    base = Config.model_validate(_good_dict())
    d = _good_dict()
    d["split_salt"] = "snapsplit-v2"
    assert Config.model_validate(d).config_hash() != base.config_hash()


def test_bad_split_fractions_rejected():
    d = _good_dict()
    d["split_fractions"] = {"train": 0.80, "val": 0.10, "test": 0.20}
    with pytest.raises(ValidationError):
        Config.model_validate(d)


def test_identity_threshold_out_of_range_rejected():
    d = _good_dict()
    d["identity_threshold"] = 1.5
    with pytest.raises(ValidationError):
        Config.model_validate(d)


def test_clustering_backend_default_is_precomputed():
    d = _good_dict()
    del d["clustering_backend"]
    assert Config.model_validate(d).clustering_backend == "precomputed"


def test_unknown_clustering_backend_rejected():
    d = _good_dict()
    d["clustering_backend"] = "blastclust"
    with pytest.raises(ValidationError):
        Config.model_validate(d)


def test_unknown_key_rejected():
    d = _good_dict()
    d["totally_made_up"] = True
    with pytest.raises(ValidationError):
        Config.model_validate(d)
