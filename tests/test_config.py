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


def test_precomputed_identity_threshold_must_be_an_rcsb_level():
    # 0.40 -> 40% is a natural-looking value but is NOT an RCSB precomputed level;
    # accepting it would silently produce all-singleton clusters (no clustering ->
    # cross-split leakage). It must be rejected for the precomputed backend.
    d = _good_dict()
    d["identity_threshold"] = 0.40
    with pytest.raises(ValidationError):
        Config.model_validate(d)


def test_all_rcsb_identity_levels_accepted_for_precomputed():
    for frac in (0.30, 0.50, 0.70, 0.90, 0.95, 1.00):
        d = _good_dict()
        d["identity_threshold"] = frac
        assert Config.model_validate(d).identity_level == round(frac * 100)


def test_nonstandard_identity_level_allowed_for_mmseqs2_backend():
    # mmseqs2 clusters at an arbitrary threshold, so a non-RCSB level is fine there
    # even though the precomputed backend rejects it.
    d = _good_dict()
    d["identity_threshold"] = 0.40
    d["clustering_backend"] = "mmseqs2"
    assert Config.model_validate(d).identity_level == 40


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


# ------------------------------- spec sharing ------------------------------ #
def test_spec_metadata_excluded_from_hash():
    # Two configs identical but for spec metadata must hash the same.
    base = Config.model_validate(_good_dict())
    d = _good_dict()
    d["spec"] = {"name": "my-split", "author": "WSobo", "description": "demo"}
    withmeta = Config.model_validate(d)
    assert withmeta.spec is not None
    assert withmeta.config_hash() == base.config_hash()


def test_to_spec_dict_roundtrips_to_same_hash():
    cfg = load_config(DEFAULT_CONFIG)
    doc = cfg.to_spec_dict(stamp_hash=True)
    assert doc["spec"]["ifsplit_spec"] == "ifsplit/config@1"
    assert doc["spec"]["expected_config_hash"] == cfg.config_hash()
    # Reloading the emitted spec yields the same output-affecting settings.
    back = Config.model_validate(doc)
    assert back.config_hash() == cfg.config_hash()


def test_spec_hash_mismatch_warns(tmp_path):
    import warnings

    import yaml

    cfg = load_config(DEFAULT_CONFIG)
    doc = cfg.to_spec_dict(stamp_hash=True)
    doc["resolution_max_A"] = 2.0  # change a setting AFTER stamping the hash
    p = tmp_path / "tampered.ifsplit.yaml"
    p.write_text(yaml.safe_dump(doc), encoding="utf-8")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        load_config(p)
    assert any("expected_config_hash" in str(x.message) for x in w)


def test_spec_matching_hash_no_warn(tmp_path):
    import warnings

    import yaml

    cfg = load_config(DEFAULT_CONFIG)
    p = tmp_path / "clean.ifsplit.yaml"
    p.write_text(yaml.safe_dump(cfg.to_spec_dict(stamp_hash=True)), encoding="utf-8")
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        load_config(p)
    assert not any("expected_config_hash" in str(x.message) for x in w)
