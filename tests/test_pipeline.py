"""Phases 3/5/6/7 tests: filter, cluster, split, manifest, determinism (offline)."""

from __future__ import annotations

import collections
from pathlib import Path

import pytest

from ifsplit.cluster import build_clusters
from ifsplit.config import Config, load_config
from ifsplit.dataset import load_dataset
from ifsplit.ligands import classify_components
from ifsplit.manifest import build_manifest, summarize_manifest, write_manifest
from ifsplit.parse import (
    DROP_CLASHSCORE,
    DROP_EM_INCLUSION,
    DROP_NO_PROTEIN,
    DROP_NO_SEQUENCE,
    DROP_NO_VALIDATION,
    DROP_RESOLUTION,
    DROP_RFREE,
    DROP_SEQUENCE_TOO_SHORT,
    DROP_TOO_LARGE,
    drop_summary,
    filter_candidates,
)
from ifsplit.schema import CandidateRecord, PolymerEntity
from ifsplit.split import assign_splits, bucket, check_no_leakage, split_for_key

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = ROOT / "config" / "default.yaml"


def _cfg(**over) -> Config:
    cfg = load_config(DEFAULT_CONFIG)
    return cfg.model_copy(update=over) if over else cfg


def _records(sample_entries, artifact_entry) -> list[CandidateRecord]:
    recs = [CandidateRecord.from_data_api(e) for e in sample_entries.values()]
    recs.append(CandidateRecord.from_data_api(artifact_entry))
    return recs


def _seq_record(entry_id: str, seq: str) -> CandidateRecord:
    """A single-protein-chain record carrying the given canonical sequence."""
    return CandidateRecord(
        entry_id=entry_id,
        methods=["X-RAY DIFFRACTION"],
        resolution_A=2.0,
        release_date="2020-01-01",
        deposited_residues=len(seq),
        assemblies={f"{entry_id}-1": len(seq)},
        polymer_entities=[
            PolymerEntity(
                entity_id=f"{entry_id}_1",
                polymer_type="Protein",
                seq_len=len(seq),
                seq=seq,
                cluster_ids={30: 1},
            )
        ],
        nonpolymer_comps=[],
        bound_components=[],
        affinity_comp_ids=[],
    )


def _unclustered_record(entry_id: str, seq: str) -> CandidateRecord:
    """A single-protein record whose chain RCSB left UNCLUSTERED at 30% (cluster_ids empty)."""
    rec = _seq_record(entry_id, seq)
    rec.polymer_entities[0].cluster_ids = {}
    return rec


def test_identical_unclustered_sequences_share_one_component():
    # Two distinct entries carrying the SAME unclustered peptide. Keyed by sequence
    # (not entity id) they must collapse into ONE component, so the sequence cannot
    # straddle splits and check_no_leakage stays clean under any salt. Entity-id
    # keying (the bug) gave two components that a salt could split apart.
    seq = "GSHMWYPQR" * 2  # short peptide RCSB would not cluster at 30%
    clusters = build_clusters(
        [_unclustered_record("1PEP", seq), _unclustered_record("2PEP", seq)], _cfg()
    )
    assert clusters.n_clusters == 1
    assert clusters.entry_to_cluster["1PEP"] == clusters.entry_to_cluster["2PEP"]
    for salt in ("s1", "s2", "s3", "s4", "s5"):
        splits = assign_splits(clusters, _cfg(split_salt=salt))
        check_no_leakage(splits, clusters)
        assert splits.entry_split["1PEP"] == splits.entry_split["2PEP"]


def _clustered_record(entry_id: str, cid: int, seq: str) -> CandidateRecord:
    """An entry whose single protein chain RCSB DID cluster, under cluster id `cid`."""
    return CandidateRecord(
        entry_id=entry_id,
        methods=["X-RAY DIFFRACTION"],
        resolution_A=2.0,
        release_date="2020-01-01",
        deposited_residues=len(seq),
        assemblies={f"{entry_id}-1": len(seq)},
        polymer_entities=[
            PolymerEntity(
                entity_id=f"{entry_id}_1",
                polymer_type="Protein",
                seq_len=len(seq),
                seq=seq,
                cluster_ids={30: cid},
            )
        ],
        nonpolymer_comps=[],
        bound_components=[],
        affinity_comp_ids=[],
    )


def test_identical_sequences_merge_despite_different_cluster_ids():
    # RCSB's cluster file is not identity-complete: byte-identical sequences can carry
    # DIFFERENT 30% cluster ids (measured on the 2026-07-22 snapshot: 69 sequences across
    # 497 entries, incl. a 621-residue chain in test AND val). Keying a clustered chain by
    # its cluster id alone (the bug) let the SAME protein straddle two splits. Exact
    # sequence identity must merge regardless of what the cluster file says.
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGHKLMNPQRSTVWY"  # 41 aa, well above the modeled gate
    clusters = build_clusters(
        [_clustered_record("1AAA", 111, seq), _clustered_record("2AAA", 222, seq)], _cfg()
    )
    assert clusters.n_clusters == 1
    assert clusters.entry_to_cluster["1AAA"] == clusters.entry_to_cluster["2AAA"]
    for salt in ("s1", "s2", "s3", "s4", "s5"):
        splits = assign_splits(clusters, _cfg(split_salt=salt))
        check_no_leakage(splits, clusters)
        assert splits.entry_split["1AAA"] == splits.entry_split["2AAA"]


def test_identity_keying_does_not_inflate_multichain_count():
    # The identity key is added to every clustered chain, so a SINGLE-chain entry now
    # carries two keys. `multichain` counts entries that BRIDGE clusters, so it must be
    # computed before identity keys fold in — otherwise every entry looks bridging.
    clusters = build_clusters(
        [
            _clustered_record("1AAA", 111, "M" + "A" * 40),
            _clustered_record("2AAA", 222, "M" + "C" * 40),
        ],
        _cfg(),
    )
    assert clusters.multichain_entries == []


def _domain_record(entry_id: str, cid: int, domains: dict[str, list[str]]) -> CandidateRecord:
    """Single-chain record in raw cluster ``cid`` carrying Pfam/InterPro ``domains``."""
    rec = _clustered_record(entry_id, cid, _seq_for(cid))
    rec.polymer_entities[0].domain_families = domains
    return rec


def test_domain_families_merge_under_pfam_and_all_but_not_structural():
    # Pfam/InterPro are HMM families over SEQUENCE, so they carry no classification lag
    # and catch homology CATH/ECOD/SCOP2 have not annotated yet. Two entries in different
    # sequence clusters sharing a Pfam family must merge under "pfam" and under "all",
    # and must NOT merge under a structural-only method that cannot see them.
    recs = [
        _domain_record("1AAA", 11, {"pfam": ["PF00042"]}),
        _domain_record("2AAA", 22, {"pfam": ["PF00042"]}),
    ]
    assert build_clusters(recs, _cfg(structural_clustering="pfam")).n_clusters == 1
    assert build_clusters(recs, _cfg(structural_clustering="all")).n_clusters == 1
    assert build_clusters(recs, _cfg(structural_clustering="scop2")).n_clusters == 2
    assert build_clusters(recs, _cfg(structural_clustering="union")).n_clusters == 2
    assert build_clusters(recs, _cfg(structural_clustering="off")).n_clusters == 2


def test_all_merges_on_any_authority_and_namespaces_keys():
    # "all" must merge on a structural OR a domain family, and must not confuse the two:
    # a CATH code and a Pfam accession that happen to collide as strings stay distinct.
    struct = [
        _fold_record("AAA1", 10, {"cath": ["1.10.490.10"]}),
        _fold_record("BBB2", 20, {"cath": ["1.10.490.10"]}),
    ]
    assert build_clusters(struct, _cfg(structural_clustering="all")).n_clusters == 1
    collide = [
        _domain_record("1AAA", 11, {"pfam": ["X1"]}),
        _fold_record("BBB2", 22, {"cath": ["X1"]}),
    ]
    # same raw string, different authority -> namespaced -> no merge
    assert build_clusters(collide, _cfg(structural_clustering="all")).n_clusters == 2


def _many_records(n_giant: int, n_tail: int) -> list[CandidateRecord]:
    """A dominant component of ``n_giant`` entries plus ``n_tail`` singleton components."""
    recs = [_clustered_record(f"G{i:03d}", 1, _seq_for(1)) for i in range(n_giant)]
    recs += [_clustered_record(f"T{i:03d}", 100 + i, _seq_for(100 + i)) for i in range(n_tail)]
    return recs


def test_maximal_holds_out_the_whole_tail_and_never_starves_val():
    # The case `balanced` cannot serve: the tail is far thinner than a 10/10 target, so
    # balanced fills test first and leaves val EMPTY. `maximal` treats the fractions as a
    # ceiling and splits whatever tail exists evenly, so both held-out sets are non-empty.
    recs = _many_records(n_giant=200, n_tail=20)
    clusters = build_clusters(recs, _cfg())
    bal = assign_splits(clusters, _cfg(split_strategy="balanced"))
    mx = assign_splits(clusters, _cfg(split_strategy="maximal"))
    counts = lambda r: collections.Counter(r.entry_split.values())  # noqa: E731
    assert counts(bal)["val"] == 0  # the starvation this strategy exists to fix
    cm = counts(mx)
    assert cm["val"] > 0 and cm["test"] > 0
    assert cm["val"] + cm["test"] == 20  # the entire tail is held out
    assert cm["train"] == 200  # the dominant component stays in train
    assert abs(cm["val"] - cm["test"]) <= 1  # filled toward the smaller side
    check_no_leakage(mx, clusters)


def test_maximal_respects_the_fraction_ceiling_when_the_tail_is_large():
    # With no merging the tail can be most of the snapshot; the ceiling stops it from
    # starving train. 100 singleton components, ceiling val+test = 20% -> 20 entries.
    recs = _many_records(n_giant=0, n_tail=100)
    clusters = build_clusters(recs, _cfg())
    mx = assign_splits(clusters, _cfg(split_strategy="maximal"))
    cm = collections.Counter(mx.entry_split.values())
    assert cm["val"] + cm["test"] <= 20
    assert cm["train"] >= 80
    check_no_leakage(mx, clusters)


def test_maximal_reports_no_gap_when_the_tail_is_thin():
    # A thin tail is an OUTCOME for `maximal`, not a missed target, so it must not be
    # reported as a shortfall the way `balanced` does.
    clusters = build_clusters(_many_records(n_giant=200, n_tail=6), _cfg())
    assert assign_splits(clusters, _cfg(split_strategy="balanced")).balance_gaps
    assert assign_splits(clusters, _cfg(split_strategy="maximal")).balance_gaps == {}


def test_maximal_never_holds_out_the_giant_when_a_pin_says_it_should():
    # The bug this guards: under `maximal` the dominant component absorbs held-out
    # components as the snapshot grows, and honoring its INHERITED pin unconditionally
    # put the giant in the holdout. Measured on the real 2026-07-22 snapshot with a
    # 2023-cutoff registry: train=1,164 / val=862 / test=214,796, reported as
    # growth_stable=true with pinned_reassignments=1. The cap must beat the pin.
    cfg = _cfg(split_strategy="maximal")
    v1_recs = _many_records(n_giant=200, n_tail=20)
    v1_clusters = build_clusters(v1_recs, cfg)
    v1 = assign_splits(v1_clusters, cfg)
    giant = v1.entry_split["G000"]
    assert giant == "train"  # precondition: the giant starts in train
    held = [k for k, s in v1.cluster_split.items() if s in ("val", "test")]
    assert held, "precondition: the v1 tail is held out"

    # Growth: a bridging entry welds a held-out tail component onto the giant, so the
    # merged component inherits that component's val/test pin. Pin the merged component
    # directly rather than replaying v1's registry: a raw cluster key is named after its
    # smallest member entity, so the bridge RENAMES the giant's key and the replayed pin
    # would not match at all (a separate, documented weakness). What is under test here is
    # what happens once the giant does carry a held-out pin, which is what the real
    # 2023-registry rebuild produced.
    bridge = _protein_record("BRDG", [1, 100])
    v2_clusters = build_clusters([*v1_recs, bridge], cfg)
    assert v2_clusters.entry_to_cluster["G000"] == v2_clusters.entry_to_cluster["T000"]
    giant_key = v2_clusters.entry_to_cluster["G000"]

    v2 = assign_splits(v2_clusters, cfg, registry={giant_key: "test"})
    check_no_leakage(v2, v2_clusters)
    counts = collections.Counter(v2.entry_split.values())
    # The invariant, at the ENTRY level: the giant is in train and the holdout stays a
    # tail. Asserting on cluster_split would pass even on the inverted split.
    assert v2.entry_split["G000"] == "train"
    assert v2.entry_split["T000"] == "train"  # absorbed into the giant, follows it
    assert counts["train"] > counts["val"] + counts["test"]
    assert counts["val"] + counts["test"] <= 0.2 * sum(counts.values())  # ceiling holds
    # And the override is reported in ENTRIES, not just as a component count of 1.
    assert v2.pinned_reassignments == 1
    assert v2.pinned_entries_reassigned >= 200  # the component-level "1" hides this


def test_distinct_unclustered_sequences_stay_separate():
    clusters = build_clusters(
        [_unclustered_record("1PEP", "GSHMWYPQRT"), _unclustered_record("2PEP", "AAAKKKDDEE")],
        _cfg(),
    )
    assert clusters.n_clusters == 2  # different sequences -> different components


_AA = "ACDEFGHIKLMNPQRSTVWY"


def _seq_for(cid: int, length: int = 100) -> str:
    """A deterministic dummy protein sequence unique to raw cluster ``cid``.

    Distinct clusters must get distinct sequences: an exact-sequence identity edge merges
    byte-identical chains regardless of their cluster id, so a single shared constant
    sequence would collapse clusters a fixture means to keep apart (and would not resemble
    real data, where different clusters never share a byte-identical chain).
    """
    tag = "".join(_AA[(cid // (20**i)) % 20] for i in range(4))
    return (tag + _AA)[:length].ljust(length, "G")


def _mixed_record(entry_id: str, cid: int, uncl_seq: str) -> CandidateRecord:
    """An entry with one CLUSTERED protein chain (cluster `cid`) + one UNCLUSTERED chain."""
    return CandidateRecord(
        entry_id=entry_id,
        methods=["X-RAY DIFFRACTION"],
        resolution_A=2.0,
        release_date="2020-01-01",
        deposited_residues=100,
        assemblies={f"{entry_id}-1": 100},
        polymer_entities=[
            PolymerEntity(
                entity_id=f"{entry_id}_1",
                polymer_type="Protein",
                seq_len=60,
                seq=_seq_for(cid, 60),
                cluster_ids={30: cid},
            ),
            PolymerEntity(
                entity_id=f"{entry_id}_2",
                polymer_type="Protein",
                seq_len=len(uncl_seq),
                seq=uncl_seq,
                cluster_ids={},
            ),
        ],
        nonpolymer_comps=[],
        bound_components=[],
        affinity_comp_ids=[],
    )


_LONG_UNCL = "MKTAYIAKQRQISFVKSHFSRQLEERLGHKLMNPQRSTVWY"  # 41 aa >= MIN_UNCLUSTERED_MERGE_LEN
_SHORT_PEP = "GSHMWYPQR"  # 9 aa, below the merge gate


def test_partial_entries_merge_via_long_unclustered_chain():
    # Two entries with DISTINCT clustered chains but the SAME long unclustered chain must
    # merge — a long shared sequence is a real leak signal, keeping them apart would leak it.
    cr = build_clusters(
        [_mixed_record("EE01", 1, _LONG_UNCL), _mixed_record("EE02", 2, _LONG_UNCL)], _cfg()
    )
    assert cr.entry_to_cluster["EE01"] == cr.entry_to_cluster["EE02"]


def test_short_shared_peptide_does_not_merge_partial_entries():
    # The SAME short peptide as a second chain must NOT union two otherwise-unrelated
    # proteins (a peptide ligand/tag is not a homology signal, and would fan out).
    cr = build_clusters(
        [_mixed_record("EE01", 1, _SHORT_PEP), _mixed_record("EE02", 2, _SHORT_PEP)], _cfg()
    )
    assert cr.entry_to_cluster["EE01"] != cr.entry_to_cluster["EE02"]


def test_promiscuous_short_peptide_forms_no_megacomponent():
    # 50 unrelated proteins each carry the SAME short peptide as a second chain. The gate
    # keeps them 50 components; WITHOUT it they collapse into ONE mega-component (== every
    # entry), which under `hash` would land wholesale in a salt-chosen split. Goes red if
    # someone removes the gate.
    recs = [_mixed_record(f"P{i:03d}", 100 + i, _SHORT_PEP) for i in range(50)]
    cr = build_clusters(recs, _cfg())
    n = sum(len(m) for m in cr.cluster_members.values())
    biggest = max(len(m) for m in cr.cluster_members.values())
    assert biggest <= 0.1 * n  # gate holds (biggest == 1); without it biggest == n


def test_polyx_unmodeled_chain_never_merges_even_when_long():
    # Full-snapshot finding (2026-07-22): unclustered fan-out is driven by UNMODELED
    # (poly-'X') chains at ALL lengths (a 72-'X' chain bridged 283 clusters). A raw-length
    # gate would let a long poly-'X' through and seed a mega-component; the modeled-residue
    # gate (0 modeled) keeps them apart. Goes red if the gate reverts to raw length.
    long_polyx = "X" * 80  # 80 residues, 0 modeled
    cr = build_clusters(
        [_mixed_record("EE01", 1, long_polyx), _mixed_record("EE02", 2, long_polyx)], _cfg()
    )
    assert cr.entry_to_cluster["EE01"] != cr.entry_to_cluster["EE02"]


def test_rebuild_diff_reports_entry_moves_into_train(tmp_path, capsys):
    # The entry-level growth signal: a same-config rebuild reports how many prior entries
    # changed split AND how many were absorbed into train (the unsafe hash-merge direction).
    import json as _json

    from ifsplit.cli import _report_rebuild_migration

    cfg = _cfg()
    (tmp_path / "train.json").write_text(_json.dumps(["B"]))
    (tmp_path / "val.json").write_text(_json.dumps([]))
    (tmp_path / "test.json").write_text(_json.dumps(["A"]))
    (tmp_path / "dataset.lock").write_text(_json.dumps({"config_hash": cfg.config_hash()}))
    # New assignment: A moved test -> train (into train), B moved train -> test (out of train)
    # — both directions contaminate, for different downstream users.
    _report_rebuild_migration(cfg, tmp_path, {"A": "train", "B": "test"})
    out = capsys.readouterr().out
    assert "CHANGED split" in out
    assert "1 INTO train" in out
    assert "1 OUT of train" in out


def test_rebuild_diff_silent_on_first_build(tmp_path, capsys):
    from ifsplit.cli import _report_rebuild_migration

    _report_rebuild_migration(_cfg(), tmp_path, {"A": "train"})  # no prior split in the dir
    assert capsys.readouterr().out == ""


# ----------------------------- Stage 3: filter ----------------------------- #
def test_filter_keeps_protein_entries(sample_entries):
    recs = [CandidateRecord.from_data_api(e) for e in sample_entries.values()]
    kept, drops = filter_candidates(recs, _cfg())
    assert {r.entry_id for r in kept} == {"1A1F", "4HHB"}
    assert drops == []


def test_filter_drops_no_protein(sample_entries):
    # Strip 1A1F down to its DNA entities only.
    rec = CandidateRecord.from_data_api(sample_entries["1A1F"])
    rec = rec.model_copy(
        update={"polymer_entities": [e for e in rec.polymer_entities if e.is_nucleic]}
    )
    kept, drops = filter_candidates([rec], _cfg())
    assert kept == []
    assert drops[0]["reason"] == DROP_NO_PROTEIN


def test_filter_drops_too_large(sample_entries):
    kept, drops = filter_candidates(
        [CandidateRecord.from_data_api(sample_entries["4HHB"])],
        _cfg(max_total_residues=100),  # 4HHB assembly has 574
    )
    assert kept == []
    assert drops[0]["reason"] == DROP_TOO_LARGE
    assert drops[0]["residues"] == 574


def test_size_cap_boundary_keeps_exactly_max_residues():
    # max_total_residues is the max KEPT (drop if > it), so an entry with exactly
    # that many residues is kept — LigandMPNN's "< 6000" is keep <= 5999, not <= 5998.
    rec = _seq_record("BND1", "A" * 100)  # assembly-1 count == 100
    kept, _ = filter_candidates([rec], _cfg(max_total_residues=100))
    assert [r.entry_id for r in kept] == ["BND1"]
    # One residue over the cap is dropped.
    kept2, drops = filter_candidates([rec], _cfg(max_total_residues=99))
    assert kept2 == []
    assert drops[0]["reason"] == DROP_TOO_LARGE


def test_filter_drops_poly_unk_sequence():
    # Every protein chain is all-X (poly-UNK): no known residue identities, so no
    # learnable inverse-folding label -> always dropped (even at the default min=0).
    kept, drops = filter_candidates([_seq_record("UNK1", "X" * 80)], _cfg())
    assert kept == []
    assert drops[0]["reason"] == DROP_NO_SEQUENCE


def test_filter_keeps_partially_modeled_sequence():
    # A chain with even a few modeled residues is usable at the default (min=0).
    kept, _ = filter_candidates([_seq_record("OK1", "X" * 70 + "ACDEFGHIKL")], _cfg())
    assert [r.entry_id for r in kept] == ["OK1"]


def test_min_modeled_residues_drops_short_chain():
    rec = _seq_record("SHRT", "ACDEFGHIKLMNPQR")  # 15 modeled residues
    kept, drops = filter_candidates([rec], _cfg(min_modeled_residues=20))
    assert kept == []
    assert drops[0]["reason"] == DROP_SEQUENCE_TOO_SHORT
    assert drops[0]["modeled"] == 15
    # Off by default (min=0): the same record is kept.
    kept2, _ = filter_candidates([rec], _cfg())
    assert [r.entry_id for r in kept2] == ["SHRT"]


def test_filter_drops_high_clashscore(sample_entries):
    # 4HHB's real clashscore is 142; a 40 cap drops it but keeps 1A1F (4.5).
    recs = [CandidateRecord.from_data_api(e) for e in sample_entries.values()]
    kept, drops = filter_candidates(recs, _cfg(max_clashscore=40.0))
    assert {r.entry_id for r in kept} == {"1A1F"}
    assert drops[0]["reason"] == DROP_CLASHSCORE
    assert drops[0]["value"] == 142.32


def test_filter_keeps_when_metric_absent(sample_entries):
    # 4HHB has no diffraction summary -> rfree is None -> an rfree cap can't drop it.
    rec = CandidateRecord.from_data_api(sample_entries["4HHB"])
    assert rec.quality.rfree is None
    kept, drops = filter_candidates([rec], _cfg(max_rfree=0.25))
    assert [r.entry_id for r in kept] == ["4HHB"]
    assert drops == []


def test_filter_drops_high_rfree(sample_entries):
    # 1A1F has DCC_Rfree 0.21; a 0.20 cap drops it.
    rec = CandidateRecord.from_data_api(sample_entries["1A1F"])
    kept, drops = filter_candidates([rec], _cfg(max_rfree=0.20))
    assert kept == []
    assert drops[0]["reason"] == DROP_RFREE


def test_require_validation_report_drops_reportless_entry():
    # _protein_record() builds a record with no validation summary at all.
    rec = _protein_record("AAA1", [10])
    assert rec.quality.has_report is False
    kept, drops = filter_candidates([rec], _cfg(require_validation_report=True))
    assert kept == []
    assert drops[0]["reason"] == DROP_NO_VALIDATION


# ------------------------ Stage 3: resolution re-filter -------------------- #
def test_filter_resolution_refilter_is_auditable():
    # Stage 3 re-derives the resolution cut (Search applied it too) so it is auditable
    # from candidates.jsonl. An entry over the cap is dropped; at the cap it is kept.
    over = _seq_record("RES1", "A" * 50).model_copy(update={"resolution_A": 3.8})
    kept, drops = filter_candidates([over], _cfg())  # default cap 3.5
    assert kept == []
    assert drops[0]["reason"] == DROP_RESOLUTION
    assert drops[0]["resolution"] == 3.8
    at_cap = _seq_record("RES2", "A" * 50).model_copy(update={"resolution_A": 3.5})
    assert [r.entry_id for r in filter_candidates([at_cap], _cfg())[0]] == ["RES2"]


def test_filter_resolution_missing_is_kept():
    rec = _seq_record("RESN", "A" * 50).model_copy(update={"resolution_A": None})
    assert [r.entry_id for r in filter_candidates([rec], _cfg())[0]] == ["RESN"]


def test_per_method_resolution_cap():
    cfg = _cfg(resolution_max_A_by_method={"ELECTRON MICROSCOPY": 3.0})
    # A 3.2 A cryo-EM entry is dropped by the tighter EM cap...
    em = _seq_record("EM1", "A" * 50).model_copy(
        update={"resolution_A": 3.2, "methods": ["ELECTRON MICROSCOPY"]}
    )
    kept, drops = filter_candidates([em], cfg)
    assert kept == []
    assert drops[0]["reason"] == DROP_RESOLUTION
    # ...but a 3.2 A X-ray entry passes (its cap is still the global 3.5).
    xr = _seq_record("XR1", "A" * 50).model_copy(
        update={"resolution_A": 3.2, "methods": ["X-RAY DIFFRACTION"]}
    )
    assert [r.entry_id for r in filter_candidates([xr], cfg)[0]] == ["XR1"]


def test_search_resolution_cap_is_loosest():
    # The Search query must pull a superset: the loosest cap across enabled methods.
    assert (
        _cfg(resolution_max_A_by_method={"ELECTRON MICROSCOPY": 3.0}).search_resolution_cap() == 3.5
    )
    assert (
        _cfg(
            resolution_max_A_by_method={"X-RAY DIFFRACTION": 4.0, "ELECTRON MICROSCOPY": 3.0}
        ).search_resolution_cap()
        == 4.0
    )


# ---------------------- Stage 3: cryo-EM map-fit floor --------------------- #
def _with_em_inclusion(entry_id: str, value: float | None) -> CandidateRecord:
    rec = _seq_record(entry_id, "A" * 50)
    return rec.model_copy(
        update={"quality": rec.quality.model_copy(update={"em_backbone_inclusion": value})}
    )


def test_em_backbone_inclusion_floor_drops_low_fit():
    cfg = _cfg(min_em_backbone_inclusion=0.7)
    kept, drops = filter_candidates([_with_em_inclusion("EMLO", 0.6)], cfg)
    assert kept == []
    assert drops[0]["reason"] == DROP_EM_INCLUSION
    assert drops[0]["value"] == 0.6
    # At/above the floor is kept.
    assert [r.entry_id for r in filter_candidates([_with_em_inclusion("EMHI", 0.85)], cfg)[0]] == [
        "EMHI"
    ]


def test_em_floor_ignores_entries_without_the_metric():
    # X-ray entries have no em_backbone_inclusion -> the floor never drops them.
    cfg = _cfg(min_em_backbone_inclusion=0.7)
    assert [r.entry_id for r in filter_candidates([_with_em_inclusion("XNOEM", None)], cfg)[0]] == [
        "XNOEM"
    ]


def test_drop_summary_counts():
    drops = [
        {"entry_id": "A", "reason": DROP_NO_PROTEIN},
        {"entry_id": "B", "reason": DROP_NO_PROTEIN},
        {"entry_id": "C", "reason": DROP_TOO_LARGE},
    ]
    assert drop_summary(drops) == {DROP_NO_PROTEIN: 2, DROP_TOO_LARGE: 1}


# ---------------------------- Stage 5: cluster ----------------------------- #
def test_cluster_groups_by_membership(sample_entries, artifact_entry):
    recs = _records(sample_entries, artifact_entry)
    kept, _ = filter_candidates(recs, _cfg())
    cr = build_clusters(kept, _cfg())
    # 4HHB (two protein clusters merge into one component); 1A1F; artifact: 3 entries.
    assert set(cr.entry_to_cluster) == {"4HHB", "1A1F", "pdb_00009xyz"}
    # canonical keys are entity ids (smallest member), not raw integers.
    assert all(":" not in k or k.startswith("singleton:") for k in cr.cluster_members)


def test_cluster_multichain_detected(sample_entries):
    # 4HHB has two different protein clusters (alpha/beta) -> multichain.
    kept, _ = filter_candidates([CandidateRecord.from_data_api(sample_entries["4HHB"])], _cfg())
    cr = build_clusters(kept, _cfg())
    assert "4HHB" in cr.multichain_entries


# ------------------------------ Stage 6: split ----------------------------- #
def test_bucket_is_deterministic_and_unit_range():
    b = bucket("4HHB_1", "snapsplit-v1")
    assert b == bucket("4HHB_1", "snapsplit-v1")
    assert 0.0 <= b < 1.0
    assert bucket("4HHB_1", "other-salt") != b


def test_split_assignment_deterministic(sample_entries, artifact_entry):
    recs = _records(sample_entries, artifact_entry)
    kept, _ = filter_candidates(recs, _cfg())
    cr = build_clusters(kept, _cfg())
    a = assign_splits(cr, _cfg())
    b = assign_splits(cr, _cfg())
    assert a.cluster_split == b.cluster_split
    assert a.entry_split == b.entry_split


def test_no_cluster_leakage_invariant(sample_entries, artifact_entry):
    recs = _records(sample_entries, artifact_entry)
    kept, _ = filter_candidates(recs, _cfg())
    cr = build_clusters(kept, _cfg())
    res = assign_splits(cr, _cfg())
    check_no_leakage(res, cr)  # raises on leakage


def test_registry_pins_assignment(sample_entries, artifact_entry):
    recs = _records(sample_entries, artifact_entry)
    cfg_hash = _cfg(split_strategy="hash")
    kept, _ = filter_candidates(recs, cfg_hash)
    cr = build_clusters(kept, cfg_hash)
    # Force every cluster to "test" via a registry; the hash is overridden.
    reg = {k: "test" for k in cr.cluster_members}
    assert set(assign_splits(cr, cfg_hash, registry=reg).cluster_split.values()) == {"test"}

    # Under "maximal" the same registry must NOT be honored: the holdout ceiling is not
    # negotiable, and a pin that would hold out more than it allows loses to the cap.
    # Honoring it unconditionally is what inverted a real growth rebuild into a
    # 214,796-entry test set.
    mx = assign_splits(cr, _cfg(split_strategy="maximal"), registry=reg)
    n = sum(mx.counts.values())
    assert mx.counts["train"] > 0
    assert mx.counts["val"] + mx.counts["test"] <= 0.2 * n
    assert mx.pinned_entries_reassigned > 0  # overrides reported, never silent


def test_test_minimums_recruit_components_no_leakage(sample_entries, artifact_entry):
    # 1A1F carries a functional metal (bound Zn). With the pure hash it may not be
    # in test; a metal floor of 1 must pull its whole component into test.
    recs = _records(sample_entries, artifact_entry)
    base = _cfg(structural_clustering="off", split_strategy="hash")
    kept, _ = filter_candidates(recs, base)
    class_map = {r.entry_id: classify_components(r, base) for r in kept}
    entry_classes = {eid: info["classes"] for eid, info in class_map.items()}
    cr = build_clusters(kept, base)
    cfg_min = _cfg(
        structural_clustering="off", split_strategy="hash", test_min_per_class={"metal": 1}
    )
    res = assign_splits(cr, cfg_min, entry_classes=entry_classes)
    # The floor is met and the structural no-leakage invariant still holds.
    metal_in_test = sum(
        1 for e, s in res.entry_split.items() if s == "test" and "metal" in entry_classes.get(e, [])
    )
    assert metal_in_test >= 1
    assert res.minimum_shortfalls == {}
    check_no_leakage(res, cr)


def test_test_minimums_report_shortfall_when_supply_short(sample_entries, artifact_entry):
    recs = _records(sample_entries, artifact_entry)
    kept, _ = filter_candidates(recs, _cfg())
    class_map = {r.entry_id: classify_components(r, _cfg()) for r in kept}
    entry_classes = {eid: info["classes"] for eid, info in class_map.items()}
    cr = build_clusters(kept, _cfg())
    # Demand far more metal entries than exist -> shortfall reported, not forced.
    cfg_min = _cfg(test_min_per_class={"metal": 999})
    res = assign_splits(cr, cfg_min, entry_classes=entry_classes)
    assert res.minimum_shortfalls.get("metal", 0) > 0
    check_no_leakage(res, cr)


def test_minimums_off_by_default_matches_pure_hash(sample_entries, artifact_entry):
    recs = _records(sample_entries, artifact_entry)
    kept, _ = filter_candidates(recs, _cfg())
    cr = build_clusters(kept, _cfg())
    base = assign_splits(cr, _cfg())  # no entry_classes, default empty minimums
    assert base.minimum_shortfalls == {}
    # Providing classes but no minimums must not change the assignment.
    class_map = {r.entry_id: classify_components(r, _cfg()) for r in kept}
    ec = {eid: info["classes"] for eid, info in class_map.items()}
    same = assign_splits(cr, _cfg(), entry_classes=ec)
    assert same.entry_split == base.entry_split


def test_fractions_roughly_respected_on_many_keys():
    # Synthetic: hash 3000 distinct keys, check broad proportions hold.
    cfg = _cfg()
    buckets = [split_for_key(f"K{i}", cfg) for i in range(3000)]
    train = buckets.count("train") / len(buckets)
    assert 0.74 < train < 0.86  # ~0.80 with sampling slack


# --------------------- Stages 6/7: manifest + loader ----------------------- #
def _full_manifest(sample_entries, artifact_entry, cfg):
    recs = _records(sample_entries, artifact_entry)
    kept, drops = filter_candidates(recs, cfg)
    class_map = {r.entry_id: classify_components(r, cfg) for r in kept}
    cr = build_clusters(kept, cfg)
    sp = assign_splits(cr, cfg)
    return build_manifest(
        cfg,
        candidates_sha256="deadbeef",
        n_candidates=len(recs),
        drops=drops,
        drop_counts=drop_summary(drops),
        clusters=cr,
        splits=sp,
        class_map=class_map,
    )


def test_manifest_is_deterministic(sample_entries, artifact_entry):
    cfg = _cfg()
    import json

    m1 = json.dumps(_full_manifest(sample_entries, artifact_entry, cfg), sort_keys=True)
    m2 = json.dumps(_full_manifest(sample_entries, artifact_entry, cfg), sort_keys=True)
    assert m1 == m2  # no wall-clock fields -> byte-identical


def test_manifest_has_all_entries(sample_entries, artifact_entry):
    m = _full_manifest(sample_entries, artifact_entry, _cfg())
    # The manifest holds counts only (per-entry lists live in train/val/test.json).
    total = sum(m["splits"]["entry_counts"].values())
    assert total == 3  # 4HHB, 1A1F, pdb_00009xyz


def test_manifest_is_lightweight(sample_entries, artifact_entry):
    import json

    m = _full_manifest(sample_entries, artifact_entry, _cfg())
    # No per-entry arrays in the manifest itself — only counts + a files index.
    assert "entries" not in m["splits"]
    assert "entry_clusters" not in m["splits"]
    assert "classes" not in m["ligands"]
    assert "tiers" not in m["ligands"]
    assert set(m["files"]["splits"]) == {"train", "val", "test"}
    # Sanity: the whole manifest is tiny (well under 10 KB for 3 entries).
    assert len(json.dumps(m)) < 10_000


def test_loader_roundtrip(tmp_path, sample_entries, artifact_entry):
    from ifsplit.manifest import write_classes, write_clusters, write_split_files

    cfg = _cfg()
    recs = _records(sample_entries, artifact_entry)
    kept, drops = filter_candidates(recs, cfg)
    class_map = {r.entry_id: classify_components(r, cfg) for r in kept}
    cr = build_clusters(kept, cfg)
    sp = assign_splits(cr, cfg)
    m = build_manifest(
        cfg,
        candidates_sha256="deadbeef",
        n_candidates=len(recs),
        drops=drops,
        drop_counts=drop_summary(drops),
        clusters=cr,
        splits=sp,
        class_map=class_map,
    )
    write_split_files(sp, class_map, tmp_path)
    write_clusters(cr.entry_to_cluster, tmp_path)
    write_classes(class_map, tmp_path)
    path = write_manifest(m, tmp_path)
    ds = load_dataset(path)
    total = len(ds.train) + len(ds.val) + len(ds.test)
    assert total == 3
    assert ds.config_hash == _cfg().config_hash()


# ------------------------- resplit (offline, no RCSB) ---------------------- #
def test_resplit_reproduces_build_offline(tmp_path, fake_client):
    import argparse

    from ifsplit.cli import _run_pipeline, cmd_resplit
    from ifsplit.enumerate import enumerate_candidates
    from ifsplit.manifest import read_lock, read_manifest

    cfg = _cfg()
    # Reference: enumerate to `ref` (writes candidates.jsonl) + run Stages 3-7 there.
    ref = tmp_path / "ref"
    records, cand_path, sha = enumerate_candidates(cfg, ref, client=fake_client)
    _run_pipeline(cfg, records, sha, ref, limit=None, registry_path=None)
    man_ref = read_manifest(ref / "manifest.json")

    # Resplit re-derives from the SAME candidates.jsonl offline (no client).
    out = tmp_path / "out"
    args = argparse.Namespace(
        config=str(DEFAULT_CONFIG), candidates=str(cand_path), out=str(out), registry=None
    )
    assert cmd_resplit(args) == 0
    man_out = read_manifest(out / "manifest.json")
    # Same snapshot + config -> identical split output; the offline sha matches the
    # sha enumerate computed (the file bytes hash identically).
    assert man_out["splits"]["entry_counts"] == man_ref["splits"]["entry_counts"]
    out_lock = read_lock(out / "dataset.lock")
    assert out_lock["candidates"]["sha256"] == sha
    # The lock records how it was produced: resplit vs a live build.
    assert out_lock["source"] == "resplit"
    assert read_lock(ref / "dataset.lock")["source"] == "build"


# ----------------------------- Phase 7: growth ----------------------------- #
def test_existing_cluster_does_not_move_when_dataset_grows(sample_entries, artifact_entry):
    cfg = _cfg()
    # Snapshot A: just the two sample entries.
    recs_a = [CandidateRecord.from_data_api(e) for e in sample_entries.values()]
    kept_a, _ = filter_candidates(recs_a, cfg)
    cr_a = build_clusters(kept_a, cfg)
    sp_a = assign_splits(cr_a, cfg)

    # Snapshot B: adds a third entry (growth).
    recs_b = _records(sample_entries, artifact_entry)
    kept_b, _ = filter_candidates(recs_b, cfg)
    cr_b = build_clusters(kept_b, cfg)
    sp_b = assign_splits(cr_b, cfg, registry=sp_a.cluster_split)

    # Every cluster present in A keeps its split in B.
    for key, split in sp_a.cluster_split.items():
        assert sp_b.cluster_split[key] == split


def test_growth_bridging_merge_is_honest_and_registry_stable():
    # The case the disjoint-growth test above misses: a later snapshot adds a bridging
    # 2-chain entry C that MERGES A's and B's previously-separate components. Assert at
    # the ENTRY level (where the invariant actually lives), for both strategies of pinning.
    a = _protein_record("AAAA", [1])
    b = _protein_record("BBBB", [2])
    c = _protein_record("CCCC", [1, 2])  # bridges clusters 1 and 2 -> merges the components

    # hash: this test is about registry pinning across a merge, which is the hash path
    # (maximal/balanced pin via the registry too, but their fill is size-driven, so a
    # 2-entry synthetic set has no holdout budget at all).
    def _cfg_g(**kw):
        return _cfg(structural_clustering="off", split_strategy="hash", **kw)

    def _v1_entry_split(s):
        c = _cfg_g(split_salt=s)
        return assign_splits(build_clusters(filter_candidates([a, b], c)[0], c), c).entry_split

    # A salt that, in the v1 (two-component) build, holds B out in test and puts A elsewhere.
    salt = None
    for s in (f"g{i}" for i in range(300)):
        es = _v1_entry_split(s)
        if es["BBBB"] == "test" and es["AAAA"] != "test":
            salt = s
            break
    assert salt is not None, "no salt placed B in test and A elsewhere"
    cfg = _cfg_g(split_salt=salt)

    v1_clusters = build_clusters(filter_candidates([a, b], cfg)[0], cfg)
    v1 = assign_splits(v1_clusters, cfg)
    assert v1_clusters.n_clusters == 2
    registry = dict(v1.cluster_split)  # {AAAA_1: <A's split>, BBBB_1: 'test'}

    v2_clusters = build_clusters(filter_candidates([a, b, c], cfg)[0], cfg)
    assert v2_clusters.n_clusters == 1  # the bridge merged the two components into one

    # Without a registry: B is absorbed into A's component and leaves test — the
    # unavoidable migration, now honest (pinned_reassignments needs a registry to count).
    v2 = assign_splits(v2_clusters, cfg)
    assert v2.entry_split["BBBB"] != "test"
    assert v2.entry_split["AAAA"] == v2.entry_split["BBBB"]  # one component -> one split
    assert v2.pinned_reassignments == 0

    # With the v1 registry: B's held-out (test) pin wins across the merge (test
    # precedence), so B STAYS in test; A's prior pin is overridden and that reassignment
    # is counted, not dropped silently.
    v2r = assign_splits(v2_clusters, cfg, registry=registry)
    check_no_leakage(v2r, v2_clusters)
    assert v2r.entry_split["BBBB"] == "test"  # held-out data preserved across the merge
    assert v2r.entry_split["AAAA"] == "test"  # forced to B's split (single component)
    assert v2r.pinned_reassignments == 1  # A's non-test pin overridden -> reported


def _clusters_from(members):
    """Build a ClusterResult from a {component_key: [entry_ids]} map (for split tests)."""
    from ifsplit.cluster import ClusterResult

    e2c, eraw = {}, {}
    for key, ents in members.items():
        for e in ents:
            e2c[e] = key
            eraw[e] = [key]
    return ClusterResult(
        30, dict(sorted(e2c.items())), dict(sorted(members.items())), dict(sorted(eraw.items()))
    )


def _members(prefix, n, start=0):
    """n components with tail sizes 1..20, so the balanced val/test fill is exercised."""
    return {
        f"{prefix}{i:05d}": [f"{prefix}{i:05d}_{j}" for j in range(1 + i % 20)]
        for i in range(start, start + n)
    }


def test_balanced_growth_stability_needs_registry():
    # A balanced split's val/test fill boundaries scale with total entries, so WITHOUT a
    # registry a grown snapshot moves prior components across splits; the registry pins them.
    cfg = _cfg(split_strategy="balanced")
    a_mem = _members("A", 1000)
    clusters_a = _clusters_from(a_mem)
    clusters_b = _clusters_from({**a_mem, **_members("B", 1000)})  # A + 1000 new components

    sp_a = assign_splits(clusters_a, cfg)
    sp_b_noreg = assign_splits(clusters_b, cfg)  # registry defaults to {}
    moved = sum(1 for k in a_mem if sp_a.cluster_split[k] != sp_b_noreg.cluster_split[k])
    assert moved > 0, "balanced must NOT be assumed growth-stable without a registry"

    sp_b_reg = assign_splits(clusters_b, cfg, registry=sp_a.cluster_split)
    moved_reg = sum(1 for k in a_mem if sp_a.cluster_split[k] != sp_b_reg.cluster_split[k])
    assert moved_reg == 0, "the registry must restore growth stability for balanced"


def test_balanced_rebuild_auto_pins_registry(tmp_path, sample_entries, artifact_entry, capsys):
    # The fix: a balanced rebuild into the same --out auto-adopts the prior registry when
    # the config matches, so the lock records it and the manifest reports growth_stable.
    from ifsplit.cli import _run_pipeline
    from ifsplit.manifest import read_lock, read_manifest

    cfg = _cfg(split_strategy="balanced")
    recs = _records(sample_entries, artifact_entry)
    out = tmp_path / "d"

    _run_pipeline(cfg, recs, "sha", out, limit=None, registry_path=None)  # first build
    assert read_lock(out / "dataset.lock")["split"]["registry_sha256"] is None  # registry-free

    capsys.readouterr()
    _run_pipeline(cfg, recs, "sha", out, limit=None, registry_path=None)  # in-place rebuild
    assert "pinning" in capsys.readouterr().out
    assert read_lock(out / "dataset.lock")["split"]["registry_sha256"] is not None
    assert read_manifest(out / "manifest.json")["splits"]["growth_stable"] is True

    _run_pipeline(cfg, recs, "sha", out, limit=None, registry_path=None, fresh=True)  # opt out
    assert read_lock(out / "dataset.lock")["split"]["registry_sha256"] is None


def test_hash_rebuild_stays_registry_free(tmp_path, sample_entries, artifact_entry):
    # No regression: hash is input-independent, so it is never auto-pinned and stays
    # registry-free (verify can still certify the split output).
    from ifsplit.cli import _run_pipeline
    from ifsplit.manifest import read_lock, read_manifest

    cfg = _cfg(split_strategy="hash")  # hash is registry-free by design
    recs = _records(sample_entries, artifact_entry)
    out = tmp_path / "d"
    _run_pipeline(cfg, recs, "sha", out, limit=None, registry_path=None)
    _run_pipeline(cfg, recs, "sha", out, limit=None, registry_path=None)  # rebuild
    assert read_lock(out / "dataset.lock")["split"]["registry_sha256"] is None
    assert read_manifest(out / "manifest.json")["splits"]["growth_stable"] is True


# ----------------------- Phase 7: novel-fold benchmark --------------------- #
def test_build_fold_benchmark_novel_is_train_unseen():
    from ifsplit.manifest import build_fold_benchmark

    labels = {"T1": ["fam.seen"], "T2": ["fam.novel"], "TR1": ["fam.seen"], "V1": ["fam.novel"]}
    entry_split = {"T1": "test", "T2": "test", "TR1": "train", "V1": "val"}
    fb = build_fold_benchmark(labels, entry_split, "scop2")
    assert fb["novel_fold_test"] == ["T2"]  # its fold is not in train; T1's is
    assert fb["summary"]["n_test_classified"] == 2
    assert fb["summary"]["n_test_novel_fold"] == 1
    assert fb["fold_groups"]["fam.seen"] == {"novel": False, "test_entries": ["T1"]}
    assert fb["fold_groups"]["fam.novel"] == {"novel": True, "test_entries": ["T2"]}
    assert fb["per_entry"]["V1"]["novel_fold"] is True  # val held out + fold unseen
    assert fb["per_entry"]["TR1"]["novel_fold"] is False  # train is never novel
    assert build_fold_benchmark(labels, entry_split, "off") is None


def test_fold_benchmark_decoupled_from_split():
    # structural_clustering OFF but fold_benchmark ON: labels are emitted WITHOUT merging,
    # so the split is byte-identical to a pure off/off build (the decoupling guarantee).
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.1.1"]}),
        _fold_record("BBB2", 20, {"cath": ["1.10.1.1"]}),  # same fold, different seq cluster
    ]
    base = _cfg(structural_clustering="off")  # off / off
    bench = _cfg(structural_clustering="off", fold_benchmark_method="cath")

    cr_base = build_clusters(filter_candidates(recs, base)[0], base)
    cr_bench = build_clusters(filter_candidates(recs, bench)[0], bench)
    assert cr_base.cluster_members == cr_bench.cluster_members  # labels never merged them
    assert cr_bench.entry_families == {}  # structural_clustering off -> no merge labels
    assert cr_bench.entry_fold_labels == {"AAA1": ["1.10.1.1"], "BBB2": ["1.10.1.1"]}
    assert assign_splits(cr_base, base).entry_split == assign_splits(cr_bench, bench).entry_split


def test_fold_benchmark_off_is_config_hash_stable():
    # Off (default) is a pure export toggle: omitted from the hash so legacy configs are
    # unchanged; turning it on changes the hash (it changes the outputs).
    assert "fold_benchmark_method" not in _cfg().canonical_dict()
    assert "fold_benchmark_method" in _cfg(fold_benchmark_method="cath").canonical_dict()
    assert _cfg().config_hash() != _cfg(fold_benchmark_method="cath").config_hash()


def test_fold_benchmark_export_end_to_end(tmp_path):
    from ifsplit.cli import _run_pipeline
    from ifsplit.dataset import load_dataset
    from ifsplit.manifest import read_manifest

    cfg = _cfg(fold_benchmark_method="cath")
    recs = [_fold_record(f"E{i:03d}", i, {"cath": [f"1.10.{i}.1"]}) for i in range(30)]
    out = tmp_path / "d"
    _run_pipeline(cfg, recs, "sha", out, limit=None, registry_path=None)

    for fname in ("folds.json", "fold_groups.json", "novel_fold_test.json"):
        assert (out / fname).exists()
    assert read_manifest(out / "manifest.json")["fold_benchmark"]["method"] == "cath"

    ds = load_dataset(out / "manifest.json")
    novel = ds.test.novel_fold_entries()
    assert set(novel) <= set(ds.test.entry_ids)
    # Every cath fold here is unique to one entry, so a test fold is never in train:
    # all classified test entries are novel-fold, and every fold group is novel.
    assert set(novel) == set(ds.test.entry_ids)
    assert all(g["novel"] for g in ds.fold_groups().values())
    for e in ds.test.entry_ids:
        assert ds.test.is_novel_fold(e) and ds.test.folds_of(e)


# --------------- union-find: structural leakage prevention ----------------- #
def _protein_record(entry_id: str, cluster30_ids: list[int]) -> CandidateRecord:
    """A record whose protein chains sit in the given raw clusters (id at 30%)."""
    pes = [
        PolymerEntity(
            entity_id=f"{entry_id}_{i + 1}",
            polymer_type="Protein",
            seq_len=100,
            seq=_seq_for(cid),
            cluster_ids={30: cid},
        )
        for i, cid in enumerate(cluster30_ids)
    ]
    return CandidateRecord(
        entry_id=entry_id,
        methods=["X-RAY DIFFRACTION"],
        resolution_A=2.0,
        release_date="2020-01-01",
        deposited_residues=100,
        assemblies={f"{entry_id}-1": 100},
        polymer_entities=pes,
        nonpolymer_comps=[],
        bound_components=[],
        affinity_comp_ids=[],
    )


def test_union_find_merges_bridged_clusters_no_leakage():
    cfg = _cfg()
    # X bridges raw clusters 1 and 2 via two chains; Y is in 1, Z is in 2.
    # Without union-find, clusters 1 and 2 could hash to different splits and Y/Z
    # would leak X's sequences across splits. With it, {1,2} is one component.
    recs = [
        _protein_record("X1AA", [1, 2]),
        _protein_record("Y2BB", [1]),
        _protein_record("Z3CC", [2]),
    ]
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    assert cr.n_raw_clusters == 2
    assert cr.n_clusters == 1  # the two raw clusters merged into one component
    res = assign_splits(cr, cfg)
    check_no_leakage(res, cr)  # would raise if 1 and 2 split apart
    assert len(set(res.entry_split.values())) == 1  # all three co-assigned


def test_independent_clusters_can_differ_and_check_passes():
    cfg = _cfg()
    # Two unrelated single-chain entries: separate components, no shared sequence.
    recs = [_protein_record("AAA1", [10]), _protein_record("BBB2", [20])]
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    assert cr.n_clusters == 2
    res = assign_splits(cr, cfg)
    check_no_leakage(res, cr)  # passes regardless of which splits they land in


def _fold_record(entry_id: str, cluster30: int, families: dict[str, list[str]]) -> CandidateRecord:
    """Single-chain record in raw cluster ``cluster30`` carrying structural ``families``."""
    pe = PolymerEntity(
        entity_id=f"{entry_id}_1",
        polymer_type="Protein",
        seq_len=100,
        seq=_seq_for(cluster30),
        cluster_ids={30: cluster30},
        structural_families=families,
    )
    return CandidateRecord(
        entry_id=entry_id,
        methods=["X-RAY DIFFRACTION"],
        resolution_A=2.0,
        release_date="2020-01-01",
        deposited_residues=100,
        assemblies={f"{entry_id}-1": 100},
        polymer_entities=[pe],
        nonpolymer_comps=[],
        bound_components=[],
        affinity_comp_ids=[],
    )


def test_structural_clustering_merges_same_fold():
    # Two entries in DIFFERENT sequence clusters (10, 20) but the same CATH
    # superfamily. Sequence-only leaves them separable (a fold-leakage risk);
    # cath clustering folds them into one leakage-safe component.
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.490.10"]}),
        _fold_record("BBB2", 20, {"cath": ["1.10.490.10"]}),
    ]
    cfg = _cfg(structural_clustering="cath")
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    assert cr.n_raw_clusters == 2
    assert cr.n_seq_only_components == 2  # sequence edges alone: two components
    assert cr.n_clusters == 1  # ... folded into one by the shared superfamily
    assert cr.structural_method == "cath"
    assert cr.n_structural_families == 1

    # Off -> prior behavior, two separate components.
    cr_off = build_clusters(kept, _cfg(structural_clustering="off"))
    assert cr_off.n_clusters == 2
    assert cr_off.structural_method == "off"


def test_structural_method_is_selectable():
    # Same ECOD family name but different CATH codes: 'cath' keeps them apart,
    # 'ecod' merges them. The classification method is the config's to choose.
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.1.1"], "ecod": ["Globin-like"]}),
        _fold_record("BBB2", 20, {"cath": ["2.20.2.2"], "ecod": ["Globin-like"]}),
    ]
    kept, _ = filter_candidates(recs, _cfg())
    assert build_clusters(kept, _cfg(structural_clustering="cath")).n_clusters == 2
    assert build_clusters(kept, _cfg(structural_clustering="ecod")).n_clusters == 1


def test_union_merges_on_any_authority():
    # "union" merges chains that agree under ANY of CATH/ECOD/SCOP2. Here they differ under
    # CATH and carry no SCOP2, but share an ECOD family — cath keeps them apart, union merges.
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.1.1"], "ecod": ["Bcl-2"]}),
        _fold_record("BBB2", 20, {"cath": ["2.20.2.2"], "ecod": ["Bcl-2"]}),
    ]
    kept, _ = filter_candidates(recs, _cfg())
    assert build_clusters(kept, _cfg(structural_clustering="cath")).n_clusters == 2
    assert build_clusters(kept, _cfg(structural_clustering="union")).n_clusters == 1
    # Namespacing: a CATH code that equals an (unrelated) ECOD name must NOT merge under union.
    recs2 = [
        _fold_record("CCC1", 30, {"cath": ["shared"]}),
        _fold_record("DDD2", 40, {"ecod": ["shared"]}),
    ]
    kept2, _ = filter_candidates(recs2, _cfg())
    assert build_clusters(kept2, _cfg(structural_clustering="union")).n_clusters == 2


def test_cath_key_is_name_stable_but_scop2_key_is_name_sensitive():
    # CATH keys on the stable superfamily code, so a display-name change does NOT change
    # the grouping key. ECOD/SCOP2 key on the free-text name (their annotation_id is
    # per-domain), so a rename DOES change the key — the documented fresh-rebuild
    # limitation (#6): a locked build still reproduces exactly via candidates.jsonl.
    from ifsplit.schema import structural_families_from_instances

    def _inst(atype, ann_id, name):
        return [
            {
                "rcsb_polymer_instance_annotation": [
                    {"type": atype, "annotation_id": ann_id, "name": name}
                ]
            }
        ]

    cath_a = structural_families_from_instances(_inst("CATH", "1.10.490.10", "Globins"))
    cath_b = structural_families_from_instances(_inst("CATH", "1.10.490.10", "Globins (renamed)"))
    assert cath_a["cath"] == cath_b["cath"]  # keyed on the code -> rename-stable

    scop_a = structural_families_from_instances(_inst("SCOP2", "8039836", "Globin-like"))
    scop_b = structural_families_from_instances(_inst("SCOP2", "8039836", "Globin fold"))
    assert scop_a["scop2"] != scop_b["scop2"]  # keyed on the name -> rename-sensitive (known)


def test_structural_clustering_keeps_split_leakage_safe():
    # A fold shared across two entries must land in ONE split, never straddle.
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.490.10"]}),
        _fold_record("BBB2", 20, {"cath": ["1.10.490.10"]}),
    ]
    cfg = _cfg(structural_clustering="cath")
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    res = assign_splits(cr, cfg)
    check_no_leakage(res, cr)
    assert len(set(res.entry_split.values())) == 1  # same fold -> same split


def test_balanced_strategy_caps_dominant_fold_and_balances_entries():
    # 400 entries share one dominant cluster (a mega-fold); 200 are singletons.
    # Per-component hashing would let the 400-entry fold balloon a split; balanced
    # caps it to train and fills val/test to their ENTRY targets from the tail.
    cfg = _cfg(split_strategy="balanced")
    recs = [_protein_record(f"D{i:04d}", [1]) for i in range(400)]
    recs += [_protein_record(f"S{i:04d}", [1000 + i]) for i in range(200)]
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    res = assign_splits(cr, cfg)
    check_no_leakage(res, cr)
    assert res.strategy == "balanced"
    assert res.capped_folds == 1  # the 400-entry fold -> train
    c = res.counts
    assert (c["train"], c["val"], c["test"]) == (480, 60, 60)  # 80/10/10 by entries
    assert not res.balance_gaps


def test_balanced_strategy_reports_thin_tail_gap():
    # Almost everything in one fold: the tail can't fill val+test to 10% each.
    # The gap is reported, not forced (never breaks leakage safety to hit a target).
    cfg = _cfg(split_strategy="balanced")
    recs = [_protein_record(f"D{i:04d}", [1]) for i in range(580)]
    recs += [_protein_record(f"S{i:04d}", [1000 + i]) for i in range(20)]
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    res = assign_splits(cr, cfg)
    check_no_leakage(res, cr)
    assert res.balance_gaps  # tail too thin -> reported
    assert "val" in res.balance_gaps


def test_test_min_does_not_recruit_capped_fold_under_balanced():
    # A dominant fold (capped to train under balanced) that ALSO carries the floored
    # class must NOT be pulled into test to meet a small floor — that would blow the
    # entry balance (the bug). The floor is met from the small-fold tail instead.
    cfg = _cfg(split_strategy="balanced", test_min_per_class={"metal": 30})
    recs = [_protein_record(f"D{i:04d}", [1]) for i in range(400)]  # one 400-entry fold
    recs += [_protein_record(f"S{i:04d}", [1000 + i]) for i in range(200)]  # size-1 tail
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    # metal in the mega-fold AND in 40 small components (a satisfiable tail supply).
    entry_classes = {f"D{i:04d}": ["metal"] for i in range(400)}
    entry_classes.update({f"S{i:04d}": ["metal"] for i in range(40)})
    res = assign_splits(cr, cfg, entry_classes=entry_classes)
    check_no_leakage(res, cr)
    mega = cr.entry_to_cluster["D0000"]
    assert res.cluster_split[mega] == "train"  # capped fold stays in train, not recruited
    assert res.counts["test"] < 200  # not ballooned by the 400-entry fold
    assert res.minimum_shortfalls == {}  # floor met from the small-fold tail
    metal_in_test = sum(
        1 for e, s in res.entry_split.items() if s == "test" and "metal" in entry_classes.get(e, [])
    )
    assert metal_in_test >= 30


def test_test_min_reports_shortfall_when_only_capped_fold_has_class():
    # If the class exists ONLY in a capped dominant fold, the floor cannot be met
    # without blowing the balance — so it is reported as a shortfall, not forced.
    cfg = _cfg(split_strategy="balanced", test_min_per_class={"metal": 50})
    recs = [_protein_record(f"D{i:04d}", [1]) for i in range(400)]
    recs += [_protein_record(f"S{i:04d}", [1000 + i]) for i in range(200)]
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    entry_classes = {f"D{i:04d}": ["metal"] for i in range(400)}  # metal only in the mega-fold
    res = assign_splits(cr, cfg, entry_classes=entry_classes)
    mega = cr.entry_to_cluster["D0000"]
    assert res.cluster_split[mega] == "train"  # not recruited despite the unmet floor
    assert res.minimum_shortfalls.get("metal") == 50  # honest shortfall, not a silent blowup


# ---------- negative leakage: the guard MUST fire on a corrupted partition ------- #
# These prove check_no_leakage is a real invariant, not a happy-path pass — it
# actually raises when a sequence cluster or a fold (super)family spans two splits.


def test_check_no_leakage_fires_on_shared_sequence_cluster():
    # X bridges raw clusters 1 & 2; Y is in 1. A partition that puts X and Y in
    # different splits leaks X's sequence via raw cluster 1 — the guard must catch it.
    cfg = _cfg()
    recs = [
        _protein_record("X1AA", [1, 2]),
        _protein_record("Y2BB", [1]),
        _protein_record("Z3CC", [2]),
    ]
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    res = assign_splits(cr, cfg)
    check_no_leakage(res, cr)  # the real (valid) partition passes
    # Corrupt it: move Y to a different split from X (they share raw cluster 1).
    res.entry_split = dict(res.entry_split)
    res.entry_split["Y2BB"] = "test" if res.entry_split["X1AA"] != "test" else "train"
    with pytest.raises(AssertionError, match="raw cluster"):
        check_no_leakage(res, cr)


def test_check_no_leakage_fires_on_fold_span_when_structural_on():
    # Two entries with DIFFERENT sequence clusters but the SAME CATH superfamily are
    # one component under structural clustering. Splitting them apart shares no
    # sequence cluster (sequence check passes) — only the fold-level guard catches
    # it, proving structural_clustering's guarantee is actually enforced.
    cfg = _cfg(structural_clustering="cath")
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.490.10"]}),
        _fold_record("BBB2", 20, {"cath": ["1.10.490.10"]}),
    ]
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    assert cr.n_clusters == 1  # merged into one component by shared fold
    res = assign_splits(cr, cfg)
    check_no_leakage(res, cr)
    res.entry_split = {"AAA1": "train", "BBB2": "test"}  # force same fold across splits
    with pytest.raises(AssertionError, match="fold leakage"):
        check_no_leakage(res, cr)


def test_fold_leakage_guard_is_a_noop_when_structural_off():
    # Same two same-fold entries, structural OFF: distinct sequence clusters, so
    # they may legitimately land in different splits. entry_families is empty, so
    # the fold guard must NOT fire (no over-merging beyond what the user asked for).
    cfg = _cfg(structural_clustering="off")
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.490.10"]}),
        _fold_record("BBB2", 20, {"cath": ["1.10.490.10"]}),
    ]
    kept, _ = filter_candidates(recs, cfg)
    cr = build_clusters(kept, cfg)
    assert cr.entry_families == {}  # off -> no fold edges recorded
    res = assign_splits(cr, cfg)
    res.entry_split = {"AAA1": "train", "BBB2": "test"}
    check_no_leakage(res, cr)  # must NOT raise — distinct folds may differ when off


def test_single_chain_only_filter(sample_entries):
    # 1A1F is a protein+DNA complex (multiple polymer entities); 4HHB has two protein
    # entities (alpha/beta). single_chain_only keeps only single-protein-entity records.
    recs = [CandidateRecord.from_data_api(e) for e in sample_entries.values()]
    kept, drops = filter_candidates(recs, _cfg(single_chain_only=True))
    assert all(len(r.polymer_entities) == 1 for r in kept)
    reasons = {d["entry_id"]: d["reason"] for d in drops}
    assert reasons.get("1A1F") == "not_single_chain"  # protein+DNA -> dropped
    # A genuine single-chain record passes.
    kept2, _ = filter_candidates([_seq_record("ONE1", "A" * 80)], _cfg(single_chain_only=True))
    assert [r.entry_id for r in kept2] == ["ONE1"]
    # Off by default: the complex is kept.
    kept3, _ = filter_candidates(recs, _cfg())
    assert "1A1F" in {r.entry_id for r in kept3}


def test_manifest_tier_reason_histogram(sample_entries, artifact_entry):
    # The tier-reason histogram summarizes every curation call (a "tier:reason" key
    # per component) so the distribution is auditable without the per-component file.
    m = _full_manifest(sample_entries, artifact_entry, _cfg())
    trc = m["ligands"]["tier_reason_counts"]
    assert isinstance(trc, dict) and trc and all(isinstance(v, int) for v in trc.values())
    assert all(":" in k for k in trc)  # keys are "tier:reason"


def test_manifest_fold_coverage_counts_distinct_folds():
    # per_split_fold_coverage counts the distinct structural families held in each
    # split, plus the unclassified count — the residual-leakage ceiling: entries no
    # fold taxonomy classifies are held out by sequence only, not by fold.
    cfg = _cfg(structural_clustering="cath")
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.1.1"]}),
        _fold_record("BBB2", 20, {"cath": ["2.20.2.2"]}),
        _fold_record("CCC3", 30, {}),  # no fold classification -> unclassified
    ]
    kept, drops = filter_candidates(recs, cfg)
    class_map = {r.entry_id: classify_components(r, cfg) for r in kept}
    cr = build_clusters(kept, cfg)
    sp = assign_splits(cr, cfg)
    m = build_manifest(
        cfg,
        candidates_sha256="deadbeef",
        n_candidates=len(recs),
        drops=drops,
        drop_counts=drop_summary(drops),
        clusters=cr,
        splits=sp,
        class_map=class_map,
    )
    cov = m["splits"]["per_split_fold_coverage"]
    assert set(cov) == {"train", "val", "test"}
    for c in cov.values():
        assert set(c) == {
            "total_entries",
            "classified_entries",
            "unclassified_entries",
            "n_distinct_folds",
        }
        assert c["unclassified_entries"] == c["total_entries"] - c["classified_entries"]
    assert sum(c["n_distinct_folds"] for c in cov.values()) == 2  # two distinct folds
    assert sum(c["classified_entries"] for c in cov.values()) == 2
    assert sum(c["total_entries"] for c in cov.values()) == 3
    assert sum(c["unclassified_entries"] for c in cov.values()) == 1  # CCC3, unclassified


def test_summarize_manifest_reports_residual_ceiling(tmp_path, capsys):
    # `stats` surfaces the unclassified fraction per split (the residual-leakage
    # ceiling) whenever fold-aware clustering is on.
    cfg = _cfg(structural_clustering="cath")
    recs = [
        _fold_record("AAA1", 10, {"cath": ["1.10.1.1"]}),
        _fold_record("CCC3", 30, {}),  # unclassified
    ]
    kept, drops = filter_candidates(recs, cfg)
    class_map = {r.entry_id: classify_components(r, cfg) for r in kept}
    cr = build_clusters(kept, cfg)
    sp = assign_splits(cr, cfg)
    m = build_manifest(
        cfg,
        candidates_sha256="deadbeef",
        n_candidates=len(recs),
        drops=drops,
        drop_counts=drop_summary(drops),
        clusters=cr,
        splits=sp,
        class_map=class_map,
    )
    path = write_manifest(m, tmp_path)
    assert summarize_manifest(path) == 0
    assert "unclassified" in capsys.readouterr().out
