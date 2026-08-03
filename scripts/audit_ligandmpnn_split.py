"""Audit LigandMPNN's *actual* train/val/test split for cross-split leakage.

LigandMPNN (Dauparas et al., Nature Methods 2025) splits the PDB by **30% sequence
identity clusters only** ("The train-test split was based on protein sequences
clustered at a 30% sequence identity cutoff"), assigns whole clusters to
train/val/test, and holds out a "nonoverlapping subset" of 317 small-molecule,
74 nucleotide and 83 metal structures as the test set. That prevents *sequence*
memorization but says nothing about **structural (fold) homology below 30%
identity**, **shared ligand context**, or **multi-chain complex bridges** — the
things a structure->sequence model actually keys on.

Leakage is measured for **both** held-out splits, because they fail differently:
  - leakage into TEST inflates the *reported* recovery number;
  - leakage into VALIDATION is worse — val loss drives early-stopping, checkpoint
    selection and hyperparameter tuning, so a fold-leaked val set means the model
    that gets *shipped* is the one that best memorized, not the one that best
    generalized. Every model-selection decision is made on a leaked signal.

This audit takes LigandMPNN's **real** partition (its published train.json /
valid.json / test_*.json entry lists) and overlays an **independent structural
authority** — RCSB's precomputed CATH / ECOD / SCOP2 classifications and current
30%-sequence clusters — to measure what leaks past their split. It is deliberately
**non-circular**: we never re-run their mmseqs2 30% clustering and compare it to
itself; we apply a *different, stronger* criterion (fold-superfamily identity from
an external classifier) to the partition they actually shipped, and count the
leakage that criterion exposes.

Reuses IF-Split's own machinery: ``ifsplit.ligands`` for ligand tiering/classes
and ``ifsplit.cluster.build_clusters`` for the leakage-safe union-find components.

Inputs (produced by ``scripts/fetch_external_split.py``):
  candidates.jsonl  canonical CandidateRecords for the split's entries
  splits.json       entry_id(UPPER) -> {"split": train|valid|test, "test_class": ...}

Usage:
  uv run python scripts/audit_ligandmpnn_split.py <candidates.jsonl> <splits.json> [out.json]
"""

from __future__ import annotations

import json
import statistics
import sys
from collections import Counter, defaultdict

from ifsplit.cluster import build_clusters
from ifsplit.config import load_config
from ifsplit.ligands import (
    CLASS_METAL,
    CLASS_NUCLEIC_ACID,
    CLASS_SMALL_MOLECULE,
    classify_components,
)
from ifsplit.schema import STRUCTURAL_METHODS, CandidateRecord, read_candidates_jsonl

TEST_CLASS_TO_LIGAND = {
    "metal": CLASS_METAL,
    "small_molecule": CLASS_SMALL_MOLECULE,
    "nucleotide": CLASS_NUCLEIC_ACID,
}


def entry_families(rec: CandidateRecord, method: str) -> set[str]:
    """Structural (super)family keys for one entry under one classification method."""
    fams: set[str] = set()
    for e in rec.polymer_entities:
        if e.is_protein:
            fams.update(e.structural_families.get(method, []))
    return fams


def pct(n: int, d: int) -> str:
    return f"{100 * n / d:5.1f}%" if d else "  n/a"


def fold_leakage(name: str, held: list, train: list, train_fams_by_m: dict) -> dict:
    """A. Does each held-out entry's fold already appear in training?"""
    print("\n" + "-" * 78)
    print(f"A. FOLD LEAKAGE into {name} — is each {name} entry's fold already in TRAINING?")
    print("-" * 78)
    print(
        f"{'method':6s} {name + ' classified':>18s} {'fold-seen in train':>20s} "
        f"{'NOVEL (held-out fold)':>22s}"
    )
    fold: dict = {}
    for m in STRUCTURAL_METHODS:
        classified = seen = 0
        novel_ids = []
        for r in held:
            fams = entry_families(r, m)
            if not fams:
                continue
            classified += 1
            if fams & train_fams_by_m[m]:
                seen += 1
            else:
                novel_ids.append(r.entry_id)
        fold[m] = {"classified": classified, "fold_seen": seen, "novel": len(novel_ids)}
        print(
            f"{m:6s} {classified:>18d} {seen:>13d} ({pct(seen, classified)}) "
            f"{len(novel_ids):>15d} ({pct(len(novel_ids), classified)})"
        )

    # Combined: fold-seen if seen under ANY method that classifies it; NOVEL only if
    # classified by >=1 method and unseen by all.
    seen = novel = unclassified = 0
    novel_ids = []
    for r in held:
        any_classified = any_seen = False
        for m in STRUCTURAL_METHODS:
            fams = entry_families(r, m)
            if fams:
                any_classified = True
                if fams & train_fams_by_m[m]:
                    any_seen = True
        if not any_classified:
            unclassified += 1
        elif any_seen:
            seen += 1
        else:
            novel += 1
            novel_ids.append(r.entry_id)
    print(
        f"\nCOMBINED (any method): of {len(held)} {name} entries — {seen} fold-seen, "
        f"{novel} genuinely NOVEL-fold, {unclassified} unclassified by all three."
    )
    print(
        f"  => the 'honest' held-out-fold {name} set is {novel} structures "
        f"({pct(novel, len(held))} of {name})."
    )
    fold["combined"] = {
        "n": len(held),
        "fold_seen": seen,
        "novel": novel,
        "unclassified": unclassified,
        "novel_ids": sorted(novel_ids),
    }
    return fold


def memorization_pressure(name: str, held: list, train: list) -> dict:
    """B. How many training structures share a held-out entry's fold?"""
    print("\n" + "-" * 78)
    print(f"B. MEMORIZATION PRESSURE — training structures sharing a {name} entry's fold")
    print("-" * 78)
    out: dict = {}
    for m in ("scop2", "cath", "ecod"):
        fam_train: dict[str, set[str]] = defaultdict(set)
        for r in train:
            for fam in entry_families(r, m):
                fam_train[fam].add(r.entry_id)
        pool_sizes = []
        for r in held:
            fams = entry_families(r, m)
            if not fams:
                continue
            pool = set().union(*[fam_train.get(f, set()) for f in fams])
            pool_sizes.append(len(pool))
        if pool_sizes:
            pool_sizes.sort()
            med = statistics.median(pool_sizes)
            mean = statistics.mean(pool_sizes)
            p90 = pool_sizes[9 * len(pool_sizes) // 10]
            zero = sum(1 for x in pool_sizes if x == 0)
            print(f"  {m:6s}: median {med:.0f}, mean {mean:.0f}, p90 {p90}  (novel/0-pool: {zero})")
            out[m] = {"median": med, "mean": mean, "p90": p90, "zero_pool": zero}
    return out


def main() -> None:
    cand_path = sys.argv[1]
    splits_path = sys.argv[2]
    out_path = sys.argv[3] if len(sys.argv) > 3 else None

    cfg = load_config("config/default.yaml")
    records = read_candidates_jsonl(cand_path)
    with open(splits_path) as fh:
        splits = json.load(fh)
    by_id = {r.entry_id.upper(): r for r in records}

    # Partition record lists by the split's ACTUAL assignment.
    split_of = {eid: meta["split"] for eid, meta in splits.items()}
    test_class_of = {eid: meta.get("test_class") for eid, meta in splits.items()}
    resolved = {eid: r for eid, r in by_id.items() if eid in split_of}
    train = [r for eid, r in resolved.items() if split_of[eid] == "train"]
    test = [r for eid, r in resolved.items() if split_of[eid] == "test"]
    val = [r for eid, r in resolved.items() if split_of[eid] == "valid"]

    n_requested = len(splits)
    n_resolved = len(resolved)
    print("=" * 78)
    print("MPNN split leakage audit  (independent criterion: RCSB CATH/ECOD/SCOP2)")
    print("=" * 78)
    print(f"entries in split lists : {n_requested}")
    print(f"resolved via RCSB      : {n_resolved}  (obsoleted: {n_requested - n_resolved})")
    print(f"  train {len(train)}   val {len(val)}   test {len(test)}")

    summary: dict = {
        "n_requested": n_requested,
        "n_resolved": n_resolved,
        "n_train": len(train),
        "n_val": len(val),
        "n_test": len(test),
    }

    cls = {eid: classify_components(r, cfg) for eid, r in resolved.items()}
    train_fams_by_m = {
        m: set().union(*[entry_families(r, m) for r in train]) if train else set()
        for m in STRUCTURAL_METHODS
    }

    # ---- A + B, for BOTH held-out splits (val leaks worse: it drives model selection) ----
    holdouts = [("VALIDATION", "val", val), ("TEST", "test", test)]
    for label, key, held in holdouts:
        if not held:
            continue
        summary[f"fold_leakage_{key}"] = fold_leakage(label, held, train, train_fams_by_m)
        summary[f"memorization_pressure_{key}"] = memorization_pressure(label, held, train)

    # ---- fold leakage by test ligand class (test only) ----
    print("\n  TEST fold leakage by ligand class (combined, any method):")
    per_class = defaultdict(lambda: {"seen": 0, "novel": 0, "unclassified": 0, "n": 0})
    for r in test:
        klass = test_class_of[r.entry_id.upper()]
        d = per_class[klass]
        d["n"] += 1
        any_classified = any(entry_families(r, m) for m in STRUCTURAL_METHODS)
        any_seen = any(entry_families(r, m) & train_fams_by_m[m] for m in STRUCTURAL_METHODS)
        if not any_classified:
            d["unclassified"] += 1
        elif any_seen:
            d["seen"] += 1
        else:
            d["novel"] += 1
    for klass, d in sorted(per_class.items(), key=lambda kv: str(kv[0])):
        print(
            f"    {klass!s:16s} n={d['n']:4d}  fold-seen={d['seen']:4d}  "
            f"novel={d['novel']:3d}  unclassified={d['unclassified']:3d}"
        )
    summary["test_fold_leakage_by_class"] = {str(k): v for k, v in per_class.items()}

    # ============================== C. independent 30%-cluster overlap, both holdouts
    print("\n" + "-" * 78)
    print("C. SEQUENCE-CLUSTER OVERLAP under RCSB's INDEPENDENT 30% clustering")
    print("   (does each held-out split survive a second, independent 30% clustering?)")
    print("-" * 78)
    level = cfg.identity_level
    clust_splits: dict[int, set[str]] = defaultdict(set)
    for eid, r in resolved.items():
        for e in r.polymer_entities:
            if e.is_protein and level in e.cluster_ids:
                clust_splits[e.cluster_ids[level]].add(split_of[eid])
    summary["rcsb_30pct_train_overlap"] = {}
    for label, key, held in holdouts:
        if not held:
            continue
        leaky = 0
        for r in held:
            s: set[str] = set()
            for e in r.polymer_entities:
                if e.is_protein and level in e.cluster_ids:
                    s |= clust_splits[e.cluster_ids[level]]
            if "train" in s:
                leaky += 1
        print(f"  {label:11s}: {leaky}/{len(held)} ({pct(leaky, len(held))}) co-cluster with TRAIN")
        summary["rcsb_30pct_train_overlap"][key] = {"leaky": leaky, "n": len(held)}
    print("  (caveat: RCSB DIAMOND clustering & 2026 membership != their Dec-2022 mmseqs2.)")

    # ============================== D. ligand-context leakage (test only)
    print("\n" + "-" * 78)
    print("D. LIGAND-CONTEXT LEAKAGE — is the TEST conditioning context seen in training?")
    print("-" * 78)
    train_comps: set[str] = set()
    train_fold_class: set = set()
    for r in train:
        c = cls[r.entry_id.upper()]
        comps = set(c["metals"]) | set(c["small_molecules"])
        train_comps |= comps
        for fam in entry_families(r, "scop2"):
            for klass in c["classes"]:
                train_fold_class.add((fam, klass))
    comp_seen = comp_total = fc_seen = fc_total = 0
    for r in test:
        c = cls[r.entry_id.upper()]
        comps = set(c["metals"]) | set(c["small_molecules"])
        for comp in comps:
            comp_total += 1
            comp_seen += comp in train_comps
        for fam in entry_families(r, "scop2"):
            for klass in c["classes"]:
                fc_total += 1
                fc_seen += (fam, klass) in train_fold_class
    print(
        f"  test functional ligand COMPONENTS whose comp-id is also functional in train: "
        f"{comp_seen}/{comp_total} ({pct(comp_seen, comp_total)})"
    )
    print(
        f"  test (SCOP2 fold x ligand-class) pairs also present in train: "
        f"{fc_seen}/{fc_total} ({pct(fc_seen, fc_total)})"
    )
    summary["ligand_context_leakage"] = {
        "comp_seen": comp_seen,
        "comp_total": comp_total,
        "fold_class_seen": fc_seen,
        "fold_class_total": fc_total,
    }

    # ============================== E. complex/fold bridge leakage, both holdouts
    print("\n" + "-" * 78)
    print("E. COMPLEX / FOLD BRIDGE LEAKAGE — IF-Split leakage-safe components vs their split")
    print("-" * 78)
    all_records = list(resolved.values())
    summary["bridge_leakage"] = {}
    for method in ("off", "scop2"):
        cr = build_clusters(all_records, cfg.model_copy(update={"structural_clustering": method}))
        comp_splits: dict[str, set[str]] = defaultdict(set)
        for eid, comp in cr.entry_to_cluster.items():
            comp_splits[comp].add(split_of[eid.upper()])
        tag = "seq+multichain only" if method == "off" else "+ SCOP2 fold union "
        summary["bridge_leakage"][method] = {}
        for label, key, held in holdouts:
            if not held:
                continue
            # No-protein (DNA/RNA-only) entries are not placed in a component; the
            # fold-bridge metric is over protein-bearing held-out entries.
            placed = [r for r in held if r.entry_id in cr.entry_to_cluster]
            bridged = sum(
                1 for r in placed if "train" in comp_splits[cr.entry_to_cluster[r.entry_id]]
            )
            print(
                f"  [{tag}] {label:11s} in a component with a TRAIN entry: "
                f"{bridged}/{len(placed)} ({pct(bridged, len(placed))})"
            )
            summary["bridge_leakage"][method][key] = {"bridged": bridged, "n": len(placed)}

    # ============================== F. test-set curation contamination (test only)
    print("\n" + "-" * 78)
    print("F. TEST-SET CONTAMINATION — IF-Split's ligand tiering on their class-stratified test")
    print("-" * 78)
    contamination = {}
    for klass in ("metal", "small_molecule", "nucleotide"):
        want = TEST_CLASS_TO_LIGAND[klass]
        members = [r for r in test if test_class_of[r.entry_id.upper()] == klass]
        good = [r for r in members if want in cls[r.entry_id.upper()]["classes"]]
        bad = [r for r in members if want not in cls[r.entry_id.upper()]["classes"]]
        reasons: Counter = Counter()
        for r in bad:
            reasons.update(v["reason"] for v in cls[r.entry_id.upper()]["tiers"].values())
        print(
            f"  {klass:15s} test n={len(members):4d}  functional={len(good):4d}  "
            f"NOT-functional={len(bad):3d} ({pct(len(bad), len(members))})"
        )
        if bad:
            print(f"      reasons: {', '.join(f'{k}:{v}' for k, v in reasons.most_common(6))}")
        contamination[klass] = {
            "n": len(members),
            "functional": len(good),
            "not_functional": len(bad),
            "not_functional_ids": sorted(r.entry_id for r in bad),
            "reasons": dict(reasons.most_common()),
        }
    summary["test_contamination"] = contamination

    print("\n" + "=" * 78)
    if out_path:
        with open(out_path, "w") as fh:
            json.dump(summary, fh, indent=2, sort_keys=True)
        print(f"summary written to {out_path}")


if __name__ == "__main__":
    main()
