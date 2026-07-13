"""Reproducible benchmark for the small-molecule / glycan tiering, over the live PDB.

Runs the tool's own Stage 1 fetch + Stage 4 classify over a representative stride
sample of the whole snapshot and reports how non-protein small molecules tier:

  - functional small-molecule targets by reason (ligand_bound / ligand_investigated /
    ligand_affinity),
  - carbohydrates demoted to `glycan` (decorative glycosylation / sugar-detergents,
    excluded from the functional corpus unless a measured affinity vouches for them),
  - additives (blacklist) and unbound.

Re-run after any change to the small-molecule rule to see, at scale, the corpus
cleanliness. Run: ``uv run python scripts/eval_sm_tiering.py``.
"""

from __future__ import annotations

from collections import Counter

from ifsplit.config import load_config
from ifsplit.ligands import classify_components
from ifsplit.rcsb import RcsbClient
from ifsplit.schema import CandidateRecord

N_SAMPLE = 3000


def main() -> None:
    cfg = load_config("config/default.yaml")
    client = RcsbClient()
    try:
        all_ids = client.search_entry_ids(cfg)
        print(f"snapshot: {len(all_ids)} entries")
        step = len(all_ids) / N_SAMPLE
        sample = sorted({all_ids[int(i * step)] for i in range(min(N_SAMPLE, len(all_ids)))})
        print(f"stride-sampling {len(sample)}; fetching...")
        records: list[CandidateRecord] = []
        for i, raw in enumerate(client.fetch_entries(sample), 1):
            records.append(CandidateRecord.from_data_api(raw))
            if i % 1000 == 0:
                print(f"  enriched {i}/{len(sample)}")
    finally:
        client.close()

    func_reason: Counter[str] = Counter()
    glycan_demoted: Counter[str] = Counter()
    n_additive = 0
    func_comps: Counter[str] = Counter()
    for rec in records:
        res = classify_components(rec, cfg)
        for comp in res.get("small_molecules", []):
            func_reason[res["tiers"][comp]["reason"]] += 1
            func_comps[comp] += 1
        for comp, t in res["tiers"].items():
            if t["reason"] == "glycan":
                glycan_demoted[comp] += 1
            elif t["reason"] == "additive":
                n_additive += 1

    n_func = sum(func_reason.values())
    print(f"\nrecords: {len(records)}")
    print(f"functional small_molecule targets: {n_func}")
    for r, n in func_reason.most_common():
        print(f"  {r:22s} {n:5d}  ({100 * n / max(1, n_func):4.1f}%)")
    n_gly = sum(glycan_demoted.values())
    print(f"\ncarbohydrates demoted to glycan (not a target unless affinity): {n_gly}")
    print(f"  top glycans: {glycan_demoted.most_common(12)}")
    print(f"additive (blacklist) comp instances: {n_additive}")
    print(f"\ntop functional small_molecule comps: {func_comps.most_common(20)}")


if __name__ == "__main__":
    main()
