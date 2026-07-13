"""Reproducible benchmark for the lone-Ni/Co metal tiering, over the live PDB.

Runs the tool's *own* Stage 1 fetch (now carrying RCSB metal annotations) and
Stage 4 classify_components over every lone Ni/Co entry in the snapshot, and
reports the distribution of tier reasons. Re-run after any change to the metal
rule to see, at scale, how many entries are:

  - rescued to functional (metal_annotated / metal_affinity / metal_investigated /
    metal_bound),
  - reported ambiguous as a non-native metal site (metal_site_nonnative), or
  - demoted with no corroboration (purification_metal_uncorroborated / metal_unbound).

This is the labelled eval scaffolding: the STRONG rescues and the metal_site_
nonnative bucket are the two levers to tune, and this makes any change measurable
instead of asserted. Run: ``uv run python scripts/eval_metal_tiering.py``.
"""

from __future__ import annotations

from collections import Counter

import httpx

from ifsplit.config import load_config
from ifsplit.ligands import classify_components, metal_comps
from ifsplit.rcsb import SEARCH_URL, RcsbClient
from ifsplit.schema import CandidateRecord

METHODS = ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"]
RES_MAX = 3.5
CUTOFF = "2026-05-30T23:59:59Z"
PURIFICATION = {"NI", "CO"}


def _snapshot_nodes() -> list[dict]:
    return [
        {"type": "terminal", "service": "text",
         "parameters": {"attribute": "exptl.method", "operator": "in", "value": METHODS}},
        {"type": "terminal", "service": "text",
         "parameters": {"attribute": "rcsb_entry_info.resolution_combined",
                        "operator": "less_or_equal", "value": RES_MAX}},
        {"type": "terminal", "service": "text",
         "parameters": {"attribute": "rcsb_accession_info.initial_release_date",
                        "operator": "less_or_equal", "value": CUTOFF}},
    ]  # fmt: skip


def search_entries(http: httpx.Client, comp: str) -> list[str]:
    ids: list[str] = []
    start = 0
    while True:
        node = {"type": "terminal", "service": "text_chem",
                "parameters": {"attribute": "rcsb_chem_comp_container_identifiers.comp_id",
                               "operator": "exact_match", "value": comp}}  # fmt: skip
        nodes = [*_snapshot_nodes(), node]
        body = {
            "query": {"type": "group", "logical_operator": "and", "nodes": nodes},
            "return_type": "entry",
            "request_options": {"paginate": {"start": start, "rows": 10000},
                                "sort": [{"sort_by": "rcsb_entry_container_identifiers.entry_id",
                                          "direction": "asc"}]},
        }  # fmt: skip
        resp = http.post(SEARCH_URL, json=body)
        if resp.status_code == 204:
            break
        resp.raise_for_status()
        page = resp.json().get("result_set", [])
        if not page:
            break
        ids.extend(hit["identifier"] for hit in page)
        start += len(page)
        if len(page) < 10000:
            break
    return ids


def main() -> None:
    cfg = load_config("config/default.yaml")
    http = httpx.Client(timeout=90, headers={"User-Agent": "IF-Split-audit/0.1"})
    union = sorted(set(search_entries(http, "NI")) | set(search_entries(http, "CO")))
    http.close()
    print(f"snapshot Ni/Co-containing entries: {len(union)}")

    rc = RcsbClient()
    records: list[CandidateRecord] = []
    try:
        for i, raw in enumerate(rc.fetch_entries(union), 1):
            records.append(CandidateRecord.from_data_api(raw))
            if i % 1000 == 0:
                print(f"  enriched {i}/{len(union)}")
    finally:
        rc.close()

    reasons: Counter[str] = Counter()
    examples: dict[str, list[str]] = {}
    n_lone = 0
    for rec in records:
        metals = {c.comp_id for c in metal_comps(rec)}
        if not (metals and metals <= PURIFICATION):
            continue
        n_lone += 1
        res = classify_components(rec, cfg)
        for cid in metals:
            reason = res["tiers"].get(cid, {}).get("reason", "?")
            reasons[reason] += 1
            examples.setdefault(reason, [])
            if len(examples[reason]) < 8:
                examples[reason].append(rec.entry_id)

    functional = {"metal_annotated", "metal_affinity", "metal_investigated", "metal_bound"}
    print(f"\nlone Ni/Co entries: {n_lone}\ntier reason distribution (per Ni/Co component):")
    for reason, n in reasons.most_common():
        kind = "functional" if reason in functional else "ambiguous/artifact"
        print(f"  {reason:34s} {n:5d}  [{kind}]  e.g. {examples[reason][:6]}")
    resc = sum(reasons[r] for r in functional)
    print(f"\n  rescued to functional: {resc}  ({100 * resc / max(1, sum(reasons.values())):.1f}%)")
    print(f"  metal_site_nonnative : {reasons['metal_site_nonnative']}")
    print(f"  uncorroborated       : {reasons['purification_metal_uncorroborated']}")


if __name__ == "__main__":
    main()
