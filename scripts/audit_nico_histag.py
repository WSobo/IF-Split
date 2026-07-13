"""Reproduce the lone-Ni/Co His-tag audit that justifies the metal-curation rule.

Run: ``uv run python scripts/audit_nico_histag.py`` (needs network + RCSB access).

The metal curation demotes a *lone* Ni/Co (an entry whose only metal is Ni or Co)
to ``ambiguous`` when nothing corroborates it. This script re-derives the numbers
that decision rests on, against the live RCSB snapshot, using the tool's own
``RcsbClient`` / ``CandidateRecord`` / ``has_histag`` plus an INDEPENDENT His-run
scanner (so the figure isn't just a re-run of the tool's own logic):

  1. What fraction of lone Ni/Co carry no detectable His-tag? (~82%, robust across
     tag definitions and between the independent scanner and the tool's has_histag.)
  2. Where a stray "96%" comes from: of the no-tag entries, ~96.7% are protein-bound
     — a DIFFERENT ratio than "fraction of lone Ni/Co with no tag".
  3. Blast radius: how many lone Ni/Co are actually re-tiered functional->ambiguous.
  4. Over-fire lower bound: how many demoted entries are recognizable Ni/Co enzymes
     (urease, Co-methionine-aminopeptidase, nitrile hydratase, ...) by title keyword.
"""

from __future__ import annotations

import re

import httpx

from ifsplit.ligands import has_histag, metal_comps
from ifsplit.rcsb import DATA_GRAPHQL_URL, SEARCH_URL, RcsbClient
from ifsplit.schema import CandidateRecord

# Snapshot filters mirror config/default.yaml.
METHODS = ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"]
RES_MAX = 3.5
CUTOFF = "2026-05-30T23:59:59Z"
PURIFICATION = {"NI", "CO"}

# Rough lower-bound keyword match for real Ni/Co metalloenzymes (not exhaustive).
ENZYME_PAT = re.compile(
    r"urease|hydrogenase|superoxide dismutase|glyoxalase|methionine aminopeptidase|"
    r"peptide deformylase|nitrile hydratase|xylose isomerase|phosphoglucomutase|"
    r"cobalamin|methyltransferase|nickel|cobalt",
    re.IGNORECASE,
)

TITLE_QUERY = (
    "query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id struct{title} "
    "polymer_entities{rcsb_polymer_entity{pdbx_description}}}}"
)


def _snapshot_nodes() -> list[dict]:
    return [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {"attribute": "exptl.method", "operator": "in", "value": METHODS},
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_entry_info.resolution_combined",
                "operator": "less_or_equal",
                "value": RES_MAX,
            },
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_accession_info.initial_release_date",
                "operator": "less_or_equal",
                "value": CUTOFF,
            },
        },
    ]


def _ligand_node(comp: str) -> dict:
    return {
        "type": "terminal",
        "service": "text_chem",
        "parameters": {
            "attribute": "rcsb_chem_comp_container_identifiers.comp_id",
            "operator": "exact_match",
            "value": comp,
        },
    }


def search_entries(http: httpx.Client, comp: str) -> list[str]:
    """Entry ids under the snapshot filters that contain nonpolymer component ``comp``."""
    ids: list[str] = []
    start = 0
    while True:
        body = {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [*_snapshot_nodes(), _ligand_node(comp)],
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {"start": start, "rows": 10000},
                "sort": [
                    {
                        "sort_by": "rcsb_entry_container_identifiers.entry_id",
                        "direction": "asc",
                    }
                ],
            },
        }
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


def max_his_run(seq: str) -> int:
    """Longest consecutive His run — independent of the tool's has_histag."""
    runs = [len(m.group(0)) for m in re.finditer(r"H+", (seq or "").upper())]
    return max(runs) if runs else 0


def fetch_titles(http: httpx.Client, ids: list[str]) -> dict[str, str]:
    """entry_id -> lowercased title + entity descriptions (for enzyme keyword match)."""
    text: dict[str, str] = {}
    for i in range(0, len(ids), 100):
        batch = ids[i : i + 100]
        resp = http.post(DATA_GRAPHQL_URL, json={"query": TITLE_QUERY, "variables": {"ids": batch}})
        if resp.status_code != 200:
            continue
        for e in (resp.json().get("data") or {}).get("entries") or []:
            if not e:
                continue
            parts = [(e.get("struct") or {}).get("title") or ""]
            for pe in e.get("polymer_entities") or []:
                parts.append((pe.get("rcsb_polymer_entity") or {}).get("pdbx_description") or "")
            text[e["rcsb_id"]] = " ".join(parts).lower()
    return text


def main() -> None:
    http = httpx.Client(timeout=90, headers={"User-Agent": "IF-Split-audit/0.1"})
    union = sorted(set(search_entries(http, "NI")) | set(search_entries(http, "CO")))
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

    lone = [r for r in records if (m := {c.comp_id for c in metal_comps(r)}) and m <= PURIFICATION]
    n = len(lone)
    print(f"lone Ni/Co entries (every metal in {{NI,CO}}): {n}")
    if not n:
        return
    prot = {r.entry_id: [e.seq for e in r.polymer_entities if e.is_protein] for r in lone}

    print("\n(1) fraction with NO detectable His-tag:")
    for k in (6, 5, 4):
        no = sum(1 for r in lone if all(max_his_run(s) < k for s in prot[r.entry_id]))
        print(f"    independent His-run >= {k}: {no}/{n} ({100 * no / n:.1f}% no-tag)")
    tool_no = [r for r in lone if not any(has_histag(s, 6, 3) for s in prot[r.entry_id])]
    print(
        f"    tool has_histag(6,term3) : {len(tool_no)}/{n} ({100 * len(tool_no) / n:.1f}% no-tag)"
    )

    bound = [r for r in tool_no if PURIFICATION & set(r.bound_components)]
    bound_pct = 100 * len(bound) / len(tool_no)
    print("\n(2) source of a stray '96%' — of the no-tag entries, how many are protein-bound:")
    print(f"    {len(bound)}/{len(tool_no)} ({bound_pct:.1f}%) — a DIFFERENT ratio")

    demoted = [
        r
        for r in tool_no
        if (PURIFICATION & set(r.bound_components))
        and not (PURIFICATION & set(r.affinity_comp_ids))
        and not (PURIFICATION & set(r.investigated_comp_ids))
    ]
    dem_pct = 100 * len(demoted) / n
    print(f"\n(3) re-tiered functional->ambiguous: {len(demoted)}/{n} ({dem_pct:.1f}% of lone)")

    titles = fetch_titles(http, [r.entry_id for r in demoted])
    real = [r for r in demoted if ENZYME_PAT.search(titles.get(r.entry_id, ""))]
    real_pct = 100 * len(real) / max(1, len(demoted))
    examples = [r.entry_id for r in real[:12]]
    print("\n(4) demoted entries that are recognizable Ni/Co enzymes (keyword lower bound):")
    print(f"    {len(real)}/{len(demoted)} ({real_pct:.1f}%)  examples: {examples}")
    http.close()


if __name__ == "__main__":
    main()
