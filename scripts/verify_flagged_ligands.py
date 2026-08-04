"""Check IF-Split's flagged ligand-class test entries against the deposited coordinates.

Stage 4's tiering is metadata-only by design, so it decides "is this ligand real?"
from RCSB's ``nonpolymer_bound_components`` (bond-based), its curated
subject-of-investigation flag, measured affinities, an additive blacklist and the
CCD carbohydrate type. Every one of those is a proxy. This script tests the proxies
the only way that settles it -- by measuring the distance from the flagged component
to the protein in the deposited model.

It is a *validation* script, not part of the build path: `build` still never
downloads coordinates (see CLAUDE.md). Only the handful of entries the tiering
already flagged are fetched here.

For each flagged entry it re-runs the tiering, then for every component that did NOT
come back ``functional`` it reports the closest protein contact and the coordinating
residues. A metal within 2.8 A of protein donor atoms is coordinated; an organic
component with tens of contacts under 4 A sits in a pocket. Either way the metadata
call of "unbound" is wrong, and the script says so.

Coordinating residues that fall inside a terminal poly-His run of the entity sequence
are marked ``[tag]``, which is what separates a genuine IMAC artifact (Ni held by the
purification tag) from a native surface site.

Usage:
  uv run python scripts/verify_flagged_ligands.py                    # audit's flagged ids
  uv run python scripts/verify_flagged_ligands.py 3I9Z 2B4L          # explicit ids
  uv run python scripts/verify_flagged_ligands.py --json out.json
  # the archived record in examples/ligandmpnn-audit/ is regenerated with:
  uv run python scripts/verify_flagged_ligands.py --all-comps \\
      --json examples/ligandmpnn-audit/flagged_ligand_verification.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import gemmi
import httpx

from ifsplit import __version__
from ifsplit.config import load_config
from ifsplit.ligands import (
    TIER_FUNCTIONAL,
    classify_components,
    is_metal_ion,
    terminal_residue_run,
)
from ifsplit.rcsb import RcsbClient
from ifsplit.schema import CandidateRecord

# The ids IF-Split's tiering flags in LigandMPNN's three ligand-class test sets, as
# reported by scripts/audit_ligandmpnn_split.py (examples/ligandmpnn-audit/).
DEFAULT_IDS = [
    # metal test set
    "1T31",
    "2CFV",
    "2NZ6",
    "3HG9",
    "3I9Z",
    "4X68",
    # small-molecule test set
    "2B4L",
    "3UEU",
    "4GNY",
    "5YFS",
    "5YFT",
    "6I67",
    # nucleotide test set
    "2ZIO",
    "3OLT",
]

CIF_URL = "https://files.rcsb.org/download/{pid}.cif.gz"
# A metal-donor bond is ~1.8-2.6 A; 2.8 A is a generous coordination cutoff. Organic
# ligands bind non-covalently, so contact (van der Waals) distance is the right test.
METAL_COORD_CUTOFF = 2.8
LIGAND_CONTACT_CUTOFF = 4.0
# A His run at least this long, within a terminal window, is tag-like -- the same
# rule ligands.py uses to detect the tag in the first place.
TAG_MIN_RUN = 4
TAG_WINDOW = 20


def fetch_cif(pid: str, cache: Path) -> gemmi.Structure:
    """Download (and cache) one deposited mmCIF."""
    cache.mkdir(parents=True, exist_ok=True)
    path = cache / f"{pid.lower()}.cif.gz"
    if not path.exists():
        with httpx.Client(timeout=120, follow_redirects=True) as client:
            resp = client.get(CIF_URL.format(pid=pid.upper()))
            resp.raise_for_status()
            path.write_bytes(resp.content)
    with gzip.open(path, "rt") as fh:
        st = gemmi.read_structure_string(fh.read())
    st.setup_entities()
    st.remove_hydrogens()
    return st


def tag_spans(seq: str) -> list[tuple[int, int]]:
    """1-based [start, end] spans of terminal poly-His runs in an entity sequence.

    Only terminal runs count: a purification tag is terminal by construction, and an
    internal His cluster is more likely a real metal site (ligands.py makes the same
    distinction).
    """
    if terminal_residue_run(seq, "H", TAG_WINDOW) < TAG_MIN_RUN:
        return []
    spans, run_start = [], None
    for i, ch in enumerate(seq.upper(), start=1):
        if ch == "H":
            run_start = i if run_start is None else run_start
            continue
        if run_start is not None:
            if i - run_start >= TAG_MIN_RUN:
                spans.append((run_start, i - 1))
            run_start = None
    if run_start is not None and len(seq) - run_start + 1 >= TAG_MIN_RUN:
        spans.append((run_start, len(seq)))
    # Keep only the runs that are actually near a terminus.
    return [(a, b) for a, b in spans if a <= TAG_WINDOW or b > len(seq) - TAG_WINDOW]


def contacts(st: gemmi.Structure, comp_id: str, cutoff: float) -> list[dict]:
    """Every protein atom within ``cutoff`` of any instance of ``comp_id``.

    Symmetry is dropped first, so this measures the *deposited* asymmetric unit only.
    That is the conservative direction: a contact found here cannot be an artifact of
    generating the biological assembly, so it refutes an "unbound" call outright.
    """
    st.cell = gemmi.UnitCell()
    st.spacegroup_hm = "P 1"
    model = st[0]
    search = gemmi.NeighborSearch(model, st.cell, cutoff).populate()
    found = []
    for chain in model:
        for res in chain:
            if res.name != comp_id:
                continue
            for lig_atom in res:
                for mark in search.find_atoms(lig_atom.pos, "\0", radius=cutoff):
                    cra = mark.to_cra(model)
                    if cra.residue.het_flag != "A":  # keep standard polymer residues
                        continue
                    d = lig_atom.pos.dist(cra.atom.pos)
                    if d <= cutoff:
                        found.append(
                            {
                                "ligand_instance": f"{chain.name}/{res.seqid.num}",
                                "ligand_atom": lig_atom.name,
                                "chain": cra.chain.name,
                                "res": f"{cra.residue.name}{cra.residue.seqid.num}",
                                "label_seq": cra.residue.label_seq,
                                "atom": cra.atom.name,
                                "distance": round(d, 2),
                            }
                        )
    return sorted(found, key=lambda c: c["distance"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("ids", nargs="*", default=None, help="entry ids (default: the audit's)")
    ap.add_argument("--config", default="config/default.yaml")
    ap.add_argument("--cache", type=Path, default=Path("data/cache/cif"))
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument(
        "--all-comps",
        action="store_true",
        help="measure every component, not only the flagged ones. Use for an archived "
        "record: a gate-order change moves components in and out of the flagged set, "
        "so a flags-only file silently loses the evidence for whatever got rescued.",
    )
    args = ap.parse_args()

    ids = [i.upper() for i in (args.ids or DEFAULT_IDS)]
    cfg = load_config(args.config)

    client = RcsbClient()
    try:
        records = {
            CandidateRecord.from_data_api(raw).entry_id.upper(): CandidateRecord.from_data_api(raw)
            for raw in client.fetch_entries(ids)
        }
    finally:
        client.close()

    report: dict[str, dict] = {}
    for pid in ids:
        rec = records.get(pid)
        if rec is None:
            print(f"{pid}: UNRESOLVED via RCSB")
            continue
        tiers = classify_components(rec, cfg)["tiers"]
        flagged = {c: v for c, v in tiers.items() if v["tier"] != TIER_FUNCTIONAL}
        measured = tiers if args.all_comps else flagged
        if not flagged:
            print(f"{pid}: nothing flagged")
            if not args.all_comps:
                continue

        st = fetch_cif(pid, args.cache)
        # label_seq_id indexes into the entity sequence, so a tag span maps straight
        # onto the modeled residues without needing author numbering.
        tag_positions: set[int] = set()
        for entity in st.entities:
            seq = gemmi.one_letter_code(entity.full_sequence).upper()
            for start, end in tag_spans(seq):
                tag_positions.update(range(start, end + 1))

        by_comp = {c.comp_id: c for c in rec.nonpolymer_comps}
        print("=" * 78)
        summary = ", ".join(f"{c} ({v['reason']})" for c, v in sorted(measured.items()))
        print(f"{pid}  {summary}")
        entry_out = {}
        for comp_id, verdict in sorted(measured.items()):
            comp = by_comp.get(comp_id)
            metal = comp is not None and is_metal_ion(comp)
            cutoff = METAL_COORD_CUTOFF if metal else LIGAND_CONTACT_CUTOFF
            found = contacts(st, comp_id, cutoff)
            residues, seen = [], set()
            for c in found:
                key = (c["ligand_instance"], c["chain"], c["res"])
                if key in seen:
                    continue
                seen.add(key)
                tag = " [tag]" if c["label_seq"] in tag_positions else ""
                residues.append(f"{c['res']}.{c['atom']} {c['distance']}A{tag}")
            kind = "coordination" if metal else "contact"
            print(f"  {comp_id:5s} {kind}s <= {cutoff}A: {len(found)}")
            if residues:
                print(f"        {'; '.join(residues[:8])}")
            else:
                print("        NONE -- no protein atom within cutoff")
            entry_out[comp_id] = {
                "reason": verdict["reason"],
                "tier": verdict["tier"],
                "is_metal": metal,
                "cutoff": cutoff,
                "n_contacts": len(found),
                "closest": found[0] if found else None,
                "contact_residues": residues,
                "any_tag_residue": any("[tag]" in r for r in residues),
            }
        report[pid] = entry_out

    if args.json:
        # Stamp the tiering code that produced the `tier`/`reason` fields. They are the
        # one part of this file that is not a property of the PDB: a Stage 4 gate-order
        # change rewrites them, while the measured distances never move.
        out = {
            "_provenance": {
                "ifsplit_version": __version__,
                "config": str(args.config),
                "all_comps": args.all_comps,
                "metal_coord_cutoff_angstrom": METAL_COORD_CUTOFF,
                "ligand_contact_cutoff_angstrom": LIGAND_CONTACT_CUTOFF,
                "note": "distances are from the deposited asymmetric unit, symmetry dropped",
            },
            "entries": report,
        }
        args.json.write_text(json.dumps(out, indent=1, sort_keys=True) + "\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    sys.exit(main())
