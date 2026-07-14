"""Reference example: go from an IF-Split target to its pocket in the coordinates.

IF-Split's job ends at the *labels*: leakage-safe train/val/test lists and a
per-ligand conditioning-target corpus (`targets.jsonl`). It never parses
coordinates. This script shows the intended last mile a researcher owns -- join a
target row to a fetched structure with *their own* parser and pull the ligand
atoms + pocket residues -- so an inverse-folding model can condition on exactly
the right ligand.

It is deliberately NOT part of the package: featurization is model-specific
(ligand atoms? SMILES? a pocket mask? all-atom context?), so this is a template
to copy and adapt, not an API to depend on.

Prereqs:
  1. Build a split:   uv run if-split build --out data/out
  2. Fetch structures: uv run if-split fetch data/out/manifest.json --split test
  3. A structure parser. This example uses gemmi (`uv pip install gemmi`); biotite
     or Biopython work equally well -- the pattern is the same.

Run:
  uv run python scripts/consume_split.py \
      --manifest data/out/manifest.json --structures data/structures --split test --limit 5
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

from ifsplit.dataset import ConditioningTarget, load_dataset
from ifsplit.download import rel_path_for
from ifsplit.manifest import read_manifest

try:
    import gemmi
except ImportError:  # pragma: no cover - example-only dependency
    gemmi = None


def _load_structure(path: Path):
    """Parse a (gzipped) mmCIF with gemmi. Swap in your own parser here."""
    with gzip.open(path, "rt") as fh:
        doc = gemmi.cif.read_string(fh.read())
    st = gemmi.make_structure_from_block(doc.sole_block())
    st.setup_entities()
    return st


def _ligand_copies(model, comp_id: str) -> list:
    """Every residue in the model whose component id matches -- i.e. all copies.

    This is the crux of the multi-ligand / multi-copy case: a structure may hold
    several copies of one ligand (and unrelated adventitious ions). We surface
    *all* instances; which one to condition on is the caller's per-instance choice.
    """
    copies = []
    for chain in model:
        for res in chain:
            if res.name == comp_id:
                copies.append((chain.name, res.seqid.num, res))
    return copies


def _pocket_residues(model, ligand_res, radius: float) -> list[tuple[str, int, str]]:
    """Polymer residues with any atom within ``radius`` Å of the ligand's atoms."""
    lig_atoms = [atom.pos for atom in ligand_res]
    pocket: dict[tuple[str, int], str] = {}
    for chain in model:
        for res in chain:
            if res.name == ligand_res.name and res.seqid.num == ligand_res.seqid.num:
                continue  # skip the ligand itself
            info = gemmi.find_tabulated_residue(res.name)
            if not (info and info.is_amino_acid()):
                continue  # keep only standard protein residues (drops water, HET, NA)
            for atom in res:
                if any(atom.pos.dist(p) <= radius for p in lig_atoms):
                    pocket[(chain.name, res.seqid.num)] = res.name
                    break
    return [(c, n, name) for (c, n), name in sorted(pocket.items())]


def _radius_from_manifest(manifest_path: str) -> float:
    cfg = read_manifest(manifest_path).get("config", {})
    return float(cfg.get("ligand_context_radius_A", 8.0))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True, help="Path to a built manifest.json")
    ap.add_argument("--structures", required=True, help="Root passed to `fetch --out`")
    ap.add_argument("--split", default="test", choices=("train", "val", "test"))
    ap.add_argument("--limit", type=int, default=5, help="How many entries to demo")
    ap.add_argument("--radius", type=float, default=None, help="Å (default: manifest config)")
    ap.add_argument(
        "--asymmetric-unit",
        dest="assembly",
        action="store_false",
        help="Read the AU instead of assembly 1 (match your `fetch` choice).",
    )
    args = ap.parse_args()

    if gemmi is None:
        raise SystemExit(
            "This example needs gemmi: `uv pip install gemmi` (or adapt to your own parser)."
        )

    radius = args.radius if args.radius is not None else _radius_from_manifest(args.manifest)
    ds = load_dataset(args.manifest)
    view = ds.split(args.split)
    root = Path(args.structures)

    print(
        f"{args.split}: {len(view)} backbones, "
        f"{len(view.conditioning_targets())} functional targets; pocket radius {radius} Å\n"
    )

    shown = 0
    # Group targets by structure so we open each mmCIF once (condition-on-all view).
    for entry_id, targets in view.targets_by_entry().items():
        if shown >= args.limit:
            break
        path = root / rel_path_for(entry_id, args.split, assembly=args.assembly)
        if not path.exists():
            print(f"{entry_id}: not fetched ({path}); run `if-split fetch` first -- skipping")
            continue

        model = _load_structure(path)[0]
        print(f"{entry_id}  ({len(targets)} target(s))")
        for t in targets:
            _describe_target(model, t, radius)
        print()
        shown += 1
    return 0


def _describe_target(model, t: ConditioningTarget, radius: float) -> None:
    if t.comp_id is None:
        # A nucleic_acid target is a protein<->DNA/RNA interface, not a HET group.
        # There's nothing to look up by comp_id: condition on the nucleic-acid
        # polymer chain(s) in the model (select DNA/RNA chains directly).
        print(f"  - {t.ligand_class:14s} [{t.tier}] nucleic-acid chain (select DNA/RNA polymer)")
        return
    copies = _ligand_copies(model, t.comp_id)
    if not copies:
        print(f"  - {t.comp_id:14s} [{t.tier}] not found in assembly coordinates")
        return
    for chain_name, seqid, res in copies:
        pocket = _pocket_residues(model, res, radius)
        print(
            f"  - {t.comp_id:6s} copy {chain_name}/{seqid:<5d} [{t.tier}, {t.reason}] "
            f"-> {len(pocket)} pocket residues within {radius} Å"
        )


if __name__ == "__main__":
    raise SystemExit(main())
