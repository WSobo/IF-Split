"""IF-Split command-line interface.

`build` runs Stage 1 (enumerate via RCSB Search + Data API -> candidates.jsonl)
and writes the snapshot lock; `verify` re-enumerates from a lock and reports
drift. Later stages (filter -> ligands -> cluster -> split -> manifest) are
stubs that raise NotImplementedError; the CLI catches those and reports cleanly.
"""

from __future__ import annotations

import argparse
import sys

from pydantic import ValidationError

from . import __version__
from .config import load_config
from .rcsb import RcsbError

DEFAULT_CONFIG = "config/default.yaml"


def cmd_build(args: argparse.Namespace) -> int:
    from .enumerate import enumerate_candidates
    from .manifest import build_lock, write_lock

    cfg = load_config(args.config)
    sf = cfg.split_fractions
    assembly = "biological (assembly 1)" if cfg.use_biological_assembly else "asymmetric unit"
    print(f"IF-Split {__version__}")
    print(f"  config file:   {args.config}")
    print(f"  config hash:   {cfg.config_hash()}")
    print(f"  dataset:       {cfg.dataset_version}")
    print(f"  snapshot_date: {cfg.snapshot_date}  (selects release_date <= this)")
    print(f"  methods:       {', '.join(cfg.experimental_methods)}")
    print(f"  resolution:    <= {cfg.resolution_max_A} A")
    print(f"  max residues:  < {cfg.max_total_residues + 1}")
    print(f"  assembly:      {assembly}")
    print(f"  identity:      {cfg.identity_threshold}")
    print(f"  clustering:    {cfg.clustering_backend}")
    print(f"  splits:        train={sf.train} val={sf.val} test={sf.test}  salt={cfg.split_salt!r}")
    if args.limit is not None:
        print(f"  limit:         {args.limit} (dev: first N by sorted entry id)")
    print()

    print("Stage 1 - enumerate candidates (Search + Data API, no coordinates):")
    records, _candidates_path, sha = enumerate_candidates(
        cfg, args.out, limit=args.limit, progress=lambda m: print(f"  {m}")
    )

    lock = build_lock(
        cfg,
        entry_ids=[r.entry_id for r in records],
        candidates_sha256=sha,
        limit=args.limit,
    )
    lock_path = write_lock(lock, args.out)
    print(f"  wrote {lock_path}")
    print()
    print(f"Stage 1 complete: {len(records)} candidates. Stages 3-7 pending (see PLAN.md §8).")
    return 0


def cmd_verify(args: argparse.Namespace) -> int:
    from .manifest import verify_lock

    return verify_lock(args.lock)


def cmd_stats(args: argparse.Namespace) -> int:
    from .manifest import summarize_manifest

    return summarize_manifest(args.manifest)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="if-split",
        description="Reproducible, ligand-aware PDB train/val/test splitter (LigandMPNN-style).",
    )
    p.add_argument("--version", action="version", version=f"if-split {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    pb = sub.add_parser("build", help="Run the pipeline (Stages 1-7) and emit manifest + lock.")
    pb.add_argument(
        "--config",
        default=DEFAULT_CONFIG,
        help=f"Path to config YAML (default: {DEFAULT_CONFIG}).",
    )
    pb.add_argument(
        "--out",
        default="data/out",
        help="Output dir for candidates.jsonl + dataset.lock (default: data/out).",
    )
    pb.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Dev only: cap to the first N candidates by sorted entry id (reproducible).",
    )
    pb.set_defaults(func=cmd_build)

    pv = sub.add_parser(
        "verify", help="Re-enumerate from a lock file and report drift vs the live PDB."
    )
    pv.add_argument("lock", help="Path to dataset.lock")
    pv.set_defaults(func=cmd_verify)

    ps = sub.add_parser(
        "stats", help="Report split sizes and per-class test counts from a manifest."
    )
    ps.add_argument("manifest", help="Path to manifest.json")
    ps.set_defaults(func=cmd_stats)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except FileNotFoundError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    except ValidationError as e:
        print(f"invalid config:\n{e}", file=sys.stderr)
        return 2
    except NotImplementedError as e:
        print(f"not implemented yet: {e}", file=sys.stderr)
        return 3
    except RcsbError as e:
        print(f"RCSB request failed: {e}", file=sys.stderr)
        return 4


if __name__ == "__main__":
    sys.exit(main())
