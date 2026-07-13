"""IF-Split command-line interface.

`build` runs the full pipeline: Stage 1 enumerate (RCSB Search + Data API ->
candidates.jsonl + dataset.lock), Stage 3 filter, Stage 4 ligand classification,
Stage 5 cluster, Stage 6 deterministic split, Stage 7 manifest + registry. No
structure coordinates are downloaded. `verify` re-enumerates from a lock and
reports drift; `stats` summarizes a manifest.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pydantic import ValidationError

from . import __version__
from .config import load_config
from .rcsb import RcsbError

DEFAULT_CONFIG = "config/default.yaml"
SPLITS_CHOICES = ("train", "val", "test")


def cmd_build(args: argparse.Namespace) -> int:
    from .cluster import build_clusters
    from .enumerate import enumerate_candidates, make_console_progress
    from .ligands import classify_components
    from .manifest import (
        build_lock,
        build_manifest,
        build_tiers_doc,
        read_registry,
        write_classes,
        write_clusters,
        write_lock,
        write_manifest,
        write_registry,
        write_split_files,
        write_tiers,
    )
    from .parse import drop_summary, filter_candidates
    from .split import assign_splits, check_no_leakage

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
    say = make_console_progress()  # timestamped + line-flushed (survives redirect)
    records, _candidates_path, sha = enumerate_candidates(
        cfg, args.out, limit=args.limit, progress=say
    )
    lock_path = write_lock(
        build_lock(
            cfg,
            entry_ids=[r.entry_id for r in records],
            candidates_sha256=sha,
            limit=args.limit,
        ),
        args.out,
    )
    print(f"  wrote {lock_path}")

    print("Stage 3 - filter (metadata only):")
    kept, drops = filter_candidates(records, cfg)
    dcounts = drop_summary(drops)
    print(f"  kept {len(kept)} / {len(records)}; dropped {len(drops)} {dcounts or ''}")

    print("Stage 4 - ligand classification + curation:")
    class_map = {r.entry_id: classify_components(r, cfg) for r in kept}
    n_artifact = sum(1 for v in class_map.values() if v["purification_artifact"])
    print(f"  classified {len(class_map)} entries; {n_artifact} purification artifact(s) flagged")

    print(f"Stage 5 - cluster ({cfg.clustering_backend} @ {cfg.identity_level}%):")
    clusters = build_clusters(kept, cfg)
    print(
        f"  {clusters.n_clusters} leakage-safe components "
        f"from {clusters.n_raw_clusters} raw clusters "
        f"({len(clusters.multichain_entries)} multi-chain merged)"
    )

    print("Stage 6 - assign splits (deterministic hash):")
    registry = read_registry(args.registry) if args.registry else {}
    entry_classes = {eid: info["classes"] for eid, info in class_map.items()}
    splits = assign_splits(clusters, cfg, registry=registry, entry_classes=entry_classes)
    check_no_leakage(splits, clusters)  # structural guarantee; raises on violation
    c = splits.counts
    print(f"  train={c['train']} val={c['val']} test={c['test']}  (no cross-split leakage)")
    if cfg.test_min_per_class:
        if splits.minimum_shortfalls:
            short = ", ".join(f"{k}:{v}" for k, v in splits.minimum_shortfalls.items())
            print(f"  test minimums: applied; SHORTFALL (not enough supply) -> {short}")
        else:
            print("  test minimums: applied; all per-class floors met")

    print("Stage 7 - manifest + registry:")
    manifest = build_manifest(
        cfg,
        candidates_sha256=sha,
        n_candidates=len(records),
        drops=drops,
        drop_counts=dcounts,
        clusters=clusters,
        splits=splits,
        class_map=class_map,
    )
    mpath = write_manifest(manifest, args.out)
    split_paths = write_split_files(splits, class_map, args.out)
    write_clusters(clusters.entry_to_cluster, args.out)
    write_classes(class_map, args.out)
    write_registry(splits.cluster_split, args.out)
    write_tiers(build_tiers_doc(class_map), args.out)
    print(f"  wrote {mpath} (provenance + counts)")
    for s in ("train", "val", "test"):
        print(f"  wrote {split_paths[s]}")
    test_class_paths = [p for k, p in split_paths.items() if k.startswith("test:")]
    for p in sorted(test_class_paths):
        print(f"  wrote {p}")
    print("  wrote clusters.json, ligands.classes.json, ligands.tiers.json, splits.registry.json")
    print()
    print(f"Build complete: {len(kept)} structures across train/val/test.")
    print(f"Run `if-split stats {mpath}`.")
    return 0


def cmd_verify(args: argparse.Namespace) -> int:
    from .manifest import verify_lock

    return verify_lock(args.lock)


def cmd_stats(args: argparse.Namespace) -> int:
    from .manifest import summarize_manifest

    return summarize_manifest(args.manifest)


def cmd_spec(args: argparse.Namespace) -> int:
    """Emit a portable, self-identifying split spec (YAML) from a build or config.

    The source may be a manifest.json (the config is embedded in it) or an existing
    config YAML. The output is a small file you can share so anyone can reproduce
    your split with `if-split build --config <that file>`.
    """
    import yaml

    from .config import Config, SpecMeta
    from .manifest import read_manifest

    src = Path(args.source)
    if src.name.endswith(".json"):
        cfg = Config.model_validate(read_manifest(src)["config"])
    else:
        cfg = load_config(src)

    # Apply optional human metadata overrides from flags.
    meta = (cfg.spec or SpecMeta()).model_dump(exclude_none=True)
    for field in ("name", "description", "author"):
        val = getattr(args, field, None)
        if val is not None:
            meta[field] = val
    cfg = cfg.model_copy(update={"spec": SpecMeta(**meta)})

    doc = cfg.to_spec_dict(stamp_hash=True)
    text = yaml.safe_dump(doc, sort_keys=False, default_flow_style=False)
    if args.out:
        Path(args.out).write_text(text, encoding="utf-8")
        print(f"wrote split spec -> {args.out}  (config_hash {cfg.config_hash()})")
    else:
        sys.stdout.write(text)
    return 0


# Confirm before a pull larger than this many structures unless --yes is given.
_FETCH_CONFIRM_THRESHOLD = 1000


def cmd_fetch(args: argparse.Namespace) -> int:
    from .download import SPLITS, StructureFetcher
    from .hydrate import hydrate, select_targets
    from .manifest import read_manifest

    if args.all and args.split:
        print("error: use either --all or --split, not both", file=sys.stderr)
        return 2
    if not args.all and not args.split:
        print(
            "error: choose a scope: --split test (repeatable) or --all.\n"
            "       (explicit by design — `fetch` can pull a lot of data)",
            file=sys.stderr,
        )
        return 2

    splits = list(SPLITS) if args.all else args.split
    unknown = [s for s in splits if s not in SPLITS]
    if unknown:
        print(f"error: unknown split(s): {', '.join(unknown)}", file=sys.stderr)
        return 2

    manifest = read_manifest(args.manifest)
    # Split id-lists (train.json etc.) live next to the manifest, NOT in the cwd —
    # resolve them relative to the manifest so the documented invocation works from
    # any directory (hydrate() below already does this for the real fetch).
    targets = select_targets(manifest, splits, base_dir=Path(args.manifest).parent)
    if not targets:
        print(f"nothing to fetch for split(s): {', '.join(splits)}")
        return 0

    assembly = not args.asymmetric_unit
    print(f"fetch: {len(targets)} structures across {', '.join(splits)} -> {args.out}")

    # Size estimate + large-pull guard (no accidental terabyte).
    unit_label = "assembly 1" if assembly else "asymmetric unit"
    with StructureFetcher(assembly=assembly, workers=args.workers) as fetcher:
        est = fetcher.estimate_bytes([e for e, _ in targets])
        if est is not None:
            print(f"  estimated download: ~{est / 1e9:.2f} GB ({unit_label})")
        if len(targets) > _FETCH_CONFIRM_THRESHOLD and not args.yes:
            print(
                f"  refusing to fetch {len(targets)} structures without --yes "
                f"(threshold {_FETCH_CONFIRM_THRESHOLD}).",
                file=sys.stderr,
            )
            return 5
        summary = hydrate(
            args.manifest,
            args.out,
            splits=splits,
            assembly=assembly,
            workers=args.workers,
            fetcher=fetcher,
            progress=lambda m: print(f"  {m}"),
        )

    print(
        f"done: {summary['fetched']} fetched, {summary['skipped']} cached, "
        f"{len(summary['failed'])} failed"
    )
    if summary["failed"]:
        for eid, reason in summary["failed"][:10]:
            print(f"  ! {eid}: {reason}", file=sys.stderr)
    for kind, path in summary["index"].items():
        print(f"  index ({kind}): {path}")
    print(f"  dataset card: {args.out}/DATASET_CARD.md")
    return 0 if not summary["failed"] else 6


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
    pb.add_argument(
        "--registry",
        default=None,
        help="Optional prior splits.registry.json to pin existing cluster->split "
        "assignments (growth stability).",
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

    psp = sub.add_parser(
        "spec",
        help="Emit a portable, shareable split spec (YAML) from a manifest or config.",
    )
    psp.add_argument("source", help="Path to a manifest.json or a config YAML.")
    psp.add_argument("--out", default=None, help="Write spec here (default: stdout).")
    psp.add_argument("--name", default=None, help="Human-readable split name.")
    psp.add_argument("--description", default=None, help="One-line description.")
    psp.add_argument("--author", default=None, help="Author/attribution.")
    psp.set_defaults(func=cmd_spec)

    pf = sub.add_parser(
        "fetch",
        help="OPTIONAL: download structures for a built manifest into an ML-ready tree.",
    )
    pf.add_argument("manifest", help="Path to manifest.json")
    pf.add_argument(
        "--out", default="data/structures", help="Output root (default: data/structures)."
    )
    pf.add_argument(
        "--split",
        action="append",
        choices=list(SPLITS_CHOICES),
        help="Split to fetch (repeatable, e.g. --split test --split val).",
    )
    pf.add_argument("--all", action="store_true", help="Fetch all splits.")
    pf.add_argument(
        "--asymmetric-unit",
        action="store_true",
        help="Fetch the asymmetric unit instead of biological assembly 1.",
    )
    pf.add_argument("--workers", type=int, default=8, help="Concurrent downloads (default: 8).")
    pf.add_argument("--yes", action="store_true", help="Proceed without confirming large pulls.")
    pf.set_defaults(func=cmd_fetch)

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
