"""IF-Split command-line interface.

`build` runs the full pipeline: Stage 1 enumerate (RCSB Search + Data API ->
candidates.jsonl + dataset.lock), Stage 3 filter, Stage 4 ligand classification,
Stage 5 cluster, Stage 6 deterministic split, Stage 7 manifest + registry. No
structure coordinates are downloaded. `resplit` re-derives Stages 3-7 from a
cached candidates.jsonl with no RCSB access — for ablating curation/clustering/
split settings on a fixed snapshot without re-enumerating. `verify` re-derives
from a lock and reports drift (against the live PDB, or offline against a local
candidates.jsonl via --candidates); `stats` summarizes a manifest.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import httpx
from pydantic import ValidationError

from . import __version__
from .config import load_config
from .rcsb import RcsbError

DEFAULT_CONFIG = "config/default.yaml"
SPLITS_CHOICES = ("train", "val", "test")

# Documented process exit codes (see `if-split --help` and the README).
EXIT_OK = 0
EXIT_BAD_INPUT = 2  # missing/malformed file, invalid config, or bad argument value
EXIT_NOT_IMPLEMENTED = 3
EXIT_NETWORK = 4  # RCSB / HTTP failure after retries
EXIT_INTERRUPTED = 130


def _print_config_header(cfg, config_path: str, *, limit: int | None = None) -> None:
    sf = cfg.split_fractions
    assembly = "biological (assembly 1)" if cfg.use_biological_assembly else "asymmetric unit"
    print(f"IF-Split {__version__}")
    print(f"  config file:   {config_path}")
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
    if limit is not None:
        print(f"  limit:         {limit} (dev: first N by sorted entry id)")
    print()


def _resolve_registry(cfg, out, registry_path, fresh) -> dict[str, str]:
    """Resolve the growth-stability registry for this build.

    Precedence: an explicit ``--registry`` path wins; ``--fresh`` forces a clean
    (unpinned) lineage; otherwise, for the ``balanced`` strategy only, auto-adopt
    ``<out>/splits.registry.json`` when a prior in-place build used the SAME config
    (its ``dataset.lock`` ``config_hash`` matches). This makes an in-place rebuild
    growth-stable by default — the fix for ``balanced``, whose fill boundaries shift
    as the snapshot grows unless prior components are pinned. ``hash`` is already
    input-independent and registry-free (so ``verify`` can still certify it), so it
    is never auto-pinned.
    """
    from .manifest import read_lock, read_registry

    if registry_path:
        return read_registry(registry_path)
    if fresh or cfg.split_strategy != "balanced":
        return {}
    reg_file = Path(out) / "splits.registry.json"
    if not reg_file.exists():
        return {}  # first build into this dir — nothing to pin
    same_config = False
    lock_file = Path(out) / "dataset.lock"
    if lock_file.exists():
        try:
            lock = read_lock(lock_file)
            same_config = isinstance(lock, dict) and lock.get("config_hash") == cfg.config_hash()
        except (OSError, ValueError):
            same_config = False
    if same_config:
        reg = read_registry(reg_file)
        if reg:
            print(
                f"  growth stability: pinning {len(reg)} prior assignments from "
                f"{reg_file.name} (same config; --fresh to start a new lineage)"
            )
        return reg
    print(
        f"  growth stability: {reg_file.name} exists but its build config differs — "
        f"starting a NEW lineage (a balanced split is only growth-stable within one "
        f"config; pass --registry <prior> to pin, or --fresh to silence)"
    )
    return {}


def _run_pipeline(
    cfg, records, sha, out, *, limit, registry_path, fresh=False, source="build"
) -> tuple[Path, int]:
    """Stages 3-7 from an in-memory candidate set (shared by ``build`` and ``resplit``).

    Filters, classifies ligands, clusters, assigns the split, checks the no-leakage
    invariant, and writes the full output tree (manifest, lock, split files, targets,
    tiers, classes, registry). Returns ``(manifest_path, n_kept)``. Coordinate-free.
    """
    from .cluster import build_clusters
    from .ligands import classify_components
    from .manifest import (
        build_fold_benchmark,
        build_lock,
        build_manifest,
        build_targets,
        build_tiers_doc,
        write_classes,
        write_clusters,
        write_fold_benchmark,
        write_lock,
        write_manifest,
        write_registry,
        write_split_files,
        write_targets,
        write_tiers,
    )
    from .parse import drop_summary, filter_candidates
    from .split import (
        assign_splits,
        check_no_leakage,
        registry_fingerprint,
        split_fingerprint,
    )

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
    registry = _resolve_registry(cfg, out, registry_path, fresh)
    growth_stable = cfg.split_strategy != "balanced" or bool(registry)
    entry_classes = {eid: info["classes"] for eid, info in class_map.items()}
    splits = assign_splits(clusters, cfg, registry=registry, entry_classes=entry_classes)
    check_no_leakage(splits, clusters)  # structural guarantee; raises on violation
    c = splits.counts
    print(f"  train={c['train']} val={c['val']} test={c['test']}  (no cross-split leakage)")
    if splits.strategy != "hash":
        note = f"  strategy={splits.strategy}: {splits.capped_folds} dominant folds -> train"
        if splits.balance_gaps:
            note += f"; WARNING fold tail too thin, val/test short by {splits.balance_gaps}"
        print(note)
    if cfg.test_min_per_class:
        if splits.minimum_shortfalls:
            short = ", ".join(f"{k}:{v}" for k, v in splits.minimum_shortfalls.items())
            print(f"  test minimums: applied; SHORTFALL (not enough supply) -> {short}")
        else:
            print("  test minimums: applied; all per-class floors met")

    fold_benchmark = build_fold_benchmark(
        clusters.entry_fold_labels, splits.entry_split, cfg.fold_benchmark_method
    )

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
        growth_stable=growth_stable,
        fold_benchmark_summary=fold_benchmark["summary"] if fold_benchmark else None,
    )
    mpath = write_manifest(manifest, out)
    # The lock is written here (after Stage 6), not at Stage 1, so it can pin the
    # split OUTPUT (split_sha256) alongside the candidate inputs — verify then
    # certifies the split reproduced, not just the candidates.
    lock_path = write_lock(
        build_lock(
            cfg,
            entry_ids=[r.entry_id for r in records],
            candidates_sha256=sha,
            limit=limit,
            split_sha256=split_fingerprint(splits.entry_split),
            registry_sha256=registry_fingerprint(registry),
            split_strategy=cfg.split_strategy,
            source=source,
        ),
        out,
    )
    print(f"  wrote {lock_path}")
    split_paths = write_split_files(splits, class_map, out)
    write_clusters(clusters.entry_to_cluster, out)
    write_classes(class_map, out)
    write_registry(splits.cluster_split, out)
    write_tiers(build_tiers_doc(class_map), out)
    targets = build_targets(class_map, splits, clusters)
    write_targets(targets, out)
    write_fold_benchmark(fold_benchmark, out)  # writes sidecars, or clears stale ones
    n_targets = sum(1 for t in targets if t["tier"] == "functional")
    print(f"  wrote {mpath} (provenance + counts)")
    for s in ("train", "val", "test"):
        print(f"  wrote {split_paths[s]}")
    for p in sorted(p for k, p in split_paths.items() if k.startswith("test:")):
        print(f"  wrote {p}")
    print("  wrote clusters.json, ligands.classes.json, ligands.tiers.json, splits.registry.json")
    print(f"  wrote targets.jsonl ({len(kept)} backbones, {n_targets} conditioning targets)")
    if fold_benchmark is not None:
        n_novel = fold_benchmark["summary"]["n_test_novel_fold"]
        print(
            f"  wrote novel_fold_test.json, folds.json, fold_groups.json "
            f"({n_novel} novel-fold test entries)"
        )
    return mpath, len(kept)


def cmd_build(args: argparse.Namespace) -> int:
    from .enumerate import enumerate_candidates, make_console_progress

    cfg = load_config(args.config)
    _print_config_header(cfg, args.config, limit=args.limit)
    if args.count:
        from .rcsb import RcsbClient

        with RcsbClient() as client:
            n = client.count_entries(cfg)
        capped = f"; a build would cap to the first {args.limit}" if args.limit else ""
        print(f"Search matches {n} entries for this snapshot{capped}.")
        print("(--count is a preview — no candidates were fetched or written.)")
        return 0
    print("Stage 1 - enumerate candidates (Search + Data API, no coordinates):")
    say = make_console_progress()  # timestamped + line-flushed (survives redirect)
    records, _candidates_path, sha = enumerate_candidates(
        cfg, args.out, limit=args.limit, progress=say
    )
    mpath, n_kept = _run_pipeline(
        cfg,
        records,
        sha,
        args.out,
        limit=args.limit,
        registry_path=args.registry,
        fresh=getattr(args, "fresh", False),
    )
    print()
    print(f"Build complete: {n_kept} structures across train/val/test.")
    print(f"Run `if-split stats {mpath}`.")
    return 0


def _warn_if_config_would_reenumerate(cfg, candidates_path, sha: str) -> None:
    """Warn if this resplit config would enumerate a SUPERSET of the cached snapshot.

    ``resplit`` re-derives Stages 3-7 from a *fixed* candidate set: it can tighten or
    change curation, clustering, and the split, but it cannot add entries a looser
    Stage-1 snapshot (later cutoff, more methods, higher resolution) would include. If a
    ``dataset.lock`` sits beside the candidates file and describes it (matching hash),
    compare the Stage-1-affecting fields and warn on a widening change; otherwise print
    a generic caveat.
    """
    from .config import Config
    from .manifest import read_lock

    lock_path = Path(candidates_path).parent / "dataset.lock"
    locked = None
    if lock_path.exists():
        try:
            lock = read_lock(lock_path)
            # A malformed lock (valid JSON of the wrong shape) must degrade to the
            # generic caveat, never crash this purely-informational guard.
            if (
                isinstance(lock, dict)
                and isinstance(lock.get("candidates"), dict)
                and lock["candidates"].get("sha256") == sha
            ):
                locked = Config.model_validate(lock["config"])
        except (KeyError, ValueError, TypeError, AttributeError):
            locked = None
    if locked is None:
        print(
            "  note: resplit re-derives from this fixed candidate snapshot; it cannot add "
            "entries a looser config (later snapshot_date, more methods, higher resolution) "
            "would enumerate — run `build` for those."
        )
        return
    widened = []
    if cfg.snapshot_date > locked.snapshot_date:
        widened.append(f"snapshot_date {locked.snapshot_date} -> {cfg.snapshot_date}")
    if set(cfg.experimental_methods) - set(locked.experimental_methods):
        widened.append("experimental_methods added")
    if cfg.search_resolution_cap() > locked.search_resolution_cap():
        widened.append(
            f"resolution {locked.search_resolution_cap()} -> {cfg.search_resolution_cap()} A"
        )
    if widened:
        print(
            f"  WARNING: this config would enumerate MORE than the cached snapshot "
            f"({'; '.join(widened)}), so the resplit result is INCOMPLETE for it. "
            f"Run `build` to enumerate the full set; resplit only re-derives/tightens an "
            f"existing snapshot."
        )


def cmd_resplit(args: argparse.Namespace) -> int:
    from .schema import read_candidates_jsonl, sha256_hex

    cfg = load_config(args.config)
    cand = Path(args.candidates)
    if not cand.exists():
        print(f"error: candidates file not found: {cand}", file=sys.stderr)
        return 2
    _print_config_header(cfg, args.config)
    print("Stage 1 - SKIPPED (offline resplit: re-deriving Stages 3-7 from cached candidates):")
    sha = sha256_hex(cand.read_bytes())
    try:
        records = read_candidates_jsonl(cand)
    except (ValueError, OSError) as exc:
        print(f"error: {cand} is not a valid candidates.jsonl ({exc})", file=sys.stderr)
        return 2
    print(f"  read {len(records)} candidates from {cand} (sha256={sha[:12]}...)")
    _warn_if_config_would_reenumerate(cfg, cand, sha)
    mpath, n_kept = _run_pipeline(
        cfg,
        records,
        sha,
        args.out,
        limit=None,
        registry_path=args.registry,
        fresh=getattr(args, "fresh", False),
        source="resplit",
    )
    print()
    print(f"Resplit complete: {n_kept} structures across train/val/test (from a fixed snapshot).")
    print(f"Run `if-split stats {mpath}`.")
    return 0


def cmd_verify(args: argparse.Namespace) -> int:
    from .manifest import verify_lock

    return verify_lock(args.lock, candidates_path=args.candidates)


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

    if args.workers < 1:
        print(f"error: --workers must be >= 1, got {args.workers}", file=sys.stderr)
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
        "--count",
        action="store_true",
        help="Preview only: print how many entries the snapshot filters match "
        "(one fast Search API call) and exit, without fetching or writing anything.",
    )
    pb.add_argument(
        "--registry",
        default=None,
        help="Optional prior splits.registry.json to pin existing cluster->split "
        "assignments (growth stability).",
    )
    pb.add_argument(
        "--fresh",
        action="store_true",
        help="Start a new split lineage: ignore any splits.registry.json already in "
        "--out. By default a balanced rebuild into the same --out auto-pins the prior "
        "assignments when the config matches (growth stability); --fresh opts out.",
    )
    pb.set_defaults(func=cmd_build)

    prs = sub.add_parser(
        "resplit",
        help="Re-derive the split (Stages 3-7) from a cached candidates.jsonl — no RCSB.",
    )
    prs.add_argument(
        "--candidates",
        required=True,
        help="Path to a candidates.jsonl from a prior build (the fixed snapshot).",
    )
    prs.add_argument(
        "--config", default=DEFAULT_CONFIG, help=f"Config YAML (default: {DEFAULT_CONFIG})."
    )
    prs.add_argument("--out", default="data/out", help="Output dir (default: data/out).")
    prs.add_argument(
        "--registry",
        default=None,
        help="Optional prior splits.registry.json to pin existing cluster->split assignments.",
    )
    prs.add_argument(
        "--fresh",
        action="store_true",
        help="Start a new split lineage: ignore any splits.registry.json already in --out.",
    )
    prs.set_defaults(func=cmd_resplit)

    pv = sub.add_parser(
        "verify", help="Re-derive from a lock and report drift; --candidates verifies offline."
    )
    pv.add_argument("lock", help="Path to dataset.lock")
    pv.add_argument(
        "--candidates",
        default=None,
        help="Verify offline against this local candidates.jsonl instead of re-enumerating "
        "from RCSB (checks it matches the lock, then re-derives the split).",
    )
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
    except KeyboardInterrupt:
        print("\ninterrupted.", file=sys.stderr)
        return EXIT_INTERRUPTED
    except FileNotFoundError as e:
        print(f"error: file not found: {e}", file=sys.stderr)
        return EXIT_BAD_INPUT
    except ValidationError as e:
        print(f"invalid config:\n{e}", file=sys.stderr)
        return EXIT_BAD_INPUT
    except json.JSONDecodeError as e:
        # JSONDecodeError subclasses ValueError, so catch it first for a precise message.
        print(f"error: malformed JSON ({e})", file=sys.stderr)
        return EXIT_BAD_INPUT
    except KeyError as e:
        print(
            f"error: malformed or old-schema file — missing expected field {e}. "
            "Was it produced by a different if-split version?",
            file=sys.stderr,
        )
        return EXIT_BAD_INPUT
    except NotImplementedError as e:
        print(f"not implemented: {e}", file=sys.stderr)
        return EXIT_NOT_IMPLEMENTED
    except RcsbError as e:
        print(f"RCSB request failed: {e}", file=sys.stderr)
        return EXIT_NETWORK
    except httpx.HTTPError as e:
        print(f"network error: {e}", file=sys.stderr)
        return EXIT_NETWORK
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return EXIT_BAD_INPUT


if __name__ == "__main__":
    sys.exit(main())
