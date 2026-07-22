# Contributing to IF-Split

Thanks for your interest in improving IF-Split! Contributions of all kinds —
bug reports, feature requests, documentation, and code — are welcome.

## Reporting bugs and requesting features

Please open an issue on the [GitHub issue tracker](https://github.com/WSobo/IF-Split/issues).
For bugs, include:

- the `if-split` version (`if-split --version`) and your Python version,
- the exact command or code you ran and the full error/traceback,
- what you expected to happen instead.

Because a split is defined by a config, attaching the `config.yaml` (or the
`spec:` block from a `*.ifsplit.yaml`) makes most issues reproducible.

## Getting help

Open a [GitHub Discussion or issue](https://github.com/WSobo/IF-Split/issues)
with the `question` label. There is no separate support channel.

## Development setup

IF-Split uses [`uv`](https://docs.astral.sh/uv/):

```bash
git clone https://github.com/WSobo/IF-Split && cd IF-Split
uv sync                      # create .venv from uv.lock (deps + dev tools)
uv run pytest                # run the test suite (offline by default)
uv run ruff check .          # lint (must pass)
uv run ruff format .         # format (line length 100)
```

The test suite is offline by default. One opt-in live-network round-trip test runs
only with `IFSPLIT_NETWORK_TESTS=1`.

## Submitting changes

1. Fork the repo and create a topic branch off `main`.
2. Make your change, adding or updating tests. Two invariants are load-bearing and
   have guarding tests — please keep them green:
   - **Determinism:** the same config produces a byte-identical `manifest.json`.
   - **No cross-split leakage:** a sequence/fold component maps to exactly one split.
3. Ensure `uv run pytest` and `uv run ruff check .` pass, and run `uv run ruff format .`.
4. Open a pull request describing the change and the motivation. CI (GitHub Actions)
   must pass.

## Code conventions

- Python ≥ 3.11, `ruff` for lint + format (line length 100).
- Keep the core `build` path **metadata-only** — it must never download structure
  coordinates (see [PLAN.md](PLAN.md) §1.5). Coordinate access belongs to the
  optional `fetch` path.
- Preserve PDB identifiers verbatim (never slice/case-fold entry or entity ids).

By contributing, you agree that your contributions are licensed under the project's
[MIT License](LICENSE).
