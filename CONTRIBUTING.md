# Contributing

Bugs and questions go to
[the issue tracker](https://github.com/mwhansen/symfn/issues); a wrong value
has its own issue template, because it is the most serious defect this
library can have. A vulnerability goes to the private route in
[SECURITY.md](SECURITY.md), which also says what does and does not count as
one here. [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) applies everywhere the
project is discussed.

## Which checks need what

- `cargo test` needs nothing beyond the Rust toolchain. The oracle tests
  (`tests/sage_oracle.rs`, `tests/lrcalc_oracle.rs`) read committed fixtures
  under `tests/fixtures/`, so they run offline.
- `scripts/preflight.sh` is the gate for a Rust change: the format check,
  both test suites, and the prose and link checks. It needs the toolchain
  and python3, no network and no Sage. Run it before every commit; the
  pre-commit hook re-checks only formatting.
- `scripts/preflight_python.sh` is the gate for anything under
  `src/python.rs`, `python/symfn/` or `docsite/` — and for the README's
  Python example, which runs in its page-examples step. It needs
  `cargo build --features python`, ruff, mypy and Sphinx.
- **Regenerating a fixture** needs the oracle that produced it — Sage or
  lrcalc — and nearly every script in `scripts/` needs Sage.
  [scripts/README.md](scripts/README.md) documents the main ones, including
  the release and packaging scripts.

## Once per clone

```sh
git config core.hooksPath .githooks
git config blame.ignoreRevsFile .git-blame-ignore-revs
```

## Where things are

- [docs/layout.md](docs/layout.md) maps every module in the tree.
- [docs/public-api.md](docs/public-api.md) maps the public surface and its
  tiers.
- [docs/record/](docs/record/) is the long-term memory, one file per
  subsystem. Before working in a subsystem, read its record file — dead
  ends are recorded with their premises exactly so they are not re-explored
  at full price.
- [CLAUDE.md](CLAUDE.md) is the session bootstrap for coding agents. It
  routes to the five rulebooks that govern prose, failure handling, the
  Python surface, validation, and the record; the rules bind human
  contributors the same way.
- [docs/release-readiness.md](docs/release-readiness.md) is the release
  plan.
