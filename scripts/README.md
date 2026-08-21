# scripts

- **`gen_lrcalc_oracle.py`** — regenerate the lrcalc oracle fixture:
  ```
  LRCALC=/path/to/lrcalc python3 gen_lrcalc_oracle.py > ../tests/fixtures/lrcalc_oracle.txt
  ```
  lrcalc is invoked as an external program and only its output is used; see
  `../NOTICE.md`. It complements the Sage oracle by reaching much larger
  shapes — that is where the single-traversal LR backend differs most from the
  coefficient-at-a-time one.

- **`check_deltaop.py`** — hold `src/deltaop.rs` to Sage and to the identities
  Sage cannot check:
  ```
  cargo run --release --example deltaop_dump -- 7 > /tmp/dop.txt
  sage -python check_deltaop.py /tmp/dop.txt
  ```
  Sage has `nabla` and nothing else, so ∇ and ∇² are compared against it
  directly while Δ', Θ are tied to it through
  `Θ_{e_k}∇e_{n−k} = Δ'_{e_{n−k−1}}e_n` and `Δ'_{e_{n−1}}e_n = ∇e_n`, plus the
  rise version of the Delta conjecture against a labelled-Dyck-path enumeration.

- **`verify_deltaop_formulas.py`** — verify every formula in
  `docs/record/macdonald-operators.md` numerically, before any of it is
  implemented:
  ```
  sage -python verify_deltaop_formulas.py
  ```
  Not a `check_*.py`: there is no symfn dump to compare against yet. It
  implements ∇, Δ_f, Δ'_f, Π and Θ_f from the papers on top of Sage's `Ht` basis
  and checks them against Sage's own `nabla`, against the published
  `Θ_{e_k}∇e_{n−k} = Δ'_{e_{n−k−1}}e_n` identity, and against a direct
  enumeration of labelled Dyck paths for both versions of the Delta conjecture.
  Prints `all formulas verified` and exits nonzero otherwise.

- **`check_jack.py`** — hold `src/jack.rs` and `src/afrac.rs` to Sage, and to
  the laws Sage cannot state:
  ```
  cargo run --release --example jack_dump -- 7 5 > /tmp/jack.txt
  sage -python check_jack.py /tmp/jack.txt
  ```
  Unlike the Delta-operator check, Sage has all three normalizations, so `P`,
  `Q`, `J`, `J → p`, the norms and the three inverse directions (`m` written
  in each normalization) are a real external oracle. Three things have
  no oracle and are checked as laws instead: [KS] Thm 1.1 (`[m_μ]J_λ ∈ ℕ[α]`
  and divisible by `u_μ`), the closed-form norms against `scalar_jack`, and
  **Stanley's open conjecture** — where a violation is a result to report, not
  a bug to fix. Two Sage warts are re-checked here rather than assumed:
  `scalar_jack` cannot *return* zero, and `zonal()` is `P^{(2)}`, not `J^{(2)}`.

- **`bench_jack.py`** — the Sage side of the Jack ladder, so both sides can be
  measured in one session and one power state:
  ```
  sage -python -u bench_jack.py 11
  ```
  ⚠️ Record the power state. This machine drifts about 1.8× on battery, and
  `docs/record/jack.md`'s own sweep straddled an AC detach mid-run.

- **`spec_jack_verify.py`** / **`spec_jack_walls.py`** / **`spec_jack_swell.py`**
  — the pre-implementation trio behind `docs/record/jack.md`, in the
  `verify_deltaop_formulas.py` role: run before any Rust existed.
  `verify` checks all three routes (Laplace–Beltrami, branching, Knop–Sahi)
  against Sage and prints `all formulas verified`; `walls` measures where Sage
  stops; `swell` simulates the proposed factored-atom arithmetic through the
  recursion and reports that there is **no denominator swell at all** — the
  measurement that licensed the design.

- **`gen_sage_oracle.sage`** — regenerate the Sage oracle fixture:
  ```
  SAGE_DISABLE_SYMFN=1 sage gen_sage_oracle.sage > ../tests/fixtures/sage_oracle.txt
  ```
  ⚠️ The variable is not optional and the script refuses without it. Sage
  reaches this library through its optional backend, so a fixture taken with it
  enabled is symfn quoting itself. The header of the script is the line-format
  key; every family's forward expansion is there and, since 2026-08-21, every
  inverse one too (`sinhlp`, `sinhlqp`, `sinht`, `sinj`, `minp`, `minq`,
  `jminp`, `jminq`, `jminj`).
- **`bench_vs_sage.py`** — benchmark symfn against Sage's own symmetric
  functions. Requires the extension module built (see repo README):
  ```
  sage -python bench_vs_sage.py
  ```
  Note: Sage memoizes symmetric-function products, so benchmarks that repeat an
  *identical* computation measure its cache, not its algorithm. Use distinct
  inputs computed once, as this script does.

- **`doc_review.py`** — needs neither Sage nor a build. Walks every rustdoc
  paragraph in `src/` one at a time and collects a comment on each, so a prose
  pass can be reviewed in one sitting rather than across forty files:
  ```
  python3 doc_review.py review --skip-tests     # one paragraph per screen
  python3 doc_review.py collect --skip-tests    # only the ones you commented on
  python3 doc_review.py files                   # paragraph counts, to scope a pass
  ```
  Enter means *nothing to say*, and that is **recorded as an answer** rather
  than as silence: a resumed pass starts at the first paragraph never put on
  screen, not at the first one without a comment. `s` skips without answering,
  for the ones worth a second look. Every answer is written before the next
  paragraph is drawn, so Ctrl-C costs at most the one in front of you.

  `extract` writes the same worksheet (`../doc-review.txt`) for editing by
  hand, which suits a long sitting with a lot to say; the two modes share a
  file and can be alternated. Answers survive a re-`extract`: a paragraph is
  matched on its text and the item it documents, never its line number, so
  unrelated edits above it do not detach the note. Answers whose paragraph
  *did* change are listed at the top of the regenerated file and nowhere else.
  `--path`, `--module-only` and `--min-words` scope it down; the worksheet is
  gitignored.

- **`bench_macdonald_cache.py`** — Sage's Macdonald `J` change of basis with the
  backend against without it, one degree at a time:
  ```
  python -u bench_macdonald_cache.py 11 2
  ```
  The two arms alternate in separate processes and the control gets
  `SAGE_DISABLE_SYMFN` in its *environment*, which is the only place it works:
  Sage fills its conversion table at import. ⚠️ Record the power state.

- **`bench_inverse.py`** — time the seven inverse expansions against Sage —
  Macdonald `s -> Htilde`, `s -> J`, `m -> P`, `m -> Q` and Jack `m -> P`,
  `m -> Q`, `m -> J` — one degree and one arm per process:
  ```
  python bench_inverse.py 8
  ```
  The unit is the degree — every partition of `n` expanded out of the basis
  the family's forward direction returns, Schur for the two `s ->` arms and
  monomial for the rest — which is what both sides amortize a transition
  matrix over. One process **per arm** is not optional: Sage shares a family's
  transition matrix between its normalizations, and two Jack arms in one
  process read 20x faster than they do alone. The two families print as two
  tables and run over different base rings. The Sage
  arm gets `SAGE_DISABLE_SYMFN` in its *environment* and the script **refuses**
  without it taking effect: Sage's Macdonald bases reach this library through
  the optional backend, and an arm with it enabled reports about 1.0x. Drives
  `cargo run --release --example bench_inverse` for the symfn side, which
  clears the memo between workloads. ⚠️ Record the power state.

- **`check_python_pointers.py`** — fail when a docstring or `#:` comment under
  `python/symfn/`, or the `///` on a `#[pyfunction]` or the `#[pymodule]` in
  `src/python.rs`, names a `docs/`, `scripts/`, `examples/` or `*.rs` path,
  which a reader of the wheel cannot open (docs/policies/python.md, P11).
  Reads source; needs nothing built.

- **`check_python_stubs.py`** — hold `symfn.pyi` to the module it describes:
  ```
  cargo build --features python
  python3 check_python_stubs.py ../target/debug/libsymfn.dylib
  ```
  Needs **no Sage**, like its siblings `check_python_boundary.py` — which
  feeds every entry point malformed input and demands a typed exception — and
  `check_python_marshalling.py` below. This one is what makes the stub file
  *be* the supported surface rather than describe it: it fails when a name is exported
  without a stub, stubbed without being exported, or when a parameter name
  differs between the two. That last case is the quiet one — PyO3 exports every
  argument as keyword-callable, so a parameter name is contract, and a stub
  saying `mu` where the module says `nu` type-checks a call that fails.

- **`check_python_marshalling.py`** — the round-trip half of the boundary
  suite (docs/policies/python.md, P1):
  ```
  cargo build --features python
  python3 check_python_marshalling.py ../target/debug/libsymfn.dylib
  ```
  Every exported callable runs once on a small valid input and its return is
  validated against the shape the stub file promises, `type() is` strict;
  coefficients at both edges of the `i128` fast path and past it round-trip
  through identity-shaped calls; and the permissive inbound spellings — list
  or tuple, padded or not, an output handed straight back — agree. Needs no
  Sage. Its first run caught two defects at `i128::MIN`
  (docs/record/python-and-sage-interop.md).
