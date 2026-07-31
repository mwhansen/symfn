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
  `Q`, `J`, `J → p` and the norms are a real external oracle. Three things have
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
  sage gen_sage_oracle.sage > ../tests/fixtures/sage_oracle.txt
  ```
- **`bench_vs_sage.py`** — benchmark symfn against Sage's own symmetric
  functions. Requires the extension module built (see repo README):
  ```
  sage -python bench_vs_sage.py
  ```
  Note: Sage memoizes symmetric-function products, so benchmarks that repeat an
  *identical* computation measure its cache, not its algorithm. Use distinct
  inputs computed once, as this script does.
