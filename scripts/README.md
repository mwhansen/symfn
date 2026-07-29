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
  `../docs/spec-macdonald-operators.md` numerically, before any of it is
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
