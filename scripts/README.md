# scripts

- **`gen_lrcalc_oracle.py`** — regenerate the lrcalc oracle fixture:
  ```
  LRCALC=/path/to/lrcalc python3 gen_lrcalc_oracle.py > ../tests/fixtures/lrcalc_oracle.txt
  ```
  lrcalc is invoked as an external program and only its output is used; see
  `../NOTICE.md`. It complements the Sage oracle by reaching much larger
  shapes — that is where the single-traversal LR backend differs most from the
  coefficient-at-a-time one.

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
