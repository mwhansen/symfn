# Third-party notices

symfn is licensed under **MIT OR Apache-2.0** (see `LICENSE-MIT`,
`LICENSE-APACHE`). It contains no third-party code. This file records the
external projects symfn *relates* to, and why none of them changes that.

## lrcalc — studied, not copied

The Littlewood-Richardson Calculator (Anders S. Buch,
<https://bitbucket.org/asbuch/lrcalc>) is GPL-licensed. symfn's LR engine
(`src/skew_lr.rs`) uses ideas that lrcalc also uses:

- expanding a skew Schur function in a **single traversal**, binning tableaux by
  the content they turn out to have, rather than recomputing per candidate;
- using the running **content vector as the ballot test** (place value `v` only
  while `cont[v] < cont[v-1]`), so the lattice-word condition and content
  validity are one check;
- **normalizing the shape** before counting — here, transposing when that
  reduces the row count, since rows drive branching.

These are published mathematical facts and algorithmic techniques, not
expression: they appear in the LR literature and in other implementations.
symfn's implementation is independent Rust, written against those ideas; no
lrcalc source was translated, adapted, or copied, and symfn neither links
against nor bundles lrcalc.

lrcalc is used in one place only, as an **external test oracle**:
`scripts/gen_lrcalc_oracle.py` invokes the `lrcalc` binary as a separate
program to generate `tests/fixtures/lrcalc_oracle.txt`. Running a GPL program
and recording its output imposes no license obligation on symfn, and the
committed fixture means the test suite needs no lrcalc installed. The fixture
holds computed mathematical values — LR coefficients are facts, not authorship.

## GMP / rug — LGPL, and optional

The `gmp` feature pulls in `rug`, which links GMP. Both are **LGPL** (GMP is
dual LGPLv3+ / GPLv2+; `rug` is LGPL-3.0+). LGPL permits use from a
differently-licensed program, so this does **not** require symfn to be GPL and
does not affect symfn's own MIT/Apache-2.0 terms.

The feature is also **off by default** — the stock build has zero dependencies
and does not link GMP at all. Downstream users who enable `gmp` and ship
binaries take on the usual LGPL relinking obligation for that dependency; that
obligation is theirs and attaches to GMP, not to symfn.

## Sage — test oracle only

`scripts/gen_sage_oracle.sage` invokes Sage (GPL) as a separate program to
generate `tests/fixtures/sage_oracle.txt`, on the same terms as lrcalc above.

## Symmetrica — no relationship in code

symfn is a clean-room successor *in spirit* to Symmetrica. It shares no code
with it and derives nothing from it.
