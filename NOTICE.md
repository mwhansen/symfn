# Third-party notices

symfn is licensed under **MIT OR Apache-2.0** (see `LICENSE-MIT`,
`LICENSE-APACHE`). It contains no third-party code. This file records the
external projects symfn *relates* to, and why none of them changes that.

## lrcalc — clean-room separated

The Littlewood-Richardson Calculator (Anders S. Buch,
<https://bitbucket.org/asbuch/lrcalc>) is GPL-licensed.

An earlier version of symfn's LR engine was written after reading lrcalc. It was
not copied, and the ideas involved are textbook — but the author had *access* to
lrcalc's source, including implementation code that ships inline in its headers
(`lriter.h`). Since infringement is established by access plus substantial
similarity, that version was discarded rather than defended.

**The shipped `src/skew_lr.rs` was produced by a two-team clean room.** The
specification (`docs/cleanroom-spec-skew-lr.md`, committed as the audit trail)
was written by the party with lrcalc access and contains only textbook
mathematics, functional and performance requirements, and symfn's own interface.
It contains no implementation technique. The implementation was then written by
an independent party with no access to lrcalc, to the discarded version, or to
the project documents describing either.

The result is *not* a reconstruction of the discarded version. It is a row-level
dynamic program over a merged frontier — partial fillings that agree on the
preceding row and the content so far collapse into one weighted state — where
the discarded version enumerated individual tableaux. It is structurally closer
to symfn's own `src/strip_lr.rs` than to anything external, and it is
substantially faster than what it replaced. (Against lrcalc it is faster on some
shapes and slower on others; see `ROADMAP.md` for the measurements. Performance
is not the point here — independence is.)

Two findings from the exercise are worth recording:

- The independent implementation **did** arrive at the juxtaposition
  construction for products (place μ and ν in a skew diagram sharing no row or
  column, so the skew Schur function factors as the product). Converging on it
  without access is evidence it is an obvious consequence of the mathematics
  rather than anyone's expression.
- It **did not** arrive at the per-cell neighbour structure that was the
  discarded version's closest resemblance to lrcalc. That construct is simply
  absent here, because a row-level formulation has no use for it.

symfn shares no code with lrcalc, and neither links against nor bundles it.

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

## Symmetrica — public domain, and no relationship in code

symfn is a successor *in spirit* to Symmetrica. It shares no code with it and
derives nothing from it.

**Symmetrica is public domain**, stated by its authors at
<https://www.algorithm.uni-bayreuth.de/en/research/SYMMETRICA/>: "Symmetrica is
public domain." So unlike lrcalc above, there is no obligation here and never
was — reading, porting, or copying Symmetrica would all have been permitted.
The separation is a design decision (real types over a single untyped `OP`
object; see `src/lib.rs`), not a legal one, and it should not be mistaken for
the lrcalc clean room, which was necessary rather than chosen.

⚠️ The dedication lives on the project webpage, **not in the distribution**. A
search of the `Symmetrica_2.0` source tree in this repository for "public
domain", "copyright", "license", "GPL" and "warranty" — case-insensitive, all
142 files — returns nothing. Anyone receiving only the source has no evidence
of its status, so cite the webpage rather than the tree.

Scope, for reference: Symmetrica also covers modular and projective
representation theory of the symmetric group, Schubert polynomials (commutative
and non-commutative), Hecke algebras of type A, and finite group operations.
symfn implements a strict subset of it — the symmetric-function core — and
none of the rest.
