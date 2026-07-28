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

## Bignum coefficients — permissive, and deliberately not GMP

The `bignum` feature (which `python` enables) pulls in **num-bigint**,
**num-rational** and **num-traits**, all dual **MIT OR Apache-2.0** — the same
terms as symfn. Nothing in the dependency graph is copyleft, so the published
wheel carries no relinking obligation and can be redistributed under symfn's own
licence.

GMP, via `rug`, was the obvious alternative and was rejected on two grounds.
Licensing: GMP is dual LGPLv3+ / GPLv2+ and `rug` is LGPL-3.0+, so statically
linking either into a distributed wheel attaches LGPL obligations to a binary
this project intends to ship under MIT/Apache-2.0. Performance: it would not
have bought anything. Coefficients past `i128` in this library are 2-5 limbs
(`examples/coeff_sizes.rs`), where every bignum implementation runs schoolbook
and GMP's Karatsuba/Toom/FFT paths never engage, and coefficient arithmetic is
only a few percent of runtime to begin with. Exactness was the requirement;
speed was not.

## Sage — test oracle only

`scripts/gen_sage_oracle.sage` invokes Sage (GPL) as a separate program to
generate `tests/fixtures/sage_oracle.txt`, on the same terms as lrcalc above.

## Symmetrica — public domain, and consulted once

symfn is a successor *in spirit* to Symmetrica. It shares no code with it.

It is not, however, wholly independent of it, and this is the one place that
matters: **Symmetrica's source was read to identify the algorithm behind
`m → s`.** `tms.c` computes m_μ in the Schur basis as the product m_μ · s_∅ and
routes it to `muir.c`, which names the rule — Muir's. Knowing the name was the
useful part; the implementation in `src/convert.rs` was then derived from the
bialternant identity rather than transcribed, and it differs in substance (see
`muir_expand`: β-numbers in a u64 mask, slots processed in decreasing value so
that a collision is provably permanent and prunes on the spot). Symmetrica's
`muir_lim_new` is a different construction over explicit vectors, and the
readable form preserved in `mms.c`'s `#ifdef UNDEF` block enumerates subsets and
permutations directly.

This is recorded because it is true, not because anything requires it —
Symmetrica is public domain, and copying it outright would have been permitted.

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
