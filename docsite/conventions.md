# Conventions

Every family in this library has rivals in the literature that differ by a
twist — a `q ↔ t` swap, a `t → 1/t`, an `α → 1/α`, a transpose of the indices.
A wrong convention does not raise. It returns a well-formed answer to a question
you did not ask, and it will agree with the right one on the small cases you are
most likely to check by hand.

So each family's documentation states its convention and gives a value that
distinguishes it, and the specializations below are the ones to run when you
need to be sure which library you are talking to.

## The classical bases

`s`, `h`, `e`, `p`, `m`, `f` are the Schur, complete homogeneous, elementary,
power-sum, monomial and forgotten bases, in the standard normalizations. The
Hall inner product makes Schur orthonormal:

```pycon
>>> from symfn import s, p
>>> s([2, 1]).scalar(s([2, 1])), s([2, 1]).scalar(s([3]))
(1, 0)
>>> p([2]).scalar(p([2]))
2
```

The second value is what tells you the power-sum basis is orthogonal but not
orthonormal: `⟨p_λ, p_λ⟩ = z_λ`.

The `ω` involution conjugates the shape, and the antipode is `ω` with a sign:

```pycon
>>> s([2, 1, 1]).omega()
s[3,1]
>>> s([2, 1]).antipode()
-s[2,1]
```

## Macdonald

`P` is monic in the monomial basis; `Q = b_λ P`; `J = c_λ P` is the integral
form, whose coefficients are polynomials rather than fractions. `Htilde` is the
modified form `H̃_μ`, whose Schur coefficients are the `(q,t)`-Kostka
polynomials in the Garsia–Haiman orientation.

The specialization that pins the whole family is `P_λ(x; q, q) = s_λ`:

```pycon
>>> from symfn import s, macdonald
>>> macdonald.P([3, 1]).at(q=5, t=5) == s([3, 1]).to("m")
True
```

The smallest value separating `P` from `Q`:

```pycon
>>> macdonald.P([1]).coefficient([1]), macdonald.Q([1]).coefficient([1])
(1, (1 - t)/(1 - q))
```

And the orientation of `K̃`:

```pycon
>>> macdonald.qt_kostka([2], [1, 1])
t
>>> macdonald.Htilde([2])
q*s[1,1] + s[2]
```

Under the `q ↔ t` mirror both of those swap, which is what makes them the
values worth checking.

## Jack

`P` is monic in the monomial basis. The parameter is `α` in the convention
where `α = 1` gives Schur functions and `α = 2` gives zonal polynomials — not
the `α → 1/α` mirror:

```pycon
>>> from symfn import s, jack
>>> jack.P([3, 1]).at(alpha=1) == s([3, 1]).to("m")
True
>>> jack.P([2]).at(alpha=2) == jack.zonal([2])
True
```

## Hall–Littlewood

`Qp` is `Q'_λ`, the family whose Schur coefficients are the Kostka–Foulkes
polynomials `K_{μλ}(t)` in the charge convention. Two specializations pin it,
and both fail under `t → 1/t`:

```pycon
>>> from symfn import s, h, hl
>>> hl.Qp([2, 1]).at(t=0) == s([2, 1])
True
>>> hl.Qp([2, 1]).at(t=1) == h([2, 1]).to("s")
True
```

`K_{λμ}(1)` is the Kostka number, which is the cross-check against the contract
layer's own count:

```pycon
>>> import symfn
>>> hl.kostka_foulkes([2, 2], [1, 1, 1, 1]).at(1)
2
>>> symfn.kostka_number([2, 2], [1, 1, 1, 1])
2
```

The absence of a constant term in `K_{(2,2),(1^4)}(t) = t^2 + t^4` is what
distinguishes charge from its cocharge rival.

The same matrix runs the other way. `to_P` rewrites a Schur-basis element in
the `P` basis, `to_Qp` in the `Q'` basis, and the two smallest values tell the
directions apart — the `t` lands on the smaller shape with a plus sign in one
and on the larger shape with a minus in the other:

```pycon
>>> hl.to_P(s([2]))
t*HLP[1,1] + HLP[2]
>>> hl.to_Qp(s([1, 1]))
HLQp[1,1] - t*HLQp[2]
>>> hl.to_P(hl.P([2, 1]))
HLP[2,1]
```

The tags are the names Sage prints — `HLP(s[2])` and `HLQp(s[1,1])` give the
same two values — and an element in one of them cannot be specialized: `at`
raises, because a `Sym` carries only the six classical bases.

## LLT

One family, two presentations, and the conventions in circulation differ by
more than a twist, so each entry point states which function it computes. `G`
is the tuple model: `G_ν(x; q) = Σ_T q^{inv(T)} x^T` over semistandard
fillings of a tuple of skew shapes, `inv` counting attacking pairs that are
out of order. `Gtilde`, `Htilde` and `H` are the ribbon model at level `k`,
summing over `k`-ribbon tableaux: `H` grades by spin, `Gtilde` and `Htilde`
by cospin, and `Htilde(mu, k)` is `Gtilde(k·mu, k)`.

Every entry point returns the **monomial** basis, not Schur, and carries one
parameter rather than two:

```pycon
>>> from symfn import llt
>>> llt.G([[1], [1]]).basis, llt.G([[1], [1]]).parameters
('m', ('q',))
>>> llt.G([[1], [1]]).at(q=2)
3*m[1,1] + m[2]
```

The spin/cospin split is the first value to check, because the two gradings
reverse — `H = q^{s*}·H̃(x; 1/q)`, and the `q` below sits where `Htilde`
carries the constant:

```pycon
>>> llt.H([1, 1], 2)
(1 + q)*m[1,1] + q*m[2]
>>> llt.Htilde([1, 1], 2)
(1 + q)*m[1,1] + m[2]
```

The raw inversion grading is **not** normalized — the floor `min_T inv(T)` is
left in, because it is real data about the shape tuple and hiding it is how
the quotient dictionary gets misread. The floor can be forced: `((1), (11))`
is the 2-quotient of `(2, 2, 2)`, and no offset choice brings its floor to
zero:

```pycon
>>> llt.min_inv([[1], [1, 1]])
1
```

Sage's dictionary, with the grading variable named `t` there and `q` here:
`Sym.llt(k).hspin()` is `llt.H`, `hcospin()` is `llt.Htilde`, and `cospin()`
on a partition is `llt.Gtilde`. On a tuple, Sage's `cospin()` divides the
floor out; divide `G` by `q^{llt.min_inv(...)}` to compare.

## Schubert

Schubert polynomials are indexed by permutations in one-line notation, in the
convention where `𝔖_{21} · 𝔖_{21} = 𝔖_{312}`:

```pycon
>>> from symfn import X
>>> X[2, 1] * X[2, 1]
X[3,1,2]
```

The rival convention indexes by inverse permutations and gives `𝔖_{231}`.
