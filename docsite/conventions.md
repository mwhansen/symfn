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

## LLT

The LLT conventions differ by more than a twist, so each entry point names
which `G` or `H` it computes, and the raw inversion grading is **not**
normalized — the floor `min_T inv(T)` is left in, because it is real data about
the shape tuple and hiding it is how the quotient dictionary gets misread.

The LLT family returns in the **monomial** basis, not Schur:

```pycon
>>> from symfn import llt
>>> llt.G([[1], [1]]).basis
'm'
```

Sage's `llt(k).cospin(tuple)` divides the floor out; divide by
`q^{symfn.llt_min_inv(...)}` to compare.

## Schubert

Schubert polynomials are indexed by permutations in one-line notation, in the
convention where `𝔖_{21} · 𝔖_{21} = 𝔖_{312}`:

```pycon
>>> from symfn import X
>>> X[2, 1] * X[2, 1]
X[3,1,2]
```

The rival convention indexes by inverse permutations and gives `𝔖_{231}`.
