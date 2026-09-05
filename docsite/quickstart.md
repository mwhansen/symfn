# Quickstart

Everything below runs in a bare interpreter with nothing but symfn installed.

## Building elements

The six classical bases are module-level factories. Each takes a partition, and
each is also subscriptable, which reads like the output:

```pycon
>>> from symfn import s, h, e, p, m, f
>>> s([2, 1])
s[2,1]
>>> s[2, 1] == s([2, 1])
True
>>> s()
1
```

Arithmetic works as you would expect, with `int` and `Fraction` scalars:

```pycon
>>> 2 * s([1]) + s([2])
2*s[1] + s[2]
>>> s([1]) ** 3
s[1,1,1] + 2*s[2,1] + s[3]
>>> s([2]) / 2
1/2*s[2]
>>> (s([2]) * s([1])).coefficient([2, 1])
1
```

## Bases do not mix silently

An element knows which basis it is written in, and combining two that disagree
raises rather than picking a basis for the output:

```pycon
>>> s([2, 1]) + h([2])
Traceback (most recent call last):
  ...
symfn.BasisError: cannot combine s with h; convert one with .to()
```

The `to` method provides explicit conversion:

```pycon
>>> h([2]).to("s")
s[2]
>>> e([1, 1]).to("m")
2*m[1,1] + m[2]
>>> s([1, 1]).to("p")
1/2*p[1,1] - 1/2*p[2]
```

Conversions into the power-sum basis are rational, so the coefficients are
`Fraction` objects. Everything else has `int` coefficients.

## Operations on the ring

```pycon
>>> s([2, 1, 1]).omega()
s[3,1]
>>> s([2, 1]).scalar(s([2, 1]))
1
>>> s([2, 1]).skew_by(s([1]))
s[1,1] + s[2]
>>> s([2]).plethysm(s([2]))
s[2,2] + s[4]
>>> s([2, 1]).internal_product(s([2, 1]))
s[1,1,1] + s[2,1] + s[3]
>>> s([2]).coproduct()
{((), (2,)): 1, ((1,), (1,)): 1, ((2,), ()): 1}
```

## Leaving the ring

```pycon
>>> s([2]).expand(2)
{(1, 1): 1, (2, 0): 1, (0, 2): 1}
>>> s([2, 1]).evaluate([1, 1, 1])
8
>>> s([2, 1]).dimension()
2
>>> (s([2, 1]) + s([3])).principal_specialization(3)
18
```

## Parameter families

Macdonald, Jack and Hall–Littlewood each have their own bases:

```pycon
>>> from symfn import macdonald, jack, hl
>>> jack.P([2])
JackP[2]
>>> macdonald.Q([1])
McdQ[1]
>>> hl.Qp([1, 1])
HLQp[1,1]
```

`to` expands one into the classical basis its family is defined in — monomial
for Macdonald and Jack, Schur for Hall–Littlewood — and that is where the
parameters become visible:

```pycon
>>> jack.P([2]).to("m")
2/(alpha + 1)*m[1,1] + m[2]
>>> macdonald.Q([1]).to("m").coefficient([1])
(1 - t)/(1 - q)
>>> hl.Qp([1, 1]).to("s")
s[1,1] + t*s[2]
```

`basis` says which family an element belongs to, and any other classical basis
is one more change of basis on top of the expansion — a ℤ-linear map on the
partitions, so the parameters ride along:

```pycon
>>> macdonald.P([2]).basis, hl.Qp([2]).basis
('McdP', 'HLQp')
>>> hl.Qp([1, 1]).to("m")
(1 + t)*m[1,1] + t*m[2]
>>> hl.Qp([1, 1]).to("h")
h[1,1] + (-1 + t)*h[2]
```

That holds for the rational-function families too, where the coefficients are
in ℚ(q,t) or ℚ(α):

```pycon
>>> macdonald.P([2]).to("s")
(-t + q)/(1 - q*t)*s[1,1] + s[2]
>>> jack.P([2]).to("s")
(1 - alpha)/(alpha + 1)*s[1,1] + s[2]
```

The second vanishes at α = 1, where Jack `P` is the Schur function, and the
first at `q = t`, where Macdonald `P` is. The power-sum basis is reached the
same way — that conversion divides by z_μ, so it is the one whose
coefficients pick up rational parts:

```pycon
>>> jack.P([2]).to("p")
1/(alpha + 1)*p[1,1] + alpha/(alpha + 1)*p[2]
```

Each family is orthogonal under its own deformed pairing rather than the
classical one — `scalar_t`, `scalar_qt` and `scalar_jack`, named as Sage
spells them:

```pycon
>>> jack.P([2]).scalar_jack(jack.P([1, 1]))
0
>>> jack.P([2]).scalar_jack(jack.Q([2]))
1
>>> macdonald.P([2]).scalar_qt(macdonald.Q([2]))
1
>>> s([1]).scalar_t(s([1]))
1/(1 - t)
```

The expansions run the other way too. `to_P`, `to_Q`, `to_J` and `to_Htilde`
rewrite a classical element *into* a family's basis — the direction a
positivity question asks in — and each takes the basis its forward sibling
expands in:

```pycon
>>> from symfn import m, s
>>> jack.to_P(m([2]))
-2/(alpha + 1)*JackP[1,1] + JackP[2]
>>> hl.to_Qp(s([1, 1]))
HLQp[1,1] - t*HLQp[2]
>>> jack.to_P(jack.P([2, 1]).to("m"))
JackP[2,1]
```

The family tags are basis codes on the same footing as the classical six —
`basis` above returned one — so the same `to` reaches them. It takes any of
the fifteen codes, converting to an element's own basis is the identity, and
the `to_*` methods above are the named conventions behind the parametric
targets:

```pycon
>>> s([1, 1]).to("HLQp")
HLQp[1,1] - t*HLQp[2]
>>> s([1, 1]).to("HLQp") == hl.to_Qp(s([1, 1]))
True
>>> hl.P([2, 1]).to("HLP")
HLP[2,1]
```

To build an element rather than compute one, multiply a shape by a scalar.
`q`, `t` and `alpha` are values, and an element takes an `int`, a `Fraction` or
any polynomial in them:

```pycon
>>> from symfn import q, t, alpha
>>> q * macdonald.Htilde([2, 1]) + t * macdonald.Htilde([3])
q*McdHt[2,1] + t*McdHt[3]
>>> (1 - alpha) * jack.P([2])
(1 - alpha)*JackP[2]
>>> ((1 - alpha) * jack.P([2])).to("m")
(2 - 2*alpha)/(alpha + 1)*m[1,1] + (1 - alpha)*m[2]
```

A classical element scales the same way and stays in its own basis, now
carrying the parameter — which is how a scaled one reaches the expansions:

```pycon
>>> q * m([2])
q*m[2]
>>> macdonald.to_P(q * m([2])).coefficient([2])
q
```

A coefficient read off an element has the same arithmetic, so a value can be
carried around and put back without leaving exact form:

```pycon
>>> c = macdonald.Q([1]).to("m").coefficient([1])
>>> c + c
(2 - 2*t)/(1 - q)
>>> c * m([1])
(1 - t)/(1 - q)*m[1]
```

Two elements add or multiply when they are in the same basis — a parametric
basis multiplies through the basis its family expands in and comes back in
its own. ⚠️ Hall–Littlewood is in `t` alone and LLT in `q` alone, with their
own one-variable symbols `t_hl` and `q_llt`; the two-variable `q` and `t`
above also scale their elements when the value involves only the element's
variable, and a value genuinely in both refuses.

`at` substitutes the parameters and hands back an ordinary `Sym`, expanding
first if it has to, which is how a family rejoins the arithmetic above:

```pycon
>>> jack.P([2]).at(alpha=1)
m[1,1] + m[2]
>>> hl.Qp([2, 1]).at(t=1)
s[2,1] + s[3]
>>> macdonald.qt_kostka([2], [1, 1])
t
```

A subset of the parameters, by name, is substituted with the rest kept —
`q = 0` is the Hall–Littlewood degeneration of the Macdonald family:

```pycon
>>> macdonald.P([2]).at(q=0)
(1 - t)*m[1,1] + m[2]
```

LLT is the exception: it has no basis of its own here, because the kernel has
no expansion *into* one, so `llt.G` and its siblings come back in the monomial
basis directly. A shape in `G`'s tuple may be skew, written as an
`(outer, inner)` pair beside the plain partitions:

```pycon
>>> from symfn import llt
>>> llt.G([([2, 1], [1]), [1]])
(3 + 3*q)*m[1,1,1] + (2 + q)*m[2,1] + m[3]
```

The Delta operators and the Macdonald eigenoperators are on the same namespace:

```pycon
>>> from symfn import s, macdonald
>>> macdonald.nabla_e(2)
(t + q)*s[1,1] + s[2]
>>> macdonald.delta_prime_e(1, 2)
(t + q)*s[1,1] + s[2]
>>> macdonald.theta_ek(1, s([1, 1])).support()
[(1, 1, 1), (2, 1)]
```

## Schubert polynomials

Schubert polynomials are indexed by permutations in one-line notation, so they
live in their own type and never mix with symmetric functions:

```pycon
>>> from symfn import X
>>> X[2, 1] * X[2, 1]
X[3,1,2]
>>> X[1, 3, 2].divided_difference(2)
1
>>> X[1, 3, 2].expand()
{(0, 1): 1, (1,): 1}
>>> X[1, 3, 2].dimension()
2
>>> X[1, 3, 2].pairing(X[3, 1, 2], 3)
1
```

`expand` and `from_polynomial` are inverse, and `stanley_schur` crosses back to
the symmetric functions:

```pycon
>>> from symfn import from_polynomial, stanley_schur
>>> from_polynomial(X[1, 3, 2].expand()) == X[1, 3, 2]
True
>>> stanley_schur([3, 2, 1])
s[2,1]
```

## Dropping to the contract layer

When you are marshalling in bulk, the compiled entry points take and return
plain lists and skip the object construction entirely:

```pycon
>>> import symfn
>>> symfn.schur_multiply([([2, 1], 1)], [([1], 1)])
[((2, 1, 1), 1), ((2, 2), 1), ((3, 1), 1)]
>>> symfn.kostka_number([3, 1], [2, 1, 1])
2
```

Both layers are supported and both are documented; the convenience layer is
built out of exactly these calls.
