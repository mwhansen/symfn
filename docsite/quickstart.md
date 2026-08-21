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
>>> (s([2]) * s([1])).coefficient([2, 1])
1
```

## Bases do not mix silently

An element carries the basis it is written in, and combining two that disagree
raises rather than picking one for you:

```pycon
>>> s([2, 1]) + h([2])
Traceback (most recent call last):
  ...
symfn.BasisError: cannot combine s with h; convert one with .to()
```

Convert explicitly, and the conversion is exact:

```pycon
>>> h([2]).to("s")
s[2]
>>> e([1, 1]).to("m")
2*m[1,1] + m[2]
>>> s([1, 1]).to("p")
1/2*p[1,1] - 1/2*p[2]
```

Conversions into the power-sum basis are rational, so they come back as
`Fraction`. Everything else lands in the integers.

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

Macdonald, Jack, Hall–Littlewood and LLT come back with their parameters
intact, as objects that print readably:

```pycon
>>> from symfn import macdonald, jack, hl
>>> jack.P([2])
2/(alpha + 1)*m[1,1] + m[2]
>>> macdonald.Q([1]).coefficient([1])
(1 - t)/(1 - q)
>>> hl.Qp([1, 1])
s[1,1] + t*s[2]
```

`at` substitutes the parameters and hands back an ordinary `Sym`, which is how
a family rejoins the arithmetic above:

```pycon
>>> jack.P([2]).at(alpha=1)
m[1,1] + m[2]
>>> hl.Qp([2, 1]).at(t=1)
s[2,1] + s[3]
>>> macdonald.qt_kostka([2], [1, 1])
t
```

Each family's element states which basis it is written in, and it is not always
the same one — Macdonald and Jack come back in the monomial basis, Hall–Littlewood
in Schur:

```pycon
>>> macdonald.P([2]).basis, hl.Qp([2]).basis
('m', 's')
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
