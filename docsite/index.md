# symfn

Exact symmetric functions for Python: the six classical bases, Hall–Littlewood,
Macdonald, LLT, Jack, and Schubert polynomials, with no dependencies and no
Sage required.

```pycon
>>> from symfn import s, h, macdonald, jack
>>> s([2, 1]) * s([1])
s[2,1,1] + s[2,2] + s[3,1]
>>> h([2, 1]).to("s")
s[2,1] + s[3]
>>> macdonald.P([2]).at(q=5, t=5) == s([2]).to("m")
True
>>> jack.P([2])
2/(alpha + 1)*m[1,1] + m[2]
```

Every value is exact or the call fails: coefficients are Python `int` of any
size and `Fraction` where a denominator exists, nothing is rounded, and a
computation that cannot be exact raises rather than approximating.

## Two layers, both supported

The **convenience layer** is what most callers want. `Sym` is a symmetric
function that knows which basis it is written in, so a basis mix-up raises
instead of quietly returning a plausible wrong answer, and the parameter
families come back as objects you can read and specialize.

The **contract layer** is the compiled module underneath: whole-object entry
points over plain lists of `(partition, coefficient)` pairs. Reach for it when
you are marshalling in bulk, or building another library on top — it is the
surface Sage itself is pinned to, and it changes slowly and deliberately.

The convenience layer is defined entirely in terms of the contract layer and
computes nothing of its own, so the two cannot disagree.

```{toctree}
:maxdepth: 2
:caption: Guide

install
quickstart
conventions
```

```{toctree}
:maxdepth: 2
:caption: Reference

api/convenience
api/families
api/contract
```
