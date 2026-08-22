# symfn

An efficient library built on a Rust core for working with symmetric
functions and other related objects: the six classical bases,
Hall–Littlewood, Macdonald, LLT, Jack, and Schubert polynomials.

```pycon
>>> from symfn import s, h, macdonald, jack
>>> s([2, 1]) * s([1])
s[2,1,1] + s[2,2] + s[3,1]
>>> h([2, 1]).to("s")
s[2,1] + s[3]
>>> macdonald.P([2]).at(q=5, t=5) == s([2]).to("m")
True
>>> jack.P([2]).to("m")
2/(alpha + 1)*m[1,1] + m[2]
```

Every value is exact or the call fails: coefficients are Python `int` of any
size and `Fraction` where a denominator exists, nothing is rounded, and a
computation that cannot be exact raises rather than approximating.

## Two layers, both supported

The **convenience layer** is what most users want and provides a more
usable high-level interface at the expense of some overhead. The
primary class is `Sym` representing a symmetric function that knows
which basis it is written in.

The **contract layer** is the low-level compiled module underneath:
whole-object entry points which use plain `(partition, coefficient)`
pairs. Use this if you are marshalling objects in bulk, or building
another library on top of this. This interface should be the most
stable and change slowly and deliberately.

The computation uses the same compiled kernel either way — the
convenience layer only shapes arguments and wraps results — so the two
layers agree on every value, and the contract layer just skips the
object construction.

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
