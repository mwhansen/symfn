# The convenience layer

The pure-Python surface: a basis-tagged element, the six basis factories, and
the Schubert type. Every method here is a composition of contract-layer calls
and computes nothing of its own, which `scripts/check_convenience.py` checks by
running each one against the contract sequence it claims to be.

## Elements

```{eval-rst}
.. autoclass:: symfn.Sym
```

## Basis factories

The six factories are instances rather than classes: `s([2, 1])` and `s[2, 1]`
both build the Schur function indexed by that partition, and `s()` is the unit.

```{eval-rst}
.. autoclass:: symfn._sym._Factory

.. autodata:: symfn.s
   :annotation:

.. autodata:: symfn.h
   :annotation:

.. autodata:: symfn.e
   :annotation:

.. autodata:: symfn.p
   :annotation:

.. autodata:: symfn.m
   :annotation:

.. autodata:: symfn.f
   :annotation:
```

## Skew Schur functions

```{eval-rst}
.. autofunction:: symfn.skew
```

## Stanley symmetric functions

```{eval-rst}
.. autofunction:: symfn.stanley_schur
```

## Schubert polynomials

```{eval-rst}
.. autoclass:: symfn.Schub

.. autoclass:: symfn._schubert._SchubFactory

.. autodata:: symfn.X
   :annotation:

.. autofunction:: symfn.from_polynomial
```

## Basis identity, and the base ring

An element carries two facts that are independent of each other: the basis it
is written in, and the ring its coefficients live in. A mismatch in either is
refused, and the two are different exceptions because only one of them is
fixed by converting.

```{eval-rst}
.. autoexception:: symfn.BasisError
.. autoexception:: symfn.BaseRingError
```
