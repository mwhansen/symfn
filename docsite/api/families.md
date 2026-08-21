# The parameter families

Macdonald, Jack, Hall–Littlewood and LLT, each as a namespace over the contract
layer, returning coefficients you can read and specialize.

Each namespace runs in **both directions**, and both speak of the same
objects. The capitalized methods — `P`, `Q`, `J`, `Htilde`, `Qp` — take a
shape and return it in the family's own basis, tagged as Sage prints it; the
`to_*` methods take a classical element and rewrite it *into* that same basis,
which is the direction a positivity question asks in. `Param.to` expands
either one into the classical basis the family is defined in — monomial for
Macdonald and Jack, Schur for Hall–Littlewood and `H̃` — and each `to_*` takes
that basis, refusing another rather than converting silently.

LLT is the exception on both counts: `G` and its siblings return the monomial
basis directly and there is no `to_*`, because LLT polynomials are not a basis
of Λ and there is nothing to rewrite into.

Every family here has a normalization that plausible rivals disagree with, so
each method names its convention and carries an example whose value rules the
rivals out. [Conventions](../conventions.md) collects the specializations that
pin them — and ⚠️ note that for the inverse direction those specializations are
not enough on their own: `α = 1` and `q = t` fix the very twists they would
have to separate, so the `to_*` methods pin themselves at free parameters.

## Macdonald

Reached as `symfn.macdonald`.

```{eval-rst}
.. autoclass:: symfn._families._Macdonald

.. autodata:: symfn.macdonald
   :annotation:
```

## Jack

Reached as `symfn.jack`.

```{eval-rst}
.. autoclass:: symfn._families._Jack

.. autodata:: symfn.jack
   :annotation:
```

## Hall–Littlewood

Reached as `symfn.hl`. The name is `hl` rather than `hall_littlewood` because
`symfn.hall_littlewood` is a contract entry point, and the flat surface keeps
its names.

```{eval-rst}
.. autoclass:: symfn._families._HallLittlewood

.. autodata:: symfn.hl
   :annotation:
```

## LLT

Reached as `symfn.llt`.

```{eval-rst}
.. autoclass:: symfn._families._LLT

.. autodata:: symfn.llt
   :annotation:
```

## Coefficient types

What a parameter family's coefficients come back as. Each holds the contract
layer's rows unchanged and adds a `repr` and an evaluation map.

```{eval-rst}
.. autoclass:: symfn.Param


.. autoclass:: symfn.Poly


.. autoclass:: symfn.QtPoly


.. autoclass:: symfn.QtFrac


.. autoclass:: symfn.QtRatio


.. autoclass:: symfn.AlphaFrac

```
