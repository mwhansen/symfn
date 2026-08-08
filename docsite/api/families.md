# The parameter families

Macdonald, Jack, Hall–Littlewood and LLT, each as a namespace over the contract
layer, returning coefficients you can read and specialize.

Every family here has a normalization that plausible rivals disagree with, so
each method names its convention and carries an example whose value rules the
rivals out. [Conventions](../conventions.md) collects the specializations that
pin them.

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


.. autoclass:: symfn.AlphaFrac

```
