# The contract layer

The compiled module. Every name here is also available flat at `symfn.*` —
`symfn.schur_multiply` and `symfn.symfn.schur_multiply` are the same function.

This is the surface that changes slowly: the Sage adapter is pinned to it, so
it does not break between releases; the convenience layer above it is free to
grow.

```{eval-rst}
.. automodule:: symfn.symfn
   :members:
   :undoc-members:
   :member-order: bysource
```
