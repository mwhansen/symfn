# The contract layer

The compiled module: whole-object entry points over plain data. Every name here
is also available flat at `symfn.*` — `symfn.schur_multiply` and
`symfn.symfn.schur_multiply` are the same function.

This is the surface that changes slowly. It is what the Sage adapter is pinned
to, so a break here is a break for every consumer; the convenience layer above
is free to grow, and this is not.

The type aliases at the top of the page — `Partition`, `Element`, `QtElement`,
`JackCell`, … — name the plain-data shapes every signature is written in; each
one says what its tuples mean.

```{eval-rst}
.. automodule:: symfn.symfn
   :members:
   :undoc-members:
   :member-order: bysource
```
