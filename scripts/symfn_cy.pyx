# cython: language_level=3, boundscheck=False, wraparound=False
"""The per-term loop of the Sage backend shim, compiled.

Everything else in `sage_backend.py` runs once per *call*; this runs once per
output *term*, and a conversion of degree 18 produces tens of thousands of them.
In pure Python that loop measured **185 ns/term** -- unpack a 3-tuple, index a
list, build an `Integer`, insert into a dict -- which came to 0.76x the entire
Rust computation it was wrapping. None of those four steps is avoidable in
Python, so the loop itself has to stop being Python.

The partitions arrive as **indices** into `symfn.partitions(degree)` rather than
as lists of parts (see `convert_indexed`), which is what makes the inner step a
C array read instead of building a tuple and hashing it. Sage's own Symmetrica
wrapper cannot do this: `symmetrica.pxi` gets lists of parts back from C and
calls `Partition(res)` on each, so it pays object construction per term where
this pays an index.
"""

from cpython.dict cimport PyDict_SetItem
from cpython.list cimport PyList_GET_ITEM, PyList_GET_SIZE
from cpython.long cimport PyLong_AsLongAndOverflow
from cpython.tuple cimport PyTuple_GET_ITEM

from sage.rings.integer cimport Integer, smallInteger


cpdef dict build_terms(list raw, dict parts_by_degree):
    """`[(degree, index, coefficient)]` -> `{Partition: Integer}`.

    Zero coefficients are dropped, matching what Sage's `_from_dict` expects.
    """
    cdef dict out = {}
    cdef Py_ssize_t k, size = PyList_GET_SIZE(raw)
    cdef object item, degree, index, coeff, table
    cdef Integer value
    cdef long small
    cdef int overflow

    # Nearly every call is homogeneous, so the degree lookup is hoisted and
    # only repeated when the degree actually changes.
    cdef object last_degree = None
    cdef list parts = None

    for k in range(size):
        item = <object>PyList_GET_ITEM(raw, k)
        coeff = <object>PyTuple_GET_ITEM(<tuple>item, 2)
        if not coeff:
            continue
        degree = <object>PyTuple_GET_ITEM(<tuple>item, 0)
        if degree is not last_degree:
            table = parts_by_degree[degree]
            parts = <list>table
            last_degree = degree
        index = <object>PyTuple_GET_ITEM(<tuple>item, 1)
        # `Integer(coeff)` parses generically; `smallInteger` is the constructor
        # Sage uses internally for values that fit a C long, which every
        # structure constant here does apart from the rare escalated one.
        small = PyLong_AsLongAndOverflow(coeff, &overflow)
        if overflow == 0:
            value = smallInteger(small)
        else:
            value = Integer(coeff)
        PyDict_SetItem(out, <object>PyList_GET_ITEM(parts, <Py_ssize_t>index), value)
    return out
