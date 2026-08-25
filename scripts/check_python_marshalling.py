"""Round-trip the marshalling layer: values cross intact, in the promised shapes.

    cargo build --features python
    python3 scripts/check_python_marshalling.py target/debug/libsymfn.dylib

The proposition is `docs/policies/python.md` P1 read as a test: **a value
handed in comes back out intact, and every value handed out has exactly the
promised shape** — supports as tuples of ints, coefficients as Python ints of
any size, every denominator positive and in lowest terms, the factored
encodings factored. This is a boundary test, not an oracle: catching a kernel
defect is `docs/policies/validation.md`'s job, discharged by the Rust suites
and the committed fixtures. What no Rust test can see is the crossing itself —
a PyO3 conversion that truncates, a `list` where a `tuple` is promised, a
denominator dumped negative — and that is all this script asserts.

Needs no Sage — it imports the built extension module and nothing else, like
its siblings `check_python_boundary.py` (the failure half of the boundary) and
`check_convenience.py` (the layer above it).

Three checks:

  * **shapes** — every exported callable runs once on a small valid input and
    its return is validated against the stub alias `symfn.pyi` declares for
    it, `type() is` strict: a `tuple` where a tuple is promised, `int`
    coefficients (`bool` is not an `int` here), partitions weakly decreasing
    with no zeros, permutations with trailing fixed points dropped. The table
    must name every export, so a new entry point cannot ship unclassified.
  * **widths** — coefficients at both edges of the `i128` fast path and past
    it survive identity-shaped calls unchanged, which is the module doc's "no
    ceiling" claim exercised on both sides of the escalation. The `(q,t)`
    inbound wall at `i128` raises a typed `ValueError` rather than wrapping.
  * **permissive inbound** — the `*Arg` halves: a support crosses in as list,
    tuple, or padded with trailing zeros, an element as any sequence of
    pairs, and an output handed straight back in is accepted — all with
    identical answers.
"""

import importlib.machinery
import importlib.util
import math
import sys


def load(path):
    """Import the built cdylib under the name its `#[pymodule]` declares.

    Cargo names it `libsymfn.dylib`/`.so`, which is neither the name nor the
    suffix an import would accept, so the loader is named explicitly rather
    than inferred — this takes the artifact `cargo build` leaves behind and
    saves a maturin step.
    """
    loader = importlib.machinery.ExtensionFileLoader("symfn", path)
    spec = importlib.util.spec_from_loader("symfn", loader)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["symfn"] = mod
    loader.exec_module(mod)
    return mod


# --- validators, one per stub alias -----------------------------------------
#
# Each takes the returned value and answers whether it has the promised shape.
# `type() is` rather than `isinstance`: the promise is `tuple` and `int`
# themselves — a subclass (`bool` under `int` is the live case) is not what a
# consumer keying a dict on the value was told to expect.


def is_int(x):
    return type(x) is int


def is_partition(x):
    """A tuple of positive ints, weakly decreasing — the outbound normal form.

    Trailing zeros are tolerated inbound and never emitted outbound; a zero
    part coming back means the normalization was skipped on some path.
    """
    return (
        type(x) is tuple
        and all(type(p) is int and p > 0 for p in x)
        and all(a >= b for a, b in zip(x, x[1:]))
    )


def is_permutation(x):
    """One-line notation, 1-based, trailing fixed points dropped."""
    if type(x) is not tuple or not all(type(p) is int for p in x):
        return False
    if sorted(x) != list(range(1, len(x) + 1)):
        return False
    return not x or x[-1] != len(x)


def is_exponents(x):
    """A tuple of non-negative ints — a monomial's exponent vector, or an
    area sequence: keys that are positional rather than partitions."""
    return type(x) is tuple and all(type(p) is int and p >= 0 for p in x)


def is_pairs(key_ok, val_ok):
    def check(x):
        return type(x) is list and all(
            type(t) is tuple and len(t) == 2 and key_ok(t[0]) and val_ok(t[1])
            for t in x
        )

    return check


def is_rational(x):
    """`(numerator, denominator)`, denominator positive and in lowest terms."""
    return (
        type(x) is tuple
        and len(x) == 2
        and is_int(x[0])
        and is_int(x[1])
        and x[1] > 0
        and math.gcd(x[0], x[1]) == 1
    )


def is_t_coeff(x):
    """`(t_exponent, coefficient)` pairs, exponents non-negative."""
    return type(x) is list and all(
        type(t) is tuple
        and len(t) == 2
        and is_int(t[0])
        and t[0] >= 0
        and is_int(t[1])
        for t in x
    )


def is_qt_coeff(x):
    """`(q_exponent, t_exponent, coefficient)` triples."""
    return type(x) is list and all(
        type(t) is tuple
        and len(t) == 3
        and is_int(t[0])
        and t[0] >= 0
        and is_int(t[1])
        and t[1] >= 0
        and is_int(t[2])
        for t in x
    )


def is_qt_factors(x):
    """Denominator factors `(q_exp, t_exp, multiplicity)` for `(1 - q^a t^b)^m`:
    multiplicity at least 1, and no `(0, 0)` factor, which would be `0^m`."""
    return type(x) is list and all(
        type(t) is tuple
        and len(t) == 3
        and is_int(t[0])
        and is_int(t[1])
        and is_int(t[2])
        and t[2] >= 1
        and (t[0], t[1]) != (0, 0)
        for t in x
    )


def is_alpha_atoms(x):
    """Factors `(u, v, m)` for `(u·α + v)^m`: primitive, in increasing order."""
    if type(x) is not list:
        return False
    for t in x:
        if not (type(t) is tuple and len(t) == 3 and all(is_int(v) for v in t)):
            return False
        if t[2] < 1 or math.gcd(t[0], t[1]) != 1:
            return False
    return all(a[:2] < b[:2] for a, b in zip(x, x[1:]))


def is_int_list(x):
    return type(x) is list and all(is_int(v) for v in x)


def is_jack_cell(x):
    """`(dense numerator, denominator atoms, scale, tail)`, the value
    `numerator / (scale · Π atoms · tail)`; a zero scale would be division by
    zero. The tail is the general denominator factor a plethysm can produce and
    nothing else does, so it is empty on almost every row."""
    return (
        type(x) is tuple
        and len(x) == 4
        and is_int_list(x[0])
        and is_alpha_atoms(x[1])
        and is_int(x[2])
        and x[2] != 0
        and is_int_list(x[3])
    )


element = is_pairs(is_partition, is_int)
rational_element = is_pairs(is_partition, is_rational)
monomials = is_pairs(is_exponents, is_int)
t_element = is_pairs(is_partition, is_t_coeff)
qt_element = is_pairs(is_partition, is_qt_coeff)
schubert_element = is_pairs(is_permutation, is_int)
coproduct_terms = is_pairs(
    lambda k: type(k) is tuple
    and len(k) == 2
    and is_partition(k[0])
    and is_partition(k[1]),
    is_int,
)
path_refinement = is_pairs(is_exponents, qt_element)


def is_pair_key(k):
    """The support of a coproduct: an ordered pair of partitions."""
    return (
        type(k) is tuple
        and len(k) == 2
        and is_partition(k[0])
        and is_partition(k[1])
    )


def pair_rows(*cell):
    """Coproduct rows whose coefficient is spread over `cell` trailing slots,
    one predicate each — how each coefficient ring encodes one value.
    """

    def ok(x):
        return type(x) is list and all(
            type(t) is tuple
            and len(t) == 1 + len(cell)
            and is_pair_key(t[0])
            and all(p(v) for p, v in zip(cell, t[1:]))
            for t in x
        )

    return ok


def expo_rows(*cell):
    """Expansion rows whose coefficient is spread over `cell` trailing slots."""

    def ok(x):
        return type(x) is list and all(
            type(t) is tuple
            and len(t) == 1 + len(cell)
            and is_exponents(t[0])
            and all(p(v) for p, v in zip(cell, t[1:]))
            for t in x
        )

    return ok


coproduct_qt_terms = pair_rows(is_qt_coeff)
coproduct_mac_terms = pair_rows(is_qt_coeff, is_qt_factors)
coproduct_jack_terms = pair_rows(
    is_int_list, is_alpha_atoms, lambda v: is_int(v) and v != 0, is_int_list
)
coproduct_ht_terms = pair_rows(
    is_qt_coeff, lambda v: type(v) is list and all(is_atom(a) for a in v)
)


def indexed_element(x):
    return type(x) is list and all(
        type(t) is tuple
        and len(t) == 3
        and is_int(t[0])
        and t[0] >= 0
        and is_int(t[1])
        and t[1] >= 0
        and is_int(t[2])
        for t in x
    )


def macdonald_element(x):
    return type(x) is list and all(
        type(t) is tuple
        and len(t) == 3
        and is_partition(t[0])
        and is_qt_coeff(t[1])
        and is_qt_factors(t[2])
        for t in x
    )


def ht_element(x):
    """`(mu, numerator, denominator atoms)` rows: a polynomial in q and t over
    `(kind, a, b, multiplicity)` atoms.

    The atom list may be empty — that is a denominator of 1, which is what a
    polynomial coefficient has. Kind is 0 or 1 and nothing else, and kind 1
    needs both exponents positive.
    """
    return type(x) is list and all(
        type(t) is tuple
        and len(t) == 3
        and is_partition(t[0])
        and is_qt_coeff(t[1])
        and type(t[2]) is list
        and all(is_atom(a) for a in t[2])
        for t in x
    )


def is_atom(a):
    """One `(kind, a, b, multiplicity)` denominator atom."""
    return (
        type(a) is tuple
        and len(a) == 4
        and all(type(v) is int and v >= 0 for v in a)
        and a[0] in (0, 1)
        and (a[1] or a[2])
        and (a[0] == 0 or (a[1] and a[2]))
        and a[3] > 0
    )


def jack_element(x):
    return type(x) is list and all(
        type(t) is tuple and len(t) == 5 and is_partition(t[0]) and is_jack_cell(t[1:])
        for t in x
    )


def mac_cell(x):
    """`(numerator, factored denominator)` — one Macdonald-family coefficient,
    which is a `macdonald_element` row without its partition.
    """
    return (
        type(x) is tuple
        and len(x) == 2
        and is_qt_coeff(x[0])
        and is_qt_factors(x[1])
    )


def ht_cell(x):
    """`(numerator, denominator atoms)` — one `H̃` coefficient."""
    return (
        type(x) is tuple
        and len(x) == 2
        and is_qt_coeff(x[0])
        and type(x[1]) is list
        and all(is_atom(a) for a in x[1])
    )


def is_rat_qt_coeff(x):
    """`(q_exponent, t_exponent, numerator, denominator)` rows — a polynomial
    in q and t whose coefficients divide, each denominator positive."""
    return type(x) is list and all(
        type(t) is tuple
        and len(t) == 4
        and is_int(t[0])
        and t[0] >= 0
        and is_int(t[1])
        and t[1] >= 0
        and is_int(t[2])
        and is_int(t[3])
        and t[3] > 0
        for t in x
    )


rat_qt_element = is_pairs(is_partition, is_rat_qt_coeff)


def rat_macdonald_element(x):
    """`macdonald_element` rows with the numerator coefficients split — what
    the conversions that divide return."""
    return type(x) is list and all(
        type(t) is tuple
        and len(t) == 3
        and is_partition(t[0])
        and is_rat_qt_coeff(t[1])
        and is_qt_factors(t[2])
        for t in x
    )


def rat_ht_element(x):
    """`ht_element` rows with the numerator coefficients split."""
    return type(x) is list and all(
        type(t) is tuple
        and len(t) == 3
        and is_partition(t[0])
        and is_rat_qt_coeff(t[1])
        and type(t[2]) is list
        and all(is_atom(a) for a in t[2])
        for t in x
    )


def zonal_terms(x):
    """`(mu, numerator, denominator)` rows, each fraction in lowest terms."""
    return type(x) is list and all(
        type(t) is tuple
        and len(t) == 3
        and is_partition(t[0])
        and is_rational(t[1:])
        for t in x
    )


def keyed_by_partition(val_ok):
    return is_pairs(is_partition, val_ok)


def list_of(item_ok):
    def check(x):
        return type(x) is list and all(item_ok(v) for v in x)

    return check


def optional(ok):
    def check(x):
        return x is None or ok(x)

    return check


def stanley_rows(x):
    """`(lambda, mu, nu, numerator, denominator atoms, scalar)` rows — the
    shapes as lists here, not tuples, as `symfn.pyi` declares for the tables."""
    return type(x) is list and all(
        type(t) is tuple
        and len(t) == 6
        and is_int_list(t[0])
        and is_int_list(t[1])
        and is_int_list(t[2])
        and is_int_list(t[3])
        and is_alpha_atoms(t[4])
        and is_int(t[5])
        and t[5] != 0
        for t in x
    )


def b_table(x):
    return type(x) is list and all(
        type(t) is tuple
        and len(t) == 5
        and is_int_list(t[0])
        and is_int_list(t[1])
        and is_int_list(t[2])
        and is_int_list(t[3])
        and is_int(t[4])
        and t[4] != 0
        for t in x
    )


def gj_tables(x):
    return type(x) is tuple and len(x) == 2 and b_table(x[0]) and b_table(x[1])


def core_and_quotient(x):
    return (
        type(x) is tuple
        and len(x) == 2
        and is_partition(x[0])
        and list_of(is_partition)(x[1])
    )


def tableaux(x):
    return list_of(list_of(is_int_list))(x)


def is_none(x):
    return x is None


# --- the shape table --------------------------------------------------------
#
# name -> (args, validator). Every exported callable must appear;
# `check_surface_is_covered` fails when one is added without a row, which is
# what keeps a new encoding from shipping unvalidated. Inputs are small — the
# whole sweep is under a second — and valid; the malformed half is
# `check_python_boundary.py`'s.

A = [([2, 1], 3), ([1, 1, 1], -1)]
B = [([1], 1)]
QT_A = [([2, 1], [(0, 0, 2), (1, 1, -1)])]
SCHUB = [([2, 1], 2), ([1, 3, 2], 1)]
SHAPES = {
    "antipode": ((A,), element),
    "big_pi": ((QT_A,), qt_element),
    "character_table": ((3,), list_of(is_int_list)),
    "character_value": (([2, 1], [1, 1, 1]), is_int),
    "chromatic_from_llt": ((3, [(0, 1)], []), qt_element),
    "class_algebra_coefficient": (([2, 1], [2, 1], [1, 1, 1]), is_int),
    "clear_caches": ((), is_none),
    "convert_indexed": ((A, "s", "m"), indexed_element),
    "convert_terms": ((A, "s", "m"), element),
    "convert_qt_terms": ((QT_A, "s", "m"), qt_element),
    "convert_macdonald_terms": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1)])], "m", "s"),
        macdonald_element,
    ),
    "convert_jack_terms": (
        ([([2, 1], [1, -2], [(1, 1, 1)], 3, [])], "m", "s"),
        jack_element,
    ),
    "convert_ht_terms": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])], "s", "m"),
        ht_element,
    ),
    "schur_multiply_qt": ((QT_A, QT_A), qt_element),
    "schur_multiply_macdonald": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1)])], [([], [(0, 0, 1)], [])]),
        macdonald_element,
    ),
    "schur_multiply_jack": (
        ([([2, 1], [1, -2], [(1, 1, 1)], 3, [])], [([], [1], [], 1, [])]),
        jack_element,
    ),
    "schur_multiply_ht": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])], [([], [(0, 0, 1)], [])]),
        ht_element,
    ),
    "omega_qt_terms": ((QT_A,), qt_element),
    "antipode_qt_terms": ((QT_A,), qt_element),
    "omega_macdonald_terms": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1)])],),
        macdonald_element,
    ),
    "antipode_macdonald_terms": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1)])],),
        macdonald_element,
    ),
    "omega_jack_terms": (([([2, 1], [1, -2], [(1, 1, 1)], 3, [])],), jack_element),
    "antipode_jack_terms": (
        ([([2, 1], [1, -2], [(1, 1, 1)], 3, [])],),
        jack_element,
    ),
    "omega_ht_terms": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])],),
        ht_element,
    ),
    "antipode_ht_terms": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])],),
        ht_element,
    ),
    "skew_by_qt": ((QT_A, [([], [(0, 0, 1)])], "s"), qt_element),
    "skew_by_macdonald": (
        (
            [([2, 1], [(0, 0, 1)], [(1, 1, 1)])],
            [([], [(0, 0, 1)], [])],
            "s",
        ),
        macdonald_element,
    ),
    "skew_by_jack": (
        ([([2, 1], [1, -2], [(1, 1, 1)], 3, [])], [([], [1], [], 1, [])], "s"),
        jack_element,
    ),
    "skew_by_ht": (
        (
            [([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])],
            [([], [(0, 0, 1)], [])],
            "s",
        ),
        ht_element,
    ),
    "coproduct": ((A,), coproduct_terms),
    "coproduct_qt": ((QT_A,), coproduct_qt_terms),
    "coproduct_macdonald": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1)])],),
        coproduct_mac_terms,
    ),
    "coproduct_jack": (
        ([([2, 1], [1, -2], [(1, 1, 1)], 3, [])],),
        coproduct_jack_terms,
    ),
    "coproduct_ht": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])],),
        coproduct_ht_terms,
    ),
    "expand_qt": ((QT_A, 3), expo_rows(is_qt_coeff)),
    "expand_macdonald": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1)])], 3),
        expo_rows(is_qt_coeff, is_qt_factors),
    ),
    "expand_jack": (
        ([([2, 1], [1, -2], [(1, 1, 1)], 3, [])], 3),
        expo_rows(
            is_int_list,
            is_alpha_atoms,
            lambda v: is_int(v) and v != 0,
            is_int_list,
        ),
    ),
    "expand_ht": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])], 3),
        expo_rows(
            is_qt_coeff,
            lambda v: type(v) is list and all(is_atom(a) for a in v),
        ),
    ),
    "evaluate_qt": ((QT_A, [1, 1, 1]), is_qt_coeff),
    "evaluate_macdonald": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1)])], [1, 1, 1]),
        mac_cell,
    ),
    "evaluate_jack": (
        ([([2, 1], [1, -2], [(1, 1, 1)], 3, [])], [1, 1, 1]),
        is_jack_cell,
    ),
    "evaluate_ht": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])], [1, 1, 1]),
        ht_cell,
    ),
    "dimension_qt": ((QT_A,), is_qt_coeff),
    "dimension_macdonald": (([([2, 1], [(0, 0, 1)], [(1, 1, 1)])],), mac_cell),
    "dimension_jack": (([([2, 1], [1, -2], [(1, 1, 1)], 3, [])],), is_jack_cell),
    "dimension_ht": (([([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])],), ht_cell),
    "principal_specialization_qt": ((QT_A, 3), is_qt_coeff),
    "principal_specialization_macdonald": (([([2, 1], [(0, 0, 1)], [(1, 1, 1)])], 3), mac_cell),
    "principal_specialization_jack": (([([2, 1], [1, -2], [(1, 1, 1)], 3, [])], 3), is_jack_cell),
    "principal_specialization_ht": (([([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])], 3), ht_cell),
    "principal_specialization_q_qt": (([([2, 1], [(0, 1, 1)])], 3), is_qt_coeff),
    "internal_product_qt": ((QT_A, QT_A), qt_element),
    "internal_product_macdonald": (
        (
            [([2, 1], [(0, 0, 1)], [(1, 1, 1)])],
            [([2, 1], [(0, 0, 1)], [])],
        ),
        macdonald_element,
    ),
    "internal_product_jack": (
        ([([2, 1], [1, -2], [(1, 1, 1)], 3, [])], [([2, 1], [1], [], 1, [])]),
        jack_element,
    ),
    "internal_product_ht": (
        (
            [([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])],
            [([2, 1], [(0, 0, 1)], [])],
        ),
        ht_element,
    ),
    "principal_specialization_at_qt": ((QT_A, 3, [(0, 1, 1)]), is_qt_coeff),
    "principal_specialization_at_macdonald": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1)])], 3, ([(1, 0, 1)], [])),
        mac_cell,
    ),
    "principal_specialization_at_jack": (
        ([([2, 1], [1, -2], [(1, 1, 1)], 3, [])], 3, ([0, 1], [], 1, [])),
        is_jack_cell,
    ),
    "principal_specialization_at_ht": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])], 3, ([(0, 1, 1)], [])),
        ht_cell,
    ),
    "plethysm_qt": ((QT_A, [([1], [(0, 1, 1)])]), qt_element),
    "plethysm_jack": (
        ([([2], [1], [], 1, [])], [([1], [1], [(1, 1, 1)], 1, [])]),
        jack_element,
    ),
    "plethysm_macdonald": (
        (
            [([2, 1], [(0, 0, 1)], [])],
            [([1], [(0, 0, 1)], [(1, 1, 1)])],
        ),
        macdonald_element,
    ),
    "plethysm_ht": (
        (
            [([2, 1], [(0, 0, 1)], [])],
            [([1], [(0, 1, 1)], [])],
        ),
        ht_element,
    ),
    "delta_conjecture_side": ((3, "rise"), list_of(qt_element)),
    "delta_ek": ((1, QT_A), qt_element),
    "delta_prime_e": ((1, 3), qt_element),
    "delta_prime_ek": ((1, QT_A), qt_element),
    "dimension": (([2, 1],), optional(is_int)),
    "elementary_to_schur": ((A,), element),
    "evaluate_schur": ((A, [1, 2]), is_int),
    "expand_alphabet": ((A, "s", 2), monomials),
    "forgotten_to_schur": ((A,), element),
    "gj_connection_tables": ((3,), gj_tables),
    "hall_inner_product": ((A, A), is_int),
    "hall_inner_product_qt": ((QT_A, QT_A), is_qt_coeff),
    "hall_inner_product_macdonald": (
        (
            [([2, 1], [(0, 0, 1)], [(1, 1, 1)])],
            [([2, 1], [(0, 0, 1)], [(1, 1, 1)])],
        ),
        mac_cell,
    ),
    "hall_inner_product_jack": (
        ([([2, 1], [1, -2], [(1, 1, 1)], 3, [])], [([2, 1], [1], [], 1, [])]),
        is_jack_cell,
    ),
    "hall_inner_product_ht": (
        (
            [([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])],
            [([2, 1], [(0, 0, 1)], [])],
        ),
        ht_cell,
    ),
    "hall_littlewood": (([2, 1],), t_element),
    "hall_littlewood_p": (([2, 1],), t_element),
    "hall_littlewood_p_table": ((3,), keyed_by_partition(t_element)),
    "schur_to_hall_littlewood_p": (([([2, 1], [(0, 1), (2, -3)])],), t_element),
    "schur_to_hall_littlewood_qp": (([([2, 1], [(0, 1), (2, -3)])],), t_element),
    "hall_littlewood_table": ((3,), keyed_by_partition(t_element)),
    "homogeneous_to_schur": ((A,), element),
    "ht_multiply": ((A, B), optional(element)),
    "ht_to_schur": ((A,), element),
    "htilde_by_llt": (([2, 1],), qt_element),
    "internal_product": ((A, A), element),
    "jack_j": (([2, 1],), jack_element),
    "jack_j_powersum": (([2, 1],), jack_element),
    "jack_norm_j": (([2, 1],), is_alpha_atoms),
    "jack_p": (([2, 1],), jack_element),
    "jack_q": (([2, 1],), jack_element),
    "jack_scalar": (
        ([([1], [0, 1], [], 1, [])], [([1], [1], [], 1, [])]),
        is_jack_cell,
    ),
    "scalar_t": (
        (
            [([2, 1], [(0, 0, 1)], [(1, 1, 1)])],
            [([2, 1], [(0, 0, 1)], [])],
        ),
        mac_cell,
    ),
    "scalar_qt": (
        (
            [([2, 1], [(0, 0, 1)], [(1, 1, 1)])],
            [([2, 1], [(0, 0, 1)], [])],
        ),
        mac_cell,
    ),
    "scalar_qt_ht": (
        (
            [([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])],
            [([2, 1], [(0, 0, 1)], [])],
        ),
        ht_cell,
    ),
    "to_power_qt": (([([2, 1], [(0, 1, 1)])], "s"), rat_qt_element),
    "to_power_macdonald": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1)])], "s"),
        rat_macdonald_element,
    ),
    "to_power_jack": (
        ([([2, 1], [1, -2], [(1, 1, 1)], 3, [])], "s"),
        jack_element,
    ),
    "to_power_ht": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])], "s"),
        rat_ht_element,
    ),
    "jack_structure_constant": (([2], [1], [2, 1]), is_jack_cell),
    "jack_table": ((3,), keyed_by_partition(jack_element)),
    "k_core_quotient": (([3, 1], 2), core_and_quotient),
    "kostka_foulkes": (([2, 1], [1, 1, 1]), is_t_coeff),
    "kostka_foulkes_column": (([1, 1, 1],), t_element),
    "kostka_foulkes_table": ((3,), list_of(list_of(is_t_coeff))),
    "kostka_number": (([2, 1], [1, 1, 1]), is_int),
    "kostka_table": ((3,), list_of(is_int_list)),
    "kronecker_coefficient": (([2, 1], [2, 1], [2, 1]), is_int),
    "llt_e_expansion": ((3, [(0, 1)], []), qt_element),
    "llt_fundamental": (([[1], [1]],), qt_element),
    "llt_g": (([[1], [1]],), qt_element),
    "llt_g_lt": (([3, 1], 2), qt_element),
    "llt_graph": ((3, [(0, 1)], []), qt_element),
    "llt_gtilde": (([3, 1], 2), qt_element),
    "llt_gtilde_table": ((4, 2), keyed_by_partition(qt_element)),
    "llt_h": (([3, 1], 2), qt_element),
    "llt_h_table": ((4, 2), keyed_by_partition(qt_element)),
    "llt_h_tilde": (([3, 1], 2), qt_element),
    "llt_kl_column": (([2], 2), qt_element),
    "llt_min_inv": (([[1], [1]],), is_int),
    "llt_schur": (([2], 2), qt_element),
    "lr_coefficient": (([3, 1], [2, 1], [1]), is_int),
    "macdonald_ht": (([2, 1],), qt_element),
    "schur_to_macdonald_ht": (([([2, 1], [(0, 0, 1), (1, 2, -3)])],), ht_element),
    "macdonald_ht_to_schur": (
        (([([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])],)),
        ht_element,
    ),
    "macdonald_ht_element_add": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])], [([2, 1], [(0, 0, 2)], [(1, 1, 1, 1)])]),
        ht_element,
    ),
    "macdonald_ht_element_scale": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1, 1)])], [(1, 0, 1)], []),
        ht_element,
    ),
    "schur_to_macdonald_j": (([([2, 1], [(0, 0, 1), (1, 2, -3)])],), macdonald_element),
    "macdonald_j": (([2, 1],), macdonald_element),
    "macdonald_p": (([2, 1],), macdonald_element),
    "macdonald_q": (([2, 1],), macdonald_element),
    "monomial_to_macdonald_p": (
        (([([2, 1], [(0, 0, 1)], [(1, 1, 1)])],)),
        macdonald_element,
    ),
    "monomial_to_macdonald_q": (
        (([([2, 1], [(0, 0, 1)], [(1, 1, 1)])],)),
        macdonald_element,
    ),
    "monomial_to_jack_p": ((([([2, 1], [1, -2], [(1, 1, 1)], 3, [])],)), jack_element),
    "monomial_to_jack_q": ((([([2, 1], [1, -2], [(1, 1, 1)], 3, [])],)), jack_element),
    "monomial_to_jack_j": ((([([2, 1], [1, -2], [(1, 1, 1)], 3, [])],)), jack_element),
    "jack_p_to_monomial": ((([([2, 1], [1, -2], [(1, 1, 1)], 3, [])],)), jack_element),
    "jack_q_to_monomial": ((([([2, 1], [1, -2], [(1, 1, 1)], 3, [])],)), jack_element),
    "jack_j_to_monomial": ((([([2, 1], [1, -2], [(1, 1, 1)], 3, [])],)), jack_element),
    "macdonald_p_to_monomial": (
        (([([2, 1], [(0, 0, 1)], [(1, 1, 1)])],)),
        macdonald_element,
    ),
    "macdonald_q_to_monomial": (
        (([([2, 1], [(0, 0, 1)], [(1, 1, 1)])],)),
        macdonald_element,
    ),
    "macdonald_j_to_monomial": (
        (([([2, 1], [(0, 0, 1)], [(1, 1, 1)])],)),
        macdonald_element,
    ),
    "macdonald_element_add": (
        (
            [([2, 1], [(0, 0, 1)], [(1, 1, 1)])],
            [([2, 1], [(0, 0, 2)], [(1, 1, 1)])],
        ),
        macdonald_element,
    ),
    "macdonald_element_scale": (
        ([([2, 1], [(0, 0, 1)], [(1, 1, 1)])], [(1, 0, 1)], []),
        macdonald_element,
    ),
    "jack_element_add": (
        ([([2, 1], [1, -2], [(1, 1, 1)], 3, [])], [([2, 1], [2], [(1, 1, 1)], 3, [])]),
        jack_element,
    ),
    "jack_element_scale": (
        ([([2, 1], [1, -2], [(1, 1, 1)], 3, [])], [0, 1], [], 1),
        jack_element,
    ),
    "hall_littlewood_p_to_schur": (([([2, 1], [(0, 1), (2, -3)])],), t_element),
    "hall_littlewood_qp_to_schur": (([([2, 1], [(0, 1), (2, -3)])],), t_element),
    "monomial_multiply": ((A, B), element),
    "monomial_to_schur": ((A,), element),
    "nabla": ((QT_A,), qt_element),
    "nabla_e": ((3,), qt_element),
    "nabla_e_by_path": ((3,), path_refinement),
    "nabla_power": ((QT_A, 2), qt_element),
    "omega": ((A,), element),
    "partitions": ((4,), list_of(is_partition)),
    "plethysm": ((A, B), element),
    "polynomial_to_schubert": (([([1], 1)],), schubert_element),
    "power_to_schur": ((A,), element),
    "principal_specialization": (([2, 1], 3), optional(is_int)),
    "principal_specialization_q": (([2, 1], 3), is_int_list),
    "qt_kostka": (([2, 1], [1, 1, 1]), is_qt_coeff),
    "qt_kostka_column": (([1, 1, 1],), qt_element),
    "qt_kostka_table": ((3,), list_of(list_of(is_qt_coeff))),
    "reduced_kronecker": (([2, 1], [2, 1], [2, 1]), is_int),
    "reduced_kronecker_product": (([2, 1], [2, 1]), element),
    "schubert_coefficient": (([2, 1], [2, 1], [3, 1, 2]), is_int),
    "schubert_dimension": (([2, 1, 4, 3],), is_int),
    "schubert_divided_difference": ((SCHUB, 1), schubert_element),
    "schubert_divided_difference_perm": ((SCHUB, [2, 1]), schubert_element),
    "schubert_expand": ((SCHUB,), monomials),
    "schubert_monomial_mass": (([2, 1], [1, 3, 2]), is_int),
    "schubert_multiply": ((SCHUB, SCHUB), schubert_element),
    "schubert_multiply_variable": ((SCHUB, 1), schubert_element),
    "schubert_pairing": ((SCHUB, SCHUB, 3), is_int),
    "schubert_scalar_product": ((SCHUB, SCHUB, 3), schubert_element),
    "schubert_to_stanley_schur": (([2, 1, 4, 3],), element),
    "schur_in_macdonald_j": ((3,), keyed_by_partition(macdonald_element)),
    "schur_multiply": ((A, B), element),
    "schur_to_elementary": ((A,), element),
    "schur_to_forgotten": ((A,), element),
    "schur_to_homogeneous": ((A,), element),
    "schur_to_ht": ((A,), element),
    "schur_to_monomial": ((A,), element),
    "schur_to_st": ((A,), element),
    "semistandard_tableaux": (([2, 1], [1, 1, 1]), tableaux),
    "skew_by": ((A, B, "s"), element),
    "skew_schur": (([2, 1], [1]), element),
    "st_multiply": ((A, B), element),
    "st_to_schur": ((A,), element),
    "stanley_table": ((2,), stanley_rows),
    "theta_ek": ((1, QT_A), qt_element),
    "to_power": ((A, "s"), rational_element),
    "zonal": (([2, 1], False), zonal_terms),
}


def exported(mod):
    return {n for n in dir(mod) if not n.startswith("_") and callable(getattr(mod, n))}


def check_surface_is_covered(mod):
    """Every exported callable has a shape row — an unclassified encoding
    cannot ship."""
    missing = sorted(exported(mod) - set(SHAPES))
    stale = sorted(set(SHAPES) - exported(mod))
    out = []
    if missing:
        out.append(f"these pyfunctions have no shape row: {', '.join(missing)}")
    if stale:
        out.append(f"these names are listed but not exported: {', '.join(stale)}")
    return out


def check_shapes(mod):
    checked = 0
    out = []
    for name, (args, ok) in sorted(SHAPES.items()):
        checked += 1
        try:
            got = getattr(mod, name)(*args)
        except BaseException as e:  # noqa: BLE001 - a raise on valid input is the finding
            out.append(f"{name}{args}: raised {type(e).__name__} on valid input: {e}")
            continue
        if not ok(got):
            shown = repr(got)
            if len(shown) > 160:
                shown = shown[:160] + "..."
            out.append(f"{name}{args}: return violates the promised shape: {shown}")
    return checked, out


# --- widths -----------------------------------------------------------------
#
# The module doc's claim: coefficients cross as Python ints of arbitrary size,
# in both directions, with an i128 fast path and a BigInt escalation underneath.
# Both edges of the fast path and values past it must come back unchanged —
# a wrapped coefficient here is exactly the silent-wrong outcome
# docs/policies/failure.md R1 forbids. The identity-shaped calls keep the
# mathematics trivial so any change in the value is the marshalling's.

WIDTHS = [
    1,
    -1,
    2**63 - 1,  # i64::MAX: the widest a older fast path would have offered
    -(2**63),
    2**64,  # past u64
    2**127 - 1,  # i128::MAX: the last value the fast path holds
    -(2**127),  # i128::MIN
    2**127,  # the first value that must take the BigInt path
    10**40,
    -(2**200),
]


def check_widths(mod):
    checked = 0
    out = []
    for v in WIDTHS:
        cases = [
            ("convert_terms", ([([1], v)], "s", "s"), [((1,), v)]),
            (
                "convert_qt_terms",
                ([([1], [(0, 1, v)])], "s", "s"),
                [((1,), [(0, 1, v)])],
            ),
            # (2) is self-conjugate under neither, but (1, 1) comes back at the
            # same width; the antipode's sign is even here, so both copy.
            ("omega_qt_terms", ([([1, 1], [(0, 1, v)])], ), [((2,), [(0, 1, v)])]),
            # Times the identity, so the coefficient crosses the product path
            # at full width rather than being multiplied down.
            (
                "schur_multiply_qt",
                ([([2, 1], [(0, 1, v)])], [([], [(0, 0, 1)])]),
                [((2, 1), [(0, 1, v)])],
            ),
            (
                "antipode_qt_terms",
                ([([1, 1], [(0, 1, v)])],),
                [((2,), [(0, 1, v)])],
            ),
            ("schur_multiply", ([([2, 1], v)], [([], 1)]), [((2, 1), v)]),
            ("st_multiply", ([([1], v)], [([], 1)]), [((1,), v)]),
            # ω fixes s_1, so the answer is the input — through the omega path.
            ("omega", ([([1], v)],), [((1,), v)]),
            # Multiplying by the identity permutation exercises the Schubert
            # marshalling both ways.
            ("schubert_multiply", ([([2, 1], v)], [([1], 1)]), [((2, 1), v)]),
            # v·s_1 = v·p_1: the rational dump's (numerator, denominator).
            ("to_power", ([([1], v)], "s"), [((1,), (v, 1))]),
            ("evaluate_schur", ([([1], v)], [1]), v),
        ]
        # ∇s_1 = s_1, so the (q,t) inbound path is identity here — but it is
        # walled at the fixed width by policy (tested below), so only in-range
        # values: i128 minus `i128::MIN`, which has no negation in the width.
        if -(2**127) < v < 2**127:
            cases.append(("nabla", ([([1], [(0, 0, v)])],), [((1,), [(0, 0, v)])]))
        for name, args, want in cases:
            checked += 1
            try:
                got = getattr(mod, name)(*args)
            except BaseException as e:  # noqa: BLE001
                out.append(f"{name} at {v}: raised {type(e).__name__}: {e}")
                continue
            if got != want:
                out.append(f"{name} at {v}: got {got}, want {want}")
    # The (q,t) families take coefficients through fixed width, and outside it
    # the outcome is a typed refusal — never a wrapped value, never a panic.
    # Both ends of the window: past i128, and i128::MIN, which fits the type
    # but has no negation in it and once panicked at the first sign flip.
    for v in (2**127, -(2**127)):
        checked += 1
        try:
            got = mod.nabla([([1], [(0, 0, v)])])
            out.append(f"nabla at {v}: returned {got} instead of raising")
        except ValueError:
            pass
        except BaseException as e:  # noqa: BLE001
            out.append(f"nabla at {v}: {type(e).__name__}, not ValueError: {e}")
    return checked, out


# --- permissive inbound -----------------------------------------------------


def check_permissive_inbound(mod):
    """The `*Arg` halves of the stub aliases: any sequence, padding tolerated.

    Every spelling of the same value must produce the identical answer — and
    the strict outbound form handed straight back in must be among them, which
    is the round trip a caller who post-processes an output actually performs.
    """
    checked = 0
    out = []
    groups = [
        # One support, four spellings.
        ("kostka_number", [
            (([2, 1], [1, 1, 1]),),
            (((2, 1), (1, 1, 1)),),
            (([2, 1, 0, 0], [1, 1, 1, 0]),),
            (((2, 1, 0), [1, 1, 1]),),
        ]),
        # An element as list-of-pairs and tuple-of-pairs, padded and not.
        ("schur_to_monomial", [
            (([([2, 1], 3), ([1, 1, 1], -1)],),),
            (((([2, 1, 0], 3), ((1, 1, 1), -1)),),),
        ]),
        # A permutation padded with trailing fixed points.
        ("schubert_multiply", [
            (([([2, 1], 1)], [([1, 3, 2], 1)]),),
            (([([2, 1, 3, 4], 1)], [((1, 3, 2, 4), 1)]),),
        ]),
        # A (q,t) element with its sequences as lists and as tuples — the
        # triples themselves stay tuples, which is what the stub promises.
        ("nabla", [
            (([([2], [(0, 0, 1)]), ([1, 1], [(1, 0, 2)])],),),
            ((((((2,)), ((0, 0, 1),)), (((1, 1)), ((1, 0, 2),))),),),
        ]),
    ]
    for name, spellings in groups:
        fn = getattr(mod, name)
        answers = []
        for (args,) in spellings:
            checked += 1
            try:
                answers.append(fn(*args))
            except BaseException as e:  # noqa: BLE001
                out.append(f"{name}{args}: raised {type(e).__name__}: {e}")
        if len(set(map(repr, answers))) > 1:
            out.append(f"{name}: spellings of one value disagree: {answers}")
    # The literal round trip: an output fed back in, unmodified.
    checked += 1
    once = mod.schur_multiply([([2, 1], 3)], [([1], 1)])
    again = mod.schur_multiply(once, [([], 1)])
    if once != again:
        out.append(f"schur_multiply: output refused as input: {once} -> {again}")
    return checked, out


def main():
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {sys.argv[0]} <path to built symfn extension>")
    mod = load(sys.argv[1])

    failures = check_surface_is_covered(mod)
    shapes, bad = check_shapes(mod)
    failures += bad
    widths, bad = check_widths(mod)
    failures += bad
    perm, bad = check_permissive_inbound(mod)
    failures += bad

    if failures:
        total = shapes + widths + perm
        print(f"FAIL: {len(failures)} problem(s) across {total} checks\n")
        for f in failures:
            print(f"  {f}")
        raise SystemExit(1)
    print(
        f"ok: {shapes} returns in the promised shapes, {widths} width "
        f"round-trips, {perm} inbound spellings"
    )


if __name__ == "__main__":
    main()
