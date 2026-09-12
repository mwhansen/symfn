"""Feed malformed input to every `#[pyfunction]` and demand a typed exception.

    cargo build --features python
    python3 scripts/check_python_boundary.py target/debug/libsymfn.dylib

The proposition: **no argument a Python caller can construct produces a
`PanicException`.** `docs/policies/failure.md` R2 puts it the other way round —
a `PanicException` surfacing in Sage is by definition a bug report, the crate's
bug or evidence that a wall needs a real mechanism, never an interface — so
what this script asserts is that the boundary has no such interface left.

Malformed here means *correctly typed and mathematically impossible*: PyO3's
own extraction succeeds, and what fails is a precondition. A `[1, 3]` where a
partition belongs, a `[1, 1]` where a permutation belongs, `k = 0` for a ribbon
level, mixed degrees where an operator needs one. Wrong *types* are PyO3's job
and are not tested here; wrong *values* are this boundary's job and are.

Needs no Sage — it imports the built extension module and nothing else. Run it
after touching `src/python.rs`; `scripts/preflight.sh` cannot, since it builds
only the default features and this needs `--features python`.

Two cases carry their own history:

  * `convert_indexed([([2, 1, 0], 1)], "Schur", "Schur")` is a **valid**
    partition, padded — the fixed-width lists Sage hands over. The identity
    conversion passed the raw list to a table keyed by normal forms, so the one
    input shape the module documents as tolerated was the one that crashed.
  * every Schubert entry point taking a term list reached
    `build_schubert(..).unwrap()` on the escalation path. The fast pass declined
    a malformed word by returning `None`, so the only way to see the panic was
    to pass a coefficient that fits — i.e. almost always.
"""

import importlib.machinery
import importlib.util
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


# A partition argument that is not a partition, and one that is valid but
# padded — the boundary must reject the first and accept the second.
BAD_PART = [1, 3]
BAD_PART_ZERO = [2, 0, 1]
BAD_TERMS = [(BAD_PART, 1)]
# A one-line word that is not a bijection, with a coefficient small enough to
# take the fixed-width pass — which is what made the old panic reachable.
BAD_SCHUB = [([1, 1], 1)]
BAD_QTSCHUR = [([2], [(0, 0, 1)]), ([1, 1, 1], [(0, 0, 1)])]

# name -> [(label, args), ...]. Every exported callable must appear here or in
# TOTAL below; `test_surface_is_covered` fails if a new one is added without a
# case, which is why the table is written out rather than introspected.
CASES = {
    # --- partitions crossing as bare arguments ---
    "lr_coefficient": [("lambda", (BAD_PART, [1], [2]))],
    "kostka_number": [("lambda", (BAD_PART, [2, 2]))],
    "character_value": [
        ("lambda", (BAD_PART, [2, 2])),
        ("off-degree", ([2, 1], [2, 2])),
    ],
    "kronecker_coefficient": [
        ("lambda", (BAD_PART, [2, 2], [2, 2])),
        ("off-degree", ([2, 1], [2, 2], [2, 2])),
    ],
    "class_algebra_coefficient": [
        ("lambda", (BAD_PART, [2, 2], [2, 2])),
        ("off-degree", ([2, 1], [2, 2], [2, 2])),
    ],
    "semistandard_tableaux": [("lambda", (BAD_PART, [2, 2]))],
    "dimension": [("lambda", (BAD_PART,))],
    "principal_specialization": [("lambda", (BAD_PART, 3))],
    "principal_specialization_q": [("lambda", (BAD_PART, 3))],
    "skew_schur": [("lambda", (BAD_PART, [1]))],
    "hall_littlewood": [("lambda", (BAD_PART,))],
    "hall_littlewood_p": [("lambda", (BAD_PART,))],
    "schur_to_hall_littlewood_p": [("lambda", ([(BAD_PART, [(0, 1)])],))],
    "schur_to_hall_littlewood_qp": [("lambda", ([(BAD_PART, [(0, 1)])],))],
    "schur_to_macdonald_ht": [("lambda", ([(BAD_PART, [(0, 0, 1)])],))],
    "schur_to_macdonald_j": [("lambda", ([(BAD_PART, [(0, 0, 1)])],))],
    "monomial_to_macdonald_p": [
        ("lambda", ([(BAD_PART, [(0, 0, 1)], [])],)),
        # (0, 0) is `1 - q^0 t^0 = 0`, so a zero denominator rather than a
        # factor; `Frac::mul_factors` asserts on it and the boundary must
        # raise first.
        ("zero denominator factor", ([([2], [(0, 0, 1)], [(0, 0, 1)])],)),
    ],
    "monomial_to_macdonald_q": [
        ("lambda", ([(BAD_PART, [(0, 0, 1)], [])],)),
        ("zero denominator factor", ([([2], [(0, 0, 1)], [(0, 0, 1)])],)),
    ],
    "monomial_to_jack_p": [
        ("lambda", ([(BAD_PART, [1], [], 1)],)),
        # (0, 0) is the zero linear form, so a zero denominator rather than a
        # factor; `AFrac::mul_factors` asserts on it and the boundary must
        # raise first.
        ("zero denominator atom", ([([2], [1], [(0, 0, 1)], 1)],)),
        # A scale of 0 is the same trap in the other field of the row.
        ("zero scale", ([([2], [1], [], 0)],)),
    ],
    "monomial_to_jack_q": [
        ("lambda", ([(BAD_PART, [1], [], 1)],)),
        ("zero denominator atom", ([([2], [1], [(0, 0, 1)], 1)],)),
    ],
    "monomial_to_jack_j": [
        ("lambda", ([(BAD_PART, [1], [], 1)],)),
        ("zero scale", ([([2], [1], [], 0)],)),
    ],
    "macdonald_ht_to_schur": [
        ("lambda", ([(BAD_PART, [(0, 0, 1)], [])],)),
        # (0, 0) of kind 0 is `1 - q^0 t^0 = 0`, a zero denominator.
        ("zero atom", ([([2], [(0, 0, 1)], [(0, 0, 0, 1)])],)),
        # `q^0 - t^b` is `1 - t^b`: it belongs to kind 0, and letting it in
        # under kind 1 would leave two spellings of one polynomial.
        ("degenerate kind 1", ([([2], [(0, 0, 1)], [(1, 0, 1, 1)])],)),
        ("unknown kind", ([([2], [(0, 0, 1)], [(2, 1, 1, 1)])],)),
    ],
    "macdonald_ht_element_add": [
        ("lambda", ([(BAD_PART, [(0, 0, 1)], [])], [])),
        ("unknown kind", ([([2], [(0, 0, 1)], [(2, 1, 1, 1)])], [])),
    ],
    "macdonald_ht_element_scale": [
        ("lambda", ([(BAD_PART, [(0, 0, 1)], [])], [(0, 0, 1)], [])),
        ("unknown kind", ([([2], [(0, 0, 1)], [])], [(0, 0, 1)], [(2, 1, 1, 1)])),
    ],
    "macdonald_element_add": [
        ("lambda", ([(BAD_PART, [(0, 0, 1)], [])], [])),
        ("lambda in the second", ([], [(BAD_PART, [(0, 0, 1)], [])])),
    ],
    "macdonald_element_scale": [
        ("lambda", ([(BAD_PART, [(0, 0, 1)], [])], [(0, 0, 1)], [])),
        ("zero denominator factor", ([([2], [(0, 0, 1)], [])], [(0, 0, 1)], [(0, 0, 1)])),
    ],
    "jack_element_add": [
        ("lambda", ([(BAD_PART, [1], [], 1)], [])),
        ("zero scale", ([([2], [1], [], 0)], [])),
    ],
    "jack_element_scale": [
        ("lambda", ([(BAD_PART, [1], [], 1)], [1], [], 1)),
        ("zero scale in the scalar", ([([2], [1], [], 1)], [1], [], 0)),
    ],
    "hall_littlewood_p_to_schur": [("lambda", ([(BAD_PART, [(0, 1)])],))],
    "hall_littlewood_qp_to_schur": [("lambda", ([(BAD_PART, [(0, 1)])],))],
    "macdonald_p_to_monomial": [
        ("lambda", ([(BAD_PART, [(0, 0, 1)], [])],)),
        ("zero denominator factor", ([([2], [(0, 0, 1)], [(0, 0, 1)])],)),
    ],
    "macdonald_q_to_monomial": [
        ("lambda", ([(BAD_PART, [(0, 0, 1)], [])],)),
        ("zero denominator factor", ([([2], [(0, 0, 1)], [(0, 0, 1)])],)),
    ],
    "macdonald_j_to_monomial": [
        ("lambda", ([(BAD_PART, [(0, 0, 1)], [])],)),
        ("zero denominator factor", ([([2], [(0, 0, 1)], [(0, 0, 1)])],)),
    ],
    "jack_p_to_monomial": [
        ("lambda", ([(BAD_PART, [1], [], 1)],)),
        ("zero denominator atom", ([([2], [1], [(0, 0, 1)], 1)],)),
        ("zero scale", ([([2], [1], [], 0)],)),
    ],
    "jack_q_to_monomial": [
        ("lambda", ([(BAD_PART, [1], [], 1)],)),
        ("zero denominator atom", ([([2], [1], [(0, 0, 1)], 1)],)),
    ],
    "jack_j_to_monomial": [
        ("lambda", ([(BAD_PART, [1], [], 1)],)),
        ("zero scale", ([([2], [1], [], 0)],)),
    ],
    "kostka_foulkes": [
        ("lambda", (BAD_PART, [2, 2])),
        ("off-degree", ([2, 1], [2, 2])),
    ],
    "kostka_foulkes_column": [("mu", (BAD_PART,))],
    "macdonald_p": [("lambda", (BAD_PART,))],
    "macdonald_q": [("lambda", (BAD_PART,))],
    "macdonald_j": [("lambda", (BAD_PART,))],
    "macdonald_ht": [("mu", (BAD_PART,))],
    "qt_kostka": [
        ("lambda", (BAD_PART, [2, 2])),
        ("off-degree", ([2, 1], [2, 2])),
    ],
    "qt_kostka_column": [("mu", (BAD_PART,))],
    "jack_p": [("lambda", (BAD_PART,))],
    "jack_q": [("lambda", (BAD_PART,))],
    "jack_j": [("lambda", (BAD_PART,))],
    "jack_j_powersum": [("lambda", (BAD_PART,))],
    "jack_norm_j": [("lambda", (BAD_PART,))],
    "jack_structure_constant": [("lambda", (BAD_PART, [1], [2, 1]))],
    "jack_scalar": [("f", ([(BAD_PART, [1])], [([1], [1])]))],
    "zonal": [("lambda", (BAD_PART, False))],
    # --- elements crossing as (partition, coefficient) lists ---
    "schur_multiply": [
        ("a", (BAD_TERMS, [([1], 1)])),
        ("interior zero", ([(BAD_PART_ZERO, 1)], [([1], 1)])),
    ],
    "monomial_multiply": [("a", (BAD_TERMS, [([1], 1)]))],
    "st_multiply": [("a", (BAD_TERMS, [([1], 1)]))],
    "ht_multiply": [("a", (BAD_TERMS, [([1], 1)]))],
    "schur_to_st": [("a", (BAD_TERMS,))],
    "st_to_schur": [("a", (BAD_TERMS,))],
    "schur_to_ht": [("a", (BAD_TERMS,))],
    "ht_to_schur": [("a", (BAD_TERMS,))],
    "reduced_kronecker_product": [("lambda", (BAD_PART, [1]))],
    "reduced_kronecker": [("lambda", (BAD_PART, [1], [1]))],
    "schur_to_homogeneous": [("a", (BAD_TERMS,))],
    "schur_to_elementary": [("a", (BAD_TERMS,))],
    "schur_to_monomial": [("a", (BAD_TERMS,))],
    "schur_to_forgotten": [("a", (BAD_TERMS,))],
    "to_power": [
        ("a", (BAD_TERMS, "Schur")),
        ("src", ([([2], 1)], "zzz")),
    ],
    "homogeneous_to_schur": [("a", (BAD_TERMS,))],
    "elementary_to_schur": [("a", (BAD_TERMS,))],
    "monomial_to_schur": [("a", (BAD_TERMS,))],
    "forgotten_to_schur": [("a", (BAD_TERMS,))],
    "power_to_schur": [("a", (BAD_TERMS,))],
    "omega": [("a", (BAD_TERMS,))],
    "antipode": [("a", (BAD_TERMS,))],
    "coproduct": [("a", (BAD_TERMS,))],
    "hall_inner_product": [("a", (BAD_TERMS, [([1], 1)]))],
    "internal_product": [("a", (BAD_TERMS, [([1], 1)]))],
    "plethysm": [("f", (BAD_TERMS, [([1], 1)]))],
    "evaluate_schur": [("a", (BAD_TERMS, [1, 2]))],
    "skew_by": [
        ("f", (BAD_TERMS, [([1], 1)], "s")),
        ("basis", ([([2], 1)], [([1], 1)], "zzz")),
    ],
    "expand_alphabet": [
        ("a", (BAD_TERMS, "Schur", 3)),
        ("src", ([([2], 1)], "zzz", 3)),
    ],
    "convert_indexed": [
        ("a", (BAD_TERMS, "Schur", "Schur")),
        ("src", ([([2], 1)], "zzz", "Schur")),
        ("dst", ([([2], 1)], "Schur", "zzz")),
    ],
    "convert_terms": [
        ("a", (BAD_TERMS, "Schur", "Schur")),
        ("src", ([([2], 1)], "zzz", "Schur")),
        ("dst", ([([2], 1)], "Schur", "zzz")),
        ("dst powersum", ([([2], 1)], "Schur", "powersum")),
        ("dst p", ([([2], 1)], "s", "p")),
    ],
    "convert_macdonald_terms": [
        ("a", ([(BAD_PART, [(0, 0, 1)], [])], "Schur", "Schur")),
        ("src", ([([2], [(0, 0, 1)], [])], "zzz", "Schur")),
        ("dst", ([([2], [(0, 0, 1)], [])], "Schur", "zzz")),
        ("dst p", ([([2], [(0, 0, 1)], [])], "s", "p")),
    ],
    "convert_jack_terms": [
        ("a", ([(BAD_PART, [1], [], 1)], "Schur", "Schur")),
        ("src", ([([2], [1], [], 1)], "zzz", "Schur")),
        ("dst", ([([2], [1], [], 1)], "Schur", "zzz")),
        ("dst p", ([([2], [1], [], 1)], "s", "p")),
    ],
    "convert_ht_terms": [
        ("a", ([(BAD_PART, [(0, 0, 1)], [])], "Schur", "Schur")),
        ("src", ([([2], [(0, 0, 1)], [])], "zzz", "Schur")),
        ("dst", ([([2], [(0, 0, 1)], [])], "Schur", "zzz")),
        ("dst p", ([([2], [(0, 0, 1)], [])], "s", "p")),
    ],
    "schur_multiply_macdonald": [
        ("a", ([(BAD_PART, [(0, 0, 1)], [])], [([1], [(0, 0, 1)], [])])),
        ("b", ([([1], [(0, 0, 1)], [])], [(BAD_PART, [(0, 0, 1)], [])])),
    ],
    "schur_multiply_jack": [
        ("a", ([(BAD_PART, [1], [], 1)], [([1], [1], [], 1)])),
        ("b", ([([1], [1], [], 1)], [(BAD_PART, [1], [], 1)])),
    ],
    "schur_multiply_ht": [
        ("a", ([(BAD_PART, [(0, 0, 1)], [])], [([1], [(0, 0, 1)], [])])),
        ("b", ([([1], [(0, 0, 1)], [])], [(BAD_PART, [(0, 0, 1)], [])])),
    ],
    "schur_multiply_qt": [
        ("a", ([(BAD_PART, [(0, 0, 1)])], [([1], [(0, 0, 1)])])),
        ("b", ([([1], [(0, 0, 1)])], [(BAD_PART, [(0, 0, 1)])])),
    ],
    "omega_qt_terms": [("a", ([(BAD_PART, [(0, 0, 1)])],))],
    "antipode_qt_terms": [("a", ([(BAD_PART, [(0, 0, 1)])],))],
    "omega_macdonald_terms": [("a", ([(BAD_PART, [(0, 0, 1)], [])],))],
    "antipode_macdonald_terms": [("a", ([(BAD_PART, [(0, 0, 1)], [])],))],
    "omega_jack_terms": [("a", ([(BAD_PART, [1], [], 1)],))],
    "antipode_jack_terms": [("a", ([(BAD_PART, [1], [], 1)],))],
    "omega_ht_terms": [("a", ([(BAD_PART, [(0, 0, 1)], [])],))],
    "antipode_ht_terms": [("a", ([(BAD_PART, [(0, 0, 1)], [])],))],
    "skew_by_qt": [
        ("f", ([(BAD_PART, [(0, 0, 1)])], [([1], [(0, 0, 1)])], "s")),
        ("g", ([([1], [(0, 0, 1)])], [(BAD_PART, [(0, 0, 1)])], "s")),
        ("basis", ([([1], [(0, 0, 1)])], [([1], [(0, 0, 1)])], "zzz")),
    ],
    "skew_by_macdonald": [
        ("f", ([(BAD_PART, [(0, 0, 1)], [])], [([1], [(0, 0, 1)], [])], "s")),
        ("basis", ([([1], [(0, 0, 1)], [])], [([1], [(0, 0, 1)], [])], "zzz")),
    ],
    "skew_by_jack": [
        ("f", ([(BAD_PART, [1], [], 1)], [([1], [1], [], 1)], "s")),
        ("basis", ([([1], [1], [], 1)], [([1], [1], [], 1)], "zzz")),
    ],
    "skew_by_ht": [
        ("f", ([(BAD_PART, [(0, 0, 1)], [])], [([1], [(0, 0, 1)], [])], "s")),
        ("basis", ([([1], [(0, 0, 1)], [])], [([1], [(0, 0, 1)], [])], "zzz")),
    ],
    "hall_inner_product_qt": [
        ("a", ([(BAD_PART, [(0, 0, 1)])], [([1], [(0, 0, 1)])])),
        ("b", ([([1], [(0, 0, 1)])], [(BAD_PART, [(0, 0, 1)])])),
    ],
    "hall_inner_product_macdonald": [
        ("a", ([(BAD_PART, [(0, 0, 1)], [])], [([1], [(0, 0, 1)], [])])),
    ],
    "hall_inner_product_jack": [
        ("a", ([(BAD_PART, [1], [], 1)], [([1], [1], [], 1)])),
    ],
    "hall_inner_product_ht": [
        ("a", ([(BAD_PART, [(0, 0, 1)], [])], [([1], [(0, 0, 1)], [])])),
    ],
    "scalar_t": [
        ("a", ([(BAD_PART, [(0, 0, 1)], [])], [([1], [(0, 0, 1)], [])])),
        ("b", ([([1], [(0, 0, 1)], [])], [(BAD_PART, [(0, 0, 1)], [])])),
    ],
    "scalar_qt": [
        ("a", ([(BAD_PART, [(0, 0, 1)], [])], [([1], [(0, 0, 1)], [])])),
        ("b", ([([1], [(0, 0, 1)], [])], [(BAD_PART, [(0, 0, 1)], [])])),
    ],
    "scalar_qt_ht": [
        ("a", ([(BAD_PART, [(0, 0, 1)], [])], [([1], [(0, 0, 1)], [])])),
        ("atom", ([([1], [(0, 0, 1)], [(2, 0, 0, 1)])], [([1], [(0, 0, 1)], [])])),
    ],
    "to_power_qt": [
        ("f", ([(BAD_PART, [(0, 0, 1)])], "s")),
        ("src", ([([2], [(0, 0, 1)])], "zzz")),
    ],
    "to_power_macdonald": [
        ("f", ([(BAD_PART, [(0, 0, 1)], [])], "s")),
        ("src", ([([2], [(0, 0, 1)], [])], "zzz")),
        ("den", ([([2], [(0, 0, 1)], [(0, 0, 1)])], "s")),
    ],
    "to_power_jack": [
        ("f", ([(BAD_PART, [1], [], 1, [])], "s")),
        ("src", ([([2], [1], [], 1, [])], "zzz")),
        ("scale", ([([2], [1], [], 0, [])], "s")),
    ],
    "to_power_ht": [
        ("f", ([(BAD_PART, [(0, 0, 1)], [])], "s")),
        ("src", ([([2], [(0, 0, 1)], [])], "zzz")),
    ],
    "coproduct_qt": [("a", ([(BAD_PART, [(0, 0, 1)])],))],
    "coproduct_macdonald": [("a", ([(BAD_PART, [(0, 0, 1)], [])],))],
    "coproduct_jack": [("a", ([(BAD_PART, [1], [], 1)],))],
    "coproduct_ht": [("a", ([(BAD_PART, [(0, 0, 1)], [])],))],
    "expand_qt": [("a", ([(BAD_PART, [(0, 0, 1)])], 2))],
    "expand_macdonald": [("a", ([(BAD_PART, [(0, 0, 1)], [])], 2))],
    "expand_jack": [("a", ([(BAD_PART, [1], [], 1)], 2))],
    "expand_ht": [("a", ([(BAD_PART, [(0, 0, 1)], [])], 2))],
    "evaluate_qt": [("a", ([(BAD_PART, [(0, 0, 1)])], [1, 1]))],
    "evaluate_macdonald": [("a", ([(BAD_PART, [(0, 0, 1)], [])], [1, 1]))],
    "evaluate_jack": [("a", ([(BAD_PART, [1], [], 1)], [1, 1]))],
    "evaluate_ht": [("a", ([(BAD_PART, [(0, 0, 1)], [])], [1, 1]))],
    "dimension_qt": [("a", ([(BAD_PART, [(0, 0, 1)])],))],
    "dimension_macdonald": [("a", ([(BAD_PART, [(0, 0, 1)], [])],))],
    "dimension_jack": [("a", ([(BAD_PART, [1], [], 1)],))],
    "dimension_ht": [("a", ([(BAD_PART, [(0, 0, 1)], [])],))],
    "principal_specialization_qt": [("a", ([(BAD_PART, [(0, 0, 1)])], 3))],
    "principal_specialization_macdonald": [("a", ([(BAD_PART, [(0, 0, 1)], [])], 3))],
    "principal_specialization_jack": [("a", ([(BAD_PART, [1], [], 1)], 3))],
    "principal_specialization_ht": [("a", ([(BAD_PART, [(0, 0, 1)], [])], 3))],
    "principal_specialization_q_qt": [
        ("a", ([(BAD_PART, [(0, 0, 1)])], 3)),
        ("q slot in use", ([([2], [(1, 0, 1)])], 3)),
    ],
    "internal_product_qt": [
        ("a", ([(BAD_PART, [(0, 0, 1)])], [([1], [(0, 0, 1)])])),
        ("b", ([([1], [(0, 0, 1)])], [(BAD_PART, [(0, 0, 1)])])),
    ],
    "internal_product_macdonald": [
        ("a", ([(BAD_PART, [(0, 0, 1)], [])], [([1], [(0, 0, 1)], [])])),
    ],
    "internal_product_jack": [
        ("a", ([(BAD_PART, [1], [], 1)], [([1], [1], [], 1)])),
    ],
    "internal_product_ht": [
        ("a", ([(BAD_PART, [(0, 0, 1)], [])], [([1], [(0, 0, 1)], [])])),
    ],
    "principal_specialization_at_qt": [
        ("a", ([(BAD_PART, [(0, 0, 1)])], 3, [(0, 0, 1)])),
        ("z", ([([2], [(0, 0, 1)])], 3, [(0, 0, "x")])),
    ],
    "principal_specialization_at_macdonald": [
        ("a", ([(BAD_PART, [(0, 0, 1)], [])], 3, ([(0, 0, 1)], []))),
        ("z", ([([2], [(0, 0, 1)], [])], 3, ([(0, 0, 1)], [(0, 0, 1)]))),
    ],
    "principal_specialization_at_jack": [
        ("a", ([(BAD_PART, [1], [], 1)], 3, ([1], [], 1))),
        ("z scale", ([([2], [1], [], 1)], 3, ([1], [], 0))),
    ],
    "principal_specialization_at_ht": [
        ("a", ([(BAD_PART, [(0, 0, 1)], [])], 3, ([(0, 0, 1)], []))),
    ],
    "plethysm_qt": [
        ("f", ([(BAD_PART, [(0, 0, 1)])], [([1], [(0, 0, 1)])])),
        ("g", ([([1], [(0, 0, 1)])], [(BAD_PART, [(0, 0, 1)])])),
    ],
    "plethysm_jack": [
        ("f", ([(BAD_PART, [1], [], 1, [])], [([1], [1], [], 1, [])])),
    ],
    "plethysm_macdonald": [
        ("f", ([(BAD_PART, [(0, 0, 1)], [])], [([1], [(0, 0, 1)], [])])),
    ],
    "plethysm_ht": [
        ("f", ([(BAD_PART, [(0, 0, 1)], [])], [([1], [(0, 0, 1)], [])])),
    ],
    "convert_qt_terms": [
        ("a", ([(BAD_PART, [(0, 0, 1)])], "Schur", "Schur")),
        ("src", ([([2], [(0, 0, 1)])], "zzz", "Schur")),
        ("dst", ([([2], [(0, 0, 1)])], "Schur", "zzz")),
        ("dst powersum", ([([2], [(0, 0, 1)])], "Schur", "powersum")),
        ("dst p", ([([2], [(0, 0, 1)])], "s", "p")),
    ],
    # --- Schubert: the term list, the index, and the rank ---
    "schubert_multiply": [("a", (BAD_SCHUB, [([1], 1)]))],
    "schubert_multiply_variable": [
        ("a", (BAD_SCHUB, 1)),
        ("i zero", ([([2, 1], 1)], 0)),
        ("i past MAX_SUPPORT", ([([2, 1], 1)], 32)),
        ("i at u32 max", ([([2, 1], 1)], 4294967295)),
    ],
    "schubert_divided_difference": [
        ("a", (BAD_SCHUB, 1)),
        ("i zero", ([([2, 1], 1)], 0)),
        ("i at u32 max", ([([2, 1], 1)], 4294967295)),
    ],
    "schubert_divided_difference_perm": [
        ("a", (BAD_SCHUB, [1, 2])),
        ("w", ([([2, 1], 1)], [1, 1])),
    ],
    "schubert_expand": [("a", (BAD_SCHUB,))],
    "schubert_pairing": [
        ("a", (BAD_SCHUB, [([1], 1)], 2)),
        ("n past MAX_SUPPORT", ([([2, 1], 1)], [([2, 1], 1)], 33)),
        ("term outside S_n", ([([3, 2, 1], 1)], [([1], 1)], 2)),
    ],
    "schubert_scalar_product": [
        ("a", (BAD_SCHUB, [([1], 1)], 2)),
        ("n past MAX_SUPPORT", ([([2, 1], 1)], [([2, 1], 1)], 33)),
    ],
    "polynomial_to_schubert": [("exponent past MAX_SUPPORT", ([([32], 1)],))],
    "schubert_dimension": [("w", ([1, 1],))],
    "schubert_coefficient": [("u", ([1, 1], [1], [2, 1]))],
    "schubert_monomial_mass": [("u", ([1, 1], [1]))],
    "schubert_to_stanley_schur": [("w", ([1, 1],))],
    # --- the cache budget: a byte count or None, nothing else ---
    "set_cache_budget": [("negative", (-1,)), ("not a count", ("1 GB",))],
    # --- the Macdonald operator algebra: homogeneity ---
    "nabla": [("inhomogeneous", (BAD_QTSCHUR,))],
    "nabla_power": [("inhomogeneous", (BAD_QTSCHUR, 2))],
    "delta_ek": [("inhomogeneous", (1, BAD_QTSCHUR))],
    "delta_prime_ek": [("inhomogeneous", (1, BAD_QTSCHUR))],
    "theta_ek": [("inhomogeneous", (1, BAD_QTSCHUR))],
    "big_pi": [("inhomogeneous", (BAD_QTSCHUR,))],
    "delta_conjecture_side": [("side", (3, "zzz"))],
    # --- LLT: the level, the abacus, the tuple, the graph ---
    "llt_gtilde": [("k zero", ([2, 1], 0)), ("abacus", ([130], 1))],
    "llt_h": [("k zero", ([2, 1], 0)), ("abacus", ([130], 1))],
    "llt_h_tilde": [("k zero", ([2, 1], 0)), ("abacus", ([130], 1))],
    "llt_g_lt": [("k zero", ([2, 1], 0)), ("abacus", ([130], 1))],
    "llt_schur": [("k zero", ([2, 1], 0)), ("abacus", ([130], 1))],
    "llt_kl_column": [("k zero", ([2, 1], 0)), ("abacus", ([130], 1))],
    "llt_h_table": [("k zero", (2, 0))],
    "monomial_in_llt_h_table": [("k zero", (2, 0)), ("abacus", (130, 1))],
    "monomial_in_llt_h_tilde_table": [("k zero", (2, 0)), ("abacus", (130, 1))],
    "llt_gtilde_table": [("k zero", (2, 0))],
    "k_core_quotient": [("k zero", ([2, 1], 0)), ("lambda", (BAD_PART, 2))],
    "htilde_by_llt": [("mu", (BAD_PART,)), ("cells", ([1] * 65,))],
    "nabla_e_by_path": [("cells", (65,))],
    "llt_g": [
        ("shapes", ([BAD_PART],)),
        ("cells", ([[65]],)),
        ("skew", ([([1], [2])],)),
    ],
    "llt_min_inv": [
        ("shapes", ([BAD_PART],)),
        ("cells", ([[65]],)),
        ("skew", ([([1], [2])],)),
    ],
    "llt_fundamental": [("shapes", ([BAD_PART],)), ("cells", ([[65]],))],
    "llt_graph": [
        ("strict misoriented", (3, [], [(1, 0)])),
        ("edge out of range", (3, [(0, 9)], [])),
        ("edge both weak and strict", (3, [(0, 1)], [(0, 1)])),
    ],
    "chromatic_from_llt": [("edge out of range", (3, [(0, 9)], []))],
    "llt_e_expansion": [
        ("edge out of range", (3, [(0, 9)], [])),
        (
            "free edges past the mask",
            (9, [(i, j) for i in range(9) for j in range(i + 1, 9)], []),
        ),
    ],
}

# Functions with no malformed *typed* input: their arguments are scalars whose
# whole range is meaningful, so there is no precondition for a caller to
# violate. They are listed rather than omitted so the completeness check stays
# honest — an entry here is a claim, not a skip.
TOTAL = {
    "clear_caches": "no arguments",
    "cache_stats": "no arguments",
    "cache_budget": "no arguments",
    "partitions": "every u32 is a degree",
    "character_table": "every u32 is a degree",
    "kostka_table": "every u32 is a degree",
    "kostka_foulkes_table": "every u32 is a degree",
    "qt_kostka_table": "every u32 is a degree",
    "schur_in_macdonald_j": "every u32 is a degree",
    "hall_littlewood_table": "every u32 is a degree",
    "hall_littlewood_p_table": "every u32 is a degree",
    "jack_table": "every u32 is a degree",
    "stanley_table": "every u32 is a degree",
    "gj_connection_tables": "every u32 is a degree",
    "nabla_e": "every u32 is a degree",
    "delta_prime_e": "every (k, n) pair has a value; k >= n is zero",
}


def exported(mod):
    return {n for n in dir(mod) if not n.startswith("_") and callable(getattr(mod, n))}


def run(mod):
    """Call every case; return the failures.

    A failure is anything that is not a typed Python exception — a
    `PanicException`, or a call that returned normally when a precondition was
    violated. `PanicException` subclasses `BaseException` rather than
    `Exception`, so the catch here is deliberately as wide as it can be.
    """
    failures = []
    checked = 0
    for name, cases in sorted(CASES.items()):
        fn = getattr(mod, name, None)
        if fn is None:
            failures.append(f"{name}: not exported by the module")
            continue
        for label, args in cases:
            checked += 1
            try:
                fn(*args)
            except (ValueError, TypeError, OverflowError):
                continue
            except BaseException as e:  # noqa: BLE001 - the point is to catch panics
                kind = type(e).__name__
                failures.append(f"{name} [{label}]: {kind}, not a typed exception: {e}")
                continue
            failures.append(f"{name} [{label}]: returned normally, no exception raised")
    return checked, failures


def check_padding_still_works(mod):
    """A padded but valid partition is data, not an error.

    Sage hands over fixed-width lists, so trailing zeros are the shape this
    boundary promises to tolerate — the validation added for `[1, 3]` must not
    have cost that. `convert_indexed`'s identity conversion is the case that
    used to crash on it.
    """
    out = []
    checks = [
        ("convert_indexed", ([([2, 1, 0], 1)], "Schur", "Schur"), [(3, 1, 1)]),
        ("kostka_number", ([2, 1, 0], [1, 1, 1, 0]), 2),
        ("lr_coefficient", ([2, 1, 0], [2, 0], [1, 0]), 1),
    ]
    for name, args, want in checks:
        try:
            got = getattr(mod, name)(*args)
        except BaseException as e:  # noqa: BLE001
            out.append(f"{name}{args}: raised {type(e).__name__}: {e}")
            continue
        if got != want:
            out.append(f"{name}{args}: got {got}, want {want}")
    return out


def check_basis_codes_alias_names(mod):
    """A one-letter basis code and its full name select the same computation.

    Every entry point that takes a basis argument parses it through one table,
    and this is the check that the table has both spellings for all six bases
    on every such function — a code missing from one match arm would raise, and
    a code mapped to the wrong basis would answer differently. Both are visible
    only from here, since the convenience layer only ever sends codes.
    """
    out = []
    pairs = [
        ("s", "Schur"),
        ("h", "homogeneous"),
        ("e", "elementary"),
        ("p", "powersum"),
        ("m", "monomial"),
        ("f", "forgotten"),
    ]
    a = [([2, 1], 1), ([1, 1, 1], 2)]
    g = [([1], 1)]
    for code, name in pairs:
        calls = [
            ("to_power", (a, code), (a, name)),
            ("expand_alphabet", (a, code, 3), (a, name, 3)),
            ("skew_by", (a, g, code), (a, g, name)),
            ("convert_terms", (a, code, "s"), (a, name, "Schur")),
            ("convert_indexed", (a, code, "s"), (a, name, "Schur")),
        ]
        if code != "p":
            calls += [
                ("convert_terms", (a, "s", code), (a, "Schur", name)),
                ("convert_indexed", (a, "s", code), (a, "Schur", name)),
            ]
        for fn, by_code, by_name in calls:
            try:
                got, want = getattr(mod, fn)(*by_code), getattr(mod, fn)(*by_name)
            except BaseException as e:  # noqa: BLE001
                out.append(f"{fn} {code!r}/{name!r}: raised {type(e).__name__}: {e}")
                continue
            if got != want:
                out.append(f"{fn}: {code!r} gave {got}, {name!r} gave {want}")
    return out


def check_theorem_zeros_still_answer(mod):
    """A zero that is a theorem must stay a zero.

    The counterpart to `CASES`, and the reason it exists: the entry points
    that raise off-degree do so because their object has no referent there,
    not because a mismatch is suspicious. These others vanish by a theorem —
    `s_1·s_1` really has no `s_3` term, `s_λ` in `n` variables really is 0 when
    `ℓ(λ) > n` — and a caller sweeping a range depends on getting the value.
    Without this check, "validate more" would eventually eat them.
    """
    out = []
    checks = [
        ("lr_coefficient", ([3], [1], [1]), 0),
        ("lr_coefficient", ([3, 1], [2, 2], [1]), 0),
        ("kostka_number", ([2, 1], [2, 2]), 0),
        ("kostka_number", ([1, 1, 1], [3]), 0),
        ("semistandard_tableaux", ([2, 1], [2, 2]), []),
        ("skew_schur", ([2, 1], [3]), []),
        ("principal_specialization", ([1, 1, 1], 2), 0),
        ("evaluate_schur", ([([1, 1, 1], 1)], [2, 3]), 0),
        # Zero is the empty numerator, not a [0] one — the AFrac normal form.
        ("jack_structure_constant", ([2], [1], [2]), ([], [], 1, [])),
    ]
    for name, args, want in checks:
        try:
            got = getattr(mod, name)(*args)
        except BaseException as e:  # noqa: BLE001
            out.append(
                f"{name}{args}: raised {type(e).__name__} — this zero is a "
                f"theorem, not a malformed question: {e}"
            )
            continue
        if got != want:
            out.append(f"{name}{args}: got {got}, want {want}")
    return out


def check_surface_is_covered(mod):
    """Every exported callable is either exercised or declared total."""
    named = set(CASES) | set(TOTAL)
    missing = sorted(exported(mod) - named)
    stale = sorted(named - exported(mod))
    out = []
    if missing:
        out.append(
            "these pyfunctions have no malformed-input case and are not "
            f"declared total in TOTAL: {', '.join(missing)}"
        )
    if stale:
        out.append(f"these names are listed but not exported: {', '.join(stale)}")
    return out


def main():
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {sys.argv[0]} <path to built symfn extension>")
    mod = load(sys.argv[1])

    failures = check_surface_is_covered(mod)
    checked, bad = run(mod)
    failures += bad
    failures += check_padding_still_works(mod)
    failures += check_basis_codes_alias_names(mod)
    failures += check_theorem_zeros_still_answer(mod)

    if failures:
        print(f"FAIL: {len(failures)} problem(s) across {checked} malformed calls\n")
        for f in failures:
            print(f"  {f}")
        raise SystemExit(1)
    print(
        f"ok: {checked} malformed calls across {len(CASES)} pyfunctions, "
        f"every one a typed exception; {len(TOTAL)} declared total"
    )


if __name__ == "__main__":
    main()
