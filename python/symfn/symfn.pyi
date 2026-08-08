"""Type stubs for the symfn extension module.

This file is the **supported surface made machine-readable**: every function
`symfn` exports appears here exactly once, and a name absent from this file is
not part of the module (`docs/policies/python.md`, P10).
`scripts/check_python_stubs.py` fails if the two ever disagree.

The prose contract — how data crosses, what a coefficient can be, the
convention each family uses — lives in the module's own docstring and in each
function's, where `help(symfn.jack_p)` reaches it. Summaries here are the
first sentence of that text, for editors that read stubs rather than the
built module.

**Arguments and returns are deliberately asymmetric.** The boundary is strict
about what it promises and permissive about what it accepts, so a partition
comes back as `tuple[int, ...]` and is accepted as any `Sequence[int]` — a
caller may hand back what it received, or pass the list it already has. The
`*Arg` aliases below are the permissive halves.

A partition is weakly decreasing and positive, with trailing zeros tolerated,
except in the Schubert entry points, where the same shape is a **permutation
in one-line notation**. The aliases carry that distinction where a Rust type
alias already drew it; where it is inline, the parameter name and the
docstring say which.
"""

from typing import Sequence

__version__: str

# A symmetric function, as `(partition, coefficient)` pairs.
Element = list[tuple[tuple[int, ...], int]]
ElementArg = Sequence[tuple[Sequence[int], int]]
# A Schubert polynomial, as `(permutation, coefficient)` pairs.
SchubertElement = list[tuple[tuple[int, ...], int]]
SchubertElementArg = Sequence[tuple[Sequence[int], int]]
# Coefficients that need a denominator, as `(numerator, denominator)`.
RationalElement = list[tuple[tuple[int, ...], tuple[int, int]]]
# `(partition, [(q_exp, t_exp, coefficient), ...])` rows.
QtElement = list[tuple[tuple[int, ...], list[tuple[int, int, int]]]]
# A Macdonald element: monomial rows, then the factored denominator.
MacdonaldElement = list[
    tuple[tuple[int, ...], list[tuple[int, int, int]], list[tuple[int, int, int]]]
]
# A Jack element: monomial rows over ℚ(α), numerators then factored atoms.
JackElement = list[
    tuple[tuple[int, ...], list[int], list[tuple[int, int, int]], int]
]

def hall_littlewood(la: Sequence[int]) -> list[tuple[tuple[int, ...], list[tuple[int, int]]]]:
    """`Q'_λ(x; t) = Σ_μ K_{μλ}(t) s_μ`, as `[(mu, [(t_exponent,
    coefficient), ...])]`.
    """
    ...

def hall_littlewood_table(n: int) -> list[tuple[tuple[int, ...], list[tuple[tuple[int, ...], list[tuple[int, int]]]]]]:
    """Every `Q'_λ` for `λ ⊢ n`, sharing the recursion's suffixes across
    the degree.
    """
    ...

def kostka_foulkes(la: Sequence[int], mu: Sequence[int]) -> list[tuple[int, int]]:
    """`K_{λμ}(t)` as `[(t_exponent, coefficient), ...]`.
    """
    ...

def kostka_foulkes_column(mu: Sequence[int]) -> list[tuple[tuple[int, ...], list[tuple[int, int]]]]:
    """Every `K_{λμ}(t)` for a fixed μ, as `[(lambda, [(t_exponent,
    coefficient)])]`.
    """
    ...

def kostka_foulkes_table(n: int) -> list[list[list[tuple[int, int]]]]:
    """The whole `K_{λμ}(t)` matrix for degree `n`, indexed as
    `partitions(n)` is.
    """
    ...

def hall_littlewood_p(la: Sequence[int]) -> list[tuple[tuple[int, ...], list[tuple[int, int]]]]:
    """`P_λ(x; t)` in the Schur basis — the other Hall–Littlewood
    normalisation.
    """
    ...

def hall_littlewood_p_table(n: int) -> list[tuple[tuple[int, ...], list[tuple[tuple[int, ...], list[tuple[int, int]]]]]]:
    """Every `P_λ` for `λ ⊢ n`, from one inversion of the Kostka–Foulkes
    matrix.
    """
    ...

def macdonald_p(la: Sequence[int]) -> MacdonaldElement:
    """Macdonald `P_λ(x; q, t)` in the monomial basis.
    """
    ...

def macdonald_q(la: Sequence[int]) -> MacdonaldElement:
    """Macdonald `Q_λ = b_λ · P_λ`. Escalates, as [`macdonald_p`] does; the
    `i128` wall underneath is n = 26 at λ = (n).
    """
    ...

def macdonald_j(la: Sequence[int]) -> MacdonaldElement:
    """Macdonald `J_λ = c_λ · P_λ`, the integral form — every coefficient
    is a polynomial, so the denominator list comes back empty.
    Escalates, as [`macdonald_p`] does; the `i128` wall underneath is n
    = 26 at λ = (n).
    """
    ...

def schur_in_macdonald_j(n: int) -> list[tuple[tuple[int, ...], MacdonaldElement]]:
    """The Schur functions of degree `n` in the Macdonald `J` basis, as
    `[(lambda, [(mu, numerator, denominator), ...]), ...]` — the
    **inverse** of the `J → s` transition, in the cell encoding
    [`MacTerms`] already carries.
    """
    ...

def jack_p(la: Sequence[int]) -> JackElement:
    """Jack `P_λ(x; α)` in the monomial basis: monic in `m_λ`, dominance-
    triangular.
    """
    ...

def jack_q(la: Sequence[int]) -> JackElement:
    """Jack `Q_λ = (H_λ/H'_λ)·P_λ`, the basis dual to `P` under `⟨·,·⟩_α`.
    """
    ...

def jack_j(la: Sequence[int]) -> JackElement:
    """Jack `J_λ = H_λ·P_λ`, the integral form.
    """
    ...

def jack_table(n: int) -> list[tuple[tuple[int, ...], JackElement]]:
    """Every `P_λ` of degree `n` — the unit of work Sage has no entry point
    for, and the one `docs/record/jack.md` measures the walls in.
    """
    ...

def jack_j_powersum(la: Sequence[int]) -> JackElement:
    """`J_λ` in the **power-sum** basis — the Jack character table, and the
    unit the Goulden–Jackson pipeline consumes.
    """
    ...

def jack_norm_j(la: Sequence[int]) -> list[tuple[int, int, int]]:
    """`⟨J_λ, J_λ⟩_α = H_λ·H'_λ`, returned **factored** as `[(u, v,
    mult)]`.
    """
    ...

def jack_scalar(f: Sequence[tuple[Sequence[int], Sequence[int]]], g: Sequence[tuple[Sequence[int], Sequence[int]]]) -> tuple[list[int], list[tuple[int, int, int]], int]:
    """`⟨f, g⟩_α` for two monomial-basis elements whose coefficients are
    **integer** polynomials in α, given densely: `[(partition, [c0, c1,
    …])]`.
    """
    ...

def jack_structure_constant(la: Sequence[int], mu: Sequence[int], nu: Sequence[int]) -> tuple[list[int], list[tuple[int, int, int]], int]:
    """`⟨J_λ J_μ, J_ν⟩_α` — **Stanley's object**, whose membership in
    `ℕ[α]` is his 1989 conjecture and still open.
    """
    ...

def stanley_table(k: int) -> list[tuple[list[int], list[int], list[int], list[int], list[tuple[int, int, int]], int]]:
    """Stanley's **whole table**: every `⟨J_λ J_μ, J_ν⟩_α` with `|λ| = |μ|
    = k`, as `(lambda, mu, nu, numerator, denominator atoms, scalar)`.
    """
    ...

def zonal(la: Sequence[int], integral_form: bool) -> list[tuple[tuple[int, ...], int, int]]:
    """The zonal polynomial, in **both** circulating normalizations, as
    exact `(numerator, denominator)` pairs.
    """
    ...

def gj_connection_tables(n: int) -> tuple[list[tuple[list[int], list[int], list[int], list[int], int]], list[tuple[list[int], list[int], list[int], list[int], int]]]:
    """The Goulden–Jackson connection tables `c^λ_{μν}(b)` and
    `h^λ_{μν}(b)` at degree `n`, as `(lambda, mu, nu, [b-coefficients],
    denominator)`.
    """
    ...

def class_algebra_coefficient(la: Sequence[int], mu: Sequence[int], nu: Sequence[int]) -> int:
    """`a^λ_{μν}`, the class-algebra connection coefficient of `S_n`, from
    characters alone.
    """
    ...

def qt_kostka(la: Sequence[int], mu: Sequence[int]) -> list[tuple[int, int, int]]:
    """The (q,t)-Kostka polynomial `K_{λμ}(q,t)`, from `J_μ = Σ_λ K_{λμ}
    S_λ(x;t)`.
    """
    ...

def qt_kostka_column(mu: Sequence[int]) -> list[tuple[tuple[int, ...], list[tuple[int, int, int]]]]:
    """Every `K_{λμ}(q,t)` for a fixed μ — one `J_μ`, which is what a
    single [`qt_kostka`] costs anyway.
    """
    ...

def qt_kostka_table(n: int) -> list[list[list[tuple[int, int, int]]]]:
    """The whole `K_{λμ}(q,t)` matrix for degree `n`, indexed as
    `partitions(n)` is — the same orientation as [`kostka_table`] and
    [`kostka_foulkes_table`], of which this is the two-variable
    analogue. `q = 0` recovers the latter.
    """
    ...

def macdonald_ht(mu: Sequence[int]) -> list[tuple[tuple[int, ...], list[tuple[int, int, int]]]]:
    """The modified Macdonald polynomial `H̃_μ(x;q,t)` in the **Schur**
    basis, as `[(lambda, [(q_exp, t_exp, coeff), ...])]`.
    """
    ...

def nabla_e(n: int) -> QtElement:
    """`∇e_n` in the Schur basis — the shuffle theorem's object.
    """
    ...

def delta_prime_e(k: int, n: int) -> QtElement:
    """`Δ'_{e_k} e_n` in the Schur basis — the Delta conjecture's object.
    """
    ...

def nabla(f: Sequence[tuple[Sequence[int], Sequence[tuple[int, int, int]]]]) -> QtElement:
    """`∇F` for an arbitrary homogeneous `F`, given in the Schur basis.
    """
    ...

def nabla_power(f: Sequence[tuple[Sequence[int], Sequence[tuple[int, int, int]]]], r: int) -> QtElement:
    """`∇^r F`, sharing one change of basis across the powers — the object
    Qiu–Zhang's 2026 theorem is about.
    """
    ...

def delta_ek(k: int, f: Sequence[tuple[Sequence[int], Sequence[tuple[int, int, int]]]]) -> QtElement:
    """`Δ_{e_k} F`, with eigenvalue `e_k[B_μ]`.
    """
    ...

def delta_prime_ek(k: int, f: Sequence[tuple[Sequence[int], Sequence[tuple[int, int, int]]]]) -> QtElement:
    """`Δ'_{e_k} F`, with eigenvalue `e_k[B_μ − 1]`.
    """
    ...

def theta_ek(k: int, f: Sequence[tuple[Sequence[int], Sequence[tuple[int, int, int]]]]) -> QtElement:
    """`Θ_{e_k} F`, which raises the degree by `k`.
    """
    ...

def big_pi(f: Sequence[tuple[Sequence[int], Sequence[tuple[int, int, int]]]]) -> QtElement:
    """`ΠF`, with eigenvalue `Π_μ`.
    """
    ...

def delta_conjecture_side(n: int, side: str) -> list[QtElement]:
    """The combinatorial side of the Delta conjecture, in the **monomial**
    basis, for every `k` at once — entry `k` of the returned list.
    """
    ...

def llt_gtilde(la: Sequence[int], k: int) -> QtElement:
    """`G̃^(k)_λ(x;q)`, the **cospin** ribbon generating function of [LLT]
    (26), in the monomial basis.
    """
    ...

def llt_h(mu: Sequence[int], k: int) -> QtElement:
    """`H^(k)_μ(x;q) = Σ_R q^{s(R)} x^{w(R)}`, the **spin** family of [LLT]
    (28).
    """
    ...

def llt_h_tilde(mu: Sequence[int], k: int) -> QtElement:
    """`H̃^(k)_μ = G̃^(k)_{kμ}` ([LLT] (27)) — Sage's
    `llt(k).hcospin()[μ]`.
    """
    ...

def llt_g_lt(la: Sequence[int], k: int) -> QtElement:
    """`Σ_R q^{2s(R)} x^{w(R)}`, the spin-generating grading of [LT] (43).
    """
    ...

def llt_h_table(n: int, k: int) -> list[tuple[tuple[int, ...], QtElement]]:
    """`H^(k)_μ` for **every** μ ⊢ n — the whole degree, which is the unit
    `docs/record/llt.md` measures the walls in.
    """
    ...

def llt_gtilde_table(n: int, k: int) -> list[tuple[tuple[int, ...], QtElement]]:
    """`G̃^(k)_λ` for **every** λ ⊢ k·n with empty k-core, from a single
    walk.
    """
    ...

def llt_schur(la: Sequence[int], k: int) -> QtElement:
    """`G̃^(k)_λ` in the **Schur** basis.
    """
    ...

def llt_g(shapes: Sequence[Sequence[int]], offsets: Sequence[int] | None = None) -> QtElement:
    """`G_ν(x;q)` for a tuple of shapes, in the monomial basis and the
    **raw** inv grading.
    """
    ...

def llt_min_inv(shapes: Sequence[Sequence[int]], offsets: Sequence[int] | None = None) -> int:
    """`min_T inv(T)` over the semistandard fillings of a tuple — the
    forced `q`-floor that [`llt_g`] does not divide out.
    """
    ...

def llt_fundamental(shapes: Sequence[Sequence[int]], offsets: Sequence[int] | None = None) -> list[tuple[tuple[int, ...], list[tuple[int, int, int]]]]:
    """The **fundamental quasisymmetric** expansion of `G_ν`, as
    `[(composition, [(q_exp, t_exp, coeff), ...]), ...]`.
    """
    ...

def llt_kl_column(la: Sequence[int], k: int) -> list[tuple[tuple[int, ...], list[tuple[int, int, int]]]]:
    """One **column** of the Schur-expansion table: `c^λ_μ` for every shape
    μ ⊢ k|λ|, in the [KMS] variable `v`, as `[(mu, [(v_exp, 0, coeff),
    ...]), ...]`.
    """
    ...

def llt_graph(n: int, weak: Sequence[tuple[int, int]], strict: Sequence[tuple[int, int]]) -> QtElement:
    """`G_Γ(x;q) = Σ_κ q^{asc(κ)} x^κ` over the colorings of a decorated
    graph.
    """
    ...

def chromatic_from_llt(n: int, weak: Sequence[tuple[int, int]], strict: Sequence[tuple[int, int]]) -> QtElement:
    """The Shareshian–Wachs chromatic quasisymmetric function `X_Γ(x;q)` of
    Γ, from its LLT polynomial by the `(q−1)`-plethysm of [CM] Prop 3.5.
    """
    ...

def llt_e_expansion(n: int, weak: Sequence[tuple[int, int]], strict: Sequence[tuple[int, int]]) -> list[tuple[tuple[int, ...], list[tuple[int, int, int]]]]:
    """The [AS] **e-expansion** of `Ĝ_Γ(x; q+1)`: `Σ_θ q^{asc(θ)} e_{λ(θ)}`
    over orientations of the free edges, as `[(partition, poly), ...]`.
    """
    ...

def htilde_by_llt(mu: Sequence[int]) -> QtElement:
    """`H̃_μ(x;q,t) = Σ_D q^{−a(D)} t^{maj(D)} G_{ν(μ,D)}(x;q)` — the [HHL]
    decomposition, in the monomial basis.
    """
    ...

def nabla_e_by_path(n: int) -> list[tuple[tuple[int, ...], QtElement]]:
    """`∇e_n = Σ_D t^{area(D)} G_D(x;q)`, as `[(area_sequence, G_D), ...]`.
    """
    ...

def k_core_quotient(la: Sequence[int], k: int) -> tuple[tuple[int, ...], list[tuple[int, ...]]]:
    """The k-core and k-quotient of λ, as `(core, [component, ...])`.
    """
    ...

def clear_caches() -> None:
    """Drop every memo cache.
    """
    ...

def schur_multiply(a: ElementArg, b: ElementArg) -> Element:
    """Multiply two Schur-basis elements (Littlewood–Richardson).
    """
    ...

def st_multiply(a: ElementArg, b: ElementArg) -> Element:
    """Multiply two `st`-basis elements — the Orellana–Zabrocki irreducible
    character basis `s̃`, whose structure constants **are** the reduced
    (stable) Kronecker coefficients.
    """
    ...

def reduced_kronecker_product(la: Sequence[int], mu: Sequence[int]) -> Element:
    """The reduced Kronecker product `s̃_λ · s̃_μ`, as one column.
    """
    ...

def reduced_kronecker(la: Sequence[int], mu: Sequence[int], nu: Sequence[int]) -> int:
    """One reduced Kronecker coefficient `ḡ^ν_{λμ}`.
    """
    ...

def schur_to_st(a: ElementArg) -> Element:
    """`s → st`: rewrite a Schur-basis element in the character basis.
    """
    ...

def st_to_schur(a: ElementArg) -> Element:
    """`st → s`: rewrite a character-basis element in the Schur basis.
    """
    ...

def ht_multiply(a: ElementArg, b: ElementArg) -> Element | None:
    """Multiply two `ht`-basis elements — the Orellana–Zabrocki **induced
    trivial** character basis `h̃`, by the double-coset matrix rule.
    """
    ...

def schur_to_ht(a: ElementArg) -> Element:
    """`s → ht`: rewrite a Schur-basis element in the induced trivial
    character basis.
    """
    ...

def ht_to_schur(a: ElementArg) -> Element:
    """`ht → s`: rewrite an induced trivial character element in the Schur
    basis.
    """
    ...

def lr_coefficient(la: Sequence[int], mu: Sequence[int], nu: Sequence[int]) -> int:
    """A single Littlewood–Richardson coefficient c^λ_{μν}.
    """
    ...

def schur_to_homogeneous(a: ElementArg) -> Element: ...

def schur_to_elementary(a: ElementArg) -> Element: ...

def schur_to_monomial(a: ElementArg) -> Element: ...

def schur_to_forgotten(a: ElementArg) -> Element: ...

def to_power(a: ElementArg, src: str) -> RationalElement:
    """A conversion into the power-sum basis, from any basis
    [`convert_indexed`] names.
    """
    ...

def homogeneous_to_schur(a: ElementArg) -> Element: ...

def elementary_to_schur(a: ElementArg) -> Element: ...

def monomial_to_schur(a: ElementArg) -> Element: ...

def forgotten_to_schur(a: ElementArg) -> Element: ...

def power_to_schur(a: ElementArg) -> Element: ...

def plethysm(f: ElementArg, g: ElementArg) -> Element:
    """Plethysm f[g] of two Schur-basis elements.
    """
    ...

def kostka_number(la: Sequence[int], mu: Sequence[int]) -> int:
    """Kostka number K_{λμ}.
    """
    ...

def evaluate_schur(a: ElementArg, xs: Sequence[int]) -> int:
    """Evaluate a Schur-basis element at the alphabet `xs`.
    """
    ...

def expand_alphabet(a: ElementArg, src: str, n: int) -> list[tuple[tuple[int, ...], int]]:
    """The expansion of an element in `n` variables, as `[(exponent vector,
    coefficient), ...]` with each vector of length `n`.
    """
    ...

def monomial_multiply(a: ElementArg, b: ElementArg) -> Element:
    """Multiply two monomial-basis elements.
    """
    ...

def semistandard_tableaux(la: Sequence[int], mu: Sequence[int]) -> list[list[list[int]]]:
    """The semistandard Young tableaux of shape λ and weight μ, as lists of
    rows.
    """
    ...

def dimension(la: Sequence[int]) -> int | None:
    """f^λ — the number of standard Young tableaux of shape λ, i.e. the
    dimension of the irreducible S_{|λ|} representation. `None` past
    `u128`.
    """
    ...

def principal_specialization(la: Sequence[int], n: int) -> int | None:
    """s_λ(1^n), the dimension of the GL_n irreducible. `None` on overflow.
    """
    ...

def principal_specialization_q(la: Sequence[int], n: int) -> list[int]:
    """s_λ(1, q, …, q^{n−1}) as a coefficient list in q, lowest degree
    first.
    """
    ...

def character_value(la: Sequence[int], mu: Sequence[int]) -> int:
    """Symmetric-group character χ^λ(μ).
    """
    ...

def internal_product(a: ElementArg, b: ElementArg) -> Element:
    """The internal (Kronecker) product of two Schur-basis elements.
    """
    ...

def partitions(n: int) -> list[tuple[int, ...]]:
    """The partitions of `n`, in the order the table functions below index
    by.
    """
    ...

def kronecker_coefficient(la: Sequence[int], mu: Sequence[int], nu: Sequence[int]) -> int:
    """A single Kronecker coefficient g^ν_{λμ}, computed **without forming
    the product**.
    """
    ...

def convert_indexed(a: ElementArg, src: str, dst: str) -> list[tuple[int, int, int]]:
    """A conversion whose output partitions are returned as **indices**
    rather than as lists: `[(degree, index, coefficient), ...]`, where
    `index` is into [`partitions`] of that degree.
    """
    ...

def convert_terms(a: ElementArg, src: str, dst: str) -> Element:
    """The same conversion as [`convert_indexed`], with each output
    partition spelled out rather than given as an index.
    """
    ...

def character_table(n: int) -> list[list[int]]:
    """The full character table of S_n: `table[i][j]` = χ^{λⁱ}(λʲ).
    """
    ...

def kostka_table(n: int) -> list[list[int]]:
    """The full Kostka table of degree `n`: `table[i][j]` = K_{λⁱ λʲ}.
    """
    ...

def omega(a: ElementArg) -> Element:
    """The ω involution on a Schur-basis element.
    """
    ...

def hall_inner_product(a: ElementArg, b: ElementArg) -> int:
    """The Hall inner product of two Schur-basis elements.
    """
    ...

def skew_schur(la: Sequence[int], mu: Sequence[int]) -> Element:
    """The skew Schur function s_{λ/μ}.
    """
    ...

def skew_by(f: ElementArg, g: ElementArg, basis: str = "s") -> Element:
    """Skew a Schur-basis element by `g`, given in `basis` — the adjoint of
    multiplication by g under the Hall inner product.
    """
    ...

def coproduct(a: ElementArg) -> list[tuple[tuple[tuple[int, ...], tuple[int, ...]], int]]:
    """The coproduct Δ, as `[((mu, nu), coefficient), ...]`.
    """
    ...

def antipode(a: ElementArg) -> Element:
    """The antipode S.
    """
    ...

def schubert_multiply(a: SchubertElementArg, b: SchubertElementArg) -> SchubertElement:
    """Multiply two Schubert polynomials.
    """
    ...

def schubert_multiply_variable(a: SchubertElementArg, i: int) -> SchubertElement:
    """`x_i · f`, the signed Monk rule. **1-based**, unlike Symmetrica's
    `mult_schubert_variable`, which is 0-based while its own
    `divdiff_schubert` is 1-based. One convention, stated.
    """
    ...

def schubert_divided_difference(a: SchubertElementArg, i: int) -> SchubertElement:
    """`∂_i f` on the Schubert basis, 1-based.
    """
    ...

def schubert_divided_difference_perm(a: SchubertElementArg, w: Sequence[int]) -> SchubertElement:
    """`∂_w f`, composing along a reduced word of `w`.
    """
    ...

def schubert_expand(a: SchubertElementArg) -> Element:
    """Expand into monomials: `(exponent vector, coefficient)` pairs.
    """
    ...

def polynomial_to_schubert(terms: ElementArg) -> SchubertElement:
    """Write a polynomial in the Schubert basis (the greedy triangular
    peel).
    """
    ...

def schubert_pairing(a: SchubertElementArg, b: SchubertElementArg, n: int) -> int:
    """The Poincaré pairing on `H*(Fl(n))`.
    """
    ...

def schubert_scalar_product(
    a: SchubertElementArg, b: SchubertElementArg, n: int
) -> SchubertElement:
    """`∂_{w0(n)}(a·b)`, Symmetrica's `scalarproduct_schubert` — what Sage
    exposes as `SchubertPolynomial.scalar_product`.
    """
    ...

def schubert_dimension(w: Sequence[int]) -> int:
    """`S_w(1,…,1)`: the number of pipe dreams, i.e. the size
    `schubert_expand` would produce. Cheap — it never builds the
    expansion.
    """
    ...

def schubert_coefficient(u: Sequence[int], v: Sequence[int], w: Sequence[int]) -> int:
    """A single structure constant `c^w_{uv}`, **without building the
    product**.
    """
    ...

def schubert_monomial_mass(u: Sequence[int], v: Sequence[int]) -> int:
    """The product's total monomial mass `S_u(1,…,1)·S_v(1,…,1)`, a count.
    """
    ...

def schubert_to_stanley_schur(w: Sequence[int]) -> Element:
    """The **Stanley symmetric function** `F_w` in the Schur basis.
    """
    ...

