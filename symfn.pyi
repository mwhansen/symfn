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

`list[int]` means a **partition**, weakly decreasing and positive with
trailing zeros tolerated, except in the Schubert entry points, where it is a
**permutation in one-line notation**. The two are structurally identical and
the aliases below carry the distinction where a Rust type alias already drew
it; where it is inline, the parameter name and the docstring say which.
"""

from typing import TypeAlias

__version__: str

# A symmetric function, as `(partition, coefficient)` pairs.
Element: TypeAlias = list[tuple[list[int], int]]
# A Schubert polynomial, as `(permutation, coefficient)` pairs.
SchubertElement: TypeAlias = list[tuple[list[int], int]]
# Coefficients that need a denominator, as `(numerator, denominator)`.
RationalElement: TypeAlias = list[tuple[list[int], tuple[int, int]]]
# `(partition, [(q_exp, t_exp, coefficient), ...])` rows.
QtElement: TypeAlias = list[tuple[list[int], list[tuple[int, int, int]]]]
# A Macdonald element: monomial rows, then the factored denominator.
MacdonaldElement: TypeAlias = list[
    tuple[list[int], list[tuple[int, int, int]], list[tuple[int, int, int]]]
]
# A Jack element: monomial rows over ℚ(α), numerators then factored atoms.
JackElement: TypeAlias = list[
    tuple[list[int], list[int], list[tuple[int, int, int]], int]
]

def hall_littlewood(la: list[int]) -> list[tuple[list[int], list[tuple[int, int]]]]:
    """`Q'_λ(x; t) = Σ_μ K_{μλ}(t) s_μ`, as `[(mu, [(t_exponent,
    coefficient), ...])]`.
    """
    ...

def hall_littlewood_table(n: int) -> list[tuple[list[int], list[tuple[list[int], list[tuple[int, int]]]]]]:
    """Every `Q'_λ` for `λ ⊢ n`, sharing the recursion's suffixes across
    the degree.
    """
    ...

def kostka_foulkes(la: list[int], mu: list[int]) -> list[tuple[int, int]]:
    """`K_{λμ}(t)` as `[(t_exponent, coefficient), ...]`.
    """
    ...

def kostka_foulkes_column(mu: list[int]) -> list[tuple[list[int], list[tuple[int, int]]]]:
    """Every `K_{λμ}(t)` for a fixed μ, as `[(lambda, [(t_exponent,
    coefficient)])]`.
    """
    ...

def kostka_foulkes_table(n: int) -> list[list[list[tuple[int, int]]]]:
    """The whole `K_{λμ}(t)` matrix for degree `n`, indexed as
    `partitions(n)` is.
    """
    ...

def hall_littlewood_p(la: list[int]) -> list[tuple[list[int], list[tuple[int, int]]]]:
    """`P_λ(x; t)` in the Schur basis — the other Hall–Littlewood
    normalisation.
    """
    ...

def hall_littlewood_p_table(n: int) -> list[tuple[list[int], list[tuple[list[int], list[tuple[int, int]]]]]]:
    """Every `P_λ` for `λ ⊢ n`, from one inversion of the Kostka–Foulkes
    matrix.
    """
    ...

def macdonald_p(la: list[int]) -> MacdonaldElement:
    """Macdonald `P_λ(x; q, t)` in the monomial basis.
    """
    ...

def macdonald_q(la: list[int]) -> MacdonaldElement:
    """Macdonald `Q_λ = b_λ · P_λ`. Escalates, as [`macdonald_p`] does; the
    `i128` wall underneath is n = 26 at λ = (n).
    """
    ...

def macdonald_j(la: list[int]) -> MacdonaldElement:
    """Macdonald `J_λ = c_λ · P_λ`, the integral form — every coefficient
    is a polynomial, so the denominator list comes back empty.
    Escalates, as [`macdonald_p`] does; the `i128` wall underneath is n
    = 26 at λ = (n).
    """
    ...

def jack_p(la: list[int]) -> JackElement:
    """Jack `P_λ(x; α)` in the monomial basis: monic in `m_λ`, dominance-
    triangular.
    """
    ...

def jack_q(la: list[int]) -> JackElement:
    """Jack `Q_λ = (H_λ/H'_λ)·P_λ`, the basis dual to `P` under `⟨·,·⟩_α`.
    """
    ...

def jack_j(la: list[int]) -> JackElement:
    """Jack `J_λ = H_λ·P_λ`, the integral form.
    """
    ...

def jack_table(n: int) -> list[tuple[list[int], JackElement]]:
    """Every `P_λ` of degree `n` — the unit of work Sage has no entry point
    for, and the one `docs/record/jack.md` measures the walls in.
    """
    ...

def jack_j_powersum(la: list[int]) -> JackElement:
    """`J_λ` in the **power-sum** basis — the Jack character table, and the
    unit the Goulden–Jackson pipeline consumes.
    """
    ...

def jack_norm_j(la: list[int]) -> list[tuple[int, int, int]]:
    """`⟨J_λ, J_λ⟩_α = H_λ·H'_λ`, returned **factored** as `[(u, v,
    mult)]`.
    """
    ...

def jack_scalar(f: list[tuple[list[int], list[int]]], g: list[tuple[list[int], list[int]]]) -> tuple[list[int], list[tuple[int, int, int]], int]:
    """`⟨f, g⟩_α` for two monomial-basis elements whose coefficients are
    **integer** polynomials in α, given densely: `[(partition, [c0, c1,
    …])]`.
    """
    ...

def jack_structure_constant(la: list[int], mu: list[int], nu: list[int]) -> tuple[list[int], list[tuple[int, int, int]], int]:
    """`⟨J_λ J_μ, J_ν⟩_α` — **Stanley's object**, whose membership in
    `ℕ[α]` is his 1989 conjecture and still open.
    """
    ...

def stanley_table(k: int) -> list[tuple[list[int], list[int], list[int], list[int], list[tuple[int, int, int]], int]]:
    """Stanley's **whole table**: every `⟨J_λ J_μ, J_ν⟩_α` with `|λ| = |μ|
    = k`, as `(lambda, mu, nu, numerator, denominator atoms, scalar)`.
    """
    ...

def zonal(la: list[int], integral_form: bool) -> list[tuple[list[int], int, int]]:
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

def class_algebra_coefficient(la: list[int], mu: list[int], nu: list[int]) -> int:
    """`a^λ_{μν}`, the class-algebra connection coefficient of `S_n`, from
    characters alone.
    """
    ...

def qt_kostka(la: list[int], mu: list[int]) -> list[tuple[int, int, int]]:
    """The (q,t)-Kostka polynomial `K_{λμ}(q,t)`, from `J_μ = Σ_λ K_{λμ}
    S_λ(x;t)`.
    """
    ...

def qt_kostka_column(mu: list[int]) -> list[tuple[list[int], list[tuple[int, int, int]]]]:
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

def macdonald_ht(mu: list[int]) -> list[tuple[list[int], list[tuple[int, int, int]]]]:
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

def nabla(f: QtElement) -> QtElement:
    """`∇F` for an arbitrary homogeneous `F`, given in the Schur basis.
    """
    ...

def nabla_power(f: QtElement, r: int) -> QtElement:
    """`∇^r F`, sharing one change of basis across the powers — the object
    Qiu–Zhang's 2026 theorem is about.
    """
    ...

def delta_ek(k: int, f: QtElement) -> QtElement:
    """`Δ_{e_k} F`, with eigenvalue `e_k[B_μ]`.
    """
    ...

def delta_prime_ek(k: int, f: QtElement) -> QtElement:
    """`Δ'_{e_k} F`, with eigenvalue `e_k[B_μ − 1]`.
    """
    ...

def theta_ek(k: int, f: QtElement) -> QtElement:
    """`Θ_{e_k} F`, which raises the degree by `k`.
    """
    ...

def big_pi(f: QtElement) -> QtElement:
    """`ΠF`, with eigenvalue `Π_μ`.
    """
    ...

def delta_conjecture_side(n: int, side: str) -> list[QtElement]:
    """The combinatorial side of the Delta conjecture, in the **monomial**
    basis, for every `k` at once — entry `k` of the returned list.
    """
    ...

def llt_gtilde(la: list[int], k: int) -> QtElement:
    """`G̃^(k)_λ(x;q)`, the **cospin** ribbon generating function of [LLT]
    (26), in the monomial basis.
    """
    ...

def llt_h(mu: list[int], k: int) -> QtElement:
    """`H^(k)_μ(x;q) = Σ_R q^{s(R)} x^{w(R)}`, the **spin** family of [LLT]
    (28).
    """
    ...

def llt_h_tilde(mu: list[int], k: int) -> QtElement:
    """`H̃^(k)_μ = G̃^(k)_{kμ}` ([LLT] (27)) — Sage's
    `llt(k).hcospin()[μ]`.
    """
    ...

def llt_g_lt(la: list[int], k: int) -> QtElement:
    """`Σ_R q^{2s(R)} x^{w(R)}`, the spin-generating grading of [LT] (43).
    """
    ...

def llt_h_table(n: int, k: int) -> list[tuple[list[int], QtElement]]:
    """`H^(k)_μ` for **every** μ ⊢ n — the whole degree, which is the unit
    `docs/record/llt.md` measures the walls in.
    """
    ...

def llt_gtilde_table(n: int, k: int) -> list[tuple[list[int], QtElement]]:
    """`G̃^(k)_λ` for **every** λ ⊢ k·n with empty k-core, from a single
    walk.
    """
    ...

def llt_schur(la: list[int], k: int) -> QtElement:
    """`G̃^(k)_λ` in the **Schur** basis.
    """
    ...

def llt_g(shapes: list[list[int]], offsets: list[int] | None = None) -> QtElement:
    """`G_ν(x;q)` for a tuple of shapes, in the monomial basis and the
    **raw** inv grading.
    """
    ...

def llt_min_inv(shapes: list[list[int]], offsets: list[int] | None = None) -> int:
    """`min_T inv(T)` over the semistandard fillings of a tuple — the
    forced `q`-floor that [`llt_g`] does not divide out.
    """
    ...

def llt_fundamental(shapes: list[list[int]], offsets: list[int] | None = None) -> list[tuple[list[int], list[tuple[int, int, int]]]]:
    """The **fundamental quasisymmetric** expansion of `G_ν`, as
    `[(composition, [(q_exp, t_exp, coeff), ...]), ...]`.
    """
    ...

def llt_kl_column(la: list[int], k: int) -> list[tuple[list[int], list[tuple[int, int, int]]]]:
    """One **column** of the Schur-expansion table: `c^λ_μ` for every shape
    μ ⊢ k|λ|, in the [KMS] variable `v`, as `[(mu, [(v_exp, 0, coeff),
    ...]), ...]`.
    """
    ...

def llt_graph(n: int, weak: list[tuple[int, int]], strict: list[tuple[int, int]]) -> QtElement:
    """`G_Γ(x;q) = Σ_κ q^{asc(κ)} x^κ` over the colorings of a decorated
    graph.
    """
    ...

def chromatic_from_llt(n: int, weak: list[tuple[int, int]], strict: list[tuple[int, int]]) -> QtElement:
    """The Shareshian–Wachs chromatic quasisymmetric function `X_Γ(x;q)` of
    Γ, from its LLT polynomial by the `(q−1)`-plethysm of [CM] Prop 3.5.
    """
    ...

def llt_e_expansion(n: int, weak: list[tuple[int, int]], strict: list[tuple[int, int]]) -> list[tuple[list[int], list[tuple[int, int, int]]]]:
    """The [AS] **e-expansion** of `Ĝ_Γ(x; q+1)`: `Σ_θ q^{asc(θ)} e_{λ(θ)}`
    over orientations of the free edges, as `[(partition, poly), ...]`.
    """
    ...

def htilde_by_llt(mu: list[int]) -> QtElement:
    """`H̃_μ(x;q,t) = Σ_D q^{−a(D)} t^{maj(D)} G_{ν(μ,D)}(x;q)` — the [HHL]
    decomposition, in the monomial basis.
    """
    ...

def nabla_e_by_path(n: int) -> list[tuple[list[int], QtElement]]:
    """`∇e_n = Σ_D t^{area(D)} G_D(x;q)`, as `[(area_sequence, G_D), ...]`.
    """
    ...

def k_core_quotient(la: list[int], k: int) -> tuple[list[int], list[list[int]]]:
    """The k-core and k-quotient of λ, as `(core, [component, ...])`.
    """
    ...

def clear_caches() -> None:
    """Drop every memo cache.
    """
    ...

def schur_multiply(a: Element, b: Element) -> Element:
    """Multiply two Schur-basis elements (Littlewood–Richardson).
    """
    ...

def lr_coefficient(la: list[int], mu: list[int], nu: list[int]) -> int:
    """A single Littlewood–Richardson coefficient c^λ_{μν}.
    """
    ...

def schur_to_homogeneous(a: Element) -> Element: ...

def schur_to_elementary(a: Element) -> Element: ...

def schur_to_monomial(a: Element) -> Element: ...

def schur_to_forgotten(a: Element) -> Element: ...

def schur_to_power(a: Element) -> RationalElement:
    """s → p. Coefficients are rational, returned as `(numerator,
    denominator)`.
    """
    ...

def homogeneous_to_schur(a: Element) -> Element: ...

def elementary_to_schur(a: Element) -> Element: ...

def monomial_to_schur(a: Element) -> Element: ...

def forgotten_to_schur(a: Element) -> Element: ...

def power_to_schur(a: Element) -> Element: ...

def plethysm(f: Element, g: Element) -> Element:
    """Plethysm f[g] of two Schur-basis elements.
    """
    ...

def kostka_number(la: list[int], mu: list[int]) -> int:
    """Kostka number K_{λμ}.
    """
    ...

def evaluate_schur(a: Element, xs: list[int]) -> int:
    """Evaluate a Schur-basis element at the alphabet `xs`.
    """
    ...

def expand_alphabet(a: Element, src: str, n: int) -> list[tuple[list[int], int]]:
    """The expansion of an element in `n` variables, as `[(exponent vector,
    coefficient), ...]` with each vector of length `n`.
    """
    ...

def monomial_multiply(a: Element, b: Element) -> Element:
    """Multiply two monomial-basis elements.
    """
    ...

def semistandard_tableaux(la: list[int], mu: list[int]) -> list[list[list[int]]]:
    """The semistandard Young tableaux of shape λ and weight μ, as lists of
    rows.
    """
    ...

def dimension(la: list[int]) -> int | None:
    """f^λ — the number of standard Young tableaux of shape λ, i.e. the
    dimension of the irreducible S_{|λ|} representation. `None` past
    `u128`.
    """
    ...

def principal_specialization(la: list[int], n: int) -> int | None:
    """s_λ(1^n), the dimension of the GL_n irreducible. `None` on overflow.
    """
    ...

def principal_specialization_q(la: list[int], n: int) -> list[int]:
    """s_λ(1, q, …, q^{n−1}) as a coefficient list in q, lowest degree
    first.
    """
    ...

def character_value(la: list[int], mu: list[int]) -> int:
    """Symmetric-group character χ^λ(μ).
    """
    ...

def internal_product(a: Element, b: Element) -> Element:
    """The internal (Kronecker) product of two Schur-basis elements.
    """
    ...

def partitions(n: int) -> list[list[int]]:
    """The partitions of `n`, in the order the table functions below index
    by.
    """
    ...

def kronecker_coefficient(la: list[int], mu: list[int], nu: list[int]) -> int:
    """A single Kronecker coefficient g^ν_{λμ}, computed **without forming
    the product**.
    """
    ...

def convert_indexed(a: Element, src: str, dst: str) -> list[tuple[int, int, int]]:
    """A conversion whose output partitions are returned as **indices**
    rather than as lists: `[(degree, index, coefficient), ...]`, where
    `index` is into [`partitions`] of that degree.
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

def omega(a: Element) -> Element:
    """The ω involution on a Schur-basis element.
    """
    ...

def hall_inner_product(a: Element, b: Element) -> int:
    """The Hall inner product of two Schur-basis elements.
    """
    ...

def skew_schur(la: list[int], mu: list[int]) -> Element:
    """The skew Schur function s_{λ/μ}.
    """
    ...

def skew_by(f: Element, g: Element, basis: str = "s") -> Element:
    """Skew a Schur-basis element by `g`, given in `basis` — the adjoint of
    multiplication by g under the Hall inner product.
    """
    ...

def coproduct(a: Element) -> list[tuple[tuple[list[int], list[int]], int]]:
    """The coproduct Δ, as `[((mu, nu), coefficient), ...]`.
    """
    ...

def antipode(a: Element) -> Element:
    """The antipode S.
    """
    ...

def schubert_multiply(a: SchubertElement, b: SchubertElement) -> SchubertElement:
    """Multiply two Schubert polynomials.
    """
    ...

def schubert_multiply_variable(a: SchubertElement, i: int) -> SchubertElement:
    """`x_i · f`, the signed Monk rule. **1-based**, unlike Symmetrica's
    `mult_schubert_variable`, which is 0-based while its own
    `divdiff_schubert` is 1-based. One convention, stated.
    """
    ...

def schubert_divided_difference(a: SchubertElement, i: int) -> SchubertElement:
    """`∂_i f` on the Schubert basis, 1-based.
    """
    ...

def schubert_divided_difference_perm(a: SchubertElement, w: list[int]) -> SchubertElement:
    """`∂_w f`, composing along a reduced word of `w`.
    """
    ...

def schubert_expand(a: SchubertElement) -> Element:
    """Expand into monomials: `(exponent vector, coefficient)` pairs.
    """
    ...

def polynomial_to_schubert(terms: Element) -> SchubertElement:
    """Write a polynomial in the Schubert basis (the greedy triangular
    peel).
    """
    ...

def schubert_pairing(a: SchubertElement, b: SchubertElement, n: int) -> int:
    """The Poincaré pairing on `H*(Fl(n))`.
    """
    ...

def schubert_dimension(w: list[int]) -> int:
    """`S_w(1,…,1)`: the number of pipe dreams, i.e. the size
    `schubert_expand` would produce. Cheap — it never builds the
    expansion.
    """
    ...

def schubert_coefficient(u: list[int], v: list[int], w: list[int]) -> int:
    """A single structure constant `c^w_{uv}`, **without building the
    product**.
    """
    ...

def schubert_monomial_mass(u: list[int], v: list[int]) -> int:
    """The product's total monomial mass `S_u(1,…,1)·S_v(1,…,1)`, a count.
    """
    ...

def schubert_to_stanley_schur(w: list[int]) -> Element:
    """The **Stanley symmetric function** `F_w` in the Schur basis.
    """
    ...

