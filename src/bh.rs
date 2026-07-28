//! The modified Macdonald polynomials `H̃_μ` by the Bergeron–Haiman Pieri
//! recursion.
//!
//! ## Reference
//!
//! F. Bergeron, M. Haiman, *Tableaux formulas for Macdonald polynomials*,
//! International Journal of Algebra and Computation **23** (2013), 833–852.
//! Cited as **[BH]**; the multi-box recursion is their Proposition 5.
//!
//! ## What this computes, and why it is here
//!
//! ```text
//!   H̃_μ = Σ_ν L_{μν}(q,t) m_ν,        L_{μν} = ⟨H̃_μ, h_ν⟩
//! ```
//!
//! and `L_{μν} ∈ ℕ[q,t]`, being `Σ_λ K̃_{λμ} K_{λν}` — a non-negative combination
//! of modified Kostka coefficients. So the (q,t)-Kostka matrix is `L` composed
//! with the ordinary `m → s` transition, which is integral and already in the
//! crate.
//!
//! This is a third route to the same numbers, and the reason to prefer it is
//! measured rather than aesthetic: the branching formula in
//! [`macdonald`](crate::macdonald) grows 5.5× per degree and the
//! Lapointe–Lascoux–Morse eigenvector route in [`macop`](crate::macop) grows
//! 3.1×, while the recursion here is over *pairs of partitions ordered by
//! containment* and shares every value it computes across every μ and ν that
//! reach it.
//!
//! ## The recursion
//!
//! Write `h_r^⊥` for the adjoint of multiplying by `h_r`, and expand its action
//! on the modified basis,
//!
//! ```text
//!   h_r^⊥ H̃_μ = Σ_{γ ⊆ μ, |γ| = |μ| − r} c⁽ʳ⁾_{μγ} H̃_γ.
//! ```
//!
//! Then with `r` the smallest part of ν and `ν̂` the rest, `⟨H̃_μ, h_ν̂ h_r⟩ =
//! ⟨h_r^⊥ H̃_μ, h_ν̂⟩` gives
//!
//! ```text
//!   L_{μν} = Σ_{γ ⊆ μ, |γ| = |ν̂|} c⁽ʳ⁾_{μγ} L_{γν̂},        L_{μ,(n)} = 1.
//! ```
//!
//! The base case is `K̃_{(n)μ} = 1`: `(n)` is dominance-maximal, so `s_{(n)}` is
//! the only Schur function contributing `m_{(n)}`.
//!
//! ## The coefficients
//!
//! One box has a closed form — a product of ratios of `(q,t)`-hook weights over
//! the cells of ν in the row and the column the box was removed from. Every cell
//! outside that row and column has the same arm and leg in both shapes, so its
//! factors cancel and it is not visited.
//!
//! More than one box has no product formula; [BH] Proposition 5 gives a
//! recursion through the **bi-exponent generator** of a skew shape,
//! `B_{μ/ν} = Σ_{(i,j) ∈ μ/ν} t^i q^j`, which peels one box at a time and so
//! bottoms out at the closed form.
//!
//! ## Arithmetic: a second factored fraction field
//!
//! The hook weights are `q^a − t^b`, which is **not** the class
//! [`Frac`](crate::Frac) holds — that one is closed under `1 − qᵃtᵇ` and a
//! weight like `t² − q³` is not of that shape for any exponent pair. [`Rat`]
//! here is the same design over the family that does close: a denominator is a
//! multiset of `(a,b)` meaning `q^a − t^b`, products and lcms stay inside it, and
//! reduction is [`QtPoly::divide_exact`](crate::qt::QtPoly::divide_exact).
//!
//! No field is needed. `q^a − t^b` is lex-monic up to sign — its leading term is
//! `q^a` when `a > 0` and `−t^b` when `a = 0`, coefficient ±1 either way — so the
//! elimination never divides a coefficient by anything but a unit, and `ℤ[q,t]`
//! carries the whole computation.

use std::collections::BTreeMap;

use crate::coeff::Ring;
use crate::partition::Partition;
use crate::qt::QtPoly;

/// `q^a − t^b`, the shape every denominator here factors into.
type Atom = (u32, u32);

fn atom<C: Ring>((a, b): Atom) -> QtPoly<C> {
    let mut p = QtPoly::term(a, 0, C::one());
    p.add_term(0, b, C::one().neg());
    p
}

/// A rational function whose denominator is a product of `q^a − t^b`.
///
/// [`Frac`](crate::Frac) over a different family — see the module docs for why a
/// second one is needed rather than a reuse.
#[derive(Clone, Debug)]
pub(crate) struct Rat<C: Ring> {
    num: QtPoly<C>,
    den: BTreeMap<Atom, u32>,
}

impl<C: Ring> Rat<C> {
    fn one() -> Self {
        Rat {
            num: <QtPoly<C> as Ring>::one(),
            den: BTreeMap::new(),
        }
    }

    fn zero() -> Self {
        Rat {
            num: QtPoly::zero(),
            den: BTreeMap::new(),
        }
    }

    fn from_poly(num: QtPoly<C>) -> Self {
        Rat {
            num,
            den: BTreeMap::new(),
        }
    }

    fn is_zero(&self) -> bool {
        self.num.is_empty()
    }

    /// `self · (q^{a₁} − t^{b₁}) / (q^{a₂} − t^{b₂})`.
    fn scale_by_ratio(&mut self, up: Atom, down: Atom) {
        if up == down {
            return;
        }
        self.num = self.num.mul(&atom(up));
        *self.den.entry(down).or_insert(0) += 1;
    }

    fn mul(&self, other: &Self) -> Self {
        let mut den = self.den.clone();
        for (k, &m) in &other.den {
            *den.entry(*k).or_insert(0) += m;
        }
        Rat {
            num: self.num.mul(&other.num),
            den,
        }
    }

    fn lift(&self, target: &BTreeMap<Atom, u32>) -> QtPoly<C> {
        let mut num = self.num.clone();
        for (&k, &m) in target {
            for _ in 0..(m - self.den.get(&k).copied().unwrap_or(0)) {
                num = num.mul(&atom(k));
            }
        }
        num
    }

    fn add_assign(&mut self, other: &Self) {
        if other.is_zero() {
            return;
        }
        if self.is_zero() {
            *self = other.clone();
            return;
        }
        if other.den != self.den {
            let mut lcm = self.den.clone();
            for (&k, &m) in &other.den {
                let e = lcm.entry(k).or_insert(0);
                *e = (*e).max(m);
            }
            self.num = self.lift(&lcm);
            self.den = lcm;
        }
        let lifted = other.lift(&self.den);
        self.num.add_assign(&lifted);
    }

    /// Divide out every denominator factor that also divides the numerator.
    fn reduce(&mut self) {
        if self.num.is_empty() {
            self.den.clear();
            return;
        }
        self.den.retain(|&k, m| {
            while *m > 0 {
                match self.num.divide_exact(&atom(k)) {
                    Some(q) => {
                        self.num = q;
                        *m -= 1;
                    }
                    None => break,
                }
            }
            *m > 0
        });
    }

    /// The numerator, if the denominator cancelled entirely.
    fn into_poly(mut self) -> Option<QtPoly<C>> {
        self.reduce();
        self.den.is_empty().then_some(self.num)
    }
}

/// `B_{μ/ν} = Σ_{(i,j) ∈ μ/ν} t^i q^j`, the bi-exponent generator of a skew
/// shape ([BH] Proposition 5).
fn bi_exponent<C: Ring>(mu: &Partition, nu: &Partition) -> QtPoly<C> {
    let mut out = QtPoly::zero();
    for i in 0..mu.len() {
        for j in nu.part(i)..mu.part(i) {
            out.add_term(j, i as u32, C::one());
        }
    }
    out
}

/// The one-box coefficient `c⁽¹⁾_{μν}`, for ν = μ with a single box removed.
///
/// A product of ratios of `(q,t)`-hook weights over the cells of ν lying in the
/// row and the column that box vacated. Every other cell has the same arm and
/// leg in both shapes, so its two weights are equal and cancel — which is why
/// this is a walk along one row and one column and not over the whole diagram.
fn one_box<C: Ring>(mu: &Partition, nu: &Partition) -> Rat<C> {
    let (mut row, mut col) = (usize::MAX, u32::MAX);
    for i in 0..mu.len() {
        if mu.part(i) != nu.part(i) {
            row = i;
            col = nu.part(i);
            break;
        }
    }
    debug_assert!(row != usize::MAX, "nu must be mu with a box removed");

    let (m, n) = (mu.parts(), nu.parts());
    let mut out = Rat::one();
    // Along the row the box left: the `t^l − q^{a+1}` weights.
    for j in 0..nu.part(row) as usize {
        let up = (crate::macdonald::arm(m, row, j) + 1, crate::macdonald::leg(m, row, j));
        let down = (crate::macdonald::arm(n, row, j) + 1, crate::macdonald::leg(n, row, j));
        // `t^l − q^{a+1}` is `−(q^{a+1} − t^l)`; the two signs cancel in the
        // ratio, so the atoms are stored unsigned.
        out.scale_by_ratio(up, down);
    }
    // Down the column: the `q^a − t^{l+1}` weights.
    for i in 0..crate::macdonald::count_above(n, col as usize) {
        let up = (crate::macdonald::arm(m, i, col as usize), crate::macdonald::leg(m, i, col as usize) + 1);
        let down = (crate::macdonald::arm(n, i, col as usize), crate::macdonald::leg(n, i, col as usize) + 1);
        out.scale_by_ratio(up, down);
    }
    out.reduce();
    out
}

/// The Pieri and `L` recursions, with their caches.
///
/// Both are keyed on pairs of partitions and both are reached from many places,
/// which is the whole reason this beats a per-shape enumeration: a value like
/// `c⁽¹⁾_{(3,1),(3)}` is shared by every μ and ν whose recursion passes through
/// it, across the entire degree.
struct Recursion<C: Ring> {
    /// `c⁽ʳ⁾_{μν}` with `r = |μ| − |ν|`, so the arity is implied by the key.
    pieri: std::collections::HashMap<(Partition, Partition), Rat<C>>,
    /// `L_{μν}`.
    ell: std::collections::HashMap<(Partition, Partition), Rat<C>>,
}

impl<C: Ring> Recursion<C> {
    fn new() -> Self {
        Recursion {
            pieri: std::collections::HashMap::new(),
            ell: std::collections::HashMap::new(),
        }
    }

    /// `c⁽ʳ⁾_{μν}`, the coefficient of `H̃_ν` in `h_r^⊥ H̃_μ`.
    ///
    /// [BH] Proposition 5 for `r ≥ 2`:
    ///
    /// ```text
    ///   c⁽ʳ⁾_{μν} = ( Σ_{ν ⋖ α ⊆ μ} c⁽ʳ⁻¹⁾_{μα} · c⁽¹⁾_{αν} · B_{α/ν} ) / B_{μ/ν}
    /// ```
    ///
    /// `B_{α/ν}` is a single box, so it is one monomial `t^i q^j`; `B_{μ/ν}`
    /// divides the sum exactly, which is the step that would need a gcd in a
    /// general fraction field and here is
    /// [`divide_exact`](crate::qt::QtPoly::divide_exact).
    fn pieri(&mut self, mu: &Partition, nu: &Partition) -> Rat<C> {
        if mu == nu {
            return Rat::one();
        }
        if !mu.contains(nu) {
            return Rat::zero();
        }
        let key = (mu.clone(), nu.clone());
        if let Some(v) = self.pieri.get(&key) {
            return v.clone();
        }
        let value = if mu.size() == nu.size() + 1 {
            one_box(mu, nu)
        } else {
            let mut acc = Rat::zero();
            for alpha in covers_within(nu, mu) {
                let outer = self.pieri(mu, &alpha);
                if outer.is_zero() {
                    continue;
                }
                let mut term = outer.mul(&one_box(&alpha, nu));
                term = term.mul(&Rat::from_poly(bi_exponent(&alpha, nu)));
                acc.add_assign(&term);
            }
            // Divide by B_{μ/ν} **before** reducing. `Rat`'s denominator is a
            // multiset of whole atoms `q^a − t^b`, and those are not
            // irreducible — `q⁴ − t²` is `(q² − t)(q² + t)` — so `reduce` can
            // cancel a proper factor of an atom against the numerator and leave
            // `B` no longer dividing it. First seen at μ = (5,2,2,2),
            // ν = (4,2,1), which is degree 11: everything below that is clean
            // either way.
            let b = bi_exponent::<C>(mu, nu);
            acc.num = acc.num.divide_exact(&b).unwrap_or_else(|| {
                panic!("B_{{mu/nu}} must divide the Pieri sum at {mu} / {nu}")
            });
            acc.reduce();
            acc
        };
        self.pieri.insert(key, value.clone());
        value
    }

    /// `L_{μν} = ⟨H̃_μ, h_ν⟩`.
    fn ell(&mut self, mu: &Partition, nu: &Partition) -> Rat<C> {
        debug_assert_eq!(mu.size(), nu.size());
        if nu.len() <= 1 {
            return Rat::one();
        }
        let key = (mu.clone(), nu.clone());
        if let Some(v) = self.ell.get(&key) {
            return v.clone();
        }
        // Peel the smallest part of ν.
        let r = nu.part(nu.len() - 1);
        let hat = Partition::new(nu.parts()[..nu.len() - 1].iter().copied());
        let mut acc = Rat::zero();
        for gamma in crate::partitions_of(hat.size()) {
            if !mu.contains(&gamma) {
                continue;
            }
            let c = self.pieri(mu, &gamma);
            if c.is_zero() {
                continue;
            }
            let l = self.ell(&gamma, &hat);
            acc.add_assign(&c.mul(&l));
        }
        acc.reduce();
        debug_assert_eq!(r, nu.part(nu.len() - 1));
        self.ell.insert(key, acc.clone());
        acc
    }
}

/// Every partition covering `nu` (one box more) and still inside `mu`.
fn covers_within(nu: &Partition, mu: &Partition) -> Vec<Partition> {
    let mut out = Vec::new();
    for i in 0..=nu.len() {
        let mut parts: Vec<u32> = nu.parts().to_vec();
        if i == parts.len() {
            parts.push(1);
        } else {
            parts[i] += 1;
        }
        // Still weakly decreasing, and still inside mu.
        if i > 0 && parts[i] > parts[i - 1] {
            continue;
        }
        let alpha = Partition::new(parts);
        if mu.contains(&alpha) {
            out.push(alpha);
        }
    }
    out
}

/// `H̃_μ = Σ_ν L_{μν} m_ν`, the modified Macdonald polynomial in the monomial
/// basis, for every μ of the degree.
///
/// The whole degree at once, because the recursions share: `L_{γν̂}` for `γ ⊆ μ`
/// is reached from every μ containing γ, and `c⁽ʳ⁾` from more places still.
pub fn htilde_monomial_table<C: Ring>(n: u32) -> Vec<(Partition, crate::sym::Monomial<QtPoly<C>>)> {
    use crate::sym::SymFn;
    let parts = crate::memo::partitions_cached(n);
    let mut rec: Recursion<C> = Recursion::new();
    parts
        .iter()
        .map(|mu| {
            let mut out = crate::sym::Monomial::zero();
            for nu in parts.iter() {
                let l = rec.ell(mu, nu);
                let poly = l
                    .into_poly()
                    .unwrap_or_else(|| panic!("L_{{{mu},{nu}}} must be a polynomial"));
                out.add_term(nu.clone(), poly);
            }
            (mu.clone(), out)
        })
        .collect()
}

/// `H̃_μ` in the **Schur** basis for a whole degree — the modified (q,t)-Kostka
/// coefficients `K̃_{λμ}`.
///
/// `m → s` is the inverse Kostka transition, which is integral, so no
/// `(q,t)`-arithmetic happens in it and the whole route stays in ℤ[q,t]. That is
/// why this is bounded on [`Ring`] where the other two routes to the same
/// numbers need [`QAlgebra`](crate::coeff::QAlgebra): neither the recursion nor
/// the basis change ever divides by an integer.
pub fn htilde_table<C: Ring>(n: u32) -> Vec<(Partition, crate::sym::Schur<QtPoly<C>>)> {
    use crate::convert::ToSchur;
    htilde_monomial_table::<C>(n)
        .into_iter()
        .map(|(mu, m)| (mu, m.to_schur()))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;
    use crate::sym::SymFn;

    /// Against the route that is already checked against Sage.
    ///
    /// [`macdonald_ht`](crate::macdonald_ht) reaches `H̃` by enumerating tableaux
    /// for `J`, inverting the `S` basis and reflecting the `t`-exponents; this
    /// reaches it by a Pieri recursion on pairs of partitions ordered by
    /// containment. They share `Partition` and `QtPoly`.
    #[test]
    fn the_recursion_agrees_with_the_branching_route() {
        for n in 1..=6u32 {
            for (mu, got) in htilde_table::<Rational>(n) {
                let want = crate::macdonald_ht::<Rational>(&mu);
                for lambda in crate::partitions_of(n) {
                    assert_eq!(got.coeff(&lambda), want.coeff(&lambda), "K~_{{{lambda},{mu}}}");
                }
            }
        }
    }

    /// Haiman: `K̃_{λμ} ∈ ℕ[q,t]`.
    ///
    /// Nothing here arranges it — the Pieri coefficients are rational functions
    /// with hook-weight denominators, and every one of them has to cancel. A
    /// negative coefficient or a surviving denominator is a bug, and
    /// `into_poly` already refuses the second.
    #[test]
    fn the_coefficients_are_non_negative_integers() {
        for n in 1..=6u32 {
            for (mu, s) in htilde_table::<Rational>(n) {
                for (lambda, k) in s.terms() {
                    for (_, c) in k.terms() {
                        assert_eq!(c.denom(), 1, "K~_{{{lambda},{mu}}} has {c:?}");
                        assert!(c.numer() > 0, "K~_{{{lambda},{mu}}} has {c:?}");
                    }
                }
            }
        }
    }

    /// `L_{μ,(n)} = 1` and, more usefully, `Σ_ν L_{μν}` is not what is checked —
    /// the monomial expansion must have `L_{μν} = Σ_λ K̃_{λμ} K_{λν}`.
    ///
    /// That identity is the definition `L = ⟨H̃, h⟩` written out, and it ties the
    /// two halves of this module together: get the recursion right and the `m → s`
    /// step wrong and it fails, and vice versa.
    #[test]
    fn the_monomial_coefficients_pair_against_the_kostka_numbers() {
        for n in 1..=6u32 {
            let parts = crate::partitions_of(n);
            let schur = htilde_table::<Rational>(n);
            for ((mu, s), (mu2, m)) in schur.iter().zip(htilde_monomial_table::<Rational>(n)) {
                assert_eq!(mu, &mu2);
                for nu in &parts {
                    let mut want = QtPoly::zero();
                    for lambda in &parts {
                        let k = crate::kostka::kostka(lambda, nu);
                        if k == 0 {
                            continue;
                        }
                        let scale = QtPoly::term(0, 0, Rational::from_int(k as i128));
                        want.add_assign(&s.coeff(lambda).mul(&scale));
                    }
                    assert_eq!(m.coeff(nu), want, "L_{{{mu},{nu}}}");
                }
            }
        }
    }
}
