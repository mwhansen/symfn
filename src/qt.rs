//! Bivariate `(q, t)` polynomials as a coefficient ring.
//!
//! This is the coefficient type the non-classical families need:
//! Hall–Littlewood lives over `ℤ[t]`, Macdonald over ℚ(q,t), Jack over ℚ(α).
//! [`QtPoly`] covers the polynomial half of that — the fraction field is a
//! separate layer on top, and is only needed once Macdonald arrives.
//!
//! ## Why this can exist at all
//!
//! `ℚ[q,t]` is **not a field**, and until the dividing paths were re-bounded on
//! [`QAlgebra`] rather than
//! [`Field`](crate::coeff::Field) they were unavailable over it. Every division
//! in this library is by z_μ — an integer — so a ring containing ℚ suffices.
//! That is what makes `s → p`, the internal product and plethysm work here.
//!
//! ## Sparse, and generic over the coefficients
//!
//! Terms are held as a **sorted `Vec`** of (exponent pair, coefficient), with
//! zeros removed, so equality is structural and a polynomial costs what it
//! uses. Hall–Littlewood and Macdonald expansions are sparse in (q, t) —
//! Kostka– Foulkes polynomials in particular have few terms relative to their
//! degree — so a dense representation would mostly store zeros.
//!
//! A `BTreeMap` is the obvious choice and the wrong one: these polynomials
//! hold a handful of terms, and at that size a B-tree pays a node allocation
//! and a pointer chase for what a `Vec` does in one cache line. The insert is a
//! memmove instead of a rebalance, and the whole polynomial is a single
//! allocation.
//!
//! The coefficient ring is a parameter for the same reason it is everywhere
//! else here: `QtPoly<i64>` is `ℤ[q,t]` for exact small work,
//! `QtPoly<Rational>` is `ℚ[q,t]`, and `QtPoly<BigInt>` (under `bignum`) has no
//! ceiling. Kostka–Foulkes coefficients are integers, Macdonald's are not, and
//! neither should force the other's representation.
//!
//! ## Plethysm acts on q and t
//!
//! [`Plethystic::frobenius`] raises the *variables*, so `p_n` sends q^a t^b to
//! q^{an} t^{bn}. This is Sage's default convention, and getting it wrong is
//! invisible over ℚ — see [`crate::plethysm`](mod@crate::plethysm).

use core::fmt;

use crate::coeff::{Plethystic, QAlgebra, Ring};

/// A polynomial in `q` and `t`, sparse in the exponent pair `(a, b)` for
/// `q^a t^b`, with coefficients in `C`.
///
/// Never stores a zero coefficient, so `PartialEq` is mathematical equality.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct QtPoly<C: Ring>(Vec<((u32, u32), C)>);

impl<C: Ring> QtPoly<C> {
    /// `c · q^a t^b`.
    pub fn term(a: u32, b: u32, c: C) -> Self {
        QtPoly(if c.is_zero() {
            Vec::new()
        } else {
            vec![((a, b), c)]
        })
    }

    /// The variable `q`.
    pub fn q() -> Self {
        Self::term(1, 0, C::one())
    }

    /// The variable `t`.
    pub fn t() -> Self {
        Self::term(0, 1, C::one())
    }

    /// Coefficient of `q^a t^b`.
    pub fn coeff(&self, a: u32, b: u32) -> C {
        match self.0.binary_search_by_key(&(a, b), |e| e.0) {
            Ok(i) => self.0[i].1.clone(),
            Err(_) => C::zero(),
        }
    }

    /// The terms, ascending by exponent pair.
    pub fn terms(&self) -> impl Iterator<Item = (&(u32, u32), &C)> {
        self.0.iter().map(|(k, c)| (k, c))
    }

    /// The number of terms. Explicit zeros are never stored, so this counts
    /// the nonzero ones.
    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// Whether this is the zero polynomial.
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// `self += c · q^a t^b`, dropping the term if it cancels.
    pub fn add_term(&mut self, a: u32, b: u32, c: C) {
        if c.is_zero() {
            return;
        }
        match self.0.binary_search_by_key(&(a, b), |e| e.0) {
            Ok(i) => {
                self.0[i].1.add_assign(&c);
                if self.0[i].1.is_zero() {
                    self.0.remove(i);
                }
            }
            Err(i) => self.0.insert(i, ((a, b), c)),
        }
    }

    /// `self += ± q^{sa} t^{sb} · other`, in one merging pass.
    ///
    /// Worth its own method because of *how* the terms arrive: `other` is
    /// already sorted, and a uniform shift is monotone for the lexicographic
    /// key — `(x,y) < (x',y')` implies `(x+sa, y+sb) < (x'+sa, y'+sb)` — so the
    /// incoming run is a sorted sequence being merged into a sorted sequence.
    /// Adding it term by term instead costs a binary search and a memmove each,
    /// which is enough to give away everything the `Vec` representation buys
    /// (`docs/record/hall-littlewood.md`).
    pub(crate) fn add_shifted(&mut self, other: &Self, shift: (u32, u32), negate: bool) {
        if other.0.is_empty() {
            return;
        }
        let (sa, sb) = shift;
        let signed = |c: &C| if negate { c.neg() } else { c.clone() };
        if self.0.is_empty() {
            self.0 = other
                .0
                .iter()
                .map(|((a, b), c)| ((a + sa, b + sb), signed(c)))
                .collect();
            return;
        }
        let mine = core::mem::take(&mut self.0);
        let mut out = Vec::with_capacity(mine.len() + other.0.len());
        let (mut i, mut j) = (0, 0);
        while i < mine.len() && j < other.0.len() {
            let key = (other.0[j].0 .0 + sa, other.0[j].0 .1 + sb);
            match mine[i].0.cmp(&key) {
                core::cmp::Ordering::Less => {
                    out.push(mine[i].clone());
                    i += 1;
                }
                core::cmp::Ordering::Greater => {
                    out.push((key, signed(&other.0[j].1)));
                    j += 1;
                }
                core::cmp::Ordering::Equal => {
                    let mut v = mine[i].1.clone();
                    v.add_assign(&signed(&other.0[j].1));
                    if !v.is_zero() {
                        out.push((key, v));
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        out.extend_from_slice(&mine[i..]);
        out.extend(
            other.0[j..]
                .iter()
                .map(|((a, b), c)| ((a + sa, b + sb), signed(c))),
        );
        self.0 = out;
    }

    /// `self += c · q^{sa} t^{sb} · other`, in one merging pass.
    ///
    /// [`add_shifted`](Self::add_shifted) with an arbitrary coefficient instead
    /// of a sign, and kept separate from it because that one multiplies
    /// nothing: negation is a `neg` per term where this is a `mul`, and the
    /// binomial paths that dominate Macdonald should not pay for a coefficient
    /// they know is ±1.
    pub(crate) fn add_scaled_shifted(&mut self, other: &Self, shift: (u32, u32), c: &C) {
        if other.0.is_empty() || c.is_zero() {
            return;
        }
        let (sa, sb) = shift;
        if self.0.is_empty() {
            self.0 = other
                .0
                .iter()
                .map(|((a, b), v)| ((a + sa, b + sb), v.mul(c)))
                .collect();
            return;
        }
        let mine = core::mem::take(&mut self.0);
        let mut out = Vec::with_capacity(mine.len() + other.0.len());
        let (mut i, mut j) = (0, 0);
        while i < mine.len() && j < other.0.len() {
            let key = (other.0[j].0 .0 + sa, other.0[j].0 .1 + sb);
            match mine[i].0.cmp(&key) {
                core::cmp::Ordering::Less => {
                    out.push(mine[i].clone());
                    i += 1;
                }
                core::cmp::Ordering::Greater => {
                    out.push((key, other.0[j].1.mul(c)));
                    j += 1;
                }
                core::cmp::Ordering::Equal => {
                    let mut v = mine[i].1.clone();
                    v.add_assign(&other.0[j].1.mul(c));
                    if !v.is_zero() {
                        out.push((key, v));
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        out.extend_from_slice(&mine[i..]);
        out.extend(
            other.0[j..]
                .iter()
                .map(|((a, b), v)| ((a + sa, b + sb), v.mul(c))),
        );
        self.0 = out;
    }

    /// The terms as a slice, ascending by exponent pair.
    pub(crate) fn raw(&self) -> &[((u32, u32), C)] {
        &self.0
    }

    /// Build from terms already ascending by exponent pair, with no zeros.
    ///
    /// The caller owes both invariants; nothing else in the type checks them.
    pub(crate) fn from_sorted(terms: Vec<((u32, u32), C)>) -> Self {
        debug_assert!(
            terms.windows(2).all(|w| w[0].0 < w[1].0),
            "terms must be strictly ascending"
        );
        debug_assert!(terms.iter().all(|(_, c)| !c.is_zero()), "no zero terms");
        QtPoly(terms)
    }

    /// Multiply by the binomial `1 − qᵃtᵇ`.
    ///
    /// `self − qᵃtᵇ·self` is a sorted run merged with a shifted copy of itself,
    /// so it is one pass over `self` twice and one allocation — no sort, and no
    /// intermediate copy.
    ///
    /// This is not a micro-optimisation of [`Ring::mul`], it is the whole
    /// workload. Instrumenting `macdonald_p` at degree 9 found **every one of
    /// its 100k `mul` calls had a two-term operand**, averaging 129 terms on
    /// the other side: [`Frac`](crate::Frac) multiplies by a binomial and never
    /// by anything else, because `from_factors`, `lift` and `denominator` build
    /// products of `1 − qᵃtᵇ` and nothing more. The general `mul` collected 258
    /// products into a `Vec` and quicksorted it, putting `quicksort` +
    /// `small_sort` high in the profile for what is a concatenation of two
    /// already-sorted runs (`docs/record/macdonald-operators.md`).
    /// # Panics
    ///
    /// Panics if `a == 0 && b == 0`. The factor would be `1 − q⁰t⁰ = 0`, so a
    /// caller that reaches it built a degenerate atom upstream — the same
    /// requirement [`Atom::unit`](crate::deltaop::Atom::unit) states.
    pub fn mul_binomial(&self, a: u32, b: u32) -> Self {
        assert!(a > 0 || b > 0, "1 - q^0 t^0 is zero");
        let n = self.0.len();
        // The shifted copy is `self` read with `(a, b)` added to every key, so
        // both runs are the same slice walked at two offsets.
        let mut out: Vec<((u32, u32), C)> = Vec::with_capacity(n + 1);
        let (mut i, mut j) = (0, 0);
        while i < n && j < n {
            let key = (self.0[j].0 .0 + a, self.0[j].0 .1 + b);
            match self.0[i].0.cmp(&key) {
                core::cmp::Ordering::Less => {
                    out.push(self.0[i].clone());
                    i += 1;
                }
                core::cmp::Ordering::Greater => {
                    out.push((key, self.0[j].1.neg()));
                    j += 1;
                }
                core::cmp::Ordering::Equal => {
                    let mut v = self.0[i].1.clone();
                    v.add_assign(&self.0[j].1.neg());
                    if !v.is_zero() {
                        out.push((key, v));
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        out.extend_from_slice(&self.0[i..]);
        out.extend(
            self.0[j..]
                .iter()
                .map(|((x, y), c)| ((x + a, y + b), c.neg())),
        );
        QtPoly(out)
    }

    /// Multiply by `qᵃ − tᵇ`, the other two-term factor this library divides
    /// by.
    ///
    /// [`mul_binomial`](Self::mul_binomial) is this for `1 − qᵃtᵇ`, and the
    /// argument for having both is the same one, measured again. `q^a·self` and
    /// `t^b·self` are the same sorted run read at two different uniform shifts,
    /// and a uniform shift preserves the lexicographic order. So the product is
    /// a **merge of two sorted runs** rather than a general product that
    /// collects `2n` pairs and sorts them.
    ///
    /// [`deltaop`](crate::deltaop) is what needs it: `w_μ` factors into this
    /// family, and lifting an accumulator to a common denominator multiplies by
    /// these atoms over and over.
    ///
    /// ⚠️ **It bought nothing on its own, and is kept anyway.** Introduced on
    /// the reasoning above — that [`Ring::mul`] would quicksort a concatenation
    /// of two sorted runs, exactly what [`mul_binomial`](Self::mul_binomial)'s
    /// notes record for Macdonald `P` — it moved `∇e_12` not at all
    /// (`docs/record/macdonald-operators.md`). The sort really was a third of
    /// that profile, but it was a
    /// *different* product: `deltaop` was lifting its accumulator to the common
    /// denominator and only then multiplying by a `K̃` entry, so the big
    /// operand was in the general `mul` and not here. Reordering those two
    /// fixed it. Recorded because the reasoning was sound, the measurement
    /// still said no, and the honest conclusion is that this is the right
    /// primitive for a cost that lives somewhere else.
    /// # Panics
    ///
    /// Panics if `a == 0 && b == 0`. The factor would be `q⁰ − t⁰ = 0`; see
    /// [`mul_binomial`](Self::mul_binomial).
    pub fn mul_diff(&self, a: u32, b: u32) -> Self {
        assert!(a > 0 || b > 0, "q^0 - t^0 is zero");
        let n = self.0.len();
        let mut out: Vec<((u32, u32), C)> = Vec::with_capacity(2 * n);
        let (mut i, mut j) = (0, 0);
        while i < n && j < n {
            let ki = (self.0[i].0 .0 + a, self.0[i].0 .1);
            let kj = (self.0[j].0 .0, self.0[j].0 .1 + b);
            match ki.cmp(&kj) {
                core::cmp::Ordering::Less => {
                    out.push((ki, self.0[i].1.clone()));
                    i += 1;
                }
                core::cmp::Ordering::Greater => {
                    out.push((kj, self.0[j].1.neg()));
                    j += 1;
                }
                core::cmp::Ordering::Equal => {
                    let mut v = self.0[i].1.clone();
                    v.add_assign(&self.0[j].1.neg());
                    if !v.is_zero() {
                        out.push((ki, v));
                    }
                    i += 1;
                    j += 1;
                }
            }
        }
        out.extend(
            self.0[i..]
                .iter()
                .map(|((x, y), c)| ((x + a, *y), c.clone())),
        );
        out.extend(self.0[j..].iter().map(|((x, y), c)| ((*x, y + b), c.neg())));
        QtPoly(out)
    }

    /// Exact division: `Some(q)` with `self == q * d`, or `None` if `d` does
    /// not divide `self` (including `d == 0`).
    ///
    /// Division, not a gcd — the quotient is assumed to exist and the routine
    /// only finds it. That is why this is affordable in a ring where
    /// gcd is not: [`Frac`](crate::Frac) exists precisely because
    /// bivariate polynomial gcd is a real algorithm, but *this* is leading-term
    /// elimination, and [`divide_by_factor`](crate::frac) is already its
    /// special case for `d = 1 − qᵃtᵇ`.
    ///
    /// ## Why the leading term is well defined
    ///
    /// Terms are sorted lexicographically on the exponent pair, which is a
    /// monomial order: a well-order that a uniform shift preserves (the same
    /// fact [`mul_binomial`](Self::mul_binomial) rests on). So `lt(qd) =
    /// lt(q)·lt(d)`, the leading term of the remainder must be divisible by
    /// `lt(d)` at every step, and subtracting `m·d` strictly lowers it. The
    /// quotient monomials therefore come out in **descending** order and are
    /// reversed once at the end rather than sorted.
    ///
    /// ## Failure is failure, never a wrong answer
    ///
    /// Every way out is `None`: an exponent that would go negative, a
    /// coefficient [`Ring::div_exact`] declines, or a remainder that is not
    /// empty when the leading terms run out. A coefficient ring that does not
    /// implement `div_exact` at all makes this conservative — it will report
    /// some genuine divisions as failures — but it cannot make it wrong.
    pub fn divide_exact(&self, d: &Self) -> Option<Self> {
        let (dkey, dcoeff) = d.0.last()?; // `None` on d == 0
        if self.0.is_empty() {
            return Some(QtPoly(Vec::new()));
        }
        // The remainder needs max-extraction and arbitrary-key subtraction, so
        // it is a map here and not the sorted `Vec` the type normally uses.
        let mut rem: std::collections::BTreeMap<(u32, u32), C> = self.0.iter().cloned().collect();
        let mut quot: Vec<((u32, u32), C)> = Vec::new();

        while let Some((k, c)) = rem.iter().next_back() {
            let (rkey, rcoeff) = (*k, c.clone());
            if rkey.0 < dkey.0 || rkey.1 < dkey.1 {
                return None; // leading monomial is not a multiple
            }
            let m = (rkey.0 - dkey.0, rkey.1 - dkey.1);
            let c = rcoeff.div_exact(dcoeff)?;
            for (k, dc) in &d.0 {
                let key = (k.0 + m.0, k.1 + m.1);
                let sub = c.mul(dc).neg();
                match rem.entry(key) {
                    std::collections::btree_map::Entry::Occupied(mut e) => {
                        e.get_mut().add_assign(&sub);
                        if e.get().is_zero() {
                            e.remove();
                        }
                    }
                    std::collections::btree_map::Entry::Vacant(e) => {
                        e.insert(sub);
                    }
                }
            }
            quot.push((m, c));
        }
        quot.reverse();
        Some(QtPoly::from_sorted(quot))
    }

    /// Multiply by `t^b`.
    ///
    /// A uniform exponent shift is injective on monomials, so nothing can
    /// collide and the map is rebuilt directly — no accumulation, no zero
    /// checks. [`Ring::mul`] against `t^b` would allocate a temporary and walk
    /// the general double loop for what is a rename of the keys.
    pub fn shift_t(&self, b: u32) -> Self {
        if b == 0 {
            return self.clone();
        }
        // Shifting t leaves the q exponent alone and is monotone in the t
        // exponent, so the sort order is preserved and no re-sort is needed.
        QtPoly(
            self.0
                .iter()
                .map(|((x, y), c)| ((*x, y + b), c.clone()))
                .collect(),
        )
    }

    /// Substitute numbers for `q` and `t`.
    ///
    /// The specialisations that matter are exactly this: Hall–Littlewood at
    /// t = 0 is Schur and at t = 1 is monomial, which is how the family gets
    /// checked against bases that already have oracles.
    pub fn eval(&self, q: &C, t: &C) -> C {
        let mut total = C::zero();
        for ((a, b), c) in &self.0 {
            let mut term = c.clone();
            for _ in 0..*a {
                term = term.mul(q);
            }
            for _ in 0..*b {
                term = term.mul(t);
            }
            total.add_assign(&term);
        }
        total
    }

    /// Highest `q` and `t` exponents present, or `None` when zero.
    pub fn degrees(&self) -> Option<(u32, u32)> {
        let a = self.0.iter().map(|e| e.0 .0).max()?;
        let b = self.0.iter().map(|e| e.0 .1).max()?;
        Some((a, b))
    }
}

impl<C: Ring> Ring for QtPoly<C> {
    fn zero() -> Self {
        QtPoly(Vec::new())
    }
    fn one() -> Self {
        Self::term(0, 0, C::one())
    }
    fn is_zero(&self) -> bool {
        self.0.is_empty()
    }
    fn add_assign(&mut self, other: &Self) {
        // The same merge `add_shifted` performs, with no shift. Adding term by
        // term instead is a binary search and a memmove each — fine for the two
        // or three terms a Hall-Littlewood coefficient holds, quadratic for the
        // hundreds a Macdonald numerator holds.
        self.add_shifted(other, (0, 0), false);
    }
    fn mul(&self, other: &Self) -> Self {
        if self.0.is_empty() || other.0.is_empty() {
            return Self::zero();
        }
        // A merge per term of the smaller operand, not one big sort.
        //
        // A uniform shift is monotone for the lexicographic key, so `q^a t^b ·
        // other` is *already sorted* — the same fact `mul_binomial` rests on.
        // The product is therefore a merge of `self.len()` sorted runs, and
        // `add_scaled_shifted` folds one in per pass.
        //
        // This replaced collect-and-sort, which was itself a fix for
        // accumulating with `add_term` (that visited keys in no useful order and
        // put 2654 samples in `memmove` against 247 in the multiplication). The
        // sort was the right answer for the workload of the time, where one
        // operand was almost always a binomial and `mul_binomial` now handles
        // that case anyway. On the Lapointe–Lascoux–Morse solve, where both
        // operands are general, sorting was **64% of the profile** —
        // `quicksort` and `small_sort` between them — for a product whose terms
        // arrive in sorted runs.
        //
        // Driving from the smaller side keeps the number of merges down; each
        // costs O(|accumulator| + |other|).
        let (driver, run) = if self.0.len() <= other.0.len() {
            (self, other)
        } else {
            (other, self)
        };
        let mut acc = QtPoly(Vec::new());
        for (k, c) in &driver.0 {
            acc.add_scaled_shifted(run, *k, c);
        }
        acc
    }
    fn neg(&self) -> Self {
        QtPoly(self.0.iter().map(|(k, c)| (*k, c.neg())).collect())
    }
    fn from_i64(n: i64) -> Self {
        Self::term(0, 0, C::from_i64(n))
    }
    fn from_u128(n: u128) -> Self {
        Self::term(0, 0, C::from_u128(n))
    }
    fn from_i128(n: i128) -> Self {
        Self::term(0, 0, C::from_i128(n))
    }

    // `as_ratio` / `from_ratio` are deliberately left declining. They would have
    // to answer in `i128`, which can only represent a *constant* polynomial, and
    // a batch mixing constants with genuine polynomials would have to fall back
    // anyway. Declining keeps `convert::integral_sweep` on its generic path,
    // which is always correct.
}

impl<C: QAlgebra> QAlgebra for QtPoly<C> {
    /// Coefficientwise, and exact: `ℚ[q,t]` contains ℚ, so dividing by an
    /// integer never needs an inverse of q or t. This is why the bound is
    /// `QAlgebra` and not `Field` — see the module docs.
    fn div_u128(&self, n: u128) -> Self {
        QtPoly(self.0.iter().map(|(k, c)| (*k, c.div_u128(n))).collect())
    }
}

impl<C: Plethystic> Plethystic for QtPoly<C> {
    /// `p_n` raises the variables: q^a t^b ↦ q^{an} t^{bn}, and the
    /// coefficients are pushed through their own Frobenius.
    ///
    /// Exponents only grow, so the map is injective on monomials and no two
    /// terms can collide — the result is built directly rather than
    /// accumulated.
    fn frobenius(&self, n: u32) -> Self {
        QtPoly(
            self.0
                .iter()
                .map(|((a, b), c)| ((a * n, b * n), c.frobenius(n)))
                .collect(),
        )
    }
}

impl<C: Ring> fmt::Display for QtPoly<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.is_empty() {
            return f.write_str("0");
        }
        let mut parts = Vec::new();
        for ((a, b), c) in &self.0 {
            let mut s = format!("{c:?}");
            if *a > 0 {
                s.push_str(&if *a == 1 {
                    "*q".into()
                } else {
                    format!("*q^{a}")
                });
            }
            if *b > 0 {
                s.push_str(&if *b == 1 {
                    "*t".into()
                } else {
                    format!("*t^{b}")
                });
            }
            parts.push(s);
        }
        f.write_str(&parts.join(" + "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;
    use crate::convert::{FromSchur, ToSchur};
    use crate::partition::Partition;
    use crate::sym::{PowerSum, Schur, SymFn};

    type P = QtPoly<i64>;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    /// A spread of polynomials to run division against — dense blocks, sparse
    /// ones with gaps, pure powers, and things with negative coefficients.
    fn division_cases() -> Vec<P> {
        let mut cases: Vec<P> = vec![
            <P as Ring>::one(),
            QtPoly::q(),
            QtPoly::t(),
            QtPoly::term(3, 4, 7),
        ];
        let mut block: P = QtPoly::zero();
        for a in 0..4 {
            for b in 0..4 {
                block.add_term(a, b, (a as i64) + (b as i64) - 3);
            }
        }
        cases.push(block);
        let mut sparse: P = QtPoly::zero();
        sparse.add_term(0, 0, 1);
        sparse.add_term(5, 0, -2);
        sparse.add_term(0, 5, 3);
        sparse.add_term(7, 9, -1);
        cases.push(sparse);
        cases
    }

    /// `divide_exact` must invert multiplication, which is the only property it
    /// promises: `(a·b) / b == a` for every pair, with both orders tried since
    /// the leading term of the divisor is what drives the elimination and the
    /// two operands have different ones.
    #[test]
    fn division_inverts_multiplication() {
        let cases = division_cases();
        for a in &cases {
            for b in &cases {
                if b.is_empty() {
                    continue;
                }
                let prod = a.mul(b);
                assert_eq!(
                    prod.divide_exact(b).as_ref(),
                    Some(a),
                    "({a}) * ({b}) / ({b})"
                );
                assert_eq!(
                    prod.divide_exact(a).as_ref(),
                    (!a.is_empty()).then_some(b),
                    "({a}) * ({b}) / ({a})"
                );
            }
        }
    }

    /// Every way out is `None`, and none of them is a wrong quotient.
    #[test]
    fn division_declines_rather_than_guessing() {
        let one = <P as Ring>::one();
        // Division by zero.
        assert_eq!(one.divide_exact(&QtPoly::zero()), None);
        // Zero divided by anything is zero, and needs no elimination at all.
        assert_eq!(QtPoly::zero().divide_exact(&one), Some(QtPoly::zero()));
        // A monomial the divisor cannot reach: q does not divide t.
        let (t, q): (P, P) = (QtPoly::t(), QtPoly::q());
        assert_eq!(t.divide_exact(&q), None);
        // Divides the leading term but leaves a remainder: (1 + q + t) / (1 + q).
        let mut n: P = QtPoly::zero();
        n.add_term(0, 0, 1);
        n.add_term(1, 0, 1);
        n.add_term(0, 1, 1);
        let mut d: P = QtPoly::zero();
        d.add_term(0, 0, 1);
        d.add_term(1, 0, 1);
        assert_eq!(n.divide_exact(&d), None);
    }

    /// The divisor this routine exists for: `v_λ = ∏_{μ≠λ}([|λ|] − [|μ|])`, the
    /// normalising factor of Lapointe–Lascoux–Morse, *Determinantal expressions
    /// for Macdonald polynomials* (IMRN **1998** no. 18, 957–978;
    /// arXiv:math/9808050), their 3.10 — over **ℤ** and not ℚ.
    ///
    /// `[|α|] = Σ_i q^{α_i} t^{n−i}` is the eigenvalue of the Macdonald operator
    /// `M₁`, and the plan for replacing the branching formula is to clear those
    /// denominators, work in `ℤ[q,t]`, and divide them back out at the end. So
    /// the question this test answers is not "does division work" but "does it
    /// work on *that*, without a field".
    ///
    /// It does. Every coefficient of `[|λ|] − [|μ|]` is ±1: two of its
    /// monomials could only collide if `λ_i = μ_i` for the same `i`, and then
    /// they cancel to nothing rather than accumulating. Lex order is
    /// multiplicative, so a product of such factors still has leading
    /// coefficient ±1, and the elimination never needs to divide a coefficient
    /// by anything but a unit. The assertion
    /// below states that directly — if it ever fails, `Ring::div_exact` over ℤ
    /// starts declining and the whole route needs ℚ.
    #[test]
    fn division_handles_the_macdonald_eigenvalue_products() {
        // The real one, not a copy: a private duplicate here would keep passing
        // if `macop` changed its indexing, which is exactly the mistake that
        // convention has already caused once.
        let eigenvalue = |p: &Partition, k: usize| crate::macop::eigenvalue_of::<i64>(p, k);
        for n in 2..=6u32 {
            let parts = crate::partitions_of(n);
            let k = n as usize;
            for lambda in &parts {
                let ev_l = eigenvalue(lambda, k);
                let mut v: P = <P as Ring>::one();
                for mu in &parts {
                    if mu == lambda {
                        continue;
                    }
                    let mut diff = ev_l.clone();
                    diff.sub_assign(&eigenvalue(mu, k));
                    assert!(
                        diff.terms().all(|(_, c)| c.abs() == 1),
                        "[|{lambda}|] - [|{mu}|] should have unit coefficients"
                    );
                    v = v.mul(&diff);
                }
                let lead = v.raw().last().expect("v is nonzero").1;
                assert_eq!(lead.abs(), 1, "v_{lambda} should be lex-monic up to sign");

                // The round trip the algorithm will actually perform.
                for x in division_cases() {
                    assert_eq!(
                        v.mul(&x).divide_exact(&v),
                        Some(x.clone()),
                        "v_{lambda} dividing its own multiple"
                    );
                }
            }
        }
    }

    /// The coefficient ring decides, and ℤ is not ℚ.
    ///
    /// `3q / 2` has no answer over `QtPoly<i64>` and a perfectly good one over
    /// `QtPoly<Rational>`. `Ring::div_exact` is the seam that reports the
    /// difference instead of truncating, which is the failure a plain `/` on
    /// integers would have produced silently.
    #[test]
    fn division_respects_the_coefficient_ring() {
        let three_q: P = QtPoly::term(1, 0, 3);
        let two: P = QtPoly::term(0, 0, 2);
        assert_eq!(three_q.divide_exact(&two), None, "3q/2 is not in Z[q,t]");

        let three_q: QtPoly<Rational> = QtPoly::term(1, 0, Rational::from_int(3));
        let two: QtPoly<Rational> = QtPoly::term(0, 0, Rational::from_int(2));
        let want = QtPoly::term(1, 0, Rational::new(3, 2));
        assert_eq!(three_q.divide_exact(&two), Some(want), "3q/2 is in Q[q,t]");
    }

    /// `mul_binomial` must agree with the general product, term for term, on a
    /// spread of shapes — including ones where the shift makes terms collide
    /// and cancel, which is the only place a merge can go wrong that a double
    /// loop cannot.
    #[test]
    fn mul_binomial_agrees_with_the_general_product() {
        let mut cases: Vec<P> = vec![
            <P as Ring>::one(),
            QtPoly::q(),
            QtPoly::t(),
            QtPoly::term(3, 4, 7),
        ];
        // 1 + q + q^2 + ... + t + qt + ...: a dense block, where every shift
        // lands on an existing term.
        let mut block: P = QtPoly::zero();
        for a in 0..4 {
            for b in 0..4 {
                block.add_term(a, b, (a as i64) + (b as i64) - 3);
            }
        }
        cases.push(block);
        // A sparse one with gaps the shift steps across.
        let mut sparse: P = QtPoly::zero();
        sparse.add_term(0, 0, 1);
        sparse.add_term(5, 0, -2);
        sparse.add_term(0, 5, 3);
        sparse.add_term(5, 5, 4);
        cases.push(sparse);

        for f in &cases {
            for (a, b) in [(1, 0), (0, 1), (1, 1), (2, 0), (5, 5), (3, 2)] {
                let mut binom: P = <P as Ring>::one();
                binom.add_term(a, b, -1);
                assert_eq!(
                    f.mul_binomial(a, b),
                    f.mul(&binom),
                    "(1 - q^{a} t^{b}) * ({f})"
                );
            }
        }
    }

    /// The cancelling case on its own: `(1 - q)·(1 + q) = 1 - q^2`, where the
    /// middle terms must vanish and leave no stored zero behind.
    #[test]
    fn mul_binomial_drops_cancelled_terms() {
        let mut one_plus_q: P = <P as Ring>::one();
        one_plus_q.add_term(1, 0, 1);
        let got = one_plus_q.mul_binomial(1, 0);
        assert_eq!(got.coeff(0, 0), 1);
        assert_eq!(got.coeff(2, 0), -1);
        assert_eq!(got.len(), 2, "the q terms must cancel: {got}");
    }

    #[test]
    fn ring_axioms_hold() {
        let q: P = QtPoly::q();
        let t: P = QtPoly::t();
        let one = <P as Ring>::one();
        // (q + t)(q - t) = q^2 - t^2
        let sum = q.add_ring(&t);
        let diff = q.add_ring(&t.neg());
        let prod = sum.mul(&diff);
        assert_eq!(prod.coeff(2, 0), 1);
        assert_eq!(prod.coeff(0, 2), -1);
        assert_eq!(prod.len(), 2, "cross terms must cancel: {prod}");
        // distributivity and identity
        assert_eq!(q.mul(&one), q);
        assert_eq!(q.mul(&sum), q.mul(&q).add_ring(&q.mul(&t)));
        // no explicit zeros are ever stored
        assert!(q.add_ring(&q.neg()).is_zero());
    }

    /// Helper: `Ring` has `add_assign` but no `add`, so tests build one.
    trait AddRing: Ring {
        fn add_ring(&self, other: &Self) -> Self {
            let mut x = self.clone();
            x.add_assign(other);
            x
        }
    }
    impl<T: Ring> AddRing for T {}

    #[test]
    fn division_by_integers_is_exact_and_needs_no_field() {
        let x: QtPoly<Rational> = QtPoly::term(2, 3, Rational::from_int(6));
        let half = x.div_u128(4);
        assert_eq!(half.coeff(2, 3), Rational::new(3, 2));
        // and it is genuinely a division, not a truncation
        assert_eq!(half.mul(&<QtPoly<Rational> as Ring>::from_i64(4)), x);
    }

    /// `frobenius` must be a ring homomorphism and the identity at n = 1 —
    /// the contract `Plethystic` states, and what plethysm relies on.
    ///
    /// Over `ℚ[q,t]` rather than `ℤ[q,t]` because `Plethystic: QAlgebra`, and
    /// that is deliberate: plethysm routes through the power-sum basis and so
    /// carries z_μ⁻¹. `ℤ[q,t]` is a perfectly good ring for *holding*
    /// Hall–Littlewood coefficients and cannot support plethysm, which the
    /// bound says out loud.
    #[test]
    fn frobenius_is_a_ring_homomorphism() {
        type R = QtPoly<Rational>;
        let r = Rational::from_int;
        let a: R = QtPoly::term(1, 2, r(3)).add_ring(&QtPoly::term(0, 1, r(-1)));
        let b: R = QtPoly::term(2, 0, r(5)).add_ring(&<R as Ring>::one());
        // psi^0 collapses every monomial to its coefficient, so it is not
        // multiplicative; the Adams operations start at n = 1.
        for n in 1..=4u32 {
            assert_eq!(
                a.mul(&b).frobenius(n),
                a.frobenius(n).mul(&b.frobenius(n)),
                "multiplicative at n = {n}"
            );
            assert_eq!(
                a.add_ring(&b).frobenius(n),
                a.frobenius(n).add_ring(&b.frobenius(n)),
                "additive at n = {n}"
            );
        }
        assert_eq!(a.frobenius(1), a, "identity at n = 1");
        assert_eq!(a.frobenius(3).coeff(3, 6), r(3), "q t^2 -> q^3 t^6");
    }

    /// Sage's convention, checked on the two-variable case that `ℚ[t]` alone
    /// cannot distinguish: `s_2[q·t·s_1] = q²t²·s_2`.
    #[test]
    fn plethysm_raises_both_variables() {
        let qt: QtPoly<Rational> = QtPoly::term(1, 1, Rational::from_int(1));
        let inner: Schur<QtPoly<Rational>> = Schur::monomial(part(&[1]), qt);
        let outer: Schur<QtPoly<Rational>> =
            Schur::monomial(part(&[2]), <QtPoly<Rational> as Ring>::one());
        let r = crate::plethysm::plethysm(&outer, &inner);
        assert_eq!(r.terms().len(), 1);
        assert_eq!(r.coeff(&part(&[2])).coeff(2, 2), Rational::from_int(1));
        assert_eq!(r.coeff(&part(&[2])).len(), 1, "no other (q,t) term");
    }

    /// The case a single monomial cannot distinguish: a **sum** in the
    /// coefficient, where the Frobenius has to act on each term and the answer
    /// mixes. Sage:
    ///
    /// ```text
    ///   sage: R.<q,t> = QQ[]; s = SymmetricFunctions(R).schur()
    ///   sage: s[2]((q+t)*s[1])
    ///   q*t*s[1, 1] + (q^2+q*t+t^2)*s[2]
    /// ```
    ///
    /// A frobenius that scaled only the leading term, or that scaled the whole
    /// polynomial by q^n t^n rather than raising each variable, still passes
    /// `plethysm_raises_both_variables` and fails here.
    #[test]
    fn plethysm_of_a_polynomial_coefficient_matches_sage() {
        type R = QtPoly<Rational>;
        let r = Rational::from_int;
        let q_plus_t: R = QtPoly::term(1, 0, r(1)).add_ring(&QtPoly::term(0, 1, r(1)));
        let inner: Schur<R> = Schur::monomial(part(&[1]), q_plus_t);
        let outer: Schur<R> = Schur::monomial(part(&[2]), <R as Ring>::one());
        let res = crate::plethysm::plethysm(&outer, &inner);

        assert_eq!(res.terms().len(), 2);
        let c11 = res.coeff(&part(&[1, 1]));
        assert_eq!(c11.len(), 1);
        assert_eq!(c11.coeff(1, 1), r(1), "q*t on s_11");

        let c2 = res.coeff(&part(&[2]));
        assert_eq!(c2.len(), 3, "q^2 + q t + t^2, got {c2}");
        assert_eq!(c2.coeff(2, 0), r(1));
        assert_eq!(c2.coeff(1, 1), r(1));
        assert_eq!(c2.coeff(0, 2), r(1));
    }

    /// The library's dividing paths must work over `ℚ[q,t]`, which is the point
    /// of the type. `s → p` carries z_μ⁻¹ and would need a `Field` if the bound
    /// had not been fixed.
    #[test]
    fn s_to_p_round_trips_over_qt() {
        for n in 0..=6u32 {
            for lambda in crate::partitions_of(n) {
                let c: QtPoly<Rational> = QtPoly::term(2, 1, Rational::from_int(3));
                let f: Schur<QtPoly<Rational>> = Schur::monomial(lambda.clone(), c);
                let p: PowerSum<QtPoly<Rational>> = PowerSum::from_schur(&f);
                assert_eq!(p.to_schur(), f, "s -> p -> s at {lambda}");
            }
        }
    }

    /// Evaluation is what makes the specialisations checkable: Hall–Littlewood
    /// at t = 0 is Schur and at t = 1 is monomial, so this is the harness those
    /// tests will use.
    #[test]
    fn evaluation_specialises() {
        // 1 + t + t^2, at t = 0 and t = 1
        let f: P = QtPoly::term(0, 0, 1)
            .add_ring(&QtPoly::term(0, 1, 1))
            .add_ring(&QtPoly::term(0, 2, 1));
        assert_eq!(f.eval(&0, &0), 1);
        assert_eq!(f.eval(&0, &1), 3);
        assert_eq!(f.eval(&0, &2), 7);
        // q and t are independent
        let g: P = QtPoly::term(1, 1, 1);
        assert_eq!(g.eval(&5, &7), 35);
        assert_eq!(g.eval(&5, &0), 0);
    }
}
