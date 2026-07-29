//! The Goulden–Jackson connection-coefficient pipeline: `c^λ_{μν}(b)` and
//! `h^λ_{μν}(b)`.
//!
//! Two open conjectures live on these two families, and no package computes
//! either table. [GJ] TAMS 348 (1996) 873–892 defines them from a single
//! Jack-Cauchy kernel, with `α = 1 + b`:
//!
//! ```text
//!   Φ(x,y,z;t,α) = Σ_θ  t^{|θ|} · J_θ(x)J_θ(y)J_θ(z) / ⟨J_θ,J_θ⟩_α        (1)
//!   Ψ            = α · t ∂/∂t log Φ                                        (2)
//!   c^λ_{μν}(b)  = z_λ(1+b)^{ℓ(λ)} · [t^n p_λ(x)p_μ(y)p_ν(z)] Φ            (5)
//!   h^λ_{μν}(b)  =                    [t^n p_λ(x)p_μ(y)p_ν(z)] Ψ           (4)
//! ```
//!
//! **Matchings-Jack conjecture**: `c^λ_{μν}(b) = Σ_δ b^{wt_λ(δ)}` over
//! matchings — in particular it lies in ℕ[b]. **b-conjecture**:
//! `h^λ_{μν}(b) = Σ_M b^{ϑ(M)}` over rooted connected bipartite maps in
//! locally orientable surfaces. Status: ℚ[b]-polynomiality **proven** for both
//! ([DF], [arXiv:1601.01501](https://arxiv.org/abs/1601.01501)); integrality of
//! `c` **proven** ([BD], [arXiv:2203.14879](https://arxiv.org/abs/2203.14879));
//! **positivity open for both**.
//!
//! So this module computes two things it must not "fix":
//!
//! - Polynomiality and `c`'s integrality are **theorems**, so a failure there
//!   is *our* bug and is reported as one.
//! - Positivity is **the open question itself**. A negative coefficient is a
//!   result to report, not a bug to debug away — the valley-Delta-conjecture
//!   posture of `spec-macdonald-operators.md` §5.9, verbatim. [`GjTables`]
//!   collects them rather than asserting.
//!
//! ## How it is computed
//!
//! `Φ_n := [t^n]Φ` is an element of `Sym^{⊗3}` in the power-sum basis, so it is
//! a map from triples of partitions to scalars. Each `J_θ` arrives in the p
//! basis from [`jack_j_powersum`](crate::jack_j_powersum), and `⟨J_θ,J_θ⟩` is
//! the closed product `H_θH'_θ` — never a pairing.
//!
//! The `log` is the standard exponential recurrence. From `Φ = exp(L)` and
//! `Φ' = L'Φ` in `t`, writing `G_k := k·L_k`,
//!
//! ```text
//!   G_n = n·Φ_n − Σ_{k=1}^{n−1} G_k · Φ_{n−k},      Ψ_n = α · G_n
//! ```
//!
//! and every product there is a **multiset union in each of the three
//! factors** — the power-sum basis is multiplicative, so no basis change of any
//! product is ever formed.
//!
//! ## The b = 0 check, and why it is a real one
//!
//! At `b = 0` (`α = 1`) the whole pipeline must collapse onto the class algebra
//! of `S_n`: `J_θ = h_θ s_θ` and `⟨J_θ,J_θ⟩ = h_θ²`, so
//!
//! ```text
//!   c^λ_{μν}(0) = z_λ · [p_λp_μp_ν] Σ_θ h_θ s_θ^{⊗3}
//!               = (n!/(z_μ z_ν)) · Σ_θ χ^θ_λ χ^θ_μ χ^θ_ν / f^θ
//!               = a^λ_{μν},
//! ```
//!
//! the connection coefficient `C_μ C_ν = Σ_λ a^λ_{μν} C_λ`. That is computable
//! from [`character`](crate::character) alone — no Jack polynomial, no `AFrac`,
//! no fraction — so agreement pins the normalization of (1) and (5) against an
//! object with an independent definition. It is the check that says the
//! transcription of [GJ]'s formulas is right.
//!
//! ## Fixed width
//!
//! Concrete over `i128`, not generic: positivity needs an order, and the
//! `b = 0` check needs `n!`-sized integers to be exact. `i128` holds
//! comfortably at the degrees this reaches — [`GjTables::peak_bits`] reports
//! how close it gets, and the `b = 0` check would fail loudly on a wrap.

use std::collections::{BTreeMap, HashMap};

use crate::afrac::AFrac;
use crate::coeff::Ring;
use crate::jack::{jack_j_powersum, jack_norm_j};
use crate::partition::Partition;
use crate::sym::SymFn;

/// A triple of partitions indexing `p_λ ⊗ p_μ ⊗ p_ν`.
pub type Key = (Partition, Partition, Partition);

/// An element of `Sym^{⊗3}` in the power-sum basis.
type Tensor = HashMap<Key, AFrac<i128>>;

/// A polynomial in `b` with a common integer denominator:
/// `(Σ num[k]·b^k) / den`.
///
/// [DF] says both families land here; [BD] says `c`'s `den` is 1. Keeping the
/// denominator rather than demanding it be 1 is what lets the weaker theorem be
/// checked separately from the stronger one.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct BPoly {
    pub num: Vec<i128>,
    pub den: u128,
}

impl BPoly {
    /// Whether every coefficient is ≥ 0 — the open half of both conjectures.
    pub fn is_positive(&self) -> bool {
        self.num.iter().all(|&c| c >= 0)
    }
    /// Whether this is an integer polynomial — proven for `c` [BD], open for
    /// `h`.
    pub fn is_integral(&self) -> bool {
        self.den == 1
    }
    /// The value at `b = 0`, as an exact ratio.
    pub fn at_zero(&self) -> (i128, u128) {
        (self.num.first().copied().unwrap_or(0), self.den)
    }
}

/// One degree of both [GJ] tables, with the collapse checks already run.
#[derive(Clone, Debug)]
pub struct GjTables {
    pub n: u32,
    /// `c^λ_{μν}(b)`, keyed `(λ, μ, ν)`. Zero entries are omitted.
    pub c: BTreeMap<Key, BPoly>,
    /// `h^λ_{μν}(b)`, keyed `(λ, μ, ν)`.
    pub h: BTreeMap<Key, BPoly>,
    /// Keys where the ℚ(α) value refused to collapse to a polynomial in `b`.
    /// **[DF] says this must be empty; a nonempty one is our bug.**
    pub not_polynomial: Vec<(&'static str, Key)>,
    /// `c` keys with a surviving denominator. **[BD] says this must be empty.**
    pub c_not_integral: Vec<Key>,
    /// Keys with a negative coefficient. **Open — a result, not a bug.**
    pub negative: Vec<(&'static str, Key, BPoly)>,
    /// Widest numerator coefficient seen anywhere, in bits — the `i128` margin.
    pub peak_bits: u32,
}

impl GjTables {
    /// The two proven laws, as one question: did anything that must not happen
    /// happen?
    pub fn laws_hold(&self) -> bool {
        self.not_polynomial.is_empty() && self.c_not_integral.is_empty()
    }
}

/// Multiply two `Sym^{⊗3}` elements: a multiset union in each factor.
fn tensor_mul(a: &Tensor, b: &Tensor) -> Tensor {
    let mut out: Tensor = HashMap::new();
    for ((l1, m1, n1), ca) in a {
        for ((l2, m2, n2), cb) in b {
            let key = (union(l1, l2), union(m1, m2), union(n1, n2));
            let v = ca.mul(cb);
            if v.is_zero() {
                continue;
            }
            match out.entry(key) {
                std::collections::hash_map::Entry::Occupied(mut e) => e.get_mut().add_assign(&v),
                std::collections::hash_map::Entry::Vacant(e) => {
                    e.insert(v);
                }
            }
        }
    }
    out.retain(|_, v| !v.is_zero());
    out
}

fn union(a: &Partition, b: &Partition) -> Partition {
    let mut parts = a.parts().to_vec();
    parts.extend_from_slice(b.parts());
    Partition::new(parts)
}

/// `[t^n]Φ` — the degree-`n` slice of the Jack-Cauchy kernel, in `p^{⊗3}`.
fn phi_slice(n: u32) -> Tensor {
    let mut out: Tensor = HashMap::new();
    if n == 0 {
        out.insert(
            (
                Partition::default(),
                Partition::default(),
                Partition::default(),
            ),
            <AFrac<i128> as Ring>::one(),
        );
        return out;
    }
    for theta in crate::partitions_of(n) {
        let j = jack_j_powersum::<i128>(&theta);
        // Fold 1/⟨J_θ,J_θ⟩ into one of the three factors, so the triple loop
        // does two multiplies instead of four.
        let inverse: BTreeMap<(u32, u32), i32> = jack_norm_j(&theta)
            .into_iter()
            .map(|(k, m)| (k, -m))
            .collect();
        let scaled: Vec<(&Partition, AFrac<i128>)> = j
            .terms()
            .iter()
            .map(|(rho, c)| {
                let mut v = c.mul_factors(&inverse);
                v.reduce();
                (rho, v)
            })
            .collect();
        for (la, ca) in &scaled {
            for (mu, cb) in j.terms() {
                let ab = ca.mul(cb);
                if ab.is_zero() {
                    continue;
                }
                for (nu, cc) in j.terms() {
                    let v = ab.mul(cc);
                    if v.is_zero() {
                        continue;
                    }
                    let key = ((*la).clone(), mu.clone(), nu.clone());
                    match out.entry(key) {
                        std::collections::hash_map::Entry::Occupied(mut e) => {
                            e.get_mut().add_assign(&v)
                        }
                        std::collections::hash_map::Entry::Vacant(e) => {
                            e.insert(v);
                        }
                    }
                }
            }
        }
    }
    for v in out.values_mut() {
        v.reduce();
    }
    out.retain(|_, v| !v.is_zero());
    out
}

/// Substitute `α = 1 + b` into a dense polynomial in α — a binomial transform.
fn shift_to_b(coeffs: &[i128]) -> Vec<i128> {
    let mut out = vec![0i128; coeffs.len()];
    for (k, &c) in coeffs.iter().enumerate() {
        if c == 0 {
            continue;
        }
        // (1 + b)^k = Σ_j C(k,j) b^j
        let mut binom: i128 = 1;
        for j in 0..=k {
            out[j] += c * binom;
            binom = binom * (k - j) as i128 / (j + 1) as i128;
        }
    }
    while out.last() == Some(&0) {
        out.pop();
    }
    out
}

/// Collapse one ℚ(α) value to a polynomial in `b`, recording what it failed.
fn collapse(v: AFrac<i128>, tag: &'static str, key: &Key, t: &mut GjTables) -> Option<BPoly> {
    let Some((coeffs, den)) = v.into_rational_poly() else {
        t.not_polynomial.push((tag, key.clone()));
        return None;
    };
    for &c in &coeffs {
        t.peak_bits = t.peak_bits.max(128 - c.unsigned_abs().leading_zeros());
    }
    let poly = BPoly {
        num: shift_to_b(&coeffs),
        den,
    };
    if poly.num.is_empty() {
        return None;
    }
    if tag == "c" && !poly.is_integral() {
        t.c_not_integral.push(key.clone());
    }
    if !poly.is_positive() {
        t.negative.push((tag, key.clone(), poly.clone()));
    }
    Some(poly)
}

/// Both [GJ] tables at degree `n`, with the [DF] / [BD] collapse checks run.
///
/// The unit of work is the whole degree, because `Ψ` needs every lower `Φ_k`.
/// Sage's unit of work for the same pipeline — a single `J → p` at n = 12 —
/// already costs 303.8 s before the triple product starts
/// (`docs/spec-jack.md` §2.2).
pub fn gj_connection_tables(n: u32) -> GjTables {
    let mut t = GjTables {
        n,
        c: BTreeMap::new(),
        h: BTreeMap::new(),
        not_polynomial: Vec::new(),
        c_not_integral: Vec::new(),
        negative: Vec::new(),
        peak_bits: 0,
    };
    if n == 0 {
        return t;
    }

    let phi: Vec<Tensor> = (0..=n).map(phi_slice).collect();

    // G_k = k·Φ_k − Σ_{j<k} G_j·Φ_{k−j}, so that Ψ_k = α·G_k.
    let mut g: Vec<Tensor> = vec![HashMap::new()];
    for k in 1..=n as usize {
        let mut acc: Tensor = phi[k]
            .iter()
            .map(|(key, v)| (key.clone(), v.scale_int(k as i64)))
            .collect();
        for (j, gj) in g.iter().enumerate().skip(1) {
            for (key, v) in tensor_mul(gj, &phi[k - j]) {
                match acc.entry(key) {
                    std::collections::hash_map::Entry::Occupied(mut e) => {
                        e.get_mut().add_assign(&v.neg())
                    }
                    std::collections::hash_map::Entry::Vacant(e) => {
                        e.insert(v.neg());
                    }
                }
            }
        }
        for v in acc.values_mut() {
            v.reduce();
        }
        acc.retain(|_, v| !v.is_zero());
        g.push(acc);
    }

    // c: multiply by z_λ·α^{ℓ(λ)}; h: multiply by α. Both then α = 1 + b.
    let alpha = |power: i32| -> BTreeMap<(u32, u32), i32> {
        let mut f = BTreeMap::new();
        if power != 0 {
            f.insert((1, 0), power);
        }
        f
    };
    for (key, v) in &phi[n as usize] {
        let scaled = v
            .mul_factors(&alpha(key.0.len() as i32))
            .mul(&AFrac::from_u128(key.0.z()));
        if let Some(poly) = collapse(scaled, "c", key, &mut t) {
            t.c.insert(key.clone(), poly);
        }
    }
    for (key, v) in &g[n as usize] {
        let scaled = v.mul_factors(&alpha(1));
        if let Some(poly) = collapse(scaled, "h", key, &mut t) {
            t.h.insert(key.clone(), poly);
        }
    }
    t
}

/// `a^λ_{μν}`, the class-algebra connection coefficient of `S_n`, from
/// characters alone.
///
/// ```text
///   C_μ · C_ν = Σ_λ a^λ_{μν} C_λ,
///   a^λ_{μν}  = (n!/(z_μ z_ν)) · Σ_θ χ^θ_λ χ^θ_μ χ^θ_ν / f^θ
/// ```
///
/// The independent definition [`gj_connection_tables`] is pinned against at
/// `b = 0`. Nothing here touches a Jack polynomial, an [`AFrac`], or a
/// fraction field, so agreement is evidence about the transcription of [GJ]'s
/// (1) and (5) and not a restatement of it.
pub fn class_algebra_coefficient(la: &Partition, mu: &Partition, nu: &Partition) -> i128 {
    let n = la.size();
    if mu.size() != n || nu.size() != n {
        return 0;
    }
    let factorial: i128 = (1..=i128::from(n)).product::<i128>().max(1);
    // Σ_θ χχχ/f^θ over a common denominator, so the sum stays in ℤ.
    let thetas = crate::partitions_of(n);
    let dims: Vec<i128> = thetas
        .iter()
        .map(|th| crate::dimension(th).expect("a partition has a dimension") as i128)
        .collect();
    let mut num = 0i128;
    let mut den = 1i128;
    for (th, &f) in thetas.iter().zip(dims.iter()) {
        let term = crate::character(th, la) * crate::character(th, mu) * crate::character(th, nu);
        // num/den += term/f
        let g = gcd(den, f);
        num = num * (f / g) + term * (den / g);
        den = den / g * f;
        let g = gcd(num.abs(), den);
        if g > 1 {
            num /= g;
            den /= g;
        }
    }
    // Over ONE denominator. `n!/(z_μ z_ν)` need not be an integer on its own —
    // it is 3/2 already for two transpositions in S_3 — so dividing in stages
    // truncates.
    let whole = num * factorial;
    let denom = den * (mu.z() as i128) * (nu.z() as i128);
    debug_assert_eq!(
        whole % denom,
        0,
        "the connection coefficient must be an integer"
    );
    whole / denom
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a.max(1)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    /// The class-algebra coefficients themselves, on cases that can be counted
    /// by hand in `S_3`.
    ///
    /// `C_{(2,1)}` is the three transpositions. `C_{(2,1)}² = 3·C_{(1,1,1)} +
    /// 3·C_{(3)}`: a product of two transpositions is the identity in 3 ways
    /// and a 3-cycle in 6, spread over the 2 elements of `C_{(3)}`.
    #[test]
    fn class_algebra_on_s3() {
        let t = part(&[2, 1]);
        assert_eq!(class_algebra_coefficient(&part(&[1, 1, 1]), &t, &t), 3);
        assert_eq!(class_algebra_coefficient(&part(&[3]), &t, &t), 3);
        assert_eq!(class_algebra_coefficient(&part(&[2, 1]), &t, &t), 0);
    }

    /// **The transcription check.** `c^λ_{μν}(0)` must be the class-algebra
    /// connection coefficient — computed here from characters, with no Jack
    /// polynomial anywhere in it.
    ///
    /// This is what says [GJ]'s (1) and (5) went in correctly, including the
    /// `z_λ(1+b)^{ℓ(λ)}` prefactor, which is the part with no independent
    /// derivation available.
    #[test]
    fn the_c_table_at_b_zero_is_the_class_algebra() {
        for n in 1..=6u32 {
            let t = gj_connection_tables(n);
            assert!(t.laws_hold(), "a proven law failed at n = {n}: {t:?}");
            let parts = crate::partitions_of(n);
            for la in &parts {
                for mu in &parts {
                    for nu in &parts {
                        let key = (la.clone(), mu.clone(), nu.clone());
                        let got = t.c.get(&key).map_or((0, 1), BPoly::at_zero);
                        let want = class_algebra_coefficient(la, mu, nu);
                        assert_eq!(got.1, 1, "c must be integral at {key:?}");
                        assert_eq!(
                            got.0, want,
                            "c^{la}_{{{mu},{nu}}}(0) = {} but a^λ_μν = {want}",
                            got.0
                        );
                    }
                }
            }
        }
    }

    /// [DF]: both families are polynomials in `b`; [BD]: `c` is integral.
    /// Both are theorems, so both are assertions.
    #[test]
    fn the_proven_laws_hold() {
        for n in 1..=6u32 {
            let t = gj_connection_tables(n);
            assert!(
                t.not_polynomial.is_empty(),
                "n = {n}: [DF] says these collapse: {:?}",
                t.not_polynomial
            );
            assert!(
                t.c_not_integral.is_empty(),
                "n = {n}: [BD] says c is integral, but {:?}",
                t.c_not_integral
            );
            assert!(
                !t.c.is_empty() && !t.h.is_empty(),
                "n = {n} produced nothing"
            );
        }
    }

    /// **Both conjectures, observed.** Positivity is the *open* question, so a
    /// failure here is a result to report and not a bug to fix — the
    /// valley-Delta posture.
    #[test]
    fn positivity_is_observed_not_assumed() {
        for n in 1..=6u32 {
            let t = gj_connection_tables(n);
            assert!(
                t.negative.is_empty(),
                "n = {n}: a NEGATIVE coefficient — this is a RESULT, verify the \
                 normalization against docs/spec-jack.md §3.6 and REPORT it, \
                 do not 'fix' it: {:?}",
                t.negative
            );
        }
    }

    /// `h` is the connected version of `c`, so at `n = 1` they coincide, and
    /// `h`'s support is contained in `c`'s.
    #[test]
    fn h_is_supported_where_c_is() {
        for n in 1..=5u32 {
            let t = gj_connection_tables(n);
            for key in t.h.keys() {
                assert!(
                    t.c.contains_key(key),
                    "n = {n}: h has {key:?} where c does not"
                );
            }
        }
    }

    /// The α → 1+b substitution, on cases with an obvious answer.
    #[test]
    fn the_binomial_shift() {
        assert_eq!(shift_to_b(&[5]), vec![5]);
        assert_eq!(shift_to_b(&[0, 1]), vec![1, 1], "α = 1 + b");
        assert_eq!(shift_to_b(&[0, 0, 1]), vec![1, 2, 1], "α² = 1 + 2b + b²");
        assert_eq!(shift_to_b(&[0, -1, 1]), vec![0, 1, 1], "α² − α = b + b²");
    }
}
