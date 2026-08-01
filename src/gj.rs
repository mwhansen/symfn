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
//!   result to report, not a bug to debug away — the same posture
//!   [`Side`](crate::dyck::Side) states for the Delta conjecture, where the
//!   rise version is a theorem and the valley version is open.
//!   [`GjTables`] collects them rather than asserting.
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
//! ## The b = 1 check
//!
//! `b = 1` (`α = 2`) is the other known specialization: the double coset algebra
//! of the hyperoctahedral group `H_n` inside `S_2n`.
//! [`double_coset_coefficient`] computes that side by **enumerating matchings**,
//! which is what makes it independent — the obvious route, via zonal spherical
//! functions, would re-use Jack at α = 2 and check nothing.
//!
//! ## What is actually open
//!
//! Most triples fall in cases already proved, so a bulk count of positive
//! coefficients is not evidence for the conjecture. [`matchings_jack_coverage`]
//! separates them; it is a bibliography and the most stale-prone thing in this
//! file.
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
    /// How many primes the modular engine needed, counting the held-back check.
    /// `0` for this engine, which works in ℚ(α) and has no modulus to run out
    /// of. Reported because that engine *grows* its range on demand, and how
    /// often it had to is a measurement rather than a constant.
    pub primes_used: u32,
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
/// (`docs/record/jack.md`).
pub fn gj_connection_tables(n: u32) -> GjTables {
    let mut t = GjTables {
        n,
        c: BTreeMap::new(),
        h: BTreeMap::new(),
        not_polynomial: Vec::new(),
        c_not_integral: Vec::new(),
        negative: Vec::new(),
        peak_bits: 0,
        primes_used: 0,
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
///
/// # Panics
///
/// If any intermediate leaves `i128` — `n!`, the common denominator, and the
/// character products all grow with `n`, and this route holds nothing back.
/// The wall is unmeasured; `docs/record/jack.md` owns it.
///
/// Off-degree inputs are an *answer*, not a panic: `a^λ_{μν} = 0` unless
/// `|λ| = |μ| = |ν|`, and a caller sweeping a range depends on getting it.
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

// ------------------------------------------------------------- b = 1 --------

/// The coset type of a pair of perfect matchings on `2n` points.
///
/// The union of two perfect matchings is a disjoint set of even cycles that
/// alternate between them; halving the cycle lengths gives a partition of `n`.
/// Matchings are `m[i] = partner of i`.
fn coset_type(a: &[usize], b: &[usize]) -> Partition {
    let mut seen = vec![false; a.len()];
    let mut parts = Vec::new();
    for start in 0..a.len() {
        if seen[start] {
            continue;
        }
        // Alternate an `a` edge and a `b` edge; one such step consumes one of
        // each, so a `k`-step cycle is a part of size `k`.
        let mut j = start;
        let mut k = 0u32;
        loop {
            seen[j] = true;
            seen[a[j]] = true;
            j = b[a[j]];
            k += 1;
            if j == start {
                break;
            }
        }
        parts.push(k);
    }
    Partition::new(parts)
}

/// A matching whose coset type against `pairs_in_order` is `lambda`.
///
/// Within each block of `2λᵢ` consecutive points, shift by one: the reference
/// pairs `(0,1),(2,3),…` and this pairs `(1,2),(3,4),…,(2λᵢ−1,0)`, so the union
/// of the block is a single cycle of length `2λᵢ`.
fn matching_of_type(lambda: &Partition) -> Vec<usize> {
    let n = lambda.size() as usize;
    let mut m = vec![0usize; 2 * n];
    let mut base = 0usize;
    for &part in lambda.parts() {
        let len = 2 * part as usize;
        for s in (1..len).step_by(2) {
            m[base + s] = base + (s + 1) % len;
            m[base + (s + 1) % len] = base + s;
        }
        base += len;
    }
    m
}

/// `(0,1),(2,3),…` — the reference matching.
fn pairs_in_order(n: usize) -> Vec<usize> {
    let mut m = vec![0usize; 2 * n];
    for i in 0..n {
        m[2 * i] = 2 * i + 1;
        m[2 * i + 1] = 2 * i;
    }
    m
}

/// Call `f` on every perfect matching of `2n` points.
fn for_each_matching(n: usize, mut f: impl FnMut(&[usize])) {
    let mut m = vec![usize::MAX; 2 * n];
    fn go(m: &mut Vec<usize>, f: &mut impl FnMut(&[usize])) {
        match m.iter().position(|&v| v == usize::MAX) {
            None => f(m),
            Some(i) => {
                for j in i + 1..m.len() {
                    if m[j] == usize::MAX {
                        m[i] = j;
                        m[j] = i;
                        go(m, f);
                        m[i] = usize::MAX;
                        m[j] = usize::MAX;
                    }
                }
            }
        }
    }
    go(&mut m, &mut f);
}

/// `b^λ_{μν}`, the double-coset connection coefficient of `(S_2n, H_n)`, by
/// counting matchings.
///
/// ```text
///   b^λ_{μν} = #{ δ : type(δ₀,δ) = μ and type(δ,δ₁) = ν },
///   where δ₀, δ₁ are fixed with type(δ₀,δ₁) = λ
/// ```
///
/// The `b = 1` analogue of [`class_algebra_coefficient`]: [GJ] specialize their
/// series to the double coset algebra of the hyperoctahedral group at `b = 1`,
/// exactly as `b = 0` gives the class algebra of `S_n`. This computes the
/// right-hand side by **enumerating the `(2n−1)!! matchings directly** — no Jack
/// polynomial, no zonal polynomial, no character. That independence is the whole
/// point: computing it from zonal spherical functions would re-use Jack at
/// α = 2 and check nothing.
///
/// Cost is `(2n−1)!!`, so 105 at n = 4 and 2,027,025 at n = 8. Fine as a pin at
/// small degree and hopeless as an engine, which is the usual shape for these.
pub fn double_coset_coefficient(la: &Partition, mu: &Partition, nu: &Partition) -> u64 {
    let n = la.size() as usize;
    if mu.size() as usize != n || nu.size() as usize != n {
        return 0;
    }
    let d0 = pairs_in_order(n);
    let d1 = matching_of_type(la);
    debug_assert_eq!(&coset_type(&d0, &d1), la, "the witness has the wrong type");
    let mut count = 0u64;
    for_each_matching(n, |d| {
        if &coset_type(&d0, d) == mu && &coset_type(d, &d1) == nu {
            count += 1;
        }
    });
    count
}

/// Every `b^λ_{μν}` at degree `n` at once, keyed as [`Key`]. Zeros omitted.
///
/// One sweep of the `(2n−1)!!` matchings per λ rather than per triple, which is
/// a factor of `p(n)²` — 1331 at n = 6 — and the difference between a pin that
/// runs in a unit test and one that does not.
pub fn double_coset_table(n: u32) -> BTreeMap<Key, u64> {
    let mut out: BTreeMap<Key, u64> = BTreeMap::new();
    if n == 0 {
        return out;
    }
    let d0 = pairs_in_order(n as usize);
    for la in crate::partitions_of(n) {
        let d1 = matching_of_type(&la);
        for_each_matching(n as usize, |d| {
            let mu = coset_type(&d0, d);
            let nu = coset_type(d, &d1);
            *out.entry((la.clone(), mu, nu)).or_insert(0) += 1;
        });
    }
    out
}

// ----------------------------------------------------- what is still open ---

/// Which triples the Matchings–Jack conjecture is still open on.
///
/// ⚠️ **This is a bibliography, and it is the part of this file most likely to
/// go stale.** It exists because a bulk count of positive coefficients mostly
/// counts already-proved cases, and reporting "2,045,553 terms verified" without
/// it overstates what the computation shows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Coverage {
    /// `λ = [1ⁿ]` or `λ = [2,1^{n−2}]` — Goulden and Jackson constructed the
    /// statistic and proved these two in the original paper.
    ProvedByGj,
    /// One of the three partitions is `(n)`.
    ///
    /// ⚠️ Kanunnikov–Vassilieva proved `μ = ν = (n)` outright. The extension to
    /// *any one* of the three being `(n)`, with Promyslov, proves **a variation
    /// involving additional labellings on matchings** — not the original
    /// statement. So this is weaker than [`Coverage::ProvedByGj`] and is
    /// reported separately rather than merged into it.
    SinglePartVariant,
    /// Neither — the conjecture is open here.
    ///
    /// This is where a computation is worth anything. The smallest such triple
    /// is `λ = μ = ν = (2,2)` at n = 4.
    Open,
}

/// Classify a triple against the literature; see [`Coverage`].
///
/// Independent of anything computed here, and deliberately conservative: a case
/// is only called covered when a cited result covers it.
pub fn matchings_jack_coverage(la: &Partition, mu: &Partition, nu: &Partition) -> Coverage {
    let n = la.size();
    let is_single = |p: &Partition| p.parts() == [n];
    let all_ones = la.parts().iter().all(|&p| p == 1);
    let hook21 = {
        let p = la.parts();
        p.first() == Some(&2) && p[1..].iter().all(|&v| v == 1)
    };
    if all_ones || hook21 {
        Coverage::ProvedByGj
    } else if is_single(la) || is_single(mu) || is_single(nu) {
        Coverage::SinglePartVariant
    } else {
        Coverage::Open
    }
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
                 normalization against docs/record/jack.md and REPORT \
                 it, do not 'fix' it: {:?}",
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

    /// **`b = 1` is the double coset algebra of `(S_2n, H_n)`** — the second pin,
    /// alongside `b = 0` and the class algebra.
    ///
    /// The right-hand side counts matchings and nothing else: no Jack
    /// polynomial, no zonal polynomial, no character. Computing it from zonal
    /// spherical functions instead would re-use Jack at α = 2 and check nothing.
    ///
    /// ⚠️ The **normalization was measured, not read from [GJ]**. The ratio came
    /// back exactly 1 on all 285 live triples through n = 5 — no factor of
    /// `z_λ`, `2^{ℓ}`, or anything else — and this test then requires it through
    /// n = 6, which is 484 triples the constant was not fitted on.
    #[test]
    fn the_b_one_slice_is_the_double_coset_algebra() {
        for n in 1..=6u32 {
            let t = gj_connection_tables(n);
            let want = double_coset_table(n);
            let parts = crate::partitions_of(n);
            let mut checked = 0;
            for la in &parts {
                for mu in &parts {
                    for nu in &parts {
                        let key = (la.clone(), mu.clone(), nu.clone());
                        let got =
                            t.c.get(&key)
                                .map_or((0i128, 1u128), |p| (p.num.iter().sum::<i128>(), p.den));
                        let b = *want.get(&key).unwrap_or(&0) as i128;
                        assert_eq!(got.1, 1, "c^{la}_{{{mu},{nu}}}(1) is not an integer");
                        assert_eq!(got.0, b, "c^{la}_{{{mu},{nu}}}(1)");
                        checked += 1;
                    }
                }
            }
            assert!(checked > 0);
        }
    }

    /// The matchings machinery itself, independent of anything [GJ].
    #[test]
    fn coset_types_are_what_they_should_be() {
        // Against itself: the union is n cycles of length 2, so type [1^n].
        for n in 1..=5usize {
            let d0 = pairs_in_order(n);
            assert_eq!(coset_type(&d0, &d0), Partition::new(vec![1; n]));
            // And every witness has the type it advertises.
            for la in crate::partitions_of(n as u32) {
                assert_eq!(coset_type(&d0, &matching_of_type(&la)), la);
            }
            // (2n-1)!! matchings, all distinct.
            let mut count = 0usize;
            for_each_matching(n, |_| count += 1);
            let want: usize = (1..=n).map(|k| 2 * k - 1).product();
            assert_eq!(count, want, "matchings of 2*{n} points");
        }
    }

    /// The coverage classifier, on the cases that matter.
    #[test]
    fn the_coverage_classifier_is_conservative() {
        let p = |v: Vec<u32>| Partition::new(v);
        // Goulden-Jackson's own two lambda cases.
        let ones = p(vec![1, 1, 1, 1]);
        let hook = p(vec![2, 1, 1]);
        assert_eq!(
            matchings_jack_coverage(&ones, &p(vec![2, 2]), &p(vec![2, 2])),
            Coverage::ProvedByGj
        );
        assert_eq!(
            matchings_jack_coverage(&hook, &p(vec![2, 2]), &p(vec![2, 2])),
            Coverage::ProvedByGj
        );
        // Any one of the three equal to (n) -- the labelled variant.
        assert_eq!(
            matchings_jack_coverage(&p(vec![2, 2]), &p(vec![4]), &p(vec![2, 2])),
            Coverage::SinglePartVariant
        );
        // The smallest genuinely open triple.
        assert_eq!(
            matchings_jack_coverage(&p(vec![2, 2]), &p(vec![2, 2]), &p(vec![2, 2])),
            Coverage::Open
        );
        // And nothing at n < 4 is open, which is why n = 4 is the smallest.
        for n in 1..=3u32 {
            for la in crate::partitions_of(n) {
                for mu in crate::partitions_of(n) {
                    for nu in crate::partitions_of(n) {
                        assert_ne!(
                            matchings_jack_coverage(&la, &mu, &nu),
                            Coverage::Open,
                            "{la} {mu} {nu}"
                        );
                    }
                }
            }
        }
    }
}
