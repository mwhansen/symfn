//! The Goulden–Jackson tables again, by modular evaluation and interpolation.
//!
//! No package computes either table, Sage included
//! (`docs/record/oracles-and-comparisons.md`).
//! A second engine for
//! [`gj_connection_tables`](crate::gj::gj_connection_tables), sharing no
//! arithmetic with it — the `qtkostka.rs` "three routes" standard. The exact
//! one works in ℚ(α) throughout; this one runs the entire pipeline at numeric α
//! over a prime field and reconstructs the answers.
//!
//! ## Why, and why this shape
//!
//! `Φ_n := [t^n]Φ` is `Σ_θ` of a rank-1 tensor over `p(n)³` entries, so
//! `p(n)⁴` coefficient operations. `examples/probe_gj.rs` measured the four
//! candidate fixes and killed three of them:
//!
//! - `J → p` carries **no** atoms — every one comes from a single
//!   `1/⟨J_θ,J_θ⟩` per θ — so hoisting them out of the inner loop looks
//!   obvious. It is worth nothing, because the numerator polynomial dominates
//!   the multiply and atom-carrying values have *smaller* numerators.
//! - One global denominator `D = lcm_θ⟨J_θ,J_θ⟩` would remove the atoms
//!   entirely, but `deg D` is 73 at n = 10 against `2n = 20` and grows like
//!   `n²`. It trades atom bookkeeping for degree blowup.
//! - Evaluating at numeric α over ℚ is cheaper per operation, but not by
//!   enough to pay for the `~n+2` points an interpolation needs.
//!
//! What survives is that last idea with a scalar that is actually cheap:
//! `Rational` runs a 128-bit gcd per operation and a prime field does not
//! (`docs/record/jack.md`).
//!
//! ## Two things make it safe rather than a guess
//!
//! **There are no poles for α > 0.** Every atom is `uα + v` with `u, v ≥ 0` and
//! not both zero, so every positive α is a legal evaluation point. That is a
//! proof about the shape of the hooks, not a sampling argument, and
//! `probe_gj` asserts it over every atom the engine produces.
//!
//! **The object interpolated is the answer, not an intermediate.** `Φ` is not a
//! polynomial in α at all; `c` and `h` are polynomials in `b` ([DF]). So the
//! whole pipeline — the `Φ_k` slices, the `G_k = kΦ_k − Σ G_jΦ_{k−j}`
//! recurrence and the `z_λα^{ℓ(λ)}` scaling — runs in residues, and only `c`
//! and `h` are reconstructed.
//!
//! ## What it gives up, and what replaces it
//!
//! The exact engine enforces [DF] for free: a value that failed to collapse out
//! of ℚ(α) would be caught by `into_rational_poly` returning `None`. Working in
//! residues there is nothing to collapse, so that check is gone. Two things
//! replace it, and both are stronger than they look:
//!
//! - the **degree bound is verified, not assumed** — interpolation runs with
//!   spare points and requires the top coefficients to vanish, which is exactly
//!   the statement that the answer is a polynomial of the expected degree;
//! - the **`b = 0` slice must be the class algebra of `S_n`**, computed from
//!   characters alone, which is now the law that has to catch an error;
//! - an **independent prime**, held back from the reconstruction and required
//!    to agree with it, because rational reconstruction returns a spurious
//!    small rational rather than failing when the true value is out of range.
//!
//! Plus the exact engine itself, at every degree it can still reach.
//!
//! ## The prime size, which is where the cost went
//!
//! Under `2^31` a product fits a `u64`, and that is why the bound is there:
//! `u128 %` is a function call on aarch64 and `u64 %` is not. See
//! `modular::Md` for the arithmetic and
//! `docs/record/jack.md` for the ladder it was measured on.
//!
//! [DF]: https://arxiv.org/abs/1601.01501

// The two wide casts carry checks at their sites; the rest are shape indices.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use std::collections::BTreeMap;

use crate::gj::{BPoly, GjTables, Key};
use crate::jack::{jack_j_powersum, jack_norm_j};
use crate::modular as mp;
use crate::partition::Partition;
use crate::sym::SymFn;

// ------------------------------------------------------ the dense layout ----

/// Everything about degree `k` that does not depend on α or the prime.
struct Degrees {
    /// `parts[k]` — the partitions of `k`, indexed.
    parts: Vec<Vec<Partition>>,
    /// `index[k]` — partition → its position in `parts[k]`.
    index: Vec<BTreeMap<Partition, usize>>,
    /// `union[j][m][a*p(m)+b]` — the index in degree `j+m` of the multiset
    /// union of `parts[j][a]` and `parts[m][b]`.
    ///
    /// The `G_j · Φ_{k−j}` products are multiset unions in each of the three
    /// tensor factors, so this table turns the whole recurrence into index
    /// arithmetic — no partition is ever built inside the loop.
    union: Vec<Vec<Vec<u32>>>,
}

impl Degrees {
    fn new(n: u32) -> Self {
        let parts: Vec<Vec<Partition>> = (0..=n).map(crate::partitions_of).collect();
        let index: Vec<BTreeMap<Partition, usize>> = parts
            .iter()
            .map(|v| v.iter().cloned().zip(0..).collect())
            .collect();
        let mut union = vec![vec![Vec::new(); n as usize + 1]; n as usize + 1];
        for j in 0..=n as usize {
            for m in 0..=n as usize - j {
                let mut t = Vec::with_capacity(parts[j].len() * parts[m].len());
                for a in &parts[j] {
                    for b in &parts[m] {
                        let mut v = a.parts().to_vec();
                        v.extend_from_slice(b.parts());
                        t.push(index[j + m][&Partition::new(v)] as u32);
                    }
                }
                union[j][m] = t;
            }
        }
        Degrees {
            parts,
            index,
            union,
        }
    }

    fn len(&self, k: usize) -> usize {
        self.parts[k].len()
    }
    /// Entries in a `Sym^{⊗3}` slab of degree `k`.
    fn cube(&self, k: usize) -> usize {
        let p = self.len(k);
        p * p * p
    }
}

/// `J_θ` in the p basis and `⟨J_θ,J_θ⟩`, held exactly so they can be evaluated
/// at any α without recomputing the Jack polynomial.
struct Exact {
    /// `jp[k][θ][ρ]` — the numerator coefficients and the integer scalar.
    jp: Vec<Vec<Vec<(Vec<i128>, u128)>>>,
    /// `norm[k][θ]` — the atoms of `H_θH'_θ`.
    norm: Vec<Vec<Vec<((u32, u32), i32)>>>,
}

impl Exact {
    fn new(n: u32, deg: &Degrees) -> Self {
        let mut jp = vec![Vec::new(); n as usize + 1];
        let mut norm = vec![Vec::new(); n as usize + 1];
        for k in 1..=n as usize {
            for theta in &deg.parts[k] {
                let f = jack_j_powersum::<i128>(theta);
                let mut row = vec![(Vec::new(), 1u128); deg.len(k)];
                for (rho, c) in f.terms() {
                    let (num, den, scale) = c.parts();
                    // Measured in `probe_gj`: J → p carries no linear atoms at
                    // any degree. It is asserted rather than assumed, because
                    // the whole evaluation shortcut rests on it.
                    assert_eq!(
                        den.count(),
                        0,
                        "J_{theta} at p_{rho} carries a linear atom; the modular \
                         engine's evaluation of it would be wrong"
                    );
                    row[deg.index[k][rho]] = (num.to_vec(), scale);
                }
                jp[k].push(row);
                norm[k].push(jack_norm_j(theta).into_iter().collect());
            }
        }
        Exact { jp, norm }
    }
}

// --------------------------------------------------------- the pipeline -----

/// One `(prime, α)` run: `c` and `h` at degree `n`, as dense residue slabs.
fn run_point(n: usize, deg: &Degrees, ex: &Exact, alpha: u64, p: mp::Md) -> (Vec<u64>, Vec<u64>) {
    // Φ_k, dense over (λ, μ, ν).
    let mut phi: Vec<Vec<u64>> = Vec::with_capacity(n + 1);
    phi.push(vec![1u64]); // Φ_0 = 1
    for k in 1..=n {
        let pk = deg.len(k);
        let mut slab = vec![0u64; deg.cube(k)];
        for (t, row) in ex.jp[k].iter().enumerate() {
            // a_θρ at this α, and the shared 1/⟨J_θ,J_θ⟩ folded into one factor.
            let a: Vec<u64> = row
                .iter()
                .map(|(num, scale)| {
                    if num.is_empty() {
                        return 0;
                    }
                    let mut acc = 0u64;
                    for c in num.iter().rev() {
                        acc = mp::add(mp::mul(acc, alpha, p), mp::from_i128(*c, p), p);
                    }
                    mp::mul(acc, mp::inv(p.from_u128(*scale), p), p)
                })
                .collect();
            let mut den = 1u64;
            for &((u, v), m) in &ex.norm[k][t] {
                // α > 0 and u, v ≥ 0 not both zero, so this is never zero.
                let f = mp::add(mp::mul(u as u64 % p.p, alpha, p), v as u64 % p.p, p);
                for _ in 0..m {
                    den = mp::mul(den, f, p);
                }
            }
            // One inversion, not one per ρ: `inv` is a 61-squaring ladder.
            let inv_den = mp::inv(den, p);
            let scaled: Vec<u64> = a.iter().map(|&x| mp::mul(x, inv_den, p)).collect();
            // The rank-1 tensor a' ⊗ a ⊗ a, accumulated.
            for i in 0..pk {
                if scaled[i] == 0 {
                    continue;
                }
                for j in 0..pk {
                    if a[j] == 0 {
                        continue;
                    }
                    let ij = mp::mul(scaled[i], a[j], p);
                    let base = (i * pk + j) * pk;
                    for l in 0..pk {
                        if a[l] == 0 {
                            continue;
                        }
                        slab[base + l] = mp::add(slab[base + l], mp::mul(ij, a[l], p), p);
                    }
                }
            }
        }
        phi.push(slab);
    }

    // G_k = k·Φ_k − Σ_{j<k} G_j·Φ_{k−j}, so that Ψ_k = α·G_k.
    let mut g: Vec<Vec<u64>> = Vec::with_capacity(n + 1);
    g.push(vec![0u64]);
    for k in 1..=n {
        let pk = deg.len(k);
        let mut acc: Vec<u64> = phi[k]
            .iter()
            .map(|&x| mp::mul(x, k as u64 % p.p, p))
            .collect();
        for j in 1..k {
            let m = k - j;
            let (pj, pm) = (deg.len(j), deg.len(m));
            let (uj, gj, pm_slab) = (&deg.union[j][m], &g[j], &phi[m]);
            for a1 in 0..pj {
                for b1 in 0..pj {
                    for c1 in 0..pj {
                        let lhs = gj[(a1 * pj + b1) * pj + c1];
                        if lhs == 0 {
                            continue;
                        }
                        for a2 in 0..pm {
                            let ta = uj[a1 * pm + a2] as usize;
                            for b2 in 0..pm {
                                let tb = uj[b1 * pm + b2] as usize;
                                let base = (ta * pk + tb) * pk;
                                for c2 in 0..pm {
                                    let rhs = pm_slab[(a2 * pm + b2) * pm + c2];
                                    if rhs == 0 {
                                        continue;
                                    }
                                    let tc = uj[c1 * pm + c2] as usize;
                                    let v = mp::mul(lhs, rhs, p);
                                    acc[base + tc] = mp::sub(acc[base + tc], v, p);
                                }
                            }
                        }
                    }
                }
            }
        }
        g.push(acc);
    }

    // c = z_λ·α^{ℓ(λ)}·Φ_n and h = α·G_n.
    let pn = deg.len(n);
    let mut c = vec![0u64; deg.cube(n)];
    for i in 0..pn {
        let lam = &deg.parts[n][i];
        let w = mp::mul(p.from_u128(lam.z()), mp::pow(alpha, lam.len() as u64, p), p);
        for j in 0..pn {
            for l in 0..pn {
                let idx = (i * pn + j) * pn + l;
                c[idx] = mp::mul(phi[n][idx], w, p);
            }
        }
    }
    let h: Vec<u64> = g[n].iter().map(|&x| mp::mul(x, alpha, p)).collect();
    (c, h)
}

/// Everything about one prime that the reconstruction pass consumes: for each
/// key, the interpolated coefficients of the answer in `b`, mod that prime.
type Rows = (Vec<Vec<u64>>, Vec<Vec<u64>>);

/// Run the whole pipeline at `points` values of α mod `p`, then interpolate.
fn rows_for_prime(n: usize, deg: &Degrees, ex: &Exact, points: usize, p: mp::Md) -> Rows {
    let xs: Vec<u64> = (1..=points as u64).collect();
    let mut cs: Vec<Vec<u64>> = Vec::with_capacity(points);
    let mut hs: Vec<Vec<u64>> = Vec::with_capacity(points);
    for &x in &xs {
        let (c, h) = run_point(n, deg, ex, x, p);
        cs.push(c);
        hs.push(h);
    }
    // Transpose: per key, the values across evaluation points.
    let size = deg.cube(n);
    let tmat = mp::lagrange_matrix(&xs, 1, p);
    let mut cpoly = Vec::with_capacity(size);
    let mut hpoly = Vec::with_capacity(size);
    let mut cv = vec![0u64; points];
    let mut hv = vec![0u64; points];
    for idx in 0..size {
        for k in 0..points {
            cv[k] = cs[k][idx];
            hv[k] = hs[k][idx];
        }
        cpoly.push(mp::apply_matrix(&tmat, &cv, p));
        hpoly.push(mp::apply_matrix(&tmat, &hv, p));
    }
    (cpoly, hpoly)
}

/// Assemble the answer from `rows`, holding the last prime back as a check.
///
/// Returns the tables and the number of keys that failed **for range reasons**
/// — the lift found no small rational, or the held-back prime disagreed. Those
/// are the failures more primes can fix, and they are counted separately from a
/// nonvanishing spare coefficient, which no amount of range would change.
fn assemble(n: u32, deg: &Degrees, rows: &[Rows]) -> (GjTables, usize) {
    let mut t = GjTables {
        n,
        c: BTreeMap::new(),
        h: BTreeMap::new(),
        not_polynomial: Vec::new(),
        c_not_integral: Vec::new(),
        negative: Vec::new(),
        peak_bits: 0,
        primes_used: rows.len() as u32,
    };
    let pn = deg.len(n as usize);
    let size = deg.cube(n as usize);
    let lift_count = rows.len() - 1;
    let lift_primes: Vec<mp::Md> = (0..lift_count).map(mp::nth_prime).collect();
    let check = mp::nth_prime(lift_count);
    let modulus: u128 = lift_primes.iter().map(|q| q.p as u128).product();
    let bound = mp::bound_for(modulus);
    let mut short = 0usize;
    for idx in 0..size {
        let (i, j, l) = (idx / (pn * pn), (idx / pn) % pn, idx % pn);
        let key: Key = (
            deg.parts[n as usize][i].clone(),
            deg.parts[n as usize][j].clone(),
            deg.parts[n as usize][l].clone(),
        );
        for (tag, which) in [("c", 0usize), ("h", 1usize)] {
            let all: Vec<&Vec<u64>> = rows
                .iter()
                .map(|r| if which == 0 { &r.0[idx] } else { &r.1[idx] })
                .collect();
            // The spare top coefficients must vanish at every prime. That
            // requirement IS the statement that the answer is a polynomial of
            // the expected degree, so it is verified rather than assumed. A
            // failure here is not a range problem — more primes would report it
            // just as loudly — so it does not count toward `short`.
            if all
                .iter()
                .any(|poly| poly.iter().skip(n as usize + 1).any(|&v| v != 0))
            {
                t.not_polynomial.push((tag, key.clone()));
                continue;
            }
            let Some(poly) = lift(&all[..lift_count], &lift_primes, modulus, bound) else {
                t.not_polynomial.push((tag, key.clone()));
                short += 1;
                continue;
            };
            // Independent verification: reduce the reconstructed rational mod
            // the held-back prime and require it to be what that prime
            // computed. Reconstruction returns a spurious small rational rather
            // than failing when the true value is too big, so a large bound is
            // not by itself evidence that it was large enough.
            if !agrees_mod(&poly, all[lift_count], check) {
                t.not_polynomial.push((tag, key.clone()));
                short += 1;
                continue;
            }
            if poly.num.is_empty() {
                continue;
            }
            for &v in &poly.num {
                t.peak_bits = t.peak_bits.max(128 - v.unsigned_abs().leading_zeros());
            }
            if tag == "c" && !poly.is_integral() {
                t.c_not_integral.push(key.clone());
            }
            if !poly.is_positive() {
                t.negative.push((tag, key.clone(), poly.clone()));
            }
            if which == 0 {
                t.c.insert(key.clone(), poly);
            } else {
                t.h.insert(key.clone(), poly);
            }
        }
    }
    (t, short)
}

/// Both \[GJ\] tables at degree `n`, by modular evaluation and interpolation.
///
/// At `n = 0` both tables come back empty and no prime is run.
///
/// Runs the whole pipeline at `points` values of α over several primes, CRTs
/// all but one to reconstruct the rationals, and holds the last back as an
/// independent check — which is what replaces the exact engine's free "the
/// denominator collapsed" law. `points` is chosen with spare capacity and the
/// surplus top coefficients are **required to vanish**; that requirement *is*
/// the statement that the answer is a polynomial of the expected degree.
///
/// ## The range grows on demand
///
/// Three lifting primes give a reconstruction bound of `2^46`, and the widest
/// coefficient measured runs about `3n` bits — 35 at n = 12. So the bound is
/// not a constant that can be checked once and forgotten: somewhere past n = 15
/// the coefficients pass it, and a fixed prime count would report that as
/// *non-polynomiality*, blaming \[DF\] for this crate's arithmetic.
///
/// Instead the failures that more primes could fix are counted, and while there
/// are any the engine adds a prime and reassembles. Each prime costs one full
/// evaluation pass, so this is paid only at the degrees that need it — and once
/// the prime budget is exhausted the remaining failures are reported as
/// findings, which by then they have earned.
///
/// # Panics
///
/// Panics if a reconstructed denominator leaves `i128` — either the common
/// denominator scaled onto one coefficient, or the gcd cancelled back out of
/// it. Panics if a `J → p` coefficient carries a linear atom, which would make
/// its evaluation at numeric α wrong.
pub fn gj_connection_tables_modular(n: u32) -> GjTables {
    /// Three to lift plus one to check: the smallest set that both reaches past
    /// a single prime's useless `2^15` bound and verifies the result.
    const MIN_PRIMES: usize = 4;
    /// `2^124`, hence a `2^61` bound — past the point where the `i128`
    /// coefficients themselves would overflow, so nothing is gained by more.
    const MAX_PRIMES: usize = 5;

    if n == 0 {
        return GjTables {
            n,
            c: BTreeMap::new(),
            h: BTreeMap::new(),
            not_polynomial: Vec::new(),
            c_not_integral: Vec::new(),
            negative: Vec::new(),
            peak_bits: 0,
            primes_used: 0,
        };
    }
    let deg = Degrees::new(n);
    let ex = Exact::new(n, &deg);

    // `c` and `h` are polynomials in b of degree ≤ n−1 at every degree measured,
    // so n+4 points leaves spare — enough that "the top ones vanish" is a real
    // test and not a tautology.
    let points = n as usize + 4;

    let mut rows: Vec<Rows> = Vec::new();
    let mut want = MIN_PRIMES;
    loop {
        while rows.len() < want {
            let p = mp::nth_prime(rows.len());
            rows.push(rows_for_prime(n as usize, &deg, &ex, points, p));
        }
        let (t, short) = assemble(n, &deg, &rows);
        if short == 0 || want >= MAX_PRIMES {
            return t;
        }
        want += 1;
    }
}

/// Does `poly`, reduced mod `p`, equal what that prime independently computed?
fn agrees_mod(poly: &BPoly, want: &[u64], p: mp::Md) -> bool {
    let inv_den = mp::inv(p.from_u128(poly.den), p);
    (0..want.len()).all(|k| {
        let got = poly
            .num
            .get(k)
            .map_or(0, |&v| mp::mul(mp::from_i128(v, p), inv_den, p));
        got == want[k]
    })
}

/// Reconstruct a dense residue polynomial as a [`BPoly`] over one denominator,
/// from its residues at two primes.
fn lift(rows: &[&Vec<u64>], primes: &[mp::Md], modulus: u128, bound: u128) -> Option<BPoly> {
    let width = rows[0].len();
    let mut nums: Vec<i128> = Vec::with_capacity(width);
    let mut dens: Vec<u128> = Vec::with_capacity(width);
    let mut residues = vec![0u64; rows.len()];
    for k in 0..width {
        for (r, row) in residues.iter_mut().zip(rows.iter()) {
            *r = row[k];
        }
        let (x, m) = mp::crt(&residues, primes);
        debug_assert_eq!(m, modulus);
        let (num, den) = mp::reconstruct(x, modulus, bound)?;
        nums.push(num);
        dens.push(den);
    }
    let mut common: u128 = 1;
    for &d in &dens {
        common = common / mp::gcd128(common, d) * d;
    }
    let num: Vec<i128> = nums
        .iter()
        .zip(dens.iter())
        // `common` is an lcm of denominators and can in principle outgrow the
        // signed half of the width; it is a *value* here, not an index, so it
        // checks rather than proving (R5).
        .map(|(&n, &d)| {
            let scale = i128::try_from(common / d)
                .unwrap_or_else(|_| panic!("the scale {common}/{d} does not fit i128"));
            n * scale
        })
        .collect();
    let mut out = BPoly { num, den: common };
    while out.num.last() == Some(&0) {
        out.num.pop();
    }
    let mut g = out.den;
    for &v in &out.num {
        g = mp::gcd128(g, v.unsigned_abs());
    }
    if g > 1 {
        out.den /= g;
        // `g` divides `out.den` and every `|num|`, so it is bounded by the
        // magnitudes it divides — but those are `u128` on the denominator side,
        // so the narrowing is checked rather than assumed.
        let g = i128::try_from(g).unwrap_or_else(|_| panic!("the gcd {g} does not fit i128"));
        for v in &mut out.num {
            *v /= g;
        }
    }
    Some(out)
}

/// Assert the exact and modular engines agree, for tests and examples.
///
/// Returns the number of entries compared, `c` and `h` together.
///
/// They share `Partition`, `AFrac` (only to *build* `J → p`) and nothing else:
/// one works in ℚ(α) with factored linear denominators throughout, the other
/// never forms a rational function at all. Agreement is evidence.
///
/// # Errors
///
/// Returns `Err` if the modular engine reports a broken law, or if the two
/// tables differ in entry count. Returns `Err` if a key is missing from the
/// modular table or disagrees with the exact one.
pub fn engines_agree(n: u32) -> Result<usize, String> {
    let exact = crate::gj::gj_connection_tables(n);
    let modular = gj_connection_tables_modular(n);
    if !modular.laws_hold() {
        return Err(format!(
            "modular engine reported a broken law at n = {n}: {:?} / {:?}",
            modular.not_polynomial, modular.c_not_integral
        ));
    }
    let mut checked = 0;
    for (name, a, b) in [("c", &exact.c, &modular.c), ("h", &exact.h, &modular.h)] {
        if a.len() != b.len() {
            return Err(format!(
                "{name} at n = {n}: exact has {} entries, modular {}",
                a.len(),
                b.len()
            ));
        }
        for (key, want) in a {
            match b.get(key) {
                Some(got) if got == want => checked += 1,
                Some(got) => {
                    return Err(format!(
                        "{name} at {key:?}: exact {want:?}, modular {got:?}"
                    ))
                }
                None => return Err(format!("{name} at {key:?}: missing from the modular table")),
            }
        }
    }
    Ok(checked)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **The two engines agree**, entry for entry, on both tables.
    ///
    /// One works in ℚ(α) with factored linear denominators throughout; the
    /// other never forms a rational function. They share no arithmetic, so this
    /// is evidence rather than tautology — V3 in `docs/policies/validation.md`.
    #[test]
    fn the_modular_engine_agrees_with_the_exact_one() {
        // The [GJ] tables are empty at n = 0 by construction, so a sweep from 0
        // asserts nothing; that edge is pinned by the_tables_are_empty_at_zero.
        for n in 1..=6u32 {
            let checked = engines_agree(n).unwrap_or_else(|e| panic!("{e}"));
            assert!(checked > 0, "n = {n} compared nothing");
        }
    }

    /// The degree bound is checked, not assumed: the spare interpolation points
    /// must come back zero.
    #[test]
    fn the_spare_coefficients_vanish() {
        // The [GJ] tables are empty at n = 0 by construction, so a sweep from 0
        // asserts nothing; that edge is pinned by the_tables_are_empty_at_zero.
        for n in 1..=5u32 {
            let t = gj_connection_tables_modular(n);
            assert!(
                t.not_polynomial.is_empty(),
                "n = {n}: an interpolated answer exceeded its degree bound or \
                 failed reconstruction: {:?}",
                t.not_polynomial
            );
            for poly in t.c.values().chain(t.h.values()) {
                assert!(
                    poly.num.len() <= n as usize,
                    "n = {n}: b-degree {} exceeds the measured n−1 bound",
                    poly.num.len() - 1
                );
            }
        }
    }

    /// `b = 0` is the class algebra of `S_n` — the law that has to catch an
    /// error once the "denominators collapsed" check is gone.
    #[test]
    fn the_modular_b_zero_slice_is_the_class_algebra() {
        // The [GJ] tables are empty at n = 0 by construction, so a sweep from 0
        // asserts nothing; that edge is pinned by the_tables_are_empty_at_zero.
        for n in 1..=6u32 {
            let t = gj_connection_tables_modular(n);
            for la in crate::partitions_of(n) {
                for mu in crate::partitions_of(n) {
                    for nu in crate::partitions_of(n) {
                        let key = (la.clone(), mu.clone(), nu.clone());
                        let got = t.c.get(&key).map_or((0, 1), BPoly::at_zero);
                        let want = crate::class_algebra_coefficient(&la, &mu, &nu);
                        assert_eq!(got, (want, 1), "c^{la}_{{{mu},{nu}}}(0)");
                    }
                }
            }
        }
    }

    /// The range really does run out, and adding a prime really does fix it.
    ///
    /// Both halves matter. The engine grows its prime count on a *count of
    /// failures*, so if that count could never be nonzero the whole mechanism
    /// would be untested scaffolding — and if a failure were not curable by one
    /// more prime, growing would be the wrong response to it.
    ///
    /// Degree 8 is chosen because its coefficients reach 18 bits: comfortably
    /// inside the `2^46` that three lifting primes give, and comfortably
    /// outside the `2^15` that one gives.
    #[test]
    fn the_reconstruction_bound_is_reached_and_then_escaped() {
        let n = 8u32;
        let deg = Degrees::new(n);
        let ex = Exact::new(n, &deg);
        let points = n as usize + 4;
        let rows: Vec<Rows> = (0..4)
            .map(|k| rows_for_prime(n as usize, &deg, &ex, points, mp::nth_prime(k)))
            .collect();

        // One lifting prime plus the check: a bound of 2^15 against 18-bit
        // answers, so the lift must come up short somewhere.
        let (thin, short) = assemble(n, &deg, &rows[..2]);
        assert!(
            short > 0,
            "a 2^15 bound was somehow enough for 18-bit values"
        );
        assert!(!thin.laws_hold(), "and it must be reported, not swallowed");

        // Three plus the check, which is what the driver settles on here.
        let (full, short) = assemble(n, &deg, &rows[..4]);
        assert_eq!(short, 0, "three primes must be enough at n = 8");
        assert!(full.laws_hold());
        assert_eq!(full.primes_used, 4);

        // And the driver's own answer is that one, having grown to it or not.
        let driver = gj_connection_tables_modular(n);
        assert_eq!(driver.c, full.c);
        assert_eq!(driver.h, full.h);
        assert_eq!(driver.primes_used, 4);
    }
}
