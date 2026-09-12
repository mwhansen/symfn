//! Exhibit a *nonzero* structure constant on the pair whose full product
//! cannot be materialized (`S_13 ℓ=25,36`, monomial mass 4.3×10¹⁶).
//!
//! Random `w` of the right length are essentially never in the support, so the
//! earlier sampling found only zeros. This walks the support instead.
//!
//! **Propose:** start from `S_u` and apply the Monk chain of `x^{code(v)}`,
//! keeping only the `K` largest-coefficient terms after each pass. Every
//! surviving term is automatically a candidate of exactly the right degree —
//! each Monk pass raises `ℓ` by one and `Σ code(v) = ℓ(v)` — and automatically
//! `≥ u` in Bruhat order, since Monk only moves up.
//!
//! **Dispose:** the beam is *not* the answer. Truncation discards terms that
//! would have canceled, so its coefficients are wrong. Every candidate is
//! re-checked with the exact pruned query `schubert_coeff`, which is the thing
//! actually verified against the full product. The beam only decides *where to
//! look*.

use std::time::Instant;
use symfn::permutation::Perm;
use symfn::schubert::{schubert_coeff, Schubert};

fn main() {
    let beam: usize = std::env::var("BEAM")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(4000);
    // CASE=132 runs the *computable* neighbor instead, where the full
    // product exists and can referee both the beam and the query. Validating
    // the pipeline there is the only way to trust its answers on S_13.0,
    // which by construction has no referee.
    let computable = std::env::var("CASE").as_deref() == Ok("132");
    let (u, v) = if computable {
        (
            Perm::new([6, 7, 10, 1, 12, 11, 8, 3, 9, 2, 5, 13, 4]).unwrap(),
            Perm::new([2, 12, 9, 11, 4, 10, 6, 3, 1, 8, 7, 5, 13]).unwrap(),
        )
    } else {
        (
            Perm::new([1, 8, 2, 12, 9, 3, 4, 6, 11, 7, 5, 13, 10]).unwrap(),
            Perm::new([1, 12, 7, 3, 2, 11, 8, 6, 13, 10, 5, 9, 4]).unwrap(),
        )
    };
    let deg = u.length() + v.length();
    println!(
        "u: ℓ={}  v: ℓ={}  target degree {deg}  beam {beam}",
        u.length(),
        v.length()
    );

    let t = Instant::now();
    let mut f: Schubert<i128> = Schubert::monomial(u, 1);
    for (i, &e) in v.code().iter().enumerate() {
        for _ in 0..e {
            f = f.mul_variable(i as u32 + 1);
            if f.terms().len() > beam {
                let mut v: Vec<_> = f.terms().iter().map(|(w, c)| (*w, *c)).collect();
                v.sort_by_key(|(_, c)| -c.abs());
                v.truncate(beam);
                let mut g = Schubert::zero();
                for (w, c) in v {
                    g.add_term(w, &c);
                }
                f = g;
            }
        }
    }
    println!(
        "beam produced {} candidates in {:.2}s",
        f.terms().len(),
        t.elapsed().as_secs_f64()
    );

    // Dispose: exact query on the candidates, best-looking first.
    let mut cands: Vec<Perm> = f
        .terms()
        .iter()
        .filter(|(w, _)| w.length() == deg && v.bruhat_le(w))
        .map(|(w, _)| *w)
        .collect();
    cands.sort_by_key(|w| std::cmp::Reverse(f.coeff(w).abs()));
    println!(
        "{} of them satisfy v ≤ w (u ≤ w holds by construction)",
        cands.len()
    );

    let mut found = 0;
    let t = Instant::now();
    for (k, w) in cands.iter().enumerate() {
        let c: i128 = schubert_coeff(&u, &v, w);
        if c != 0 {
            found += 1;
            println!("  NONZERO  c^{w}_uv = {c}   (candidate #{k})");
            if found == 5 {
                break;
            }
        }
        if k >= 400 {
            break;
        }
    }
    println!(
        "{found} nonzero found, {:.2}s of exact queries",
        t.elapsed().as_secs_f64()
    );

    if computable {
        // Referee: build the whole product and check the query against it on
        // every candidate, zeros included -- a pruning bug shows up as a
        // spurious zero, so the zeros are the part worth checking.
        let a: Schubert<i128> = Schubert::monomial(u, 1);
        let b: Schubert<i128> = Schubert::monomial(v, 1);
        let t = Instant::now();
        let full = a.mul_e2(&b);
        println!(
            "referee: full product {} terms in {:.2}s",
            full.terms().len(),
            t.elapsed().as_secs_f64()
        );
        let (mut ok, mut bad, mut nz) = (0, 0, 0);
        for w in cands.iter().take(400) {
            let got: i128 = schubert_coeff(&u, &v, w);
            let want = full.coeff(w);
            if got == want {
                ok += 1;
                if want != 0 {
                    nz += 1;
                }
            } else {
                bad += 1;
                if bad <= 3 {
                    println!("  MISMATCH c^{w} query={got} full={want}");
                }
            }
        }
        println!("query vs full product on {ok} candidates: {bad} mismatches ({nz} nonzero)");
    }
}
