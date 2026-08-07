//! Dump the Macdonald operators in the Schur basis for
//! `scripts/check_deltaop.py`.
//!
//! ```text
//!   cargo run --release --example deltaop_dump -- 8 > /tmp/dop.txt
//! ```
//!
//! One line per (operator, input, λ) with a nonzero coefficient:
//!
//! ```text
//!   op | arg | lambda | a,b:c;a,b:c;…
//! ```
//!
//! where `a,b:c` is `c·qᵃtᵇ`. Every coefficient is a **polynomial** — that is
//! the point of the operators, and a surviving denominator panics upstream
//! rather than being printed.
//!
//! Only `nabla` has an external oracle (Sage has no Δ, Δ' or Θ), so the other
//! rows exist to be checked against the *identities* in the Python script.

use symfn::coeff::Rational;
use symfn::sym::SymFn;
use symfn::{QtPoly, Ring, Schur};

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);

    for n in 1..=top {
        // nabla on e_n, h_n, and every Schur function of the degree.
        emit("nabla", &format!("e{n}"), &symfn::nabla_e::<Rational>(n));
        for lambda in symfn::partitions_of(n) {
            let s: Schur<QtPoly<Rational>> =
                Schur::monomial(lambda.clone(), <QtPoly<Rational> as Ring>::one());
            emit(
                "nabla",
                &format!("s{}", join(lambda.parts())),
                &symfn::nabla(&s),
            );
        }
        // nabla^2, the object [QZ] 2026 is about.
        for lambda in symfn::partitions_of(n) {
            let s: Schur<QtPoly<Rational>> =
                Schur::monomial(lambda.clone(), <QtPoly<Rational> as Ring>::one());
            emit(
                "nabla2",
                &format!("s{}", join(lambda.parts())),
                &symfn::nabla_power(&s, 2),
            );
        }
        // the Delta conjecture family
        for k in 0..n {
            emit(
                "deltaprime",
                &format!("e{k},e{n}"),
                &symfn::delta_prime_e::<Rational>(k, n),
            );
        }
        // The labeled-Dyck-path generating function, for the one slice Sage
        // can independently confirm: distinct labels and k = n-1, where the
        // z-selection is empty and this is sum q^dinv t^area over parking
        // functions.
        {
            let ones = symfn::Partition::new(std::iter::repeat_n(1, n as usize));
            let g = symfn::side_at_content::<Rational>(&ones, n - 1, symfn::Side::Rise);
            println!("dyck | pf{n} |  | {}", poly(&g));
        }
        // Theta, on the composite the theorem constrains
        for k in 1..n {
            let inner = symfn::nabla_e::<Rational>(n - k);
            emit(
                "theta",
                &format!("e{k},nabla_e{}", n - k),
                &symfn::theta(&symfn::deltaop::elementary(k), &inner),
            );
        }
    }
}

fn emit(op: &str, arg: &str, f: &Schur<QtPoly<Rational>>) {
    for (lambda, c) in f.terms() {
        println!("{op} | {arg} | {} | {}", join(lambda.parts()), poly(c));
    }
}

fn join(v: &[u32]) -> String {
    v.iter().map(u32::to_string).collect::<Vec<_>>().join(",")
}

/// Coefficients must be **integers**: every operator here maps ℤ[q,t]-Schur
/// combinations to ℤ[q,t] ones, and a surviving `1/2` would be a bug the
/// comparison could not see.
fn poly(p: &QtPoly<Rational>) -> String {
    p.terms()
        .map(|(&(a, b), c)| {
            assert_eq!(
                c.denom(),
                1,
                "coefficient {c:?} of q^{a}t^{b} is not an integer"
            );
            format!("{a},{b}:{}", c.numer())
        })
        .collect::<Vec<_>>()
        .join(";")
}
