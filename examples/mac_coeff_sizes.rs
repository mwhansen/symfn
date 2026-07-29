//! How wide do Macdonald numerator coefficients get, and where does `i128` stop
//! being safe?
//!
//! ```text
//!   cargo run --release --features bignum --example mac_coeff_sizes -- 10
//! ```
//!
//! Runs each degree over `Frac<i128>` and over `Frac<BigInt>` and compares. A
//! wrapped `i128` is silent — `impl_ring_for_int` uses plain `*` and `+` — so
//! the only honest way to find the ceiling is to compute the same thing twice
//! in two widths and look for the first disagreement.

#[cfg(not(feature = "bignum"))]
fn main() {
    eprintln!("needs --features bignum");
}

#[cfg(feature = "bignum")]
fn main() {
    use num_bigint::BigInt;
    use symfn::sym::SymFn;
    use symfn::{macdonald_p, Frac, Monomial};

    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    println!(
        "{:>3}  {:>26}  {:>8}  {}",
        "n", "max |numerator coeff|", "bits", "i128 exact?"
    );
    for n in 1..=top {
        let mut widest = BigInt::from(0);
        let mut agree = true;
        for lambda in symfn::partitions_of(n) {
            let small: Monomial<Frac<i128>> = macdonald_p(&lambda);
            let big: Monomial<Frac<BigInt>> = macdonald_p(&lambda);
            for (_mu, c) in big.terms() {
                let (num, _) = c.parts();
                for (_, v) in num.terms() {
                    let a = if *v < BigInt::from(0) {
                        -v.clone()
                    } else {
                        v.clone()
                    };
                    if a > widest {
                        widest = a;
                    }
                }
            }
            // Compare the two runs term by term, via the expanded denominator so
            // the factored forms need not match structurally.
            for (mu, cb) in big.terms() {
                let cs = small.coeff(mu);
                let (sn, _) = cs.parts();
                let (bn, _) = cb.parts();
                let sn_big: Vec<_> = sn.terms().map(|(k, v)| (*k, BigInt::from(*v))).collect();
                let bn_v: Vec<_> = bn.terms().map(|(k, v)| (*k, v.clone())).collect();
                if sn_big != bn_v || cs.denominator().len() != cb.denominator().len() {
                    agree = false;
                }
            }
        }
        let bits = widest.bits();
        println!(
            "{n:>3}  {:>26}  {bits:>8}  {}",
            widest.to_string(),
            if agree {
                "yes"
            } else {
                "NO — i128 has wrapped"
            }
        );
    }
}
