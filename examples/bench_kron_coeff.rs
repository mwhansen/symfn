//! Single Kronecker coefficient: the character sum against the product route.
//!
//! `ops::kronecker` reaches one coefficient by building the whole internal
//! product — `s → p`, a diagonal multiply, `p → s` back — and reading a term
//! off. `ops::kronecker_coeff` never forms a symmetric function: it sums
//! `χ^λ(ρ)χ^μ(ρ)χ^ν(ρ)/z_ρ` over ρ ⊢ n.
//!
//! Two things this is meant to establish, in order:
//!
//! 1. **Where the crossover is.** The product route amortises across ν, so it
//!    should win at small n where `p → s` is cheap. Guessing which side wins
//!    where is exactly the mistake `docs/record/README.md` records three times.
//! 2. **How the two scale once both are exact.** Both routes are run over
//!    `BigRational` here, so the comparison is like-for-like at every degree —
//!    which it was not before `Partition::z_in` and `Partition::div_by_z`
//!    removed the `u128` cap at degree 34.
//!
//! Run it:
//!
//! ```text
//!   cargo run --release --features bignum --example bench_kron_coeff
//! ```

use num_rational::BigRational;
use num_traits::ToPrimitive;
use std::time::Instant;
use symfn::ops::{kronecker, kronecker_coeff};
use symfn::partition::Partition;

/// Above this the product route is left out: it is `p(n)²` character lookups,
/// and past the low 30s that is minutes per point. Not a wall — a budget, and
/// stated so the empty cells do not read as "cannot".
const PRODUCT_BUDGET_TO: u32 = 32;

fn main() {
    println!("   n  coefficient   product route   character sum   ratio");
    println!("  --  -----------  --------------  --------------  ------");

    for n in [8u32, 12, 16, 20, 24, 28, 32, 36, 40, 44] {
        // A shape family with enough rows to be non-trivial, scaled with n.
        let lambda = Partition::new(vec![n - 5, 3, 2]);
        let mu = Partition::new(vec![n - 6, 4, 2]);
        let nu = Partition::new(vec![n - 4, 3, 1]);

        let t = Instant::now();
        let got = kronecker_coeff(&lambda, &mu, &nu);
        let char_time = t.elapsed();

        if n <= PRODUCT_BUDGET_TO {
            let t = Instant::now();
            let want: BigRational = kronecker(&lambda, &mu, &nu);
            let prod_time = t.elapsed();
            assert!(want.is_integer(), "g is not an integer at n = {n}");
            assert_eq!(
                want.to_integer().to_i128(),
                Some(got),
                "routes disagree at n = {n}: g^{nu}_{{{lambda},{mu}}}"
            );
            let ratio = prod_time.as_secs_f64() / char_time.as_secs_f64();
            println!("  {n:2}  {got:11}  {prod_time:>14.2?}  {char_time:>14.2?}  {ratio:>5.2}x");
        } else {
            println!(
                "  {n:2}  {got:11}  {:>14}  {char_time:>14.2?}  {:>6}",
                "-", "-"
            );
        }
    }

    println!();
    println!("  Both routes run over BigRational, so every row is exact and the");
    println!("  comparison is like-for-like. '-' is the {PRODUCT_BUDGET_TO}-degree budget above,");
    println!("  not a limit of the product route.");
}
