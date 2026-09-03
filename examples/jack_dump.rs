//! Dump Jack expansions in a form `scripts/check_jack.py` can hold to Sage.
//!
//! ```text
//!   cargo run --release --example jack_dump -- 7 > /tmp/jack.txt
//!   sage -python scripts/check_jack.py /tmp/jack.txt
//! ```
//!
//! One line per coefficient:
//!
//! ```text
//!   kind | lambda | mu | numerator | denominator | scale
//!   p|3,1|2,1,1|0,1,2|1:2:1;2:1:1|4
//! ```
//!
//! meaning the `mu` coefficient of `kind_lambda` is
//! `(0 + 1·α + 2·α²) / (4 · (α+2)^1 · (2α+1)^1)`. Kinds are `p`, `q`, `j` in the
//! monomial basis and `jp` for `J` in the **power-sum** basis, which is the
//! Jack-character table and the unit the [GJ] pipeline consumes. `mp`, `mq`
//! and `mj` are the inverse direction — `m_λ` written in each of the three
//! normalizations — where λ indexes the *monomial* and μ the Jack shape.
//!
//! Everything runs over `i128`, on purpose: the denominators are then a
//! property of the representation rather than something a ℚ could absorb, so a
//! mismatch in `scale` is a real mismatch.

use symfn::afrac::AFrac;
use symfn::sym::{Monomial, PowerSum, SymFn};
use symfn::{Partition, Ring};

fn cell(c: &AFrac<i128>) -> String {
    let (num, den, scale) = c.parts();
    let n: Vec<String> = num.iter().map(|x| x.to_string()).collect();
    let d: Vec<String> = den.map(|(&(u, v), &m)| format!("{u}:{v}:{m}")).collect();
    format!("{}|{}|{scale}", n.join(","), d.join(";"))
}

fn shape(p: &Partition) -> String {
    p.parts()
        .iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn dump_monomial(kind: &str, lambda: &Partition, f: &Monomial<AFrac<i128>>) {
    for (mu, c) in f.terms() {
        println!("{kind}|{}|{}|{}", shape(lambda), shape(mu), cell(c));
    }
}

fn dump_map(
    kind: &str,
    lambda: &Partition,
    f: &std::collections::BTreeMap<Partition, AFrac<i128>>,
) {
    for (mu, c) in f {
        println!("{kind}|{}|{}|{}", shape(lambda), shape(mu), cell(c));
    }
}

fn dump_powersum(kind: &str, lambda: &Partition, f: &PowerSum<AFrac<i128>>) {
    for (mu, c) in f.terms() {
        println!("{kind}|{}|{}|{}", shape(lambda), shape(mu), cell(c));
    }
}

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);
    // The tableau route is exponential; keep it to shapes it can actually run.
    let tab_top: u32 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(5);

    for n in 1..=top {
        for lambda in symfn::partitions_of(n) {
            let p: Monomial<AFrac<i128>> = symfn::jack_p(&lambda);
            dump_monomial("p", &lambda, &p);
            dump_monomial("q", &lambda, &symfn::jack_q(&lambda));
            dump_monomial("j", &lambda, &symfn::jack_j(&lambda));
            dump_powersum("jp", &lambda, &symfn::jack_j_powersum(&lambda));
            let m: Monomial<AFrac<i128>> =
                Monomial::monomial(lambda.clone(), <AFrac<i128> as Ring>::one());
            dump_map("mp", &lambda, &symfn::monomial_to_jack_p(&m));
            dump_map("mq", &lambda, &symfn::monomial_to_jack_q(&m));
            dump_map("mj", &lambda, &symfn::monomial_to_jack_j(&m));
            // The branching cross-check travels too, so the oracle sees both
            // engines rather than only whichever one `jack_p` dispatches to.
            let b: Monomial<AFrac<i128>> = symfn::jack::jack_p_branching(&lambda);
            assert_eq!(b, p, "E1 vs E2 disagree at {lambda}");
            if n <= tab_top {
                let t: Monomial<AFrac<i128>> = symfn::jack::jack_j_tableaux(&lambda);
                assert_eq!(t, symfn::jack_j(&lambda), "E3 disagrees at {lambda}");
            }
        }
    }

    // The closed-form norms — the table Sage prices like a full expansion
    // (>360 s at n = 12). Computed as a factor multiset and expanded only to
    // cross the wire; that it *is* factored is what `jack_norm_j`'s type says
    // and what `norms_and_the_dual_basis` checks.
    for n in 1..=top {
        for lambda in symfn::partitions_of(n) {
            let f = AFrac::<i128>::from_factors(&symfn::jack_norm_j(&lambda));
            println!("normj|{}||{}", shape(&lambda), cell(&f));
        }
    }

    // Stanley's structure constants, the open-conjecture object.
    for (na, nb) in [(2u32, 2u32), (2, 3), (3, 3)] {
        for la in symfn::partitions_of(na) {
            for mu in symfn::partitions_of(nb) {
                for nu in symfn::partitions_of(na + nb) {
                    let g: AFrac<i128> = symfn::jack_structure_constant(&la, &mu, &nu);
                    if g.is_zero() {
                        continue;
                    }
                    println!(
                        "g|{}+{}|{}|{}",
                        shape(&la),
                        shape(&mu),
                        shape(&nu),
                        cell(&g)
                    );
                }
            }
        }
    }
}
