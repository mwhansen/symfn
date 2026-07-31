//! A CLI mirroring `lrcalc`'s interface, for head-to-head comparison.
//!
//! Deliberately argument- and output-compatible with lrcalc so the two can be
//! run on identical inputs and diffed:
//!
//! ```text
//!   lr_cli mult 3 2 1 - 2 1
//!   lr_cli skew 3 2 1 / 2 1
//!   lr_cli coef 3 2 1 - 2 1 - 2 1
//! ```
//!
//! Output is `COEFF  (p1,p2,...)` per line, one per nonzero term, exactly as
//! lrcalc prints it (lrcalc emits terms in hash order; `scripts/compare_lrcalc.py`
//! sorts both sides before comparing). `coef` prints a single integer.
//!
//! Running as a separate process per case is the point: it makes the comparison
//! fair. Both sides pay process startup and both format their full output, and
//! neither gets to reuse a warm memo cache across cases.

use symfn::{expand_skew, AutoLr, LrBackend, Partition};

fn parse(args: &[String]) -> Partition {
    Partition::new(args.iter().map(|a| a.parse::<u32>().expect("integer part")))
}

/// Split on a literal separator token, e.g. `-` or `/`.
fn split_on(args: &[String], sep: &str) -> Vec<Vec<String>> {
    let mut out = vec![Vec::new()];
    for a in args {
        if a == sep {
            out.push(Vec::new());
        } else {
            out.last_mut().unwrap().push(a.clone());
        }
    }
    out
}

fn print_terms(terms: &[(Partition, u128)]) {
    use std::fmt::Write as _;
    // One buffer, no per-term allocation: at tens of thousands of terms the
    // old format!-per-term path cost ~360ns a line, which is the same order
    // as the whole DP on mid-sized products.
    let mut s = String::with_capacity(terms.len() * 24);
    for (lambda, c) in terms {
        write!(s, "{c}  (").expect("write to String is infallible");
        let mut sep = "";
        for p in lambda.parts() {
            write!(s, "{sep}{p}").expect("write to String is infallible");
            sep = ",";
        }
        s.push_str(")\n");
    }
    print!("{s}");
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() {
        eprintln!("usage: lr_cli mult A - B | skew A / B | coef C - A - B");
        std::process::exit(2);
    }
    let rest = &args[1..];
    match args[0].as_str() {
        "mult" => {
            let parts = split_on(rest, "-");
            assert_eq!(parts.len(), 2, "mult takes: A - B");
            print_terms(&AutoLr.schur_product(&parse(&parts[0]), &parse(&parts[1])));
        }
        "skew" => {
            let parts = split_on(rest, "/");
            assert_eq!(parts.len(), 2, "skew takes: A / B");
            print_terms(&expand_skew(&parse(&parts[0]), &parse(&parts[1])));
        }
        "coef" => {
            let parts = split_on(rest, "-");
            assert_eq!(parts.len(), 3, "coef takes: C - A - B");
            println!(
                "{}",
                AutoLr.lr_coeff(&parse(&parts[0]), &parse(&parts[1]), &parse(&parts[2]))
            );
        }
        other => {
            eprintln!("unknown operation {other:?}");
            std::process::exit(2);
        }
    }
}
