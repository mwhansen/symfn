//! Time the four walks that compute one product `s_a · s_b` through
//! [`SkewLr`]'s layer, so the choice among them can be calibrated.
//!
//! A product is the skew expansion of a disconnected shape: one factor sits
//! above and to the right with its single canonical filling, the other is the
//! block the layer actually enumerates, with the first as its ballot offset.
//! Either factor can be the enumerated one, and the walk can run on the
//! diagram or on its conjugate — and conjugating the juxtaposed shape swaps
//! the roles, so the four walks enumerate `b`, `a`, `a'` and `b'`
//! respectively. All four produce the same coefficients; they differ in how
//! many row fillings the layer commits and how large it grows.
//!
//!   cargo run --release --example calibrate_orientation [reps] [a1,a2 b1,b2 ...]
//!
//! One line per walk: the pair, the walk, the enumerated block's rows, cells
//! and widest row, min time over `reps` in-process runs (order rotated so no
//! walk always inherits a warm allocator), productions, peak states, terms,
//! and — for the walk the library currently takes — a `<-` marker with its
//! ratio to the fastest. Every walk's expansion is checked equal to the
//! first's. Cases default to a spread of asymmetric pairs and squares.

use std::time::Instant;
use symfn::skew_lr::{take_peak_layer_states, take_productions};
use symfn::{clear_caches, expand_skew, LrBackend, Partition, SkewLr};

/// The disconnected shape whose Schur function is s_top · s_bottom.
fn juxtapose(top: &Partition, bottom: &Partition) -> (Partition, Partition) {
    let shift = bottom.part(0);
    let mut outer: Vec<u32> = top.parts().iter().map(|&p| p + shift).collect();
    outer.extend_from_slice(bottom.parts());
    let inner = vec![shift; if shift == 0 { 0 } else { top.len() }];
    (Partition::new(outer), Partition::new(inner))
}

fn parse(s: &str) -> Partition {
    Partition::new(s.split(',').map(|x| x.parse::<u32>().expect("shape part")))
}

struct Walk {
    /// Which factor is enumerated, as text.
    name: &'static str,
    top: Partition,
    bottom: Partition,
    conj: bool,
    /// The block the layer enumerates: `bottom` for a direct walk, `top'`
    /// for the conjugate one.
    block: Partition,
}

fn main() {
    let mut args = std::env::args().skip(1).peekable();
    let reps: usize = match args.peek() {
        Some(a) if a.parse::<usize>().is_ok() => args.next().unwrap().parse().unwrap(),
        _ => 3,
    };
    let explicit: Vec<Partition> = args.map(|a| parse(&a)).collect();
    let pairs: Vec<(Partition, Partition)> = if explicit.is_empty() {
        [
            // Asymmetric, four-row factors (the layer's own regime).
            ("16,13,10,7", "8,6,4,2"),
            ("20,16,12,8", "8,6,4,2"),
            ("18,15,12,9", "9,7,5,3"),
            ("24,20,16,12", "10,8,6,4"),
            ("20,16,12,8", "10,8,6,4"),
            ("20,16,12,8", "3,2,1"),
            ("16,13,10,7", "5,4,3,2,1"),
            // Deep against shallow.
            ("12,11,10,9,8,7", "6,5,4"),
            ("12,11,10,9,8,7", "6,5,4,3"),
            ("10,9,8,7,6,5", "7,6,5,4"),
            ("14,12,10,8,6", "7,5,3,1"),
            ("20,16,12,8", "6,5,4,3,2,1"),
            ("10,9,8,7,6,5,4,3,2,1", "6,4,2"),
            // Squares: only direct-versus-conjugate can differ.
            ("16,13,10,7", "16,13,10,7"),
            ("12,10,8,6", "12,10,8,6"),
            ("8,7,6,5,4,3", "8,7,6,5,4,3"),
            ("10,9,8,7,6,5", "10,9,8,7,6,5"),
            ("12,12,12,12,12,12", "12,12,12,12,12,12"),
        ]
        .iter()
        .map(|(a, b)| (parse(a), parse(b)))
        .collect()
    } else {
        explicit
            .chunks(2)
            .map(|c| (c[0].clone(), c[1].clone()))
            .collect()
    };

    println!(
        "{:<28} {:<7} {:>4} {:>5} {:>4} {:>10} {:>12} {:>10} {:>8}",
        "product", "walk", "rows", "cells", "wide", "time", "productions", "peak", "terms"
    );
    for (a, b) in &pairs {
        let walks = [
            Walk {
                name: "b",
                top: a.clone(),
                bottom: b.clone(),
                conj: false,
                block: b.clone(),
            },
            Walk {
                name: "a",
                top: b.clone(),
                bottom: a.clone(),
                conj: false,
                block: a.clone(),
            },
            Walk {
                name: "a'",
                top: a.clone(),
                bottom: b.clone(),
                conj: true,
                block: a.conjugate(),
            },
            Walk {
                name: "b'",
                top: b.clone(),
                bottom: a.clone(),
                conj: true,
                block: b.conjugate(),
            },
        ];
        // The library's own choice (`skew_lr::product_walk`): the larger
        // factor by cells, then lexicographically, on top — so the smaller is
        // enumerated — and the transpose only for a near-square. The rule is
        // restated here and checked below against the production count of
        // the walk the library actually ran, so a drift between the two is
        // an assertion failure rather than a mislabeled row.
        let a_on_top = (a.size(), a.parts()) >= (b.size(), b.parts());
        let (top, bottom) = if a_on_top { (a, b) } else { (b, a) };
        let (outer, _) = juxtapose(top, bottom);
        let cells = a.size() + b.size();
        let rule_conj = top.len() == bottom.len()
            && 10 * bottom.size() >= 7 * top.size()
            && cells >= 72
            && outer.len() >= 8
            && outer.part(0) as usize > outer.len();
        let current = match (a_on_top, rule_conj) {
            (true, false) => "b",
            (true, true) => "a'",
            (false, false) => "a",
            (false, true) => "b'",
        };
        std::env::remove_var("SKEW_ORIENT");
        clear_caches();
        take_productions();
        let library = SkewLr.schur_product(a, b);
        let library_prod = take_productions();

        let mut t_min = vec![f64::INFINITY; walks.len()];
        let mut stats = vec![(0u64, 0usize, 0usize); walks.len()];
        for r in 0..reps {
            // Rotate so each walk leads equally often and none always
            // inherits the previous walk's warm allocator.
            for k in 0..walks.len() {
                let i = (k + r) % walks.len();
                let w = &walks[i];
                let (outer, inner) = juxtapose(&w.top, &w.bottom);
                std::env::set_var("SKEW_ORIENT", if w.conj { "conj" } else { "direct" });
                clear_caches();
                take_productions();
                take_peak_layer_states();
                let t = Instant::now();
                let v = expand_skew(&outer, &inner);
                let dt = t.elapsed().as_secs_f64();
                stats[i] = (take_productions(), take_peak_layer_states(), v.len());
                t_min[i] = t_min[i].min(dt);
                assert_eq!(v, library, "walk {} disagrees on s{a}·s{b}", w.name);
            }
        }
        std::env::remove_var("SKEW_ORIENT");
        let best = t_min.iter().cloned().fold(f64::INFINITY, f64::min);
        let cur_i = walks.iter().position(|w| w.name == current).unwrap();
        assert_eq!(
            stats[cur_i].0, library_prod,
            "the restated rule does not match the walk the library ran on s{a}·s{b}"
        );
        for (i, w) in walks.iter().enumerate() {
            let (t, (prod, peak, terms)) = (t_min[i], stats[i]);
            let mark = if i == cur_i {
                format!("<- {:.2}x", t / best)
            } else {
                String::new()
            };
            let star = if (t - best).abs() < f64::EPSILON {
                "*"
            } else {
                " "
            };
            println!(
                "{:<28} {:<7} {:>4} {:>5} {:>4} {:>9.2}ms {:>12} {:>10} {:>8} {star} {mark}",
                if w.name == "b" {
                    format!("s{a}·s{b}")
                } else {
                    String::new()
                },
                w.name,
                w.block.len(),
                w.block.size(),
                w.block.part(0),
                t * 1e3,
                prod,
                peak,
                terms,
            );
        }
        println!();
    }
}
