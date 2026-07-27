//! Attack the 2-D state: per-output counting for a three-row ν.
//!
//! State after row j is (λ¹ⱼ, aⱼ, bⱼ) — the current λ¹ part, plus cells added
//! so far by strips 1 and 2. a and b are bounded by ν₁ and ν₂ rather than by
//! |μ|+ν₁, which is the only reason this is worth trying at all.
//!
//! Reports operations per output term. The bar: the frontier spends ~3.7µs per
//! term on `[20,16,12]²`, so anything past a few thousand ops per term is dead
//! regardless of constant factors.

use std::io::Write;
use std::time::Instant;

use symfn::{AutoLr, LrBackend, Partition};

macro_rules! p {
    ($($a:tt)*) => {{ println!($($a)*); std::io::stdout().flush().ok(); }};
}

/// Dense generation-stamped state table: O(1) insert, no per-row clearing.
struct Table {
    gen: Vec<u32>,
    val: Vec<u128>,
    key: Vec<u32>,
    touched: Vec<u32>,
    era: u32,
    stride_a: usize,
    stride_l: usize,
}

impl Table {
    fn new(max_l1: usize, n1: usize, n2: usize) -> Table {
        let stride_a = n2 + 1;
        let stride_l = (n1 + 1) * stride_a;
        let size = (max_l1 + 1) * stride_l;
        Table {
            gen: vec![0; size],
            val: vec![0; size],
            key: vec![0; size],
            touched: Vec::with_capacity(256),
            era: 0,
            stride_a,
            stride_l,
        }
    }
    fn clear(&mut self) {
        self.era += 1;
        self.touched.clear();
    }
    fn idx(&self, l1: usize, a: usize, b: usize) -> usize {
        l1 * self.stride_l + a * self.stride_a + b
    }
    fn add(&mut self, l1: usize, a: usize, b: usize, ways: u128) {
        let i = self.idx(l1, a, b);
        if self.gen[i] != self.era {
            self.gen[i] = self.era;
            self.val[i] = 0;
            self.key[i] = i as u32;
            self.touched.push(i as u32);
        }
        self.val[i] += ways;
    }
}

/// Count pairs (λ¹, λ²) for one λ. Returns (coefficient, transitions examined).
fn count_fibre3(
    mu: &[u32],
    lam: &[u32],
    nu: &[u32],
    cur: &mut Table,
    nxt: &mut Table,
    live: &mut Vec<(usize, usize, usize, u128)>,
) -> (u128, u64) {
    let at = |v: &[u32], i: usize| -> i64 { v.get(i).copied().map_or(0, i64::from) };
    let rows = lam.len();
    let (n1, n2) = (i64::from(nu[0]), i64::from(nu[1]));
    let max_l1 = at(lam, 0) as usize;

    cur.clear();
    cur.add(max_l1, 0, 0, 1); // λ¹ prev unbounded ⇒ max_l1 is an exact stand-in
    let mut ops = 0u64;

    let (mut big_lam, mut big_mu) = (0i64, 0i64);
    for j in 0..rows {
        big_lam += at(lam, j);
        big_mu += at(mu, j);
        let l1_lo = at(mu, j).max(at(lam, j + 2));
        let l1_hi = if j == 0 { at(lam, 0) } else { at(mu, j - 1).min(at(lam, j)) };

        live.clear();
        for &t in &cur.touched {
            let i = t as usize;
            let b = i % cur.stride_a;
            let a = (i / cur.stride_a) % (n1 as usize + 1);
            let l1 = i / cur.stride_l;
            live.push((l1, a, b, cur.val[i]));
        }
        nxt.clear();
        for &(l1_prev, a_prev, b_prev, ways) in live.iter() {
            let (a_prev, b_prev) = (a_prev as i64, b_prev as i64);
            for l1 in l1_lo..=l1_hi {
                let a = a_prev + l1 - at(mu, j);
                if a > n1 {
                    break;
                }
                let l2_lo = l1.max(at(lam, j + 1));
                let l2_hi = (l1_prev as i64).min(at(lam, j));
                for l2 in l2_lo..=l2_hi {
                    ops += 1;
                    let b = b_prev + l2 - l1;
                    if b > n2 || b > a_prev {
                        break;
                    }
                    let c = (big_lam - big_mu) - a - b;
                    if c < 0 || c > b_prev {
                        continue;
                    }
                    nxt.add(l1 as usize, a as usize, b as usize, ways);
                }
            }
        }
        std::mem::swap(cur, nxt);
        if cur.touched.is_empty() {
            return (0, ops);
        }
    }
    let mut total = 0u128;
    for &t in &cur.touched {
        let i = t as usize;
        let b = i % cur.stride_a;
        let a = (i / cur.stride_a) % (n1 as usize + 1);
        if a as i64 == n1 && b as i64 == n2 {
            total += cur.val[i];
        }
    }
    (total, ops)
}

/// Enumerate candidate λ and count each fibre — self-contained, no oracle.
/// λ ⊇ μ with λⱼ ≤ μⱼ₋₃ (through three strips) and at most three new rows.
#[allow(clippy::too_many_arguments)]
fn walk3(
    j: usize,
    rows: usize,
    left: u32,
    prev: u32,
    cand: &mut Vec<u32>,
    mu: &[u32],
    nu: &[u32],
    out: &mut Vec<(Partition, u128)>,
    ops: &mut u64,
    cands: &mut u64,
    ta: &mut Table,
    tb: &mut Table,
    live: &mut Vec<(usize, usize, usize, u128)>,
) {
    let at = |i: usize| mu.get(i).copied().unwrap_or(0);
    if j == rows {
        if left != 0 {
            return;
        }
        *cands += 1;
        let end = cand.iter().rposition(|&x| x > 0).map_or(0, |i| i + 1);
        let (c, o) = count_fibre3(mu, &cand[..end], nu, ta, tb, live);
        *ops += o;
        if c > 0 {
            out.push((Partition::new(cand[..end].iter().copied()), c));
        }
        return;
    }
    let lo = at(j);
    let hi = if j < 3 { prev.min(left + lo) } else { at(j - 3).min(prev).min(left + lo) };
    if hi < lo {
        return;
    }
    for v in lo..=hi {
        cand[j] = v;
        walk3(j + 1, rows, left - (v - lo), v, cand, mu, nu, out, ops, cands, ta, tb, live);
    }
    cand[j] = 0;
}

fn main() {
    for shape in ["4,3,2", "6,5,4", "8,6,4", "10,8,6", "12,10,8", "14,12,10", "20,16,12", "24,20,16", "30,24,18"] {
        let v: Vec<u32> = shape.split(',').map(|x| x.parse().unwrap()).collect();
        let mu = Partition::new(v.iter().copied());
        let nu = mu.clone();

        symfn::clear_caches();
        let t = Instant::now();
        let want = AutoLr.schur_product(&mu, &nu);
        let t_dp = t.elapsed();

        let mu_v: Vec<u32> = mu.parts().to_vec();
        let nu_v: Vec<u32> = nu.parts().to_vec();
        // Widest candidate: every cell of ν landing in row 0.
        let max_l1 = (mu.part(0) + nu.size()) as usize;
        let mut ta = Table::new(max_l1, nu_v[0] as usize, nu_v[1] as usize);
        let mut tb = Table::new(max_l1, nu_v[0] as usize, nu_v[1] as usize);
        let mut live = Vec::new();
        let t = Instant::now();
        let (mut ops, mut cands) = (0u64, 0u64);
        let mut got: Vec<(Partition, u128)> = Vec::new();
        let nrows = mu_v.len() + 3;
        let mut cand = vec![0u32; nrows];
        let total: u32 = nu_v.iter().sum();
        walk3(0, nrows, total, u32::MAX, &mut cand, &mu_v, &nu_v,
              &mut got, &mut ops, &mut cands, &mut ta, &mut tb, &mut live);
        got.sort_by(|x, y| x.0.cmp(&y.0));
        let t_c = t.elapsed();
        let bad = usize::from(got != want);
        p!(
            "s[{shape}]²  {:>7} terms  frontier {:>10.4?}  counting {:>10.4?}  {:>9} ops  {:>7.0} ops/term  {}",
            want.len(),
            t_dp,
            t_c,
            ops,
            ops as f64 / want.len() as f64,
            if bad == 0 { format!("ok ({} cands)", cands) } else { "MISMATCH".into() }
        );
    }
}
