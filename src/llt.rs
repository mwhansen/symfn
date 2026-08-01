//! LLT polynomials: the ribbon model, the tuple model, and the dictionary
//! between them.
//!
//! ## One family, two presentations
//!
//! **The tuple model** ([HHL] Def 3.2). A tuple `ν = (ν⁽¹⁾, …, ν⁽ʳ⁾)` of skew
//! shapes, each carrying an integer content offset; a cell's content is
//! `col − row + offset`. Reading order sorts cells by (content, component,
//! row). Two cells in *distinct* components **attack** when they share a
//! content (earlier component first) or sit on adjacent contents (later
//! component first). Then
//!
//! ```text
//!   G_ν(x; q) = Σ_{T ∈ SSYT(ν)} q^{inv(T)} x^T
//! ```
//!
//! with `inv` counting attacking pairs that are out of order.
//!
//! **The ribbon model** ([LLT] §6). Fix a level `k`. A k-ribbon tableau of
//! shape λ peels λ by horizontal k-ribbon strips; `spin` adds `(h−1)/2` per
//! ribbon and `cospin` is `s*(λ) − spin`. Four normalizations circulate, and
//! this module offers all four because the literature uses all four:
//!
//! ```text
//!   G_LT,λ  = Σ_R q^{2s(R)} x^{w(R)}                  [LT] (43)   `llt_g_lt`
//!   G̃^(k)_λ = Σ_R q^{s̃(R)} x^{w(R)}                   [LLT] (26)  `llt_gtilde`
//!   H̃^(k)_μ = G̃^(k)_{kμ}                              [LLT] (27)  `llt_h_tilde`
//!   H^(k)_μ = Σ_R q^{s(R)} x^{w(R)} = q^{s*} H̃(x;1/q) [LLT] (28)  `llt_h`
//! ```
//!
//! ## The convention minefield
//!
//! The two models agree, but not on the nose, and every trap below is silent —
//! it returns a wrong-by-a-twist answer rather than an error. Each is pinned by
//! a test.
//!
//! - **The min-inv floor.** The dictionary is
//!   `q^{−min inv} G_{k-quotient(λ)} = G̃^(k)_λ`, offsets zero. The floor is
//!   *forced*: λ = (2,2,2) at k = 2 has 2-quotient `((1), (1,1))` and
//!   `min inv = 1`, which no offset choice removes. [`llt_g`] therefore
//!   exposes the **raw** inv grading and [`llt_min_inv`] beside it, rather
//!   than quietly normalizing. λ = (2,2) is *not* a witness — its quotient
//!   `((1), (1))` has floor 0 — and reading it as one breaks the convention
//!   gate, which depends on that zero. Both halves are pinned by
//!   `the_min_inv_floor_is_forced`.
//! - **The tilde collision.** [HHL] writes `G̃` for the *spin*-flavored
//!   function; [LLT] writes `G̃` for the *cospin* one. Reading [HHL]'s
//!   remark-shaped `q^e G(1/q)` as the quotient dictionary gives the wrong
//!   answer; the direct floored equality above is the law.
//! - **Tuples lose absolute spin.** λ = (1,1,1,1) at k = 2 has one ribbon
//!   tableau with `s* = 1`, while its 2-quotient is two single cells with
//!   `max inv = min inv = 0`. No statistic of the tuple recovers `s*`, so the
//!   spin-graded `H` is exposed on the **partition-plus-level** side only.
//! - **The shuffle-world reversal** ([DA] Remark 2.2).
//!   [`SkewTuple::from_area`] lists the runs of a Dyck path in **reverse row
//!   order**; with any other order [HHL]'s attacking inversions stop agreeing
//!   with [HRW] `dinv` pointwise.
//! - **Two graphs per path.** [`DecoratedGraph::unit_interval`] (the
//!   cell-below-path rule of [CM]/[AS]) and [`DecoratedGraph::from_area`] (the
//!   dinv-faithful one) are *different graphs*: `a = (0,0)` has no cell below
//!   the path but its two rows form a primary dinv pair. Both are real objects
//!   and each has its own job here.
//!
//! ## The routes
//!
//! **R1**, [`llt_g`] — standard fillings bucketed by descent set, [HHL] (82).
//! The reference engine, and the one that emits the fundamental
//! quasisymmetric expansion ([`llt_fundamental`]) as a by-product. Cost is
//! `#SYT(ν)` for *all* monomials at once, where direct semistandard
//! enumeration pays one walk per content.
//!
//! **R2**, [`llt_gtilde`] and friends — β-number ribbon strips on the abacus.
//! The engine for the ribbon side, decomposed by runner: a bead moves `+k`, so
//! it never leaves its residue class, and a horizontal strip is a choice of
//! *prefix of a maximal block* on each runner independently. That is what
//! keeps this out of the `C(rows, m)` subset enumeration the naive reading
//! suggests — for λ = (kn) the subset count is `C(kn, n)` where the block count
//! is one.
//!
//! **R3**, [`llt_kl_column`] — [KMS] Fock-space straightening. Emits a whole
//! *column* of the Schur-expansion table (one λ, every shape μ), which is the
//! transpose of what tableau enumeration produces, and whose entries are
//! parabolic affine Kazhdan–Lusztig polynomials ([LT] Thm 4.2).
//!
//! ## Range
//!
//! Over a fixed-width `C` this family refuses rather than wrapping past its
//! wall (`docs/policies/failure.md`, R3), and **this is the one `(q,t)` family
//! in the crate whose wall a caller reaches cheaply** — so it is a measured
//! degree rather than a projection, and it is the one that earned an
//! escalation ladder. Through the Python boundary the single-shape entry
//! points escalate: the fixed-width pass reports, and the same generic code
//! re-runs over `BigInt`, so a caller there has no wall at all. A Rust caller
//! choosing `C = i128` still meets the degrees below.
//!
//! [`llt_h`] at μ = 1ⁿ overflows `i128` at the degrees below — each in about a
//! second, so the arithmetic wall arrives first and there is no runtime
//! obstacle in front of it:
//!
//! ```text
//!   k = 2   n = 124        k = 4   n = 71
//!   k = 3   n = 87
//! ```
//!
//! The wall falls as `k` rises: a larger ribbon level packs more coefficient
//! into the same degree. At k = 3 the widest coefficient is 127 bits at n = 86
//! and the next degree overflows in the accumulation, so this bounds the
//! *answers*, not merely an intermediate.
//!
//! The table entry points are stopped by runtime long before this:
//! [`llt_h_table`] gains ~2.6 bits per degree and would reach 127 bits near
//! n ≈ 54, and [`llt_gtilde_table`] ~2.5 near n ≈ 55, against tables that stop
//! finishing around n = 17–20. Range here is therefore a statement about the
//! *entry point*, not about the family.
//!
//! Degrees, slopes and the harness are in `docs/record/failure-and-overflow.md`
//! (`examples/probe_qt_walls.rs`).
//!
//! ## References
//!
//! - **[LLT]** Lascoux, Leclerc, Thibon, *Ribbon tableaux, Hall–Littlewood
//!   functions, quantum affine algebras and unipotent varieties*,
//!   [arXiv:q-alg/9512031](https://arxiv.org/abs/q-alg/9512031) — spin (15),
//!   cospin (24)–(25), `G̃` (26), `H̃` (27), `H` (28); Thm 6.6; Ex 6.8.
//! - **[LT]** Leclerc, Thibon, *Littlewood–Richardson coefficients and
//!   Kazhdan–Lusztig polynomials*,
//!   [arXiv:math/9809122](https://arxiv.org/abs/math/9809122) — the
//!   spin-generating `G` (their (43)), Ex 4.1, Lemma 6.5 (β-set strips).
//! - **[KMS]** Kashiwara, Miwa, Stern, *Decomposition of q-deformed Fock
//!   spaces*, [arXiv:q-alg/9508006](https://arxiv.org/abs/q-alg/9508006) — the
//!   **normative** straightening rules (43)/(45). [LLT] §7's printing of the
//!   same rules carries two misprints.
//! - **[HHL]** Haglund, Haiman, Loehr, *A combinatorial formula for Macdonald
//!   polynomials*, [arXiv:math/0409538](https://arxiv.org/abs/math/0409538) —
//!   Def 3.2 (the tuple model), the standardization identity (82), and
//!   `H̃_μ = Σ_D q^{−a} t^{maj} G_{ν(μ,D)}`.
//! - **[DA]** D'Adderio, *e-positivity of vertical strip LLT polynomials*,
//!   [arXiv:1906.02633](https://arxiv.org/abs/1906.02633) — Remark 2.2 (the
//!   tuple reversal), Ex 5.6.
//! - **[AS]** Alexandersson, Sulzgruber, *A combinatorial expansion of
//!   vertical-strip LLT polynomials in the basis of elementary symmetric
//!   functions*, [arXiv:2004.09198](https://arxiv.org/abs/2004.09198) — the
//!   orientation/e-expansion formula, Ex 6.1.
//! - **[CM]** Carlsson, Mellit, *A proof of the shuffle conjecture*,
//!   [arXiv:1508.06239](https://arxiv.org/abs/1508.06239) — Prop 3.5, the
//!   chromatic bridge.
//! - **[AP]** Alexandersson, Panova, *LLT polynomials, chromatic quasisymmetric
//!   functions and graphs with cycles*,
//!   [arXiv:1705.10353](https://arxiv.org/abs/1705.10353) — Def 16 (the
//!   decorated-graph presentation), Lemma 47 (the chromatic bridge, with
//!   [CM] Prop 3.5), and Conj 25 on unicellular e-positivity.
//!
//! `docs/record/llt.md` carries the measured record: the ladders, the
//! profiling, and the formula-by-formula verification against Sage.
//!
//! [LLT]: https://arxiv.org/abs/q-alg/9512031
//! [LT]: https://arxiv.org/abs/math/9809122
//! [KMS]: https://arxiv.org/abs/q-alg/9508006
//! [HHL]: https://arxiv.org/abs/math/0409538
//! [DA]: https://arxiv.org/abs/1906.02633
//! [AP]: https://arxiv.org/abs/1705.10353
//! [AS]: https://arxiv.org/abs/2004.09198
//! [CM]: https://arxiv.org/abs/1508.06239
//! [HRW]: https://arxiv.org/abs/1509.07058

// Every `as` in this module converts a *cell coordinate, component index, or
// q-exponent* — bounded by |λ|, by the number of components, and by the ribbon
// level respectively, all `u32` where they are stored. None of them carries a
// coefficient, and none can: coefficients here are `QtPoly<C>` over a generic
// `C: Ring`, and a generic parameter cannot be `as`-cast at all. That is what
// makes a module-level allow safe where a per-site one would normally be
// required — a value-carrying narrowing cannot be written in this module
// without first introducing a concrete integer coefficient
// (`docs/policies/failure.md`, R5).
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use std::collections::HashMap;

use crate::coeff::{QAlgebra, Ring};
use crate::partition::Partition;
use crate::qt::QtPoly;
use crate::sym::{Monomial, PowerSum, Schur, SymFn};

// ============================================================== the tuple model

/// One cell of a [`SkewTuple`].
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Cell {
    comp: u32,
    row: i32,
    col: i32,
    /// `col − row + offset`, the classical content.
    content: i32,
}

/// A tuple of skew shapes with content offsets — the \[HHL\] Def 3.2 object.
///
/// Components are ordered and the order matters: attacking pairs are
/// asymmetric in the component index, so permuting components changes `G_ν`.
/// Cell coordinates are signed because the ribbon components of the \[HHL\]
/// Macdonald decomposition walk left out of the first column.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SkewTuple {
    cells: Vec<Cell>,
    /// Content offset per component, kept so [`conjugate`](Self::conjugate) can
    /// rebuild.
    offsets: Vec<i32>,
    /// Where each cell sits in reading order (content ↑, component ↑, row ↑).
    reading: Vec<u32>,
    /// Attacking pairs `(a, b)`, `a` earlier than `b` in reading order.
    attack: Vec<(u32, u32)>,
    /// `T[a] ≤ T[b]`: the right neighbour inside a component.
    weak: Vec<(u32, u32)>,
    /// `T[a] < T[b]`: the neighbour one row down inside a component.
    strict: Vec<(u32, u32)>,
    /// Bit `b` of `attack_mask[a]` iff `(a, b)` attacks. The standard-filling
    /// walk's `inv` delta is one `popcount` against the assigned set, where the
    /// adjacency-list form was a pointer chase and a loop — and that walk is
    /// nearly all of the by-path shuffle refinement, so the difference is the
    /// module's headline number.
    attack_mask: Vec<u64>,
    /// Bit `p` of `pred_mask[c]` iff `p` is an immediate predecessor of `c`
    /// under weak+strict. A cell is available exactly when
    /// `pred_mask[c] & !assigned == 0`, which replaces the in-degree counters
    /// and the `avail` stack with two bit ops.
    pred_mask: Vec<u64>,
}

impl SkewTuple {
    /// From explicit cell sets: `(cells, offset)` per component, cells as
    /// `(row, col)` in any integer coordinates.
    ///
    /// The general constructor — semistandardness is read off cell *adjacency*
    /// (right neighbour weak, next row down strict), so any set of cells works
    /// and skew shapes are just the common case.
    ///
    /// # Panics
    ///
    /// If the components hold more than [`MAX_CELLS`] cells in total — the
    /// width of the `u64` attack masks, a representation limit rather than a
    /// mathematical one.
    pub fn from_cells(comps: &[(Vec<(i32, i32)>, i32)]) -> Self {
        let mut cells = Vec::new();
        let mut offsets = Vec::with_capacity(comps.len());
        for (k, (cs, off)) in comps.iter().enumerate() {
            offsets.push(*off);
            for &(row, col) in cs {
                cells.push(Cell {
                    comp: k as u32,
                    row,
                    col,
                    content: col - row + off,
                });
            }
        }
        let n = cells.len();

        // Reading order. Ties in content are broken by component then row, and
        // that tie-break is what makes [HHL]'s standardization an identity
        // rather than a convention — see `syt_buckets`.
        let mut order: Vec<u32> = (0..n as u32).collect();
        order.sort_by_key(|&i| {
            let c = cells[i as usize];
            (c.content, c.comp, c.row)
        });
        let mut reading = vec![0u32; n];
        for (pos, &i) in order.iter().enumerate() {
            reading[i as usize] = pos as u32;
        }

        // Attacking pairs, between distinct components only.
        //
        // The bit masks below cap a tuple at 64 cells. That is the same cap
        // `syt_buckets`'s descent mask already imposed, so it is not a new
        // limit — but it is now enforced at construction, and says so.
        assert!(
            n <= MAX_CELLS,
            "a SkewTuple holds at most {MAX_CELLS} cells, got {n}; this is a \
             representation limit, not a mathematical one"
        );
        let mut attack = Vec::new();
        let mut attack_mask = vec![0u64; n];
        for a in 0..n {
            for b in 0..n {
                let (x, y) = (cells[a], cells[b]);
                let attacks = (y.content == x.content && x.comp < y.comp)
                    || (y.content == x.content + 1 && y.comp < x.comp);
                if attacks {
                    attack.push((a as u32, b as u32));
                    attack_mask[a] |= 1u64 << b;
                }
            }
        }

        // Semistandardness covers, by adjacency inside a component.
        let index: HashMap<(u32, i32, i32), u32> = cells
            .iter()
            .enumerate()
            .map(|(i, c)| ((c.comp, c.row, c.col), i as u32))
            .collect();
        let mut weak = Vec::new();
        let mut strict = Vec::new();
        let mut pred_mask = vec![0u64; n];
        for (i, c) in cells.iter().enumerate() {
            if let Some(&j) = index.get(&(c.comp, c.row, c.col + 1)) {
                weak.push((i as u32, j));
                pred_mask[j as usize] |= 1u64 << i;
            }
            if let Some(&j) = index.get(&(c.comp, c.row + 1, c.col)) {
                strict.push((i as u32, j));
                pred_mask[j as usize] |= 1u64 << i;
            }
        }

        SkewTuple {
            cells,
            offsets,
            reading,
            attack,
            weak,
            strict,
            attack_mask,
            pred_mask,
        }
    }

    /// From skew shapes `outer/inner`, one content offset each.
    ///
    /// # Panics
    ///
    /// If `skews` and `offsets` differ in length, if any `inner ⊄ outer`, or if
    /// the shapes hold more than [`MAX_CELLS`] cells between them.
    pub fn from_skews(skews: &[(Partition, Partition)], offsets: &[i32]) -> Self {
        assert_eq!(skews.len(), offsets.len(), "one offset per component");
        let comps: Vec<(Vec<(i32, i32)>, i32)> = skews
            .iter()
            .zip(offsets)
            .map(|((outer, inner), &off)| {
                assert!(
                    outer.contains(inner),
                    "a skew shape needs inner ⊆ outer: {inner} ⊄ {outer}"
                );
                let cells = (0..outer.len())
                    .flat_map(|r| (inner.part(r)..outer.part(r)).map(move |c| (r as i32, c as i32)))
                    .collect();
                (cells, off)
            })
            .collect();
        Self::from_cells(&comps)
    }

    /// From straight shapes, one content offset each.
    pub fn from_partitions(shapes: &[Partition], offsets: &[i32]) -> Self {
        let skews: Vec<(Partition, Partition)> = shapes
            .iter()
            .map(|s| (s.clone(), Partition::default()))
            .collect();
        Self::from_skews(&skews, offsets)
    }

    /// The k-quotient of λ with **offsets zero** — the verified normal form of
    /// the quotient dictionary.
    ///
    /// `None` when λ has no k-ribbon tableaux (nonempty k-core), which is
    /// exactly when the dictionary has nothing to say. Note what the caller
    /// still owes: `G` of this tuple carries a forced `q^{min inv}` floor that
    /// `G̃^(k)_λ` does not — see [`llt_min_inv`] and the module docs.
    pub fn quotient(lambda: &Partition, k: u32) -> Option<Self> {
        if !lambda.has_empty_k_core(k) {
            return None;
        }
        let shapes = lambda.k_quotient(k);
        let offsets = vec![0i32; shapes.len()];
        Some(Self::from_partitions(&shapes, &offsets))
    }

    /// The vertical-strip tuple of a Dyck path given by its area sequence.
    ///
    /// Components are the maximal `+1`-runs of the area sequence, each a single
    /// column; row `i`'s cell sits at content `−a_i`. Components are listed in
    /// **reverse row order** (\[DA\] Remark 2.2) — the one order under which
    /// \[HHL\] attacking inversions and \[HRW\] `dinv` agree pointwise, not
    /// merely in total.
    pub fn from_area(area: &[u32]) -> Self {
        let n = area.len();
        let mut comps: Vec<(Vec<(i32, i32)>, i32)> = Vec::new();
        let mut start = 0usize;
        for i in 1..=n {
            if i == n || area[i] != area[i - 1] + 1 {
                let cells = (start..i).map(|r| ((r - start) as i32, 0i32)).collect();
                comps.push((cells, -(area[start] as i32)));
                start = i;
            }
        }
        comps.reverse();
        Self::from_cells(&comps)
    }

    /// The conjugate tuple: transpose every component, negate every offset,
    /// reverse the component order.
    ///
    /// The shape the ω-duality law wants,
    /// `G_{ν'}(x; q) = q^{A(ν)} ω G_ν(x; 1/q)` with `A(ν)` the attacking-pair
    /// count. Offset negation is the law that keeps contents negated
    /// cell-for-cell; it is verified at offset 0 (see the ω-duality test) and
    /// is the only extension of that consistent with `content ↦ −content`.
    pub fn conjugate(&self) -> Self {
        let ncomps = self.offsets.len();
        let mut comps: Vec<(Vec<(i32, i32)>, i32)> = (0..ncomps)
            .map(|k| (Vec::new(), -self.offsets[k]))
            .collect();
        for c in &self.cells {
            comps[c.comp as usize].0.push((c.col, c.row));
        }
        comps.reverse();
        Self::from_cells(&comps)
    }

    /// `A(ν)`, the number of attacking pairs — the twist in the ω-duality law.
    pub fn attacking_pairs(&self) -> usize {
        self.attack.len()
    }

    /// The number of cells, i.e. the degree of `G_ν`.
    pub fn size(&self) -> usize {
        self.cells.len()
    }

    /// The number of components.
    pub fn components(&self) -> usize {
        self.offsets.len()
    }

    /// Standard fillings bucketed by descent set: `(mask, counts by inv)`.
    ///
    /// \[HHL\] (82) is `G_ν = Σ_{S ∈ SYT(ν)} q^{inv(S)} Q_{n,D(S)}`, so these
    /// buckets *are* the fundamental quasisymmetric expansion. Bit `i−1` of
    /// `mask` is set iff `i` is a descent, meaning `i+1` sits earlier than `i`
    /// in reading order.
    ///
    /// ## Why this is written in bit masks
    ///
    /// This walk is the module's hot spot: sampling `nabla_e_by_path(9)` puts
    /// nearly every sample in it, with the area enumeration, the tuple
    /// construction and the allocator sharing what is left
    /// (`docs/record/llt.md`). So the inner loop is where the module's headline
    /// number lives, and it is written accordingly:
    ///
    /// - **Availability** is `pred_mask[c] & !assigned == 0`. The textbook form
    ///   keeps in-degree counters and an `avail` stack, and pays a decrement
    ///   loop plus a stack splice on the way down *and* on the way back up; two
    ///   bit ops replace all of it, and nothing is allocated or undone.
    /// - **The `inv` delta** is `(attack_mask[c] & assigned).count_ones()`.
    ///    Assigning the value `s` to `c` can only complete attacking pairs `(c,
    ///    b)` with `b` already filled — an unfilled `b` gets a larger value and
    ///    is counted from its own side — so the delta is a popcount where the
    ///    adjacency-list form was a pointer chase and a loop.
    /// - **The descent bit** compares against the previous value's position
    ///   only, so it rides down the recursion as a `u32` instead of an array.
    ///
    /// The leaf goes through [`Sink`], monomorphized, so that the flat-table
    /// and hash-map accumulators below cost the same as writing either one
    /// inline.
    fn syt_buckets(&self) -> Vec<(u64, Vec<u128>)> {
        self.syt_buckets_within(FLAT_TABLE_BUDGET)
    }

    /// [`syt_buckets`](Self::syt_buckets) with the flat table's size budget
    /// exposed, so a test can drive the same tuple down both sinks.
    fn syt_buckets_within(&self, budget: usize) -> Vec<(u64, Vec<u128>)> {
        let n = self.cells.len();
        if n == 0 {
            return vec![(0u64, vec![1u128])];
        }
        let width = self.attack.len() + 1; // inv ≤ #attacking pairs
        let masks = 1usize << (n - 1);
        // A hash per leaf was 6% of the R1 profile, and there are `#SYT` leaves.
        // Direct-indexing the descent mask removes it outright, but the table is
        // `2^{n−1} × width`, so it is only taken when that fits a fixed budget —
        // past which the map is the honest fallback rather than a memory cliff.
        // Saturating is the *correct* comparison here, not a concession: a
        // product that overflows `usize` is one that exceeds the budget, and
        // saturation routes it to the fallback exactly as the true value would
        // (R4). `checked_mul` would say the same thing more loudly for no gain.
        let mut out = if masks.saturating_mul(width) <= budget {
            let mut sink = FlatSink {
                table: vec![0u128; masks * width],
                width,
            };
            self.walk(n, &mut sink);
            sink.table
                .chunks(width)
                .enumerate()
                .filter(|(_, row)| row.iter().any(|&c| c != 0))
                .map(|(mask, row)| (mask as u64, row.to_vec()))
                .collect()
        } else {
            let mut sink = MapSink {
                buckets: Default::default(),
            };
            self.walk(n, &mut sink);
            sink.buckets.into_iter().collect::<Vec<_>>()
        };
        out.sort_unstable_by_key(|e| e.0);
        out
    }

    fn walk<S: Sink>(&self, n: usize, sink: &mut S) {
        let all = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
        self.walk_rec(1, n, all, 0, 0, 0, 0, sink);
    }

    #[allow(clippy::too_many_arguments)]
    fn walk_rec<S: Sink>(
        &self,
        step: usize,
        n: usize,
        all: u64,
        assigned: u64,
        prev_pos: u32,
        inv: u32,
        mask: u64,
        sink: &mut S,
    ) {
        if step > n {
            sink.hit(mask, inv);
            return;
        }
        let mut cand = all & !assigned;
        while cand != 0 {
            let c = cand.trailing_zeros();
            cand &= cand - 1;
            let cell = c as usize;
            if self.pred_mask[cell] & !assigned != 0 {
                continue;
            }
            let delta = (self.attack_mask[cell] & assigned).count_ones();
            let pos = self.reading[cell];
            let bit = if step >= 2 && pos < prev_pos {
                1u64 << (step - 2)
            } else {
                0
            };
            self.walk_rec(
                step + 1,
                n,
                all,
                assigned | 1u64 << c,
                pos,
                inv + delta,
                mask | bit,
                sink,
            );
        }
    }

    /// `(min inv, max inv)` over standard fillings, without building the
    /// buckets.
    ///
    /// The floor and ceiling are wanted far more often than the whole
    /// expansion — the quotient dictionary needs the floor on every shape it
    /// sweeps — and the walk that produces them need not accumulate anything.
    fn inv_range(&self) -> (u32, u32) {
        if self.cells.is_empty() {
            return (0, 0);
        }
        let mut sink = RangeSink {
            lo: u32::MAX,
            hi: 0,
        };
        self.walk(self.cells.len(), &mut sink);
        (sink.lo.min(sink.hi), sink.hi)
    }
}

/// Where the standard-filling walk deposits a leaf. Monomorphized per
/// accumulator, so the walk pays no dispatch for having three of them.
trait Sink {
    fn hit(&mut self, mask: u64, inv: u32);
}

/// Direct-indexed `(mask, inv)` counts — one add, no hashing.
struct FlatSink {
    table: Vec<u128>,
    width: usize,
}

impl Sink for FlatSink {
    #[inline]
    fn hit(&mut self, mask: u64, inv: u32) {
        self.table[mask as usize * self.width + inv as usize] += 1;
    }
}

/// The largest flat descent-mask table [`SkewTuple::syt_buckets`] will build,
/// in `u128` entries. Past it the map fallback takes over.
const FLAT_TABLE_BUDGET: usize = 1 << 20;

/// The fallback for tuples whose descent-mask table would not fit.
///
/// [`hit`](Sink::hit) runs once per standard filling, so the per-key hash the
/// flat table exists to avoid is paid on every leaf here — and the key is a
/// bare `u64`, which is what [`crate::fasthash`] is for — SipHash was a third
/// of the runtime on the tuples that land here (`docs/record/llt.md`).
struct MapSink {
    buckets: crate::fasthash::Map<u64, Vec<u128>>,
}

impl Sink for MapSink {
    fn hit(&mut self, mask: u64, inv: u32) {
        let counts = self.buckets.entry(mask).or_default();
        if counts.len() <= inv as usize {
            counts.resize(inv as usize + 1, 0);
        }
        counts[inv as usize] += 1;
    }
}

/// Just the extremes, for [`llt_min_inv`] / [`llt_max_inv`].
struct RangeSink {
    lo: u32,
    hi: u32,
}

impl Sink for RangeSink {
    #[inline]
    fn hit(&mut self, _mask: u64, inv: u32) {
        self.lo = self.lo.min(inv);
        self.hi = self.hi.max(inv);
    }
}

/// The mask of descent positions a monomial `x^μ` tolerates: bit `s−1` for each
/// proper partial sum `s` of μ.
///
/// `[x^μ] Q_{n,D} = 1` exactly when `D` is contained in this set, so the
/// monomial expansion is a subset filter over the buckets.
fn partial_sum_mask(mu: &Partition) -> u64 {
    let mut mask = 0u64;
    let mut acc = 0u32;
    for &p in mu.parts().iter().take(mu.len().saturating_sub(1)) {
        acc += p;
        mask |= 1u64 << (acc - 1);
    }
    mask
}

/// `G_ν(x; q)` in the monomial basis, in the **raw** inv grading.
///
/// The floor is not divided out: `min_T inv(T)` can be positive and callers
/// that want the cospin normalization must divide by `q^{llt_min_inv(ν)}`
/// themselves. That is deliberate — the floor is real data about ν (see the
/// module docs), and hiding it is how the quotient dictionary gets misread.
pub fn llt_g<C: Ring>(nu: &SkewTuple) -> Monomial<QtPoly<C>> {
    let buckets = nu.syt_buckets();
    let mut out = Monomial::zero();
    for mu in crate::partitions_of(nu.size() as u32) {
        let allowed = partial_sum_mask(&mu);
        let mut poly = QtPoly::zero();
        for (mask, counts) in &buckets {
            if mask & !allowed != 0 {
                continue;
            }
            for (inv, c) in counts.iter().enumerate() {
                if *c != 0 {
                    poly.add_term(inv as u32, 0, C::from_u128(*c));
                }
            }
        }
        out.add_term(mu, poly);
    }
    out
}

/// `min_T inv(T)` over semistandard fillings — the forced `q`-floor of the
/// tuple model.
///
/// Read off the standard fillings, which is legitimate because standardization
/// preserves `inv` (\[HHL\] (82)): every semistandard filling has a standard
/// one with the same `inv`, and standard fillings are themselves semistandard,
/// so the two minima coincide.
pub fn llt_min_inv(nu: &SkewTuple) -> u32 {
    nu.inv_range().0
}

/// `max_T inv(T)`, by the same argument as [`llt_min_inv`].
pub fn llt_max_inv(nu: &SkewTuple) -> u32 {
    nu.inv_range().1
}

/// The **fundamental quasisymmetric** expansion of `G_ν`, as
/// `(composition, coefficient)` pairs.
///
/// The buckets of \[HHL\] (82) read directly: a descent set `D ⊆ [n−1]` is the
/// composition of `n` whose partial sums are `D`. The crate has no QSym type
/// and does not need one for this — the compositions carry their own meaning
/// and nothing here multiplies them.
///
/// Nobody ships this expansion, which is the point of exposing it.
pub fn llt_fundamental<C: Ring>(nu: &SkewTuple) -> Vec<(Vec<u32>, QtPoly<C>)> {
    let n = nu.size() as u32;
    let mut out = Vec::new();
    for (mask, counts) in nu.syt_buckets() {
        let mut comp = Vec::new();
        let mut prev = 0u32;
        for i in 1..n {
            if mask >> (i - 1) & 1 == 1 {
                comp.push(i - prev);
                prev = i;
            }
        }
        if n > 0 {
            comp.push(n - prev);
        }
        let mut poly = QtPoly::zero();
        for (inv, c) in counts.iter().enumerate() {
            if *c != 0 {
                poly.add_term(inv as u32, 0, C::from_u128(*c));
            }
        }
        out.push((comp, poly));
    }
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

// ============================================================= the ribbon model

/// An abacus state: bit `b` is set iff `b` is a β-number of the **conjugate**
/// shape.
///
/// `u128` because the whole-degree work at `k·n ≈ 60` needs β-numbers up to
/// `ℓ(λ) + λ₁`, and because every operation the strip enumeration performs —
/// occupancy, the `+k` move, the crossing count — is then a bit op.
type Abacus = u128;

/// The highest β-number this representation holds.
const ABACUS_BITS: u32 = 128;

/// The abacus of λ with `rows` beads.
fn abacus_of(lambda: &Partition, rows: usize) -> Abacus {
    let conj = lambda.conjugate();
    let beta = conj.beta_numbers(rows);
    let top = beta.first().copied().unwrap_or(0);
    assert!(
        top < ABACUS_BITS,
        "β-number {top} does not fit a {ABACUS_BITS}-bit abacus"
    );
    beta.iter().fold(0u128, |acc, &b| acc | 1u128 << b)
}

/// The shape whose conjugate's β-numbers are the set bits of `beta`.
fn shape_of(beta: Abacus) -> Partition {
    let rows = beta.count_ones() as usize;
    let mut desc = Vec::with_capacity(rows);
    let mut b = beta;
    while b != 0 {
        let top = ABACUS_BITS - 1 - b.leading_zeros();
        desc.push(top);
        b &= !(1u128 << top);
    }
    Partition::from_beta_numbers(&desc).conjugate()
}

/// Is the shape of `beta` contained in the shape of `target`?
///
/// Containment of partitions is elementwise domination of β-numbers at equal
/// bead counts, and conjugation preserves containment, so this reads the two
/// bitmasks without materializing either partition. The strip walk calls it
/// once per candidate strip, which makes it the pruning hot path.
///
/// **Only the differing bits are visited.** For equal-size bead sets,
/// `β(ν)_j ≤ β(λ)_j` for every `j` is the same as
/// `|ν ∩ [x,∞)| ≤ |λ ∩ [x,∞)|` for every threshold `x`, and that count
/// difference only changes at a position where the two disagree. So walking
/// `beta ^ target` from the top with a running `+1` for target-only and `−1`
/// for beta-only, and failing the moment it goes negative, decides containment
/// in `popcount(beta ^ target)` steps — two per moved bead — where marching
/// both bead lists in lockstep took one step per bead. At k = 3, n = 13 that is
/// 2 iterations instead of 39 for a weight-1 strip.
fn abacus_contained(beta: Abacus, target: Abacus) -> bool {
    debug_assert_eq!(beta.count_ones(), target.count_ones());
    let mut d = beta ^ target;
    let mut acc = 0i32;
    while d != 0 {
        let b = ABACUS_BITS - 1 - d.leading_zeros();
        d &= !(1u128 << b);
        if target >> b & 1 == 1 {
            acc += 1;
        } else {
            acc -= 1;
            if acc < 0 {
                return false;
            }
        }
    }
    true
}

/// Reusable buffers for the strip enumeration.
///
/// [`for_each_strip_up`] runs once per (state, weight) in the ribbon walk, and
/// allocating its runner, block and suffix-sum vectors per call put `malloc`
/// and `free` between them above `strip_rec` itself in the R2 profile
/// (`docs/record/llt.md`). Blocks are `(offset, len)` into one flat bead vector
/// rather than a `Vec<Vec<u32>>`, so a call touches three buffers and allocates
/// nothing.
#[derive(Default)]
struct StripScratch {
    /// Bead positions grouped by runner, ascending within each runner.
    beads: Vec<u32>,
    /// `(offset into beads, length)` per maximal block.
    blocks: Vec<(u32, u32)>,
    /// Suffix sums of block lengths, for the weight pruning.
    room: Vec<u32>,
    /// `bits ≡ r (mod k)` per runner, rebuilt only when `k` changes — which is
    /// once per walk, since a walk has one level.
    runner_masks: Vec<Abacus>,
    /// The `k` `runner_masks` was built for. `0` is never a valid level, so a
    /// default scratch always rebuilds.
    runner_k: u32,
}

/// Every horizontal k-ribbon strip of weight `m` that can be **added** to
/// `beta`, as `(new state, spin_LT)`.
///
/// `spin_LT` is the integral spin `Σ (h−1) = 2·s` of \[LT\], which is what
/// keeps the walk in ℤ; the halvings to \[LLT\]'s `s` and `s̃` happen once at
/// the end, in [`regrade`].
///
/// The enumeration is decomposed by runner (\[LT\] Lemma 6.5 read on the
/// abacus). A bead moves `+k`, so it stays in its residue class; within a
/// runner, moving a set of beads up one slot each keeps them distinct exactly
/// when the set is a **prefix of a maximal block** of occupied slots — move any
/// bead with an occupied slot above it and the two collide, and a maximal
/// block's top always has room. So a strip is one prefix length per block, and
/// the weights add. The naive reading of the same lemma enumerates all
/// `C(rows, m)` subsets and filters; for λ = (kn) that is `C(kn, n)` candidates
/// where this is one.
///
/// Crossings — the unmoved beads a moving bead passes — come out as `k−1`
/// popcounts: an unmoved bead lies strictly between `B` and `B+k` iff it sits
/// at `B+i` for some `0 < i < k`.
// Reached only from the tests since the R2 profiling work moved the live path
// elsewhere, and kept for two reasons: it is a second, independent enumerator
// for the strip lemma, and its doc carries the block-prefix derivation that
// keeps this out of the `C(rows, m)` subset enumeration.
#[allow(dead_code)]
fn for_each_strip_up(
    beta: Abacus,
    k: u32,
    m: u32,
    sc: &mut StripScratch,
    visit: &mut impl FnMut(Abacus, u32),
) {
    collect_blocks(beta, k, sc);
    // How many beads the remaining blocks can still supply, for pruning.
    sc.room.clear();
    sc.room.resize(sc.blocks.len() + 1, 0);
    for i in (0..sc.blocks.len()).rev() {
        sc.room[i] = sc.room[i + 1] + sc.blocks[i].1;
    }
    strip_rec(beta, k, &sc.beads, &sc.blocks, &sc.room, 0, m, 0, m, visit);
}

/// The maximal blocks of `beta`, per runner, into `sc`.
///
/// Driven off `beta & runner_mask[r]`, so it visits the `rows` occupied
/// positions rather than every position up to the top bead. The abacus runs at
/// roughly half density, and this was 7% of the R2 profile.
fn collect_blocks(beta: Abacus, k: u32, sc: &mut StripScratch) {
    if sc.runner_k != k {
        sc.runner_masks.clear();
        for r in 0..k {
            let mut mask: Abacus = 0;
            let mut g = r;
            while g < ABACUS_BITS {
                mask |= 1u128 << g;
                g += k;
            }
            sc.runner_masks.push(mask);
        }
        sc.runner_k = k;
    }
    sc.beads.clear();
    sc.blocks.clear();
    for r in 0..k as usize {
        let mut b = beta & sc.runner_masks[r];
        let mut start = sc.beads.len() as u32;
        let mut len = 0u32;
        let mut prev: Option<u32> = None;
        while b != 0 {
            let g = b.trailing_zeros();
            b &= b - 1;
            if prev.is_some_and(|p| g != p + k) {
                sc.blocks.push((start, len));
                start = sc.beads.len() as u32;
                len = 0;
            }
            sc.beads.push(g);
            len += 1;
            prev = Some(g);
        }
        if len > 0 {
            sc.blocks.push((start, len));
        }
    }
}

/// Every horizontal k-ribbon strip from `beta` of **every** weight `1 … max_w`,
/// as `(weight, new state, spin_LT)`.
///
/// The same block tree as [`for_each_strip_up`], but the weight falls out at
/// the leaf instead of being fixed in advance. That matters because the ribbon
/// walk wants every weight from each state: asking one weight at a time
/// re-collected the blocks and re-walked all the shared internal nodes, once
/// per weight. `strip_rec` was 70% of the R2 profile once the allocator was
/// dealt with, and this is where most of it came from.
fn for_each_strip_any(
    beta: Abacus,
    k: u32,
    max_w: u32,
    sc: &mut StripScratch,
    visit: &mut impl FnMut(u32, Abacus, u32),
) {
    collect_blocks(beta, k, sc);
    strip_any_rec(beta, k, &sc.beads, &sc.blocks, 0, 0, max_w, 0, visit);
}

#[allow(clippy::too_many_arguments)]
fn strip_any_rec(
    beta: Abacus,
    k: u32,
    beads: &[u32],
    blocks: &[(u32, u32)],
    bi: usize,
    w: u32,
    max_w: u32,
    moved: Abacus,
    visit: &mut impl FnMut(u32, Abacus, u32),
) {
    if bi == blocks.len() {
        if w >= 1 {
            let unmoved = beta & !moved;
            let shifted = moved << k;
            debug_assert_eq!(shifted & unmoved, 0, "a block prefix never collides");
            let mut cross = 0u32;
            for i in 1..k {
                cross += ((moved << i) & unmoved).count_ones();
            }
            visit(w, unmoved | shifted, (k - 1) * w - cross);
        }
        return;
    }
    let (off, len) = blocks[bi];
    let top = len.min(max_w - w);
    let mut mv = moved;
    for j in 0..=top {
        if j > 0 {
            mv |= 1u128 << beads[(off + len - j) as usize];
        }
        strip_any_rec(beta, k, beads, blocks, bi + 1, w + j, max_w, mv, visit);
    }
}

#[allow(clippy::too_many_arguments)]
// The recursion behind `for_each_strip_up`, and dead exactly when it is.
#[allow(dead_code)]
fn strip_rec(
    beta: Abacus,
    k: u32,
    beads: &[u32],
    blocks: &[(u32, u32)],
    room: &[u32],
    bi: usize,
    remaining: u32,
    moved: Abacus,
    m: u32,
    visit: &mut impl FnMut(Abacus, u32),
) {
    if remaining > room[bi] {
        return;
    }
    if bi == blocks.len() {
        if remaining == 0 {
            let unmoved = beta & !moved;
            let shifted = moved << k;
            debug_assert_eq!(shifted & unmoved, 0, "a block prefix never collides");
            let mut cross = 0u32;
            for i in 1..k {
                cross += ((moved << i) & unmoved).count_ones();
            }
            visit(unmoved | shifted, (k - 1) * m - cross);
        }
        return;
    }
    let (off, len) = blocks[bi];
    let top = len.min(remaining);
    let mut mv = moved;
    for j in 0..=top {
        if j > 0 {
            // The top `j` beads of a block are its last `j` entries, and the
            // loop adds them one at a time from the top down.
            mv |= 1u128 << beads[(off + len - j) as usize];
        }
        strip_rec(
            beta,
            k,
            beads,
            blocks,
            room,
            bi + 1,
            remaining - j,
            mv,
            m,
            visit,
        );
    }
}

/// One node of the weight walk: which shapes are reachable, and with what spin.
///
/// [`crate::fasthash`]'s hasher rather than SipHash: keys are a single small
/// structured integer and there are a great many of them, which is the exact
/// case that module was extracted for.
type States<C> = crate::fasthash::Map<Abacus, QtPoly<C>>;

/// The k-ribbon tableau walk, in the \[LT\] spin grading.
///
/// Builds up from the empty shape by horizontal strips of **weakly decreasing**
/// weight, so the weight word of a leaf is a partition and the recursion tree
/// is the partition trie — every weight *prefix* is walked once and shared by
/// all the partitions extending it. That sharing is the whole point of the
/// whole-degree entry points: `μ = (n)` is 94–100% of a degree in the incumbent
/// (`docs/record/llt.md`) precisely because per-shape morphisms
/// recompute what a trie shares.
///
/// `target` prunes to chains through subshapes of one λ, which is exact —
/// a chain to λ never leaves `⊆ λ`.
struct RibbonWalk<'a, C: Ring> {
    k: u32,
    target: Option<Abacus>,
    weight: Vec<u32>,
    out: HashMap<Partition, Monomial<QtPoly<C>>>,
    /// The shape to keep, when only one is wanted.
    keep: Option<&'a Partition>,
}

impl<C: Ring> RibbonWalk<'_, C> {
    fn step(&mut self, states: &States<C>, remaining: u32, max_w: u32, sc: &mut StripScratch) {
        if remaining == 0 {
            let weight = Partition::from_sorted(self.weight.clone());
            for (state, poly) in states {
                let shape = shape_of(*state);
                if let Some(keep) = self.keep {
                    if &shape != keep {
                        continue;
                    }
                }
                self.out
                    .entry(shape)
                    .or_insert_with(Monomial::zero)
                    .add_term(weight.clone(), poly.clone());
            }
            return;
        }
        // Copied out so the visit closure below does not borrow `self`, which
        // the scratch buffer and the recursive call both need.
        let (k, target) = (self.k, self.target);
        // Every next weight from one walk per state, bucketed by weight, rather
        // than one walk per (state, weight) — see `for_each_strip_any`. The cost
        // is holding `wmax` layer maps at once instead of one; the levels
        // where `wmax` is large are the shallow ones, which have few states, so
        // the peak is set by the deep narrow levels either way.
        let wmax = max_w.min(remaining);
        let mut next: Vec<States<C>> = vec![States::default(); wmax as usize + 1];
        for (state, poly) in states {
            for_each_strip_any(*state, k, wmax, sc, &mut |w, nb, spin| {
                if let Some(t) = target {
                    if !abacus_contained(nb, t) {
                        return;
                    }
                }
                next[w as usize]
                    .entry(nb)
                    .or_insert_with(QtPoly::zero)
                    .add_shifted(poly, (spin, 0), false);
            });
        }
        for w in (1..=wmax).rev() {
            if next[w as usize].is_empty() {
                continue;
            }
            self.weight.push(w);
            let level = std::mem::take(&mut next[w as usize]);
            self.step(&level, remaining - w, w, sc);
            self.weight.pop();
        }
    }
}

/// The highest β-number a walk can reach, asserted to fit.
///
/// `β_j = conj(ν)_j + rows − 1 − j`, so the top β is `ℓ(ν) + rows − 1`, and a
/// move adds `k` on top of that. The two entry points know different things
/// about `ℓ(ν)`:
///
/// - the **pruned** walk stays inside one λ, so `ℓ(ν) ≤ ℓ(λ)` and the bound is
///   `ℓ(λ) + λ₁ + k` — for a one-row λ that is barely more than `λ₁`;
/// - the **unpruned** walk visits every shape of the degree, where the only
///   bound on `ℓ(ν)` is the size.
///
/// Worth separating rather than taking the loose bound for both: `H^(5)_{(13)}`
/// is a 65-cell one-row shape needing 71 bits, and the loose bound rejects it
/// at 135.
fn assert_abacus_fits(top_beta: u32, what: &str) {
    assert!(
        top_beta < ABACUS_BITS,
        "{what} needs β-numbers up to {top_beta}, past this abacus's \
         {ABACUS_BITS} bits"
    );
}

/// The top β-number the **pruned** walk on λ at level `k` will need.
///
/// Public because the Python boundary refuses past this wall rather than
/// letting the assert above become a `PanicException`
/// (`docs/policies/failure.md`, R2). Exposing the expression rather than
/// restating it there is the point: a bound copied to a second site is a bound
/// that drifts, and this one is already subtle enough to need the note above
/// about why the two walks cannot share it.
pub fn abacus_reach(lambda: &Partition, k: u32) -> u32 {
    lambda.len() as u32 + lambda.part(0).max(1) + k
}

/// The same for the **unpruned** whole-degree walk, which visits every shape of
/// size `k·n` and so can only use the size as its bound on `ℓ(ν)`.
pub fn abacus_reach_table(n: u32, k: u32) -> u32 {
    let size = k * n;
    size + size + k
}

/// The width the two functions above are measured against.
pub const ABACUS_REACH_LIMIT: u32 = ABACUS_BITS;

/// The most cells a [`SkewTuple`] can hold, a representation limit of the
/// `u64` attack masks rather than a mathematical one.
pub const MAX_CELLS: usize = 64;

/// The edges [`llt_e_expansion`] orients itself: the weak ones, as unordered
/// pairs, that are not also strict.
///
/// The sum it drives is `2^{free.len()}` terms, so the count is what a caller
/// needs to know before committing — and what the Python boundary refuses on,
/// rather than letting the assert below reach Sage as a `PanicException`.
pub fn free_edges(g: &DecoratedGraph) -> Vec<(u32, u32)> {
    let mut free: Vec<(u32, u32)> = g
        .weak
        .iter()
        .map(|&(u, v)| (u.min(v), u.max(v)))
        .filter(|e| !g.strict.contains(e))
        .collect();
    free.sort_unstable();
    free.dedup();
    free
}

/// The ceiling on [`free_edges`], set by the `u32` orientation mask.
pub const MAX_FREE_EDGES: usize = 32;

/// `Σ_R q^{spin_LT(R)} x^{w(R)}` over k-ribbon tableaux of shape λ — the \[LT\]
/// (43) grading `q^{2s}`.
///
/// The rawest of the four normalizations, and the one the Fock route
/// ([`llt_kl_column`]) is pinned against. Zero when λ has no k-ribbon tableaux.
///
/// # Panics
///
/// If `k == 0`, and if the abacus this shape needs is wider than
/// [`ABACUS_REACH_LIMIT`] — [`abacus_reach`] is that requirement, exposed so a
/// caller can ask before committing rather than discover it here.
///
/// A λ with no k-ribbon tableaux is an *answer* (zero), not a panic.
pub fn llt_g_lt<C: Ring>(lambda: &Partition, k: u32) -> Monomial<QtPoly<C>> {
    assert!(k >= 1, "a ribbon level needs k ≥ 1");
    let n = lambda.size();
    if !n.is_multiple_of(k) || !lambda.has_empty_k_core(k) {
        return Monomial::zero();
    }
    if n == 0 {
        return Monomial::monomial(Partition::default(), <QtPoly<C> as Ring>::one());
    }
    let rows = lambda.part(0).max(1) as usize;
    assert_abacus_fits(
        abacus_reach(lambda, k),
        &format!("the shape {lambda} at level {k}"),
    );
    let target = abacus_of(lambda, rows);
    let empty = abacus_of(&Partition::default(), rows);
    let mut walk = RibbonWalk::<C> {
        k,
        target: Some(target),
        weight: Vec::new(),
        out: HashMap::new(),
        keep: Some(lambda),
    };
    let mut start: States<C> = States::default();
    start.insert(empty, <QtPoly<C> as Ring>::one());
    let r = n / k;
    walk.step(&start, r, r, &mut StripScratch::default());
    walk.out.remove(lambda).unwrap_or_else(Monomial::zero)
}

/// Which of \[LLT\]'s halvings to apply to a `spin_LT`-graded expansion.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Grading {
    /// `s̃ = (s* − spin)/2` — \[LLT\] (24)–(26)'s cospin.
    Cospin,
    /// `s = spin/2` — \[LLT\] (28)'s spin.
    Spin,
}

/// Halve the `spin_LT` grading into \[LLT\]'s `s` or `s̃`.
///
/// Integrality is an assertion, not a hope: cospin is an integer on every
/// shape (\[LLT\] (25)), and spin is an integer on the shapes `kμ` that `H`
/// is defined on but *not* in general — a shape with odd total `spin_LT` has
/// half-integral spin and \[LLT\] (28) carries a `q^{1/2}`. Asking for
/// [`Grading::Spin`] there is a caller error and says so.
fn regrade<C: Ring>(f: &Monomial<QtPoly<C>>, how: Grading) -> Monomial<QtPoly<C>> {
    let smax = f
        .terms()
        .values()
        .filter_map(|p| p.degrees())
        .map(|(a, _)| a)
        .max()
        .unwrap_or(0);
    let mut out = Monomial::zero();
    for (mu, p) in f.terms() {
        let mut np = QtPoly::zero();
        for (&(a, b), c) in p.terms() {
            let e = match how {
                Grading::Cospin => smax - a,
                Grading::Spin => a,
            };
            assert!(
                e % 2 == 0,
                "{how:?} of {mu} is half-integral (spin_LT {a}, s* {smax}): this \
                 shape's grading carries a q^(1/2)"
            );
            np.add_term(e / 2, b, c.clone());
        }
        out.add_term(mu.clone(), np);
    }
    out
}

/// `G̃^(k)_λ(x; q)`, the **cospin** generating function of \[LLT\] (26), in the
/// monomial basis.
///
/// Zero when λ has no k-ribbon tableaux. `k = 1` gives `s_λ`; a shape `kμ`
/// gives `H̃^(k)_μ` (\[LLT\] (27), and [`llt_h_tilde`]).
pub fn llt_gtilde<C: Ring>(lambda: &Partition, k: u32) -> Monomial<QtPoly<C>> {
    let raw = llt_g_lt::<C>(lambda, k);
    if raw.is_zero() {
        return raw;
    }
    regrade(&raw, Grading::Cospin)
}

/// `H̃^(k)_μ = G̃^(k)_{kμ}` (\[LLT\] (27)).
pub fn llt_h_tilde<C: Ring>(mu: &Partition, k: u32) -> Monomial<QtPoly<C>> {
    llt_gtilde(&scale(mu, k), k)
}

/// `H^(k)_μ = Σ_R q^{s(R)} x^{w(R)} = q^{s*} H̃^(k)_μ(x; 1/q)` (\[LLT\] (28)).
///
/// Exposed on the **partition-plus-level** side only, and that is a
/// mathematical fact rather than an API choice: the k-quotient tuple of a shape
/// does not determine `s*` (λ = (1,1,1,1) at k = 2 is the smallest witness), so
/// there is no honest `H` of a bare tuple.
///
/// `H^(1)_μ = s_μ`, and for `k` past \[LLT\] Thm 6.6's bound this is the
/// Hall–Littlewood `Q'_μ` — the k-interpolation between the two.
pub fn llt_h<C: Ring>(mu: &Partition, k: u32) -> Monomial<QtPoly<C>> {
    let raw = llt_g_lt::<C>(&scale(mu, k), k);
    if raw.is_zero() {
        return raw;
    }
    regrade(&raw, Grading::Spin)
}

/// `kμ = (kμ₁, kμ₂, …)`.
fn scale(mu: &Partition, k: u32) -> Partition {
    Partition::from_sorted(mu.parts().iter().map(|&p| k * p).collect())
}

/// `H^(k)_μ` for **every** μ ⊢ n — the unit the incumbent's walls are measured
/// in (`docs/record/llt.md`).
pub fn llt_h_table<C: Ring>(n: u32, k: u32) -> Vec<(Partition, Monomial<QtPoly<C>>)> {
    crate::partitions_of(n)
        .into_iter()
        .map(|mu| {
            let f = llt_h::<C>(&mu, k);
            (mu, f)
        })
        .collect()
}

/// `G̃^(k)_λ` for **every** λ ⊢ kn with empty k-core, from a single walk.
///
/// The unpruned form of [`llt_gtilde`]: one pass over the weight trie visits
/// every shape of the degree at once, where the pruned walk pays per shape.
/// Which is cheaper depends on the degree — the shapes of `k·n` outnumber the
/// `p(n)` shapes `kμ` that [`llt_h_table`] wants, but they share every chain —
/// so both exist and `the_two_walks_agree` holds them together.
///
/// # Panics
///
/// If `k == 0`, or if the whole-degree abacus exceeds
/// [`ABACUS_REACH_LIMIT`]. This walk is unpruned and so can only bound `ℓ(ν)`
/// by the degree, which is why [`abacus_reach_table`] is a separate and much
/// coarser bound than [`abacus_reach`].
pub fn llt_gtilde_table<C: Ring>(n: u32, k: u32) -> Vec<(Partition, Monomial<QtPoly<C>>)> {
    assert!(k >= 1, "a ribbon level needs k ≥ 1");
    if n == 0 {
        return vec![(
            Partition::default(),
            Monomial::monomial(Partition::default(), <QtPoly<C> as Ring>::one()),
        )];
    }
    let size = k * n;
    let rows = size as usize;
    assert_abacus_fits(
        abacus_reach_table(n, k),
        &format!("every shape of size {size} at level {k}"),
    );
    let empty = abacus_of(&Partition::default(), rows);
    let mut walk = RibbonWalk::<C> {
        k,
        target: None,
        weight: Vec::new(),
        out: HashMap::new(),
        keep: None,
    };
    let mut start: States<C> = States::default();
    start.insert(empty, <QtPoly<C> as Ring>::one());
    walk.step(&start, n, n, &mut StripScratch::default());
    let mut out: Vec<(Partition, Monomial<QtPoly<C>>)> = walk
        .out
        .into_iter()
        .map(|(shape, f)| {
            let g = regrade(&f, Grading::Cospin);
            (shape, g)
        })
        .collect();
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

/// `G̃^(k)_λ` in the **Schur** basis.
///
/// The row of the Schur-expansion table at λ. [`llt_kl_column`] produces the
/// transposed object — one Schur index, every shape — and the two are worth
/// having separately because the sweeps want different slices.
pub fn llt_schur<C: Ring>(lambda: &Partition, k: u32) -> Schur<QtPoly<C>> {
    use crate::convert::ToSchur;
    llt_gtilde::<C>(lambda, k).to_schur()
}

// ================================================= the shuffle-side refinement

/// `∇e_n = Σ_D t^{area(D)} G_D(x; q)`: the by-path Schur refinement of the
/// shuffle theorem, as `(area sequence, G_D)` pairs.
///
/// Every `G_D` is Schur-positive, so this is `∇e_n` written as a positive sum
/// of positive pieces — the decomposition `docs/record/dyck-paths.md` calls
/// "the win left on the table" and that no package emits.
///
/// The two obstructions recorded there are gone. Standardizing labelled paths
/// has no dinv-invariant tie-break, but standardizing *tuple fillings* does,
/// and \[HHL\] (82) makes it an identity rather than a convention — so each
/// `G_D` costs `#SYT` of its tuple instead of one enumeration per content. The
/// valley side of the Delta conjecture stays out: `Val` is not an LLT
/// statistic, and [`crate::dyck`] keeps it.
pub fn nabla_e_by_path<C: Ring>(n: u32) -> Vec<(Vec<u32>, Monomial<QtPoly<C>>)> {
    let mut out = Vec::new();
    for_each_area(n as usize, &mut |area| {
        let g = llt_g::<C>(&SkewTuple::from_area(area));
        let areasum: u32 = area.iter().sum();
        let mut scaled = Monomial::zero();
        for (mu, poly) in g.terms() {
            scaled.add_term(mu.clone(), poly.shift_t(areasum));
        }
        out.push((area.to_vec(), scaled));
    });
    out
}

/// Every area sequence of size `n`: `a₁ = 0` and `a_i ≤ a_{i−1} + 1`.
///
/// The same enumeration as [`crate::dyck`]'s, kept local because this module
/// needs the sequence and that one needs it labelled.
fn for_each_area(n: usize, visit: &mut impl FnMut(&[u32])) {
    fn rec(n: usize, area: &mut Vec<u32>, visit: &mut impl FnMut(&[u32])) {
        if area.len() == n {
            visit(area);
            return;
        }
        let top = if area.is_empty() {
            0
        } else {
            area[area.len() - 1] + 1
        };
        for a in 0..=top {
            area.push(a);
            rec(n, area, visit);
            area.pop();
        }
    }
    let mut area = Vec::with_capacity(n);
    rec(n, &mut area, visit);
}

// ======================================================== the graph presentation

/// A decorated graph: the coloring presentation of a vertical-strip LLT
/// polynomial (\[AP\] Def 16, \[DA\] §2.2).
///
/// Vertices are `0 … n−1`. A **weak** edge is an ordered pair `(u, v)` whose
/// *ascents* are counted: a coloring κ contributes `q` for it when
/// `κ(u) < κ(v)`. A **strict** edge `(u, v)` — always with `u < v` — instead
/// *constrains*: `κ(u) < κ(v)` or the coloring does not count at all. Strict
/// edges never contribute to the statistic.
///
/// The two constructors build genuinely different graphs from the same Dyck
/// path and the module docs say why; do not reach for one expecting the other.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct DecoratedGraph {
    n: u32,
    weak: Vec<(u32, u32)>,
    strict: Vec<(u32, u32)>,
}

impl DecoratedGraph {
    /// From explicit edge sets. Strict edges must run `u < v`, and the two sets
    /// must be disjoint as unordered pairs.
    ///
    /// Disjointness is checked because it is the invariant [`llt_graph`] relies
    /// on to skip the strict edges when counting ascents. An edge that is both
    /// would be constrained *and* counted, and the statistic would silently
    /// gain a `q` per strict edge.
    ///
    /// # Panics
    ///
    /// If a strict edge is not oriented `u < v`, if any endpoint is outside
    /// `0..n`, or if an edge appears in both sets as an unordered pair.
    pub fn new(n: u32, weak: &[(u32, u32)], strict: &[(u32, u32)]) -> Self {
        assert!(
            strict.iter().all(|&(u, v)| u < v),
            "a strict edge is oriented u < v: κ(u) < κ(v) is the constraint"
        );
        assert!(
            weak.iter().chain(strict).all(|&(u, v)| u < n && v < n),
            "edges must land inside 0..n"
        );
        let unordered = |&(u, v): &(u32, u32)| (u.min(v), u.max(v));
        assert!(
            weak.iter()
                .map(unordered)
                .all(|e| !strict.iter().map(unordered).any(|s| s == e)),
            "an edge is either weak (counted) or strict (constraining), never both"
        );
        DecoratedGraph {
            n,
            weak: weak.to_vec(),
            strict: strict.to_vec(),
        }
    }

    /// The **dinv-faithful** graph of a Dyck path.
    ///
    /// Weak edges are the \[HRW\] dinv pairs, oriented so that an ascent *is* a
    /// dinv: primary (`a_x = a_y`, `x < y`) as `(x, y)`, secondary
    /// (`a_x = a_y + 1`, `x < y`) as `(y, x)`. Strict edges are consecutive
    /// same-run rows, where labels must strictly increase up a column. The two
    /// sets are disjoint: a same-run pair is never a dinv pair.
    ///
    /// This is *not* [`unit_interval`](Self::unit_interval). `a = (0,0)` has no
    /// cell below its path, and so no edge there, but its two rows do form a
    /// primary dinv pair.
    pub fn from_area(area: &[u32]) -> Self {
        let n = area.len();
        let mut weak = Vec::new();
        for x in 0..n {
            for y in (x + 1)..n {
                if area[x] == area[y] {
                    weak.push((x as u32, y as u32));
                } else if area[x] == area[y] + 1 {
                    weak.push((y as u32, x as u32));
                }
            }
        }
        let strict: Vec<(u32, u32)> = (1..n)
            .filter(|&i| area[i] == area[i - 1] + 1)
            .map(|i| ((i - 1) as u32, i as u32))
            .collect();
        DecoratedGraph {
            n: n as u32,
            weak,
            strict,
        }
    }

    /// The **unit interval graph** of a Dyck path by the cell-below-path rule
    /// of \[CM\] §3 / \[AS\] §2.6: an edge `{j, i}`, `j < i`, exactly when
    /// `j ≥ i − a_i`.
    ///
    /// No strict edges — this is the unicellular case. It is the presentation
    /// the chromatic bridge consumes, and the generator of the unit-interval
    /// corpus; see [`from_area`](Self::from_area) for the other graph.
    pub fn unit_interval(area: &[u32]) -> Self {
        let n = area.len();
        let mut weak = Vec::new();
        for i in 0..n {
            for j in (i.saturating_sub(area[i] as usize))..i {
                weak.push((j as u32, i as u32));
            }
        }
        DecoratedGraph {
            n: n as u32,
            weak,
            strict: Vec::new(),
        }
    }

    /// The number of vertices; a coloring assigns a value to each of them.
    pub fn order(&self) -> u32 {
        self.n
    }

    /// Every coloring of content μ obeying the strict constraints, visited with
    /// its ascent count.
    fn for_each_coloring(&self, mu: &Partition, visit: &mut impl FnMut(&[u32])) {
        let n = self.n as usize;
        let nvals = mu.len();
        let mut counts: Vec<u32> = mu.parts().to_vec();
        let mut kappa: Vec<u32> = Vec::with_capacity(n);
        self.color_rec(n, nvals, &mut counts, &mut kappa, visit);
    }

    fn color_rec(
        &self,
        n: usize,
        nvals: usize,
        counts: &mut [u32],
        kappa: &mut Vec<u32>,
        visit: &mut impl FnMut(&[u32]),
    ) {
        let v = kappa.len();
        if v == n {
            visit(kappa);
            return;
        }
        for j in 0..nvals {
            if counts[j] == 0 {
                continue;
            }
            let c = j as u32 + 1;
            // Strict edges run u < v, so the constraint on the vertex being
            // assigned is always against an already-coloured endpoint.
            let ok = self
                .strict
                .iter()
                .all(|&(a, b)| b as usize != v || kappa[a as usize] < c);
            if !ok {
                continue;
            }
            counts[j] -= 1;
            kappa.push(c);
            self.color_rec(n, nvals, counts, kappa, visit);
            kappa.pop();
            counts[j] += 1;
        }
    }

    fn ascents(&self, kappa: &[u32]) -> u32 {
        self.weak
            .iter()
            .filter(|&&(a, b)| kappa[a as usize] < kappa[b as usize])
            .count() as u32
    }
}

/// `G_Γ(x; q) = Σ_κ q^{asc(κ)} x^κ` over the colorings of a decorated graph,
/// in the monomial basis.
///
/// The third presentation of the vertical-strip LLTs, and the one the chromatic
/// bridge and the \[AS\] e-expansion both speak. Equal to
/// `llt_g(SkewTuple::from_area(a))` when the graph is
/// [`DecoratedGraph::from_area(a)`](DecoratedGraph::from_area).
pub fn llt_graph<C: Ring>(g: &DecoratedGraph) -> Monomial<QtPoly<C>> {
    let mut out = Monomial::zero();
    for mu in crate::partitions_of(g.n) {
        let mut poly = QtPoly::zero();
        g.for_each_coloring(&mu, &mut |kappa| {
            poly.add_term(g.ascents(kappa), 0, C::one());
        });
        out.add_term(mu, poly);
    }
    out
}

/// `X_Γ(x; q) = (q−1)^{−n} G_Γ[x(q−1); q]` — the Shareshian–Wachs chromatic
/// quasisymmetric function of Γ, from its LLT polynomial (\[CM\] Prop 3.5,
/// \[AP\] Lemma 47).
///
/// The plethysm is a p-basis twist, `p_r ↦ (q^r − 1) p_r`, and the division by
/// `(q−1)^n` is exact — a remainder would mean the input was not a coloring
/// generating function of the right degree.
///
/// Γ should carry the \[CM\] presentation — natural orientation, no strict
/// edges, as [`DecoratedGraph::unit_interval`] builds — for the output to be
/// the chromatic function of the graph rather than of a decorated relative of
/// it. Isolated vertices are part of the graph and must be counted in `n`;
/// dropping them is a classic way to get a plausible wrong answer.
///
/// # Panics
///
/// If `(q−1)ⁿ` does not divide the twisted expansion exactly. That is not a
/// capacity wall — it says the input was not a coloring generating function of
/// degree `n`, which for the presentations named above cannot happen, so
/// reaching it means Γ was built some other way.
pub fn chromatic_from_llt<C: QAlgebra>(g: &DecoratedGraph) -> Monomial<QtPoly<C>> {
    use crate::convert::{FromSchur, ToSchur};
    let n = g.n;
    let f: Monomial<QtPoly<C>> = llt_graph(g);
    let p: PowerSum<QtPoly<C>> = PowerSum::from_schur(&f.to_schur());
    let mut twisted: PowerSum<QtPoly<C>> = PowerSum::zero();
    for (rho, c) in p.terms() {
        let mut factor = <QtPoly<C> as Ring>::one();
        for &r in rho.parts() {
            // q^r − 1
            let mut term: QtPoly<C> = QtPoly::term(r, 0, C::one());
            term.add_term(0, 0, C::one().neg());
            factor = factor.mul(&term);
        }
        twisted.add_term(rho.clone(), c.mul(&factor));
    }
    let mut divisor = <QtPoly<C> as Ring>::one();
    let mut qm1: QtPoly<C> = QtPoly::term(1, 0, C::one());
    qm1.add_term(0, 0, C::one().neg());
    for _ in 0..n {
        divisor = divisor.mul(&qm1);
    }
    let mut out: PowerSum<QtPoly<C>> = PowerSum::zero();
    for (rho, c) in twisted.terms() {
        let q = c
            .divide_exact(&divisor)
            .expect("(q−1)^n must divide the twisted expansion exactly");
        out.add_term(rho.clone(), q);
    }
    Monomial::from_schur(&out.to_schur())
}

/// The \[AS\] **e-expansion** of `Ĝ_Γ(x; q+1)`: `Σ_θ q^{asc(θ)} e_{λ(θ)}` over
/// orientations of the free (non-strict) edges.
///
/// Blocks come from the highest vertex reachable along strict and *ascending*
/// edges; the block sizes are the partition. By \[DA\]'s theorem the
/// coefficients are non-negative, so this is **certified positive output**
/// rather than a conjecture to check — which is why the test asserts positivity
/// here and only *records* it for the Schur side of [`nabla_e_by_path`].
///
/// Weak edges are read as unordered pairs `{min, max}`: \[AS\]'s formula
/// orients the edges itself, so the input's orientation is not consulted. Cost
/// is `2^{#free edges}`.
///
/// # Panics
///
/// If [`free_edges`] returns [`MAX_FREE_EDGES`] or more, the ceiling set by the
/// `u32` orientation mask. Since the sum is `2^{free}` terms, a caller wanting
/// to know before committing asks `free_edges(g).len()` — which is why that
/// function is public.
pub fn llt_e_expansion<C: Ring>(g: &DecoratedGraph) -> Vec<(Partition, QtPoly<C>)> {
    let n = g.n as usize;
    let strict: Vec<(u32, u32)> = g.strict.clone();
    let free = free_edges(g);
    assert!(
        free.len() < MAX_FREE_EDGES,
        "the orientation sum is 2^{} terms",
        free.len()
    );

    let mut acc: HashMap<Partition, QtPoly<C>> = HashMap::new();
    for mask in 0u32..(1u32 << free.len()) {
        let mut adj: Vec<Vec<u32>> = vec![Vec::new(); n];
        for &(u, v) in &strict {
            adj[u as usize].push(v);
        }
        let mut asc = 0u32;
        for (i, &(u, v)) in free.iter().enumerate() {
            if mask >> i & 1 == 1 {
                adj[u as usize].push(v);
                asc += 1;
            }
        }
        // Highest reachable vertex along strict + ascending edges.
        let mut blocks: HashMap<u32, u32> = HashMap::new();
        for u in 0..n {
            let mut seen = vec![false; n];
            let mut stack = vec![u as u32];
            seen[u] = true;
            let mut best = u as u32;
            while let Some(w) = stack.pop() {
                best = best.max(w);
                for &z in &adj[w as usize] {
                    if !seen[z as usize] {
                        seen[z as usize] = true;
                        stack.push(z);
                    }
                }
            }
            *blocks.entry(best).or_insert(0) += 1;
        }
        let lambda = Partition::new(blocks.values().copied());
        acc.entry(lambda)
            .or_insert_with(QtPoly::zero)
            .add_term(asc, 0, C::one());
    }
    let mut out: Vec<(Partition, QtPoly<C>)> = acc.into_iter().collect();
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

// =============================================== the [HHL] Macdonald assembly

/// `H̃_μ(x; q, t) = Σ_D q^{−a(D)} t^{maj(D)} G_{ν(μ,D)}(x; q)` — the \[HHL\]
/// decomposition, in the monomial basis.
///
/// The fourth route to `H̃` in this crate, and the only one that is positively
/// graded at every intermediate step: Macdonald positivity *is* LLT
/// positivity, and this is where that becomes a computation rather than a
/// slogan. The other three live in [`crate::qtkostka`] and must agree with it.
///
/// `D` ranges over subsets of the non-bottom cells of μ, so this is
/// `2^{|μ|−μ₁}` LLT evaluations — a reference route, not a fast one.
///
/// `a(D)` is negative in the exponent, and [`QtPoly`] exponents are unsigned,
/// so the accumulation carries a uniform `q^{Σ arms}` and divides it out at the
/// end. The division is exact; a remainder would mean the offset was too small.
///
/// # Panics
///
/// If that division leaves a remainder — a bug in the offset, not an input the
/// caller can pick. Also on the [`MAX_CELLS`] wall, since each `D` builds a
/// [`SkewTuple`] holding `|μ|` cells.
pub fn htilde_by_llt<C: Ring>(mu: &Partition) -> Monomial<QtPoly<C>> {
    let k = mu.part(0);
    if k == 0 {
        return Monomial::monomial(Partition::default(), <QtPoly<C> as Ring>::one());
    }
    let conj = mu.conjugate();
    // Candidate descent cells (i, j), 1-based, every cell not in the bottom row
    // of its column.
    let cand: Vec<(u32, u32)> = (1..=k)
        .flat_map(|j| (2..=conj.part((j - 1) as usize)).map(move |i| (i, j)))
        .collect();
    let arm = |(i, j): (u32, u32)| mu.part((i - 1) as usize) - j;
    let leg = |(i, j): (u32, u32)| conj.part((j - 1) as usize) - i;
    let offset: u32 = cand.iter().map(|&c| arm(c)).sum();

    let mut acc: Monomial<QtPoly<C>> = Monomial::zero();
    for mask in 0u64..(1u64 << cand.len()) {
        let d: Vec<(u32, u32)> = (0..cand.len())
            .filter(|i| mask >> i & 1 == 1)
            .map(|i| cand[i])
            .collect();
        let a: u32 = d.iter().map(|&c| arm(c)).sum();
        let maj: u32 = d.iter().map(|&c| leg(c) + 1).sum();

        let comps: Vec<(Vec<(i32, i32)>, i32)> = (1..=k)
            .map(|j| {
                let size = conj.part((j - 1) as usize);
                let des: Vec<u32> = d
                    .iter()
                    .filter(|&&(_, jj)| jj == j)
                    .map(|&(i, _)| i)
                    .collect();
                (ribbon_from_descents(size, &des), 0)
            })
            .collect();
        let g = llt_g::<C>(&SkewTuple::from_cells(&comps));
        for (nu, poly) in g.terms() {
            acc.terms_mut()
                .entry(nu.clone())
                .or_insert_with(<QtPoly<C> as Ring>::zero)
                .add_shifted(poly, (offset - a, maj), false);
        }
    }
    acc.terms_mut().retain(|_, c| !c.is_zero());

    let divisor: QtPoly<C> = QtPoly::term(offset, 0, C::one());
    let mut out = Monomial::zero();
    for (nu, poly) in acc.terms() {
        out.add_term(
            nu.clone(),
            poly.divide_exact(&divisor)
                .expect("the q-offset must divide out exactly"),
        );
    }
    out
}

/// The ribbon component of the \[HHL\] decomposition: `size` cells on
/// consecutive contents, turning a corner at every descent.
///
/// Cell contents run `−1, −2, …` downward in this crate's convention (`col −
/// row`), which is [HHL]'s `row − col` negated; the descent set says which
/// steps go up rather than left.
fn ribbon_from_descents(size: u32, des: &[u32]) -> Vec<(i32, i32)> {
    if size == 0 {
        return Vec::new();
    }
    let mut cells = vec![(1i32, 0i32)];
    for c in 2..=size {
        let (r, cc) = *cells.last().expect("cells is seeded with the first cell");
        if des.contains(&c) {
            cells.push((r + 1, cc));
        } else {
            cells.push((r, cc - 1));
        }
    }
    cells
}

// ================================================ R3: the Fock straightening

/// A wedge: a sequence of integers, strictly decreasing once straightened.
type Wedge = Vec<i32>;

/// Normal-ordered wedges with their coefficients in the \[KMS\] variable `v`.
type Straightened<C> = crate::fasthash::Map<Wedge, QtPoly<C>>;

/// `−v · p`, in one pass: shift the `v` exponents by one and negate.
///
/// `p.mul(&QtPoly::q()).neg()` is the same value through two allocations and
/// the general merge, and this is the inner loop of the straightening ladder —
/// which a sampling profile put at 63% allocator before this and the in-place
/// wedge walk.
fn neg_v_times<C: Ring>(p: &QtPoly<C>) -> QtPoly<C> {
    let mut out = <QtPoly<C> as Ring>::zero();
    out.add_shifted(p, (1, 0), true);
    out
}

/// `(v² − 1) · p`, in one pass, for the same reason.
fn v2_minus_one_times<C: Ring>(p: &QtPoly<C>) -> QtPoly<C> {
    let mut out = <QtPoly<C> as Ring>::zero();
    out.add_shifted(p, (2, 0), false);
    out.add_shifted(p, (0, 0), true);
    out
}

/// \[KMS\] (43)/(45): normal-order a wedge, accumulating into `out`.
///
/// The rules, with `i = (m−l) mod k`:
///
/// ```text
///   u_l ∧ u_m = −u_m ∧ u_l                                  when l ≡ m (mod k)
///   u_l ∧ u_m = −v u_m∧u_l + (v²−1)( u_{m−i}∧u_{l+i}
///               − v u_{m−k}∧u_{l+k} + v² u_{m−k−i}∧u_{l+k+i} + ⋯ )   otherwise
/// ```
///
/// offsets `i, k, k+i, 2k, 2k+i, …` truncated at `2a < m − l`. \[LLT\] §7
/// prints the same rules with two misprints; \[KMS\] is what this implements,
/// and the dictionary to the ribbon side is `q = −v` — pinned by the multi-term
/// Kazhdan–Lusztig entries, which single monomials cannot separate from the
/// mirrored laws.
///
/// **Written in place.** Each branch changes only the two positions `i, i+1`,
/// so it mutates, recurses and restores rather than cloning the wedge. Cloning
/// per branch made this route **63% allocator** in a sampling profile —
/// `malloc` and `free` together outweighed the straightening itself by three to
/// one — and the offset ladder is generated on the fly for the same reason: the
/// candidate offsets `i, k, k+i, 2k, …` are already sorted and already
/// distinct, so building, sorting and filtering a `Vec` of them per ascent
/// bought nothing.
fn straighten<C: Ring>(w: &mut [i32], coeff: &QtPoly<C>, k: u32, out: &mut Straightened<C>) {
    if coeff.is_zero() {
        return;
    }
    for i in 0..w.len().saturating_sub(1) {
        if w[i] == w[i + 1] {
            return; // a repeated index kills the wedge
        }
        if w[i] < w[i + 1] {
            let (l, m) = (w[i], w[i + 1]);
            let diff = (m - l) as u32;
            w[i] = m;
            w[i + 1] = l;
            if diff.is_multiple_of(k) {
                straighten(w, &coeff.neg(), k, out);
            } else {
                straighten(w, &neg_v_times(coeff), k, out);
                let mut term = v2_minus_one_times(coeff);
                let im = diff % k;
                let mut a = im;
                // The offsets i, k, k+i, 2k, 2k+i, … in order, without a Vec:
                // already sorted and already distinct, so building one bought
                // nothing.
                while 2 * a < diff {
                    w[i] = m - a as i32;
                    w[i + 1] = l + a as i32;
                    straighten(w, &term, k, out);
                    term = neg_v_times(&term); // the (−v)^t ladder
                    a = if a.is_multiple_of(k) {
                        a + im
                    } else {
                        a - im + k
                    };
                }
            }
            w[i] = l;
            w[i + 1] = m;
            return;
        }
    }
    // Most leaves land on a wedge already seen — that accumulation is the whole
    // point of the route — so probe before paying for the key.
    if let Some(slot) = out.get_mut(&*w) {
        slot.add_assign(coeff);
    } else {
        out.insert(w.to_vec(), coeff.clone());
    }
}

/// `V_m`, the `h_m`-shaped boson: add `k` to a multiset of `m` positions, then
/// straighten.
///
/// The multiset walk is incremental for the same reason [`straighten`] is: one
/// buffer, `+k` on the way down and `−k` on the way back up, instead of a fresh
/// wedge per multiset.
fn v_op<C: Ring>(wedges: &Straightened<C>, m: u32, k: u32) -> Straightened<C> {
    let mut out: Straightened<C> = Straightened::default();
    let mut buf: Wedge = Vec::new();
    for (wedge, coeff) in wedges {
        buf.clear();
        buf.extend_from_slice(wedge);
        boson_rec(&mut buf, 0, m as usize, k, coeff, &mut out);
    }
    out.retain(|_, c| !c.is_zero());
    out
}

/// Every `m`-multiset of positions, applied to `buf` in place.
fn boson_rec<C: Ring>(
    buf: &mut [i32],
    lo: usize,
    m: usize,
    k: u32,
    coeff: &QtPoly<C>,
    out: &mut Straightened<C>,
) {
    if m == 0 {
        straighten(buf, coeff, k, out);
        return;
    }
    for j in lo..buf.len() {
        buf[j] += k as i32;
        boson_rec(buf, j, m - 1, k, coeff, out);
        buf[j] -= k as i32;
    }
}

/// One **column** of the Schur-expansion table: `⟨μ + ρ| S_λ |ρ⟩` for every
/// shape μ ⊢ k|λ|, in the \[KMS\] variable `v`.
///
/// The transpose of what tableau enumeration gives. One straightening run at
/// fixed λ produces the coefficient of `s_λ` in `G_LT,μ` for *every* μ at once,
/// and those coefficients are parabolic affine Kazhdan–Lusztig polynomials
/// (\[LT\] Thm 4.2, Varagnolo–Vasserot) — computed by exact straightening with
/// no Hecke algebra in sight.
///
/// The variable is `v`, occupying [`QtPoly`]'s `q` slot; the ribbon side's
/// grading is recovered at **`q = −v`**. That dictionary is the one thing here
/// no single monomial can pin, and the module's test range includes the
/// multi-term degree-4 entries that do.
///
/// # Panics
///
/// If `k == 0`, and on the same abacus wall as [`llt_g_lt`]
/// ([`abacus_reach_table`] against [`ABACUS_REACH_LIMIT`]).
pub fn llt_kl_column<C: Ring>(lambda: &Partition, k: u32) -> Vec<(Partition, QtPoly<C>)> {
    use crate::convert::FromSchur;
    assert!(k >= 1, "a ribbon level needs k ≥ 1");
    let deg = lambda.size();
    if deg == 0 {
        return vec![(Partition::default(), <QtPoly<C> as Ring>::one())];
    }
    let rows = (k * deg) as usize;
    let rho: Wedge = (0..rows as i32).rev().collect();

    // κ: s_λ in the h-basis, one determinant row per λ.
    let s: Schur<i64> = Schur::monomial(lambda.clone(), 1);
    let h = crate::sym::Homogeneous::from_schur(&s);

    let mut acc: HashMap<Wedge, QtPoly<C>> = HashMap::new();
    for (nu, c) in h.terms() {
        let mut wedges: Straightened<C> = Straightened::default();
        wedges.insert(rho.clone(), QtPoly::term(0, 0, C::from_i64(*c)));
        for &part in nu.parts() {
            wedges = v_op(&wedges, part, k);
        }
        for (w, cc) in wedges {
            acc.entry(w)
                .or_insert_with(<QtPoly<C> as Ring>::zero)
                .add_assign(&cc);
        }
    }

    let mut out = Vec::new();
    for mu in crate::partitions_of(k * deg) {
        if mu.len() > rows {
            continue;
        }
        let key: Wedge = (0..rows).map(|j| mu.part(j) as i32 + rho[j]).collect();
        if let Some(c) = acc.get(&key) {
            if !c.is_zero() {
                out.push((mu, c.clone()));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;
    use crate::convert::{FromSchur, ToSchur};
    use crate::sym::Elementary;

    type Q = QtPoly<i64>;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    fn q_pow(e: u32) -> Q {
        QtPoly::term(e, 0, 1)
    }

    fn poly(terms: &[(u32, i64)]) -> Q {
        let mut p = QtPoly::zero();
        for &(e, c) in terms {
            p.add_term(e, 0, c);
        }
        p
    }

    fn tuple(shapes: &[&[u32]]) -> SkewTuple {
        let shapes: Vec<Partition> = shapes.iter().map(|s| part(s)).collect();
        let offsets = vec![0i32; shapes.len()];
        SkewTuple::from_partitions(&shapes, &offsets)
    }

    /// The map fallback is only reached by tuples too large to run in a test
    /// suite, so the budget is driven to zero instead: same tuples, same
    /// answers, the other sink. Without this the fallback's only coverage is
    /// the cases nobody can afford to check.
    ///
    /// Compared on trimmed rows, because the two sinks disagree on trailing
    /// zeros and always have: the flat table pads every row to `A(ν) + 1`,
    /// while the map grows a row only to the largest `inv` it saw. Both
    /// consumers skip zero counts, so the length is not part of the value —
    /// which is exactly the claim being pinned here.
    #[test]
    fn both_syt_bucket_sinks_agree_up_to_trailing_zeros() {
        fn trim(b: Vec<(u64, Vec<u128>)>) -> Vec<(u64, Vec<u128>)> {
            b.into_iter()
                .map(|(mask, mut counts)| {
                    while counts.last() == Some(&0) {
                        counts.pop();
                    }
                    (mask, counts)
                })
                .collect()
        }
        for shapes in [
            &[&[2u32, 1][..], &[2, 1][..]][..],
            &[&[3, 2][..], &[2][..]][..],
            &[&[2, 2][..], &[2, 1][..], &[1][..]][..],
            &[&[3, 1, 1][..]][..],
        ] {
            let nu = tuple(shapes);
            assert_eq!(
                trim(nu.syt_buckets()),
                trim(nu.syt_buckets_within(0)),
                "sinks disagree on {shapes:?}"
            );
        }
    }

    /// Direct semistandard enumeration, sharing nothing with the descent
    /// buckets: fill cells in index order, check the covers, read off `x^μ`.
    fn g_direct(nu: &SkewTuple) -> Monomial<Q> {
        let n = nu.size();
        let mut out = Monomial::zero();
        for mu in crate::partitions_of(n as u32) {
            let nvals = mu.len();
            let mut counts: Vec<u32> = mu.parts().to_vec();
            let mut filling: Vec<u32> = Vec::with_capacity(n);
            let mut acc = QtPoly::zero();
            fill_rec(nu, n, nvals, &mut counts, &mut filling, &mut acc);
            out.add_term(mu, acc);
        }
        out
    }

    fn fill_rec(
        nu: &SkewTuple,
        n: usize,
        nvals: usize,
        counts: &mut [u32],
        filling: &mut Vec<u32>,
        acc: &mut Q,
    ) {
        let cell = filling.len();
        if cell == n {
            let inv = nu
                .attack
                .iter()
                .filter(|&&(a, b)| filling[a as usize] > filling[b as usize])
                .count() as u32;
            acc.add_term(inv, 0, 1);
            return;
        }
        for j in 0..nvals {
            if counts[j] == 0 {
                continue;
            }
            let v = j as u32 + 1;
            // Each disjunct is one way the edge is *violated*, so the negation
            // sits outside the whole search rather than on each clause: "no
            // weak edge objects to v here". De Morgan of the per-clause form.
            let ok_weak = !nu.weak.iter().any(|&(x, y)| {
                let (x, y) = (x as usize, y as usize);
                (y == cell && filling[x] > v) || (x == cell && y < cell && v > filling[y])
            });
            let ok_strict = !nu.strict.iter().any(|&(x, y)| {
                let (x, y) = (x as usize, y as usize);
                (y == cell && filling[x] >= v) || (x == cell && y < cell && v >= filling[y])
            });
            if ok_weak && ok_strict {
                counts[j] -= 1;
                filling.push(v);
                fill_rec(nu, n, nvals, counts, filling, acc);
                filling.pop();
                counts[j] += 1;
            }
        }
    }

    /// Divide out the `q^{min inv}` floor — the cospin normalization of the
    /// tuple model.
    fn floored(f: &Monomial<Q>, floor: u32) -> Monomial<Q> {
        let divisor = q_pow(floor);
        let mut out = Monomial::zero();
        for (mu, p) in f.terms() {
            out.add_term(mu.clone(), p.divide_exact(&divisor).expect("floor divides"));
        }
        out
    }

    /// **The convention gate.** One attacking pair, and the line that dies if
    /// reading order, the attack rule, or the inv orientation drifts.
    #[test]
    fn the_convention_gate() {
        let g: Monomial<Q> = llt_g(&tuple(&[&[1], &[1]]));
        assert_eq!(g.coeff(&part(&[2])), poly(&[(0, 1)]), "m_2");
        assert_eq!(g.coeff(&part(&[1, 1])), poly(&[(0, 1), (1, 1)]), "m_11");
        assert_eq!(g.terms().len(), 2);
    }

    /// R1 against direct semistandard enumeration, on straight and skew tuples.
    /// This is the identity that makes \[HHL\] (82) an identity rather than a
    /// convention, so it is checked on everything else's behalf.
    #[test]
    fn standardization_agrees_with_direct_enumeration() {
        let straight: &[&[&[u32]]] = &[
            &[&[2], &[1]],
            &[&[1, 1], &[1]],
            &[&[2, 1], &[1]],
            &[&[1], &[1], &[1]],
            &[&[2], &[2]],
            &[&[2, 1], &[2]],
        ];
        for shapes in straight {
            let nu = tuple(shapes);
            assert_eq!(llt_g::<i64>(&nu), g_direct(&nu), "{shapes:?}");
        }
        let skews: &[&[(&[u32], &[u32])]] = &[
            &[(&[2, 1], &[1]), (&[1], &[])],
            &[(&[2, 2], &[1]), (&[2], &[1])],
        ];
        for sk in skews {
            let shapes: Vec<(Partition, Partition)> =
                sk.iter().map(|(o, i)| (part(o), part(i))).collect();
            let nu = SkewTuple::from_skews(&shapes, &vec![0i32; sk.len()]);
            assert_eq!(llt_g::<i64>(&nu), g_direct(&nu), "{sk:?}");
        }
    }

    /// Contents are the model's content: shifting one component's offset must
    /// move the answer, or the offsets are not being read at all.
    #[test]
    fn content_offsets_are_load_bearing() {
        let shapes = [part(&[1]), part(&[1])];
        let a: Monomial<Q> = llt_g(&SkewTuple::from_partitions(&shapes, &[0, 0]));
        let b: Monomial<Q> = llt_g(&SkewTuple::from_partitions(&shapes, &[0, 5]));
        assert_ne!(a, b);
    }

    /// `q = 1` collapses the tuple model to a product of skew Schur functions.
    /// Checked against [`crate::skew_lr`], which shares no code with this.
    #[test]
    fn at_q_one_it_is_a_product_of_skew_schurs() {
        let cases: &[&[(&[u32], &[u32])]] = &[
            &[(&[2, 1], &[]), (&[1], &[])],
            &[(&[2, 2], &[1]), (&[2, 1], &[])],
            &[(&[3, 1], &[1]), (&[2], &[]), (&[1, 1], &[])],
        ];
        for sk in cases {
            let shapes: Vec<(Partition, Partition)> =
                sk.iter().map(|(o, i)| (part(o), part(i))).collect();
            let nu = SkewTuple::from_skews(&shapes, &vec![0i32; sk.len()]);
            let got: Schur<i64> = {
                let g: Monomial<Q> = llt_g(&nu);
                let mut at_one: Monomial<i64> = Monomial::zero();
                for (mu, p) in g.terms() {
                    at_one.add_term(mu.clone(), p.eval(&1, &1));
                }
                at_one.to_schur()
            };
            let mut want: Schur<i64> = Schur::monomial(Partition::default(), 1);
            for (o, i) in &shapes {
                let mut factor: Schur<i64> = Schur::zero();
                for (nu, c) in crate::skew_lr::expand_skew_shared(o, i).iter() {
                    factor.add_term(nu.clone(), *c as i64);
                }
                want = want.mul(&factor);
            }
            assert_eq!(got, want, "{sk:?} at q = 1");
        }
    }

    /// **The quotient dictionary**, the law that ties the two models: `q^{−min
    /// inv} G_{k-quotient(λ)} = G̃^(k)_λ`, offsets zero. Swept over every
    /// empty-core shape of the small degrees, and the sweep must *contain*
    /// shapes with a nonzero floor or it is not testing the floor.
    #[test]
    fn the_quotient_dictionary_is_the_floored_tuple() {
        let mut floors = 0;
        let mut shapes = 0;
        for k in 2..=3u32 {
            for n in [k, 2 * k, 3 * k] {
                for lambda in crate::partitions_of(n) {
                    let Some(nu) = SkewTuple::quotient(&lambda, k) else {
                        continue;
                    };
                    if nu.size() == 0 {
                        continue;
                    }
                    shapes += 1;
                    let floor = llt_min_inv(&nu);
                    if floor > 0 {
                        floors += 1;
                    }
                    let got = floored(&llt_g::<i64>(&nu), floor);
                    let want: Monomial<Q> = llt_gtilde(&lambda, k);
                    assert_eq!(got, want, "λ = {lambda}, k = {k}, floor {floor}");
                }
            }
        }
        assert!(shapes >= 40, "the sweep shrank to {shapes} shapes");
        assert!(
            floors > 0,
            "no shape in the sweep has a nonzero min-inv floor, so the floor is \
             untested — λ = (2,2,2) at k = 2 is the smallest witness"
        );
    }

    /// The floor is *forced*, and the smallest witness deserves its own line:
    /// λ = (2,2,2) at k = 2 has 2-quotient `((1), (1,1))` and `min inv = 1`.
    ///
    /// ⚠️ λ = (2,2) is the near miss, and was written down as the witness
    /// twice before being checked: its 2-quotient is `((1), (1))`, whose `G`
    /// is `m₂ + (1+q)m₁₁` — floor 0, which is what the convention gate above
    /// depends on, so the two claims were never consistent. Both directions
    /// are asserted below. The floored shapes at k = 2 through |λ| = 8 are
    /// (2,2,2), (4,2,2), (3,3,2), (2,2,2,2) and (2,2,2,1,1); the law itself
    /// is unaffected.
    #[test]
    fn the_min_inv_floor_is_forced() {
        let nu = SkewTuple::quotient(&part(&[2, 2, 2]), 2).expect("empty 2-core");
        assert_eq!(
            part(&[2, 2, 2]).k_quotient(2),
            vec![part(&[1]), part(&[1, 1])]
        );
        assert_eq!(llt_min_inv(&nu), 1);
        // And the spec's stated witness really does have no floor.
        let flat = SkewTuple::quotient(&part(&[2, 2]), 2).expect("empty 2-core");
        assert_eq!(llt_min_inv(&flat), 0);
    }

    /// Tuples lose absolute spin: λ = (1,1,1,1) at k = 2 has `s* = 1` while its
    /// quotient's inv range is a single point. No tuple statistic can recover
    /// the spin grading, which is why `llt_h` takes a partition and a level.
    #[test]
    fn tuples_lose_absolute_spin() {
        let lambda = part(&[1, 1, 1, 1]);
        let raw: Monomial<Q> = llt_g_lt(&lambda, 2);
        let smax = raw
            .terms()
            .values()
            .filter_map(|p| p.degrees())
            .map(|(a, _)| a)
            .max()
            .unwrap();
        assert_eq!(smax, 2, "spin_LT of the one tableau, i.e. s* = 1");
        let nu = SkewTuple::quotient(&lambda, 2).unwrap();
        assert_eq!(llt_min_inv(&nu), llt_max_inv(&nu));
    }

    /// **\[LLT\] Ex 6.8(i)**: `G̃_{(3,3,3,2,1)}` at k = 3, in m and in s.
    #[test]
    fn llt_example_68i() {
        let g: Monomial<Q> = llt_gtilde(&part(&[3, 3, 3, 2, 1]), 3);
        assert_eq!(g.coeff(&part(&[3, 1])), poly(&[(0, 1)]));
        assert_eq!(g.coeff(&part(&[2, 2])), poly(&[(0, 1), (1, 1)]));
        assert_eq!(g.coeff(&part(&[2, 1, 1])), poly(&[(0, 2), (1, 2), (2, 1)]));
        assert_eq!(
            g.coeff(&part(&[1, 1, 1, 1])),
            poly(&[(0, 3), (1, 5), (2, 3), (3, 1)])
        );
        assert_eq!(g.terms().len(), 4);

        let s = g.to_schur();
        assert_eq!(s.coeff(&part(&[3, 1])), poly(&[(0, 1)]));
        assert_eq!(s.coeff(&part(&[2, 2])), poly(&[(1, 1)]));
        assert_eq!(s.coeff(&part(&[2, 1, 1])), poly(&[(1, 1), (2, 1)]));
        assert_eq!(s.coeff(&part(&[1, 1, 1, 1])), poly(&[(3, 1)]));
        assert_eq!(s.terms().len(), 4);
    }

    /// **\[LLT\] Ex 6.8(ii)**: `H^(2)_{3211}` in the Schur basis.
    #[test]
    fn llt_example_68ii() {
        let h: Monomial<Q> = llt_h(&part(&[3, 2, 1, 1]), 2);
        let s = h.to_schur();
        let want: &[(&[u32], &[(u32, i64)])] = &[
            (&[3, 2, 1, 1], &[(0, 1)]),
            (&[3, 2, 2], &[(1, 1)]),
            (&[3, 3, 1], &[(1, 1)]),
            (&[4, 1, 1, 1], &[(1, 1)]),
            (&[4, 2, 1], &[(1, 1), (2, 1)]),
            (&[4, 3], &[(2, 1)]),
            (&[5, 1, 1], &[(2, 1)]),
            (&[5, 2], &[(3, 1)]),
        ];
        for (nu, terms) in want {
            assert_eq!(s.coeff(&part(nu)), poly(terms), "s_{nu:?}");
        }
        assert_eq!(s.terms().len(), want.len());
    }

    /// **\[LT\] Ex 4.1**: the same shape in the spin-generating (43) grading.
    #[test]
    fn lt_example_41() {
        let g: Monomial<Q> = llt_g_lt(&part(&[3, 3, 3, 2, 1]), 3);
        let s = g.to_schur();
        assert_eq!(s.coeff(&part(&[3, 1])), poly(&[(7, 1)]));
        assert_eq!(s.coeff(&part(&[2, 2])), poly(&[(5, 1)]));
        assert_eq!(s.coeff(&part(&[2, 1, 1])), poly(&[(3, 1), (5, 1)]));
        assert_eq!(s.coeff(&part(&[1, 1, 1, 1])), poly(&[(1, 1)]));
    }

    /// `k = 1`: one-ribbons are cells, spin is zero, `G̃^(1)_λ = s_λ`.
    #[test]
    fn level_one_is_the_schur_function() {
        for n in 0..=6u32 {
            for lambda in crate::partitions_of(n) {
                let got: Schur<Q> = llt_schur(&lambda, 1);
                let want: Schur<Q> = Schur::monomial(lambda.clone(), <Q as Ring>::one());
                assert_eq!(got, want, "G̃^(1)_{lambda}");
            }
        }
    }

    /// **\[LLT\] (28)**, internally: `H = q^{s*} G̃(1/q)` on the shapes where
    /// the spin grading is integral.
    #[test]
    fn h_is_the_cospin_flip() {
        for k in 2..=3u32 {
            for n in 0..=4u32 {
                for mu in crate::partitions_of(n) {
                    let cos: Monomial<Q> = llt_h_tilde(&mu, k);
                    let h: Monomial<Q> = llt_h(&mu, k);
                    let sstar = cos
                        .terms()
                        .values()
                        .filter_map(|p| p.degrees())
                        .map(|(a, _)| a)
                        .max()
                        .unwrap_or(0);
                    let mut flipped: Monomial<Q> = Monomial::zero();
                    for (nu, p) in cos.terms() {
                        let mut np = QtPoly::zero();
                        for (&(a, b), c) in p.terms() {
                            np.add_term(sstar - a, b, *c);
                        }
                        flipped.add_term(nu.clone(), np);
                    }
                    assert_eq!(flipped, h, "H^({k})_{mu} = q^s* H̃(1/q)");
                }
            }
        }
    }

    /// **\[LLT\] Thm 6.6**: at large level the spin ribbon function *is*
    /// Hall–Littlewood `Q'_μ`. Checked against [`crate::hl`] in-crate — two
    /// entirely different recursions landing on the same polynomial.
    #[test]
    fn large_level_is_hall_littlewood() {
        for n in 0..=4u32 {
            for mu in crate::partitions_of(n) {
                let k = mu.len().max(1) as u32;
                let got: Schur<Q> = llt_h::<i64>(&mu, k).to_schur();
                // `hl` grades in t; this module grades in q. Same polynomial,
                // different slot.
                let want_t: Schur<Q> = crate::hl::hall_littlewood(&mu);
                let mut want: Schur<Q> = Schur::zero();
                for (nu, p) in want_t.terms() {
                    let mut np = QtPoly::zero();
                    for (&(a, b), c) in p.terms() {
                        assert_eq!(a, 0, "Hall–Littlewood lives in t alone");
                        np.add_term(b, 0, *c);
                    }
                    want.add_term(nu.clone(), np);
                }
                assert_eq!(got, want, "H^({k})_{mu} vs Q'_{mu}");
            }
        }
    }

    /// **ω-duality**: `G_{ν'}(x; q) = q^{A(ν)} ω G_ν(x; 1/q)`.
    #[test]
    fn omega_duality() {
        for shapes in [
            &[&[1u32][..], &[1]][..],
            &[&[2], &[1]],
            &[&[2, 1], &[1]],
            &[&[2], &[2]],
        ] {
            let nu = tuple(shapes);
            let a = nu.attacking_pairs() as u32;
            let lhs: Schur<Q> = llt_g::<i64>(&nu.conjugate()).to_schur();
            // ω on the Schur basis transposes the index; the q-flip is
            // q^A · p(1/q), which is exactly reversing each coefficient.
            let g: Schur<Q> = llt_g::<i64>(&nu).to_schur();
            let mut rhs: Schur<Q> = Schur::zero();
            for (lambda, p) in g.terms() {
                let mut np = QtPoly::zero();
                for (&(e, b), c) in p.terms() {
                    assert!(e <= a, "the flip must stay polynomial");
                    np.add_term(a - e, b, *c);
                }
                rhs.add_term(lambda.conjugate(), np);
            }
            assert_eq!(lhs, rhs, "{shapes:?}");
        }
    }

    /// The three vertical-strip presentations must agree **pointwise per
    /// path**, not merely in total: the reverse-run tuple, the dinv model, and
    /// the coloring model on the dinv-faithful graph.
    #[test]
    fn the_three_vertical_strip_models_agree() {
        for n in 0..=6usize {
            for_each_area(n, &mut |area| {
                let from_tuple: Monomial<Q> = llt_g(&SkewTuple::from_area(area));
                let from_graph: Monomial<Q> = llt_graph(&DecoratedGraph::from_area(area));
                assert_eq!(from_tuple, from_graph, "area {area:?}");
                assert_eq!(from_tuple, dinv_model(area), "area {area:?}");
            });
        }
    }

    /// `Σ_labellings q^{dinv} x^ℓ` straight from \[HRW\], sharing nothing with
    /// either model under test.
    fn dinv_model(area: &[u32]) -> Monomial<Q> {
        let n = area.len();
        let mut out = Monomial::zero();
        for mu in crate::partitions_of(n as u32) {
            let nvals = mu.len();
            let mut counts: Vec<u32> = mu.parts().to_vec();
            let mut labels: Vec<u32> = Vec::with_capacity(n);
            let mut acc = QtPoly::zero();
            label_rec(area, nvals, &mut counts, &mut labels, &mut acc);
            out.add_term(mu, acc);
        }
        out
    }

    fn label_rec(
        area: &[u32],
        nvals: usize,
        counts: &mut [u32],
        labels: &mut Vec<u32>,
        acc: &mut Q,
    ) {
        let i = labels.len();
        if i == area.len() {
            let mut dinv = 0u32;
            for x in 0..area.len() {
                for y in (x + 1)..area.len() {
                    if (area[x] == area[y] && labels[x] < labels[y])
                        || (area[x] == area[y] + 1 && labels[x] > labels[y])
                    {
                        dinv += 1;
                    }
                }
            }
            acc.add_term(dinv, 0, 1);
            return;
        }
        let rise = i > 0 && area[i] == area[i - 1] + 1;
        for j in 0..nvals {
            if counts[j] == 0 {
                continue;
            }
            let v = j as u32 + 1;
            if rise && v <= labels[i - 1] {
                continue;
            }
            counts[j] -= 1;
            labels.push(v);
            label_rec(area, nvals, counts, labels, acc);
            labels.pop();
            counts[j] += 1;
        }
    }

    /// **The shuffle refinement**: `Σ_D t^{area} G_D = ∇e_n`, against
    /// [`crate::deltaop`]'s operator side.
    #[test]
    fn the_shuffle_refinement_is_nabla_e() {
        for n in 0..=5u32 {
            let mut total: Monomial<QtPoly<Rational>> = Monomial::zero();
            for (_, g) in nabla_e_by_path::<Rational>(n) {
                total = total.add(&g);
            }
            let want = crate::deltaop::nabla_e::<Rational>(n);
            assert_eq!(total.to_schur(), want, "∇e_{n} by path");
        }
    }

    /// Per-path Schur positivity is **recorded, never repaired**: it is a
    /// theorem for tuples of partitions (\[LT\]/Varagnolo–Vasserot) and rests
    /// on an unpublished preprint (\[GH\]) in general, so a negative
    /// coefficient here is a finding of the first order, not a bug to paper
    /// over.
    #[test]
    fn every_path_piece_is_schur_positive() {
        for n in 0..=5u32 {
            for (area, g) in nabla_e_by_path::<i64>(n) {
                for (lambda, p) in g.to_schur().terms() {
                    for (_, c) in p.terms() {
                        assert!(
                            *c >= 0,
                            "negative Schur coefficient {c} on s_{lambda} for area \
                             {area:?}: report this, do not repair it"
                        );
                    }
                }
            }
        }
    }

    /// **The chromatic bridge.** Shareshian–Wachs define `X_Γ(x;q)` by summing
    /// `q^{asc}` over *proper* colorings; the LLT side sums over *all* of them
    /// and undoes the difference with a `(q−1)`-plethysm. Both are computed
    /// here and they must meet.
    #[test]
    fn the_chromatic_bridge_is_shareshian_wachs() {
        for n in 0..=5usize {
            for_each_area(n, &mut |area| {
                let g = DecoratedGraph::unit_interval(area);
                let got: Monomial<QtPoly<Rational>> = chromatic_from_llt(&g);
                assert_eq!(got, proper_colorings(&g), "area {area:?}");
            });
        }
    }

    /// The graph must keep its isolated vertices. A two-vertex edgeless graph
    /// is not the empty graph, and collapsing it is a silent way to a plausible
    /// wrong answer.
    #[test]
    fn isolated_vertices_survive() {
        let g = DecoratedGraph::new(2, &[], &[]);
        let x: Monomial<QtPoly<Rational>> = chromatic_from_llt(&g);
        assert_eq!(x.coeff(&part(&[2])), <QtPoly<Rational> as Ring>::one());
        assert_eq!(
            x.coeff(&part(&[1, 1]))
                .eval(&Rational::from_int(1), &Rational::from_int(1)),
            Rational::from_int(2)
        );
    }

    /// `Σ_{proper κ} q^{asc(κ)} x^κ`, straight from Shareshian–Wachs.
    fn proper_colorings(g: &DecoratedGraph) -> Monomial<QtPoly<Rational>> {
        let n = g.n as usize;
        let edges: Vec<(u32, u32)> = g.weak.iter().map(|&(u, v)| (u.min(v), u.max(v))).collect();
        let mut out = Monomial::zero();
        for mu in crate::partitions_of(g.n) {
            let nvals = mu.len();
            let mut counts: Vec<u32> = mu.parts().to_vec();
            let mut kappa: Vec<u32> = Vec::with_capacity(n);
            let mut acc: QtPoly<Rational> = QtPoly::zero();
            proper_rec(n, nvals, &edges, &mut counts, &mut kappa, &mut acc);
            out.add_term(mu, acc);
        }
        out
    }

    fn proper_rec(
        n: usize,
        nvals: usize,
        edges: &[(u32, u32)],
        counts: &mut [u32],
        kappa: &mut Vec<u32>,
        acc: &mut QtPoly<Rational>,
    ) {
        let v = kappa.len();
        if v == n {
            let asc = edges
                .iter()
                .filter(|&&(a, b)| kappa[a as usize] < kappa[b as usize])
                .count() as u32;
            acc.add_term(asc, 0, Rational::from_int(1));
            return;
        }
        for j in 0..nvals {
            if counts[j] == 0 {
                continue;
            }
            let c = j as u32 + 1;
            if edges
                .iter()
                .any(|&(a, b)| b as usize == v && kappa[a as usize] == c)
            {
                continue;
            }
            counts[j] -= 1;
            kappa.push(c);
            proper_rec(n, nvals, edges, counts, kappa, acc);
            kappa.pop();
            counts[j] += 1;
        }
    }

    /// `p(q + d, t)`, by the binomial theorem. The \[AS\] formula lands on
    /// `G(x; q+1)`, so comparing it with `G` needs one of these in each
    /// direction.
    fn shift_q<C: Ring>(p: &QtPoly<C>, d: i64) -> QtPoly<C> {
        let mut out = QtPoly::zero();
        for (&(a, b), c) in p.terms() {
            let mut binom: i128 = 1;
            for i in (0..=a).rev() {
                // binom = C(a, i)
                let scale = binom * (d as i128).pow(a - i);
                out.add_term(i, b, C::from_i128(scale).mul(c));
                if i > 0 {
                    binom = binom * i as i128 / (a - i + 1) as i128;
                }
            }
        }
        out
    }

    /// **The \[AS\] orientation formula**, on the unicellular corpus: the
    /// e-expansion of `G(x; q+1)` must reproduce `G` after `q ↦ q−1`. And by
    /// \[DA\]'s theorem the coefficients are non-negative — a theorem, so a
    /// violation here is a bug in this crate, the *reverse* of the posture in
    /// `every_path_piece_is_schur_positive`.
    #[test]
    fn the_as_orientation_formula_reproduces_g() {
        for n in 0..=5usize {
            for_each_area(n, &mut |area| {
                let g = DecoratedGraph::unit_interval(area);
                let e = llt_e_expansion::<i64>(&g);
                for (lambda, p) in &e {
                    for (_, c) in p.terms() {
                        assert!(*c >= 0, "[DA] positivity fails on e_{lambda}: our bug");
                    }
                }
                let mut got: Elementary<Q> = Elementary::zero();
                for (lambda, p) in e {
                    got.add_term(lambda, shift_q(&p, -1));
                }
                let want: Elementary<Q> = Elementary::from_schur(&llt_graph::<i64>(&g).to_schur());
                assert_eq!(got, want, "area {area:?}");
            });
        }
    }

    /// **\[AS\] Ex 6.1** and **\[DA\] Ex 5.6**, the two published Schröder-path
    /// fixtures. Words are a test constructor, not API: graphs are the API.
    #[test]
    fn published_schroder_fixtures() {
        let g = graph_from_word("nndee");
        let s = llt_graph::<i64>(&g).to_schur();
        assert_eq!(s.coeff(&part(&[1, 1, 1])), poly(&[(2, 1)]), "[AS] Ex 6.1");
        assert_eq!(s.coeff(&part(&[2, 1])), poly(&[(1, 1)]), "[AS] Ex 6.1");
        assert_eq!(s.terms().len(), 2);

        let g = graph_from_word("ndndee");
        let e: Elementary<Q> = Elementary::from_schur(&llt_graph::<i64>(&g).to_schur());
        // q e_1 e_3 + q(q−1) e_4
        assert_eq!(e.coeff(&part(&[3, 1])), poly(&[(1, 1)]), "[DA] Ex 5.6");
        assert_eq!(
            e.coeff(&part(&[4])),
            poly(&[(1, -1), (2, 1)]),
            "[DA] Ex 5.6"
        );
        assert_eq!(e.terms().len(), 2);
    }

    /// A Schröder path word in `{n, e, d}` to its decorated graph, \[AS\] §2.6:
    /// vertices along the diagonal, an edge when the cell is below the path,
    /// and a diagonal step makes its endpoint pair strict.
    fn graph_from_word(word: &str) -> DecoratedGraph {
        let (mut x, mut y) = (0i64, 0i64);
        let mut colheight: Vec<i64> = Vec::new();
        let mut diag: Vec<(usize, usize)> = Vec::new();
        for st in word.chars() {
            match st {
                'n' => y += 1,
                'e' => {
                    colheight.push(y);
                    x += 1;
                }
                'd' => {
                    colheight.push(y + 1);
                    x += 1;
                    y += 1;
                    diag.push((x as usize, y as usize));
                }
                c => panic!("unexpected step {c}"),
            }
        }
        let n = x as usize;
        let mut weak = Vec::new();
        for u in 1..=n {
            for v in (u + 1)..=n {
                if (v as i64) <= colheight[u - 1] {
                    weak.push(((u - 1) as u32, (v - 1) as u32));
                }
            }
        }
        let strict: Vec<(u32, u32)> = diag
            .iter()
            .map(|&(u, v)| ((u - 1) as u32, (v - 1) as u32))
            .collect();
        // A diagonal step's endpoint pair is also below the path, so the cell
        // rule already emitted it; the strict decoration replaces it rather than
        // joining it.
        weak.retain(|e| !strict.contains(e));
        DecoratedGraph::new(n as u32, &weak, &strict)
    }

    /// **The \[HHL\] assembly**: the fourth route to `H̃_μ`, against the three
    /// in [`crate::qtkostka`].
    #[test]
    fn the_hhl_assembly_is_htilde() {
        for n in 0..=4u32 {
            for mu in crate::partitions_of(n) {
                let got = htilde_by_llt::<i64>(&mu).to_schur();
                let want = crate::qtkostka::macdonald_ht::<i64>(&mu);
                assert_eq!(got, want, "H̃_{mu} by LLT decomposition");
            }
        }
    }

    /// **R3**: the Fock column against the ribbon side, at `q = −v`.
    ///
    /// The range is chosen to *pin the variable dictionary*, not for coverage.
    /// A single monomial cannot tell `q = −v` from its mirrors, so the sweep
    /// has to reach a multi-term Kazhdan–Lusztig entry, and the final assertion
    /// fails if it does not.
    ///
    /// ⚠️ The prototype swept k = 2 through degree 4 and believed that pinned
    /// the dictionary. It does not: **every** Schur coefficient of `G_LT` at
    /// k = 2 is a single monomial through |μ| = 10. The first multi-term
    /// entries are at **k = 3** — `s_{21}` in `G_LT,(3,3,3)` is `q² + q⁴`, and
    /// \[LT\] Ex 4.1's `q³ + q⁵` on `s_{211}` is the k = 3 shape (3,3,3,2,1).
    /// So the k = 2-only sweep did not in fact pin `q = −v`; this test does,
    /// and the `multi > 0` assertion below is what keeps it honest.
    #[test]
    fn fock_straightening_is_the_ribbon_column_at_q_eq_minus_v() {
        let mut multi = 0;
        for (k, top) in [(2u32, 4u32), (3, 3)] {
            for lambda in (1..=top).flat_map(crate::partitions_of) {
                let got: HashMap<Partition, Q> =
                    llt_kl_column::<i64>(&lambda, k).into_iter().collect();
                for mu in crate::partitions_of(k * lambda.size()) {
                    let s: Schur<Q> = llt_g_lt::<i64>(&mu, k).to_schur();
                    let c = s.coeff(&lambda);
                    if c.terms().count() > 1 {
                        multi += 1;
                    }
                    // q ↦ −v: negate the odd-degree coefficients.
                    let mut want = QtPoly::zero();
                    for (&(a, b), cc) in c.terms() {
                        let v = if a % 2 == 1 { cc.neg() } else { *cc };
                        want.add_term(a, b, v);
                    }
                    let have = got.get(&mu).cloned().unwrap_or_else(QtPoly::zero);
                    assert_eq!(have, want, "c^{lambda}_{mu} at k={k}, q = −v");
                }
            }
        }
        assert!(
            multi > 0,
            "no multi-term KL entry in range: the q = −v dictionary is untested"
        );
    }

    /// The block-prefix strip enumeration against the naive one: all
    /// `C(rows, m)` subsets, filtered. Same lemma, shared nothing.
    #[test]
    fn strip_blocks_agree_with_naive_subsets() {
        for k in 2..=4u32 {
            for n in 0..=8u32 {
                for lambda in crate::partitions_of(n) {
                    // `abacus_of` works on the *conjugate*, so the bead count is
                    // bounded below by λ₁, not by ℓ(λ).
                    let rows = lambda.part(0).max(1) as usize + 4;
                    if lambda.len() + rows + 4 >= 128 {
                        continue;
                    }
                    let beta = abacus_of(&lambda, rows);
                    for m in 0..=4u32 {
                        let mut got: Vec<(Abacus, u32)> = Vec::new();
                        let mut sc = StripScratch::default();
                        for_each_strip_up(beta, k, m, &mut sc, &mut |nb, spin| {
                            got.push((nb, spin))
                        });
                        got.sort_unstable();
                        let mut want = naive_strips(beta, rows, k, m);
                        want.sort_unstable();
                        assert_eq!(got, want, "λ={lambda} k={k} m={m} rows={rows}");
                    }
                }
            }
        }
    }

    /// The bit-trick containment test against the definition it replaced:
    /// build both shapes and ask [`Partition::contains`].
    #[test]
    fn abacus_containment_is_partition_containment() {
        for n in 0..=8u32 {
            for lambda in crate::partitions_of(n) {
                let rows = lambda.part(0).max(1) as usize + 3;
                let target = abacus_of(&lambda, rows);
                for m in 0..=n {
                    for nu in crate::partitions_of(m) {
                        if nu.len() > rows || nu.part(0) as usize > rows {
                            continue;
                        }
                        let beta = abacus_of(&nu, rows);
                        assert_eq!(
                            abacus_contained(beta, target),
                            lambda.contains(&nu),
                            "{nu} ⊆ {lambda}?"
                        );
                    }
                }
            }
        }
    }

    /// The all-weights-at-once walk the ribbon engine actually runs, against
    /// the fixed-weight one the test above pins to the naive subsets.
    ///
    /// Without this the chain is broken:
    /// `strip_blocks_agree_with_naive_subsets` checks `for_each_strip_up`,
    /// which after the R2 profiling work is used only by that test — production
    /// goes through `for_each_strip_any`.
    #[test]
    fn the_all_weights_walk_agrees_with_the_fixed_weight_one() {
        for k in 1..=4u32 {
            for n in 0..=8u32 {
                for lambda in crate::partitions_of(n) {
                    let rows = lambda.part(0).max(1) as usize + 4;
                    if lambda.len() + rows + 4 >= 128 {
                        continue;
                    }
                    let beta = abacus_of(&lambda, rows);
                    let max_w = 5u32;
                    let mut sc = StripScratch::default();
                    let mut got: Vec<(u32, Abacus, u32)> = Vec::new();
                    for_each_strip_any(beta, k, max_w, &mut sc, &mut |w, nb, spin| {
                        got.push((w, nb, spin))
                    });
                    got.sort_unstable();
                    let mut want: Vec<(u32, Abacus, u32)> = Vec::new();
                    for m in 1..=max_w {
                        for_each_strip_up(beta, k, m, &mut sc, &mut |nb, spin| {
                            want.push((m, nb, spin))
                        });
                    }
                    want.sort_unstable();
                    assert_eq!(got, want, "λ={lambda} k={k} rows={rows}");
                }
            }
        }
    }

    /// Every m-subset of the beads, moved `+k`, keeping the results distinct.
    fn naive_strips(beta: Abacus, rows: usize, k: u32, m: u32) -> Vec<(Abacus, u32)> {
        let beads: Vec<u32> = (0..128u32).filter(|&b| beta >> b & 1 == 1).collect();
        assert_eq!(beads.len(), rows);
        let mut out = Vec::new();
        for mask in 0u64..(1u64 << rows) {
            if (mask.count_ones()) != m {
                continue;
            }
            let mut moved: Abacus = 0;
            for (i, &b) in beads.iter().enumerate() {
                if mask >> i & 1 == 1 {
                    moved |= 1u128 << b;
                }
            }
            let unmoved = beta & !moved;
            let shifted = moved << k;
            if shifted & unmoved != 0 {
                continue;
            }
            let mut cross = 0u32;
            for i in 1..k {
                cross += ((moved << i) & unmoved).count_ones();
            }
            out.push((unmoved | shifted, (k - 1) * m - cross));
        }
        out
    }

    /// The pruned single-shape walk and the unpruned whole-degree walk are
    /// different traversals of the same tree and must land on the same numbers.
    #[test]
    fn the_two_walks_agree() {
        for k in 1..=3u32 {
            for n in 0..=4u32 {
                let table: Vec<(Partition, Monomial<Q>)> = llt_gtilde_table(n, k);
                for (lambda, got) in &table {
                    assert_eq!(*got, llt_gtilde::<i64>(lambda, k), "λ={lambda} k={k}");
                }
                // and nothing with an empty k-core is missing
                for lambda in crate::partitions_of(k * n) {
                    if lambda.has_empty_k_core(k) {
                        assert!(
                            table.iter().any(|(l, _)| *l == lambda),
                            "{lambda} missing from the k={k} degree-{n} table"
                        );
                    }
                }
            }
        }
    }

    /// The `H` table must agree term for term with the one-at-a-time calls.
    #[test]
    fn h_table_agrees_with_individual_calls() {
        for k in 1..=3u32 {
            for n in 0..=4u32 {
                for (mu, got) in llt_h_table::<i64>(n, k) {
                    assert_eq!(got, llt_h::<i64>(&mu, k), "H^({k})_{mu}");
                }
            }
        }
    }

    /// **Fixed-width honesty**: the same computation at two coefficient widths.
    /// A silent `i64` wrap would show up as a disagreement with `i128`.
    #[test]
    fn fixed_width_ladders_agree() {
        for k in 1..=3u32 {
            for n in 0..=4u32 {
                for mu in crate::partitions_of(n) {
                    let narrow: Monomial<QtPoly<i64>> = llt_h(&mu, k);
                    let wide: Monomial<QtPoly<i128>> = llt_h(&mu, k);
                    assert_eq!(narrow.terms().len(), wide.terms().len(), "{mu} at k={k}");
                    for (nu, p) in narrow.terms() {
                        let w = wide.coeff(nu);
                        for (&(a, b), c) in p.terms() {
                            assert_eq!(i128::from(*c), w.coeff(a, b), "{mu}/{nu} at k={k}");
                        }
                    }
                }
            }
        }
    }

    /// The fundamental expansion is a refinement of the monomial one: summing
    /// the compositions that refine μ must give `[m_μ] G_ν`.
    #[test]
    fn the_fundamental_expansion_refines_the_monomial_one() {
        for shapes in [
            &[&[2u32][..], &[1]][..],
            &[&[1], &[1], &[1]],
            &[&[2, 1], &[2]],
        ] {
            let nu = tuple(shapes);
            let f = llt_fundamental::<i64>(&nu);
            let m: Monomial<Q> = llt_g(&nu);
            let total: u32 = f.iter().map(|(c, _)| c.iter().sum::<u32>()).max().unwrap();
            assert_eq!(total, nu.size() as u32, "compositions of n");
            for mu in crate::partitions_of(nu.size() as u32) {
                let allowed = partial_sum_mask(&mu);
                let mut acc = QtPoly::zero();
                for (comp, p) in &f {
                    let mut mask = 0u64;
                    let mut acc_s = 0u32;
                    for &c in &comp[..comp.len() - 1] {
                        acc_s += c;
                        mask |= 1u64 << (acc_s - 1);
                    }
                    if mask & !allowed == 0 {
                        acc.add_assign(p);
                    }
                }
                assert_eq!(acc, m.coeff(&mu), "[m_{mu}]");
            }
        }
    }
}
