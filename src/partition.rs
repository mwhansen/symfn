//! Integer partitions — a "real type" with its invariant enforced at construction.
//!
//! Symmetrica represented partitions as raw vectors you simply had to trust.
//! Here a `Partition` is *always* weakly decreasing with positive parts, because
//! the only ways to build one either normalize or validate. Downstream code can
//! rely on that invariant instead of re-checking it.

use core::fmt;

/// A partition λ = (λ₁ ≥ λ₂ ≥ … ≥ λ_k > 0), stored as its nonzero parts.
///
/// The empty partition (the unique partition of 0) is the empty vector.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub struct Partition(Vec<u32>);

/// Why a slice failed to be a valid partition in strict construction.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum PartitionError {
    NotWeaklyDecreasing,
    ZeroPart,
}

impl Partition {
    /// Build from parts, **normalizing**: drops zeros and sorts weakly
    /// decreasing. Convenient for callers that don't want to pre-sort.
    pub fn new<I: IntoIterator<Item = u32>>(parts: I) -> Self {
        let mut v: Vec<u32> = parts.into_iter().filter(|&x| x != 0).collect();
        v.sort_unstable_by(|a, b| b.cmp(a));
        Partition(v)
    }

    /// Build from parts that are *claimed* already valid, validating the
    /// invariant and returning an error rather than silently fixing it.
    pub fn try_new<I: IntoIterator<Item = u32>>(parts: I) -> Result<Self, PartitionError> {
        let v: Vec<u32> = parts.into_iter().collect();
        for w in v.windows(2) {
            if w[0] < w[1] {
                return Err(PartitionError::NotWeaklyDecreasing);
            }
        }
        if v.iter().any(|&x| x == 0) {
            return Err(PartitionError::ZeroPart);
        }
        Ok(Partition(v))
    }

    /// Internal: build from a slice already known to be weakly decreasing and
    /// positive (e.g. the partition generator). Debug-asserts the invariant.
    pub(crate) fn from_sorted(v: Vec<u32>) -> Self {
        debug_assert!(v.windows(2).all(|w| w[0] >= w[1]) && v.iter().all(|&x| x != 0));
        Partition(v)
    }

    /// The parts λ₁ … λ_k.
    #[inline]
    pub fn parts(&self) -> &[u32] {
        &self.0
    }

    /// Number of (nonzero) parts, the *length* ℓ(λ).
    #[inline]
    pub fn len(&self) -> usize {
        self.0.len()
    }

    #[inline]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// |λ| = Σ λ_i, the integer being partitioned.
    #[inline]
    pub fn size(&self) -> u32 {
        self.0.iter().sum()
    }

    /// λ_i with the convention λ_i = 0 for i ≥ ℓ(λ).
    #[inline]
    pub fn part(&self, i: usize) -> u32 {
        self.0.get(i).copied().unwrap_or(0)
    }

    /// Young-diagram containment: does `self` contain `other` (μ ⊆ λ)?
    /// True iff other_i ≤ self_i for all i.
    pub fn contains(&self, other: &Partition) -> bool {
        (0..other.len()).all(|i| self.part(i) >= other.part(i))
    }

    /// The conjugate (transpose) partition λ', where λ'_j = #{ i : λ_i ≥ j }.
    pub fn conjugate(&self) -> Partition {
        if self.is_empty() {
            return Partition::default();
        }
        let width = self.part(0) as usize; // largest part = number of columns
        let mut conj = vec![0u32; width];
        for &p in &self.0 {
            for col in 0..p as usize {
                conj[col] += 1;
            }
        }
        // Column counts are automatically weakly decreasing and positive.
        Partition::from_sorted(conj)
    }

    /// β-numbers with `rows` beads: `β_j = λ_j + rows − 1 − j`, for
    /// `j = 0 … rows−1`.
    ///
    /// The abacus encoding. A partition with at most `rows` parts is the same
    /// data as the strictly decreasing sequence `β`, and the two directions are
    /// [`beta_numbers`](Self::beta_numbers) / [`from_beta_numbers`].
    ///
    /// Padding matters and is harmless: raising `rows` by one shifts every β by
    /// one and adds a bead at position 0, so *every* `rows ≥ ℓ(λ)` encodes the
    /// same partition. Callers that move beads (ribbon strips, k-quotients)
    /// choose `rows` for the range they need, not for λ.
    pub fn beta_numbers(&self, rows: usize) -> Vec<u32> {
        assert!(
            rows >= self.len(),
            "an abacus needs at least ℓ(λ) beads to hold λ"
        );
        (0..rows)
            .map(|j| self.part(j) + (rows - 1 - j) as u32)
            .collect()
    }

    /// The partition encoded by a **strictly decreasing** β-number sequence.
    ///
    /// Inverse to [`beta_numbers`](Self::beta_numbers) with `rows = beta.len()`.
    pub fn from_beta_numbers(beta: &[u32]) -> Partition {
        debug_assert!(
            beta.windows(2).all(|w| w[0] > w[1]),
            "β-numbers must be strictly decreasing"
        );
        let l = beta.len();
        let parts: Vec<u32> = (0..l)
            .map(|j| beta[j] - (l - 1 - j) as u32)
            .filter(|&x| x > 0)
            .collect();
        Partition::from_sorted(parts)
    }

    /// The number of beads this crate uses for k-abacus work: the smallest
    /// multiple of `k` that strictly exceeds `ℓ(λ)`.
    ///
    /// A multiple of `k` so that each runner is filled to the bottom, and the
    /// residue of a bead — hence *which* quotient component it lands in — does
    /// not depend on the padding. Any larger multiple of `k` gives the same
    /// answer ([`k_quotient_is_independent_of_the_padding`] pins that).
    ///
    /// [`k_quotient_is_independent_of_the_padding`]: self
    fn k_rows(&self, k: u32) -> usize {
        let k = k as usize;
        let l = self.len().max(1);
        k * ((l + k) / k)
    }

    /// The **k-core** of λ: what is left after peeling k-rim-hooks as long as
    /// any can be peeled.
    ///
    /// On the abacus this is one move: slide every bead as far down its own
    /// runner as it will go. Runner residues are invariant under `β ↦ β ± k`,
    /// so a runner holding `c` beads ends with them at `r, r+k, …, r+(c−1)k`.
    pub fn k_core(&self, k: u32) -> Partition {
        assert!(k >= 1, "a k-core needs k ≥ 1");
        let rows = self.k_rows(k);
        let beta = self.beta_numbers(rows);
        let mut counts = vec![0usize; k as usize];
        for &b in &beta {
            counts[(b % k) as usize] += 1;
        }
        let mut packed: Vec<u32> = Vec::with_capacity(rows);
        for (r, &c) in counts.iter().enumerate() {
            for p in 0..c {
                packed.push(k * p as u32 + r as u32);
            }
        }
        packed.sort_unstable_by(|a, b| b.cmp(a));
        Partition::from_beta_numbers(&packed)
    }

    /// The **k-quotient** of λ: the `k` partitions read off the abacus runners,
    /// component `r` holding the beads with `β ≡ r (mod k)`.
    ///
    /// Together with [`k_core`](Self::k_core) this is the Littlewood
    /// decomposition: `|λ| = |k-core| + k · Σ_r |quotient_r|`. The component
    /// **order** (runner 0 first) is load-bearing downstream — LLT's tuple
    /// model is not symmetric in its components — and is what
    /// `llt::SkewTuple::quotient` is pinned against.
    pub fn k_quotient(&self, k: u32) -> Vec<Partition> {
        assert!(k >= 1, "a k-quotient needs k ≥ 1");
        let rows = self.k_rows(k);
        let beta = self.beta_numbers(rows);
        (0..k)
            .map(|r| {
                let mut pos: Vec<u32> = beta
                    .iter()
                    .filter(|&&b| b % k == r)
                    .map(|&b| b / k)
                    .collect();
                pos.sort_unstable_by(|a, b| b.cmp(a));
                Partition::from_beta_numbers(&pos)
            })
            .collect()
    }

    /// Does λ admit k-ribbon tableaux? Equivalently, is its k-core empty?
    ///
    /// The existence criterion for everything in [`crate::llt`]'s ribbon model.
    pub fn has_empty_k_core(&self, k: u32) -> bool {
        self.k_core(k).is_empty()
    }

    /// The order z_λ = ∏_i i^{m_i} · m_i! of the centralizer of a permutation of
    /// cycle type λ (m_i = multiplicity of the part i). Used for the power-sum
    /// normalization ⟨p_λ, p_λ⟩ = z_λ and for s ↔ p conversions.
    pub fn z(&self) -> u128 {
        fn factorial(m: u32) -> u128 {
            (1..=m as u128).product::<u128>().max(1)
        }
        let mut result: u128 = 1;
        let mut i = 0;
        while i < self.0.len() {
            let val = self.0[i] as u128;
            let mut mult = 0u32;
            while i < self.0.len() && self.0[i] as u128 == val {
                mult += 1;
                i += 1;
            }
            result *= val.pow(mult) * factorial(mult);
        }
        result
    }
}

impl fmt::Display for Partition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("[")?;
        for (i, p) in self.0.iter().enumerate() {
            if i > 0 {
                f.write_str(",")?;
            }
            write!(f, "{p}")?;
        }
        f.write_str("]")
    }
}

/// All partitions of `n`, generated in a canonical order (weakly decreasing,
/// parts bounded above so no duplicates). `partitions_of(0)` yields the empty
/// partition.
pub fn partitions_of(n: u32) -> Vec<Partition> {
    fn rec(n: u32, max: u32, cur: &mut Vec<u32>, out: &mut Vec<Partition>) {
        if n == 0 {
            out.push(Partition::from_sorted(cur.clone()));
            return;
        }
        let hi = n.min(max);
        for k in (1..=hi).rev() {
            cur.push(k);
            rec(n - k, k, cur, out);
            cur.pop();
        }
    }
    let mut out = Vec::new();
    let mut cur = Vec::new();
    rec(n, n.max(1), &mut cur, &mut out);
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalizes() {
        let p = Partition::new([1, 3, 0, 2, 2]);
        assert_eq!(p.parts(), &[3, 2, 2, 1]);
        assert_eq!(p.size(), 8);
        assert_eq!(p.len(), 4);
    }

    #[test]
    fn strict_validation() {
        assert!(Partition::try_new([3, 2, 2, 1]).is_ok());
        assert_eq!(
            Partition::try_new([2, 3]),
            Err(PartitionError::NotWeaklyDecreasing)
        );
        assert_eq!(Partition::try_new([2, 1, 0]), Err(PartitionError::ZeroPart));
    }

    #[test]
    fn containment() {
        let lam = Partition::new([3, 2, 1]);
        let mu = Partition::new([2, 1]);
        assert!(lam.contains(&mu));
        assert!(!mu.contains(&lam));
        assert!(lam.contains(&Partition::default())); // everything contains ∅
    }

    #[test]
    fn conjugate_is_involutive_and_correct() {
        assert_eq!(Partition::new([3, 1]).conjugate().parts(), &[2, 1, 1]);
        assert_eq!(Partition::new([4, 2, 1]).conjugate().parts(), &[3, 2, 1, 1]);
        for parts in [&[3, 1][..], &[4, 2, 1], &[2, 2], &[5]] {
            let p = Partition::new(parts.iter().copied());
            assert_eq!(p.conjugate().conjugate(), p);
        }
        assert_eq!(Partition::default().conjugate(), Partition::default());
    }

    #[test]
    fn centralizer_order_z() {
        assert_eq!(Partition::new([2, 1, 1]).z(), 4); // 2·(1²·2!) = 4
        assert_eq!(Partition::new([1, 1]).z(), 2);
        assert_eq!(Partition::new([2, 2]).z(), 8); // 2²·2! = 8
        assert_eq!(Partition::new([3]).z(), 3);
        assert_eq!(Partition::default().z(), 1);
    }

    #[test]
    fn beta_numbers_round_trip_at_every_padding() {
        for parts in [&[3, 1][..], &[4, 2, 1], &[2, 2], &[5], &[]] {
            let p = Partition::new(parts.iter().copied());
            for extra in 0..4 {
                let rows = p.len() + extra;
                let beta = p.beta_numbers(rows);
                assert!(
                    beta.windows(2).all(|w| w[0] > w[1]),
                    "β strictly decreasing"
                );
                assert_eq!(
                    Partition::from_beta_numbers(&beta),
                    p,
                    "{p} at {rows} beads"
                );
            }
        }
    }

    /// The Littlewood decomposition: `|λ| = |k-core| + k · Σ |quotient|`.
    #[test]
    fn core_and_quotient_split_the_size() {
        for n in 0..=9u32 {
            for lambda in partitions_of(n) {
                for k in 1..=4u32 {
                    let core = lambda.k_core(k);
                    let quot: u32 = lambda.k_quotient(k).iter().map(Partition::size).sum();
                    assert_eq!(
                        core.size() + k * quot,
                        n,
                        "Littlewood decomposition of {lambda} at k={k}"
                    );
                    // A core has no removable k-rim-hook, so it is its own core.
                    assert_eq!(core.k_core(k), core, "{lambda} core is a fixed point");
                }
            }
        }
    }

    /// The quotient must not depend on how many spare beads the abacus carries
    /// — only on the residues, which a multiple-of-k padding preserves.
    #[test]
    fn k_quotient_is_independent_of_the_padding() {
        for n in 0..=8u32 {
            for lambda in partitions_of(n) {
                for k in 1..=4u32 {
                    let want = lambda.k_quotient(k);
                    for extra in 1..=3usize {
                        let rows = lambda.len().max(1);
                        let rows =
                            k as usize * ((rows + k as usize) / k as usize) + extra * k as usize;
                        let beta = lambda.beta_numbers(rows);
                        let got: Vec<Partition> = (0..k)
                            .map(|r| {
                                let mut pos: Vec<u32> = beta
                                    .iter()
                                    .filter(|&&b| b % k == r)
                                    .map(|&b| b / k)
                                    .collect();
                                pos.sort_unstable_by(|a, b| b.cmp(a));
                                Partition::from_beta_numbers(&pos)
                            })
                            .collect();
                        assert_eq!(got, want, "{lambda} at k={k}, {rows} beads");
                    }
                }
            }
        }
    }

    /// k = 1 sees every cell as its own hook: empty core, quotient = λ itself.
    #[test]
    fn one_cores_are_empty() {
        for n in 0..=7u32 {
            for lambda in partitions_of(n) {
                assert!(lambda.has_empty_k_core(1), "{lambda}");
                assert_eq!(lambda.k_quotient(1), vec![lambda.clone()]);
            }
        }
    }

    #[test]
    fn partition_counts() {
        // p(n): 1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42
        let expected = [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42];
        for (n, &e) in expected.iter().enumerate() {
            assert_eq!(partitions_of(n as u32).len(), e, "p({n})");
        }
    }
}
