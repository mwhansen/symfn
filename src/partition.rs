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
    fn partition_counts() {
        // p(n): 1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42
        let expected = [1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42];
        for (n, &e) in expected.iter().enumerate() {
            assert_eq!(partitions_of(n as u32).len(), e, "p({n})");
        }
    }
}
