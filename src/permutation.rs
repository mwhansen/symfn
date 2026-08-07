//! Permutations of ℕ₊ moving finitely many points — the index type for
//! [`Schubert`](crate::schubert::Schubert).
//!
//! A `Perm` is stored in one-line notation with **trailing fixed points
//! stripped at construction**, so `[2,1]`, `[2,1,3]` and `[2,1,3,4]` are all
//! the same value. That is the same move [`Partition`](crate::partition::Partition)
//! makes for weakly-decreasing, and for the same reason: downstream code gets
//! to rely on the normal form instead of re-deriving it.
//!
//! It is worth naming what this replaces. Schubert polynomials are stable
//! under `S_n ↪ S_{n+1}` — `S_w` does not change when you pad `w` with a fixed
//! point — so an implementation must decide *somewhere* that padded inputs are
//! equal. Symmetrica decided it in the comparator: `comp_permutation` compares
//! mixed-length lists as if fixed-point-padded, which means the invariant
//! lives in a function every caller must remember to route through, and the
//! identity permutation is spelled `[1,2]`. Here the invariant lives in the
//! constructor, the identity is the empty `Perm`, and `==` is `==`.
//!
//! # Index conventions
//!
//! **Everything here is 1-based**, matching the mathematics. This is a
//! deliberate stand: Symmetrica's own Schubert module is internally
//! inconsistent about it (`mult_schubert_variable` is 0-based,
//! `divdiff_schubert` is 1-based), and that inconsistency is a documented
//! source of silently-wrong results. [`Perm::at`] takes `i ≥ 1`, descent
//! positions are `1`-based, and transpositions name positions, not values.
//!
//! Positions past the stored prefix are not special-cased by callers: `at(i)`
//! returns `i` there, and the cover scans pad internally, because "the cover
//! that lands one past the end" is a real cover and forgetting it loses terms.

// Every `as` here converts a *position or a letter of a one-line word*: both
// are bounded by the permutation's length, which `Perm` stores as `u8` entries
// and rejects above that on construction. No coefficient passes through this
// module at all.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use core::fmt;

/// A permutation of ℕ₊ fixing all but finitely many points.
///
/// Stored as `w(1), …, w(m)` where `m` is the largest non-fixed point, so the
/// stored slice is always a permutation of `1..=m` with `w(m) ≠ m`. The empty
/// vector is the identity.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default)]
pub struct Perm {
    /// `w(1), …, w(m)` as `u8`, **zero-filled** past `len`.
    ///
    /// Inline rather than a `Vec` because this is the key of the engine's hot
    /// map: with heap storage every Bruhat cover allocated and every key
    /// comparison chased a pointer, and a sampling profile of `stair7²` put
    /// ~47% of the run in the allocator and `BTreeMap`
    /// (`docs/record/schubert.md`). Inline makes `Perm` `Copy`, makes
    /// comparison a 32-byte block compare, and makes the cover scans
    /// allocation-free.
    ///
    /// The zero fill is what keeps `Ord` identical to the old `Vec<u32>`
    /// lexicographic order — including the case where one permutation's
    /// one-line notation is a prefix of another's, since `0` sorts below every
    /// real value exactly as "ran out of elements" did.
    data: [u8; MAX_SUPPORT],
    len: u8,
}

/// Largest permutation this representation holds.
///
/// A product of two `S_n` Schuberts can reach support `2n − 1` (the `x₁`
/// degree of each factor is at most `n − 1`), so `S₁₃` inputs need ~25. 32
/// covers the measured ladder with room; [`Perm::new`] rejects anything larger
/// rather than silently truncating.
pub const MAX_SUPPORT: usize = 32;

/// Why a slice failed to be a valid one-line notation.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum PermError {
    /// A value was `0`; one-line notation is over `1..=n`.
    ZeroValue,
    /// A value repeated, or exceeded the length of the input.
    NotABijection,
    /// More than [`MAX_SUPPORT`] points.
    TooLarge,
}

impl Perm {
    /// The identity.
    pub fn identity() -> Self {
        Perm {
            data: [0; MAX_SUPPORT],
            len: 0,
        }
    }

    /// The stored prefix.
    #[inline]
    fn slice(&self) -> &[u8] {
        &self.data[..self.len as usize]
    }

    /// Build from a `1..=n` array, stripping trailing fixed points. The tail
    /// past the resulting length is re-zeroed, since `Ord` depends on it.
    #[inline]
    fn from_array(mut data: [u8; MAX_SUPPORT], mut len: usize) -> Self {
        while len > 0 && data[len - 1] as usize == len {
            len -= 1;
        }
        for slot in data.iter_mut().skip(len) {
            *slot = 0;
        }
        Perm {
            data,
            len: len as u8,
        }
    }

    /// Build from one-line notation, **validating** that the input is a
    /// permutation of `1..=n` and then stripping trailing fixed points.
    ///
    /// Trailing fixed points in the *input* are fine and expected — that is
    /// how stability arrives from a caller that padded to a fixed `n`.
    ///
    /// # Errors
    ///
    /// Returns [`PermError::TooLarge`] above [`MAX_SUPPORT`] points.
    /// Returns [`PermError::ZeroValue`] if a value is `0`.
    /// Returns [`PermError::NotABijection`] if a value repeats or exceeds the
    /// length of the input.
    pub fn new<I: IntoIterator<Item = u32>>(one_line: I) -> Result<Self, PermError> {
        let v: Vec<u32> = one_line.into_iter().collect();
        let n = v.len();
        if n > MAX_SUPPORT {
            return Err(PermError::TooLarge);
        }
        let mut seen = vec![false; n + 1];
        for &x in &v {
            if x == 0 {
                return Err(PermError::ZeroValue);
            }
            let x = x as usize;
            if x > n || seen[x] {
                return Err(PermError::NotABijection);
            }
            seen[x] = true;
        }
        Ok(Perm::from_valid(v))
    }

    /// Internal: strip trailing fixed points from an already-validated
    /// one-line vector.
    pub(crate) fn from_valid(v: Vec<u32>) -> Self {
        assert!(v.len() <= MAX_SUPPORT, "permutation exceeds MAX_SUPPORT");
        let mut data = [0u8; MAX_SUPPORT];
        for (slot, &x) in data.iter_mut().zip(v.iter()) {
            *slot = x as u8;
        }
        Perm::from_array(data, v.len())
    }

    /// `w(i)` for `i ≥ 1`, returning `i` beyond the stored prefix.
    ///
    /// # Panics
    ///
    /// Panics if `i == 0`. Positions are 1-based, so `w(0)` names no cell, and
    /// the check is unconditional rather than a `debug_assert` so that the
    /// message names the requirement. Without it `i - 1` underflows, which
    /// `overflow-checks` catches as "attempt to subtract with overflow" — loud,
    /// but pointing at arithmetic rather than at the convention the caller
    /// broke.
    ///
    /// It costs nothing measurable even though this is the innermost accessor
    /// on the Schubert path: the branch is never taken, and folds away wherever
    /// `i` is a loop index. The measurement is in `docs/record/schubert.md`.
    #[inline]
    pub fn at(&self, i: u32) -> u32 {
        assert!(i >= 1, "positions are 1-based; w(0) names no cell");
        match self.slice().get((i - 1) as usize) {
            Some(&x) => x as u32,
            None => i,
        }
    }

    /// The stored one-line prefix `w(1), …, w(m)`.
    #[inline]
    pub fn one_line(&self) -> &[u8] {
        self.slice()
    }

    /// `m`: the largest non-fixed point, `0` for the identity.
    ///
    /// Named `support_len` rather than `len` because `len` on a permutation
    /// reads like `ℓ(w)`, which is [`Perm::length`] and a different number.
    #[inline]
    pub fn support_len(&self) -> u32 {
        self.len as u32
    }

    /// Whether `w` moves no point — the identity, which is the empty
    /// `Perm`.
    #[inline]
    pub fn is_identity(&self) -> bool {
        self.len == 0
    }

    /// One-line notation padded to length `n` (no-op if already longer).
    pub fn padded(&self, n: u32) -> Vec<u32> {
        (1..=n.max(self.support_len()))
            .map(|i| self.at(i))
            .collect()
    }

    /// `ℓ(w)`, the number of inversions.
    ///
    /// Fixed points past the prefix contribute nothing: every stored value is
    /// `≤ m` and every value past it equals its own position `> m`, so no
    /// inversion can straddle the boundary.
    pub fn length(&self) -> u32 {
        let v = self.slice();
        let mut n = 0;
        for i in 0..v.len() {
            for j in (i + 1)..v.len() {
                if v[i] > v[j] {
                    n += 1;
                }
            }
        }
        n
    }

    /// The Lehmer code `c_i = #{ j > i : w(j) < w(i) }`, trailing zeros
    /// stripped.
    ///
    /// `|code| = ℓ(w)` and `x^{code(w)}` is the lex-**minimal** monomial of
    /// `S_w` — the triangularity the from-polynomial peel rests on.
    pub fn code(&self) -> Vec<u32> {
        let v = self.slice();
        let mut c = vec![0u32; v.len()];
        for i in 0..v.len() {
            for j in (i + 1)..v.len() {
                if v[j] < v[i] {
                    c[i] += 1;
                }
            }
        }
        while c.last() == Some(&0) {
            c.pop();
        }
        c
    }

    /// Inverse of [`Perm::code`]: `w(i)` is the `(c_i + 1)`-th smallest value
    /// not already used.
    ///
    /// **Total** — every `u32` slice is the code of a permutation, so there is
    /// no error case. The subtlety is the length: `c_i ≤ m − i` is forced, so
    /// the permutation needs `m = max(|c|, maxᵢ(i + c_i))` points, which is
    /// generally *more* than `|c|` because [`Perm::code`] strips trailing
    /// zeros. Bounding by `|c|` instead is the bug this signature invites —
    /// it rejects `(2,1)`, the code of `[3,2,1]`.
    ///
    /// # Panics
    ///
    /// Panics if the code needs more than [`MAX_SUPPORT`] points, which is the
    /// width of the inline one-line array.
    pub fn from_code(code: &[u32]) -> Self {
        let mut m = code.len() as u32;
        for (i, &c) in code.iter().enumerate() {
            m = m.max((i as u32) + 1 + c);
        }
        let mut avail: Vec<u32> = (1..=m).collect();
        let mut out = Vec::with_capacity(m as usize);
        for i in 0..m as usize {
            let c = code.get(i).copied().unwrap_or(0) as usize;
            out.push(avail.remove(c));
        }
        Perm::from_valid(out)
    }

    /// Descent positions: `i ≥ 1` with `w(i) > w(i+1)`.
    ///
    /// Only the stored prefix can contain one — position `m` compares
    /// `w(m) ≤ m` against `m+1`, which is never a descent.
    pub fn descents(&self) -> Vec<u32> {
        let v = self.slice();
        (1..v.len())
            .filter(|&i| v[i - 1] > v[i])
            .map(|i| i as u32)
            .collect()
    }

    /// The last descent, without building the descent list.
    ///
    /// [`Perm::descents`] heap-allocates a `Vec` and every caller on the
    /// transition path throws all of it away but the last entry.
    #[inline]
    pub fn last_descent(&self) -> Option<u32> {
        let v = self.slice();
        (1..v.len())
            .rev()
            .find(|&i| v[i - 1] > v[i])
            .map(|i| i as u32)
    }

    /// [`Perm::covers_left`] without allocating: hands each cover to `f`.
    ///
    /// The allocating form builds a `Vec<(u32, Perm)>` which callers then
    /// re-collect into a `Vec<Perm>` — two allocations per node to carry data
    /// that is consumed immediately.
    ///
    /// # Panics
    ///
    /// Panics if `i == 0` — positions are 1-based — or if `i` exceeds
    /// [`MAX_SUPPORT`].
    #[inline]
    pub fn for_each_cover_left(&self, i: u32, mut f: impl FnMut(u32, Perm)) {
        let n = self.support_len().max(i) as usize;
        debug_assert!(n <= MAX_SUPPORT);
        let mut base = [0u8; MAX_SUPPORT];
        for (k, slot) in base.iter_mut().enumerate().take(n) {
            *slot = self.at(k as u32 + 1) as u8;
        }
        let wi = base[(i - 1) as usize];
        let mut running_max = 0u8;
        for j in (1..i as usize).rev() {
            let wj = base[j - 1];
            if wj < wi && wj > running_max {
                running_max = wj;
                let mut u = base;
                u.swap((i - 1) as usize, j - 1);
                f(j as u32, Perm::from_array(u, n));
            }
        }
    }

    /// Dominant: the Lehmer code is weakly decreasing. Then `S_w = x^{code(w)}`
    /// is a single monomial.
    pub fn is_dominant(&self) -> bool {
        self.code().windows(2).all(|w| w[0] >= w[1])
    }

    /// Grassmannian: at most one descent. Returns the descent position, or
    /// `Some(0)` for the identity (no descent), `None` if there are two or
    /// more.
    ///
    /// For such `w`, `S_w = s_λ(x_1, …, x_k)` with `k` the descent and `λ` the
    /// code reversed — the case that dispatches into the LR engine.
    pub fn is_grassmannian(&self) -> Option<u32> {
        let d = self.descents();
        match d.len() {
            0 => Some(0),
            1 => Some(d[0]),
            _ => None,
        }
    }

    /// `w⁻¹`.
    pub fn inverse(&self) -> Self {
        let v = self.slice();
        let mut data = [0u8; MAX_SUPPORT];
        for (i, &x) in v.iter().enumerate() {
            data[(x - 1) as usize] = (i + 1) as u8;
        }
        Perm::from_array(data, v.len())
    }

    /// `w · t_{ij}`: swap the *positions* `i` and `j` (1-based).
    ///
    /// # Panics
    ///
    /// Panics if `i` or `j` is `0` — positions are 1-based — or if the larger
    /// of them exceeds [`MAX_SUPPORT`], which is the width of the inline
    /// one-line array and so the widest symmetric group this type represents.
    pub fn transpose(&self, i: u32, j: u32) -> Self {
        assert!(
            i >= 1 && j >= 1,
            "positions are 1-based; t_{{0j}} names no cell"
        );
        let n = i.max(j).max(self.support_len()) as usize;
        assert!(n <= MAX_SUPPORT, "transposition exceeds MAX_SUPPORT");
        let mut data = [0u8; MAX_SUPPORT];
        for (k, slot) in data.iter_mut().enumerate().take(n) {
            *slot = self.at(k as u32 + 1) as u8;
        }
        data.swap((i - 1) as usize, (j - 1) as usize);
        Perm::from_array(data, n)
    }

    /// The Bruhat covers `w ⋖ w·t_{ij}` with `j > i`: positions `j` where
    /// `w(i) < w(j)` and nothing strictly between them sits in between in
    /// value.
    ///
    /// This is the running-minimum scan Symmetrica uses in both
    /// `mult_schubert_variable` and `algorithmus2`, and it is correct to
    /// ignore values below `w(i)` while updating the running minimum: they can
    /// never block a cover.
    ///
    /// **At most one cover lands past the stored prefix.** Beyond `m` the
    /// values are `m+1, m+2, …` in order, so `j = m+1` may be a cover but
    /// `j = m+2` is always blocked by `j = m+1`. Scanning to `m+1` is
    /// therefore exhaustive, not a heuristic cutoff — the thing that would
    /// silently drop terms if it were guessed.
    /// # Panics
    ///
    /// Panics if `i == 0` — positions are 1-based — or if the scan, which
    /// reaches one past the stored prefix, would exceed [`MAX_SUPPORT`].
    pub fn covers_right(&self, i: u32) -> Vec<(u32, Perm)> {
        assert!(i >= 1, "positions are 1-based; w(0) names no cell");
        let n = (self.support_len().max(i) + 1) as usize;
        assert!(n <= MAX_SUPPORT, "cover scan exceeds MAX_SUPPORT");
        let mut base = [0u8; MAX_SUPPORT];
        for (k, slot) in base.iter_mut().enumerate().take(n) {
            *slot = self.at(k as u32 + 1) as u8;
        }
        let wi = base[(i - 1) as usize];
        let mut out = Vec::new();
        let mut running_min = u8::MAX;
        for j in (i as usize + 1)..=n {
            let wj = base[j - 1];
            if wj > wi && wj < running_min {
                running_min = wj;
                let mut u = base;
                u.swap((i - 1) as usize, j - 1);
                out.push((j as u32, Perm::from_array(u, n)));
            }
        }
        out
    }

    /// The Bruhat covers `w ⋖ w·t_{ji}` with `j < i` — the left scan, a
    /// running **maximum** walking leftwards from `i-1`.
    ///
    /// Bounded by construction (there is no room to the left of position 1),
    /// which is why this side needs no padding argument while
    /// [`Perm::covers_right`] does.
    /// # Panics
    ///
    /// Panics if `i == 0` — positions are 1-based — or if `i` exceeds
    /// [`MAX_SUPPORT`].
    pub fn covers_left(&self, i: u32) -> Vec<(u32, Perm)> {
        assert!(i >= 1, "positions are 1-based; w(0) names no cell");
        let n = self.support_len().max(i) as usize;
        assert!(n <= MAX_SUPPORT, "cover scan exceeds MAX_SUPPORT");
        let mut base = [0u8; MAX_SUPPORT];
        for (k, slot) in base.iter_mut().enumerate().take(n) {
            *slot = self.at(k as u32 + 1) as u8;
        }
        let wi = base[(i - 1) as usize];
        let mut out = Vec::new();
        let mut running_max = 0u8;
        for j in (1..i as usize).rev() {
            let wj = base[j - 1];
            if wj < wi && wj > running_max {
                running_max = wj;
                let mut u = base;
                u.swap((i - 1) as usize, j - 1);
                out.push((j as u32, Perm::from_array(u, n)));
            }
        }
        out
    }

    /// Bruhat order: `self ≤ other`, by the tableau criterion — for every `i`,
    /// the sorted prefix `{w(1),…,w(i)}` is entrywise ≤ the sorted prefix of
    /// `other`.
    ///
    /// This is the pruning predicate for a single-coefficient query. The
    /// signed Monk rule moves every term *up* the Bruhat order, since both its
    /// sums run over covers. So every permutation reachable from `z` is `≥ z`.
    /// If `z ≰ w`, no descendant of `z` can equal `w`, and `z` may be dropped
    /// without affecting `c^w`. Dropping is safe even though terms cancel,
    /// since a dropped term contributes to `w` neither directly nor through a
    /// cancellation with something that does.
    pub fn bruhat_le(&self, other: &Perm) -> bool {
        let n = self.support_len().max(other.support_len()) as usize;
        let (mut a, mut b) = ([0u8; MAX_SUPPORT], [0u8; MAX_SUPPORT]);
        for i in 0..n {
            // insert self.at(i+1) / other.at(i+1) into sorted prefixes
            let (x, y) = (self.at(i as u32 + 1) as u8, other.at(i as u32 + 1) as u8);
            let mut j = i;
            while j > 0 && a[j - 1] > x {
                a[j] = a[j - 1];
                j -= 1;
            }
            a[j] = x;
            let mut j = i;
            while j > 0 && b[j - 1] > y {
                b[j] = b[j - 1];
                j -= 1;
            }
            b[j] = y;
            for k in 0..=i {
                if a[k] > b[k] {
                    return false;
                }
            }
        }
        true
    }

    /// A reduced word: positions `i₁, …, i_ℓ` with `w = s_{i₁} ⋯ s_{i_ℓ}` and
    /// `ℓ = ℓ(w)`.
    ///
    /// Produced by repeatedly cancelling the first descent, so the word is
    /// canonical rather than arbitrary — but nothing may depend on *which*
    /// reduced word it is, since `∂` and the nil-Hecke generators satisfy the
    /// braid relations. The tests check that independence directly instead of
    /// trusting it.
    ///
    /// Not on any hot path: `∂_w` composes one-letter steps, and Macdonald's
    /// reduced-word checksum enumerates *all* words, which is a test-only job
    /// for Sage.
    pub fn reduced_word(&self) -> Vec<u32> {
        let mut v = self.slice().to_vec();
        let mut word = Vec::new();
        while let Some(i) = (1..v.len()).find(|&i| v[i - 1] > v[i]) {
            v.swap(i - 1, i);
            word.push(i as u32);
        }
        debug_assert_eq!(word.len() as u32, self.length());
        word.reverse();
        word
    }

    /// The transition step of Lascoux–Schützenberger, in the form Symmetrica's
    /// `newtrans` uses: `r` = the **last** descent, `s` = the **largest**
    /// `s > r` with `w(s) < w(r)`, `v = w · t_{rs}` (so `ℓ(v) = ℓ(w) − 1`).
    ///
    /// Returns `(r, v, [v·t_{qr} : q < r, a cover])`, giving
    ///
    /// ```text
    ///     S_w = x_r · S_v + Σ S_{v t_{qr}}
    /// ```
    ///
    /// `None` for the identity, which has no descent.
    pub fn transition(&self) -> Option<(u32, Perm, Vec<Perm>)> {
        let r = self.last_descent()?;
        let wr = self.at(r);
        // largest s > r with w(s) < w(r); bounded by the prefix, since past it
        // w(s) = s > m >= w(r).
        let s = (r + 1..=self.support_len()).rfind(|&s| self.at(s) < wr)?;
        let v = self.transpose(r, s);
        debug_assert_eq!(v.length() + 1, self.length());
        let ups = v.covers_left(r).into_iter().map(|(_, p)| p).collect();
        Some((r, v, ups))
    }
}

impl fmt::Display for Perm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_identity() {
            return write!(f, "id");
        }
        write!(f, "[")?;
        for (i, x) in self.slice().iter().enumerate() {
            if i > 0 {
                write!(f, ",")?;
            }
            write!(f, "{x}")?;
        }
        write!(f, "]")
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    fn p(v: &[u32]) -> Perm {
        Perm::new(v.iter().copied()).unwrap()
    }

    #[test]
    fn strips_trailing_fixed_points() {
        assert_eq!(p(&[2, 1]), p(&[2, 1, 3]));
        assert_eq!(p(&[2, 1]), p(&[2, 1, 3, 4, 5]));
        assert_eq!(p(&[1, 2, 3]), Perm::identity());
        assert!(p(&[1, 2, 3]).is_identity());
        // an interior fixed point is NOT stripped
        assert_eq!(p(&[3, 2, 1]).one_line(), &[3, 2, 1]);
    }

    #[test]
    fn rejects_non_permutations() {
        assert_eq!(Perm::new([1, 1]), Err(PermError::NotABijection));
        assert_eq!(Perm::new([0, 1]), Err(PermError::ZeroValue));
        assert_eq!(Perm::new([1, 3]), Err(PermError::NotABijection));
    }

    #[test]
    fn at_is_one_based_and_total() {
        let w = p(&[2, 1]);
        assert_eq!(w.at(1), 2);
        assert_eq!(w.at(2), 1);
        assert_eq!(w.at(3), 3);
        assert_eq!(w.at(97), 97);
    }

    #[test]
    fn length_and_code_agree() {
        for n in 0..=6u32 {
            for w in all_perms(n) {
                assert_eq!(w.code().iter().sum::<u32>(), w.length(), "{w}");
            }
        }
    }

    #[test]
    fn code_round_trips() {
        for n in 0..=7u32 {
            for w in all_perms(n) {
                assert_eq!(Perm::from_code(&w.code()), w, "{w}");
            }
        }
    }

    #[test]
    fn code_round_trips_through_padding() {
        // the stability invariant: a padded code is the same permutation
        let w = p(&[1, 4, 2, 3]);
        let mut c = w.code();
        c.extend_from_slice(&[0, 0, 0]);
        assert_eq!(Perm::from_code(&c), w);
    }

    #[test]
    fn inverse_is_an_involution() {
        for n in 0..=6u32 {
            for w in all_perms(n) {
                assert_eq!(w.inverse().inverse(), w, "{w}");
                assert_eq!(w.inverse().length(), w.length(), "{w}");
                for i in 1..=n {
                    assert_eq!(w.inverse().at(w.at(i)), i, "{w} at {i}");
                }
            }
        }
    }

    #[test]
    fn dominant_and_grassmannian_families() {
        assert!(p(&[3, 2, 1]).is_dominant()); // code (2,1)
        assert!(!p(&[1, 4, 2, 3]).is_dominant()); // code (0,2)
        assert_eq!(p(&[1, 3, 2]).is_grassmannian(), Some(2));
        assert_eq!(p(&[2, 4, 1, 3]).is_grassmannian(), Some(2));
        assert_eq!(p(&[2, 1, 4, 3]).is_grassmannian(), None);
        assert_eq!(Perm::identity().is_grassmannian(), Some(0));
        // w0 is dominant for every n
        for n in 2..=7u32 {
            let w0: Vec<u32> = (1..=n).rev().collect();
            assert!(p(&w0).is_dominant());
            assert_eq!(p(&w0).length(), n * (n - 1) / 2);
        }
    }

    #[test]
    fn covers_really_are_covers() {
        for n in 0..=6u32 {
            for w in all_perms(n) {
                for i in 1..=(n + 1) {
                    for (_, u) in w.covers_right(i) {
                        assert_eq!(u.length(), w.length() + 1, "{w} right {i} -> {u}");
                    }
                    for (_, u) in w.covers_left(i) {
                        assert_eq!(u.length(), w.length() + 1, "{w} left {i} -> {u}");
                    }
                }
            }
        }
    }

    /// The cover scans must find *every* cover, including the one landing past
    /// the stored prefix — that is the term a naive `for j in i+1..=m` loses.
    #[test]
    fn cover_scans_are_exhaustive_against_brute_force() {
        for n in 0..=5u32 {
            for w in all_perms(n) {
                for i in 1..=(n + 1) {
                    // brute force with generous headroom
                    let bound = n + 4;
                    let mut want_r: Vec<Perm> = Vec::new();
                    for j in (i + 1)..=bound {
                        let u = w.transpose(i, j);
                        if u.length() == w.length() + 1 {
                            want_r.push(u);
                        }
                    }
                    let got_r: Vec<Perm> = w.covers_right(i).into_iter().map(|x| x.1).collect();
                    assert_eq!(got_r, want_r, "right covers of {w} at {i}");

                    let mut want_l: Vec<Perm> = Vec::new();
                    for j in (1..i).rev() {
                        let u = w.transpose(j, i);
                        if u.length() == w.length() + 1 {
                            want_l.push(u);
                        }
                    }
                    let got_l: Vec<Perm> = w.covers_left(i).into_iter().map(|x| x.1).collect();
                    assert_eq!(got_l, want_l, "left covers of {w} at {i}");
                }
            }
        }
    }

    #[test]
    fn transition_drops_length_by_one() {
        for n in 0..=6u32 {
            for w in all_perms(n) {
                match w.transition() {
                    None => assert!(w.is_identity(), "{w}"),
                    Some((r, v, ups)) => {
                        assert_eq!(v.length() + 1, w.length(), "{w}");
                        assert!(r >= 1);
                        for u in &ups {
                            assert_eq!(u.length(), v.length() + 1, "{w} -> {u}");
                        }
                    }
                }
            }
        }
    }

    pub(crate) fn all_perms(n: u32) -> Vec<Perm> {
        let mut out = Vec::new();
        let mut cur: Vec<u32> = (1..=n).collect();
        permute(&mut cur, 0, &mut out);
        out
    }

    fn permute(cur: &mut Vec<u32>, k: usize, out: &mut Vec<Perm>) {
        if k == cur.len() {
            out.push(Perm::new(cur.iter().copied()).unwrap());
            return;
        }
        for i in k..cur.len() {
            cur.swap(k, i);
            permute(cur, k + 1, out);
            cur.swap(k, i);
        }
    }
}
