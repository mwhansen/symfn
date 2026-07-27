//! A cheap non-cryptographic hasher for the hot frontier maps.
//!
//! Extracted from `skew_lr`, which discovered both of the tunings below, so that
//! the other dynamic programs over integer keys — `convert::p_expand`'s β-mask
//! sweep in particular — get them too rather than each rediscovering SipHash's
//! cost separately.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

/// A multiply-xor-rotate hasher.
///
/// The frontier map is the hot data structure — several million short integer
/// keys per expansion — and SipHash's per-key setup dominates there. Keys are
/// small and structured, so a cheap mixing step is a large net win; measured at
/// roughly 1.3× on `s[8,7,6,5,4,3]²`. Not cryptographic, and nothing here is
/// exposed to adversarial input.
#[derive(Default)]
pub(crate) struct MixHasher(u64);

const MIX: u64 = 0x9e37_79b9_7f4a_7c15;

impl MixHasher {
    #[inline]
    fn add(&mut self, w: u64) {
        self.0 = (self.0.rotate_left(5) ^ w).wrapping_mul(MIX);
    }
}

impl Hasher for MixHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut chunks = bytes.chunks_exact(8);
        for c in &mut chunks {
            self.add(u64::from_le_bytes(c.try_into().unwrap()));
        }
        let rest = chunks.remainder();
        if !rest.is_empty() {
            let mut buf = [0u8; 8];
            buf[..rest.len()].copy_from_slice(rest);
            self.add(u64::from_le_bytes(buf));
        }
    }
    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.add(i as u64);
    }
    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.add(i as u64);
    }
    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.add(i);
    }
    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.add(i as u64);
    }
    #[inline]
    fn finish(&self) -> u64 {
        // Avalanche (the murmur3/splitmix-style finalizer). The per-chunk mix
        // above only moves entropy *upward* (multiply and a small left
        // rotation), which was fine when a key spanned many chunks and the
        // rotation wrapped, but byte-packed keys fit in 2–4 chunks and left
        // the low bits — exactly the ones hashbrown takes the bucket index
        // from — barely mixed. That showed up as probe-chain clustering worth
        // ~2.7× on `[20,16,12]²`; these two multiply-xorshift rounds fix the
        // dispersion for a few cycles per key.
        let mut x = self.0;
        x ^= x >> 33;
        x = x.wrapping_mul(0xff51_afd7_ed55_8ccd);
        x ^= x >> 33;
        x = x.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
        x ^ (x >> 33)
    }
}

pub(crate) type Map<K, V> = HashMap<K, V, BuildHasherDefault<MixHasher>>;
