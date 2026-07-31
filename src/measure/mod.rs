//! Heap accounting, shared by the benchmarks, the budget tests, and `heapstat`.
//!
//! ## Why this exists rather than a profiler
//!
//! The three quantities that get called "memory" move independently, and only
//! two of them are reproducible:
//!
//! * **peak live bytes** — what a computation holds at its high-water mark,
//! * **total bytes** — everything ever handed out ("churn"),
//! * **peak RSS** — pages held from the OS, which depends on the allocator's
//!   whole allocation *history* and therefore drifts between identical runs.
//!
//! Measured here: a single-threaded workload reports the first two **bit-exactly
//! across runs** (the degree-10 (q,t)-Kostka table gives 16.5 MB / 786.9 MB /
//! 729 802 allocations every time), while its RSS moves run to run. That is what
//! makes the first two usable as *test assertions* and RSS not. Workloads that
//! fan out across threads vary by about 0.1% in allocation count and 1% in peak,
//! because the shard split depends on thread scheduling — hence the tolerance
//! band in [`Budget`].
//!
//! See `docs/record/memory.md` for the rules these
//! numbers feed into, in particular why halving churn is not on its own a
//! memory improvement.
//!
//! ## Using it
//!
//! A binary that wants to measure installs the allocator once:
//!
//! ```ignore
//! #[global_allocator]
//! static A: symfn::measure::Counting = symfn::measure::Counting::new();
//! ```
//!
//! and then wraps the work:
//!
//! ```ignore
//! let (result, stats) = symfn::measure::measure(|| qt_kostka_table::<i128>(10));
//! println!("{stats}");
//! ```
//!
//! Nothing here costs anything in a binary that does not install `Counting` —
//! it is a type and some statics, and the library itself never calls into it.
//!
//! Accounting is **process-wide**, so one measurement runs at a time; that is
//! deliberate, since the paths worth measuring spawn threads and a thread-local
//! counter would silently miss their allocations.

pub mod workloads;

use std::alloc::{GlobalAlloc, Layout, System};
use std::fmt;
use std::sync::atomic::{AtomicIsize, AtomicUsize, Ordering::Relaxed};

/// Live bytes, **signed** — and it genuinely goes negative. [`reset`] zeroes it
/// while memory allocated before the measurement is still held (the test
/// harness's own capture buffer is the reliable example), and every one of
/// those blocks is freed against a counter that no longer counts it. Unsigned,
/// that underflow wrapped to ~2^64 and `PEAK` latched it: a memory budget
/// reading whatever the first stale free happened to produce. `overflow-checks`
/// turned the wrap into a panic, which is how it was found.
static LIVE: AtomicIsize = AtomicIsize::new(0);
/// High-water mark of [`LIVE`], signed for the same reason and floored at 0 on
/// the way out ([`snapshot`]).
static PEAK: AtomicIsize = AtomicIsize::new(0);
static TOTAL: AtomicUsize = AtomicUsize::new(0);
static COUNT: AtomicUsize = AtomicUsize::new(0);
/// Allocation count and bytes by size class; bucket `k` covers `[2^k, 2^{k+1})`.
static BUCKET_COUNT: [AtomicUsize; 48] = [const { AtomicUsize::new(0) }; 48];
static BUCKET_BYTES: [AtomicUsize; 48] = [const { AtomicUsize::new(0) }; 48];

/// A [`GlobalAlloc`] that forwards to the system allocator and counts.
///
/// Counting is four relaxed atomic adds per allocation, which perturbs timing
/// enough that a binary with this installed should not be used for wall-clock
/// benchmarks — measure memory and time in separate runs.
pub struct Counting;

impl Default for Counting {
    fn default() -> Self {
        Self::new()
    }
}

impl Counting {
    pub const fn new() -> Self {
        Counting
    }
}

#[inline]
fn record(size: usize) {
    TOTAL.fetch_add(size, Relaxed);
    COUNT.fetch_add(1, Relaxed);
    let k = ((usize::BITS - size.max(1).leading_zeros()) as usize - 1).min(47);
    BUCKET_COUNT[k].fetch_add(1, Relaxed);
    BUCKET_BYTES[k].fetch_add(size, Relaxed);
}

#[inline]
fn grew(by: usize) {
    // `by` is a `Layout` size, so it is at most `isize::MAX` by that type's own
    // invariant, and the sum is bounded by the address space.
    let by = by as isize;
    let now = LIVE.fetch_add(by, Relaxed) + by;
    PEAK.fetch_max(now, Relaxed);
}

#[inline]
fn shrank(by: usize) {
    LIVE.fetch_sub(by as isize, Relaxed);
}

unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        grew(l.size());
        record(l.size());
        unsafe { System.alloc(l) }
    }
    unsafe fn alloc_zeroed(&self, l: Layout) -> *mut u8 {
        grew(l.size());
        record(l.size());
        unsafe { System.alloc_zeroed(l) }
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        shrank(l.size());
        unsafe { System.dealloc(p, l) }
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, new: usize) -> *mut u8 {
        if new > l.size() {
            grew(new - l.size());
            record(new - l.size());
        } else {
            shrank(l.size() - new);
            COUNT.fetch_add(1, Relaxed);
        }
        unsafe { System.realloc(p, l, new) }
    }
}

/// What one measured region allocated.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct Stats {
    /// High-water mark of live bytes — the number a memory budget is about.
    pub peak: usize,
    /// Every byte handed out, including what was freed again.
    pub total: usize,
    /// Number of allocation calls. Drives fragmentation, and so the gap between
    /// `peak` and RSS, far more than `total` does.
    pub allocs: usize,
}

impl Stats {
    /// Total over peak: how many times the live set was churned through.
    ///
    /// High churn is only a *memory* problem when the sizes are diverse (so the
    /// allocator cannot recycle a block into the next request) or when the
    /// buffers are retained. Uniform, promptly-freed churn is a CPU cost.
    pub fn churn(&self) -> f64 {
        self.total as f64 / (self.peak.max(1) as f64)
    }
}

impl fmt::Display for Stats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        const MB: f64 = 1048576.0;
        write!(
            f,
            "peak {:>8.1} MB   total {:>9.1} MB   allocs {:>11}   churn {:>4.1}x",
            self.peak as f64 / MB,
            self.total as f64 / MB,
            self.allocs,
            self.churn(),
        )
    }
}

/// Run `f` with the counters zeroed, and report what it allocated.
///
/// Clears the memo tables first, so a measurement never depends on what ran
/// before it — without that, the second workload in a process measures a warm
/// cache and reports a fraction of its true cost.
///
/// Process-wide: do not call this on two threads at once.
pub fn measure<T>(f: impl FnOnce() -> T) -> (T, Stats) {
    crate::clear_caches();
    reset();
    let out = f();
    (out, snapshot())
}

/// Zero every counter.
pub fn reset() {
    LIVE.store(0, Relaxed);
    PEAK.store(0, Relaxed);
    TOTAL.store(0, Relaxed);
    COUNT.store(0, Relaxed);
    for k in 0..48 {
        BUCKET_COUNT[k].store(0, Relaxed);
        BUCKET_BYTES[k].store(0, Relaxed);
    }
}

/// Read the counters without disturbing them.
pub fn snapshot() -> Stats {
    Stats {
        // Floored at 0: a measurement that only ever *freed* pre-existing memory
        // has a negative high-water mark, and "it allocated nothing" is the
        // honest reading of that. Non-negative, so the cast cannot change value.
        peak: PEAK.load(Relaxed).max(0) as usize,
        total: TOTAL.load(Relaxed),
        allocs: COUNT.load(Relaxed),
    }
}

/// Allocation count and bytes per power-of-two size class, smallest first, as
/// `(size, count, bytes)` and skipping empty classes.
///
/// This is what attributes churn to a *specific buffer*: a spike in the 1–8 KB
/// classes is polynomial arithmetic, a spike at 32 bytes is one `Partition` per
/// output term.
pub fn histogram() -> Vec<(usize, usize, usize)> {
    (0..48)
        .filter_map(|k| {
            let n = BUCKET_COUNT[k].load(Relaxed);
            (n > 0).then(|| (1usize << k, n, BUCKET_BYTES[k].load(Relaxed)))
        })
        .collect()
}

/// A ceiling on what a workload may allocate, checked by `tests/memory.rs`.
///
/// Budgets are **upper bounds with headroom**, not golden values: they catch a
/// regression that changes the shape of the memory use, and stay quiet for the
/// noise. `tolerance` exists because workloads that fan out across threads are
/// not bit-reproducible — the shard split follows thread scheduling, worth about
/// 0.1% on allocation count and 1% on peak. Single-threaded workloads repeat
/// exactly and can be pinned tight.
#[derive(Clone, Copy, Debug)]
pub struct Budget {
    pub name: &'static str,
    /// Ceiling on peak live bytes.
    pub peak: usize,
    /// Ceiling on allocation count.
    pub allocs: usize,
    /// Fraction over the ceiling that still passes, e.g. `0.05` for 5%.
    pub tolerance: f64,
}

impl Budget {
    /// `Err(explanation)` if `stats` breaks the budget.
    ///
    /// The message carries the numbers to paste back in, because the correct
    /// response to a failure is often "this change legitimately costs more" and
    /// the point is to make *deciding that* explicit rather than laborious.
    pub fn check(&self, stats: &Stats) -> Result<(), String> {
        let scale = 1.0 + self.tolerance;
        let peak_max = (self.peak as f64 * scale) as usize;
        let allocs_max = (self.allocs as f64 * scale) as usize;
        if stats.peak <= peak_max && stats.allocs <= allocs_max {
            return Ok(());
        }
        const MB: f64 = 1048576.0;
        Err(format!(
            "memory budget exceeded for `{}`\n  \
             peak   {:.1} MB  (budget {:.1} MB, +{:.0}% = {:.1} MB)\n  \
             allocs {}  (budget {}, +{:.0}% = {})\n  \
             If the increase is intended, update the budget to:\n    \
             Budget {{ name: \"{}\", peak: {}, allocs: {}, tolerance: {} }}",
            self.name,
            stats.peak as f64 / MB,
            self.peak as f64 / MB,
            self.tolerance * 100.0,
            peak_max as f64 / MB,
            stats.allocs,
            self.allocs,
            self.tolerance * 100.0,
            allocs_max,
            self.name,
            stats.peak,
            stats.allocs,
            self.tolerance,
        ))
    }
}
