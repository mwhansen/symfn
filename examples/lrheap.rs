//! Peak *heap* accounting for one skew expansion, via a counting allocator.
//! RSS conflates live bytes with allocator retention; this separates them.
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicUsize, Ordering};
use symfn::{clear_caches, skew_lr::{expand_skew, take_peak_frontier_states}, Partition};

static LIVE: AtomicUsize = AtomicUsize::new(0);
static PEAK: AtomicUsize = AtomicUsize::new(0);

struct Counting;
unsafe impl GlobalAlloc for Counting {
    unsafe fn alloc(&self, l: Layout) -> *mut u8 {
        let n = LIVE.fetch_add(l.size(), Ordering::Relaxed) + l.size();
        PEAK.fetch_max(n, Ordering::Relaxed);
        unsafe { System.alloc(l) }
    }
    unsafe fn dealloc(&self, p: *mut u8, l: Layout) {
        LIVE.fetch_sub(l.size(), Ordering::Relaxed);
        unsafe { System.dealloc(p, l) }
    }
    unsafe fn realloc(&self, p: *mut u8, l: Layout, new: usize) -> *mut u8 {
        if new > l.size() {
            let n = LIVE.fetch_add(new - l.size(), Ordering::Relaxed) + (new - l.size());
            PEAK.fetch_max(n, Ordering::Relaxed);
        } else {
            LIVE.fetch_sub(l.size() - new, Ordering::Relaxed);
        }
        unsafe { System.realloc(p, l, new) }
    }
}
#[global_allocator]
static A: Counting = Counting;

fn p(v: &[u32]) -> Partition { Partition::new(v.iter().copied()) }

fn main() {
    let which = std::env::args().nth(1).unwrap_or_else(|| "0".into()).parse::<usize>().unwrap();
    let shapes: [&[u32]; 4] = [&[10,8,6,4], &[12,10,8,6], &[8,7,6,5,4,3], &[16,13,10,7]];
    let mu = p(shapes[which]);
    let w = mu.part(0);
    let outer: Vec<u32> = mu.parts().iter().map(|x| x + w).chain(mu.parts().iter().copied()).collect();
    let inner: Vec<u32> = mu.parts().iter().map(|_| w).chain(std::iter::repeat(0).take(mu.len())).collect();
    clear_caches();
    let _ = take_peak_frontier_states();
    LIVE.store(0, Ordering::Relaxed);
    PEAK.store(0, Ordering::Relaxed);
    let r = expand_skew(&p(&outer), &Partition::new(inner.into_iter()));
    let states = take_peak_frontier_states();
    let peak = PEAK.load(Ordering::Relaxed);
    eprintln!("{mu}^2  {} terms  peak states {states}  peak heap {:.1} MB  = {:.0} bytes/state",
        r.len(), peak as f64 / 1048576.0, peak as f64 / states as f64);
}
