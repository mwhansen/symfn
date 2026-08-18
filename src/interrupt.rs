//! Cancellation: how a caller stops a computation that is already running.
//!
//! A long call into this library is one uninterruptible unit of work to the
//! program that made it. That is invisible from Rust, where the caller chose
//! to block, and it is a defect from an embedder that handles signals itself:
//! CPython records SIGINT and runs the handler at its next bytecode boundary,
//! so a Ctrl-C during a multi-minute conversion is delivered only once the
//! conversion returns. The user sees a session that ignores Ctrl-C
//! (`docs/record/python-and-sage-interop.md`).
//!
//! The kernel cannot ask CPython anything — it does not depend on it
//! (`docs/policies/python.md`, P3). So the embedder installs a checker and the
//! kernel calls it:
//!
//! ```text
//!   set_checker(f)   ->  poll() consults f     ->  panics with Interrupted
//!   catch_interrupt(|| .. compute ..)  ->  Ok(answer)  or  Err(Interrupted)
//! ```
//!
//! [`poll`] is called from loops whose trip count grows with the input. Where
//! nothing is installed it is one relaxed load and a predictable branch, which
//! is what keeps it affordable in a hot loop.
//!
//! ## Why a panic and not a return value
//!
//! The same constraint that produced the overflow counter in
//! [`guard`](crate::guard) produces this: `Ring::mul` returns `Self`, so
//! generic code has no per-operation failure channel to thread a cancellation
//! through. Every alternative reaches the same place — either every loop that
//! can run long returns `Option`, and every caller of those loops does too, or
//! the stack is unwound. Unwinding is the smaller change and the one that
//! reaches the generic interiors at all.
//!
//! This is the one panic in the crate that is **not** a violated contract
//! (`docs/policies/failure.md`, R2 and the cancellation row of its table). It
//! is raised only when a checker the embedder installed says so, it carries a
//! payload no other panic uses, and the embedder that installed the checker is
//! the one that catches it. A caller who installs nothing cannot reach it.
//!
//! Two things in the crate already behave correctly under it, and neither was
//! written for this:
//!
//! - `escalate` in [`python`](crate::python) tests for `None`, so a panic
//!   passes through it rather than being read as overflow. An interrupted
//!   fixed-width pass does not silently restart over `BigInt`.
//! - `src/memo.rs` runs `compute` outside its guard and clears lock poison
//!   rather than propagating it, so an unwind through a cache leaves no
//!   half-built entry and does not disable the table for later calls.
//!
//! ## Which thread may poll
//!
//! Only the thread that entered the library. A cancellation unwinds the thread
//! that saw it, and a worker's unwind is not the caller's: the parallel
//! Littlewood–Richardson fill joins its scoped workers with an `expect` that
//! reads any worker panic as a bug, so a cancellation raised inside one would
//! be reported as that bug instead of as itself. The checker is also the
//! embedder's code — under CPython it takes the GIL, which a worker cannot
//! assume it may do. So [`poll`] goes on sequential driver loops, and a
//! parallel section is cancelled at the boundary that dispatched it.
//!
//! ## What a checker owes
//!
//! It runs at a point the kernel chose, with the crate's locks released. It
//! must not call back into this library, and it must be cheap enough to run
//! repeatedly — [`poll`] throttles the calls but states no rate. Returning
//! `true` more than once is harmless: the first one already unwound the thread
//! that saw it.

use std::cell::Cell;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::RwLock;

/// The payload of a cancellation panic, and the error [`catch_interrupt`]
/// reports.
///
/// It carries nothing. What a caller wants to know — which signal arrived,
/// what to raise — is the embedder's, and it is the embedder that held that
/// information before it installed the checker.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Interrupted;

/// Set once at install so [`poll`] can decline in one relaxed load, rather
/// than taking the lock to find the checker absent.
static ARMED: AtomicBool = AtomicBool::new(false);

static CHECKER: RwLock<Option<fn() -> bool>> = RwLock::new(None);

// Polls left before the next call to the checker. Per thread because the
// counter is only a rate limit: two threads sharing one would make each check
// the other's business and add a contended write to every poll.
thread_local! {
    static UNTIL_CHECK: Cell<u32> = const { Cell::new(STRIDE) };
}

/// How many polls pass between two calls to the checker.
///
/// The checker is the expensive half — under CPython it reacquires the GIL and
/// runs pending signal handlers — so [`poll`] is placed for a useful bound on
/// *latency* and this divides down the cost of getting one. It is a count and
/// not a duration because a clock read is itself the thing being avoided.
///
/// Public because a test asserting how promptly a cancellation lands has to
/// state the bound in the same unit the throttle uses; a number repeated in a
/// test is a number that drifts.
pub const STRIDE: u32 = 64;

/// Install the checker consulted by [`poll`], replacing any previous one.
///
/// `check` returns `true` to cancel. It is a plain function pointer rather
/// than a closure so that storing it needs no allocation and no type
/// parameter on a `static`; an embedder that needs state reaches it the same
/// way it reached this call.
pub fn set_checker(check: fn() -> bool) {
    *CHECKER.write().unwrap_or_else(|e| e.into_inner()) = Some(check);
    ARMED.store(true, Ordering::Relaxed);
}

/// Remove the checker, returning [`poll`] to a no-op.
///
/// Disarms before dropping the pointer, so a poll racing this call either
/// sees the old checker or declines, and never finds `ARMED` set with nothing
/// behind it.
pub fn clear_checker() {
    ARMED.store(false, Ordering::Relaxed);
    *CHECKER.write().unwrap_or_else(|e| e.into_inner()) = None;
}

/// Ask the installed checker, at most once per [`STRIDE`] calls, whether the
/// computation should stop — and unwind the thread if it should.
///
/// Call this from loops whose trip count grows with the input, at a boundary
/// crossed often enough to bound how long a cancellation waits. It is not a
/// safepoint in any stronger sense: it does no bookkeeping and leaves nothing
/// for a caller to clean up.
///
/// # Panics
///
/// Panics with [`Interrupted`] when the checker returns `true`. That panic is
/// this module's return path and is meant to be caught by
/// [`catch_interrupt`], not by an intermediate frame.
#[inline]
pub fn poll() {
    if !ARMED.load(Ordering::Relaxed) {
        return;
    }
    let due = UNTIL_CHECK.with(|c| {
        let n = c.get().wrapping_sub(1);
        c.set(if n == 0 { STRIDE } else { n });
        n == 0
    });
    if due {
        consult();
    }
}

/// Split from [`poll`] so the common path — decline, or decrement — inlines
/// without the call to the checker behind it.
#[cold]
#[inline(never)]
fn consult() {
    let check = *CHECKER.read().unwrap_or_else(|e| e.into_inner());
    if let Some(check) = check {
        if check() {
            std::panic::panic_any(Interrupted);
        }
    }
}

/// Run `f`, converting a cancellation into `Err(Interrupted)`.
///
/// Any other panic propagates unchanged, which is the whole distinction this
/// module rests on: a violated contract is still a crash, and only the panic
/// this module raised is an ordinary outcome.
///
/// Nested scopes are sound but redundant — the inner one reports first and the
/// outer sees an ordinary `Err`, so a computation is cancelled once rather
/// than at every level that happened to wrap it.
pub fn catch_interrupt<T>(f: impl FnOnce() -> T) -> Result<T, Interrupted> {
    match catch_unwind(AssertUnwindSafe(f)) {
        Ok(v) => Ok(v),
        Err(payload) => {
            if payload.downcast_ref::<Interrupted>().is_some() {
                Err(Interrupted)
            } else {
                std::panic::resume_unwind(payload)
            }
        }
    }
}

/// Whether the panic now unwinding is a cancellation, given its payload.
///
/// A panic hook is handed the payload and cannot catch anything, so it cannot
/// use [`catch_interrupt`] to tell the two apart. An embedder that silences
/// the default hook's report for cancellations — which are expected, and whose
/// report would otherwise reach the user's terminal on every Ctrl-C — asks
/// here and delegates everything else to the hook it replaced.
pub fn is_interrupt(payload: &(dyn std::any::Any + Send)) -> bool {
    payload.downcast_ref::<Interrupted>().is_some()
}
