//! Cancellation reaches a real computation, and the library survives it.
//!
//! ```text
//!   cargo test --test interrupt
//! ```
//!
//! Two things are asserted here, and the second is why this file is large.
//! The first is the mechanism: [`poll`](symfn::interrupt::poll) driven
//! directly, with nothing else in the way. The second is that the poll sites
//! are actually *on* the paths that run long, and that unwinding through the
//! crate's caches leaves it usable — properties of the whole tree rather than
//! of one module, and both the way this feature would quietly stop working. A
//! poll site is easy to lose in a refactor, and a memo that kept a half-built
//! entry would be wrong only on the call *after* the one that was cancelled.
//!
//! **The mechanism's tests are here rather than beside it in
//! `src/interrupt.rs` because the checker is process-wide.** `cargo test`
//! threads a module's tests alongside every other test in the same binary, so
//! a checker installed by a unit test cancelled unrelated work in `hopf` and
//! `hl`: the interrupt tests passed and three others failed. A separate
//! integration binary is a separate process, which is the isolation a
//! process-wide switch actually needs. Within this binary the tests still hold
//! a mutex and clear the checker on the way out, for the same reason.

use std::sync::atomic::{AtomicU32, Ordering};
use std::sync::Mutex;

use symfn::interrupt::{self, Interrupted};
use symfn::partition::Partition;
use symfn::sym::{Schur, SymFn};
use symfn::{Rational, Ring};

static SERIAL: Mutex<()> = Mutex::new(());

/// Polls seen since the last [`arm`]. The checker is a plain function pointer
/// and carries no state, so its state lives here.
static SEEN: AtomicU32 = AtomicU32::new(0);

/// Cancel once the computation has made some progress, rather than at the
/// first poll: a checker that fires immediately would pass even if the only
/// poll site reached were the one in the entry point itself.
///
/// Counted in *checker calls*, which `poll` throttles to one per stride — so
/// this is a few hundred poll sites crossed, not four.
const CANCEL_AFTER: u32 = 4;

fn cancel_after_progress() -> bool {
    SEEN.fetch_add(1, Ordering::Relaxed) >= CANCEL_AFTER
}

fn never() -> bool {
    false
}

/// How many checker calls [`cancel_at_depth`] lets pass, set per pass by
/// `answers_survive_cancellation_at_many_depths`.
static DEPTH: AtomicU32 = AtomicU32::new(1);

fn cancel_at_depth() -> bool {
    SEEN.fetch_add(1, Ordering::Relaxed) >= DEPTH.load(Ordering::Relaxed)
}

/// Install the checker and silence the default panic report, returning the
/// hook to put back. A cancellation is a panic, and an unsilenced one would
/// print a Rust backtrace notice into the test output on every run.
fn arm(check: fn() -> bool) -> Box<dyn Fn(&std::panic::PanicHookInfo<'_>) + Sync + Send + 'static> {
    SEEN.store(0, Ordering::Relaxed);
    let previous = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    interrupt::set_checker(check);
    previous
}

fn disarm(hook: Box<dyn Fn(&std::panic::PanicHookInfo<'_>) + Sync + Send + 'static>) {
    interrupt::clear_checker();
    std::panic::set_hook(hook);
}

fn s(v: &[u32]) -> Schur<Rational> {
    Schur::monomial(Partition::new(v.iter().copied()), Rational::one())
}

/// The case the mechanism was built for: the plethysm whose p → s tail is the
/// crate's longest reachable call (`docs/record/plethysm.md`).
#[test]
fn a_plethysm_in_flight_is_cancelled() {
    let _g = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let hook = arm(cancel_after_progress);
    let r = interrupt::catch_interrupt(|| symfn::plethysm(&s(&[4]), &s(&[4])));
    disarm(hook);
    assert_eq!(
        r,
        Err(Interrupted),
        "s_4[s_4] ran to completion uncancelled"
    );
}

/// Products are the other path every consumer reaches, and they share none of
/// plethysm's poll sites.
///
/// Many terms rather than one large shape, because the two are cancelled at
/// different sites: a wide element reaches the per-pair poll in `mul_with`,
/// while a single deep shape reaches only the Littlewood–Richardson fill's.
#[test]
fn a_schur_product_in_flight_is_cancelled() {
    let _g = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let mut a: Schur<Rational> = Schur::zero();
    for lambda in symfn::partitions_of(8) {
        a.add_term(lambda, Rational::one());
    }
    let hook = arm(cancel_after_progress);
    let r = interrupt::catch_interrupt(|| a.mul(&a));
    disarm(hook);
    assert_eq!(r, Err(Interrupted), "a Schur product ran uncancelled");
}

/// The answer after a cancellation is the answer, not a poisoned cache's idea
/// of it. `memo` computes outside its guard and clears lock poison rather than
/// propagating it (`src/memo.rs`), so an unwind through a table leaves nothing
/// behind — this is what asserts that stays true.
#[test]
fn a_cancelled_computation_leaves_the_caches_usable() {
    let _g = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let want = symfn::plethysm(&s(&[3]), &s(&[2]));

    let hook = arm(cancel_after_progress);
    let cancelled = interrupt::catch_interrupt(|| symfn::plethysm(&s(&[4]), &s(&[4])));
    disarm(hook);
    assert_eq!(cancelled, Err(Interrupted), "nothing was cancelled");

    let got = symfn::plethysm(&s(&[3]), &s(&[2]));
    assert_eq!(
        got, want,
        "s_3[s_2] changed after an unrelated cancellation"
    );
}

/// A checker that declines must cost only time, never an answer. Both sides
/// run with poll sites live, so this also pins that a poll cannot perturb a
/// result on its own.
#[test]
fn a_declining_checker_changes_no_answer() {
    let _g = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let want = symfn::plethysm(&s(&[3]), &s(&[3]));

    let hook = arm(never);
    let got = interrupt::catch_interrupt(|| symfn::plethysm(&s(&[3]), &s(&[3])));
    disarm(hook);

    assert_eq!(got, Ok(want), "an installed checker changed s_3[s_3]");
}

/// Cancel at many different depths, then demand every answer back.
///
/// The reasoning that says a cancellation cannot corrupt a cache is a claim
/// about every store in the crate: each one writes a value that is already
/// finished, computed outside the guard it is written under, so an unwind
/// leaves a table either untouched or holding correct entries
/// (`src/memo.rs`; `bold_guarded` in `src/character_basis.rs`). That is an
/// argument about code as it stands today, and the thing it protects — a wrong
/// answer on the call *after* a Ctrl-C — is silent. So it also gets an
/// experiment: cancel at a different depth on each pass and check the answers
/// against values taken before any of it.
///
/// The depths are a fixed sequence rather than a random one; what matters is
/// that they differ and that a failure names the pass that produced it.
#[test]
fn answers_survive_cancellation_at_many_depths() {
    let _g = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    interrupt::clear_checker();

    let cases: Vec<(Schur<Rational>, Schur<Rational>)> = vec![
        (s(&[3]), s(&[2])),
        (s(&[2, 1]), s(&[2])),
        (s(&[2]), s(&[2, 1])),
        (s(&[4]), s(&[2])),
    ];
    let want: Vec<Schur<Rational>> = cases.iter().map(|(f, g)| symfn::plethysm(f, g)).collect();

    for pass in 1..=40u32 {
        // A cheap deterministic spread over the depths a cancellation can land
        // at, coprime stride so successive passes do not repeat quickly.
        DEPTH.store(1 + (pass * 7) % 23, Ordering::Relaxed);
        let hook = arm(cancel_at_depth);
        let cancelled = interrupt::catch_interrupt(|| symfn::plethysm(&s(&[4]), &s(&[4])));
        disarm(hook);
        assert_eq!(cancelled, Err(Interrupted), "pass {pass} was not cancelled");

        for (i, (f, g)) in cases.iter().enumerate() {
            assert_eq!(
                symfn::plethysm(f, g),
                want[i],
                "case {i} changed after a cancellation on pass {pass}"
            );
        }
    }
}

/// With nothing installed the poll sites are inert, which is the state every
/// Rust caller of this crate is in.
#[test]
fn without_a_checker_nothing_is_cancelled() {
    let _g = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    interrupt::clear_checker();
    let r = interrupt::catch_interrupt(|| symfn::plethysm(&s(&[3]), &s(&[3])));
    assert!(r.is_ok(), "a computation was cancelled with no checker set");
}

// --- the mechanism itself, driven through `poll` rather than through the
// mathematics ---------------------------------------------------------------

fn always() -> bool {
    true
}

#[test]
fn poll_without_a_checker_is_a_no_op() {
    let _g = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    interrupt::clear_checker();
    for _ in 0..4 * interrupt::STRIDE {
        interrupt::poll();
    }
}

#[test]
fn a_checker_that_declines_never_cancels() {
    let _g = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let hook = arm(never);
    let r = interrupt::catch_interrupt(|| {
        for _ in 0..4 * interrupt::STRIDE {
            interrupt::poll();
        }
        7
    });
    disarm(hook);
    assert_eq!(r, Ok(7), "a checker returning false cancelled anyway");
}

#[test]
fn a_checker_that_cancels_is_reported_within_one_stride() {
    let _g = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let hook = arm(always);
    let mut polls = 0u32;
    let r = interrupt::catch_interrupt(|| loop {
        interrupt::poll();
        polls += 1;
        assert!(polls <= interrupt::STRIDE, "no cancellation after a stride");
    });
    disarm(hook);
    assert_eq!(r, Err(Interrupted));
}

/// The distinction the module rests on: only its own panic is an outcome.
#[test]
fn an_unrelated_panic_is_not_caught() {
    let _g = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    interrupt::clear_checker();
    let previous = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let escaped = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        let _ = interrupt::catch_interrupt(|| panic!("a violated contract"));
    }));
    std::panic::set_hook(previous);
    assert!(escaped.is_err(), "catch_interrupt swallowed a real panic");
}

#[test]
fn clearing_the_checker_stops_cancellation() {
    let _g = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let hook = arm(always);
    interrupt::clear_checker();
    let r = interrupt::catch_interrupt(|| {
        for _ in 0..4 * interrupt::STRIDE {
            interrupt::poll();
        }
        "finished"
    });
    std::panic::set_hook(hook);
    assert_eq!(r, Ok("finished"), "a cleared checker still cancelled");
}

#[test]
fn is_interrupt_separates_the_two_payloads() {
    let _g = SERIAL.lock().unwrap_or_else(|e| e.into_inner());
    let previous = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));
    let cancelled = std::panic::catch_unwind(|| {
        std::panic::panic_any(Interrupted);
    })
    .unwrap_err();
    let contract = std::panic::catch_unwind(|| panic!("a violated contract")).unwrap_err();
    std::panic::set_hook(previous);
    assert!(
        interrupt::is_interrupt(&*cancelled),
        "cancellation read as a contract"
    );
    assert!(
        !interrupt::is_interrupt(&*contract),
        "contract read as a cancellation"
    );
}
