//! The determinism contract for randomised procedures
//! (`docs/spec/architecture.md` §7.6).
//!
//! A resampling run must produce byte-identical output regardless of thread
//! count. The mechanism is `Rng::derive(repetition)`: each repetition gets
//! its own sub-stream keyed by its index, so the numbers a repetition sees do
//! not depend on which thread picked it up or in what order.

use reed_core::rng::Rng;

/// The shape every resampling driver uses: derive a sub-stream per
/// repetition, do the work, collect in repetition order.
fn per_repetition(seed: u64, reps: usize, work: impl Fn(&mut Rng) -> u64 + Sync) -> Vec<u64> {
    let base = Rng::from_seed(seed);
    (0..reps).map(|i| work(&mut base.derive(i as u64))).collect()
}

#[test]
fn serial_and_parallel_runs_agree() {
    let seed = 2024;
    let reps = 256;

    // What a single-threaded driver does.
    let serial = per_repetition(seed, reps, |r| {
        let mut v: Vec<u32> = (0..32).collect();
        r.shuffle(&mut v);
        v.iter().fold(0u64, |acc, &x| acc.wrapping_mul(31).wrapping_add(u64::from(x)))
    });

    // What a parallel driver does: workers claim indices out of order and
    // write into their slot. `reed-core` stays dependency-free (§14.6), so
    // this uses `std::thread` rather than rayon.
    let base = Rng::from_seed(seed);
    let mut slots: Vec<Option<u64>> = vec![None; reps];
    let workers = 4;
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..workers)
            .map(|w| {
                let base = base.clone();
                scope.spawn(move || {
                    // Stride by worker count: each worker visits a different
                    // subset, and the interleaving differs from the serial run.
                    let mut mine = Vec::new();
                    let mut i = w;
                    while i < reps {
                        let mut r = base.derive(i as u64);
                        let mut v: Vec<u32> = (0..32).collect();
                        r.shuffle(&mut v);
                        let acc = v
                            .iter()
                            .fold(0u64, |acc, &x| acc.wrapping_mul(31).wrapping_add(u64::from(x)));
                        mine.push((i, acc));
                        i += workers;
                    }
                    mine
                })
            })
            .collect();
        for h in handles {
            for (i, v) in h.join().expect("worker panicked") {
                slots[i] = Some(v);
            }
        }
    });
    let parallel: Vec<u64> = slots.into_iter().map(|v| v.expect("slot filled")).collect();

    assert_eq!(serial, parallel, "thread count changed the resampling output");
}

#[test]
fn a_run_is_reproducible_across_invocations() {
    let run = || per_repetition(7, 64, |r| r.next_u64());
    assert_eq!(run(), run());
}

#[test]
fn different_seeds_give_different_runs() {
    let a = per_repetition(1, 32, |r| r.next_u64());
    let b = per_repetition(2, 32, |r| r.next_u64());
    assert_ne!(a, b);
}

/// Repetition `i` must be identical whether it is reached first or last. This
/// is the property that makes the parallel result order-independent rather
/// than merely order-collected.
#[test]
fn a_repetition_does_not_depend_on_what_ran_before_it() {
    let base = Rng::from_seed(11);
    let first = {
        let mut r = base.derive(5);
        (0..8).map(|_| r.next_u64()).collect::<Vec<_>>()
    };
    // Burn a different number of draws in between.
    let mut noise = base.derive(3);
    for _ in 0..1000 {
        noise.next_u64();
    }
    let later = {
        let mut r = base.derive(5);
        (0..8).map(|_| r.next_u64()).collect::<Vec<_>>()
    };
    assert_eq!(first, later);
}

/// Two repetitions of the same run must not draw the same stream: a bootstrap
/// where every replicate used identical indices would have zero variance.
#[test]
fn repetitions_get_distinct_streams() {
    let base = Rng::from_seed(3);
    let mut seen = std::collections::HashSet::new();
    for i in 0..512 {
        let mut r = base.derive(i);
        assert!(seen.insert(r.next_u64()), "repetition {i} repeated a stream");
    }
}

/// A bootstrap replicate must actually vary: `resample_n` with replacement
/// over a small set should produce different index multisets across
/// repetitions.
#[test]
fn bootstrap_replicates_vary() {
    let base = Rng::from_seed(19);
    let draws: Vec<Vec<usize>> = (0..16)
        .map(|i| base.derive(i).resample_n(10, 20))
        .collect();
    let unique: std::collections::HashSet<_> = draws.iter().collect();
    assert_eq!(unique.len(), draws.len(), "two replicates drew identical samples");
}
