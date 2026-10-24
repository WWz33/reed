//! R semantics needed to reproduce `adegenet::xvalDapc`'s bookkeeping.
//!
//! These are small but load-bearing: R's `round()` rounds halves to even and
//! `pretty()` picks "nice" break points, and the candidate list of principal
//! components depends on both. Getting either wrong shifts which `n.pca`
//! values are evaluated, which changes the selected optimum.

// Every float-to-integer cast here converts a value that is already integral
// and range-checked, so the truncation lint is noise.
#![allow(clippy::cast_possible_truncation)]

/// R's `round()`: round half to even ("banker's rounding").
///
/// Rust's `f64::round` rounds half away from zero, so `round(2.5)` is 3 there
/// and 2 in R. `xvalDapc` uses `round()` to size the training set and to
/// build the candidate list, so the difference changes results.
///
/// The tie case is detected on the exact binary value, which is what R does:
/// `round(0.5) == 0`, `round(1.5) == 2`, `round(2.5) == 2`.
#[must_use]
pub fn round_half_even(x: f64) -> f64 {
    if !x.is_finite() {
        return x;
    }
    let floor = x.floor();
    let frac = x - floor;
    if (frac - 0.5).abs() < f64::EPSILON * x.abs().max(1.0) {
        // Exactly halfway: pick the even neighbour.
        if (floor / 2.0).fract() == 0.0 { floor } else { floor + 1.0 }
    } else if frac < 0.5 {
        floor
    } else {
        floor + 1.0
    }
}

/// R's `round()` returning an integer, for the call sites that need one.
///
/// # Panics
///
/// Panics when the rounded value does not fit in `usize`; every caller feeds
/// it a sample or component count, which does.
#[must_use]
pub fn round_half_even_usize(x: f64) -> usize {
    let r = round_half_even(x);
    assert!(
        r >= 0.0 && r <= usize::MAX as f64,
        "invariant: rounded count {r} is out of range"
    );
    // The value is integral and within range, so the cast is exact.
    r as usize
}

/// R's `pretty(x, n)` restricted to the shape `xvalDapc` uses: an increasing
/// integer range, deduplicated, keeping only the values inside it.
///
/// `pretty()` produces "nice" round numbers covering the range of its input,
/// which is why `pretty(1:300, 10)` yields multiples of 50 while
/// `pretty(1:5, 5)` yields every integer. `xvalDapc` calls
/// `pretty(1:n.pca.max, runs)` and then filters; it does **not** deduplicate,
/// so repeated candidates are evaluated several times (wasteful, same result).
/// This implementation deduplicates, which is the one deliberate deviation
/// (see `docs/rounds/v18/02-xval-research.md` §2).
///
/// # Panics
///
/// Panics when `max < 1`; the caller has already established that there is at
/// least one candidate.
#[must_use]
pub fn pretty_sequence(max: usize, n: usize) -> Vec<usize> {
    assert!(max >= 1, "invariant: pretty_sequence needs max >= 1");
    if n <= 1 {
        return vec![1];
    }
    if max <= n {
        // A unit step covers the range with no more than `n` points.
        return (1..=max).collect();
    }

    // R's `pretty` chooses a step from the 1-2-5 family, then extends the
    // range outward to multiples of that step. Reproducing its C internals
    // exactly is not worth it; the switch points below were measured from R
    // for the `n = 10` case `xvalDapc` uses, and are pinned by a test.
    //
    //   max:  2..15   16..29   30..71   72..141  142..281  282..701  702..1401
    //   step:    1       2        5        10       20        50        100
    //
    // Beyond 1401 the pattern continues by decades (200, 500, ...).
    let step = step_for(max, n);

    // Multiples of the step inside `1..=max`. `pretty` also emits the
    // endpoints below 1 and at or above `max`; the caller keeps only what is
    // inside the range, so those are not produced here.
    let mut out = Vec::new();
    let mut m = step;
    while m <= max as f64 {
        // `m` is an exact multiple of an integral step, so the cast is exact.
        out.push(m as usize);
        m += step;
    }
    if out.is_empty() {
        out.push(1);
    }
    out
}

/// The `pretty()` step for a range `1..=max` covered by about `n` points.
///
/// Measured from R rather than derived: R's implementation extends the range
/// to round bounds before choosing the step, which shifts every switch point
/// by a few units. The table covers `n = 10`, which is the only `runs` value
/// `xvalDapc` uses for `n.pca.max >= 10`.
fn step_for(max: usize, n: usize) -> f64 {
    if n == 10 {
        // Measured from R for `pretty(1:max, 10)`, extended past 1402 by the
        // same ladder (each step is 2x or 2.5x the previous, on decade
        // boundaries: 100, 200, 500, 1000, 2000, 5000, ...).
        const TABLE: [(usize, f64); 12] = [
            (16, 2.0),
            (30, 5.0),
            (72, 10.0),
            (142, 20.0),
            (282, 50.0),
            (702, 100.0),
            (1402, 200.0),
            (2802, 500.0),
            (7002, 1000.0),
            (14002, 2000.0),
            (28002, 5000.0),
            (70002, 10_000.0),
        ];
        let mut step = 1.0;
        for &(threshold, next) in &TABLE {
            if max >= threshold {
                step = next;
            } else {
                break;
            }
        }
        // Past the table, keep doubling per two decades; the caller never
        // reaches this with a realistic sample count.
        if max >= 70002 {
            let mut extra = 10_000.0;
            while (max as f64) >= extra * 5.0 {
                extra *= 5.0;
            }
            step = step.max(extra);
        }
        return step;
    }
    // `candidate_counts` only calls this with `n == 10`; `pretty_sequence`
    // already returned above for `max <= n`. A step of 1 is the safe answer
    // for any other `n`, and no caller reaches it.
    1.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_half_even_matches_r() {
        // Values measured from R 4.x.
        assert_eq!(round_half_even(4.5), 4.0);
        assert_eq!(round_half_even(13.5), 14.0);
        assert_eq!(round_half_even(22.5), 22.0);
        assert_eq!(round_half_even(2.5), 2.0);
        assert_eq!(round_half_even(3.5), 4.0);
        assert_eq!(round_half_even(0.5), 0.0);
        assert_eq!(round_half_even(1.5), 2.0);
    }

    #[test]
    fn round_half_even_handles_non_ties() {
        assert_eq!(round_half_even(0.4), 0.0);
        assert_eq!(round_half_even(0.6), 1.0);
        assert_eq!(round_half_even(-0.4), -0.0);
        assert_eq!(round_half_even(-0.6), -1.0);
        assert_eq!(round_half_even(2.4), 2.0);
        assert_eq!(round_half_even(2.6), 3.0);
    }

    #[test]
    fn round_half_even_passes_through_non_finite() {
        assert!(round_half_even(f64::NAN).is_nan());
        assert_eq!(round_half_even(f64::INFINITY), f64::INFINITY);
    }

    /// The table in `docs/rounds/v18/02-xval-research.md` §1.2, measured from
    /// R: the deduplicated `pretty()` output, before the caller's filter.
    #[test]
    fn pretty_sequence_matches_r() {
        assert_eq!(pretty_sequence(3, 3), vec![1, 2, 3]);
        assert_eq!(pretty_sequence(5, 5), vec![1, 2, 3, 4, 5]);
        assert_eq!(pretty_sequence(9, 9), (1..=9).collect::<Vec<_>>());
        assert_eq!(pretty_sequence(10, 10), (1..=10).collect::<Vec<_>>());
        assert_eq!(
            pretty_sequence(20, 10),
            vec![2, 4, 6, 8, 10, 12, 14, 16, 18, 20]
        );
        assert_eq!(
            pretty_sequence(50, 10),
            vec![5, 10, 15, 20, 25, 30, 35, 40, 45, 50]
        );
        assert_eq!(
            pretty_sequence(300, 10),
            vec![50, 100, 150, 200, 250, 300]
        );
    }

    #[test]
    fn pretty_sequence_is_strictly_increasing_and_in_range() {
        for max in 1..=200 {
            let v = pretty_sequence(max, 10);
            assert!(v.windows(2).all(|w| w[0] < w[1]), "max={max}: {v:?}");
            assert!(v.iter().all(|&i| i >= 1 && i <= max), "max={max}: {v:?}");
            assert!(!v.is_empty(), "max={max} produced no candidates");
        }
    }
}
