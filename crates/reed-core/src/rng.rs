//! Deterministic, splittable pseudo-random numbers for resampling.
//!
//! Every randomised procedure in reed (cross-validation, bootstrap, permutation
//! tests) must satisfy the contract in `docs/spec/architecture.md` §7.6:
//! *same version, same input, same spec ⇒ byte-identical output regardless of
//! thread count*. A shared mutable generator cannot provide that — the stream
//! a worker sees would depend on scheduling. So the generator here is
//! **splittable**: a run derives one independent sub-stream per repetition,
//! and repetition `i` draws the same numbers on 1 thread as on 64.
//!
//! ```
//! use reed_core::rng::Rng;
//! let base = Rng::from_seed(42);
//! let mut a = base.derive(7);   // repetition 7
//! let mut b = base.derive(7);   // same stream, independent state
//! assert_eq!(a.next_u64(), b.next_u64());
//! ```
//!
//! # Algorithm
//!
//! PCG32 (O'Neill 2014, `pcg32_random_r`), a 64-bit LCG with a 32-bit output
//! permutation. It is small, has no dependency, passes the usual statistical
//! batteries, and — unlike a bare xorshift — its output is not linear in the
//! state, which matters when the same generator drives thousands of
//! permutations of the same data.
//!
//! Bounded draws use rejection sampling rather than `% n`: a modulo bias of
//! even one part in `2^32` is a systematic error in a permutation test, and
//! these are cheap to avoid.

// The generator works in `u64` but PCG's output word is `u32` and bounded
// draws return `usize`; every cast below is range-exact by construction (see
// the comments at each site), so the truncation lint is noise here.
#![allow(clippy::cast_possible_truncation)]

/// A PCG32 generator with an explicit stream.
///
/// The type is deliberately `Clone` but not `Copy`: copying a generator and
/// then advancing both produces correlated output, which is almost always a
/// bug. Clone it only to snapshot a known state.
#[derive(Clone, Debug)]
pub struct Rng {
    state: u64,
    /// The odd increment that selects the stream (PCG's `inc`).
    inc: u64,
}

impl Rng {
    /// The stream id [`Rng::from_seed`] uses.
    ///
    /// Any value works; this one matches the reference PCG implementation's
    /// increment. `from_seed(seed)` and `from_seed_stream(seed, DEFAULT_STREAM)`
    /// are the same generator — there is no separate "no stream" case.
    pub const DEFAULT_STREAM: u64 = 1_442_695_040_888_963_407;

    /// Creates a generator from a seed on the default stream.
    #[must_use]
    pub fn from_seed(seed: u64) -> Self {
        Self::from_seed_stream(seed, Self::DEFAULT_STREAM)
    }

    /// Creates a generator from a seed and an explicit stream id.
    ///
    /// Stream ids are the axis along which a run is split across repetitions.
    /// Note that stream `0` is a normal stream, *not* the default one: use
    /// [`Rng::from_seed`] for the default, or pass
    /// [`Rng::DEFAULT_STREAM`] explicitly.
    #[must_use]
    pub fn from_seed_stream(seed: u64, stream: u64) -> Self {
        let mut r = Self {
            state: 0,
            // PCG requires an odd increment; shifting left makes room for it.
            inc: (stream << 1) | 1,
        };
        // PCG's required initialisation sequence: advance once, mix in the
        // seed, advance again.
        r.next_u32();
        r.state = r.state.wrapping_add(seed);
        r.next_u32();
        r
    }

    /// Derives an independent sub-stream for repetition `index`.
    ///
    /// This is what makes parallel resampling reproducible: the caller passes
    /// the repetition number, not a worker id, so the draw does not depend on
    /// which thread happens to run it.
    #[must_use]
    pub fn derive(&self, index: u64) -> Self {
        // Mix the parent state with the index so that two different parents
        // never hand out the same sub-stream, and neighbouring indices do not
        // produce neighbouring states.
        let mixed = splitmix64(self.state ^ index.wrapping_mul(0x9E37_79B9_7F4A_7C15));
        Self::from_seed_stream(mixed, index)
    }

    /// The next 32 random bits.
    pub fn next_u32(&mut self) -> u32 {
        let old = self.state;
        self.state = old
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(self.inc);
        // The shift leaves at most 32 significant bits.
        let xorshifted = (((old >> 18) ^ old) >> 27) as u32;
        let rot = (old >> 59) as u32;
        xorshifted.rotate_right(rot)
    }

    /// The next 64 random bits (two 32-bit draws, high word first).
    pub fn next_u64(&mut self) -> u64 {
        let hi = u64::from(self.next_u32());
        let lo = u64::from(self.next_u32());
        (hi << 32) | lo
    }

    /// A uniform draw in `[0, 1)`.
    ///
    /// Uses 53 bits, the most an `f64` represents exactly, so every value in
    /// the interval is equally likely at that resolution.
    pub fn next_f64(&mut self) -> f64 {
        // `>> 11` keeps the top 53 bits; dividing by 2^53 maps to [0, 1).
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// A uniform integer in `[0, n)`, without modulo bias.
    ///
    /// Rejection sampling: values below `2^64 mod n` are discarded and another
    /// is drawn. The expected number of iterations is below 2 for every `n`
    /// (worst case just above `2^63`, where half the draws are rejected), so
    /// the loop always terminates in practice.
    ///
    /// # Panics
    ///
    /// Panics when `n == 0`: there is no value to return, and silently
    /// returning 0 would hide a caller bug.
    pub fn below(&mut self, n: usize) -> usize {
        assert!(n > 0, "invariant: below() needs a non-empty range");
        let n = n as u64;
        // Reject the top `2^64 % n` values so the remainder is uniform.
        let threshold = n.wrapping_neg() % n;
        loop {
            let v = self.next_u64();
            if v >= threshold {
                // `v % n < n`, and `n` was a usize, so this is exact.
                return (v % n) as usize;
            }
        }
    }

    /// Shuffles a slice in place (Fisher–Yates).
    pub fn shuffle<T>(&mut self, items: &mut [T]) {
        for i in (1..items.len()).rev() {
            let j = self.below(i + 1);
            items.swap(i, j);
        }
    }

    /// Draws `k` distinct indices from `0..n`, in ascending order.
    ///
    /// This is the sampler `poppr::genotype_curve` uses to pick loci. The cost
    /// is `O(min(k, n-k))` draws and memory, not `O(n)`: picking 10 loci out
    /// of a million should not consume a million random numbers per
    /// repetition. The small-`k` side uses Floyd's algorithm; the small
    /// complement side draws the `n-k` excluded indices and inverts.
    ///
    /// # Panics
    ///
    /// Panics when `k > n`: the request is unsatisfiable.
    pub fn sample_n(&mut self, n: usize, k: usize) -> Vec<usize> {
        assert!(k <= n, "invariant: cannot draw {k} distinct of {n}");
        if k == 0 {
            return Vec::new();
        }
        if k == n {
            return (0..n).collect();
        }
        // Drawing the complement is cheaper when k is more than half of n.
        if k > n - k {
            let mut excluded = vec![false; n];
            for i in self.sample_n(n, n - k) {
                excluded[i] = true;
            }
            return (0..n).filter(|&i| !excluded[i]).collect();
        }

        // Floyd's algorithm (Knuth, TAOCP vol. 2, §3.4.2, exercise 15):
        // for each of the last k positions j, draw t from 0..=j; if t is
        // already taken, take j instead. Every k-subset is equally likely and
        // the loop runs exactly k times.
        let mut picked = Vec::with_capacity(k);
        let mut seen = std::collections::HashSet::with_capacity(k);
        for j in (n - k)..n {
            let t = self.below(j + 1);
            let chosen = if seen.contains(&t) { j } else { t };
            seen.insert(chosen);
            picked.push(chosen);
        }
        picked.sort_unstable();
        picked
    }

    /// Draws `k` indices from `0..n` **with replacement**.
    pub fn resample_n(&mut self, n: usize, k: usize) -> Vec<usize> {
        (0..k).map(|_| self.below(n)).collect()
    }

    /// Draws a multinomial count vector: `k` independent draws into `weights.len()`
    /// categories with the given (unnormalised) weights.
    ///
    /// This is the sampler `poppr`'s parametric-bootstrap shuffles use.
    ///
    /// # Panics
    ///
    /// Panics when `weights` is empty, when any weight is negative or `NaN`,
    /// or when the total is not a positive finite number.
    pub fn multinomial(&mut self, k: usize, weights: &[f64]) -> Vec<usize> {
        assert!(
            !weights.is_empty(),
            "invariant: multinomial needs categories"
        );
        // A negative weight makes the cumulative scan run backwards: the
        // category is never reached and the others absorb its share. Measured
        // on [5, -1, 5] the split came out 55.6/44.4 instead of 50/50, with no
        // error. Reject rather than clamp, since a negative weight means the
        // caller computed it wrong.
        assert!(
            weights.iter().all(|&w| w >= 0.0),
            "invariant: multinomial weights must be non-negative"
        );
        let total: f64 = weights.iter().sum();
        assert!(
            total.is_finite() && total > 0.0,
            "invariant: multinomial weights must sum to a positive finite number"
        );
        let mut counts = vec![0usize; weights.len()];
        for _ in 0..k {
            let target = self.next_f64() * total;
            // Cumulative search. Category counts are small (alleles per
            // locus), so a linear scan beats building a table.
            let mut acc = 0.0;
            let mut chosen = weights.len() - 1;
            for (i, &w) in weights.iter().enumerate() {
                acc += w;
                if target < acc {
                    chosen = i;
                    break;
                }
            }
            counts[chosen] += 1;
        }
        counts
    }
}

/// SplitMix64's finaliser: a strong 64-bit mixer used to decorrelate stream
/// ids from the parent state.
fn splitmix64(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn same_seed_same_sequence() {
        let mut a = Rng::from_seed(7);
        let mut b = Rng::from_seed(7);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn different_seeds_differ() {
        let mut a = Rng::from_seed(1);
        let mut b = Rng::from_seed(2);
        let da: Vec<u64> = (0..8).map(|_| a.next_u64()).collect();
        let db: Vec<u64> = (0..8).map(|_| b.next_u64()).collect();
        assert_ne!(da, db);
    }

    #[test]
    fn derive_is_reproducible() {
        let base = Rng::from_seed(42);
        let mut a = base.derive(7);
        let mut b = base.derive(7);
        for _ in 0..50 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    /// The property the whole design exists for: repetition `i` must draw the
    /// same numbers no matter what order the repetitions are visited in, so
    /// the result cannot depend on thread scheduling.
    #[test]
    fn derived_streams_are_independent_of_visit_order() {
        let base = Rng::from_seed(99);
        let forward: Vec<Vec<u64>> = (0..16u64)
            .map(|i| {
                let mut r = base.derive(i);
                (0..4).map(|_| r.next_u64()).collect()
            })
            .collect();
        // Visit in reverse, as a different thread interleaving would.
        let mut backward: Vec<Vec<u64>> = vec![Vec::new(); 16];
        for i in (0..16usize).rev() {
            let mut r = base.derive(i as u64);
            backward[i] = (0..4).map(|_| r.next_u64()).collect();
        }
        assert_eq!(forward, backward);
    }

    #[test]
    fn distinct_stream_ids_give_distinct_sequences() {
        let base = Rng::from_seed(5);
        let mut seen = std::collections::HashSet::new();
        for i in 0..64u64 {
            let mut r = base.derive(i);
            assert!(seen.insert(r.next_u64()), "stream {i} repeated a value");
        }
    }

    #[test]
    fn f64_stays_in_unit_interval() {
        let mut r = Rng::from_seed(11);
        for _ in 0..10_000 {
            let v = r.next_f64();
            assert!((0.0..1.0).contains(&v), "{v}");
        }
    }

    #[test]
    fn below_stays_in_range_and_is_roughly_uniform() {
        let mut r = Rng::from_seed(3);
        let n = 7;
        let mut counts = vec![0usize; n];
        let draws = 70_000;
        for _ in 0..draws {
            let v = r.below(n);
            assert!(v < n);
            counts[v] += 1;
        }
        // Each bucket should hold ~10 000; allow 5% slack.
        let expected = draws / n;
        for (i, &c) in counts.iter().enumerate() {
            let dev = (c as f64 - expected as f64).abs() / expected as f64;
            assert!(dev < 0.05, "bucket {i}: {c} vs expected {expected}");
        }
    }

    #[test]
    fn below_one_always_returns_zero() {
        let mut r = Rng::from_seed(1);
        for _ in 0..100 {
            assert_eq!(r.below(1), 0);
        }
    }

    #[test]
    #[should_panic(expected = "non-empty range")]
    fn below_zero_panics() {
        let mut r = Rng::from_seed(1);
        let _ = r.below(0);
    }

    #[test]
    fn shuffle_is_a_permutation() {
        let mut r = Rng::from_seed(13);
        let mut v: Vec<usize> = (0..100).collect();
        r.shuffle(&mut v);
        let mut sorted = v.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..100).collect::<Vec<_>>());
        assert_ne!(v, (0..100).collect::<Vec<_>>(), "shuffle did nothing");
    }

    #[test]
    fn shuffle_is_reproducible() {
        let mut a = Rng::from_seed(21);
        let mut b = Rng::from_seed(21);
        let mut va: Vec<usize> = (0..50).collect();
        let mut vb: Vec<usize> = (0..50).collect();
        a.shuffle(&mut va);
        b.shuffle(&mut vb);
        assert_eq!(va, vb);
    }

    #[test]
    fn shuffle_handles_empty_and_singleton() {
        let mut r = Rng::from_seed(1);
        let mut empty: Vec<usize> = Vec::new();
        r.shuffle(&mut empty);
        let mut one = vec![42];
        r.shuffle(&mut one);
        assert_eq!(one, vec![42]);
    }

    #[test]
    fn sample_n_returns_distinct_sorted_indices() {
        let mut r = Rng::from_seed(17);
        let picked = r.sample_n(100, 10);
        assert_eq!(picked.len(), 10);
        let unique: std::collections::HashSet<_> = picked.iter().collect();
        assert_eq!(unique.len(), 10, "duplicate index in {picked:?}");
        assert!(picked.iter().all(|&i| i < 100));
    }

    #[test]
    fn sample_n_is_uniform_over_positions() {
        // Every index must appear at the same rate; a biased Algorithm S shows
        // up as the last indices being over- or under-drawn.
        let mut r = Rng::from_seed(23);
        let n = 10;
        let k = 3;
        let trials = 30_000;
        let mut counts = vec![0usize; n];
        for _ in 0..trials {
            for i in r.sample_n(n, k) {
                counts[i] += 1;
            }
        }
        let expected = trials * k / n;
        for (i, &c) in counts.iter().enumerate() {
            let dev = (c as f64 - expected as f64).abs() / expected as f64;
            assert!(dev < 0.05, "index {i}: {c} vs expected {expected}");
        }
    }

    #[test]
    fn sample_n_edge_cases() {
        let mut r = Rng::from_seed(1);
        assert_eq!(r.sample_n(5, 0).len(), 0);
        assert_eq!(r.sample_n(5, 5), vec![0, 1, 2, 3, 4]);
    }

    #[test]
    #[should_panic(expected = "cannot draw")]
    fn sample_n_more_than_n_panics() {
        let mut r = Rng::from_seed(1);
        let _ = r.sample_n(3, 4);
    }

    #[test]
    fn resample_n_allows_repeats() {
        let mut r = Rng::from_seed(29);
        let picked = r.resample_n(5, 40);
        assert_eq!(picked.len(), 40);
        assert!(picked.iter().all(|&i| i < 5));
        let unique: std::collections::HashSet<_> = picked.iter().collect();
        assert!(unique.len() < 40, "40 draws from 5 with replacement should repeat");
    }

    #[test]
    fn multinomial_sums_to_the_draw_count() {
        let mut r = Rng::from_seed(31);
        let counts = r.multinomial(100, &[1.0, 2.0, 3.0]);
        assert_eq!(counts.iter().sum::<usize>(), 100);
        assert_eq!(counts.len(), 3);
    }

    #[test]
    fn multinomial_follows_the_weights() {
        let mut r = Rng::from_seed(37);
        let weights = [0.1, 0.2, 0.7];
        let draws = 50_000;
        let counts = r.multinomial(draws, &weights);
        for (i, &w) in weights.iter().enumerate() {
            let expected = draws as f64 * w;
            let dev = (counts[i] as f64 - expected).abs() / expected;
            assert!(dev < 0.05, "category {i}: {} vs expected {expected}", counts[i]);
        }
    }

    #[test]
    fn multinomial_handles_a_zero_weight() {
        let mut r = Rng::from_seed(41);
        let counts = r.multinomial(1000, &[0.0, 1.0]);
        assert_eq!(counts[0], 0, "a zero-weight category must never be drawn");
        assert_eq!(counts[1], 1000);
    }

    #[test]
    #[should_panic(expected = "non-negative")]
    fn multinomial_rejects_a_negative_weight() {
        // Before the guard this ran silently and split 55.6/44.4 instead of
        // 50/50, because the cumulative scan ran backwards past the negative
        // category.
        let mut r = Rng::from_seed(3);
        let _ = r.multinomial(1000, &[5.0, -1.0, 5.0]);
    }

    #[test]
    #[should_panic(expected = "positive finite")]
    fn multinomial_rejects_zero_total() {
        let mut r = Rng::from_seed(1);
        let _ = r.multinomial(10, &[0.0, 0.0]);
    }

    /// Pins the generator to O'Neill's PCG32, cross-checked against an
    /// independent implementation of the published algorithm. Without this,
    /// a refactor could silently change every number in every downstream
    /// resampling result while all the statistical tests still pass.
    #[test]
    fn matches_the_pcg32_reference_sequence() {
        let mut r = Rng::from_seed(42);
        let got: Vec<u32> = (0..5).map(|_| r.next_u32()).collect();
        assert_eq!(
            got,
            vec![492_690_617, 1_919_685_028, 3_561_993_920, 683_038_915, 1_183_706_632]
        );
    }

    /// Streams derived from the same parent must be statistically
    /// independent. Without this, a change to `derive`'s mixer could make
    /// neighbouring repetitions correlated, and every existing test would
    /// still pass — the damage would only show up as inflated false-positive
    /// rates in the permutation tests that consume this.
    #[test]
    fn derived_streams_are_uncorrelated() {
        let base = Rng::from_seed(123);
        let n = 20_000;
        let first: Vec<f64> = (0..n).map(|i| base.derive(i as u64).next_f64()).collect();
        let second: Vec<f64> = (0..n)
            .map(|i| {
                let mut r = base.derive(i as u64 + 1);
                r.next_f64()
            })
            .collect();
        let r = correlation(&first, &second);
        assert!(r.abs() < 0.02, "neighbouring streams correlate: r = {r:.4}");
    }

    /// Positions within a stream must not correlate with each other either.
    #[test]
    fn positions_within_derived_streams_are_uncorrelated() {
        let base = Rng::from_seed(321);
        let n = 20_000;
        let p0: Vec<f64> = (0..n).map(|i| base.derive(i as u64).next_f64()).collect();
        let p1: Vec<f64> = (0..n)
            .map(|i| {
                let mut r = base.derive(i as u64);
                r.next_f64();
                r.next_f64()
            })
            .collect();
        let r = correlation(&p0, &p1);
        assert!(r.abs() < 0.02, "positions correlate: r = {r:.4}");
    }

    /// `sample_n` must not walk every position when `k` is small.
    ///
    /// The test uses a tiny `n` on purpose: with Floyd's algorithm the draw
    /// count is exactly `k` (plus rejection retries), so a 5-of-20 call must
    /// leave the generator where 5 draws do. Under the old Algorithm S it
    /// walked all 20 positions. Wall-clock assertions are flaky on CI; this
    /// one is exact.
    #[test]
    fn sample_n_costs_min_k_not_n() {
        let n = 20;
        let k = 5;
        let mut used = Rng::from_seed(77);
        let _ = used.sample_n(n, k);

        // Replay k draws on a fresh generator. `below(20)` almost never
        // rejects, so the two states should line up after exactly k `u64`
        // draws (= 2k `u32` draws).
        let mut reference = Rng::from_seed(77);
        for _ in 0..k {
            let _ = reference.next_u64();
        }
        assert_eq!(
            used.next_u32(),
            reference.next_u32(),
            "sample_n({n}, {k}) did not consume exactly {k} draws; it walked the range"
        );
    }

    /// Every `k` must produce exactly `k` distinct indices in range, ascending,
    /// including both complement paths.
    #[test]
    fn sample_n_is_correct_for_every_k() {
        let mut r = Rng::from_seed(79);
        for n in [1usize, 2, 5, 50] {
            for k in 0..=n {
                let picked = r.sample_n(n, k);
                assert_eq!(picked.len(), k, "n={n} k={k}");
                let unique: std::collections::HashSet<_> = picked.iter().collect();
                assert_eq!(unique.len(), k, "n={n} k={k} has duplicates");
                assert!(picked.iter().all(|&i| i < n), "n={n} k={k} out of range");
                assert!(
                    picked.windows(2).all(|w| w[0] < w[1]),
                    "n={n} k={k} not ascending"
                );
            }
        }
    }

    /// Pearson correlation, used only by the tests above.
    fn correlation(a: &[f64], b: &[f64]) -> f64 {
        let n = a.len() as f64;
        let ma = a.iter().sum::<f64>() / n;
        let mb = b.iter().sum::<f64>() / n;
        let cov: f64 = a.iter().zip(b).map(|(x, y)| (x - ma) * (y - mb)).sum::<f64>() / n;
        let sa = (a.iter().map(|x| (x - ma).powi(2)).sum::<f64>() / n).sqrt();
        let sb = (b.iter().map(|x| (x - mb).powi(2)).sum::<f64>() / n).sqrt();
        cov / (sa * sb)
    }

    #[test]
    fn next_u64_uses_both_halves() {
        // A generator that only advanced one 32-bit half would show a constant
        // high or low word.
        let mut r = Rng::from_seed(43);
        let vals: Vec<u64> = (0..64).map(|_| r.next_u64()).collect();
        let hi_unique: std::collections::HashSet<_> = vals.iter().map(|v| v >> 32).collect();
        let lo_unique: std::collections::HashSet<_> = vals.iter().map(|v| v & 0xFFFF_FFFF).collect();
        assert!(hi_unique.len() > 32, "high word barely varies");
        assert!(lo_unique.len() > 32, "low word barely varies");
    }
}
