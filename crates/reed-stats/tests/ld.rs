//! LD canonical tests: the integer-plane r² against a naive per-sample
//! f64 Pearson correlation (the reference), hand-computed cases, decay
//! binning, and ordered-carry pruning semantics.

use reed_stats::ld::{LdDecayAcc, LdPruner, pair_sums};

/// Builds a site's four planes from per-sample genotypes
/// (`(slot0, slot1)` with `None` = missing).
fn planes(gts: &[(Option<u8>, Option<u8>)]) -> [&'static [u64]; 4] {
    let n = gts.len().div_ceil(64);
    let mut p: Vec<Vec<u64>> = vec![vec![0; n]; 4];
    for (i, (s0, s1)) in gts.iter().enumerate() {
        let bit = |i: usize| 1_u64 << (i % 64);
        match s0 {
            Some(a) if *a > 0 => p[0][i / 64] |= bit(i),
            Some(_) => {}
            None => p[1][i / 64] |= bit(i),
        }
        match s1 {
            Some(a) if *a > 0 => p[2][i / 64] |= bit(i),
            Some(_) => {}
            None => p[3][i / 64] |= bit(i),
        }
    }
    // Leak into 'static — tests only, bounded by the process lifetime.
    let boxed: Vec<Box<[u64]>> = p.into_iter().map(|v| v.into_boxed_slice()).collect();
    let leaked: Vec<&'static [u64]> = boxed
        .into_iter()
        .map(|b| -> &'static [u64] { Box::leak(b) })
        .collect();
    [leaked[0], leaked[1], leaked[2], leaked[3]]
}

/// Reference: per-sample dosage Pearson r² over the complete-case samples.
fn ref_r2(x: &[(Option<u8>, Option<u8>)], y: &[(Option<u8>, Option<u8>)]) -> Option<f64> {
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for (gx, gy) in x.iter().zip(y) {
        if let (Some(a), Some(b), Some(c), Some(d)) = (gx.0, gx.1, gy.0, gy.1) {
            xs.push(f64::from(a) + f64::from(b));
            ys.push(f64::from(c) + f64::from(d));
        }
    }
    let n = xs.len() as f64;
    if n < 2.0 {
        return None;
    }
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let mut cov = 0.0;
    let mut vx = 0.0;
    let mut vy = 0.0;
    for (a, b) in xs.iter().zip(&ys) {
        cov += (a - mx) * (b - my);
        vx += (a - mx) * (a - mx);
        vy += (b - my) * (b - my);
    }
    if vx <= 0.0 || vy <= 0.0 {
        return None;
    }
    Some(cov * cov / (vx * vy))
}

fn gts(spec: &[(u8, u8)]) -> Vec<(Option<u8>, Option<u8>)> {
    let allele = |c: u8| -> Option<u8> {
        match c {
            b'.' => None,
            b'0' => Some(0),
            _ => Some(1),
        }
    };
    spec.iter().map(|&(a, b)| (allele(a), allele(b))).collect()
}

#[test]
fn integer_r2_matches_reference_on_random_sites() {
    // Deterministic pseudo-random genotypes (xorshift).
    let mut seed = 0x1234_5678_9abc_def0_u64;
    let mut next = || {
        seed ^= seed >> 12;
        seed ^= seed << 25;
        seed ^= seed >> 27;
        seed
    };
    for case in 0..200 {
        let n = 5 + (next() % 60) as usize;
        let mut x_spec = Vec::new();
        let mut y_spec = Vec::new();
        for _ in 0..n {
            let r = next();
            x_spec.push((
                b"0011.."[((r >> 1) % 6) as usize],
                b"0011.."[((r >> 7) % 6) as usize],
            ));
            y_spec.push((
                b"0011.."[((r >> 13) % 6) as usize],
                b"0011.."[((r >> 19) % 6) as usize],
            ));
        }
        let (x, y) = (gts(&x_spec), gts(&y_spec));
        let sums = pair_sums(planes(&x), planes(&y), x.len());
        let mine = sums.r2();
        let expected = ref_r2(&x, &y);
        match (mine, expected) {
            (None, None) => {}
            (Some(a), Some(b)) => {
                // The exact-i128 path is the ground truth; the reference's
                // f64 mean-subtraction can produce spurious ~1e-33 values
                // where the true covariance is exactly zero. Tolerate
                // either side being numerically zero.
                assert!(
                    (a - b).abs() <= 1e-6 * b.abs().max(1e-300).max(a.abs())
                        || a.abs() < 1e-30
                        || b.abs() < 1e-30,
                    "case {case}: {a} vs {b}"
                );
            }
            (a, b) => panic!("case {case}: definedness mismatch {a:?} vs {b:?}"),
        }
        // The covariance also agrees with the reference.
        let (Some(mine_c), Some(ref_c)) = (sums.covariance(), ref_cov(&x, &y)) else {
            continue;
        };
        // Same artifact family as r²: the reference's f64 mean-subtraction
        // can be nonzero where the exact-i128 covariance is true zero.
        assert!(
            (mine_c - ref_c).abs() <= 1e-6 * ref_c.abs().max(1e-300).max(mine_c.abs())
                || mine_c.abs() < 1e-12
                || ref_c.abs() < 1e-12,
            "cov case {case}: {mine_c} vs {ref_c}"
        );
    }
}

fn ref_cov(x: &[(Option<u8>, Option<u8>)], y: &[(Option<u8>, Option<u8>)]) -> Option<f64> {
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for (gx, gy) in x.iter().zip(y) {
        if let (Some(a), Some(b), Some(c), Some(d)) = (gx.0, gx.1, gy.0, gy.1) {
            xs.push(f64::from(a) + f64::from(b));
            ys.push(f64::from(c) + f64::from(d));
        }
    }
    let n = xs.len() as f64;
    if n == 0.0 {
        return None;
    }
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    Some(xs.iter().zip(&ys).map(|(a, b)| (a - mx) * (b - my)).sum::<f64>() / n)
}

#[test]
fn hand_computed_r2_values() {
    // Perfect correlation: x == y → r² = 1.
    let spec = [(b'0', b'0'), (b'0', b'1'), (b'1', b'1'), (b'1', b'0')];
    let (x, y) = (gts(&spec), gts(&spec));
    assert_eq!(pair_sums(planes(&x), planes(&y), x.len()).r2(), Some(1.0));
    // Independent dosage: hand-built uncorrelated vectors.
    let x_spec = [(b'0', b'0'), (b'0', b'1'), (b'1', b'1'), (b'1', b'0')];
    let y_spec = [(b'0', b'1'), (b'0', b'1'), (b'0', b'1'), (b'0', b'1')];
    let (x, y) = (gts(&x_spec), gts(&y_spec));
    // y is constant (dosage 1 for all) → zero variance → None.
    assert_eq!(pair_sums(planes(&x), planes(&y), x.len()).r2(), None);
    // A perfect inverse: y = 2 − x.
    let x_spec = [(b'0', b'0'), (b'0', b'0'), (b'1', b'1'), (b'1', b'1')];
    let y_spec = [(b'1', b'1'), (b'1', b'1'), (b'0', b'0'), (b'0', b'0')];
    let (x, y) = (gts(&x_spec), gts(&y_spec));
    assert_eq!(pair_sums(planes(&x), planes(&y), x.len()).r2(), Some(1.0));
}

#[test]
fn missing_samples_restrict_the_valid_set() {
    // The same correlation with one site's sample 1 missing in x.
    let x_spec = [(b'0', b'0'), (b'.', b'.'), (b'1', b'1'), (b'1', b'0')];
    let y_spec = [(b'0', b'0'), (b'0', b'1'), (b'1', b'1'), (b'1', b'0')];
    let (x, y) = (gts(&x_spec), gts(&y_spec));
    let sums = pair_sums(planes(&x), planes(&y), x.len());
    assert_eq!(sums.n, 3);
    let expected = ref_r2(&x, &y).expect("invariant: variance exists");
    assert!((sums.r2().expect("invariant") - expected).abs() < 1e-12);
}

#[test]
fn decay_bins_aggregate_means() {
    let mut d = LdDecayAcc::new(vec![0, 10, 100]);
    d.push(5, 0.5);
    d.push(7, 0.3);
    d.push(50, 0.9);
    d.push(100, 0.1);
    d.push(999, 0.7);
    let means = d.mean_r2();
    // Bins: [0,10) = {5,7}; [10,100) = {50}; [100,∞) = {100,999}.
    assert_eq!(d.counts(), &[2, 1, 2]);
    assert!((means[0].expect("invariant") - 0.4).abs() < 1e-12);
    assert!((means[1].expect("invariant") - 0.9).abs() < 1e-12);
    assert!((means[2].expect("invariant") - 0.4).abs() < 1e-12);

    // Merge is additive.
    let mut a = LdDecayAcc::new(vec![0, 10]);
    a.push(1, 0.5);
    let mut b = LdDecayAcc::new(vec![0, 10]);
    b.push(2, 1.5);
    a.merge(&b);
    assert_eq!(a.counts(), &[2, 0]);
    assert!((a.mean_r2()[0].expect("invariant") - 1.0).abs() < 1e-12);
}

#[test]
fn ordered_carry_pruning_keeps_first_of_high_ld() {
    // Three identical sites: the first is kept, the rest pruned.
    let spec = vec![
        (b'0', b'0'),
        (b'0', b'1'),
        (b'1', b'1'),
        (b'1', b'0'),
        (b'0', b'0'),
        (b'1', b'1'),
    ];
    let g = gts(&spec);
    let p = planes(&g);
    let mut pruner = LdPruner::new(0.8, 10, p[0].len(), g.len());
    assert!(pruner.push(0, p));
    assert!(!pruner.push(1, p), "identical site → pruned");
    assert!(!pruner.push(2, p));
    assert_eq!(pruner.kept_sites(), vec![0]);
}

#[test]
fn pruner_window_forgets_beyond_max_sites() {
    // Identical sites, max_sites = 1: after keeping site 0 and pruning
    // site 1, site 2 compares only against the most recent KEPT site
    // (still 0 — pruned sites never enter the ring), so it prunes too.
    let spec = vec![(b'0', b'1'), (b'1', b'0'), (b'0', b'0'), (b'1', b'1')];
    let g = gts(&spec);
    let p = planes(&g);
    let mut pruner = LdPruner::new(0.8, 1, p[0].len(), g.len());
    assert!(pruner.push(0, p));
    assert!(!pruner.push(1, p));
    assert!(!pruner.push(2, p));

    // A constant site (zero variance → no r²) is kept and, with
    // max_sites = 1, EVICTS site 0 from the ring: the window holds only
    // the most recent kept site.
    let constant = gts(&[(b'0', b'0'), (b'0', b'0'), (b'0', b'0'), (b'0', b'0')]);
    let q = planes(&constant);
    assert!(pruner.push(3, q), "constant site: zero variance → no r² → kept");
    assert_eq!(pruner.kept_sites(), vec![3]);
    // The next constant site also never prunes (zero variance → no r²).
    assert!(pruner.push(4, q));
    assert_eq!(pruner.kept_sites(), vec![4]);
}
