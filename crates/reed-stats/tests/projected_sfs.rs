//! Projected-SFS canonical tests: the ratio recursion against an
//! independent exact-binomial path, hand-computed projections, and the
//! accumulator laws on the weighted bins.

use reed_core::math::sum::NeumaierSum;
use reed_stats::sfs::{ProjectionCache, ProjectedSfsAcc, projection_weights};

/// Exact hypergeometric P(k | n, m, alt) via f64 lgamma-summed binomial
/// logs (an independent path from the ratio recursion; exact to ~1e−12 via
/// Lanczos lgamma, tolerant at 1e−9, and overflow-free).
fn exact_p(k: u32, n: u32, m: u32, alt: u32) -> f64 {
    fn ln_choose(nn: f64, kk: f64) -> f64 {
        // ln C(nn, kk) via lgamma.
        crate_lgamma(nn + 1.0) - crate_lgamma(kk + 1.0) - crate_lgamma(nn - kk + 1.0)
    }
    let ln_num = ln_choose(f64::from(alt), f64::from(k))
        + ln_choose(f64::from(m - alt), f64::from(n - k));
    let ln_den = ln_choose(f64::from(m), f64::from(n));
    (ln_num - ln_den).exp()
}

/// ln Γ(x) — Lanczos approximation (g=7, n=9), ~1e−13 relative error.
#[must_use]
fn crate_lgamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const COEF: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        std::f64::consts::PI
            / ((std::f64::consts::PI * x).sin() * crate_lgamma(1.0 - x))
    } else {
        let x1 = x - 1.0;
        let mut a = COEF[0];
        let t = x1 + G + 0.5;
        for (i, &c) in COEF.iter().enumerate().skip(1) {
            a += c / (x1 + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (t).ln() * (x1 + 0.5) - t + (a).ln()
    }
}

#[test]
fn ratio_recursion_matches_exact_binomials() {
    for (m, alt, n) in [
        (8_u32, 1_u32, 4_u32),
        (8, 4, 4),
        (8, 8, 4),
        (696, 55, 600),
        (696, 340, 600),
        (50, 7, 30),
        (1000, 999, 500),
        (12, 3, 12), // m == n: the identity projection
    ] {
        let w = projection_weights(m, alt, n);
        let total: f64 = w.iter().map(|(_, w)| w).sum();
        assert!((total - 1.0).abs() < 1e-12, "m={m} alt={alt} n={n}: sum={total}");
        // Retained ks stay within the true support and the retained mass
        // is ≈ 1 (the truncation is at the 1e−15×peak noise floor).
        let (k_lo, k_hi) = (n.saturating_sub(m - alt), n.min(alt));
        for (k, _) in &w {
            assert!(*k >= k_lo && *k <= k_hi, "k={k} outside [{k_lo},{k_hi}]");
        }
        assert!((total - 1.0).abs() < 1e-9, "retained mass: {total}");
        for (k, weight) in &w {
            let expected = exact_p(*k, n, m, alt);
            assert!(
                (weight - expected).abs() <= 1e-9 * expected.abs().max(1e-300),
                "m={m} alt={alt} n={n} k={k}: {weight} vs {expected}"
            );
        }
    }
}

#[test]
fn non_segregating_sites_carry_no_spectrum_mass() {
    // A fixed site (alt = 0) projects entirely to k = 0: the returned
    // distribution is the certain event, and the caller stores none of it
    // (k = 0 is not a segregating bin).
    assert_eq!(projection_weights(8, 0, 4), vec![(0, 1.0)]);
    // A fully-alt site (alt = m) projects to k = n: same shape.
    assert_eq!(projection_weights(8, 8, 4), vec![(4, 1.0)]);
}

#[test]
fn identity_projection_is_the_exact_spectrum() {
    // m == n: every site contributes all its mass to k = alt.
    let w = projection_weights(12, 5, 12);
    assert_eq!(w.len(), 1);
    assert_eq!(w[0], (5, 1.0));
}

#[test]
fn hand_computed_projection() {
    // m=4, alt=2, n=2: P(0) = C(2,0)·C(2,2)/C(4,2) = 1/6;
    // P(1) = C(2,1)·C(2,1)/C(4,2) = 4/6; P(2) = 1/6. The retained
    // segregating weight is 4/6 — no renormalization.
    let w = projection_weights(4, 2, 2);
    assert_eq!(w, vec![(0, 1.0 / 6.0), (1, 4.0 / 6.0), (2, 1.0 / 6.0)]);
}

#[test]
fn accumulator_projects_and_merges() {
    let mut cache = ProjectionCache::new(4);
    let mut a = ProjectedSfsAcc::new(0, 4, false).expect("invariant: n ≥ 2");
    let mut b = ProjectedSfsAcc::new(0, 4, false).expect("invariant: n ≥ 2");
    let mut whole = ProjectedSfsAcc::new(0, 4, false).expect("invariant: n ≥ 2");

    // Site A: m=4, alt=1 (k=1 → weight 1). Site B: m=4, alt=2.
    let wa = cache.weights(4, 1);
    a.push_weighted(Some(&wa));
    whole.push_weighted(Some(&wa));
    let wb = cache.weights(4, 2);
    b.push_weighted(Some(&wb));
    whole.push_weighted(Some(&wb));
    // Site C: m=2 < n → skipped.
    let wc = cache.weights(2, 1);
    a.push_weighted(None);
    whole.push_weighted(None);

    a.merge(&b);
    assert_eq!(a.site_tally(), (2, 1));
    assert_eq!(whole.site_tally(), (2, 1));
    for k in 1..4 {
        assert!((a.xi(k) - whole.xi(k)).abs() <= 4.0 * f64::EPSILON, "k={k}");
    }
    // Site A contributes ξ_1 = 1; site B (m=4, alt=2, n=4) has support
    // k ∈ [2, 2] with P = 1 → ξ_2 = 1.
    assert!((a.xi(1) - 1.0).abs() < 1e-12);
    assert!((a.xi(2) - 1.0).abs() < 1e-12);
    assert!((a.xi(3) - 0.0).abs() < 1e-12);
}

#[test]
fn cache_returns_identical_weights() {
    let mut cache = ProjectionCache::new(600);
    let first = cache.weights(696, 55);
    for _ in 0..10 {
        assert!(std::sync::Arc::ptr_eq(&cache.weights(696, 55), &first));
    }
    assert_eq!(cache.n(), 600);
}

#[test]
fn projected_tajima_on_known_spectrum() {
    // The neutral spectrum ξ_k = 1/k has π = S/a₁ exactly → D = 0.
    let mut acc = ProjectedSfsAcc::new(0, 10, false).expect("invariant: n ≥ 2");
    for k in 1..10 {
        acc.push_weighted(Some(&[(k, 1.0 / f64::from(k))]));
    }
    let d = acc.tajima_d(0).expect("invariant: variance > 0");
    assert!(d.abs() < 1e-9, "neutral spectrum: D={d}");
    // Pure het excess (all mass at high k) → negative D.
    let mut excess = ProjectedSfsAcc::new(0, 10, false).expect("invariant: n ≥ 2");
    excess.push_weighted(Some(&[(9, 10.0)]));
    let d_excess = excess.tajima_d(0).expect("invariant: variance > 0");
    assert!(d_excess < 0.0, "het excess: D={d_excess}");
}

#[allow(dead_code)]
fn _neumaier_visible(n: NeumaierSum) -> NeumaierSum {
    n
}
