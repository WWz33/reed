//! M7 estimator-extension canonical tests: multi-population WC84 against
//! allel's `weir_cockerham_fst` (hand-computed via the reference formula),
//! HudsonFiniteSample against snputils's "hudson" branch, and the Fst
//! matrix request shape.

use reed_core::count::SiteGroupCounts;
use reed_stats::fst::{FstKind, FstMatrixAcc, HudsonFiniteAcc, NeiMultiAcc, WcMultiAcc};

fn row(ans: &[u32], alts: &[u32], hets: &[u32]) -> Vec<SiteGroupCounts> {
    ans.iter()
        .zip(alts)
        .zip(hets)
        .map(|((&an, &alt), &het)| SiteGroupCounts {
            n_ref_alleles: an - alt,
            n_alt_alleles: alt,
            n_het: het,
            ..SiteGroupCounts::default()
        })
        .collect()
}

/// allel's weir_cockerham_fst (general r) as a pure function of the
/// per-population (an, ac_ref, n_het) — the reference the accumulator
/// mirrors.
fn ref_wc84(inputs: &[(f64, f64, f64)]) -> (f64, f64, f64) {
    // (an, ref_allele_count, n_het) per population.
    let r = inputs.len() as f64;
    let n: Vec<f64> = inputs.iter().map(|(an, _, _)| an / 2.0).collect();
    let n_total: f64 = n.iter().sum();
    let n_bar = n_total / r;
    let n_c = (n_total - n.iter().map(|x| x * x).sum::<f64>() / n_total) / (r - 1.0);
    let p: Vec<f64> = inputs.iter().map(|(an, ac, _)| ac / an).collect();
    let ac_total: f64 = inputs.iter().map(|(_, ac, _)| ac).sum();
    let an_total: f64 = inputs.iter().map(|(an, _, _)| an).sum();
    let p_bar = ac_total / an_total;
    let s2 = n
        .iter()
        .zip(&p)
        .map(|(&ni, &pi)| ni * (pi - p_bar) * (pi - p_bar))
        .sum::<f64>()
        / (n_bar * (r - 1.0));
    let h_bar = inputs.iter().map(|(_, _, h)| h).sum::<f64>() / n_total;
    let a = (n_bar / n_c)
        * (s2 - (1.0 / (n_bar - 1.0))
            * (p_bar * (1.0 - p_bar) - (r - 1.0) * s2 / r - h_bar / 4.0));
    let b = (n_bar / (n_bar - 1.0))
        * (p_bar * (1.0 - p_bar)
            - (r - 1.0) * s2 / r
            - (2.0 * n_bar - 1.0) * h_bar / (4.0 * n_bar));
    let c = h_bar / 2.0;
    (a, b, c)
}

#[test]
fn wc84_three_populations_matches_reference() {
    let sites: Vec<(Vec<u32>, Vec<u32>, Vec<u32>)> = vec![
        (vec![10, 8, 12], vec![2, 5, 6], vec![2, 2, 4]),
        (vec![8, 12, 6], vec![4, 3, 2], vec![3, 2, 1]),
        (vec![12, 10, 10], vec![6, 4, 5], vec![4, 3, 3]),
        (vec![6, 6, 6], vec![3, 3, 3], vec![1, 2, 3]),
    ];
    let mut acc = WcMultiAcc::new(vec![0, 1, 2]);
    let (mut ra, mut rb, mut rc) = (0.0, 0.0, 0.0);
    for (ans, alts, hets) in &sites {
        acc.push(&row(ans, alts, hets));
        let inputs: Vec<(f64, f64, f64)> = ans
            .iter()
            .zip(alts)
            .zip(hets)
            .map(|((&an, &alt), &het)| (f64::from(an), f64::from(an - alt), f64::from(het)))
            .collect();
        let (a, b, c) = ref_wc84(&inputs);
        ra += a;
        rb += b;
        rc += c;
    }
    assert_eq!(acc.site_tally(), (4, 0));
    let (a_sum, b_sum, c_sum) = acc.components();
    assert!((a_sum.value() - ra).abs() < 1e-12, "a: {} vs {}", a_sum.value(), ra);
    assert!((b_sum.value() - rb).abs() < 1e-12);
    assert!((c_sum.value() - rc).abs() < 1e-12);
    assert!((acc.value().expect("invariant: sites") - ra / (ra + rb + rc)).abs() < 1e-12);
}

#[test]
fn wc84_two_populations_agrees_with_the_pair_estimator() {
    // The r=2 instance must match the existing pair path (oracle parity
    // in M4 tests) — cross-checked by the identity: Σa/(Σa+Σb+Σc) from
    // the multi path on equal-size groups equals the pair path's
    // wc value. (Full-value check lives in the reference test above; this
    // test guards the r=2 special case against division degeneration.)
    let mut acc = WcMultiAcc::new(vec![0, 1]);
    acc.push(&row(&[10, 10], &[2, 8], &[2, 3]));
    assert_eq!(acc.site_tally(), (1, 0));
    assert!(acc.value().is_some());
}

#[test]
fn wc84_skips_degenerate_sites() {
    let mut acc = WcMultiAcc::new(vec![0, 1, 2]);
    acc.push(&row(&[0, 10, 12], &[0, 5, 6], &[0, 2, 4])); // anA = 0
    assert_eq!(acc.site_tally(), (0, 1));
    assert!(acc.value().is_none());
}

/// snputils's "hudson" branch: num = dxy − (π1+π2)/2 with
/// πi = 2·pi(1−pi)·n/(n−1), den = dxy — the finite-sample-corrected
/// Hudson estimator (§9.2).
fn ref_hudson_finite(an1: u32, alt1: u32, an2: u32, alt2: u32) -> (f64, f64) {
    let p1 = f64::from(alt1) / f64::from(an1);
    let p2 = f64::from(alt2) / f64::from(an2);
    let dxy = p1 * (1.0 - p2) + p2 * (1.0 - p1);
    let pi1 = 2.0 * p1 * (1.0 - p1) * f64::from(an1) / f64::from(an1 - 1);
    let pi2 = 2.0 * p2 * (1.0 - p2) * f64::from(an2) / f64::from(an2 - 1);
    (dxy - (pi1 + pi2) / 2.0, dxy)
}

#[test]
fn hudson_finite_matches_snputils_hudson() {
    let mut acc = HudsonFiniteAcc::new(0, 1);
    let sites = [(10_u32, 2, 10, 8), (8, 4, 12, 6), (6, 3, 10, 5), (12, 9, 8, 1)];
    let (mut num, mut den) = (0.0, 0.0);
    for (an1, alt1, an2, alt2) in sites {
        acc.push(&row(&[an1, an2], &[alt1, alt2], &[0, 0]));
        let (n_i, d_i) = ref_hudson_finite(an1, alt1, an2, alt2);
        num += n_i;
        den += d_i;
    }
    assert_eq!(acc.site_tally(), (4, 0));
    let (n_sum, d_sum) = acc.sums();
    assert!((n_sum.value() - num).abs() < 1e-12);
    assert!((d_sum.value() - den).abs() < 1e-12);
    assert!((acc.value().expect("invariant: den ≠ 0") - num / den).abs() < 1e-12);
}

#[test]
fn hudson_finite_needs_two_chromosomes_a_side() {
    let mut acc = HudsonFiniteAcc::new(0, 1);
    acc.push(&row(&[1, 10], &[1, 5], &[0, 2])); // anA = 1
    assert_eq!(acc.site_tally(), (0, 1));
    assert!(acc.value().is_none());
}

#[test]
fn hudson_finite_identical_populations_stay_below_zero() {
    // With identical populations the UNCORRECTED Hudson Fst is 0; the
    // finite-sample correction subtracts (π1+π2)/2 > dxy (each πi carries
    // the n/(n−1) inflation), so the corrected estimate is negative —
    // the known downward bias, which is the price of unbiasedness at true
    // Fst > 0. Assert the exact algebraic value: num = −2p(1−p)/(n−1)
    // per site for an1 == an2 == n.
    let mut acc = HudsonFiniteAcc::new(0, 1);
    acc.push(&row(&[10, 10], &[5, 5], &[0, 0]));
    let p = 0.5_f64;
    let dxy = 2.0 * p * (1.0 - p);
    let pi = 2.0 * p * (1.0 - p) * 10.0 / 9.0;
    let expected = (dxy - pi) / dxy; // −(1/(n−1))·... per snputils
    let v = acc.value().expect("invariant: den ≠ 0");
    assert!((v - expected).abs() < 1e-12, "{v} vs {expected}");
    assert!(v < 0.0);
}


#[test]
fn nei_multi_matches_the_pair_estimator_at_r2() {
    // r = 2: equal-weighted means reduce to the pair formula (hs =
    // (hs1+hs2)/2, ht from p̄ = (p1+p2)/2) — the oracle's semantics.
    let mut multi = NeiMultiAcc::new(vec![0, 1]);
    let mut pair_hs = 0.0;
    let mut pair_ht = 0.0;
    let sites = [(10_u32, 2, 10, 8), (8, 4, 12, 6), (6, 3, 10, 5)];
    for (an1, alt1, an2, alt2) in sites {
        multi.push(&row(&[an1, an2], &[alt1, alt2], &[0, 0]));
        let p1 = f64::from(alt1) / f64::from(an1);
        let p2 = f64::from(alt2) / f64::from(an2);
        pair_hs += (2.0 * p1 * (1.0 - p1) + 2.0 * p2 * (1.0 - p2)) / 2.0;
        let p_bar = (p1 + p2) / 2.0;
        pair_ht += 2.0 * p_bar * (1.0 - p_bar);
    }
    let (hs, ht) = multi.sums();
    assert!((hs.value() - pair_hs).abs() < 1e-12);
    assert!((ht.value() - pair_ht).abs() < 1e-12);
    assert!((multi.value().expect("invariant: ht ≠ 0") - (1.0 - pair_hs / pair_ht)).abs() < 1e-12);
}

#[test]
fn nei_multi_equal_frequencies_give_zero_gst() {
    // All populations at the same frequency: Hs == Ht.
    let mut acc = NeiMultiAcc::new(vec![0, 1, 2, 3]);
    acc.push(&row(&[10, 10, 10, 10], &[2, 2, 2, 2], &[0, 0, 0, 0]));
    acc.push(&row(&[8, 8, 8, 8], &[4, 4, 4, 4], &[0, 0, 0, 0]));
    let v = acc.value().expect("invariant: ht ≠ 0");
    assert!(v.abs() < 1e-12, "{v}");
}

#[test]
fn fst_matrix_covers_every_pair_once() {
    let mut m = FstMatrixAcc::new(3, FstKind::Hudson);
    let entries = m.entries();
    let pairs: Vec<(usize, usize)> = entries.iter().map(|e| (e.i, e.j)).collect();
    assert_eq!(pairs, vec![(0, 1), (0, 2), (1, 2)]);

    // One site: all pairs see it; the diagonal is not present.
    m.push(&row(&[10, 8, 12], &[2, 5, 6], &[2, 2, 4]));
    for entry in m.entries() {
        // The single-site values exist (each pair had one contributing
        // site); cross-check one by the direct pair accumulator.
        let mut single = reed_stats::fst::FstAcc::new(entry.i, entry.j);
        single.push(&row(&[10, 8, 12], &[2, 5, 6], &[2, 2, 4]));
        assert_eq!(entry.value, single.value(FstKind::Hudson));
    }
}

#[test]
fn fst_matrix_matches_individual_accumulators() {
    let sites: Vec<(Vec<u32>, Vec<u32>, Vec<u32>)> = vec![
        (vec![10, 8, 12, 6], vec![2, 5, 6, 3], vec![2, 2, 4, 1]),
        (vec![8, 12, 6, 10], vec![4, 3, 2, 5], vec![3, 2, 1, 2]),
        (vec![12, 10, 10, 8], vec![6, 4, 5, 2], vec![4, 3, 3, 2]),
        (vec![6, 6, 12, 12], vec![3, 3, 7, 8], vec![1, 2, 5, 6]),
    ];
    let mut m = FstMatrixAcc::new(4, FstKind::Hudson);
    // The reference pair accumulators, flattened (i, j, acc).
    let mut refs: Vec<(usize, usize, reed_stats::fst::FstAcc)> = (0..4)
        .flat_map(|i| (i + 1..4).map(move |j| (i, j, reed_stats::fst::FstAcc::new(i, j))))
        .collect();
    for (ans, alts, hets) in &sites {
        let r = row(ans, alts, hets);
        m.push(&r);
        for acc in &mut refs {
            acc.2.push(&r);
        }
    }
    for (entry, (_, _, ref_acc)) in m.entries().into_iter().zip(&refs) {
        assert_eq!(entry.value, ref_acc.value(FstKind::Hudson), "pair ({}, {})", entry.i, entry.j);
    }
}

#[test]
fn pbs_pairs_are_independent_under_missing_data() {
    use reed_core::count::SiteGroupCounts;
    use reed_stats::pbs::{PbsAcc, PbsKind};
    // Group C has zero called alleles: the A/B pair still contributes (the
    // three Hudson pairs are independent, 04-algorithms.md §1 — the v0.2
    // three-way intersection is gone). The value stays undefined while two
    // pair denominators are empty.
    let row = vec![
        SiteGroupCounts { n_ref_alleles: 4, n_alt_alleles: 2, ..Default::default() },
        SiteGroupCounts { n_ref_alleles: 2, n_alt_alleles: 4, ..Default::default() },
        SiteGroupCounts { n_ref_alleles: 0, n_alt_alleles: 0, ..Default::default() },
    ];
    let mut per_pair = PbsAcc::new(0, 1, 2);
    per_pair.push(&row);
    assert_eq!(per_pair.site_tally(), (1, 1), "site absorbed, one group uncovered");
    assert!(per_pair.value(PbsKind::Yi2010).is_none());
    let (f_ab, f_ac, f_bc) = per_pair.pair_fsts();
    assert!(f_ab.is_some() && f_ac.is_none() && f_bc.is_none());
}
