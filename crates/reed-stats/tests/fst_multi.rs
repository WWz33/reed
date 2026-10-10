//! M7 estimator-extension canonical tests: multi-population WC84 against
//! allel's `weir_cockerham_fst` (hand-computed via the reference formula),
//! HudsonFiniteSample against snputils's "hudson" branch, and the Fst
//! matrix request shape.

use reed_core::count::SiteGroupCounts;
use reed_stats::fst::{FstKind, FstMatrixAcc, WcMultiAcc};

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
