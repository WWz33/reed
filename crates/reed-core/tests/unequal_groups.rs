//! Golden tests on randomized fixtures with **unequal group sizes and missing
//! data** — the regime where per-site normalisation (allel semantics) and
//! pooled ratios diverge. Reference values generated with scikit-allel
//! 1.3.13 (`default_rng(42)` / `default_rng(7)`, scripts preserved in the
//! DESIGN.md §12 workflow); embedded verbatim.
//!
//! These tests were added after a three-agent audit caught the original
//! Hudson/WC84 implementations agreeing with allel only in the degenerate
//! equal-coverage regime.

use reed_core::estimators::fst::{fst, FstEstimator};
use reed_core::estimators::pbs::{pbs, PbsEstimator};
use reed_core::{Bits, ContigId, PopulationView, SiteCompleteness, WindowBuffer};

const EPS: f64 = 1e-12;

/// Build a buffer from a flat genotype string (`-1` allowed for missing).
fn build(n_samples: u32, sites: &[&str]) -> WindowBuffer {
    let mut pos = Vec::new();
    let mut is_variant = Vec::new();
    let mut genotypes = Vec::new();
    let mut missing = Vec::new();
    let mut phased = Vec::new();
    for (i, site) in sites.iter().enumerate() {
        pos.push(i as u32);
        let mut variant = false;
        let alleles: Vec<i8> = site
            .split_whitespace()
            .map(|x| {
                let v: i8 = x.parse().unwrap();
                if v > 0 {
                    variant = true;
                }
                v
            })
            .collect();
        assert_eq!(alleles.len(), n_samples as usize * 2);
        genotypes.extend_from_slice(&alleles);
        for g in alleles.chunks(2) {
            missing.push(g[0] == -1 || g[1] == -1);
            phased.push(false);
        }
        is_variant.push(variant);
    }
    WindowBuffer::new(
        SiteCompleteness::WithInvariants,
        ContigId(0),
        100,
        0,
        100,
        n_samples,
        pos,
        Bits::from_bools(&is_variant),
        genotypes,
        Bits::from_bools(&missing),
        Bits::from_bools(&phased),
    )
    .unwrap()
}

fn assert_close(a: Option<f64>, b: f64) {
    let a = a.expect("expected Some(value)");
    assert!((a - b).abs() < EPS, "expected {b}, got {a}");
}

/// numpy default_rng(42).integers(0,2,(3,10,2)), 15% missing; A={0,1,2},
/// B={3..9} (3 vs 7 samples).
const RNG42_SITES: [&str; 3] = [
    "0 1 1 0 0 1 -1 1 0 0 1 1 1 1 1 1 1 0 1 -1",
    "1 -1 0 1 1 1 0 1 -1 -1 0 0 0 1 1 0 1 1 -1 1",
    "0 1 1 0 -1 1 0 1 1 1 1 0 0 0 -1 -1 1 0 1 1",
];

/// numpy default_rng(7).integers(0,2,(6,12,2)), 12% missing; A={0..3},
/// B={4..8}, C={9..11} (4/5/3 samples).
const RNG7_SITES: [&str; 6] = [
    "1 1 1 1 1 1 1 0 0 0 0 1 1 0 0 1 0 1 -1 0 1 0 0 0",
    "-1 0 -1 0 0 1 1 1 1 1 1 1 1 1 0 1 0 0 -1 0 1 1 0 0",
    "-1 0 0 1 1 0 1 1 1 1 0 1 0 0 0 0 1 0 0 0 1 1 1 0",
    "1 0 -1 0 1 1 1 0 1 0 1 1 0 1 1 -1 1 1 0 -1 0 0 0 1",
    "1 1 1 1 1 0 -1 1 0 0 1 -1 1 0 0 0 0 -1 -1 -1 1 1 0 1",
    "-1 1 0 1 -1 1 1 0 1 0 0 0 0 -1 1 0 0 -1 1 0 -1 1 0 0",
];

#[test]
fn hudson_fst_matches_allel_under_unequal_coverage() {
    let w = build(10, &RNG42_SITES);
    let pops = PopulationView::from_assignment(10, &[0, 0, 0, 1, 1, 1, 1, 1, 1, 1], 2).unwrap();
    let f = fst(&w, &pops, 0, 1, FstEstimator::Hudson).unwrap();
    // allel: num=-0.08636363636363636, den=1.456060606060606
    assert!((f.hudson_num.unwrap() - (-0.08636363636363636)).abs() < EPS);
    assert!((f.hudson_den.unwrap() - 1.456060606060606).abs() < EPS);
    assert_close(f.value(), -0.05931321540062435);
    // The pooled-ratio audit trail must also be self-consistent:
    // diffs/comparable ≤ 1 per term.
    assert!(f.between.diffs <= f.between.comparable);

    // Regression guard for the audit field previously left at zero.
    assert!(f.within_a.missing + f.within_b.missing + f.between.missing > 0);
}

#[test]
fn weir_cockerham_matches_allel_under_unequal_coverage() {
    let w = build(10, &RNG42_SITES);
    let pops = PopulationView::from_assignment(10, &[0, 0, 0, 1, 1, 1, 1, 1, 1, 1], 2).unwrap();
    let f = fst(&w, &pops, 0, 1, FstEstimator::WeirCockerham84).unwrap();
    let wc = f.wc.unwrap();
    // allel: a=-0.0877475667204709, b=-0.12224304293099507, c=1.6408730158730158
    assert!((wc.a - (-0.0877475667204709)).abs() < EPS);
    assert!((wc.b - (-0.12224304293099507)).abs() < EPS);
    assert!((wc.c - 1.6408730158730158).abs() < EPS);
    assert_eq!(wc.n_sites, 3);
}

#[test]
fn pbs_matches_allel_under_unequal_coverage() {
    let w = build(12, &RNG7_SITES);
    let assignment: Vec<u32> = (0..12)
        .map(|i| match i {
            0..=3 => 0,
            4..=8 => 1,
            _ => 2,
        })
        .collect();
    let pops = PopulationView::from_assignment(12, &assignment, 3).unwrap();

    let r = pbs(&w, &pops, 0, 1, 2, PbsEstimator::Yi2010).unwrap();
    // allel: pbs(..., normed=False) = 0.06629918897511186
    assert!((r.pbs.unwrap() - 0.06629918897511186).abs() < 1e-12);
    assert!(r.fst_ab.unwrap() > 0.0);

    let r = pbs(&w, &pops, 0, 1, 2, PbsEstimator::PbsN1).unwrap();
    // allel: pbs(..., normed=True) = 0.057258307386842185
    assert!((r.pbs.unwrap() - 0.057258307386842185).abs() < 1e-12);

    // Between-pair integer audit trail is now retained.
    assert!(r.between_ab.comparable + r.between_ac.comparable + r.between_bc.comparable > 0);
}
