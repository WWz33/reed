//! Hand-computed golden tests for the core estimators.
//!
//! Fixture: 4 diploid samples, groups A = {0, 1}, B = {2, 3}, window [0, 10)
//! on a contig of length 100.
//!
//! | site | pos | S0   | S1   | S2   | S3   |
//! |------|-----|------|------|------|------|
//! | 0    | 0   | 0/0  | 0/1  | 1/1  | 0/1  |
//! | 1    | 4   | 0/0  | 0/0  | 0/1  | 0/0  |
//! | 2    | 8   | 0/0  | 0/0  | 0/0  | 0/0  |

use reed_oracle::estimators::allele_counts::allele_counts;
use reed_oracle::estimators::dxy::{dxy, DxyEstimator};
use reed_oracle::estimators::fst::{fst, FstEstimator};
use reed_oracle::estimators::pi::{pi, PiEstimator};
use reed_oracle::{Ascertainment, Bits, ContigId, PopulationView, SiteCompleteness, WindowBuffer, WindowError};

const EPS: f64 = 1e-12;

fn assert_close(a: Option<f64>, b: f64) {
    let a = a.expect("expected Some(value)");
    assert!((a - b).abs() < EPS, "expected {b}, got {a}");
}

fn build(manifest: SiteCompleteness, sites: &[(u32, [[i8; 2]; 4])]) -> WindowBuffer {
    let n_samples = 4u32;
    let mut pos = Vec::new();
    let mut is_variant = Vec::new();
    let mut genotypes = Vec::new();
    let mut missing = Vec::new();
    let mut phased = Vec::new();
    for &(p, gts) in sites {
        pos.push(p);
        let mut variant = false;
        for g in gts {
            for a in g {
                genotypes.push(a);
                if a > 0 {
                    variant = true;
                }
            }
            missing.push(g[0] == -1 || g[1] == -1);
            phased.push(false);
        }
        is_variant.push(variant);
    }
    WindowBuffer::new(
        manifest,
        ContigId(0),
        100,
        0,
        10,
        n_samples,
        pos,
        Bits::from_bools(&is_variant),
        genotypes,
        Bits::from_bools(&missing),
        Bits::from_bools(&phased),
    )
    .unwrap()
}

fn fixture(manifest: SiteCompleteness) -> WindowBuffer {
    build(
        manifest,
        &[
            (0, [[0, 0], [0, 1], [1, 1], [0, 1]]),
            (4, [[0, 0], [0, 0], [0, 1], [0, 0]]),
            (8, [[0, 0], [0, 0], [0, 0], [0, 0]]),
        ],
    )
}

fn pops() -> PopulationView {
    PopulationView::from_assignment(4, &[0, 0, 1, 1], 2).unwrap()
}

#[test]
fn pi_per_pair_comparable_counts() {
    let w = fixture(SiteCompleteness::WithInvariants);
    let p = pops();

    let a = pi(&w, &p, 0, PiEstimator::PerPairComparable).unwrap();
    assert_eq!(a.counts.diffs, 3);
    assert_eq!(a.counts.comparable, 18);
    assert_eq!(a.counts.missing, 0);
    assert_eq!((a.n_variant_sites, a.n_invariant_sites), (2, 1));
    assert_close(a.value(), 3.0 / 18.0);

    let b = pi(&w, &p, 1, PiEstimator::PerPairComparable).unwrap();
    assert_eq!(b.counts.diffs, 6);
    assert_eq!(b.counts.comparable, 18);
    assert_close(b.value(), 6.0 / 18.0);
}

#[test]
fn pi_accessible_bp_uses_window_length() {
    let w = fixture(SiteCompleteness::WithInvariants);
    let p = pops();

    let a = pi(&w, &p, 0, PiEstimator::AccessibleBp).unwrap();
    assert_close(a.value(), 0.5 / 10.0); // site means: 3/6, 0, 0

    let b = pi(&w, &p, 1, PiEstimator::AccessibleBp).unwrap();
    assert_close(b.value(), 1.0 / 10.0); // site means: 3/6, 3/6, 0
}

#[test]
fn dxy_per_pair_and_accessible() {
    let w = fixture(SiteCompleteness::WithInvariants);
    let p = pops();

    let d = dxy(&w, &p, 0, 1, DxyEstimator::PerPairComparable).unwrap();
    assert_eq!(d.counts.diffs, 14);
    assert_eq!(d.counts.comparable, 48);
    assert_close(d.value(), 14.0 / 48.0);

    let d = dxy(&w, &p, 0, 1, DxyEstimator::AccessibleBp).unwrap();
    assert_close(d.value(), (10.0 / 16.0 + 4.0 / 16.0) / 10.0);
}

#[test]
fn hudson_fst_ratio_of_sums() {
    let w = fixture(SiteCompleteness::WithInvariants);
    let p = pops();

    let f = fst(&w, &p, 0, 1, FstEstimator::Hudson).unwrap();
    assert_eq!(f.within_a.diffs, 3);
    assert_eq!(f.within_a.comparable, 18);
    assert_eq!(f.within_b.diffs, 6);
    assert_eq!(f.within_b.comparable, 18);
    assert_eq!(f.between.diffs, 14);
    assert_eq!(f.between.comparable, 48);
    assert_close(f.value(), 1.0 - (9.0 / 36.0) / (14.0 / 48.0));
}

#[test]
fn missing_data_is_pairwise_deleted() {
    // Group {0, 1}, one site: S0 = 0/0, S1 = ./. .
    let w = build(
        SiteCompleteness::WithInvariants,
        &[(0, [[0, 0], [-1, -1], [0, 0], [0, 0]])],
    );
    let p = PopulationView::from_assignment(4, &[0, 0, 1, 1], 2).unwrap();

    let r = pi(&w, &p, 0, PiEstimator::PerPairComparable).unwrap();
    assert_eq!(r.counts.diffs, 0);
    assert_eq!(r.counts.comparable, 1); // only S0's two chromosomes
    assert_eq!(r.counts.missing, 5); // C(4,2) - 1
    assert_close(r.value(), 0.0);
}

#[test]
fn manifest_gates_pi_but_only_flags_fst() {
    let p = pops();

    let unknown = fixture(SiteCompleteness::VariantOnly(Ascertainment::Unknown));
    assert!(pi(&unknown, &p, 0, PiEstimator::PerPairComparable).is_err());
    assert!(pi(&unknown, &p, 0, PiEstimator::AccessibleBp).is_err());
    let f = fst(&unknown, &p, 0, 1, FstEstimator::Hudson).unwrap();
    assert!(f.bias.is_some());

    let complete = fixture(SiteCompleteness::VariantOnly(Ascertainment::Complete));
    let r = pi(&complete, &p, 0, PiEstimator::PerPairComparable).unwrap();
    assert!(r.bias.is_some()); // biased but computable
    assert_close(r.value(), 3.0 / 18.0);
    let r = pi(&complete, &p, 0, PiEstimator::AccessibleBp).unwrap();
    assert!(r.bias.is_none());

    let full = fixture(SiteCompleteness::WithInvariants);
    let r = pi(&full, &p, 0, PiEstimator::PerPairComparable).unwrap();
    assert!(r.bias.is_none());
}

#[test]
fn aggregate_allele_counts() {
    let w = fixture(SiteCompleteness::WithInvariants);
    let p = pops();
    let counts = allele_counts(&w, &p).unwrap();
    assert_eq!(counts.len(), 2);
    assert_eq!(
        (counts[0].ref_count, counts[0].alt_count, counts[0].missing),
        (11, 1, 0)
    );
    assert_eq!(
        (counts[1].ref_count, counts[1].alt_count, counts[1].missing),
        (8, 4, 0)
    );
}

#[test]
fn constructor_rejects_invalid_buffers() {
    // Zero contig length.
    let e = WindowBuffer::new(
        SiteCompleteness::WithInvariants,
        ContigId(0),
        0,
        0,
        10,
        1,
        vec![0],
        Bits::new(1),
        vec![0, 0],
        Bits::new(1),
        Bits::new(1),
    );
    assert_eq!(e.unwrap_err(), WindowError::ZeroContigLength);

    // Missing-bit mismatch: sentinel says missing, bit says not.
    let e = WindowBuffer::new(
        SiteCompleteness::WithInvariants,
        ContigId(0),
        100,
        0,
        10,
        1,
        vec![0],
        Bits::new(1),
        vec![-1, 0],
        Bits::from_bools(&[false]),
        Bits::new(1),
    );
    assert_eq!(e.unwrap_err(), WindowError::MissingBitMismatch { index: 0 });

    // Multiallelic allele index rejected.
    let e = WindowBuffer::new(
        SiteCompleteness::WithInvariants,
        ContigId(0),
        100,
        0,
        10,
        1,
        vec![0],
        Bits::from_bools(&[true]),
        vec![0, 2],
        Bits::from_bools(&[false]),
        Bits::new(1),
    );
    assert_eq!(
        e.unwrap_err(),
        WindowError::AlleleOutOfRange { index: 0, allele: 2 }
    );

    // Position outside window.
    let e = WindowBuffer::new(
        SiteCompleteness::WithInvariants,
        ContigId(0),
        100,
        0,
        10,
        1,
        vec![10],
        Bits::new(1),
        vec![0, 0],
        Bits::new(1),
        Bits::new(1),
    );
    assert_eq!(
        e.unwrap_err(),
        WindowError::PositionOutOfWindow {
            index: 0,
            pos: 10,
            start: 0,
            end: 10
        }
    );
}
