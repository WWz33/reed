//! End-to-end reader tests on hand-written fixtures.
//!
//! Positions (0-based): 2, 6, 11, 15.

use std::path::PathBuf;

use reed_core::estimators::pi::{pi, PiEstimator};
use reed_core::{Ascertainment, PopulationView, SiteCompleteness};
use reed_io::{for_each_window, read_windows, ReadError, ReadOptions};

fn fixture(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn opts(window_size: u32, step: u32) -> ReadOptions {
    ReadOptions {
        window_size,
        step,
        manifest: SiteCompleteness::VariantOnly(Ascertainment::Complete),
    }
}

#[test]
fn disjoint_windows_from_plain_vcf() {
    let windows = read_windows(&fixture("tiny.vcf"), &opts(10, 10)).unwrap();
    assert_eq!(windows.len(), 2);

    let w0 = &windows[0];
    assert_eq!((w0.start, w0.end), (0, 10));
    assert_eq!(w0.pos, vec![2, 6]);
    assert!(w0.is_variant.get(0) && w0.is_variant.get(1));
    assert_eq!(w0.n_samples, 3);

    let pops = PopulationView::from_assignment(3, &[0, 0, 0], 1).unwrap();

    let r = pi(w0, &pops, 0, PiEstimator::PerPairComparable).unwrap();
    // site pos2: 0/0,0/1,1/1 -> r3 a3, diffs 9, comps 15
    // site pos6: 0/0,./.,0/1 -> r3 a1, diffs 3, comps 6, missing pairs 9
    assert_eq!(r.counts.diffs, 12);
    assert_eq!(r.counts.comparable, 21);
    assert_eq!(r.counts.missing, 9);
    assert!((r.value().unwrap() - 12.0 / 21.0).abs() < 1e-12);

    let r = pi(w0, &pops, 0, PiEstimator::AccessibleBp).unwrap();
    assert!((r.value().unwrap() - (9.0 / 15.0 + 3.0 / 6.0) / 10.0).abs() < 1e-12);

    let w1 = &windows[1];
    assert_eq!((w1.start, w1.end), (10, 20));
    assert_eq!(w1.pos, vec![11, 15]);
    assert!(w1.is_variant.get(0) && !w1.is_variant.get(1)); // pos15 is monomorphic

    let r = pi(w1, &pops, 0, PiEstimator::PerPairComparable).unwrap();
    // pos11: 0/1,0/1,0/0 -> r4 a2, diffs 8, comps 15; pos15: all 0/0 -> 0, 15
    assert_eq!(r.counts.diffs, 8);
    assert_eq!(r.counts.comparable, 30);
    assert_eq!(r.counts.missing, 0);
    assert!((r.value().unwrap() - 8.0 / 30.0).abs() < 1e-12);
}

#[test]
fn sliding_windows_overlap_and_clamp() {
    let windows = read_windows(&fixture("tiny.vcf"), &opts(10, 4)).unwrap();
    // first record at pos 2 -> first grid start 0; starts 0, 4, 8, 12.
    assert_eq!(windows.len(), 4);
    let starts: Vec<u64> = windows.iter().map(|w| w.start).collect();
    assert_eq!(starts, vec![0, 4, 8, 12]);
    let sites: Vec<Vec<u32>> = windows.iter().map(|w| w.pos.clone()).collect();
    assert_eq!(sites, vec![vec![2, 6], vec![6, 11], vec![11, 15], vec![15]]);
    // Clamped trailing window: contig length 20, so [12, 20) in this fixture.
    assert_eq!(windows[3].end, 20);
}

#[test]
fn for_each_window_streams_with_summary() {
    let mut seen = 0u64;
    let summary = for_each_window(&fixture("tiny.vcf"), &opts(10, 10), |w| {
        seen += 1;
        assert!(w.window_bp() > 0);
        Ok(())
    })
    .unwrap();
    assert_eq!(summary, reed_io::StreamSummary { windows: 2, sites: 4 });
    assert_eq!(seen, 2);
}

#[test]
fn contig_without_length_is_an_error() {
    let err = read_windows(&fixture("tiny_nolength.vcf"), &opts(10, 10)).unwrap_err();
    match err {
        ReadError::MissingContigLength { contig } => assert_eq!(contig, "Chr01"),
        other => panic!("unexpected error: {other:?}"),
    }
}

#[test]
fn invalid_options_are_rejected() {
    let err = read_windows(&fixture("tiny.vcf"), &opts(10, 20)).unwrap_err();
    assert!(matches!(err, ReadError::User(_)));
}

#[test]
fn dot_sample_column_decodes_as_missing() {
    // A `.` sample column is spec-legal "no genotype" and must decode as ./.,
    // not abort the stream (edge-audit finding 1).
    let windows = read_windows(&fixture("dot_sample.vcf"), &opts(10, 10)).unwrap();
    assert_eq!(windows.len(), 1);
    let w = &windows[0];
    assert_eq!(w.pos, vec![2, 6]);
    // S2 = `.` at pos 2 -> missing genotype (-1,-1).
    assert_eq!(w.genotype(0, 1), (-1, -1));
    assert!(w.missing.get(1)); // site 0, sample 1
    // S1 = `.` carries no alleles; the site is variant via S3 = 0/1.
    use reed_core::estimators::pi::{pi, PiEstimator};
    let pops = PopulationView::from_assignment(3, &[0, 0, 0], 1).unwrap();
    let r = pi(w, &pops, 0, PiEstimator::PerPairComparable).unwrap();
    // pos2: S1=0/0, S2=., S3=0/1 -> called alleles 0,0,0,1: diffs 3, comps 6;
    // pos6: 0/0,0/1,1/1 -> 0,0,0,1,1,1: diffs 9, comps 15.
    assert_eq!(r.counts.diffs, 12);
    assert_eq!(r.counts.comparable, 21);
}

#[test]
fn contig_revisit_is_rejected() {
    // Chr01 records, then Chr02, then Chr01 again: spec violation.
    let err = read_windows(&fixture("revisit.vcf"), &opts(10, 10)).unwrap_err();
    assert!(matches!(err, ReadError::ContigRevisit { .. }));
}
