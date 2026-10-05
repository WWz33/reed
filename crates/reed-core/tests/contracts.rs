//! v0.2 contract tests: audit-trail completeness, coordinate metadata,
//! availability gates, and population/window validation.
//!
//! Expected numbers are hand-computed; gate behaviour follows the §4.2 matrix
//! plus the per-estimator 判定 recorded in DESIGN §12.3/§16.1.

use reed_core::estimators::allele_counts::allele_counts;
use reed_core::estimators::fstats::{d_stat, f2};
use reed_core::estimators::fst::{fst, FstEstimator};
use reed_core::estimators::hwe::{heterozygosity, hwe_exact};
use reed_core::estimators::ld::{ld_r2_pairs, ld_r2_summary};
use reed_core::estimators::pbs::{pbs, PbsEstimator};
use reed_core::estimators::pi::{pi, PiEstimator};
use reed_core::estimators::spectrum::{sfs, tajima_d, watterson_theta};
use reed_core::{
    Ascertainment, Bits, ContigId, PopError, PopulationView, SiteCompleteness, WindowBuffer,
};

const EPS: f64 = 1e-12;

fn build(
    manifest: SiteCompleteness,
    n_samples: u32,
    sites: &[(u32, Vec<[i8; 2]>)],
) -> WindowBuffer {
    let mut pos = Vec::new();
    let mut is_variant = Vec::new();
    let mut genotypes = Vec::new();
    let mut missing = Vec::new();
    let mut phased = Vec::new();
    for (p, gts) in sites {
        pos.push(*p);
        let mut variant = false;
        for g in gts {
            for a in g {
                genotypes.push(*a);
                if *a > 0 {
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
    .expect("invariant: valid test fixture")
}

/// The shared 4-sample fixture: groups A={0,1}, B={2,3}.
fn fixture() -> WindowBuffer {
    build(
        SiteCompleteness::WithInvariants,
        4,
        &[(0, vec![[0, 0], [0, 1], [1, 1], [0, 1]])],
    )
}

fn pops4() -> PopulationView {
    PopulationView::from_assignment(4, &[0, 0, 1, 1], 2).expect("invariant: valid test fixture")
}

#[test]
fn pbs_audit_trail_covers_skipped_and_partial_sites() {
    // One site where C is fully missing: the estimator skips it, but the
    // pooled audit trail must still account for every pair (§12.3-1).
    let w = build(
        SiteCompleteness::WithInvariants,
        6,
        &[(0, vec![[0, 0], [0, 1], [1, 1], [0, 0], [-1, -1], [-1, -1]])],
    );
    let pops = PopulationView::from_assignment(6, &[0, 0, 1, 1, 2, 2], 3).expect("invariant: valid test fixture");
    let r = pbs(&w, &pops, 0, 1, 2, PbsEstimator::Yi2010).expect("invariant: valid test fixture");
    assert_eq!((r.n_sites, r.n_skipped), (0, 1));
    // A: fully called (4 chroms, 6 pairs) — nothing missing.
    assert_eq!((r.within_a.comparable, r.within_a.missing), (6, 0));
    // C: fully missing — all 6 within pairs and all 16 cross pairs lost.
    assert_eq!((r.within_c.comparable, r.within_c.missing), (0, 6));
    assert_eq!((r.between_ac.comparable, r.between_ac.missing), (0, 16));
    assert_eq!((r.between_bc.comparable, r.between_bc.missing), (0, 16));
    // A×B fully comparable; diffs = 3·1 + 1·3? No — ref_a=3, alt_a=1,
    // ref_b=2, alt_b=2 → 3·2 + 1·2 = 8.
    assert_eq!((r.between_ab.comparable, r.between_ab.missing), (16, 0));
    assert_eq!(r.between_ab.diffs, 8);

    // The audit trail must agree with a standalone fst(A, B) call on the same
    // window (fst.rs accumulates over all sites too).
    let f = fst(&w, &pops, 0, 1, FstEstimator::Hudson).expect("invariant: valid test fixture");
    assert_eq!(f.within_a, r.within_a);
    assert_eq!(f.between, r.between_ab);

    // Partially-missing A at an otherwise-called site: alleles called = 3,
    // so 3 of A's 6 within pairs are lost.
    let w2 = build(
        SiteCompleteness::WithInvariants,
        6,
        &[(0, vec![[0, -1], [0, 1], [1, 1], [0, 0], [0, 1], [1, 0]])],
    );
    let r2 = pbs(&w2, &pops, 0, 1, 2, PbsEstimator::Yi2010).expect("invariant: valid test fixture");
    assert_eq!((r2.within_a.comparable, r2.within_a.missing), (3, 3));
}

#[test]
fn watterson_reports_no_value_for_empty_or_all_missing_windows() {
    let pops = PopulationView::from_assignment(4, &[0, 0, 0, 0], 1).expect("invariant: valid test fixture");

    // Empty grid window: no sites at all.
    let empty = build(SiteCompleteness::WithInvariants, 4, &[]);
    let t = watterson_theta(&empty, &pops, 0).expect("invariant: valid test fixture");
    assert_eq!(t.n_sites, 0);
    assert_eq!(t.value(), None);
    assert_eq!(t.total_theta(), None);

    // All-missing window: a site is present but no complete-case site is.
    let all_missing = build(
        SiteCompleteness::WithInvariants,
        4,
        &[(0, vec![[-1, -1], [-1, -1], [-1, -1], [-1, -1]])],
    );
    let t2 = watterson_theta(&all_missing, &pops, 0).expect("invariant: valid test fixture");
    assert_eq!((t2.n_sites, t2.n_skipped_missing), (0, 1));
    assert_eq!(t2.value(), None);
    assert_eq!(t2.total_theta(), None);

    // Fully-observed monomorphic window: theta is genuinely zero.
    let mono = build(
        SiteCompleteness::WithInvariants,
        4,
        &[(0, vec![[0, 0], [0, 0], [0, 0], [0, 0]])],
    );
    let t3 = watterson_theta(&mono, &pops, 0).expect("invariant: valid test fixture");
    assert_eq!(t3.value(), Some(0.0));
    assert_eq!(t3.total_theta(), Some(0.0));
}

#[test]
fn watterson_total_theta_matches_tajima_theta_w() {
    let w = build(
        SiteCompleteness::WithInvariants,
        4,
        &[
            (0, vec![[0, 0], [0, 1], [1, 1], [0, 1]]),
            (4, vec![[0, 0], [0, 0], [0, 1], [0, 0]]),
            (8, vec![[0, 0], [0, 0], [0, 0], [0, 0]]),
        ],
    );
    let pops = PopulationView::from_assignment(4, &[0, 0, 0, 0], 1).expect("invariant: valid test fixture");
    let t = watterson_theta(&w, &pops, 0).expect("invariant: valid test fixture");
    // S=2, n=8: total theta = 2 / a1(8) — the quantity Tajima's D compares
    // against π; per-bp divides by the 10 bp window.
    assert!((t.total_theta().expect("invariant: valid test fixture") - 0.7713498622589533).abs() < EPS);
    assert!((t.theta_per_bp().expect("invariant: valid test fixture") - 0.07713498622589533).abs() < EPS);
    assert_eq!(t.value(), t.theta_per_bp());
}

#[test]
fn pbs_and_ld_summary_expose_value() {
    // PBS needs three distinct groups.
    let w = build(
        SiteCompleteness::WithInvariants,
        6,
        &[(0, vec![[0, 0], [0, 1], [1, 1], [0, 0], [0, 1], [1, 1]])],
    );
    let p6 = PopulationView::from_assignment(6, &[0, 0, 1, 1, 2, 2], 3).expect("invariant: valid test fixture");
    let r = pbs(&w, &p6, 0, 1, 2, PbsEstimator::PbsN1).expect("invariant: valid test fixture");
    assert!(r.value().is_some());
    assert_eq!(r.value(), r.pbs);

    let l = ld_r2_summary(&fixture(), &pops4(), 0).expect("invariant: valid test fixture");
    assert_eq!(l.value(), l.mean_r2);
}

#[test]
fn window_level_results_carry_coordinates() {
    let w = fixture();
    let p = pops4();
    let coords = (w.contig, w.start, w.end);

    let a = pi(&w, &p, 0, PiEstimator::PerPairComparable).expect("invariant: valid test fixture");
    assert_eq!((a.contig, a.start, a.end), coords);

    let f2r = f2(&w, &p, 0, 1).expect("invariant: valid test fixture");
    assert_eq!((f2r.contig, f2r.start, f2r.end), coords);

    let l = ld_r2_summary(&w, &p, 0).expect("invariant: valid test fixture");
    assert_eq!((l.contig, l.start, l.end), coords);

    let t = watterson_theta(&w, &p, 0).expect("invariant: valid test fixture");
    assert_eq!((t.contig, t.start, t.end), coords);

    let s = sfs(&w, &p, 0, false).expect("invariant: valid test fixture");
    assert_eq!((s.contig, s.start, s.end), coords);

    let counts = allele_counts(&w, &p).expect("invariant: valid test fixture");
    assert_eq!((counts[0].contig, counts[0].start, counts[0].end), coords);

    // Per-site results carry the contig alongside pos.
    let hw = hwe_exact(&w, &p, 0).expect("invariant: valid test fixture");
    assert!(!hw.is_empty());
    assert!(hw.iter().all(|x| x.contig == w.contig));
    let ht = heterozygosity(&w, &p, 0).expect("invariant: valid test fixture");
    assert!(ht.iter().all(|x| x.contig == w.contig));
    let lp = ld_r2_pairs(&w, &p, 0).expect("invariant: valid test fixture");
    assert!(lp.iter().all(|x| x.contig == w.contig));
}

#[test]
fn availability_gates_follow_sealed_definitions() {
    let site = vec![[0, 0], [0, 1], [1, 1], [0, 1]];
    let unknown = build(
        SiteCompleteness::VariantOnly(Ascertainment::Unknown),
        4,
        &[(0, site.clone())],
    );
    let complete = build(
        SiteCompleteness::VariantOnly(Ascertainment::Complete),
        4,
        &[(0, site.clone())],
    );
    let filtered = build(
        SiteCompleteness::VariantOnly(Ascertainment::Filtered),
        4,
        &[(0, site.clone())],
    );
    let full = fixture();
    let p = pops4();

    // θ_W shares the AccessibleBp denominator: refused under VO(Unknown).
    assert!(watterson_theta(&unknown, &p, 0).is_err());
    assert!(watterson_theta(&filtered, &p, 0).expect("invariant: valid test fixture").bias.is_some());
    assert!(watterson_theta(&complete, &p, 0).expect("invariant: valid test fixture").bias.is_none());
    assert!(watterson_theta(&full, &p, 0).expect("invariant: valid test fixture").bias.is_none());

    // SFS / Tajima's D / LD are invariant to monomorphic-site presence:
    // VO(Complete) is exact; Filtered/Unknown carry the note.
    assert!(sfs(&unknown, &p, 0, false).expect("invariant: valid test fixture").bias.is_some());
    assert!(sfs(&complete, &p, 0, false).expect("invariant: valid test fixture").bias.is_none());
    assert!(tajima_d(&unknown, &p, 0, 1).expect("invariant: valid test fixture").bias.is_some());
    assert!(tajima_d(&complete, &p, 0, 1).expect("invariant: valid test fixture").bias.is_none());
    assert!(ld_r2_summary(&unknown, &p, 0).expect("invariant: valid test fixture").bias.is_some());
    assert!(ld_r2_summary(&complete, &p, 0).expect("invariant: valid test fixture").bias.is_none());

    // f2/f3/f4 are site means that include monomorphic sites under all-sites
    // input: any variant-only input changes the quantity.
    assert!(f2(&unknown, &p, 0, 1).expect("invariant: valid test fixture").bias.is_some());
    assert!(f2(&complete, &p, 0, 1).expect("invariant: valid test fixture").bias.is_some());
    assert!(f2(&full, &p, 0, 1).expect("invariant: valid test fixture").bias.is_none());

    // …but D is a ratio of sums: invariant, like SFS/Tajima/LD. Needs four
    // distinct groups.
    let p4 = PopulationView::from_assignment(4, &[0, 1, 2, 3], 4).expect("invariant: valid test fixture");
    assert!(d_stat(&complete, &p4, 0, 1, 2, 3).expect("invariant: valid test fixture").bias.is_none());
    assert!(d_stat(&unknown, &p4, 0, 1, 2, 3).expect("invariant: valid test fixture").bias.is_some());
    assert!(d_stat(&full, &p4, 0, 1, 2, 3).expect("invariant: valid test fixture").bias.is_none());

    // Site-local (HWE, heterozygosity) and census (allele counts) estimators
    // are exempt from gating: they run under any manifest.
    assert!(hwe_exact(&unknown, &p, 0).is_ok());
    assert!(heterozygosity(&unknown, &p, 0).is_ok());
    assert!(allele_counts(&unknown, &p).is_ok());
}

#[test]
fn population_mismatch_is_an_error_not_a_panic() {
    let w = fixture();
    // A view built for a different sample count.
    let bad = PopulationView::from_assignment(5, &[0, 0, 1, 1, 1], 2).expect("invariant: valid test fixture");
    assert_eq!(
        bad.validate_against(&w),
        Err(PopError::SampleCountMismatch {
            view_samples: 5,
            window_samples: 4
        })
    );
    assert!(pops4().validate_against(&w).is_ok());

    // Estimator entry points surface the mismatch as Unavailable instead of
    // panicking deep inside allele counting (§12.2-8).
    assert!(pi(&w, &bad, 0, PiEstimator::PerPairComparable).is_err());
    assert!(allele_counts(&w, &bad).is_err());
    assert!(hwe_exact(&w, &bad, 0).is_err());
    assert!(ld_r2_summary(&w, &bad, 0).is_err());
    assert!(sfs(&w, &bad, 0, false).is_err());
}

#[test]
fn group_out_of_range_is_an_error() {
    let w = fixture();
    let p = pops4();
    assert!(pi(&w, &p, 7, PiEstimator::PerPairComparable).is_err());
    assert!(f2(&w, &p, 0, 7).is_err());
    assert!(p.get(7).is_none());
    assert_eq!(p.get(1), Some(&[2u32, 3][..]));
}
