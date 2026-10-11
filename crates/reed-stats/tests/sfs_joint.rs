//! Joint (two-population) SFS canonical tests: the complete-case joint
//! spectrum with folding.

use reed_core::count::SiteGroupCounts;
use reed_stats::sfs::JointSfsAcc;

#[test]
fn joint_sfs_counts_and_folds() {
    // 2+2 chromosomes per group; complete-case sites only.
    let mut acc = JointSfsAcc::new(0, 1, 4, 4, false);
    let mk = |a1: u32, a2: u32| {
        vec![
            SiteGroupCounts {
                n_ref_alleles: 4 - a1,
                n_alt_alleles: a1,
                ..Default::default()
            },
            SiteGroupCounts {
                n_ref_alleles: 4 - a2,
                n_alt_alleles: a2,
                ..Default::default()
            },
        ]
    };
    acc.push(&mk(1, 2));
    acc.push(&mk(2, 2));
    acc.push(&mk(0, 0));
    // A site with missing data in group 2 skips (one full missing
    // genotype: ref+alt = 2 < 4 → two missing alleles).
    let skipped = [
        SiteGroupCounts {
            n_ref_alleles: 3,
            n_alt_alleles: 1,
            ..Default::default()
        },
        SiteGroupCounts {
            n_ref_alleles: 2,
            n_alt_alleles: 0,
            n_half_called: 1,
            ..Default::default()
        },
    ];
    acc.push(&skipped);
    assert_eq!(acc.site_tally(), (3, 1));
    assert_eq!(acc.at(1, 2), 1.0);
    assert_eq!(acc.at(2, 2), 1.0);
    assert_eq!(acc.at(0, 0), 1.0);
    assert_eq!(acc.at(1, 1), 0.0);

    // Folded: (3, 2) → (1, 2) for n=4; (0, 0) stays (0, 0).
    let mut folded = JointSfsAcc::new(0, 1, 4, 4, true);
    folded.push(&mk(3, 2));
    folded.push(&mk(0, 0));
    assert!((folded.at(1, 2) - 1.0).abs() < 1e-12);
    assert!((folded.at(0, 0) - 1.0).abs() < 1e-12);
    assert!((folded.at(3, 2) - 0.0).abs() < 1e-12);
}
