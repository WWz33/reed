//! Individual F (`--het`) canonical tests: hand-computed plink `--het`
//! moments over plane blocks.

use reed_core::genotype::{GenotypePolicies, PlaneBlock, PolicyCounters, write_genotype};
use reed_core::sample::SampleTable;
use reed_stats::het::IndividualFAcc;

/// Builds a PlaneBlock from explicit diploid genotypes
/// (`(slot0, slot1)`, `None` = missing).
fn block_from(genotypes: &[[(Option<u16>, Option<u16>); 4]]) -> PlaneBlock {
    let mut block = PlaneBlock::new(genotypes.len(), 4, false);
    let policies = GenotypePolicies::default();
    let mut counters = PolicyCounters::default();
    for (site_idx, site) in genotypes.iter().enumerate() {
        for (sample_idx, &(a, b)) in site.iter().enumerate() {
            write_genotype(
                &mut block,
                site_idx,
                sample_idx,
                [a, b],
                reed_core::genotype::Ploidy::Diploid,
                false,
                &policies,
                &mut counters,
            )
            .expect("invariant: test genotypes are biallelic diploid");
        }
    }
    block
}

#[test]
fn individual_f_matches_plink_moment() {
    // Site A: p_alt = 2/8 = 0.25 → E_het share 2p(1−p) = 0.375.
    // Site B: p_alt = 4/8 = 0.5 → 0.5.
    // Genotypes (s0, s1, s2, s3):
    //   A: 0/0, 0/1, 1/1, 1/1
    //   B: 0/1, 1/1, 0/1, 0/1
    // s0: O_hom = 1 (A 0/0), N = 2, Σ2p(1−p) = 0.375 + 0.5
    //   F = (1 − (2 − 0.875)) / (2 − (2 − 0.875)) = (1 − 1.125)/0.875
    // s1: O_hom = 1 (B 1/1), N = 2, same Σ → F = −0.125/0.875
    // s2: O_hom = 1 (A 1/1), N = 2 → F = −0.125/0.875
    // s3: O_hom = 2 (A 1/1, B 0/1? no: B 0/1 is het) → O_hom = 1...
    let block = block_from(&[
        [
            (Some(0), Some(0)),
            (Some(0), Some(1)),
            (Some(1), Some(1)),
            (Some(1), Some(1)),
        ],
        [
            (Some(0), Some(1)),
            (Some(1), Some(1)),
            (Some(0), Some(1)),
            (Some(0), Some(1)),
        ],
    ]);
    let samples = SampleTable::from_names(["s0", "s1", "s2", "s3"]).unwrap();
    let groups = reed_core::group::GroupTable::from_assignment(
        &samples,
        &[("s0", "P"), ("s1", "P"), ("s2", "P"), ("s3", "P")],
    )
    .unwrap();
    let mask = groups.masks()[0].clone();

    let mut acc = IndividualFAcc::new(4);
    // Site A: alt = 2+1+2+2 = 7? Recount: 0/0 → 0; 0/1 → 1; 1/1 → 2; 1/1 → 2.
    // Total alt = 5, AN = 8 → p = 0.625. Recompute the expectation below.
    acc.push_site(&block, 0, &mask, 5.0 / 8.0);
    acc.push_site(&block, 1, &mask, 4.0 / 8.0);

    let f = acc.f_values();
    // s0: O_hom = 2 (A 0/0, B het? B s0 is 0/1 het) → O_hom = 1. N = 2.
    // Σ2p(1−p) = 2(0.625·0.375) + 2(0.5·0.5) = 0.46875 + 0.5 = 0.96875.
    // E_hom = 2 − 0.96875 = 1.03125. F = (1 − 1.03125)/(2 − 1.03125).
    let e0 = 2.0 - (2.0 * 0.625 * 0.375 + 2.0 * 0.5 * 0.5);
    let expected_s0 = (1.0 - e0) / (2.0 - e0);
    assert!(
        (f[0].expect("invariant: denom") - expected_s0).abs() < 1e-12,
        "{f:?}"
    );
    // s1: A 0/1 het, B 1/1 hom → O_hom = 1, N = 2 → same F.
    assert!((f[1].expect("invariant") - expected_s0).abs() < 1e-12);
    // s2: A 1/1 hom, B 0/1 het → O_hom = 1 → same F.
    assert!((f[2].expect("invariant") - expected_s0).abs() < 1e-12);
    // s3: A 1/1 hom, B 0/1 het → O_hom = 1 → same F.
    assert!((f[3].expect("invariant") - expected_s0).abs() < 1e-12);
}

#[test]
fn individual_f_merge_is_additive() {
    let block = block_from(&[
        [
            (Some(0), Some(0)),
            (Some(0), Some(1)),
            (Some(1), Some(1)),
            (Some(1), Some(1)),
        ],
        [
            (Some(0), Some(1)),
            (Some(1), Some(1)),
            (Some(0), Some(1)),
            (Some(0), Some(1)),
        ],
    ]);
    let samples = SampleTable::from_names(["s0", "s1", "s2", "s3"]).unwrap();
    let groups = reed_core::group::GroupTable::from_assignment(
        &samples,
        &[("s0", "P"), ("s1", "P"), ("s2", "P"), ("s3", "P")],
    )
    .unwrap();
    let mask = groups.masks()[0].clone();

    let mut whole = IndividualFAcc::new(4);
    whole.push_site(&block, 0, &mask, 5.0 / 8.0);
    whole.push_site(&block, 1, &mask, 0.5);

    let mut a = IndividualFAcc::new(4);
    let mut b = IndividualFAcc::new(4);
    a.push_site(&block, 0, &mask, 5.0 / 8.0);
    b.push_site(&block, 1, &mask, 0.5);
    a.merge(&b);

    assert_eq!(a.o_hom, whole.o_hom);
    assert_eq!(a.n_called, whole.n_called);
    for (x, y) in a.f_values().iter().zip(whole.f_values()) {
        assert_eq!(x, &y);
    }
}
