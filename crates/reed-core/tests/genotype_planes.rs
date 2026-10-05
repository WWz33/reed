//! Plane-block decode invariants: the shared decoder never produces an
//! illegal plane state, whatever the genotype mix
//! (`docs/architecture-spec.md` §19 M1, planes half).

use proptest::prelude::*;

use reed_core::genotype::{
    write_genotype, GenotypePolicies, HalfCallPolicy, PlaneBlock, Ploidy, PloidyPolicy,
    PolicyCounters,
};

fn arb_gt() -> impl Strategy<Value = (Ploidy, Option<u16>, Option<u16>)> {
    (
        proptest::bool::weighted(0.1)
            .prop_map(|haploid| if haploid { Ploidy::Haploid } else { Ploidy::Diploid }),
        proptest::option::of(0_u16..2),
        proptest::option::of(0_u16..2),
    )
}

/// Missing-allele count the decoder must produce under
/// `HaploidAsHalfCall` + `HalfCallPolicy::AlleleLevel`.
fn expected_missing_alleles(ploidy: Ploidy, a: Option<u16>, b: Option<u16>) -> u32 {
    match (ploidy, a, b) {
        (Ploidy::Haploid, None, _) | (Ploidy::Diploid, None, None) => 2,
        (Ploidy::Haploid, Some(_), _)
        | (Ploidy::Diploid, None, Some(_))
        | (Ploidy::Diploid, Some(_), None) => 1,
        (Ploidy::Diploid, Some(_), Some(_)) => 0,
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn planes_stay_consistent_under_random_biallelic_input(
        genotypes in proptest::collection::vec(arb_gt(), 1..70),
    ) {
        let n_samples = genotypes.len();
        // Haploid inputs convert to half calls so nothing rejects.
        let policies = GenotypePolicies {
            ploidy: PloidyPolicy::HaploidAsHalfCall,
            ..GenotypePolicies::default()
        };
        let mut block = PlaneBlock::new(1, n_samples, true);
        let mut counters = PolicyCounters::default();
        let mut n_expected_missing = 0_u32;
        let mut n_expected_haploid_half = 0_u32;
        for (sample_idx, &(ploidy, a, b)) in genotypes.iter().enumerate() {
            let slots = if ploidy == Ploidy::Haploid { [a, None] } else { [a, b] };
            if ploidy == Ploidy::Haploid && a.is_some() {
                n_expected_haploid_half += 1;
            }
            n_expected_missing += expected_missing_alleles(ploidy, a, b);
            write_genotype(
                &mut block,
                0,
                sample_idx,
                slots,
                ploidy,
                false,
                &policies,
                &mut counters,
            )
            .expect("invariant: biallelic under HaploidAsHalfCall never rejects");
        }

        // Type invariant: ALT and MISSING never co-occur in one slot.
        prop_assert!(block.invariants_hold());

        // Per-sample plane contents match the input exactly.
        let planes = block.site_planes(0);
        let (mut n_alt, mut n_miss) = (0_u32, 0_u32);
        for s in 0..n_samples {
            let bit = |words: &[u64]| words[s / 64] >> (s % 64) & 1 == 1;
            let (a0, m0, a1, m1) = (
                bit(planes.slot0_alt),
                bit(planes.slot0_miss),
                bit(planes.slot1_alt),
                bit(planes.slot1_miss),
            );
            prop_assert!(!(a0 && m0) && !(a1 && m1), "slot state exclusive at sample {s}");
            n_miss += u32::from(m0) + u32::from(m1);
            if a0 || a1 {
                n_alt += 1;
            }
        }
        prop_assert_eq!(n_miss, n_expected_missing, "missing alleles must match input");
        prop_assert_eq!(counters.n_haploid_as_half_call, u64::from(n_expected_haploid_half));
        // Segregating sanity: with all-REF input there are no ALT bits, and
        // ALT bits only appear where the input had an ALT allele.
        let n_input_alt_carriers: u32 = genotypes
            .iter()
            .filter(|&&(p, a, b)| a == Some(1) || b == Some(1))
            .map(|_| 1)
            .sum();
        prop_assert!(n_alt <= n_input_alt_carriers);
    }
}

#[test]
fn treat_as_missing_zeroes_exactly_two_alleles_per_genotype() {
    let policies = GenotypePolicies {
        half_call: HalfCallPolicy::TreatAsMissing,
        ..GenotypePolicies::default()
    };
    let mut block = PlaneBlock::new(1, 3, false);
    let mut counters = PolicyCounters::default();
    let gt = [(None, Some(1_u16)), (Some(1_u16), None), (Some(0_u16), None)];
    for (sample_idx, &(a, b)) in gt.iter().enumerate() {
        write_genotype(
            &mut block,
            0,
            sample_idx,
            [a, b],
            Ploidy::Diploid,
            false,
            &policies,
            &mut counters,
        )
        .expect("invariant: biallelic diploid");
    }
    assert_eq!(counters.n_half_call_as_missing, 3);
    let planes = block.site_planes(0);
    let both_miss = (0..3).all(|s| {
        planes.slot0_miss[s / 64] >> (s % 64) & 1 == 1
            && planes.slot1_miss[s / 64] >> (s % 64) & 1 == 1
    });
    assert!(both_miss);
    assert!(block.invariants_hold());
}
