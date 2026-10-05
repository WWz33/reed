//! M1 acceptance: the bit-plane counting kernel must reproduce the legacy
//! `WindowBuffer::allele_counts` / `genotype_counts` exactly, and the
//! dispatched (possibly AVX-512) kernel must equal the scalar canonical path
//! (`docs/architecture-spec.md` §19 M1).

use proptest::prelude::*;

use reed_core::bits::Bits;
use reed_core::count::{count_site_group, count_site_group_scalar, CountTable};
use reed_core::genotype::{write_genotype, GenotypePolicies, Ploidy, PolicyCounters};
use reed_core::id::SampleId;
use reed_core::manifest::{Ascertainment, SiteCompleteness};
use reed_core::sample::SampleMask;
use reed_core::window::{ContigId, WindowBuffer};

/// One site's genotypes: `(slot0, slot1)` per sample, `None` = missing.
type SiteGt = Vec<(Option<u16>, Option<u16>)>;

fn arb_genotypes(n_samples: usize) -> impl Strategy<Value = SiteGt> {
    proptest::collection::vec(
        (proptest::option::of(0_u16..2), proptest::option::of(0_u16..2)),
        n_samples,
    )
}

/// Legacy window buffer carrying the same calls.
fn make_legacy(per_site: &[SiteGt]) -> WindowBuffer {
    let n_samples = per_site[0].len();
    let mut pos = Vec::new();
    let mut gt: Vec<i8> = Vec::new();
    let mut missing: Vec<bool> = Vec::new();
    let mut phased: Vec<bool> = Vec::new();
    let mut variant: Vec<bool> = Vec::new();
    for (site_idx, site) in per_site.iter().enumerate() {
        pos.push(100 + 2 * site_idx as u32);
        variant.push(site.iter().any(|&(a, b)| a.unwrap_or(0) > 0 || b.unwrap_or(0) > 0));
        for &(a, b) in site {
            let (ai, bi) = (a.map_or(-1, |v| v as i8), b.map_or(-1, |v| v as i8));
            gt.push(ai);
            gt.push(bi);
            missing.push(ai == -1 || bi == -1);
            phased.push(false);
        }
    }
    WindowBuffer::new(
        SiteCompleteness::VariantOnly(Ascertainment::Complete),
        ContigId(0),
        20_000,
        0,
        10_000,
        n_samples as u32,
        pos,
        Bits::from_bools(&variant),
        gt,
        Bits::from_bools(&missing),
        Bits::from_bools(&phased),
    )
    .expect("invariant: test buffer is valid")
}

/// Plane block carrying the same calls (diploid decode, default policies).
fn make_block(per_site: &[SiteGt]) -> reed_core::genotype::PlaneBlock {
    let mut block = reed_core::genotype::PlaneBlock::new(per_site.len(), per_site[0].len(), false);
    let policies = GenotypePolicies::default();
    let mut counters = PolicyCounters::default();
    for (site_idx, site) in per_site.iter().enumerate() {
        for (sample_idx, &(a, b)) in site.iter().enumerate() {
            write_genotype(
                &mut block,
                site_idx,
                sample_idx,
                [a, b],
                Ploidy::Diploid,
                false,
                &policies,
                &mut counters,
            )
            .expect("invariant: test genotypes are biallelic diploid");
        }
    }
    block
}

fn assert_counts_match(
    legacy: &WindowBuffer,
    block: &reed_core::genotype::PlaneBlock,
    site: usize,
    samples: &[u32],
) {
    let site_planes = block.site_planes(site);
    let n_group = samples.len() as u32;
    // Reference words for the group, built by hand (no SampleMask needed to
    // state the invariant being tested — but we use one, as production does).
    let mask = SampleMask::from_ids(
        block.n_samples(),
        samples.iter().map(|&s| SampleId(s)),
    )
    .expect("invariant: ids in range");

    let mut counts = reed_core::count::SiteGroupCounts::default();
    count_site_group(&site_planes, mask.words(), n_group, &mut counts);

    let allele = legacy.allele_counts(site, samples);
    assert_eq!(counts.n_ref_alleles, allele.ref_count as u32, "site {site} ref");
    assert_eq!(counts.n_alt_alleles, allele.alt_count as u32, "site {site} alt");
    assert_eq!(counts.n_missing_alleles(n_group), allele.missing as u32, "site {site} missing");

    let gt = legacy.genotype_counts(site, samples);
    assert_eq!(counts.n_hom_ref, gt.hom_ref as u32, "site {site} hom_ref");
    assert_eq!(counts.n_het, gt.het as u32, "site {site} het");
    assert_eq!(counts.n_hom_alt, gt.hom_alt as u32, "site {site} hom_alt");
    assert_eq!(counts.n_full_called(), gt.n_called as u32, "site {site} n_called");

    // Scalar canonical path must agree with the dispatched one (SIMD ≡ scalar).
    let mut scalar = reed_core::count::SiteGroupCounts::default();
    count_site_group_scalar(&site_planes, mask.words(), n_group, &mut scalar);
    assert_eq!(counts, scalar, "site {site} dispatched vs scalar");
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn kernel_matches_legacy_on_random_genotypes(
        sites in proptest::collection::vec(arb_genotypes(17), 1..4),
        group_a in prop::collection::vec(any::<bool>(), 17),
        group_b in prop::collection::vec(any::<bool>(), 17),
    ) {
        let per_site: Vec<SiteGt> = sites;
        let legacy = make_legacy(&per_site);
        let block = make_block(&per_site);
        let pick = |flags: &[bool]| -> Vec<u32> {
            flags.iter().enumerate().filter(|&(_, f)| *f).map(|(i, _)| i as u32).collect()
        };
        for site in 0..per_site.len() {
            assert_counts_match(&legacy, &block, site, &pick(&group_a));
            assert_counts_match(&legacy, &block, site, &pick(&group_b));
            // Plus the full sample set (includes every edge of the word mask).
            let all: Vec<u32> = (0..17).collect();
            assert_counts_match(&legacy, &block, site, &all);
        }
    }

    #[test]
    fn kernel_matches_legacy_across_word_widths(
        sites in proptest::collection::vec(arb_genotypes(200), 1..3),
    ) {
        // 200 samples → 4 words: exercises multi-word masks and the
        // vectorised prefix plus scalar tail together.
        let per_site: Vec<SiteGt> = sites;
        let legacy = make_legacy(&per_site);
        let block = make_block(&per_site);
        let samples: Vec<u32> = (0..200).step_by(3).collect();
        for site in 0..per_site.len() {
            assert_counts_match(&legacy, &block, site, &samples);
        }
    }
}

#[test]
fn count_table_fill_matches_legacy_over_groups() {
    // 3 sites × 6 samples, two groups; CountTable::fill must agree with the
    // legacy per-group loop.
    let mut rng_seed = 0x9E3779B9_u64;
    let mut next = move || {
        rng_seed ^= rng_seed << 13;
        rng_seed ^= rng_seed >> 7;
        rng_seed ^= rng_seed << 17;
        rng_seed
    };
    let per_site: Vec<SiteGt> = (0..3)
        .map(|_| {
            (0..6)
                .map(|_| {
                    let r = next();
                    let slot = |bits: u64| match bits % 3 {
                        0 => Some(0_u16),
                        1 => Some(1_u16),
                        _ => None,
                    };
                    (slot(r & 0xF), slot((r >> 4) & 0xF))
                })
                .collect()
        })
        .collect();
    let legacy = make_legacy(&per_site);
    let block = make_block(&per_site);

    let groups: Vec<Vec<u32>> = vec![vec![0, 1, 2], vec![3, 4, 5], (0..6).collect()];
    let masks: Vec<SampleMask> = groups
        .iter()
        .map(|g| SampleMask::from_ids(6, g.iter().copied().map(SampleId)).unwrap())
        .collect();
    let mut table = CountTable::new(3, groups.len());
    table.fill(0..3, &block, &masks);

    for (group_idx, samples) in groups.iter().enumerate() {
        for site in 0..3 {
            let allele = legacy.allele_counts(site, samples);
            let gt = legacy.genotype_counts(site, samples);
            let counts = table.get(site, group_idx).unwrap();
            assert_eq!(counts.n_ref_alleles, allele.ref_count as u32);
            assert_eq!(counts.n_alt_alleles, allele.alt_count as u32);
            assert_eq!(counts.n_hom_ref, gt.hom_ref as u32);
            assert_eq!(counts.n_het, gt.het as u32);
            assert_eq!(counts.n_hom_alt, gt.hom_alt as u32);
        }
    }
}
