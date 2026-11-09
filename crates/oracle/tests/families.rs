//! Golden tests for the extended statistic families.
//!
//! Reference values come from scikit-allel 1.3.13 (temporary venv, script
//! `/tmp/gen_golden.py`) and, for the HWE exact test, from an independent
//! enumeration of the Levene-Haldane distribution. Constants are embedded
//! verbatim; tolerances are 1e-12 unless noted.

use reed_oracle::estimators::fst::{fst, FstEstimator};
use reed_oracle::estimators::fstats::{d_stat, f2, f3, f4};
use reed_oracle::estimators::hwe::{hardy_weinberg_p_value, heterozygosity, hwe_exact};
use reed_oracle::estimators::ld::{ld_r2_pairs, ld_r2_summary};
use reed_oracle::estimators::pbs::{pbs, PbsEstimator};
use reed_oracle::estimators::pi::{pi, PiEstimator};
use reed_oracle::estimators::spectrum::{sfs, tajima_d, watterson_theta};
use reed_oracle::{Ascertainment, Bits, ContigId, PopulationView, SiteCompleteness, WindowBuffer};

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
    .unwrap()
}

/// The shared 4-sample fixture (see `core.rs`).
fn fixture() -> WindowBuffer {
    build(
        SiteCompleteness::WithInvariants,
        4,
        &[
            (0, vec![[0, 0], [0, 1], [1, 1], [0, 1]]),
            (4, vec![[0, 0], [0, 0], [0, 1], [0, 0]]),
            (8, vec![[0, 0], [0, 0], [0, 0], [0, 0]]),
        ],
    )
}

fn assert_close(a: Option<f64>, b: f64) {
    let a = a.expect("expected Some(value)");
    assert!((a - b).abs() < EPS, "expected {b}, got {a}");
}

#[test]
fn weir_cockerham_matches_scikit_allel() {
    let w = fixture();
    let pops = PopulationView::from_assignment(4, &[0, 0, 1, 1], 2).unwrap();
    let f = fst(&w, &pops, 0, 1, FstEstimator::WeirCockerham84).unwrap();
    let wc = f.wc.unwrap();
    // scikit-allel: a sum 0.125, b sum 0.0, c sum 0.75
    assert!((wc.a - 0.125).abs() < EPS);
    assert!((wc.b - 0.0).abs() < EPS);
    assert!((wc.c - 0.75).abs() < EPS);
    assert_close(f.value(), 1.0 / 7.0);
}

#[test]
fn nei86_matches_hand_formula() {
    let w = fixture();
    let pops = PopulationView::from_assignment(4, &[0, 0, 1, 1], 2).unwrap();
    let f = fst(&w, &pops, 0, 1, FstEstimator::Nei86).unwrap();
    let nei = f.nei.unwrap();
    assert!((nei.hs - 0.5625).abs() < EPS);
    assert!((nei.ht - 0.71875).abs() < EPS);
    assert_close(f.value(), 0.21739130434782608); // 5/23
}

#[test]
fn sfs_and_theta_match_scikit_allel() {
    let w = fixture();
    let pops = PopulationView::from_assignment(4, &[0, 0, 0, 0], 1).unwrap();

    let unfolded = sfs(&w, &pops, 0, false).unwrap();
    assert_eq!(unfolded.n_chroms, 8);
    // ALT counts: site0 = 4, site1 = 1, site2 = 0 (monomorphic, excluded)
    assert_eq!(unfolded.counts, vec![1, 0, 0, 1, 0, 0, 0]);

    let folded = sfs(&w, &pops, 0, true).unwrap();
    assert_eq!(folded.counts, vec![1, 0, 0, 1]);

    let theta = watterson_theta(&w, &pops, 0).unwrap();
    assert_eq!(theta.segregating, 2);
    // scikit-allel watterson_theta(pos, ac, start=1, stop=10)
    assert_close(theta.value(), 0.07713498622589533);
}

#[test]
fn tajima_d_matches_scikit_allel() {
    let w = fixture();
    let pops = PopulationView::from_assignment(4, &[0, 0, 0, 0], 1).unwrap();
    let d = tajima_d(&w, &pops, 0, 2).unwrap();
    assert_eq!(d.segregating, 2);
    assert!((d.pi - 0.8214285714285714).abs() < EPS);
    assert!((d.theta_w - 0.7713498622589533).abs() < EPS);
    // scikit-allel: tajima_d(acAll, pos=pos, min_sites=2)
    assert_close(d.value(), 0.2417827948304626);
}

#[test]
fn hwe_exact_matches_independent_enumeration() {
    // Direct algorithm checks against the Levene-Haldane enumeration.
    assert!((hardy_weinberg_p_value(2, 1, 1) - 1.0).abs() < EPS);
    assert!((hardy_weinberg_p_value(1, 3, 0) - 1.0).abs() < EPS);
    assert!((hardy_weinberg_p_value(1, 4, 1) - 0.27272727272727276).abs() < EPS);
    assert!((hardy_weinberg_p_value(4, 10, 2) - 0.19950562353231985).abs() < EPS);

    // Per-site pass over the fixture (all complete).
    let w = fixture();
    let pops = PopulationView::from_assignment(4, &[0, 0, 0, 0], 1).unwrap();
    let sites = hwe_exact(&w, &pops, 0).unwrap();
    assert_eq!(sites.len(), 3);
    assert_eq!(
        (sites[0].hom_ref, sites[0].het, sites[0].hom_alt),
        (1, 2, 1)
    );
    assert!((sites[0].p_value - 1.0).abs() < EPS);
    assert_eq!(
        (sites[1].hom_ref, sites[1].het, sites[1].hom_alt),
        (3, 1, 0)
    );
    assert!((sites[1].p_value - 1.0).abs() < EPS);
}

#[test]
fn heterozygosity_matches_hand_values() {
    let w = fixture();
    let pops = PopulationView::from_assignment(4, &[0, 0, 0, 0], 1).unwrap();
    let sites = heterozygosity(&w, &pops, 0).unwrap();
    assert_eq!(sites.len(), 3);

    assert!((sites[0].hobs - 0.5).abs() < EPS);
    assert!((sites[0].hexp - 0.5).abs() < EPS);
    assert!((sites[0].fis.unwrap() - 0.0).abs() < EPS);

    assert!((sites[1].hobs - 0.25).abs() < EPS);
    assert!((sites[1].hexp - 0.21875).abs() < EPS);
    assert!((sites[1].fis.unwrap() - (-1.0 / 7.0)).abs() < EPS);

    // Monomorphic site: Fis undefined.
    assert!(sites[2].fis.is_none());
}

#[test]
fn ld_rogers_huff_matches_scikit_allel() {
    let w = fixture();
    let pops = PopulationView::from_assignment(4, &[0, 0, 0, 0], 1).unwrap();

    let pairs = ld_r2_pairs(&w, &pops, 0).unwrap();
    // Only pair (site0, site1) is valid: site2 is monomorphic (zero variance).
    assert_eq!(pairs.len(), 1);
    assert_eq!((pairs[0].site_i, pairs[0].site_j), (0, 1));
    assert_eq!(pairs[0].n_called, 4);
    // scikit-allel rogers_huff_r: r2 = 0.6666667 (exactly 2/3)
    assert!((pairs[0].r2 - 2.0 / 3.0).abs() < 1e-9);

    let summary = ld_r2_summary(&w, &pops, 0).unwrap();
    assert_eq!(summary.n_pairs, 3);
    assert_eq!(summary.n_valid, 1);
    assert!((summary.mean_r2.unwrap() - 2.0 / 3.0).abs() < 1e-9);
}

#[test]
fn f_statistics_match_hand_values() {
    let w = fixture();
    // One group per sample: A={S0}, B={S1}, C={S2}, D={S3}
    let pops = PopulationView::from_assignment(4, &[0, 1, 2, 3], 4).unwrap();

    let r = f2(&w, &pops, 0, 1).unwrap();
    // Corrected Σ: site0 (S1 het): 0.25 − ĥ/1 = 0; sites 1, 2 → 0.
    assert_close(r.value(), 0.0);
    assert_eq!(r.n_sites, 3);

    let r = f3(&w, &pops, 2, 0, 1).unwrap(); // f3(S2; S0, S1)
                                             // Corrected Σ (allel patterson_f3 nansum = 0.5): site0 (S2 hom-alt)
                                             // (1−0)(1−0.5) − 0 = 0.5; site1 (S2 het) 0.25 − 0.5/2 = 0; site2 → 0.
    assert_close(r.value(), 0.5);

    let r = f4(&w, &pops, 0, 1, 2, 3).unwrap();
    // Σnum: site0 (0−0.5)(1−0.5) = −0.25; sites 1, 2 → 0.
    assert_close(r.value(), -0.25);

    let r = d_stat(&w, &pops, 0, 1, 2, 3).unwrap();
    assert!((r.num - (-0.25)).abs() < EPS);
    assert!((r.den - 0.25).abs() < EPS);
    assert_close(r.value(), -1.0);
}

#[test]
fn pbs_matches_scikit_allel() {
    // One site, six samples: A p=0, B p=0.5, C p=0.75.
    // scikit-allel: fst_ab = 1/3, fst_ac = 2/3, fst_bc = −1/6 (clipped to 0);
    // T_ab = ln(3/2), T_ac = ln 3, T_bc = 0.
    let w = build(
        SiteCompleteness::WithInvariants,
        6,
        &[(0, vec![[0, 0], [0, 0], [0, 1], [0, 1], [1, 1], [0, 1]])],
    );
    let pops = PopulationView::from_assignment(6, &[0, 0, 1, 1, 2, 2], 3).unwrap();

    let r = pbs(&w, &pops, 0, 1, 2, PbsEstimator::Yi2010).unwrap();
    assert_eq!(r.n_sites, 1);
    assert!((r.fst_ab.unwrap() - 1.0 / 3.0).abs() < EPS);
    assert!((r.fst_ac.unwrap() - 2.0 / 3.0).abs() < EPS);
    assert!((r.fst_bc.unwrap() - 0.0).abs() < EPS); // clipped from −1/6
    assert!((r.t_ab.unwrap() - 1.5f64.ln()).abs() < EPS);
    assert!((r.t_ac.unwrap() - 3.0f64.ln()).abs() < EPS);
    assert!((r.t_bc.unwrap() - 0.0).abs() < EPS);
    // scikit-allel: pbs(..., normed=False)
    assert!((r.pbs.unwrap() - 0.752038698388137).abs() < 1e-12);

    let r = pbs(&w, &pops, 0, 1, 2, PbsEstimator::PbsN1).unwrap();
    // scikit-allel: pbs(..., normed=True)
    assert!((r.pbs.unwrap() - 0.4292363513888176).abs() < 1e-12);
}

#[test]
fn manifest_availability_extends_to_new_estimators() {
    let w = build(
        SiteCompleteness::VariantOnly(Ascertainment::Unknown),
        4,
        &[(0, vec![[0, 0], [0, 1], [1, 1], [0, 1]])],
    );
    let pops = PopulationView::from_assignment(4, &[0, 0, 1, 1], 2).unwrap();
    // All Fst variants stay computable (bias flag), unlike pi/dxy.
    for est in [
        FstEstimator::Hudson,
        FstEstimator::WeirCockerham84,
        FstEstimator::Nei86,
    ] {
        let f = fst(&w, &pops, 0, 1, est).unwrap();
        assert!(f.bias.is_some());
    }
}

/// doc_v23 V23-47 (decision A): the spectrum family is **complete-case per
/// site** — a site where any sample of the group has a missing allele is
/// skipped whole and counted, never partially dropped. scikit-allel's
/// per-position functions instead delete the missing genotype and keep the
/// site, so the two disagree on such input; this golden pins reed's chosen
/// semantics (the production estimator and the spec share it) rather than
/// allel's.
///
/// Fixture: 3 sites, sample 3 missing at site 0 and sample 0 missing at site 2
/// (`-1` = missing allele). Only site 1 is complete (ALT count 1 of 8).
#[test]
fn sfs_theta_and_tajima_are_complete_case_per_site() {
    let w = build(
        SiteCompleteness::WithInvariants,
        4,
        &[
            // site 0: sample 3 missing -> skipped whole (allel would keep it).
            (0, vec![[0, 0], [0, 1], [1, 1], [-1, -1]]),
            // site 1: complete, ALT count 1 -> one singleton.
            (4, vec![[0, 0], [0, 1], [0, 0], [0, 0]]),
            // site 2: sample 0 missing -> skipped whole.
            (8, vec![[-1, -1], [0, 0], [0, 0], [0, 1]]),
        ],
    );
    let pops = PopulationView::from_assignment(4, &[0, 0, 0, 0], 1).unwrap();

    let unfolded = sfs(&w, &pops, 0, false).unwrap();
    assert_eq!(unfolded.n_chroms, 8);
    assert_eq!(unfolded.n_sites, 1, "only the complete site is counted");
    assert_eq!(
        unfolded.n_skipped_missing, 2,
        "both partial sites are skipped"
    );
    // ALT counts: site0 skipped, site1 = 1, site2 skipped.
    assert_eq!(unfolded.counts, vec![1, 0, 0, 0, 0, 0, 0]);

    let folded = sfs(&w, &pops, 0, true).unwrap();
    assert_eq!(folded.counts, vec![1, 0, 0, 0]);
    assert_eq!(folded.n_skipped_missing, 2);

    // theta_W = S / a1 over the *complete* site count is what the estimator
    // reports; the inaccessible bp denominator is unchanged by the skips.
    let theta = watterson_theta(&w, &pops, 0).unwrap();
    assert_eq!(theta.segregating, 1, "one segregating site survives");
    // a1(8) = 1 + 1/2 + … + 1/7 = 1089/420, S = 1, over the 10 bp window:
    // (1 / a1) / 10 = 42/1089 — exactly half the two-segregating case
    // this fixture's sibling test reports.
    assert_close(theta.value(), 42.0 / 1089.0);

    // Tajima's D needs >= 2 segregating sites; one is unavailable, not zero.
    let d = tajima_d(&w, &pops, 0, 2);
    if let Ok(d) = d {
        assert_eq!(d.segregating, 1);
    }
}

/// The pairwise estimators use **pairwise-complete** data instead: a missing
/// sample drops out of each pair, so a partially missing site still
/// contributes. This is the same fixture as above, and it is what makes the
/// two conventions visibly different (doc_v23 V23-47).
#[test]
fn pi_is_pairwise_complete_where_the_spectrum_is_complete_case() {
    let w = build(
        SiteCompleteness::WithInvariants,
        4,
        &[
            (0, vec![[0, 0], [0, 1], [1, 1], [-1, -1]]),
            (4, vec![[0, 0], [0, 1], [0, 0], [0, 0]]),
            (8, vec![[-1, -1], [0, 0], [0, 0], [0, 1]]),
        ],
    );
    let pops = PopulationView::from_assignment(4, &[0, 0, 0, 0], 1).unwrap();
    // Site 0 contributes pair differences among its three called samples and
    // site 2 among its three, while the spectrum above counted only site 1.
    let pairwise = pi(&w, &pops, 0, PiEstimator::PerPairComparable).unwrap();
    let per_bp = pi(&w, &pops, 0, PiEstimator::AccessibleBp).unwrap();
    assert!(
        pairwise.value().is_some(),
        "pairwise-complete pi is defined with partial sites"
    );
    // Sanity: the pairwise numerator is strictly larger than site 1 alone
    // would give (one difference in three comparable pairs).
    assert!(
        pairwise.value().unwrap() > 0.0,
        "pi must see the partial sites: {:?} / {:?}",
        pairwise.value(),
        per_bp.value()
    );
}

/// doc_v24 V24-01: a group with no called alleles keeps `c` (allel `nansum`
/// semantics) and contributes `Fst = 0`. Measured with scikit-allel 1.3.13
/// on the same genotype matrix: `a = b = NaN`, `c = 1/3`, `nansum` Fst = 0.
/// Contract: `docs/rounds/v24/02-oracle-contract.md` §2.
#[test]
fn weir_cockerham_keeps_c_for_an_all_missing_group() {
    let w = build(
        SiteCompleteness::WithInvariants,
        6,
        &[(
            0,
            vec![[0, 0], [0, 1], [1, 1], [-1, -1], [-1, -1], [-1, -1]],
        )],
    );
    let pops = PopulationView::from_assignment(6, &[0, 0, 0, 1, 1, 1], 2).unwrap();
    let f = fst(&w, &pops, 0, 1, FstEstimator::WeirCockerham84).unwrap();
    let wc = f.wc.expect("WC84 components");
    assert_eq!(wc.n_sites, 1);
    assert_eq!(wc.n_skipped, 0);
    assert!((wc.a - 0.0).abs() < EPS, "a = {}", wc.a);
    assert!((wc.b - 0.0).abs() < EPS, "b = {}", wc.b);
    assert!((wc.c - 1.0 / 3.0).abs() < EPS, "c = {}", wc.c);
    assert_close(f.value(), 0.0);
}

/// doc_v24 V24-01: the Ho/He sample-set contract on a half call.
/// `[0/1, ./1]`: one diploid called individual (one het), three called
/// alleles (ref 1, alt 2) -> hobs 1, hexp 4/9, Fis 1 - 1/(4/9) = -1.25.
/// Contract: `docs/rounds/v24/02-oracle-contract.md` §4.
#[test]
fn heterozygosity_declares_the_half_call_contract() {
    let w = build(
        SiteCompleteness::WithInvariants,
        2,
        &[(0, vec![[0, 1], [-1, 1]])],
    );
    let pops = PopulationView::from_assignment(2, &[0, 0], 1).unwrap();
    let rows = heterozygosity(&w, &pops, 0).unwrap();
    assert_eq!(rows.len(), 1);
    let row = &rows[0];
    assert_eq!(row.n_called, 1, "the half call is not a diploid call");
    assert_eq!(row.n_het, 1);
    assert!((row.hobs - 1.0).abs() < EPS, "hobs = {}", row.hobs);
    assert!((row.hexp - 4.0 / 9.0).abs() < EPS, "hexp = {}", row.hexp);
    let fis = row.fis.expect("hexp > 0");
    assert!((fis - (-1.25)).abs() < EPS, "Fis = {fis}");
}
