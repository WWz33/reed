//! Oracle parity (M4 acceptance gate): the reed-stats accumulators against
//! the frozen v0.2 estimator loop on random windows — integer components
//! exact, floats within 1e-12 (`docs/architecture-spec.md` milestone M4).

use reed_core::count::SiteGroupCounts;
use reed_oracle::estimators::dxy::{dxy, DxyEstimator};
use reed_oracle::estimators::fst::{fst, FstEstimator};
use reed_oracle::estimators::pbs::{pbs, PbsEstimator};
use reed_oracle::estimators::pi::{pi, PiEstimator};
use reed_oracle::estimators::spectrum::{sfs, tajima_d, watterson_theta};
use reed_oracle::{Bits, ContigId, PopulationView, SiteCompleteness, WindowBuffer};
use reed_stats::fst::FstKind;
use reed_stats::pbs::PbsKind;
use reed_stats::request::{build_accs, StatAcc, StatKind};

/// Floats agree within 1e-12 (absolute or relative, whichever is larger).
fn assert_close(actual: Option<f64>, expected: Option<f64>, what: &str) {
    match (actual, expected) {
        (None, None) => {}
        (Some(x), Some(y)) => {
            let tol = 1e-12 * y.abs().max(1.0);
            assert!(
                (x - y).abs() <= tol,
                "{what}: got {x}, expected {y} (tol {tol})"
            );
        }
        _ => panic!("{what}: definedness mismatch: {actual:?} vs {expected:?}"),
    }
}

/// 6 diploid samples, groups A = {0,1}, B = {2,3}, C = {4,5}.
const N_SAMPLES: u32 = 6;
const WINDOW: u64 = 200;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        // xorshift64*; determinism is the point, cryptography is not.
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// One allele: 10% missing, else ref/alt split.
    fn allele(&mut self) -> i8 {
        match self.next() % 20 {
            0 | 1 => -1,
            2..=10 => 0,
            _ => 1,
        }
    }
}

fn hand_counts(gts: &[[i8; 2]]) -> SiteGroupCounts {
    let mut c = SiteGroupCounts::default();
    for &[a, b] in gts {
        match (a, b) {
            (-1, -1) => {}
            (-1, x) | (x, -1) => {
                c.n_half_called += 1;
                // The one called allele still counts (oracle allele_counts).
                if x == 0 {
                    c.n_ref_alleles += 1;
                } else {
                    c.n_alt_alleles += 1;
                }
            }
            (0, 0) => {
                c.n_hom_ref += 1;
                c.n_ref_alleles += 2;
            }
            (0, 1) | (1, 0) => {
                c.n_het += 1;
                c.n_ref_alleles += 1;
                c.n_alt_alleles += 1;
            }
            (1, 1) => {
                c.n_hom_alt += 1;
                c.n_alt_alleles += 2;
            }
            _ => unreachable!("alleles are -1/0/1"),
        }
    }
    c
}

/// Builds one random window and returns (oracle buffer, count rows,
/// site-wide variant flags).
fn random_window(seed: u64, n_sites: usize) -> (WindowBuffer, Vec<Vec<SiteGroupCounts>>) {
    let mut rng = Rng(seed);
    let mut pos = Vec::with_capacity(n_sites);
    let mut is_variant = Vec::with_capacity(n_sites);
    let mut genotypes = Vec::with_capacity(n_sites * N_SAMPLES as usize * 2);
    let mut missing = Vec::with_capacity(n_sites * N_SAMPLES as usize);
    let mut phased = Vec::with_capacity(n_sites * N_SAMPLES as usize);
    let mut rows = Vec::with_capacity(n_sites);
    let mut site_pos = 0_u32;
    for _ in 0..n_sites {
        site_pos += u32::try_from(rng.next() % 12 + 1).expect("static range");
        pos.push(site_pos);
        let mut site_variant = false;
        let mut row = Vec::with_capacity(3);
        for _group in 0..3 {
            let mut gts = [[0_i8; 2]; 2];
            for gt in gts.iter_mut() {
                let alleles = [rng.allele(), rng.allele()];
                for &allele in &alleles {
                    genotypes.push(allele);
                    if allele > 0 {
                        site_variant = true;
                    }
                }
                missing.push(alleles[0] == -1 || alleles[1] == -1);
                phased.push(false);
                *gt = alleles;
            }
            row.push(hand_counts(&gts));
        }
        is_variant.push(site_variant);
        rows.push(row);
    }
    let window = WindowBuffer::new(
        SiteCompleteness::WithInvariants,
        ContigId(0),
        WINDOW as u32,
        0,
        WINDOW,
        N_SAMPLES,
        pos,
        Bits::from_bools(&is_variant),
        genotypes,
        Bits::from_bools(&missing),
        Bits::from_bools(&phased),
    )
    .expect("invariant: window fixture is valid");
    (window, rows)
}

fn kinds() -> Vec<StatKind> {
    vec![
        StatKind::Pi { group: 0, per_bp: false },
        StatKind::Pi { group: 1, per_bp: true },
        StatKind::Dxy { a: 0, b: 1, per_bp: false },
        StatKind::Dxy { a: 1, b: 2, per_bp: true },
        StatKind::Fst { a: 0, b: 1, kind: FstKind::Hudson },
        StatKind::Fst { a: 0, b: 1, kind: FstKind::WeirCockerham84 },
        StatKind::Fst { a: 1, b: 2, kind: FstKind::Hudson },
        StatKind::Pbs { a: 0, b: 1, c: 2, kind: PbsKind::Yi2010 },
        StatKind::Pbs { a: 0, b: 1, c: 2, kind: PbsKind::PbsN1 },
        StatKind::Sfs { group: 0, folded: false },
        StatKind::Sfs { group: 1, folded: true },
        StatKind::ThetaW { group: 2, per_bp: false },
        StatKind::ThetaW { group: 2, per_bp: true },
        StatKind::TajimaD { group: 0, min_sites: 0 },
    ]
}

#[test]
fn accumulators_match_oracle_on_random_windows() {
    let pops = PopulationView::from_assignment(N_SAMPLES as usize, &[0, 0, 1, 1, 2, 2], 3)
        .expect("invariant: assignment covers all samples");
    for case in 0..128_u64 {
        let n_sites = 1 + (case % 15) as usize;
        let (window, rows) = random_window(0x9E37_79B9_7F4A_7C15 ^ (case * 0x100), n_sites);
        let kinds = kinds();
        let mut accs: Vec<StatAcc> = build_accs(&kinds, &[2, 2, 2]).expect("static request");
        for row in &rows {
            let any_alt = row.iter().any(|g| g.n_alt_alleles > 0);
            for acc in &mut accs {
                acc.push(row, any_alt);
            }
        }
/// Index each accumulator by the request it was built for.
fn find<'a>(
    accs: &'a [StatAcc],
    kinds: &'a [StatKind],
    pred: impl Fn(&StatKind) -> bool,
) -> &'a StatAcc {
    accs.iter()
        .zip(kinds)
        .find(|(_, kind)| pred(kind))
        .map(|(acc, _)| acc)
        .expect("invariant: request present")
}

        for group in 0_u32..3 {
            // π, both conventions.
            if group == 0 {
                let per_pair = pi(&window, &pops, group, PiEstimator::PerPairComparable).unwrap();
                let acc = find(&accs, &kinds, |kind| matches!(kind, StatKind::Pi { group: 0, per_bp: false }));
                let StatAcc::Pi(pi_acc) = acc else { unreachable!("built for Pi") };
                assert_eq!(pi_acc.counts().diffs, per_pair.counts.diffs, "case {case} π diffs");
                assert_eq!(
                    pi_acc.counts().comparable, per_pair.counts.comparable,
                    "case {case} π comparable"
                );
                assert_eq!(pi_acc.counts().missing, per_pair.counts.missing, "case {case} π missing");
                assert_close(pi_acc.value_per_pair(), per_pair.value(), "π per-pair");
            }
            if group == 1 {
                let per_bp = pi(&window, &pops, group, PiEstimator::AccessibleBp).unwrap();
                let acc = find(&accs, &kinds, |kind| matches!(kind, StatKind::Pi { group: 1, per_bp: true }));
                let StatAcc::Pi(pi_acc) = acc else { unreachable!("built for Pi") };
                assert_close(
                    pi_acc.value_accessible_bp(WINDOW),
                    per_bp.value(),
                    "π accessible bp",
                );
            }

            // SFS / θ_W / Tajima's D (complete-case family).
            if group == 0 {
                let spectrum = sfs(&window, &pops, group, false).unwrap();
                let acc = find(&accs, &kinds, |kind| matches!(kind, StatKind::Sfs { folded: false, .. }));
                let StatAcc::Sfs(sfs_acc) = acc else { unreachable!("built for Sfs") };
                assert_eq!(sfs_acc.bins(), spectrum.counts.as_slice(), "case {case} SFS bins");
                // θ and Tajima's D no longer follow the v0.2 oracle: doc_v2
                // 04-algorithms.md §1 retires the complete-case / site-count
                // semantics in favour of pixy's mutation-count form
                // (k − 1 per site), which the gold/popgen pixy parity locks
                // instead (crates/reed/tests/e2e_pixy.rs).
                let _ = tajima_d(&window, &pops, group, 0);
            }
            if group == 1 {
                let spectrum = sfs(&window, &pops, group, true).unwrap();
                let acc = find(&accs, &kinds, |kind| matches!(kind, StatKind::Sfs { folded: true, .. }));
                let StatAcc::Sfs(sfs_acc) = acc else { unreachable!("built for Sfs") };
                assert_eq!(sfs_acc.bins(), spectrum.counts.as_slice(), "case {case} folded bins");
            }
            if group == 2 {
                // θ_W is pixy's per-allele-class form now (see the Tajima
                // note above); the oracle's site-count θ is retired.
                let _ = watterson_theta(&window, &pops, group);
            }
        }

        // dxy, all three Fst estimators, both PBS flavours.
        let d = dxy(&window, &pops, 0, 1, DxyEstimator::PerPairComparable).unwrap();
        let acc = find(&accs, &kinds, |kind| matches!(kind, StatKind::Dxy { a: 0, b: 1, per_bp: false }));
        let StatAcc::Dxy(dxy_acc) = acc else { unreachable!("built for Dxy") };
        assert_eq!(dxy_acc.counts().diffs, d.counts.diffs, "case {case} dxy diffs");
        assert_eq!(dxy_acc.counts().comparable, d.counts.comparable, "case {case} dxy comps");
        assert_close(dxy_acc.value_per_pair(), d.value(), "dxy per-pair");

        for (a, b, estimator, kind) in [
            (0_usize, 1_usize, FstEstimator::Hudson, FstKind::Hudson),
            (0, 1, FstEstimator::WeirCockerham84, FstKind::WeirCockerham84),
            (1, 2, FstEstimator::Hudson, FstKind::Hudson),
        ] {
            let oracle = fst(&window, &pops, u32::try_from(a).expect("static"), u32::try_from(b).expect("static"), estimator).unwrap();
            let acc = find(&accs, &kinds, |kind| {
                matches!(kind, StatKind::Fst { a: ka, b: kb, kind: _ } if *ka == a && *kb == b)
            });
            let StatAcc::Fst(fst_acc) = acc else { unreachable!("built for Fst") };
            let (within_a, within_b, between) = fst_acc.audit();
            assert_eq!(within_a.diffs, oracle.within_a.diffs, "case {case} Fst within_a");
            assert_eq!(
                within_a.comparable, oracle.within_a.comparable,
                "case {case} Fst within_a comps"
            );
            assert_eq!(within_b.comparable, oracle.within_b.comparable);
            assert_eq!(between.comparable, oracle.between.comparable);
            assert_close(fst_acc.value(kind), oracle.value(), estimator.id());
        }

        for kind in [PbsKind::Yi2010, PbsKind::PbsN1] {
            let estimator = match kind {
                PbsKind::Yi2010 => PbsEstimator::Yi2010,
                PbsKind::PbsN1 => PbsEstimator::PbsN1,
            };
            let oracle = pbs(&window, &pops, 0, 1, 2, estimator).unwrap();
            let acc = find(&accs, &kinds, |k| matches!(k, StatKind::Pbs { kind: pk, .. } if *pk == kind));
            let StatAcc::Pbs(pbs_acc) = acc else { unreachable!("built for Pbs") };
            let (n_sites, n_skipped) = pbs_acc.site_tally();
            assert_eq!(n_sites, oracle.n_sites, "case {case} PBS sites");
            assert_eq!(n_skipped, oracle.n_skipped, "case {case} PBS skipped");
            assert_close(pbs_acc.value(kind), oracle.value(), "PBS value");
        }
    }
}
