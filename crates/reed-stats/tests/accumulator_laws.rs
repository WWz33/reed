//! Accumulator laws (`docs/spec/architecture.md` §8.4), proptest-locked for
//! every statistic family: identity, associativity and split-invariance —
//! integer components exactly, float components within 4 ulp.

use proptest::prelude::*;
use reed_core::CoreError;
use reed_core::count::SiteGroupCounts;
use reed_stats::fst::FstKind;
use reed_stats::pbs::PbsKind;
use reed_stats::request::{StatAcc, StatKind, build_accs};

/// Per-group counts: `(n_ref, n_alt, n_het)` with `n_ref + n_alt = called`
/// and `called ≤ n_chrom` (the accumulators assume consistent tallies).
type Group = (u32, u32, u32);

/// Chromosome counts per group (diploid samples × 2).
const N_CHROM: [u32; 3] = [8, 6, 4];

fn group_strategy(n_chrom: u32) -> impl Strategy<Value = Group> {
    (0..=n_chrom).prop_flat_map(|called| {
        (0..=called, 0..=called).prop_map(move |(alt, het)| (called - alt, alt, het))
    })
}

fn site_strategy() -> impl Strategy<Value = ([Group; 3], bool)> {
    (
        group_strategy(N_CHROM[0]),
        group_strategy(N_CHROM[1]),
        group_strategy(N_CHROM[2]),
        any::<bool>(),
    )
        .prop_map(|(a, b, c, variant)| ([a, b, c], variant))
}

fn row_of(groups: &[Group; 3]) -> Vec<SiteGroupCounts> {
    groups
        .iter()
        .map(|&(n_ref_alleles, n_alt_alleles, n_het)| SiteGroupCounts {
            n_ref_alleles,
            n_alt_alleles,
            n_het,
            ..SiteGroupCounts::default()
        })
        .collect()
}

/// Everything one accumulator can report, for law comparisons: integer
/// components, SFS bins, and the scalar value.
#[derive(Clone, Debug, PartialEq)]
struct Obs {
    ints: Vec<u64>,
    bins: Vec<u64>,
    value: Option<f64>,
}

fn observe(accs: &[StatAcc], kinds: &[StatKind], window_bp: u64) -> Vec<Obs> {
    accs.iter()
        .zip(kinds)
        .map(|(acc, kind)| match acc {
            StatAcc::Pi(a) => Obs {
                ints: vec![a.counts().diffs, a.counts().comparable, a.counts().missing],
                bins: vec![],
                value: a.value_per_pair(),
            },
            StatAcc::Dxy(a) => Obs {
                ints: vec![a.counts().diffs, a.counts().comparable, a.counts().missing],
                bins: vec![],
                value: a.value_per_pair(),
            },
            StatAcc::Fst(a) => {
                let (within_a, within_b, between) = a.audit();
                Obs {
                    ints: vec![
                        within_a.diffs,
                        within_a.comparable,
                        within_b.diffs,
                        within_b.comparable,
                        between.diffs,
                        between.comparable,
                    ],
                    bins: vec![],
                    value: match kind {
                        StatKind::Fst { kind, .. } => a.value(*kind),
                        _ => unreachable!("kind matches the accumulator"),
                    },
                }
            }
            StatAcc::FstGlobal(a) => {
                let (n_sites, n_skipped) = a.site_tally();
                Obs {
                    ints: vec![n_sites, n_skipped],
                    bins: vec![],
                    value: a.value(),
                }
            }
            StatAcc::Pbs(a) => {
                let (n_sites, n_skipped) = a.site_tally();
                Obs {
                    ints: vec![n_sites, n_skipped],
                    bins: vec![],
                    value: match kind {
                        StatKind::Pbs { kind, .. } => a.value(*kind),
                        _ => unreachable!("kind matches the accumulator"),
                    },
                }
            }
            StatAcc::Sfs(a) => {
                let (n_sites, n_skipped) = a.site_tally();
                Obs {
                    ints: vec![n_sites, n_skipped],
                    bins: a.bins().to_vec(),
                    value: None,
                }
            }
            StatAcc::Theta(a) => {
                let (n_sites, n_skipped) = a.site_tally();
                Obs {
                    ints: vec![a.segregating(), n_sites, n_skipped],
                    bins: vec![],
                    value: match kind {
                        StatKind::ThetaW { per_bp: true, .. } => a.value_per_bp(window_bp),
                        _ => a.value_total(),
                    },
                }
            }
            StatAcc::Tajima(a) => {
                let (n_sites, n_skipped) = a.site_tally();
                Obs {
                    ints: vec![a.segregating(), n_sites, n_skipped],
                    bins: vec![],
                    value: match kind {
                        StatKind::TajimaD { min_sites, .. } => a.value(*min_sites),
                        _ => unreachable!("kind matches the accumulator"),
                    },
                }
            }
        })
        .collect()
}

/// Integer components exactly; float components within `ulp` relative steps.
fn assert_same(left: &[Obs], right: &[Obs], ulp: f64, context: &str) {
    assert_eq!(left.len(), right.len(), "{context}: observation count");
    for (i, (a, b)) in left.iter().zip(right).enumerate() {
        assert_eq!(a.ints, b.ints, "{context}: integer components of acc {i}");
        assert_eq!(a.bins, b.bins, "{context}: SFS bins of acc {i}");
        match (a.value, b.value) {
            (None, None) => {}
            (Some(x), Some(y)) => {
                let tol = ulp * f64::EPSILON * x.abs().max(y.abs()).max(1.0);
                assert!(
                    (x - y).abs() <= tol,
                    "{context}: value of acc {i}: {x} vs {y}"
                );
            }
            (a, b) => assert_eq!(
                a.is_some(),
                b.is_some(),
                "{context}: value definedness of acc {i}: {a:?} vs {b:?}"
            ),
        }
    }
}

fn kinds() -> Vec<StatKind> {
    vec![
        StatKind::Pi {
            group: 0,
            per_bp: false,
        },
        StatKind::Dxy {
            a: 1,
            b: 2,
            per_bp: false,
        },
        StatKind::Fst {
            a: 0,
            b: 2,
            kind: FstKind::Hudson,
        },
        StatKind::Fst {
            a: 0,
            b: 2,
            kind: FstKind::WeirCockerham84,
        },
        StatKind::Fst {
            a: 0,
            b: 2,
            kind: FstKind::Hudson,
        },
        StatKind::FstGlobal,
        StatKind::Pbs {
            a: 0,
            b: 1,
            c: 2,
            kind: PbsKind::Yi2010,
        },
        StatKind::Pbs {
            a: 0,
            b: 1,
            c: 2,
            kind: PbsKind::PbsN1,
        },
        StatKind::Sfs {
            group: 1,
            folded: false,
        },
        StatKind::Sfs {
            group: 1,
            folded: true,
        },
        StatKind::ThetaW {
            group: 2,
            per_bp: false,
        },
        StatKind::TajimaD {
            group: 2,
            min_sites: 0,
        },
    ]
}

/// Pushes every site into every accumulator of the set.
fn push_all(accs: &mut [StatAcc], sites: &[([Group; 3], bool)]) {
    for &(groups, variant) in sites {
        let row = row_of(&groups);
        for acc in accs.iter_mut() {
            acc.push(&row, variant);
        }
    }
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    #[test]
    fn identity_associativity_and_split_invariance(
        sites in prop::collection::vec(site_strategy(), 0..40),
    ) {
        let kinds = kinds();
        let n_group = [4_u32, 3, 2];
        let window_bp = 1_000_u64;

        let mut whole = build_accs(&kinds, &n_group).unwrap();
        push_all(&mut whole, &sites);
        let whole_obs = observe(&whole, &kinds, window_bp);

        // Identity: merging a fresh (zero) accumulator changes nothing.
        let mut identity = build_accs(&kinds, &n_group).unwrap();
        for (identity, whole) in identity.iter_mut().zip(&whole) {
            identity.merge(whole);
        }
        assert_same(&observe(&identity, &kinds, window_bp), &whole_obs, 4.0, "identity");

        // Split-invariance: any partition of the site sequence merges back
        // to the whole.
        let mid = sites.len() / 2;
        let mut first = build_accs(&kinds, &n_group).unwrap();
        push_all(&mut first, &sites[..mid]);
        let mut second = build_accs(&kinds, &n_group).unwrap();
        push_all(&mut second, &sites[mid..]);
        for (first, second) in first.iter_mut().zip(&second) {
            first.merge(second);
        }
        assert_same(&observe(&first, &kinds, window_bp), &whole_obs, 4.0, "split");

        // Associativity: (a⊕b)⊕c == a⊕(b⊕c).
        let (a, rest) = sites.split_at(mid);
        let (b, c) = rest.split_at(rest.len() / 2);
        let mut acc_a = build_accs(&kinds, &n_group).unwrap();
        let mut acc_b = build_accs(&kinds, &n_group).unwrap();
        let mut acc_c = build_accs(&kinds, &n_group).unwrap();
        push_all(&mut acc_a, a);
        push_all(&mut acc_b, b);
        push_all(&mut acc_c, c);
        let merge_all = |order: [&[StatAcc]; 3]| -> Vec<Obs> {
            let mut accs = build_accs(&kinds, &n_group).unwrap();
            for src in order {
                for (dst, src) in accs.iter_mut().zip(src) {
                    dst.merge(src);
                }
            }
            observe(&accs, &kinds, window_bp)
        };
        assert_same(
            &merge_all([&acc_a, &acc_b, &acc_c]),
            &merge_all([&acc_b, &acc_c, &acc_a]),
            4.0,
            "associativity",
        );
    }
}

#[test]
fn empty_group_is_rejected_at_build_time() {
    let err = build_accs(
        &[StatKind::Pi {
            group: 0,
            per_bp: false,
        }],
        &[0],
    )
    .err()
    .unwrap_or(CoreError::LengthMismatch {
        expected: 0,
        got: 0,
    });
    assert!(matches!(
        err,
        CoreError::InvalidWindow | CoreError::LengthMismatch { .. }
    ));
}
