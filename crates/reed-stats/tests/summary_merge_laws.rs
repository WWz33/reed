//! Accumulator laws for the live QC summary accumulator
//! (`docs/spec/architecture.md` §8.4): identity, associativity, split
//! invariance. The target is `summary::SummaryAcc` as re-exported from
//! `summary::report` (the product path; doc_v15 D15-37).

use proptest::prelude::*;

use reed_stats::summary::SummaryAcc;

fn arb_field() -> impl Strategy<Value = u64> {
    0_u64..10_000
}

fn arb_acc() -> impl Strategy<Value = SummaryAcc> {
    (
        arb_field(),
        arb_field(),
        arb_field(),
        arb_field(),
        arb_field(),
        arb_field(),
        arb_field(),
        arb_field(),
    )
        .prop_map(
            |(
                n_records,
                n_snv,
                n_mnv,
                n_indel,
                n_symbolic,
                n_ref,
                n_multiallelic,
                n_multiallelic_snv,
            )| SummaryAcc {
                n_records,
                n_snv,
                n_mnv,
                n_indel,
                n_symbolic,
                n_ref,
                n_multiallelic,
                n_multiallelic_snv,
            },
        )
}

proptest! {
    #[test]
    fn identity_holds(acc in arb_acc()) {
        let mut merged = acc;
        merged.merge(&SummaryAcc::default());
        prop_assert_eq!(merged, acc);
    }

    #[test]
    fn associativity_holds(a in arb_acc(), b in arb_acc(), c in arb_acc()) {
        let mut ab = a;
        ab.merge(&b);
        ab.merge(&c);
        let mut bc = b;
        bc.merge(&c);
        let mut abc = a;
        abc.merge(&bc);
        prop_assert_eq!(ab, abc);
    }

    #[test]
    fn split_invariance_holds(parts in prop::collection::vec(arb_acc(), 0..20)) {
        let mut whole = SummaryAcc::default();
        for p in &parts {
            whole.merge(p);
        }
        let (left, right) = parts.split_at(parts.len() / 2);
        let mut halves = SummaryAcc::default();
        for p in left {
            halves.merge(p);
        }
        for p in right {
            halves.merge(p);
        }
        prop_assert_eq!(whole, halves);
    }
}
