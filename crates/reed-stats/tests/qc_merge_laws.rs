//! Accumulator laws for every reed-qc accumulator (`docs/architecture-spec.md`
//! §8.4): identity, associativity, split invariance.

use proptest::prelude::*;

use reed_core::site::VariantClass;
use reed_stats::summary::SummaryAcc;

fn arb_class() -> impl Strategy<Value = VariantClass> {
    prop_oneof![
        Just(VariantClass::Snv),
        Just(VariantClass::Mnv),
        Just(VariantClass::Indel),
        Just(VariantClass::Symbolic),
        Just(VariantClass::Ref),
        Just(VariantClass::Other),
    ]
}

fn arb_acc() -> impl Strategy<Value = SummaryAcc> {
    (prop::collection::vec((arb_class(), 0_u8..3), 0..50)).prop_map(|records| {
        let mut acc = SummaryAcc::default();
        for (class, n_alt) in records {
            acc.push(class, n_alt);
        }
        acc
    })
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
    fn split_invariance_holds(records in prop::collection::vec((arb_class(), 0_u8..3), 0..200)) {
        let mut whole = SummaryAcc::default();
        for (class, n_alt) in &records {
            whole.push(*class, *n_alt);
        }
        let (left, right) = records.split_at(records.len() / 2);
        let mut halves = SummaryAcc::default();
        for (class, n_alt) in left {
            halves.push(*class, *n_alt);
        }
        for (class, n_alt) in right {
            halves.push(*class, *n_alt);
        }
        prop_assert_eq!(whole, halves);
    }
}
