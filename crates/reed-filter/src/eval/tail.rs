//! Post-mask site rules (§11.2 Tail stage, fixed order: missing → MAF →
//! MAC → HWE; HWE is last because it is the most expensive).

use reed_core::count::SiteCounts;

use crate::plan::TailRule;

/// Evaluates the tail rules in order; increments `hits[i]` for every rule
/// the site passes and returns whether it survives.
pub fn eval_tail(
    values: &crate::expr::SiteValues<'_>,
    counts: &SiteCounts,
    n_samples_out: usize,
    rules: &[TailRule],
    hits: &mut [u64],
) -> bool {
    // A short counter slice would silently skip rules (caller invariant).
    assert_eq!(rules.len(), hits.len(), "invariant: tail counters align with the rules");
    // The call-rate denominator: every kept sample carries two allele slots
    // (doc_v4 02-cli-semantics.md §3.1 — the numerator is `counts.an`, so a
    // half call contributes one called allele).
    let total_alleles = 2 * n_samples_out as u64;
    for (rule, hit) in rules.iter().zip(hits.iter_mut()) {
        let ok = match rule {
            TailRule::Expr { include, expr, .. } => {
                let truth = crate::expr::eval(expr, values);
                if *include {
                    truth
                } else {
                    !truth
                }
            }
            TailRule::MinAf(bound) => non_reference_af(values, counts)
                .is_some_and(|freq| freq >= *bound),
            TailRule::MaxAf(bound) => non_reference_af(values, counts)
                .is_some_and(|freq| freq <= *bound),
            TailRule::MinCallRate(rate) => {
                n_samples_out > 0 && rate.keeps(counts.an, total_alleles)
            }
            TailRule::MaxMissingCount(cap) => {
                total_alleles.saturating_sub(counts.an) <= *cap
            }
            TailRule::MinMaf(m) => counts.maf() >= *m,
            TailRule::MaxMaf(m) => counts.maf() <= *m,
            TailRule::MinMac(mac) => {
                let minor = counts.ac.min(counts.an.saturating_sub(counts.ac));
                minor >= u64::from(*mac)
            }
            TailRule::MaxMac(mac) => {
                let minor = counts.ac.min(counts.an.saturating_sub(counts.ac));
                minor <= u64::from(*mac)
            }
            TailRule::MinHweP(p) => counts.hwe_p() >= *p,
        };
        if !ok {
            return false;
        }
        *hit += 1;
    }
    true
}

/// The non-reference allele frequency (`-q/-Q`, bcftools default `nref`
/// type): the INFO AC/AN values when both are present (bcftools semantics —
/// they may be stale), else the genotype-derived counts. `None` when the
/// site has no called alleles at all (such sites fail both bounds).
fn non_reference_af(values: &crate::expr::SiteValues<'_>, counts: &SiteCounts) -> Option<f64> {
    if let Some(freq) = crate::eval::head::non_reference_af_info(values.info) {
        return Some(freq);
    }
    (counts.an > 0).then(|| counts.ac as f64 / counts.an as f64)
}

#[cfg(test)]
mod tests {
    use super::*;
    use reed_core::genotype::{GenotypePolicies, MultiallelicPolicy, PloidyPolicy};
    use reed_core::count::SiteCounts;

    fn counts(ac: u64, an: u64, hom_ref: u64, het: u64, hom_alt: u64) -> SiteCounts {
        SiteCounts { ac, an, ns: hom_ref + het + hom_alt, n_missing: 0, hom_ref, het, hom_alt }
    }

    fn values(counts: &SiteCounts) -> crate::expr::SiteValues<'static> {
        crate::expr::SiteValues::from_counts(crate::expr::GtCounts {
            ac: counts.ac,
            an: counts.an,
            n_missing: counts.n_missing,
            n_samples: 4,
        })
    }

    #[test]
    fn minor_allele_bounds_are_symmetric() {
        // 6 REF / 4 ALT over 5 samples: MAF = 0.4, MAC = 4.
        let c = counts(4, 10, 3, 1, 1);
        let v = values(&c);
        let run = |rules: &[TailRule]| {
            let mut hits = vec![0_u64; rules.len()];
            eval_tail(&v, &c, 5, rules, &mut hits)
        };
        assert!(run(&[TailRule::MinMaf(0.4)]));
        assert!(run(&[TailRule::MaxMaf(0.4)]));
        assert!(!run(&[TailRule::MinMaf(0.45)]));
        assert!(!run(&[TailRule::MaxMaf(0.35)]));
        assert!(run(&[TailRule::MinMac(4)]));
        assert!(run(&[TailRule::MaxMac(4)]));
        assert!(!run(&[TailRule::MinMac(5)]));
        assert!(!run(&[TailRule::MaxMac(3)]));
        let _ = (PloidyPolicy::HaploidAsHalfCall, MultiallelicPolicy::CollapseAlt, GenotypePolicies::default());
    }

    /// doc_v4 02 §7.6: the call-rate decisions are exact at the boundary —
    /// `X = 0.8, N = 10, missing = 2` keeps (the float form `1 − 0.8` would
    /// drop it), `X = 0.7, N = 10, missing = 3` keeps.
    #[test]
    fn call_rate_boundaries_are_exact() {
        use crate::spec::CallRate;
        // 16 of 20 alleles called: exactly 0.8.
        let c = counts(2, 16, 8, 0, 2);
        let v = values(&c);
        let run = |rules: &[TailRule]| {
            let mut hits = vec![0_u64; rules.len()];
            eval_tail(&v, &c, 10, rules, &mut hits)
        };
        assert!(run(&[TailRule::MinCallRate(CallRate::parse("0.8").expect("valid"))]));
        assert!(!run(&[TailRule::MinCallRate(CallRate::parse("0.9").expect("valid"))]));
        assert!(run(&[TailRule::MinCallRate(CallRate::parse("0.75").expect("valid"))]));
        // 14 of 20 called against 0.7 (1 − 0.7 = 0.30000000000000004: a
        // decision written through that float would drop this site).
        let c = counts(2, 14, 6, 2, 0);
        let v = values(&c);
        assert!(eval_tail(
            &v,
            &c,
            10,
            &[TailRule::MinCallRate(CallRate::parse("0.7").expect("valid"))],
            &mut [0]
        ));
        // --max-missing-count: 2N − an ≤ C. 17 of 20 called → 3 missing.
        let c = counts(2, 17, 7, 3, 0);
        let v = values(&c);
        assert!(eval_tail(&v, &c, 10, &[TailRule::MaxMissingCount(3)], &mut [0]));
        assert!(!eval_tail(&v, &c, 10, &[TailRule::MaxMissingCount(2)], &mut [0]));
    }

    #[test]
    fn bounds_compose_as_a_band() {
        // MAF in [0.2, 0.3] keeps the site; outside drops it.
        let rules = [TailRule::MinMaf(0.2), TailRule::MaxMaf(0.3)];
        let inside = counts(4, 20, 8, 0, 2); // MAF = 0.2
        let mut hits = vec![0_u64; rules.len()];
        assert!(eval_tail(&values(&inside), &inside, 10, &rules, &mut hits));
        assert_eq!(hits, vec![1, 1]);
        let outside = counts(6, 20, 7, 0, 3); // MAF = 0.3
        let mut hits = vec![0_u64; rules.len()];
        assert!(eval_tail(&values(&outside), &outside, 10, &rules, &mut hits));
        let too_common = counts(10, 20, 5, 0, 5); // MAF = 0.5
        let mut hits = vec![0_u64; rules.len()];
        assert!(!eval_tail(&values(&too_common), &too_common, 10, &rules, &mut hits));
    }
}
