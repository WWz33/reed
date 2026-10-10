//! Post-mask site rules (§11.2 Tail stage, fixed order: missing → MAF →
//! MAC → HWE; HWE is last because it is the most expensive).

use reed_core::count::SiteCounts;
use reed_core::math::hwe::HweCache;

use crate::plan::TailRule;

/// Evaluates the tail rules in order; increments `hits[i]` for every rule
/// the site passes and returns whether it survives.
pub fn eval_tail(
    values: &crate::expr::SiteValues<'_>,
    counts: &SiteCounts,
    n_samples_out: usize,
    rules: &[TailRule],
    hits: &mut [u64],
    hwe: &mut HweCache,
) -> bool {
    // A short counter slice would silently skip rules (caller invariant).
    assert_eq!(
        rules.len(),
        hits.len(),
        "invariant: tail counters align with the rules"
    );
    // The call-rate denominator: every kept sample carries two allele slots
    // (doc_v4 02-cli-semantics.md §3.1 — the numerator is `counts.an`, so a
    // half call contributes one called allele).
    let total_alleles = 2 * n_samples_out as u64;
    for (rule, hit) in rules.iter().zip(hits.iter_mut()) {
        let ok = match rule {
            TailRule::Expr { include, expr, .. } => {
                let truth = crate::expr::eval(expr, values);
                if *include { truth } else { !truth }
            }
            TailRule::MinAf(bound, kind) => {
                af_of(values, counts, *kind).is_some_and(|freq| freq >= *bound)
            }
            TailRule::MaxAf(bound, kind) => {
                af_of(values, counts, *kind).is_some_and(|freq| freq <= *bound)
            }
            TailRule::MaxMissing(rate) => {
                n_samples_out > 0
                    && rate.tolerates(total_alleles.saturating_sub(counts.an), total_alleles)
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
            // Cached exact test (A6-16); the degenerate sentinel stays 1.0.
            TailRule::MinHweP(p) => {
                counts
                    .hwe_test_input()
                    .map_or(1.0, |(n, mac, het)| hwe.p_exact(n, mac, het))
                    >= *p
            }
        };
        if !ok {
            return false;
        }
        *hit += 1;
    }
    true
}

/// The allele frequency a `-q/-Q` threshold applies to (bcftools `:type`;
/// docs/rounds/v5/02-cli.md §3.2; `vcfview.c` `subset_vcf` + `bcf_calc_ac`
/// with `BCF_UN_INFO|BCF_UN_FMT`, read 2026-10-10):
///
/// * the per-allele counts come from INFO/AC + INFO/AN when both exist,
///   else from the genotype tallies (`counts`);
/// * `nref` = Σ(ALT counts)/an, `alt1` = AC[1]/an;
/// * `minor`/`major` = min/max over every allele's count (REF included);
/// * `nonmajor` = (an − max count)/an.
///
/// The comparison runs in f64 exactly like bcftools; a site with an == 0
/// never passes. `None` when the site has no called alleles at all (such
/// sites fail both bounds).
fn af_of(
    values: &crate::expr::SiteValues<'_>,
    counts: &SiteCounts,
    kind: crate::spec::AfKind,
) -> Option<f64> {
    use crate::spec::AfKind;
    // INFO/AC lists the ALT counts only; prepend the REF count as
    // `an − ΣALT` (bcftools' `bcf_calc_ac` fills `ac[0]` the same way).
    let info = crate::eval::head::allele_counts_info(values.info);
    let (counts_per_allele, an) = match info {
        Some((alts, an)) => {
            let alt_sum: u64 = alts.iter().sum();
            if alt_sum > an {
                return None;
            }
            let mut counts = Vec::with_capacity(alts.len() + 1);
            counts.push(an - alt_sum);
            counts.extend(alts);
            (counts, an)
        }
        None => {
            if counts.an == 0 {
                return None;
            }
            (vec![counts.an - counts.ac, counts.ac], counts.an)
        }
    };
    if an == 0 {
        return None;
    }
    let an = an as f64;
    let alt_sum: u64 = counts_per_allele.iter().skip(1).sum();
    let max_count = *counts_per_allele.iter().max()?;
    let min_count = *counts_per_allele.iter().min()?;
    Some(match kind {
        AfKind::Nref => alt_sum as f64 / an,
        AfKind::Alt1 => *counts_per_allele.get(1)? as f64 / an,
        AfKind::Minor => min_count as f64 / an,
        AfKind::Major => max_count as f64 / an,
        AfKind::Nonmajor => (an - max_count as f64) / an,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use reed_core::count::SiteCounts;
    use reed_core::genotype::{GenotypePolicies, MultiallelicPolicy, PloidyPolicy};

    fn counts(ac: u64, an: u64, hom_ref: u64, het: u64, hom_alt: u64) -> SiteCounts {
        SiteCounts {
            ac,
            an,
            ns: hom_ref + het + hom_alt,
            n_missing: 0,
            hom_ref,
            het,
            hom_alt,
        }
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
            eval_tail(&v, &c, 5, rules, &mut hits, &mut HweCache::new())
        };
        assert!(run(&[TailRule::MinMaf(0.4)]));
        assert!(run(&[TailRule::MaxMaf(0.4)]));
        assert!(!run(&[TailRule::MinMaf(0.45)]));
        assert!(!run(&[TailRule::MaxMaf(0.35)]));
        assert!(run(&[TailRule::MinMac(4)]));
        assert!(run(&[TailRule::MaxMac(4)]));
        assert!(!run(&[TailRule::MinMac(5)]));
        assert!(!run(&[TailRule::MaxMac(3)]));
        let _ = (
            PloidyPolicy::HaploidAsHalfCall,
            MultiallelicPolicy::CollapseAlt,
            GenotypePolicies::default(),
        );
    }

    /// doc_v4 02 回退版 §7.6: the decisions are exact at the boundary —
    /// `X = 0.2, N = 10` (20 allele slots) with 2 missing keeps (equality)
    /// and 3 missing drops; `X = 0.3` with 3 missing keeps (0.3 is not a
    /// binary fraction).
    #[test]
    fn missing_rate_boundaries_are_exact() {
        use crate::spec::MissingRate;
        // 18 of 20 alleles called → 2 missing of 20: exactly 0.2.
        let c = counts(2, 18, 8, 2, 0);
        let v = values(&c);
        let run = |rules: &[TailRule]| {
            let mut hits = vec![0_u64; rules.len()];
            eval_tail(&v, &c, 10, rules, &mut hits, &mut HweCache::new())
        };
        assert!(run(&[TailRule::MaxMissing(
            MissingRate::parse("0.2").expect("valid")
        )]));
        assert!(!run(&[TailRule::MaxMissing(
            MissingRate::parse("0.05").expect("valid")
        )]));
        // 17 of 20 called → 3 missing (0.15): 0.3 keeps, 0.1 drops.
        let c = counts(2, 17, 7, 3, 0);
        let v = values(&c);
        assert!(eval_tail(
            &v,
            &c,
            10,
            &[TailRule::MaxMissing(
                MissingRate::parse("0.3").expect("valid")
            )],
            &mut [0],
            &mut HweCache::new()
        ));
        assert!(!eval_tail(
            &v,
            &c,
            10,
            &[TailRule::MaxMissing(
                MissingRate::parse("0.1").expect("valid")
            )],
            &mut [0],
            &mut HweCache::new()
        ));
    }

    #[test]
    fn bounds_compose_as_a_band() {
        // MAF in [0.2, 0.3] keeps the site; outside drops it.
        let rules = [TailRule::MinMaf(0.2), TailRule::MaxMaf(0.3)];
        let inside = counts(4, 20, 8, 0, 2); // MAF = 0.2
        let mut hits = vec![0_u64; rules.len()];
        assert!(eval_tail(
            &values(&inside),
            &inside,
            10,
            &rules,
            &mut hits,
            &mut HweCache::new()
        ));
        assert_eq!(hits, vec![1, 1]);
        let outside = counts(6, 20, 7, 0, 3); // MAF = 0.3
        let mut hits = vec![0_u64; rules.len()];
        assert!(eval_tail(
            &values(&outside),
            &outside,
            10,
            &rules,
            &mut hits,
            &mut HweCache::new()
        ));
        let too_common = counts(10, 20, 5, 0, 5); // MAF = 0.5
        let mut hits = vec![0_u64; rules.len()];
        assert!(!eval_tail(
            &values(&too_common),
            &too_common,
            10,
            &rules,
            &mut hits,
            &mut HweCache::new()
        ));
    }
}
