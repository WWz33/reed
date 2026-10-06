//! Post-mask site rules (§11.2 Tail stage, fixed order: missing → MAF →
//! MAC → HWE; HWE is last because it is the most expensive).

use reed_core::estimators::site_counts::SiteCounts;

use crate::plan::TailRule;

/// Evaluates the tail rules in order; increments `hits[i]` for every rule
/// the site passes and returns whether it survives.
pub fn eval_tail(counts: &SiteCounts, n_samples_out: usize, rules: &[TailRule], hits: &mut [u64]) -> bool {
    for (rule, hit) in rules.iter().zip(hits.iter_mut()) {
        let ok = match rule {
            TailRule::MaxMissing(m) => {
                n_samples_out > 0
                    && f64::from_u64(counts.n_missing) / f64::from_u64(n_samples_out as u64) <= *m
            }
            TailRule::MinMaf(m) => counts.maf() >= *m,
            TailRule::MinMac(mac) => {
                let minor = counts.ac.min(counts.an.saturating_sub(counts.ac));
                minor >= u64::from(*mac)
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

trait FromU64 {
    fn from_u64(v: u64) -> Self;
}
impl FromU64 for f64 {
    fn from_u64(v: u64) -> Self {
        v as f64
    }
}
