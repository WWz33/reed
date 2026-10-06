//! Fixed schedule expansion: FilterSpec -> Plan.
//!
//! Not a query optimizer: the pass dependency graph is fixed (sample
//! reduction -> mask -> recompute -> site predicates). Compilation does four
//! things: it decides which passes are active, derives ascertainment,
//! computes the provenance hash, and produces the human-readable summary
//! for the output header.

use sha2::{Digest, Sha256};

use crate::spec::FilterSpec;

/// The expanded execution plan.
#[derive(Clone, Debug)]
pub struct Plan {
    /// The spec as given (rendering source for provenance).
    pub spec: FilterSpec,
    /// Whether pass0 (the genome-wide sample missing-rate reduction) is active.
    pub reduce_sample_missing: bool,
    /// sha256 (hex) of the canonical spec rendering, written to the output
    /// header `##reed_filter=`.
    pub provenance_hash: String,
    /// Human-readable one-line summary, written to the header
    /// `##reed_filter_spec=`.
    pub spec_summary: String,
}

/// Expand. Pure function: the same spec always yields the same Plan.
pub fn compile(spec: &FilterSpec) -> Plan {
    let canonical = canonical_string(spec);
    let mut hasher = Sha256::new();
    hasher.update(canonical.as_bytes());
    let provenance_hash = format!("{:x}", hasher.finalize());
    Plan {
        spec: spec.clone(),
        reduce_sample_missing: spec.cohort.is_some(),
        provenance_hash,
        spec_summary: summarize(spec),
    }
}

/// Canonical rendering of the spec.
///
/// derive(Debug) prints fields in declaration order, which is deterministic
/// for this type (no HashMap, no random ordering); it is stable within one
/// reed-filter version. Hashes are not comparable across versions:
/// provenance exists to audit one artifact, not as a content address.
fn canonical_string(spec: &FilterSpec) -> String {
    format!("{spec:?}")
}

/// Human-readable one-line summary (gold run-log style).
fn summarize(spec: &FilterSpec) -> String {
    let s = &spec.sites;
    let fmt_opt_f32 = |v: &Option<f32>| v.map_or("-".into(), |x| format!("{x}"));
    let fmt_opt_f64 = |v: &Option<f64>| v.map_or("-".into(), |x| format!("{x}"));
    let cohort = match &spec.cohort {
        Some(crate::spec::legacy::CohortRule::MaxMissing(t)) => format!("miss<={t}"),
        None => "off".into(),
    };
    let regions = if spec.regions.is_empty() {
        "all".to_string()
    } else {
        format!("inc={} exc={}", spec.regions.include.len(), spec.regions.exclude.len())
    };
    let gt = match &spec.genotypes {
        Some(m) => format!(
            "dp=[{},{}] gq>={}",
            m.dp_min.map_or("-".into(), |v| v.to_string()),
            m.dp_max.map_or("-".into(), |v| v.to_string()),
            m.gq_min.map_or("-".into(), |v| v.to_string())
        ),
        None => "off".into(),
    };
    format!(
        "cohort={cohort} regions={regions} snps={} pass={} qual>={} info(qd>={},fs<={},mq>={},sor<={}) site_miss<={} maf>={} hwe>={} gt[{gt}]",
        s.biallelic_snps,
        s.require_pass,
        fmt_opt_f32(&s.min_qual),
        fmt_opt_f64(&s.info.min_qd),
        fmt_opt_f64(&s.info.max_fs),
        fmt_opt_f64(&s.info.min_mq),
        fmt_opt_f64(&s.info.max_sor),
        fmt_opt_f64(&s.max_missing),
        fmt_opt_f64(&s.min_maf),
        fmt_opt_f64(&s.min_hwe_p),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::legacy::{CohortRule, Region};

    #[test]
    fn provenance_is_deterministic() {
        let mut spec = FilterSpec::default();
        spec.regions.include.push(Region::parse("Chr01").unwrap());
        let p1 = compile(&spec);
        let p2 = compile(&spec);
        assert_eq!(p1.provenance_hash, p2.provenance_hash);
        assert!(p1.provenance_hash.len() == 64);
        assert!(p1.reduce_sample_missing);

        spec.cohort = None;
        let p3 = compile(&spec);
        assert!(!p3.reduce_sample_missing);
        assert_ne!(p1.provenance_hash, p3.provenance_hash);
    }

    #[test]
    fn summary_reflects_spec() {
        let spec = FilterSpec::default();
        let s = compile(&spec).spec_summary;
        assert!(s.contains("qual>=30"));
        assert!(s.contains("maf>=0.05"));
        assert!(s.contains("hwe>=-")); // off by default
        assert!(s.contains("regions=all"));

        let mut spec2 = FilterSpec { cohort: Some(CohortRule::MaxMissing(0.1)), ..FilterSpec::default() };
        spec2.regions.include.push(Region::parse("Chr01").unwrap());
        let s2 = compile(&spec2).spec_summary;
        assert!(s2.contains("miss<=0.1"));
        assert!(s2.contains("inc=1"));
    }
}
