//! Gold-conf compatibility (§11.1): `gold/config/filters.conf` KEY=VALUE
//! read-in, existing only to serve the gold parity run. Equivalence with
//! the TOML path is locked by test (`spec_validation::compat_matches_toml`).

use std::collections::HashMap;
use std::path::Path;

use super::{
    FilterSpecV2, GenotypeMaskRule, InfoPredicate, OnMissing, RecomputeSpec,
};
use crate::FilterError;

/// Reads `filters.conf` and builds the v2 spec. Unrecognized keys (INPUT,
/// tool paths, THREADS, runner concerns) are ignored.
///
/// # Errors
///
/// [`FilterError::Spec`] on malformed values.
pub fn from_gold_conf(path: &Path) -> Result<FilterSpecV2, FilterError> {
    let text = std::fs::read_to_string(path)?;
    let mut kv = HashMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') || !line.contains('=') {
            continue;
        }
        let (k, v) = line.split_once('=').expect("contains '=' checked above");
        kv.insert(k.trim().to_string(), shell_value(v.trim()));
    }
    from_gold_kv(&kv)
}

/// The value-building half of [`from_gold_conf`].
///
/// # Errors
///
/// [`FilterError::Spec`] on malformed values.
pub fn from_gold_kv(kv: &HashMap<String, String>) -> Result<FilterSpecV2, FilterError> {
    let get = |k: &str| kv.get(k).map(String::as_str);
    let num = |k: &str| -> Result<Option<f64>, FilterError> {
        get(k)
            .map(|v| {
                v.parse::<f64>()
                    .map_err(|_| FilterError::Spec(format!("{k} is not a number: {v:?}")))
            })
            .transpose()
    };
    // The raw (unparsed) text of a key, for keep_if strings that must keep
    // the conf's literal form ("4.0", not the reformatted "4").
    let raw_num = |k: &str| -> Result<Option<String>, FilterError> {
        get(k)
            .map(|v| {
                v.parse::<f64>()
                    .map(|_| v.to_string())
                    .map_err(|_| FilterError::Spec(format!("{k} is not a number: {v:?}")))
            })
            .transpose()
    };

    let mut spec = FilterSpecV2 {
        recompute: RecomputeSpec::all(),
        ..FilterSpecV2::default()
    };

    if let Some(v) = get("SAMPLE_MISS_MAX") {
        let t: f64 = v
            .parse()
            .map_err(|_| FilterError::Spec(format!("SAMPLE_MISS_MAX is not a number: {v:?}")))?;
        spec.cohort = Some(super::CohortRule::MaxMissing(t));
    }
    let truthy = |v: Option<&str>| matches!(v, Some("1" | "true" | "yes"));
    if truthy(get("KEEP_BIALLELIC_SNPS")) {
        spec.site.variant_class.push("snv".into());
        spec.site.max_alt_alleles = 1;
    }
    if truthy(get("REQUIRE_PASS")) {
        spec.site.require_filter.push("PASS".into());
    }
    if let Some(q) = num("MIN_QUAL")? {
        spec.site.min_qual = Some(q as f32);
    }
    for (key, tag, keep_if) in [
        ("INFO_QD_MIN", "QD", ">="),
        ("INFO_FS_MAX", "FS", "<="),
        ("INFO_MQ_MIN", "MQ", ">="),
        ("INFO_SOR_MAX", "SOR", "<="),
    ] {
        if let Some(v) = raw_num(key)? {
            spec.site.info.push(InfoPredicate {
                tag: tag.into(),
                keep_if: format!("{keep_if} {v}"),
                on_missing: OnMissing::Pass,
                quantifier: None,
            });
        }
    }
    if let Some(m) = num("SITE_MISS_MAX")? {
        spec.site_after_mask.max_missing = Some(m);
    }
    if let Some(m) = num("MAF_MIN")? {
        spec.site_after_mask.min_maf = Some(m);
    }
    if let Some(h) = num("HWE_MIN")? {
        spec.site_after_mask.min_hwe_p = (h > 0.0).then_some(h);
    }
    let int_of = |k: &str| -> Result<Option<u32>, FilterError> {
        get(k)
            .map(|v| {
                v.parse::<u32>()
                    .map_err(|_| FilterError::Spec(format!("{k} is not a number: {v:?}")))
            })
            .transpose()
    };
    for (key, field, op) in [
        ("GT_DP_MIN", "DP", "<"),
        ("GT_DP_MAX", "DP", ">"),
        ("GT_GQ_MIN", "GQ", "<"),
    ] {
        if let Some(v) = int_of(key)? {
            spec.genotype.mask_if.push(GenotypeMaskRule {
                field: field.into(),
                when: format!("{op} {v}"),
                on_missing: OnMissing::Keep,
            });
        }
    }
    Ok(spec)
}

/// Strips the quoting of a shell assignment and the `${VAR:-default}`
/// default form (shared with the legacy reader).
#[must_use]
pub fn shell_value(v: &str) -> String {
    let v = v.trim();
    let v = v.strip_prefix('"').unwrap_or(v);
    let v = v.strip_suffix('"').unwrap_or(v);
    if let Some(inner) = v.strip_prefix("${")
        && let Some(inner) = inner.strip_suffix('}')
    {
        if let Some((_name, default)) = inner.split_once(":-") {
            return default.to_string();
        }
        return String::new();
    }
    v.to_string()
}
