//! The sealed filter specification (value types).
//!
//! Contract: `FilterSpec` is the single source of semantic truth. The three
//! external executors (bcftools/vcftools/plink2) and reed each translate it
//! into their own execution, and [`FilterSpec::from_conf`] parses
//! `gold/config/filters.conf` (KEY=VALUE) directly, so the parity run
//! aligns at construction and parameter drift is impossible.
//!
//! Population stratification is out of scope by decision
//! (docs/filter-architecture.md §8): the filter recognizes only the cohort;
//! population grouping belongs to the statistics layer
//! (reed-core::PopulationView).

use std::path::Path;

use crate::FilterError;

/// The filter specification (sealed: the rule set and its semantics are
/// inseparable).
#[derive(Clone, Debug, PartialEq)]
pub struct FilterSpec {
    /// Sample-level rule (the pass0 global reduction); None = no sample filtering.
    pub cohort: Option<CohortRule>,
    /// Region rule (site predicate with targets semantics; never needs an index).
    pub regions: RegionRule,
    /// Site-level rules.
    pub sites: SiteRules,
    /// Genotype-level mask (sets `./.` for DP/GQ below threshold); None = no masking.
    pub genotypes: Option<GtMask>,
    /// INFO tags recomputed and overwritten after masking.
    pub recompute: Recompute,
}

impl Default for FilterSpec {
    /// Built-in defaults are the defaults of gold/config/filters.conf
    /// (profile-driven; see the conf comments).
    fn default() -> Self {
        FilterSpec {
            cohort: Some(CohortRule::MaxMissing(0.20)),
            regions: RegionRule::default(),
            sites: SiteRules {
                biallelic_snps: true,
                require_pass: true,
                min_qual: Some(30.0),
                info: InfoRules {
                    min_qd: Some(4.0),
                    max_fs: Some(50.0),
                    min_mq: Some(40.0),
                    max_sor: Some(3.0),
                },
                max_missing: Some(0.10),
                min_maf: Some(0.05),
                // Selfing species (soybean): HWE does not apply to this data,
                // so it is off by default (91.5% of sites measured p < 1e-6).
                min_hwe_p: None,
            },
            genotypes: Some(GtMask {
                dp_min: Some(5),
                dp_max: Some(100),
                gq_min: Some(20),
            }),
            recompute: Recompute::ALL,
        }
    }
}

/// Sample-level rule. v0.1 has MaxMissing only; a keep-list file follows in v0.2.
#[derive(Clone, Debug, PartialEq)]
pub enum CohortRule {
    /// Drop samples whose per-sample missing rate (raw counts, genome-wide
    /// reduction) exceeds the threshold.
    MaxMissing(f64),
}

/// Region rule (bcftools `--targets` streaming semantics: a site is kept
/// when it falls in any include region and in no exclude region).
///
/// The rule never needs an index; without one it degrades to a per-site
/// predicate. Index-based skipping and region parallelism are engine-level
/// execution optimizations (docs/filter-architecture.md §4/§7) and do not
/// change this semantics.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct RegionRule {
    /// Empty = keep everything.
    pub include: Vec<Region>,
    /// Excluded regions (a site in any exclude region is dropped).
    pub exclude: Vec<Region>,
}

impl RegionRule {
    /// True when include and exclude are both empty (keep everything; the
    /// predicate is always true).
    pub fn is_empty(&self) -> bool {
        self.include.is_empty() && self.exclude.is_empty()
    }

    /// Whether the site is kept.
    pub fn contains(&self, contig: &str, pos: u64) -> bool {
        let included = self.include.is_empty()
            || self
                .include
                .iter()
                .any(|r| r.contains(contig, pos));
        included && !self.exclude.iter().any(|r| r.contains(contig, pos))
    }
}

/// A genomic interval (1-based closed; a None start or end is open on that side).
#[derive(Clone, Debug, PartialEq)]
pub struct Region {
    /// Contig name (exact, case-sensitive match).
    pub contig: String,
    /// Inclusive 1-based lower bound (None = open).
    pub start: Option<u64>,
    /// Inclusive 1-based upper bound (None = open).
    pub end: Option<u64>,
}

impl Region {
    /// Parses `Chr01` / `Chr01:100-200` / `Chr01:100-` / `Chr01:-200`.
    pub fn parse(s: &str) -> Result<Self, FilterError> {
        let s = s.trim();
        let (contig, interval) = match s.split_once(':') {
            Some((c, iv)) => (c, Some(iv)),
            None => (s, None),
        };
        if contig.is_empty() {
            return Err(FilterError::Spec(format!("empty contig in region {s:?}")));
        }
        let (start, end) = match interval {
            None => (None, None),
            Some(iv) => {
                let (a, b) = iv
                    .split_once('-')
                    .ok_or_else(|| FilterError::Spec(format!("region is missing '-': {s:?}")))?;
                let p = |x: &str, what: &str| -> Result<Option<u64>, FilterError> {
                    if x.is_empty() {
                        return Ok(None);
                    }
                    x.parse::<u64>()
                        .map(Some)
                        .map_err(|_| FilterError::Spec(format!("region {what} is not a number: {s:?}")))
                };
                (p(a, "start")?, p(b, "end")?)
            }
        };
        Ok(Region {
            contig: contig.to_string(),
            start,
            end,
        })
    }

    /// 1-based closed-interval containment test.
    pub fn contains(&self, contig: &str, pos: u64) -> bool {
        self.contig == contig
            && self.start.is_none_or(|s| pos >= s)
            && self.end.is_none_or(|e| pos <= e)
    }
}

/// Site-level rules. Evaluation order matches the gold pipeline (bcftools
/// branch): type -> FILTER -> QUAL -> INFO (raw values) -> [mask] ->
/// missing/MAF/HWE (values recomputed after masking).
#[derive(Clone, Debug, PartialEq)]
pub struct SiteRules {
    /// Keep biallelic SNPs only (REF/ALT single base, ACGT, ALT != ".";
    /// matches bcftools -v snps intersected with plink2 --snps-only
    /// just-acgt).
    pub biallelic_snps: bool,
    /// The FILTER column must be exactly PASS (matches bcftools -f PASS;
    /// "." counts as failing).
    pub require_pass: bool,
    /// QUAL lower bound (a missing QUAL fails).
    pub min_qual: Option<f32>,
    /// INFO hard filter (evaluated against the raw INFO; a sub-rule whose
    /// tag is missing passes, matching the bcftools -e missing-comparison
    /// semantics).
    pub info: InfoRules,
    /// Site missing-rate cap (post-mask counts; denominator = kept sample count).
    pub max_missing: Option<f64>,
    /// MAF lower bound (recomputed from GT after masking; the input INFO/AC
    /// is stale in practice and never trusted).
    pub min_maf: Option<f64>,
    /// HWE exact-test p-value lower bound (recomputed after masking;
    /// None = off, the default for selfing species).
    pub min_hwe_p: Option<f64>,
}

/// INFO hard filter (the GATK SNP_FILTER failure condition:
/// QD < 4 || FS > 50 || MQ < 40 || SOR > 3).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct InfoRules {
    /// QD lower bound (None = sub-rule off).
    pub min_qd: Option<f64>,
    /// FS upper bound (None = sub-rule off).
    pub max_fs: Option<f64>,
    /// MQ lower bound (None = sub-rule off).
    pub min_mq: Option<f64>,
    /// SOR upper bound (None = sub-rule off).
    pub max_sor: Option<f64>,
}

impl InfoRules {
    /// Whether any INFO sub-rule is active.
    pub fn any(&self) -> bool {
        self.min_qd.is_some() || self.max_fs.is_some() || self.min_mq.is_some() || self.max_sor.is_some()
    }
}

/// Genotype-level mask: samples below threshold become `./.` (matches
/// bcftools `filter -S .`).
///
/// The definition lives in reed-io's record layer (the mask rule is expressed
/// purely in DP/GQ values and travels with the record transport); the filter
/// spec only references it. Samples with missing DP/GQ are not masked (a
/// bcftools numeric comparison is false on missing). The mask changes GT only
/// and keeps the original DP/GQ values (traceable, better than plink2
/// dropping FORMAT).
pub use reed_io::record::GtMask;

/// The set of INFO tags recomputed and overwritten after masking (reed-io's
/// record layer owns the definition; [`Recompute::ALL`] is the conf default).
pub use reed_io::record::Recompute;

impl FilterSpec {
    /// Starts a typed builder from [`FilterSpec::default`].
    pub fn builder() -> FilterSpecBuilder {
        FilterSpecBuilder {
            spec: FilterSpec::default(),
        }
    }

    /// Parses `gold/config/filters.conf` (KEY=VALUE).
    ///
    /// Both assignment forms are accepted: `KEY="value"` and
    /// `KEY="${VAR:-default}"` (shell default). Unrecognized keys (INPUT,
    /// tool paths, THREADS, other runner concerns) are ignored.
    pub fn from_conf(path: &Path) -> Result<Self, FilterError> {
        let text = std::fs::read_to_string(path)?;
        let mut kv = std::collections::HashMap::new();
        for line in text.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') || !line.contains('=') {
                continue;
            }
            let (k, v) = line.split_once('=').expect("contains '=' checked above");
            kv.insert(k.trim().to_string(), shell_value(v.trim()));
        }
        Self::from_kv(&kv)
    }

    /// Builds from key-value pairs (the parsing half of [`Self::from_conf`];
    /// missing keys keep their defaults).
    pub fn from_kv(kv: &std::collections::HashMap<String, String>) -> Result<Self, FilterError> {
        let mut spec = FilterSpec::default();
        let get = |k: &str| kv.get(k).map(|v| v.to_string());
        let parse_f64 = |k: &str| -> Result<Option<f64>, FilterError> {
            get(k)
                .map(|v| {
                    v.parse::<f64>()
                        .map_err(|_| FilterError::Spec(format!("{k} is not a number: {v:?}")))
                })
                .transpose()
        };

        if let Some(v) = get("SAMPLE_MISS_MAX") {
            let t: f64 = v
                .parse()
                .map_err(|_| FilterError::Spec(format!("SAMPLE_MISS_MAX is not a number: {v:?}")))?;
            spec.cohort = Some(CohortRule::MaxMissing(t));
        }
        if let Some(v) = get("REGION") {
            let mut inc = Vec::new();
            for part in v.split(',') {
                let part = part.trim();
                if part.is_empty() {
                    continue;
                }
                inc.push(Region::parse(part)?);
            }
            spec.regions.include = inc;
        }
        let flag = |k: &str| -> Result<Option<bool>, FilterError> {
            get(k).map(|v| match v.as_str() {
                "1" | "true" | "yes" => Ok(true),
                "0" | "false" | "no" => Ok(false),
                other => Err(FilterError::Spec(format!("{k} is not a boolean: {other:?}"))),
            }).transpose()
        };
        if let Some(b) = flag("KEEP_BIALLELIC_SNPS")? {
            spec.sites.biallelic_snps = b;
        }
        if let Some(b) = flag("REQUIRE_PASS")? {
            spec.sites.require_pass = b;
        }
        if let Some(q) = parse_f64("MIN_QUAL")? {
            spec.sites.min_qual = Some(q as f32);
        }
        if let Some(v) = parse_f64("INFO_QD_MIN")? {
            spec.sites.info.min_qd = Some(v);
        }
        if let Some(v) = parse_f64("INFO_FS_MAX")? {
            spec.sites.info.max_fs = Some(v);
        }
        if let Some(v) = parse_f64("INFO_MQ_MIN")? {
            spec.sites.info.min_mq = Some(v);
        }
        if let Some(v) = parse_f64("INFO_SOR_MAX")? {
            spec.sites.info.max_sor = Some(v);
        }
        if let Some(v) = parse_f64("SITE_MISS_MAX")? {
            spec.sites.max_missing = Some(v);
        }
        if let Some(v) = parse_f64("MAF_MIN")? {
            spec.sites.min_maf = Some(v);
        }
        // HWE_MIN=0 means "off" (the conf default); map it to None.
        if let Some(v) = parse_f64("HWE_MIN")? {
            spec.sites.min_hwe_p = if v > 0.0 { Some(v) } else { None };
        }
        let dp_min = get("GT_DP_MIN").map(|v| {
            v.parse::<u32>()
                .map_err(|_| FilterError::Spec(format!("GT_DP_MIN is not a number: {v:?}")))
        }).transpose()?;
        let dp_max = get("GT_DP_MAX").map(|v| {
            v.parse::<u32>()
                .map_err(|_| FilterError::Spec(format!("GT_DP_MAX is not a number: {v:?}")))
        }).transpose()?;
        let gq_min = get("GT_GQ_MIN").map(|v| {
            v.parse::<u32>()
                .map_err(|_| FilterError::Spec(format!("GT_GQ_MIN is not a number: {v:?}")))
        }).transpose()?;
        if dp_min.is_some() || dp_max.is_some() || gq_min.is_some() {
            spec.genotypes = Some(GtMask { dp_min, dp_max, gq_min });
        }
        Ok(spec)
    }

    /// Expands the spec into the fixed execution schedule ([`crate::plan`]).
    pub fn compile(&self) -> crate::plan::legacy::Plan {
        crate::plan::legacy::compile(self)
    }
}

/// Strips the quoting of a shell assignment and the `${VAR:-default}`
/// default form.
fn shell_value(v: &str) -> String {
    let v = v.trim();
    let v = v.strip_prefix('"').unwrap_or(v);
    let v = v.strip_suffix('"').unwrap_or(v);
    if let Some(inner) = v.strip_prefix("${")
        && let Some(inner) = inner.strip_suffix('}')
    {
        if let Some((_name, default)) = inner.split_once(":-") {
            return default.to_string();
        }
        return String::new(); // ${VAR} with no default yields an empty string (environment unset)
    }
    v.to_string()
}

/// Typed builder (construction-time validation; v0.1 sets fields directly
/// and validation lives in build).
pub struct FilterSpecBuilder {
    spec: FilterSpec,
}

impl FilterSpecBuilder {
    /// Sets the sample-level rule.
    pub fn cohort(mut self, rule: CohortRule) -> Self {
        self.spec.cohort = Some(rule);
        self
    }
    /// Disables sample-level filtering.
    pub fn no_cohort(mut self) -> Self {
        self.spec.cohort = None;
        self
    }
    /// Sets the region rule.
    pub fn regions(mut self, rule: RegionRule) -> Self {
        self.spec.regions = rule;
        self
    }
    /// Sets the site-level rules.
    pub fn sites(mut self, rules: SiteRules) -> Self {
        self.spec.sites = rules;
        self
    }
    /// Sets the genotype mask (None = no masking).
    pub fn genotypes(mut self, mask: Option<GtMask>) -> Self {
        self.spec.genotypes = mask;
        self
    }
    /// Sets the INFO tags recomputed after masking.
    pub fn recompute(mut self, recompute: Recompute) -> Self {
        self.spec.recompute = recompute;
        self
    }
    /// Validates the spec ranges and returns the spec.
    pub fn build(self) -> Result<FilterSpec, FilterError> {
        let spec = self.spec;
        if let Some(CohortRule::MaxMissing(t)) = spec.cohort.as_ref()
            && !(0.0..=1.0).contains(t)
        {
            return Err(FilterError::Spec(format!(
                "cohort max_missing out of range [0,1]: {t}"
            )));
        }
        if let Some(m) = spec.sites.max_missing && !(0.0..=1.0).contains(&m) {
            return Err(FilterError::Spec(format!("max_missing out of range [0,1]: {m}")));
        }
        if let Some(m) = spec.sites.min_maf && !(0.0..=0.5).contains(&m) {
            return Err(FilterError::Spec(format!("min_maf out of range [0,0.5]: {m}")));
        }
        // QUAL is a Phred-scaled score: negative values are meaningless.
        if let Some(q) = spec.sites.min_qual && q < 0.0 {
            return Err(FilterError::Spec(format!("min_qual out of range [0,inf): {q}")));
        }
        if let Some(p) = spec.sites.min_hwe_p {
            // An exact-test p-value threshold lives in [0,1].
            if !(0.0..=1.0).contains(&p) {
                return Err(FilterError::Spec(format!(
                    "min_hwe_p out of range [0,1]: {p}"
                )));
            }
        }
        Ok(spec)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn region_parse_and_contains() {
        let r = Region::parse("Chr01:100-200").unwrap();
        assert_eq!(r.contig, "Chr01");
        assert!(r.contains("Chr01", 100));
        assert!(r.contains("Chr01", 200));
        assert!(!r.contains("Chr01", 99));
        assert!(!r.contains("Chr01", 201));
        assert!(!r.contains("Chr02", 150));

        let whole = Region::parse("Chr02").unwrap();
        assert!(whole.contains("Chr02", 1));
        assert!(whole.contains("Chr02", u64::MAX));
        assert!(!whole.contains("Chr1", 1)); // contig matching is exact (case and numbering)

        let open = Region::parse("Chr03:500-").unwrap();
        assert!(open.contains("Chr03", 500));
        assert!(open.contains("Chr03", 9_999_999));
        assert!(!open.contains("Chr03", 499));
    }

    #[test]
    fn region_rule_include_exclude() {
        let mut rule = RegionRule::default();
        assert!(rule.is_empty());
        assert!(rule.contains("chrX", 7)); // an empty rule keeps everything

        rule.include.push(Region::parse("Chr01").unwrap());
        assert!(!rule.is_empty());
        assert!(rule.contains("Chr01", 1));
        assert!(!rule.contains("Chr02", 1));

        rule.exclude.push(Region::parse("Chr01:1000-2000").unwrap());
        assert!(!rule.contains("Chr01", 1500));
        assert!(rule.contains("Chr01", 999));
    }

    #[test]
    fn shell_value_forms() {
        assert_eq!(shell_value("\"${MIN_QUAL:-30}\""), "30");
        assert_eq!(shell_value("\"${REGION:-}\""), "");
        assert_eq!(shell_value("\"0.10\""), "0.10");
        assert_eq!(shell_value("30"), "30");
    }

    #[test]
    fn gt_mask_semantics() {
        let m = GtMask { dp_min: Some(5), dp_max: Some(100), gq_min: Some(20) };
        // missing DP/GQ never masks (a bcftools missing comparison is false)
        assert!(!m.fails(None, None));
        assert!(!m.fails(None, Some(99)));
        // within thresholds
        assert!(!m.fails(Some(50), Some(99)));
        // each failure branch
        assert!(m.fails(Some(4), Some(99)));
        assert!(m.fails(Some(101), Some(99)));
        assert!(m.fails(Some(50), Some(19)));
    }

    #[test]
    fn builder_validation() {
        assert!(FilterSpec::builder().no_cohort().build().is_ok());
        let bad = FilterSpec::builder()
            .cohort(CohortRule::MaxMissing(1.5))
            .build();
        assert!(bad.is_err());
    }

    #[test]
    fn builder_rejects_negative_min_qual() {
        let mut spec = FilterSpec::default();
        spec.sites.min_qual = Some(-1.0);
        let err = FilterSpecBuilder { spec }.build().unwrap_err();
        assert!(matches!(err, FilterError::Spec(_)), "{err}");
        assert!(err.to_string().contains("min_qual"), "{err}");
        // A non-negative QUAL threshold is accepted.
        assert!(FilterSpec::builder().build().is_ok());
    }

    #[test]
    fn builder_rejects_min_hwe_p_outside_unit_interval() {
        let mut spec = FilterSpec::default();
        spec.sites.min_hwe_p = Some(2.0);
        let err = FilterSpecBuilder { spec }.build().unwrap_err();
        assert!(matches!(err, FilterError::Spec(_)), "{err}");
        assert!(err.to_string().contains("min_hwe_p"), "{err}");
        // The valid range, including the boundary, is accepted.
        assert!(FilterSpec::builder().build().is_ok());
        let mut spec = FilterSpec::default();
        spec.sites.min_hwe_p = Some(1.0);
        assert!(FilterSpecBuilder { spec }.build().is_ok());
    }
}
