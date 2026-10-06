//! The plan compiler (§11.3): `compile(spec, header, budget) -> FilterPlan`
//! is the only input of execution — validation (tags declared, numeric
//! types, Number ∈ {1, A}), rule expansion, strategy choice (§11.4),
//! ascertainment derivation (§11.7) and provenance (§11.8) all happen here.

pub mod ascertainment;
pub mod legacy;
pub mod provenance;
pub mod strategy;

use std::collections::HashMap;

use noodles::vcf;
use noodles::vcf::header::record::value::map::format;
use noodles::vcf::header::record::value::map::info::{Number, Type};
use reed_core::manifest::v2::Ascertainment;

pub use ascertainment::derive_ascertainment;
pub use legacy::Plan as LegacyPlan;
pub use strategy::ExecStrategy;

use crate::spec::{
    FilterSpecV2, OnMissing, Quantifier, RecomputeSpec, RegionRule,
};
use crate::FilterError;

/// A compiled head rule (fixed columns + INFO only).
#[derive(Clone, Debug)]
pub enum HeadRule {

    /// Streaming targets: keep when inside any include and no exclude.
    Regions {
        /// Parsed include regions.
        include: Vec<(String, Option<u64>, Option<u64>)>,
        /// Parsed exclude regions.
        exclude: Vec<(String, Option<u64>, Option<u64>)>,
    },
    /// The variant class must be one of these (`snv | mnv | indel |
    /// symbolic | ref`), evaluated with reed-core's `classify_variant`.
    VariantClass(Vec<String>),
    /// The record's ALT count must not exceed this.
    MaxAltAlleles(u32),
    /// The FILTER column must be one of these (`.` never passes).
    RequireFilter(Vec<Vec<u8>>),
    /// QUAL lower bound (missing QUAL fails, bcftools semantics).
    MinQual(f32),
    /// One numeric INFO predicate.
    Info(InfoRule),
}

impl HeadRule {
    /// The stable funnel id (§11.5).
    #[must_use]
    pub fn id(&self) -> String {
        match self {
            Self::Regions { .. } => "site.regions".to_string(),
            Self::VariantClass(_) => "site.variant_class".to_string(),
            Self::MaxAltAlleles(_) => "site.max_alt_alleles".to_string(),
            Self::RequireFilter(_) => "site.filter".to_string(),
            Self::MinQual(_) => "site.qual".to_string(),
            Self::Info(rule) => {
                let tag = String::from_utf8_lossy(&rule.tag);
                format!("site.info.{tag}")
            }
        }
    }
}

/// A compiled INFO predicate.
#[derive(Clone, Debug)]
pub struct InfoRule {
    /// Tag name (bytes, for fragment matching).
    pub tag: Vec<u8>,
    /// The parsed comparison.
    pub cmp: crate::spec::predicate::Cmp,
    /// Missing-tag behaviour.
    pub on_missing: OnMissing,
    /// `Number=A` quantifier (required when multiallelics can survive).
    pub quantifier: Quantifier,
    /// Whether the header declares the tag Number=A.
    pub number_a: bool,
}

/// A compiled genotype mask rule.
#[derive(Clone, Debug)]
pub struct MaskRule {
    /// FORMAT field name (the scan key).
    pub field: Vec<u8>,
    /// The parsed comparison.
    pub cmp: crate::spec::predicate::Cmp,
    /// Missing-value behaviour (`Keep` = no mask, `Mask` = mask).
    pub on_missing: OnMissing,
    /// Position of this rule's field in `int_fields` (the scan order).
    pub int_column: usize,
}

/// A compiled post-mask rule.
#[derive(Clone, Copy, Debug)]
pub enum TailRule {
    /// Site missing-rate cap (denominator: kept samples).
    MaxMissing(f64),
    /// MAF lower bound.
    MinMaf(f64),
    /// Minor-allele-count lower bound.
    MinMac(u32),
    /// HWE exact-test p lower bound.
    MinHweP(f64),
}

impl TailRule {
    /// The stable funnel id (§11.5).
    #[must_use]
    pub fn id(self) -> &'static str {
        match self {
            Self::MaxMissing(_) => "after_mask.missing",
            Self::MinMaf(_) => "after_mask.maf",
            Self::MinMac(_) => "after_mask.mac",
            Self::MinHweP(_) => "after_mask.hwe",
        }
    }
}

/// The recomputed-tag flags (compiled from [`RecomputeSpec`]).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RecomputeFlags {
    /// AC.
    pub ac: bool,
    /// AN.
    pub an: bool,
    /// NS.
    pub ns: bool,
    /// MAF.
    pub maf: bool,
    /// HWE.
    pub hwe: bool,
}

/// The compiled plan: the only input of execution (§11.3).
#[derive(Clone, Debug)]
pub struct FilterPlan {
    /// Head rules in evaluation order (§11.2: the order is semantics).
    pub head: Vec<HeadRule>,
    /// Genotype mask rules in declaration order (first hit attributes).
    pub mask: Vec<MaskRule>,
    /// Post-mask rules in fixed order (missing → MAF → MAC → HWE).
    pub tail: Vec<TailRule>,
    /// The chosen execution strategy.
    pub strategy: ExecStrategy,
    /// FORMAT fields the mask needs scanned, in scan order.
    pub int_fields: Vec<Vec<u8>>,
    /// Stale INFO tags (bytes for the editor).
    pub stale_tags: Vec<Vec<u8>>,
    /// Recompute flags.
    pub recompute: RecomputeFlags,
    /// The derived ascertainment (cohort filled at runtime, §11.7).
    pub ascertainment: Ascertainment,
    /// Provenance: canonical spec JSON + sha256 + command line.
    pub provenance: provenance::Provenance,
}

/// Compiles the v2 spec against the input header.
///
/// # Errors
///
/// [`FilterError::Spec`] on validation failures: unknown tags, non-numeric
/// types, `Number` outside {1, A}, a missing quantifier where one is
/// required, thresholds out of range, or the M7 `midp` HWE mode.
pub fn compile(
    spec: &FilterSpecV2,
    header: &vcf::Header,
    command_line: &str,
) -> Result<FilterPlan, FilterError> {
    let info_defs: HashMap<String, (Type, Number)> = header
        .infos()
        .iter()
        .map(|(k, m)| (k.to_string(), (m.ty(), m.number())))
        .collect();
    let format_defs: HashMap<String, (format::Type, format::Number)> = header
        .formats()
        .iter()
        .map(|(k, m)| (k.to_string(), (m.ty(), m.number())))
        .collect();

    // ---- head rules (§11.2 order: targets → class/alt → FILTER → QUAL →
    // INFO in declaration order) ----
    let mut head = Vec::new();
    let regions = compile_regions(&spec.targets)?;
    if let Some(r) = regions {
        head.push(HeadRule::Regions { include: r.0, exclude: r.1 });
    }
    if !spec.site.variant_class.is_empty() {
        head.push(HeadRule::VariantClass(spec.site.variant_class.clone()));
    }
    if spec.site.max_alt_alleles > 0 {
        head.push(HeadRule::MaxAltAlleles(spec.site.max_alt_alleles));
    }
    if !spec.site.require_filter.is_empty() {
        head.push(HeadRule::RequireFilter(
            spec.site.require_filter.iter().map(String::as_bytes).map(<[u8]>::to_vec).collect(),
        ));
    }
    if let Some(q) = spec.site.min_qual {
        if q < 0.0 {
            return Err(FilterError::Spec(format!("min_qual out of range [0,inf): {q}")));
        }
        head.push(HeadRule::MinQual(q));
    }
    let multiallelics_possible = spec.policy.multiallelic
        == crate::spec::MultiallelicPolicy::Allow
        || spec.site.max_alt_alleles > 1;
    for rule in &spec.site.info {
        let (ty, number) = info_defs
            .get(&rule.tag)
            .copied()
            .ok_or_else(|| FilterError::Spec(format!("INFO tag {} is not declared in the header", rule.tag)))?;
        validate_numeric(&rule.tag, ty)?;
        let number_a = matches!(number, Number::A);
        if number != Number::Count(1) && !number_a {
            return Err(FilterError::Spec(format!(
                "INFO tag {} has an unsupported Number (want 1 or A)",
                rule.tag
            )));
        }
        let quantifier = if number_a {
            rule.quantifier.ok_or_else(|| {
                FilterError::Spec(format!(
                    "INFO tag {} is Number=A: declare quantifier = \"any\" | \"all\"",
                    rule.tag
                ))
            })?
        } else {
            // A Number=1 tag with a declared quantifier would silently
            // change nothing; keep the declaration honest.
            // A declared quantifier on a Number=1 tag is harmless (the
            // value set has one element); accepted silently.
            Quantifier::Any
        };
        let cmp = crate::spec::predicate::Cmp::parse(&rule.keep_if)?;
        // Missing-tag behaviour: bcftools semantics default to pass.
        let on_missing = rule.on_missing;
        if on_missing == OnMissing::Keep || on_missing == OnMissing::Mask {
            return Err(FilterError::Spec(format!(
                "INFO rule {} uses a genotype on_missing value (pass|fail only)",
                rule.tag
            )));
        }
        head.push(HeadRule::Info(InfoRule {
            tag: rule.tag.as_bytes().to_vec(),
            cmp,
            on_missing,
            quantifier,
            number_a,
        }));
    }

    // ---- mask rules ----
    let mut mask = Vec::new();
    let mut int_fields: Vec<Vec<u8>> = Vec::new();
    for rule in &spec.genotype.mask_if {
        let (ty, _number) = format_defs
            .get(&rule.field)
            .copied()
            .ok_or_else(|| FilterError::Spec(format!("FORMAT field {} is not declared in the header", rule.field)))?;
        validate_numeric(&rule.field, ty)?;
        let cmp = crate::spec::predicate::Cmp::parse(&rule.when)?;
        // The mask default is keep (§11.1): an omitted on_missing
        // deserializes to the enum default (Pass) and normalizes to Keep.
        let on_missing = match rule.on_missing {
            OnMissing::Keep | OnMissing::Pass => OnMissing::Keep,
            OnMissing::Mask | OnMissing::Fail => OnMissing::Mask,
        };
        let int_column = match int_fields.iter().position(|f| *f == rule.field.as_bytes()) {
            Some(i) => i,
            None => {
                int_fields.push(rule.field.as_bytes().to_vec());
                int_fields.len() - 1
            }
        };
        mask.push(MaskRule { field: rule.field.as_bytes().to_vec(), cmp, on_missing, int_column });
    }

    // ---- tail rules (fixed order §11.2) ----
    let mut tail = Vec::new();
    if let Some(m) = spec.site_after_mask.max_missing {
        if !(0.0..=1.0).contains(&m) {
            return Err(FilterError::Spec(format!("max_missing out of range [0,1]: {m}")));
        }
        tail.push(TailRule::MaxMissing(m));
    }
    if let Some(m) = spec.site_after_mask.min_maf {
        if !(0.0..=0.5).contains(&m) {
            return Err(FilterError::Spec(format!("min_maf out of range [0,0.5]: {m}")));
        }
        tail.push(TailRule::MinMaf(m));
    }
    if spec.site_after_mask.min_mac > 0 {
        tail.push(TailRule::MinMac(spec.site_after_mask.min_mac));
    }
    if let Some(p) = spec.site_after_mask.min_hwe_p {
        if spec.site_after_mask.hwe_mode == crate::spec::HweMode::MidP {
            return Err(FilterError::Spec(
                "hwe_mode = midp arrives with M7 (exact only in this version)".into(),
            ));
        }
        if !(0.0..=1.0).contains(&p) {
            return Err(FilterError::Spec(format!("min_hwe_p out of range [0,1]: {p}")));
        }
        tail.push(TailRule::MinHweP(p));
    }

    // ---- strategy (§11.4; SketchThenEmit selection is M5b) ----
    let strategy = match &spec.cohort {
        Some(crate::spec::CohortRule::MaxMissing(_)) => ExecStrategy::TwoPassReparse,
        Some(crate::spec::CohortRule::KeepList(_)) | Some(crate::spec::CohortRule::RemoveList(_))
        | None => ExecStrategy::SinglePass,
    };

    // ---- recompute + stale ----
    let recompute = RecomputeFlags {
        ac: spec.recompute.wants("AC"),
        an: spec.recompute.wants("AN"),
        ns: spec.recompute.wants("NS"),
        maf: spec.recompute.wants("MAF"),
        hwe: spec.recompute.wants("HWE"),
    };
    let stale_tags = spec.recompute.stale_bytes();

    // ---- ascertainment (§11.7) + provenance (§11.8) ----
    let provenance = provenance::Provenance::new(spec, command_line)?;
    let plan_hash = provenance.sha256.clone();
    let has_site_rules = !head.is_empty() || !tail.is_empty();
    let ascertainment = derive_ascertainment(
        has_site_rules,
        &spec.site.variant_class,
        spec.site_after_mask.min_maf,
        spec.site_after_mask.min_mac,
        &mask,
        &plan_hash,
    );

    Ok(FilterPlan {
        head,
        mask,
        tail,
        strategy,
        int_fields,
        stale_tags,
        recompute,
        ascertainment,
        provenance,
    })
}

fn validate_numeric(tag: &str, ty: impl std::fmt::Debug) -> Result<(), FilterError> {
    let name = format!("{ty:?}");
    if name == "Integer" || name == "Float" {
        Ok(())
    } else {
        Err(FilterError::Spec(format!(
            "tag {tag} is {name}-typed: numeric rules need Integer or Float"
        )))
    }
}

type RegionTriple = (String, Option<u64>, Option<u64>);

/// Parsed include/exclude region lists.
type RegionLists = (Vec<RegionTriple>, Vec<RegionTriple>);

fn compile_regions(rule: &RegionRule) -> Result<Option<RegionLists>, FilterError> {
    if rule.is_empty() {
        return Ok(None);
    }
    let parse = |list: &[String]| -> Result<Vec<RegionTriple>, FilterError> {
        list.iter()
            .map(|s| {
                let r = crate::spec::region::Region::parse(s)?;
                Ok((r.contig, r.start, r.end))
            })
            .collect()
    };
    Ok(Some((parse(&rule.include)?, parse(&rule.exclude)?)))
}
