//! The plan compiler (§11.3): `compile(spec, header, budget) -> FilterPlan`
//! is the only input of execution — validation (tags declared, numeric
//! types, Number ∈ {1, A}), rule expansion, strategy choice (§11.4),
//! ascertainment derivation (§11.7) and provenance (§11.8) all happen here.

pub mod provenance;
pub mod strategy;


use reed_core::tag_table::{TagNumber, TagTable, TagType};

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
    /// The variant class must not be any of these (bcftools `-V`).
    ExcludeVariantClass(Vec<String>),
    /// The record's ALT count must be at least this.
    MinAltAlleles(u32),
    /// The record's ALT count must not exceed this.
    MaxAltAlleles(u32),
    /// The FILTER column must be one of these (`.` never passes).
    RequireFilter(Vec<Vec<u8>>),
    /// QUAL lower bound (missing QUAL fails, bcftools semantics).
    MinQual(f32),
    /// One numeric INFO predicate.
    Info(InfoRule),
    /// A bcftools expression (02-cli.md §4): `include` must hold, `exclude`
    /// must not.
    Expr {
        /// Funnel id (`site.include#N` / `site.exclude#N`).
        id: String,
        /// Whether the expression is an `-i` (true = keep) or `-e`.
        include: bool,
        /// The compiled predicate.
        expr: crate::expr::Expr,
    },
}

impl HeadRule {
    /// The stable funnel id (§11.5).
    #[must_use]
    pub fn id(&self) -> String {
        match self {
            Self::Regions { .. } => "site.regions".to_string(),
            Self::VariantClass(_) => "site.variant_class".to_string(),
            Self::ExcludeVariantClass(_) => "site.exclude_variant_class".to_string(),
            Self::MinAltAlleles(_) => "site.min_alt_alleles".to_string(),
            Self::MaxAltAlleles(_) => "site.max_alt_alleles".to_string(),
            Self::RequireFilter(_) => "site.filter".to_string(),
            Self::MinQual(_) => "site.qual".to_string(),
            Self::Info(rule) => {
                let tag = String::from_utf8_lossy(&rule.tag);
                format!("site.info.{tag}")
            }
            Self::Expr { id, .. } => id.clone(),
        }
    }
}

/// One FORMAT tag a `--set-gts .` expression reads.
#[derive(Clone, Debug)]
pub struct FmtExprField {
    /// FORMAT key.
    pub tag: Vec<u8>,
    /// Integer column index when the tag is Integer-typed.
    pub int_column: usize,
    /// Float column index when the tag is Float-typed.
    pub float_column: Option<usize>,
}

/// Whether an expression reads a genotype-derived site field (forbidden in
/// `--set-gts .` expressions, whose mask decision must not depend on the
/// post-mask counts).
fn expr_reads_gt_fields(expr: &crate::expr::Expr) -> bool {
    use crate::expr::Expr;
    match expr {
        Expr::Field(f) => f.is_genotype_derived(),
        Expr::Bare(name) => crate::expr::Field::from_name(name)
            .is_some_and(crate::expr::Field::is_genotype_derived),
        Expr::Cmp(_, a, b) | Expr::Logic(_, a, b) => {
            expr_reads_gt_fields(a) || expr_reads_gt_fields(b)
        }
        Expr::Not(e) | Expr::Subscript(e, _) => expr_reads_gt_fields(e),
        _ => false,
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

/// A compiled genotype mask rule: either a single-field comparison or (with
/// `--set-gts .`) an expression evaluated per sample.
#[derive(Clone, Debug)]
pub struct MaskRule {
    /// A `FMT/` expression rule (`--set-gts .`, 02-cli.md §3.3): when set,
    /// the comparison fields below are unused.
    pub expr: Option<crate::expr::Expr>,
    /// Whether the expression is an `-i` (mask when false) or `-e` (mask
    /// when true); unused for field rules.
    pub include: bool,
    /// FORMAT field name (the scan key).
    pub field: Vec<u8>,
    /// The parsed comparison.
    pub cmp: crate::spec::predicate::Cmp,
    /// Missing-value behaviour (`Keep` = no mask, `Mask` = mask).
    pub on_missing: OnMissing,
    /// Position of this rule's field in `int_fields` (the scan order).
    pub int_column: usize,
    /// Position in `float_fields` when the header types the field Float.
    pub float_column: Option<usize>,
}

/// A compiled post-mask rule.
#[derive(Clone, Debug)]
pub enum TailRule {
    /// A genotype-derived bcftools expression (evaluated after masking).
    Expr {
        /// Funnel id (`after_mask.include#N` / `after_mask.exclude#N`).
        id: String,
        /// Whether the expression is an `-i` (true = keep) or `-e`.
        include: bool,
        /// The compiled predicate.
        expr: crate::expr::Expr,
    },
    /// Site missing-rate cap (denominator: kept samples).
    MaxMissing(f64),
    /// MAF lower bound.
    MinMaf(f64),
    /// MAF upper bound.
    MaxMaf(f64),
    /// Minor-allele-count lower bound.
    MinMac(u32),
    /// Minor-allele-count upper bound.
    MaxMac(u32),
    /// HWE exact-test p lower bound.
    MinHweP(f64),
}

impl TailRule {
    /// The stable funnel id (§11.5).
    #[must_use]
    pub fn id(&self) -> String {
        match self {
            Self::Expr { id, .. } => id.clone(),
            Self::MaxMissing(_) => "after_mask.missing".to_string(),
            Self::MinMaf(_) => "after_mask.maf".to_string(),
            Self::MaxMaf(_) => "after_mask.max_maf".to_string(),
            Self::MinMac(_) => "after_mask.mac".to_string(),
            Self::MaxMac(_) => "after_mask.max_mac".to_string(),
            Self::MinHweP(_) => "after_mask.hwe".to_string(),
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
    /// Float-typed FORMAT fields the mask needs scanned, in scan order (the
    /// decoder types its columns; Float-typed rules read these).
    pub float_fields: Vec<Vec<u8>>,
    /// FORMAT tags the `--set-gts .` expressions read, with the column that
    /// carries each value (the mask stage builds one value table per
    /// sample).
    pub fmt_expr_fields: Vec<FmtExprField>,
    /// Decode policies the scanner executes (ploidy from the spec; ALT
    /// indices above 1 always collapse, counted, because the site-level
    /// `multiallelic` policy is the biallelic gate below).
    pub policies: reed_core::genotype::GenotypePolicies,
    /// Stale INFO tags (bytes for the editor).
    pub stale_tags: Vec<Vec<u8>>,
    /// Recompute flags.
    pub recompute: RecomputeFlags,
    /// `-I/--no-update`: never add the automatic `AC/AN` update.
    pub no_update: bool,
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
    tags: &TagTable,
    command_line: &str,
) -> Result<FilterPlan, FilterError> {

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
    if !spec.site.exclude_variant_class.is_empty() {
        head.push(HeadRule::ExcludeVariantClass(spec.site.exclude_variant_class.clone()));
    }
    if spec.site.min_alt_alleles > 0 {
        head.push(HeadRule::MinAltAlleles(spec.site.min_alt_alleles));
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
    // ---- expressions (02-cli.md §4): non-genotype-dependent ones run here,
    // genotype-dependent ones after the mask stage ----
    let info_tags: Vec<String> = tags.info_names();
    let format_tags: Vec<String> = tags.format_names();
    for (index, text) in spec.site.include.iter().enumerate() {
        let expr = crate::expr::parse_and_check(text, &info_tags)
            .map_err(|e| FilterError::Spec(format!("-i {text:?}: {e}")))?;
        if expr.references_genotypes() {
            return Err(FilterError::Spec(format!(
                "-i {text:?} reads genotype-derived values and belongs in [site_after_mask]"
            )));
        }
        head.push(HeadRule::Expr { id: format!("site.include#{}", index + 1), include: true, expr });
    }
    for (index, text) in spec.site.exclude.iter().enumerate() {
        let expr = crate::expr::parse_and_check(text, &info_tags)
            .map_err(|e| FilterError::Spec(format!("-e {text:?}: {e}")))?;
        if expr.references_genotypes() {
            return Err(FilterError::Spec(format!(
                "-e {text:?} reads genotype-derived values and belongs in [site_after_mask]"
            )));
        }
        head.push(HeadRule::Expr { id: format!("site.exclude#{}", index + 1), include: false, expr });
    }
    let multiallelics_possible = spec.policy.multiallelic
        == crate::spec::MultiallelicPolicy::Allow
        || spec.site.max_alt_alleles > 1;
    for rule in &spec.site.info {
        let def = tags
            .info(&rule.tag)
            .ok_or_else(|| FilterError::Spec(format!("INFO tag {} is not declared in the header", rule.tag)))?;
        validate_numeric(&rule.tag, def.ty)?;
        let number_a = def.number == TagNumber::A;
        if def.number != TagNumber::Count(1) && !number_a {
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
    let mut fmt_expr_fields: Vec<FmtExprField> = Vec::new();
    let mut mask = Vec::new();
    let mut int_fields: Vec<Vec<u8>> = Vec::new();
    let mut float_fields: Vec<Vec<u8>> = Vec::new();
    for rule in &spec.genotype.mask_if {
        let def = tags
            .format(&rule.field)
            .ok_or_else(|| FilterError::Spec(format!("FORMAT field {} is not declared in the header", rule.field)))?;
        validate_numeric(&rule.field, def.ty)?;
        let is_float = def.ty == TagType::Float;
        let cmp = crate::spec::predicate::Cmp::parse(&rule.when)?;
        // The mask default is keep (§11.1): an omitted on_missing
        // deserializes to the enum default (Pass) and normalizes to Keep.
        let on_missing = match rule.on_missing {
            OnMissing::Keep | OnMissing::Pass => OnMissing::Keep,
            OnMissing::Mask | OnMissing::Fail => OnMissing::Mask,
        };
        let column_of = |fields: &mut Vec<Vec<u8>>| match fields.iter().position(|f| *f == rule.field.as_bytes()) {
            Some(i) => i,
            None => {
                fields.push(rule.field.as_bytes().to_vec());
                fields.len() - 1
            }
        };
        let (int_column, float_column) = if is_float {
            (0, Some(column_of(&mut float_fields)))
        } else {
            (column_of(&mut int_fields), None)
        };
        mask.push(MaskRule {
            expr: None,
            include: true,
            field: rule.field.as_bytes().to_vec(),
            cmp,
            on_missing,
            int_column,
            float_column,
        });
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
    if let Some(m) = spec.site_after_mask.max_maf {
        if !(0.0..=0.5).contains(&m) {
            return Err(FilterError::Spec(format!("max_maf out of range [0,0.5]: {m}")));
        }
        tail.push(TailRule::MaxMaf(m));
    }
    if spec.site_after_mask.min_mac > 0 {
        tail.push(TailRule::MinMac(spec.site_after_mask.min_mac));
    }
    if spec.site_after_mask.max_mac > 0 {
        tail.push(TailRule::MaxMac(spec.site_after_mask.max_mac));
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

    // `--set-gts .` turns the `FMT/` expressions into per-sample mask rules
    // (bcftools `filter -S .`); genotype-derived site fields stay tail rules.
    let set_gts = match spec.genotype.set_gts.as_deref() {
        None => false,
        Some(".") => true,
        Some(other) => {
            return Err(FilterError::Spec(format!(
                "--set-gts {other:?}: only `.` (the missing genotype) is supported"
            )))
        }
    };
    for (index, include, text) in spec
        .site_after_mask
        .include
        .iter()
        .enumerate()
        .map(|(i, t)| (i, true, t))
        .chain(
            spec.site_after_mask
                .exclude
                .iter()
                .enumerate()
                .map(|(i, t)| (i, false, t)),
        )
    {
        let flag = if include { "-i" } else { "-e" };
        let expr = crate::expr::parse_and_check_fmt(text, &info_tags, &format_tags)
            .map_err(|e| FilterError::Spec(format!("{flag} {text:?} (after mask): {e}")))?;
        if !expr.references_genotypes() {
            return Err(FilterError::Spec(format!(
                "{flag} {text:?} does not read genotype-derived values and belongs in [site]"
            )));
        }
        let mut fmt_tags = Vec::new();
        expr.fmt_tags(&mut fmt_tags);
        if fmt_tags.is_empty() {
            tail.push(TailRule::Expr {
                id: if include {
                    format!("after_mask.include#{}", index + 1)
                } else {
                    format!("after_mask.exclude#{}", index + 1)
                },
                include,
                expr,
            });
            continue;
        }
        if !set_gts {
            return Err(FilterError::Spec(format!(
                "{flag} {text:?} reads FORMAT values: they are masks, not site rules — add --set-gts ."
            )));
        }
        if expr_reads_gt_fields(&expr) {
            return Err(FilterError::Spec(format!(
                "{flag} {text:?} mixes FORMAT values with genotype-derived fields; \
                 --set-gts . expressions may read FORMAT and site fields only"
            )));
        }
        for tag in &fmt_tags {
            let def = tags
                .format(tag)
                .ok_or_else(|| FilterError::Spec(format!("FORMAT tag {tag} is not declared in the header")))?;
            if def.number != TagNumber::Count(1) {
                return Err(FilterError::Spec(format!(
                    "FORMAT tag {tag} is multi-valued (Number={:?}): subscripts are not supported yet", def.number
                )));
            }
            validate_numeric(tag, def.ty)?;
            let is_float = def.ty == TagType::Float;
            let column_of = |fields: &mut Vec<Vec<u8>>| match fields.iter().position(|f| f.as_slice() == tag.as_bytes()) {
                Some(i) => i,
                None => {
                    fields.push(tag.as_bytes().to_vec());
                    fields.len() - 1
                }
            };
            let (int_column, float_column) = if is_float {
                (0, Some(column_of(&mut float_fields)))
            } else {
                (column_of(&mut int_fields), None)
            };
            fmt_expr_fields.push(FmtExprField {
                tag: tag.as_bytes().to_vec(),
                int_column,
                float_column,
            });
        }
        mask.push(MaskRule {
            expr: Some(expr),
            include,
            field: Vec::new(),
            cmp: crate::spec::Cmp::parse(">= 0").expect("invariant: a constant comparison"),
            on_missing: crate::spec::OnMissing::Keep,
            int_column: 0,
            float_column: None,
        });
    }

    // ---- strategy (§11.4; SketchThenEmit selection is M5b) ----
    // The sketch path evaluates tail rules in the mid stage, where the
    // record bytes are gone; a genotype-derived expression therefore forces
    // the single-pass path (02-cli.md §4: expressions run on the record).
    let has_tail_expr = tail.iter().any(|r| matches!(r, TailRule::Expr { .. }));
    if let Some(m) = spec.mind
        && !(0.0..=1.0).contains(&m)
    {
        return Err(FilterError::Spec(format!("--mind out of range [0,1]: {m}")));
    }
    let strategy = if has_tail_expr {
        ExecStrategy::SinglePass
    } else if spec.mind.is_some() {
        ExecStrategy::TwoPassReparse
    } else {
        ExecStrategy::SinglePass
    };

    // ---- recompute + stale ----
    let no_update = spec.recompute.no_update;
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

    let policies = reed_core::genotype::GenotypePolicies {
        ploidy: match spec.policy.ploidy {
            crate::spec::PloidyPolicy::RequireDiploid => {
                reed_core::genotype::PloidyPolicy::RequireDiploid
            }
            crate::spec::PloidyPolicy::AllowHaploid => {
                reed_core::genotype::PloidyPolicy::HaploidAsHalfCall
            }
        },
        multiallelic: reed_core::genotype::MultiallelicPolicy::CollapseAlt,
        ..reed_core::genotype::GenotypePolicies::default()
    };

    Ok(FilterPlan {
        head,
        mask,
        tail,
        strategy,
        int_fields,
        float_fields,
        fmt_expr_fields,
        policies,
        stale_tags,
        recompute,
        no_update,
        provenance,
    })
}

fn validate_numeric(tag: &str, ty: TagType) -> Result<(), FilterError> {
    if ty.is_numeric() {
        Ok(())
    } else {
        Err(FilterError::Spec(format!(
            "tag {tag} is {ty:?}-typed: numeric rules need Integer or Float"
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
