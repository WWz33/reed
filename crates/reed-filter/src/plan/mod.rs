//! The plan compiler (§11.3): `compile(spec, header, budget) -> FilterPlan`
//! is the only input of execution — validation (tags declared, numeric
//! types, Number ∈ {1, A}), rule expansion, strategy choice (§11.4),
//! ascertainment derivation (§11.7) and provenance (§11.8) all happen here.

// Casts in this file are bounded by construction: The `-q/-Q` threshold is deliberately the f32 bcftools parses.
#![allow(clippy::cast_possible_truncation)]

pub mod provenance;
pub mod strategy;

use reed_core::tag_table::{TagNumber, TagTable, TagType};

pub use strategy::ExecStrategy;

use crate::FilterError;
use crate::spec::region::RegionTriple;
use crate::spec::{FilterSpec, OnMissing, Quantifier, RegionRule};

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
    /// The default biallelic gate: records with more than one ALT are
    /// dropped unless the spec allows multiallelic input (§11.6,
    /// [`FilterSpec::multiallelics_allowed`]). Pushed by `compile` when
    /// `biallelic_gate` is set and no explicit `max_alt_alleles` covers it.
    BiallelicGate,
    /// The FILTER column must be one of these (`.` never passes).
    RequireFilter(Vec<Vec<u8>>),
    /// QUAL lower bound (missing QUAL fails, bcftools semantics).
    MinQual(f32),
    /// Non-reference allele frequency lower bound on the INFO AC/AN values
    /// (bcftools `-q` when the header declares both tags); records without
    /// usable INFO pass and are re-checked by the tail rule.
    MinAfInfo(f64),
    /// See [`HeadRule::MinAfInfo`].
    MaxAfInfo(f64),
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
            Self::BiallelicGate => "site.biallelic_gate".to_string(),
            Self::RequireFilter(_) => "site.filter".to_string(),
            Self::MinQual(_) => "site.qual".to_string(),
            Self::MinAfInfo(_) => "site.min_af".to_string(),
            Self::MaxAfInfo(_) => "site.max_af".to_string(),
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

/// Collects the built-in fields an expression reads (04-cleanup §2.3: the
/// stale-INFO warning needs the names).
fn expr_fields(expr: &crate::expr::Expr, out: &mut Vec<crate::expr::Field>) {
    use crate::expr::Expr;
    match expr {
        Expr::Field(f) => out.push(*f),
        Expr::Bare(name) => {
            if let Some(f) = crate::expr::Field::from_name(name) {
                out.push(f);
            }
        }
        Expr::Cmp(_, a, b) | Expr::Logic(_, a, b) => {
            expr_fields(a, out);
            expr_fields(b, out);
        }
        Expr::Not(e) | Expr::Subscript(e, _) => expr_fields(e, out),
        _ => {}
    }
}

/// Warns when an expression reads a built-in that the header also declares as
/// an INFO tag: the INFO value may be stale (this dataset's trap), while the
/// recomputed `--maf/--mac` path always uses the genotypes.
fn stale_info_hits(expr: &crate::expr::Expr, tags: &TagTable) -> Vec<&'static str> {
    let mut fields = Vec::new();
    expr_fields(expr, &mut fields);
    fields
        .into_iter()
        .filter(|field| {
            field.is_genotype_derived()
                && !matches!(
                    field,
                    crate::expr::Field::NMissing | crate::expr::Field::FMissing
                )
                && tags.info(field.name()).is_some()
        })
        .map(crate::expr::Field::name)
        .collect()
}

fn warn_stale_info(expr: &crate::expr::Expr, text: &str, tags: &TagTable) {
    for name in stale_info_hits(expr, tags) {
        tracing::warn!(
            "{text:?} reads {name}, which the header declares as \
             INFO/{name}; the INFO value may be stale — use --maf/--mac (recomputed from \
             genotypes) or remove the tag with -x INFO/{name}"
        );
    }
}

/// `-q/-Q` reads `INFO/AC,AN` when the header declares both and no cohort
/// change forces the genotype path (doc_v12 D12-6; cli.md:78 promises this
/// warning). The dataset trap is the same as the expressions': the INFO
/// counts may describe a different cohort than the genotypes do.
fn warn_stale_af(tags: &TagTable) {
    let _ = tags;
    tracing::warn!(
        "-q/-Q reads INFO/AC,INFO/AN, which may be stale; use --maf/--mac for \
         genotype-recomputed values or remove the tags with -x INFO/AC,INFO/AN"
    );
}

/// Whether an expression reads a genotype-derived site field (forbidden in
/// `--set-gts .` expressions, whose mask decision must not depend on the
/// post-mask counts).
fn expr_reads_gt_fields(expr: &crate::expr::Expr) -> bool {
    use crate::expr::Expr;
    match expr {
        Expr::Field(f) => f.is_genotype_derived(),
        Expr::Bare(name) => {
            crate::expr::Field::from_name(name).is_some_and(crate::expr::Field::is_genotype_derived)
        }
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
    /// Allele-frequency lower bound (bcftools `-q`), with the `:type` the
    /// threshold applies to.
    MinAf(f64, crate::spec::AfKind),
    /// Allele-frequency upper bound (bcftools `-Q`), with its `:type`.
    MaxAf(f64, crate::spec::AfKind),
    /// Site missing-rate cap (`--geno`: `2N − an ≤ X · 2N`; exact integer
    /// judgment, doc_v4 02 回退版 §4).
    MaxMissing(crate::spec::MissingRate),
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
    /// One genotype-class site rule (bcftools `-g/--genotype`); a record
    /// carrying several must satisfy every one of them.
    GtClass(crate::spec::GtClassRule),
}

impl TailRule {
    /// The stable funnel id (§11.5).
    #[must_use]
    pub fn id(&self) -> String {
        match self {
            Self::Expr { id, .. } => id.clone(),
            Self::MinAf(..) => "after_mask.min_af".to_string(),
            Self::MaxAf(..) => "after_mask.max_af".to_string(),
            Self::MaxMissing(_) => crate::funnel::AFTER_MASK_MISSING.to_string(),
            Self::MinMaf(_) => "after_mask.maf".to_string(),
            Self::MaxMaf(_) => "after_mask.max_maf".to_string(),
            Self::MinMac(_) => "after_mask.mac".to_string(),
            Self::MaxMac(_) => "after_mask.max_mac".to_string(),
            Self::MinHweP(_) => "after_mask.hwe".to_string(),
            Self::GtClass(rule) => format!("after_mask.gt_class.{}", rule.id()),
        }
    }
}

/// The recomputed-tag flags (compiled from [`RecomputeSpec`](crate::spec::RecomputeSpec)).
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
    /// indices above 1 always collapse and are counted — the decoder never
    /// rejects an allele index. Records with several ALTs are handled
    /// upstream: the default [`HeadRule::BiallelicGate`] drops them unless
    /// the spec allows multiallelic input).
    pub policies: reed_core::genotype::GenotypePolicies,
    /// Stale INFO tags (bytes for the editor).
    pub stale_tags: Vec<Vec<u8>>,
    /// The subset of [`Self::stale_tags`] the driver's automatic AC/AN
    /// update added (not the user's `-x`). A record without a GT key keeps
    /// these values instead of being zeroed (doc_v12 D12-2); user-requested
    /// removals are unconditional.
    pub auto_stale_tags: Vec<Vec<u8>>,
    /// Recompute flags.
    pub recompute: RecomputeFlags,
    /// `-I/--no-update`: never add the automatic `AC/AN` update.
    pub no_update: bool,
    /// Whether the record's `INFO/AC,AN` still describe the cohort the tail
    /// rules decide on. The driver clears this when the sample set changes
    /// (a static subset, `--mind`); the tail stage then ignores the INFO
    /// values and uses the genotype counts, like bcftools
    /// (`vcfview.c:341-350`, `filter.c:1471`; doc_v12 D12-6).
    pub info_trusted: bool,
    /// The header's declared INFO/FORMAT tags: the runtime evaluator needs
    /// the declaration to read a Flag as 1/0 by presence and a multi-valued
    /// tag as a vector (doc_v12 D12-7). `Arc` so a worker's evaluation
    /// context can borrow it without cloning the map.
    pub tags: std::sync::Arc<TagTable>,
    /// Provenance: canonical spec JSON + sha256 + command line.
    pub provenance: provenance::Provenance,
}

impl FilterPlan {
    /// Whether the tail stage reads the record's own bytes.
    ///
    /// The sketch path evaluates the tail rules in the mid stage, where the
    /// record bytes are gone: `SiteValues::from_counts` carries `info = "."`,
    /// `qual = None` and dummy line fields (`expr/mod.rs`). A tail rule that
    /// reads any of those therefore cannot be evaluated there — the strategy
    /// selector (the driver re-picks after `--mind` resolves, doc_v12 D12-1)
    /// must not choose [`ExecStrategy::SketchThenEmit`] while this is true.
    ///
    /// Two kinds of rule qualify:
    /// - `TailRule::Expr` — an expression may read any line field;
    /// - `TailRule::MinAf`/`MaxAf` — `af_of` prefers the record's INFO/AC,AN
    ///   when they are present, so the mid stage (no INFO) would silently
    ///   fall back to the genotype counts and reach a different verdict.
    ///
    /// This is the single source of truth for the guard: [`compile`] uses it
    /// to pick the default strategy, and the driver consults it before any
    /// override.
    #[must_use]
    pub fn tail_reads_record(&self) -> bool {
        self.tail.iter().any(|r| {
            matches!(
                r,
                TailRule::Expr { .. } | TailRule::MinAf(..) | TailRule::MaxAf(..)
            )
        })
    }
}

/// Compiles the spec against the input header.
///
/// `biallelic_gate` is the filter driver's default (§11.6, plink2-style):
/// push [`HeadRule::BiallelicGate`] unless the spec allows multiallelic
/// input. The stats driver passes `false` — its multiallelic handling is
/// per reporter, not a record gate.
///
/// # Errors
///
/// [`FilterError::Spec`] on validation failures: unknown tags, non-numeric
/// types, `Number` outside {1, A}, a missing quantifier where one is
/// required, thresholds out of range, or the M7 `midp` HWE mode.
pub fn compile(
    spec: &FilterSpec,
    tags: &TagTable,
    command_line: &str,
    biallelic_gate: bool,
) -> Result<FilterPlan, FilterError> {
    // ---- head rules (§11.2 order: targets → class/alt → FILTER → QUAL →
    // INFO in declaration order) ----
    let mut head = Vec::new();
    let regions = compile_regions(&spec.targets)?;
    if let Some(r) = regions {
        head.push(HeadRule::Regions {
            include: r.0,
            exclude: r.1,
        });
    }
    if !spec.site.variant_class.is_empty() {
        head.push(HeadRule::VariantClass(spec.site.variant_class.clone()));
    }
    if !spec.site.exclude_variant_class.is_empty() {
        head.push(HeadRule::ExcludeVariantClass(
            spec.site.exclude_variant_class.clone(),
        ));
    }
    if spec.site.min_alt_alleles > 0 {
        head.push(HeadRule::MinAltAlleles(spec.site.min_alt_alleles));
    }
    if let Some(max) = spec.site.max_alt_alleles {
        head.push(HeadRule::MaxAltAlleles(max));
    } else if biallelic_gate && !spec.multiallelics_allowed() {
        // The default gate (§11.6, plink2-style): without `multiallelic =
        // "allow"` or an explicit `max_alt_alleles > 1`, multiallelic
        // records do not reach the decoder — the collapsed per-record
        // counts (and the recomputed AC/MAF/HWE) are biallelic-valid only.
        head.push(HeadRule::BiallelicGate);
    }
    if !spec.site.require_filter.is_empty() {
        head.push(HeadRule::RequireFilter(
            spec.site
                .require_filter
                .iter()
                .map(String::as_bytes)
                .map(<[u8]>::to_vec)
                .collect(),
        ));
    }
    if let Some(q) = spec.site.min_qual {
        if q < 0.0 {
            return Err(FilterError::Spec(format!(
                "min_qual out of range [0,inf): {q}"
            )));
        }
        head.push(HeadRule::MinQual(q));
    }
    // `-q/-Q` (bcftools): the INFO AC/AN values decide at the site stage when
    // the header declares both AND the threshold is the `nref` form; the
    // other `:type`s are genotype-derived and live in the tail only.
    //
    // The INFO preference holds only while AC/AN still describe the cohort
    // being decided on. Once the sample set changes (a static subset, a
    // genotype mask, `--mind`, or `--set-gts .`), bcftools recomputes from
    // the subset GT before applying `-q/-Q` (vcfview.c:341-350) and gates
    // its own INFO preference on `!sample_mask_set` (filter.c:1471); reed
    // mirrors that by not pushing the INFO short-circuit at all, so the
    // decision runs through `af_of`, whose `info_trusted` flag is false in
    // every execution path once the cohort changed (doc_v12 D12-6).
    let info_af = tags.info("AC").is_some() && tags.info("AN").is_some();
    let cohort_changed = spec.sample_select.is_some()
        || spec.mind.is_some()
        || !spec.genotype.mask_if.is_empty()
        || spec.genotype.set_gts.is_some();
    if info_af && !cohort_changed {
        warn_stale_af(tags);
    }
    for bound in spec.site_after_mask.min_af.iter().copied() {
        if !(0.0..=1.0).contains(&bound) {
            return Err(FilterError::Spec(format!(
                "min_af out of range [0,1]: {bound}"
            )));
        }
        if info_af
            && !cohort_changed
            && spec
                .site_after_mask
                .min_af_kind
                .unwrap_or(crate::spec::AfKind::Nref)
                == crate::spec::AfKind::Nref
        {
            head.push(HeadRule::MinAfInfo(af_threshold(bound)));
        }
    }
    for bound in spec.site_after_mask.max_af.iter().copied() {
        if !(0.0..=1.0).contains(&bound) {
            return Err(FilterError::Spec(format!(
                "max_af out of range [0,1]: {bound}"
            )));
        }
        if info_af
            && !cohort_changed
            && spec
                .site_after_mask
                .max_af_kind
                .unwrap_or(crate::spec::AfKind::Nref)
                == crate::spec::AfKind::Nref
        {
            head.push(HeadRule::MaxAfInfo(af_threshold(bound)));
        }
    }
    // ---- expressions (02-cli.md §4): non-genotype-dependent ones run here,
    // genotype-dependent ones after the mask stage ----
    let format_tags: Vec<String> = tags.format_names();
    for (index, text) in spec.site.include.iter().enumerate() {
        let expr = crate::expr::parse_and_check(text, tags)
            .map_err(|e| FilterError::Spec(format!("-i {text:?}: {e}")))?;
        warn_stale_info(&expr, text, tags);
        if expr.references_genotypes() {
            return Err(FilterError::Spec(format!(
                "-i {text:?} reads genotype-derived values and belongs in [site_after_mask]"
            )));
        }
        head.push(HeadRule::Expr {
            id: format!("site.include#{}", index + 1),
            include: true,
            expr,
        });
    }
    for (index, text) in spec.site.exclude.iter().enumerate() {
        let expr = crate::expr::parse_and_check(text, tags)
            .map_err(|e| FilterError::Spec(format!("-e {text:?}: {e}")))?;
        warn_stale_info(&expr, text, tags);
        if expr.references_genotypes() {
            return Err(FilterError::Spec(format!(
                "-e {text:?} reads genotype-derived values and belongs in [site_after_mask]"
            )));
        }
        head.push(HeadRule::Expr {
            id: format!("site.exclude#{}", index + 1),
            include: false,
            expr,
        });
    }
    for rule in &spec.site.info {
        let def = tags.info(&rule.tag).ok_or_else(|| {
            FilterError::Spec(format!(
                "INFO tag {} is not declared in the header",
                rule.tag
            ))
        })?;
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
        let def = tags.format(&rule.field).ok_or_else(|| {
            FilterError::Spec(format!(
                "FORMAT field {} is not declared in the header",
                rule.field
            ))
        })?;
        validate_numeric(&rule.field, def.ty)?;
        let is_float = def.ty == TagType::Float;
        let cmp = crate::spec::predicate::Cmp::parse(&rule.when)?;
        // The mask default is keep (§11.1): an omitted on_missing
        // deserializes to the enum default (Pass) and normalizes to Keep.
        let on_missing = match rule.on_missing {
            OnMissing::Keep | OnMissing::Pass => OnMissing::Keep,
            OnMissing::Mask | OnMissing::Fail => OnMissing::Mask,
        };
        let column_of = |fields: &mut Vec<Vec<u8>>| match fields
            .iter()
            .position(|f| *f == rule.field.as_bytes())
        {
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
    // bcftools applies -q/-Q (site-level AF bounds) before the plink2-style
    // missing/MAF blocks.
    if let Some(f) = spec.site_after_mask.min_af {
        if !(0.0..=1.0).contains(&f) {
            return Err(FilterError::Spec(format!("min_af out of range [0,1]: {f}")));
        }
        tail.push(TailRule::MinAf(
            af_threshold(f),
            spec.site_after_mask
                .min_af_kind
                .unwrap_or(crate::spec::AfKind::Nref),
        ));
    }
    if let Some(f) = spec.site_after_mask.max_af {
        if !(0.0..=1.0).contains(&f) {
            return Err(FilterError::Spec(format!("max_af out of range [0,1]: {f}")));
        }
        tail.push(TailRule::MaxAf(
            af_threshold(f),
            spec.site_after_mask
                .max_af_kind
                .unwrap_or(crate::spec::AfKind::Nref),
        ));
    }
    if let Some(rate) = spec.site_after_mask.max_missing {
        // The range is validated where the rate is parsed
        // (spec/missing_rate.rs), so a compiled spec can only carry
        // `q ≤ 10^9`.
        tail.push(TailRule::MaxMissing(rate));
    }
    if let Some(m) = spec.site_after_mask.min_maf {
        if !(0.0..=0.5).contains(&m) {
            return Err(FilterError::Spec(format!(
                "min_maf out of range [0,0.5]: {m}"
            )));
        }
        tail.push(TailRule::MinMaf(m));
    }
    if let Some(m) = spec.site_after_mask.max_maf {
        if !(0.0..=0.5).contains(&m) {
            return Err(FilterError::Spec(format!(
                "max_maf out of range [0,0.5]: {m}"
            )));
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
            return Err(FilterError::Spec(format!(
                "min_hwe_p out of range [0,1]: {p}"
            )));
        }
        tail.push(TailRule::MinHweP(p));
    }
    // `-g/--genotype` (repeated flags AND, unlike bcftools' last-wins):
    // exact duplicates collapse; requiring a class and its absence is
    // unsatisfiable — reject at compile time rather than silently drop
    // every record.
    let mut gt_rules: Vec<crate::spec::GtClassRule> =
        Vec::with_capacity(spec.site_after_mask.gt_class.len());
    for &rule in &spec.site_after_mask.gt_class {
        if gt_rules.contains(&rule) {
            continue; // an exact duplicate: AND is idempotent
        }
        if let Some(prior) = gt_rules.iter().find(|r| r.class() == rule.class()) {
            return Err(FilterError::Spec(format!(
                "-g/--genotype rules {:?} and {:?} are contradictory \
                 (a site cannot have and lack a {} call)",
                prior.id(),
                rule.id(),
                rule.class(),
            )));
        }
        gt_rules.push(rule);
    }
    tail.extend(gt_rules.into_iter().map(TailRule::GtClass));

    // `--set-gts .` turns the `FMT/` expressions into per-sample mask rules
    // (bcftools `filter -S .`); genotype-derived site fields stay tail rules.
    let set_gts = match spec.genotype.set_gts.as_deref() {
        None => false,
        Some(".") => true,
        Some(other) => {
            return Err(FilterError::Spec(format!(
                "--set-gts {other:?}: only `.` (the missing genotype) is supported"
            )));
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
        let expr = crate::expr::parse_and_check_fmt(text, tags, &format_tags)
            .map_err(|e| FilterError::Spec(format!("{flag} {text:?} (after mask): {e}")))?;
        warn_stale_info(&expr, text, tags);
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
            let def = tags.format(tag).ok_or_else(|| {
                FilterError::Spec(format!("FORMAT tag {tag} is not declared in the header"))
            })?;
            if def.number != TagNumber::Count(1) {
                return Err(FilterError::Spec(format!(
                    "FORMAT tag {tag} is multi-valued (Number={:?}): subscripts are not supported yet",
                    def.number
                )));
            }
            validate_numeric(tag, def.ty)?;
            let is_float = def.ty == TagType::Float;
            let column_of = |fields: &mut Vec<Vec<u8>>| match fields
                .iter()
                .position(|f| f.as_slice() == tag.as_bytes())
            {
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
    // record bytes are gone; a tail rule that reads them therefore forces
    // the single-pass path (02-cli.md §4: expressions run on the record).
    // The predicate mirrors `FilterPlan::tail_reads_record` (the plan does
    // not exist yet here), so the driver's later override asks the same
    // question (doc_v12 D12-1).
    let tail_reads_record = tail.iter().any(|r| {
        matches!(
            r,
            TailRule::Expr { .. } | TailRule::MinAf(..) | TailRule::MaxAf(..)
        )
    });
    let strategy = if tail_reads_record {
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
        auto_stale_tags: Vec::new(),
        recompute,
        no_update,
        info_trusted: true,
        tags: std::sync::Arc::new(tags.clone()),
        provenance,
    })
}

/// bcftools parses `-q/-Q` with `sscanf("%f")` into a **float**, so its
/// threshold is the f32 nearest the written value; at an exact boundary
/// (e.g. AC=34, AN=680 against 0.05) that rounds up and drops the site.
/// Quantising here reproduces the comparison bit for bit.
fn af_threshold(value: f64) -> f64 {
    f64::from(value as f32)
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

/// Parsed include/exclude region lists.
type RegionLists = (Vec<RegionTriple>, Vec<RegionTriple>);

fn compile_regions(rule: &RegionRule) -> Result<Option<RegionLists>, FilterError> {
    if rule.is_empty() {
        return Ok(None);
    }
    use crate::spec::region::parse_regions;
    Ok(Some((
        parse_regions(&rule.include)?,
        parse_regions(&rule.exclude)?,
    )))
}

#[cfg(test)]
mod gt_class_tests {
    use super::*;
    use crate::spec::GtClassRule as R;

    fn compile_gt(rules: Vec<R>) -> Result<FilterPlan, FilterError> {
        let spec = FilterSpec {
            site_after_mask: crate::spec::TailRules {
                gt_class: rules,
                ..crate::spec::TailRules::default()
            },
            ..FilterSpec::default()
        };
        compile(&spec, &TagTable::new(), "", true)
    }

    /// Repeated flags AND: exact duplicates collapse to one rule, so the
    /// funnel never counts the same predicate twice.
    #[test]
    fn duplicate_gt_rules_dedupe() {
        let plan = compile_gt(vec![R::Het, R::Het]).expect("invariant: compiles");
        let gt: Vec<_> = plan
            .tail
            .iter()
            .filter(|r| matches!(r, TailRule::GtClass(_)))
            .collect();
        assert_eq!(gt.len(), 1);
        assert_eq!(gt[0].id(), "after_mask.gt_class.het");
        let plan = compile_gt(vec![R::Het, R::Miss, R::Het, R::Miss])
            .expect("invariant: compiles");
        let ids: Vec<_> = plan.tail.iter().map(TailRule::id).collect();
        assert_eq!(ids, ["after_mask.gt_class.het", "after_mask.gt_class.miss"]);
    }

    /// Requiring a class and its absence is unsatisfiable: rejected at
    /// compile time, not silently true for no record.
    #[test]
    fn opposite_polarities_are_a_compile_error() {
        for rules in [
            vec![R::Hom, R::NoHom],
            vec![R::NoHet, R::Het],
            vec![R::Miss, R::NoMiss],
            vec![R::Het, R::Het, R::NoHet],
        ] {
            let err = compile_gt(rules).expect_err("must reject");
            let FilterError::Spec(text) = &err else {
                panic!("expected a spec error, got {err:?}");
            };
            assert!(text.contains("contradictory"), "{text}");
        }
    }
}

#[cfg(test)]
mod stale_info_tests {
    use super::*;
    use crate::expr;

    /// 04 §2.3: the warning lists exactly the INFO-declared genotype-derived
    /// fields an expression reads.
    #[test]
    fn stale_info_warning_targets_only_declared_tags() {
        let mut tags = TagTable::new();
        tags.add_info("AC", TagType::Integer, TagNumber::A);
        tags.add_info("DP", TagType::Integer, TagNumber::Count(1));

        for (text, expected) in [
            ("AC<5", vec!["AC"]),
            ("AC<5 || MAF>0.01", vec!["AC"]),
            ("Maf > 0.01", vec![]),
            ("DP>10", vec![]),
        ] {
            let expr = expr::parse(text).expect("invariant: parses");
            assert_eq!(stale_info_hits(&expr, &tags), expected, "{text}");
        }
    }
}
