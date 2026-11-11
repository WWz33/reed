//! The plan compiler (§11.3): `compile(spec, header, budget) -> FilterPlan`
//! is the only input of execution — validation (tags declared, numeric
//! types, Number ∈ {1, A}), rule expansion, strategy choice (§11.4),
//! ascertainment derivation (§11.7) and provenance (§11.8) all happen here.

// Casts in this file are bounded by construction: The `-q/-Q` threshold is deliberately the f32 bcftools parses.
#![allow(clippy::cast_possible_truncation)]

pub use reed_core::provenance;
mod compile;
pub use compile::compile;
pub mod strategy;

use reed_core::tag_table::{TagTable, TagType};

pub use strategy::ExecStrategy;

use crate::FilterError;
use crate::spec::region::RegionTriple;
use crate::spec::{OnMissing, Quantifier, RegionRule};
// The module's tests build specs directly; the compiler itself now lives in
// `plan/compile/` and imports it there.
#[cfg(test)]
use crate::spec::FilterSpec;
#[cfg(test)]
use reed_core::tag_table::TagNumber;

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
    /// Every REF/ALT allele must be one of A/C/G/T (plink2
    /// `--snps-only just-acgt`; doc_v21 S29). Applied after `VariantClass(snv)`.
    JustAcgt,
    /// The record's ALT count must be at least this.
    MinAltAlleles(u32),
    /// The record's ALT count must not exceed this.
    MaxAltAlleles(u32),
    /// The default biallelic gate: records with more than one ALT are
    /// dropped unless the spec allows multiallelic input (§11.6,
    /// [`FilterSpec::multiallelics_allowed`]). Pushed by `compile` when
    /// `biallelic_gate` is set and no explicit `max_alt_alleles` covers it.
    BiallelicGate,
    /// The FILTER column must contain one of these as a `;`-separated
    /// value (bcftools `-f`; the unannotated `.` column passes only when
    /// `.` is listed).
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
            Self::JustAcgt => "site.just_acgt".to_string(),
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

/// Which `-q/-Q` bounds read `INFO/AC,AN` — `(min, max)`. A bound does when
/// the header declares both tags, its threshold is the `nref` form (the
/// `:type` variants are genotype-derived), and no cohort change forces the
/// genotype path (doc_v12 D12-6; cli.md §3.2 promises the warning 取到时 —
/// when the INFO value is taken, not whenever the header declares the
/// trap: every GATK header declares AC+AN). The dataset trap is the same
/// as the expressions': the INFO counts may describe a different cohort
/// than the genotypes do. The `HeadRule` pushes and the warning both read
/// this one decision, so they cannot drift apart.
fn af_info_sides(
    tail: &crate::spec::TailRules,
    info_af: bool,
    cohort_changed: bool,
) -> (bool, bool) {
    let nref_form = |bound: Option<f64>, kind: Option<crate::spec::AfKind>| {
        bound.is_some() && kind.unwrap_or(crate::spec::AfKind::Nref) == crate::spec::AfKind::Nref
    };
    if !info_af || cohort_changed {
        return (false, false);
    }
    (
        nref_form(tail.min_af, tail.min_af_kind),
        nref_form(tail.max_af, tail.max_af_kind),
    )
}

fn warn_stale_af() {
    tracing::warn!(
        "-q/-Q reads INFO/AC,INFO/AN, which may be stale; use --maf/--mac for \
         genotype-recomputed values or remove the tags with -x INFO/AC,INFO/AN"
    );
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

/// One planned integer FORMAT column: the key and its explicit fallback.
///
/// The fallback is a policy the consumer declares, never inferred from the
/// key's name: `--min-dp/--max-dp` set `fallback = Some("AD")` (a record
/// whose FORMAT has no `DP` key uses the sample's AD sum), while an explicit
/// `FMT/DP` expression keeps `fallback = None` and therefore stays literal
/// (bcftools parity, doc_v22 D22-1).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IntField {
    /// The FORMAT key read into this column.
    pub name: Vec<u8>,
    /// Key summed when `name` is absent from a record's FORMAT.
    pub fallback: Option<Vec<u8>>,
}

impl From<&IntField> for reed_io::vcf::format_layout::IntSpec {
    fn from(field: &IntField) -> Self {
        Self::owned(field.name.clone(), field.fallback.clone())
    }
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
    /// Integer FORMAT columns the mask needs scanned, in scan order (each
    /// carries its own fallback policy — doc_v22 D22-1).
    pub int_fields: Vec<IntField>,
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

    /// Whether any stage reads the record's INFO column: typed head rules,
    /// `-i/-e` expressions (head or `--set-gts .`) and the AC/AN readers all
    /// need the same [`InfoIndex`](reed_io::vcf::info::InfoIndex).
    ///
    /// The filter pipeline asks this once per plan: when false, the line
    /// classifier keeps its allocation-free duplicate check and no index is
    /// built at all (doc_v28 V28-13).
    #[must_use]
    pub fn reads_info(&self) -> bool {
        let head = self.head.iter().any(|r| {
            matches!(
                r,
                HeadRule::Info(_)
                    | HeadRule::Expr { .. }
                    | HeadRule::MinAfInfo(_)
                    | HeadRule::MaxAfInfo(_)
            )
        });
        let tail = self.tail.iter().any(|r| {
            matches!(
                r,
                TailRule::Expr { .. } | TailRule::MinAf(..) | TailRule::MaxAf(..)
            )
        });
        head || tail || self.mask.iter().any(|r| r.expr.is_some())
    }
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

fn compile_regions(
    rule: &RegionRule,
    contigs: &[String],
) -> Result<Option<RegionLists>, FilterError> {
    if rule.is_empty() {
        return Ok(None);
    }
    use crate::spec::region::parse_regions_with;
    Ok(Some((
        parse_regions_with(&rule.include, contigs)?,
        parse_regions_with(&rule.exclude, contigs)?,
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
        let plan = compile_gt(vec![R::Het, R::Miss, R::Het, R::Miss]).expect("invariant: compiles");
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

    /// cli.md §3.2: the `-q/-Q` stale-INFO warning fires only when a
    /// threshold actually takes the INFO path — not whenever the header
    /// declares AC+AN (every GATK header does).
    #[test]
    fn af_warning_fires_only_when_the_info_path_is_taken() {
        fn tail(
            min_af: Option<f64>,
            min_kind: Option<crate::spec::AfKind>,
            max_af: Option<f64>,
            max_kind: Option<crate::spec::AfKind>,
        ) -> crate::spec::TailRules {
            crate::spec::TailRules {
                min_af,
                min_af_kind: min_kind,
                max_af,
                max_af_kind: max_kind,
                ..Default::default()
            }
        }
        // A `nref` threshold, an AC+AN header, the cohort intact: taken —
        // the default kind is `nref`, and either bound alone takes it.
        assert_eq!(
            af_info_sides(&tail(Some(0.05), None, None, None), true, false),
            (true, false)
        );
        assert_eq!(
            af_info_sides(&tail(None, None, Some(0.95), None), true, false),
            (false, true)
        );
        assert_eq!(
            af_info_sides(
                &tail(
                    Some(0.05),
                    Some(crate::spec::AfKind::Nref),
                    Some(0.95),
                    Some(crate::spec::AfKind::Nref)
                ),
                true,
                false
            ),
            (true, true)
        );
        // No threshold at all (a plain `-e 'QD<2 || ...'` run): not taken.
        assert_eq!(
            af_info_sides(&tail(None, None, None, None), true, false),
            (false, false)
        );
        // A typed threshold never reads INFO, on either bound.
        assert_eq!(
            af_info_sides(
                &tail(Some(0.05), Some(crate::spec::AfKind::Alt1), None, None),
                true,
                false
            ),
            (false, false)
        );
        assert_eq!(
            af_info_sides(
                &tail(None, None, Some(0.95), Some(crate::spec::AfKind::Alt1)),
                true,
                false
            ),
            (false, false)
        );
        // A cohort change (`--set-gts .` and friends) forces the GT path.
        assert_eq!(
            af_info_sides(&tail(Some(0.05), None, Some(0.95), None), true, true),
            (false, false)
        );
        // A header without both AC and AN.
        assert_eq!(
            af_info_sides(&tail(Some(0.05), None, Some(0.95), None), false, false),
            (false, false)
        );
    }
    /// doc_v23 V23-25: an empty `sample_select` is not a cohort change.
    /// `{"sample_select": {}}` used to switch `-q` from the INFO AC/AN
    /// short-circuit to genotype-recomputed values.
    #[test]
    fn empty_sample_select_keeps_the_info_af_short_circuit() {
        let mut tags = TagTable::new();
        tags.add_info("AC", TagType::Integer, TagNumber::A);
        tags.add_info("AN", TagType::Integer, TagNumber::Count(1));
        let request = |select: Option<crate::spec::SampleSelect>| FilterSpec {
            sample_select: select,
            site_after_mask: crate::spec::TailRules {
                min_af: Some(0.1),
                ..Default::default()
            },
            ..Default::default()
        };

        let empty = compile(
            &request(Some(crate::spec::SampleSelect {
                keep: None,
                remove: None,
            })),
            &tags,
            "",
            true,
        )
        .expect("invariant: compiles");
        assert!(
            empty
                .head
                .iter()
                .any(|rule| matches!(rule, HeadRule::MinAfInfo(_))),
            "an empty selection must keep the INFO short-circuit"
        );

        let real = compile(
            &request(Some(crate::spec::SampleSelect {
                keep: Some(crate::spec::SampleList::Names(vec!["s0".into()])),
                remove: None,
            })),
            &tags,
            "",
            true,
        )
        .expect("invariant: compiles");
        assert!(
            !real
                .head
                .iter()
                .any(|rule| matches!(rule, HeadRule::MinAfInfo(_))),
            "a real selection must recompute from the subset"
        );
    }
}

/// doc_v22 D22-1: the depth fallback is an explicit policy on the mask rule,
/// never inferred from the name `DP`; an explicit `FMT/DP` expression stays
/// literal and needs the tag declared (bcftools parity).
#[cfg(test)]
mod depth_policy_tests {
    use super::*;
    use crate::spec::{GenotypeMaskRule, GenotypeRules, TailRules};

    fn tags(with_dp: bool) -> TagTable {
        let mut table = TagTable::new();
        table.add_format("GT", TagType::Character, TagNumber::Count(1));
        table.add_format("AD", TagType::Integer, TagNumber::R);
        if with_dp {
            table.add_format("DP", TagType::Integer, TagNumber::Count(1));
        }
        table
    }

    fn mask(field: &str, when: &str, depth: bool) -> FilterSpec {
        FilterSpec {
            genotype: GenotypeRules {
                mask_if: vec![GenotypeMaskRule {
                    field: field.into(),
                    when: when.into(),
                    on_missing: OnMissing::Keep,
                    depth,
                }],
                ..GenotypeRules::default()
            },
            ..FilterSpec::default()
        }
    }

    fn fmt_dp_expr() -> FilterSpec {
        FilterSpec {
            genotype: GenotypeRules {
                set_gts: Some(".".into()),
                ..GenotypeRules::default()
            },
            site_after_mask: TailRules {
                include: vec!["FMT/DP>10".into()],
                ..TailRules::default()
            },
            ..FilterSpec::default()
        }
    }

    fn typed_tags(dp: Option<TagType>, ad: TagType) -> TagTable {
        let mut table = TagTable::new();
        table.add_format("GT", TagType::Character, TagNumber::Count(1));
        table.add_format("AD", ad, TagNumber::R);
        if let Some(ty) = dp {
            table.add_format("DP", ty, TagNumber::Count(1));
        }
        table
    }

    #[test]
    fn float_ad_fallback_is_rejected_even_when_dp_is_declared() {
        let table = typed_tags(Some(TagType::Integer), TagType::Float);
        let err = compile(&mask("DP", "< 10", true), &table, "", true)
            .expect_err("a Float AD cannot be summed as depth");
        let FilterError::Spec(text) = &err else {
            panic!("expected a spec error, got {err:?}");
        };
        assert!(text.contains("AD must be Integer"), "{text}");
    }

    #[test]
    fn float_dp_depth_mask_is_rejected() {
        let table = typed_tags(Some(TagType::Float), TagType::Integer);
        let err = compile(&mask("DP", "< 10", true), &table, "", true)
            .expect_err("a Float DP is not a depth column");
        let FilterError::Spec(text) = &err else {
            panic!("expected a spec error, got {err:?}");
        };
        assert!(text.contains("DP is Float"), "{text}");
    }

    #[test]
    fn depth_mask_without_dp_column_uses_ad_fallback() {
        let plan = compile(&mask("DP", "< 10", true), &tags(false), "", true)
            .expect("invariant: an AD-only header supports the depth mask");
        assert_eq!(
            plan.int_fields,
            vec![IntField {
                name: b"DP".to_vec(),
                fallback: Some(b"AD".to_vec()),
            }]
        );
        assert_eq!(plan.mask[0].int_column, 0);
    }

    #[test]
    fn explicit_fmt_dp_expression_requires_the_tag_and_has_no_fallback() {
        let err = compile(&fmt_dp_expr(), &tags(false), "", true).expect_err("must reject");
        let FilterError::Spec(text) = &err else {
            panic!("expected a spec error, got {err:?}");
        };
        assert!(text.contains("DP"), "{text}");
        let plan = compile(&fmt_dp_expr(), &tags(true), "", true).expect("invariant: DP declared");
        let dp: Vec<_> = plan.int_fields.iter().filter(|f| f.name == b"DP").collect();
        assert_eq!(dp.len(), 1);
        assert_eq!(dp[0].fallback, None, "an expression is literal DP");
    }

    #[test]
    fn mask_and_expression_dp_are_two_columns() {
        let mut spec = mask("DP", "< 10", true);
        spec.genotype.set_gts = Some(".".into());
        spec.site_after_mask.include = vec!["FMT/DP>10".into()];
        let plan = compile(&spec, &tags(true), "", true).expect("invariant: compiles");
        let dp: Vec<_> = plan.int_fields.iter().filter(|f| f.name == b"DP").collect();
        assert_eq!(dp.len(), 2, "the mask and the expression must not share");
        assert!(dp.iter().any(|f| f.fallback.as_deref() == Some(b"AD")));
        assert!(dp.iter().any(|f| f.fallback.is_none()));
    }

    #[test]
    fn depth_is_only_valid_for_dp() {
        let err = compile(&mask("GQ", "< 10", true), &tags(true), "", true)
            .expect_err("depth = true on GQ is a spec error");
        let FilterError::Spec(text) = &err else {
            panic!("expected a spec error, got {err:?}");
        };
        assert!(text.contains("only valid for FORMAT/DP"), "{text}");
    }
}

#[cfg(test)]
mod region_contig_tests {
    use super::*;

    /// doc_v24 V24-04: the plan compiler resolves a colon-containing contig
    /// from the header vocabulary instead of splitting at the first `:`.
    #[test]
    fn indexed_region_uses_the_header_contig_vocabulary() {
        let mut tags = TagTable::new();
        tags.add_contig("HLA-A*01:01");
        let spec = FilterSpec {
            targets: crate::spec::RegionRule {
                indexed: vec!["HLA-A*01:01".into()],
                ..Default::default()
            },
            ..Default::default()
        };
        let plan = compile(&spec, &tags, "", true).expect("invariant: compiles");
        let mut seen = false;
        for rule in &plan.head {
            if let HeadRule::Regions { include, .. } = rule {
                seen |= include
                    .iter()
                    .any(|(c, s, e)| c == "HLA-A*01:01" && s.is_none() && e.is_none());
            }
        }
        assert!(
            seen,
            "the whole HLA contig must be one region: {:?}",
            plan.head
        );
    }

    /// The same text with an explicit interval; without the vocabulary it
    /// would parse as contig `HLA-A*01` position 1.
    #[test]
    fn indexed_region_with_a_range_on_a_colon_contig() {
        let mut tags = TagTable::new();
        tags.add_contig("HLA-A*01:01");
        let spec = FilterSpec {
            targets: crate::spec::RegionRule {
                indexed: vec!["HLA-A*01:01:10-20".into()],
                ..Default::default()
            },
            ..Default::default()
        };
        let plan = compile(&spec, &tags, "", true).expect("invariant: compiles");
        let mut seen = false;
        for rule in &plan.head {
            if let HeadRule::Regions { include, .. } = rule {
                seen |= include
                    .iter()
                    .any(|(c, s, e)| c == "HLA-A*01:01" && *s == Some(10) && *e == Some(20));
            }
        }
        assert!(seen, "expected HLA-A*01:01:10-20: {:?}", plan.head);
    }
}

#[cfg(test)]
mod routing_tests {
    use super::*;
    use crate::spec::{SiteRules, TailRules};

    fn tags() -> TagTable {
        let mut table = TagTable::new();
        table.add_format("GT", TagType::Character, TagNumber::Count(1));
        table.add_info("AC", TagType::Integer, TagNumber::A);
        table.add_info("DP", TagType::Integer, TagNumber::Count(1));
        table
    }

    /// doc_v25 V25-04: `INFO/AC` is the same built-in as the bare `AC`, so
    /// the two spellings must route the same way. `INFO/AC` used to fall
    /// through `references_genotypes` and reach the head stage, where INFO
    /// was assumed to describe the output cohort.
    #[test]
    fn info_spelling_routes_like_the_bare_builtin() {
        let table = tags();
        // doc_v25 review: an explicit `INFO/` read is a plain INFO value and
        // stays in the head stage (bcftools evaluates `-i/-e` before it
        // subsets the samples, `vcfview.c:322` vs `:341`).
        let spec = FilterSpec {
            site: SiteRules {
                include: vec!["INFO/AC>1".into()],
                ..SiteRules::default()
            },
            ..FilterSpec::default()
        };
        let plan = compile(&spec, &table, "", false).expect("INFO/AC is a head rule");
        assert_eq!(plan.head.len(), 1);
        assert!(plan.tail.is_empty());
        // The bare built-in still needs the decoded genotypes.
        let spec = FilterSpec {
            site: SiteRules {
                include: vec!["AC>1".into()],
                ..SiteRules::default()
            },
            ..FilterSpec::default()
        };
        let err = compile(&spec, &table, "", false)
            .expect_err("the bare genotype-derived built-in belongs after the mask");
        assert!(err.to_string().contains("site_after_mask"), "{err}");
        // A plain INFO tag is not genotype-derived: it stays in the head.
        let spec = FilterSpec {
            site: SiteRules {
                include: vec!["INFO/DP>10".into()],
                ..SiteRules::default()
            },
            ..FilterSpec::default()
        };
        let plan = compile(&spec, &table, "", false).expect("INFO/DP is a head rule");
        assert_eq!(plan.head.len(), 1);
        assert!(plan.tail.is_empty());
    }

    /// The after-mask stage still accepts the bare built-in (the routed
    /// half of the same rule); an explicit `INFO/` read now belongs in
    /// `[site]` and is rejected there with a pointer to it.
    #[test]
    fn info_spelling_compiles_in_the_after_mask_stage() {
        let spec = FilterSpec {
            site_after_mask: TailRules {
                include: vec!["AC>1".into()],
                ..TailRules::default()
            },
            ..FilterSpec::default()
        };
        let plan = compile(&spec, &tags(), "", false).expect("after-mask expression");
        assert!(plan.head.is_empty());
        assert_eq!(plan.tail.len(), 1);

        let spec = FilterSpec {
            site_after_mask: TailRules {
                include: vec!["INFO/AC>1".into()],
                ..TailRules::default()
            },
            ..FilterSpec::default()
        };
        let err =
            compile(&spec, &tags(), "", false).expect_err("an INFO-only read belongs in [site]");
        assert!(err.to_string().contains("[site]"), "{err}");
    }
}
