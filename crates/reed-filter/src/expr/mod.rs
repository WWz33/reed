//! bcftools expression subset (02-cli.md §4): parser, validator and
//! evaluator for `-i/-e`.
//!
//! The subset is deliberately explicit: anything outside it is reported as
//! "not supported yet" by [`check`], never executed with different
//! semantics. Supported: the site fields `QUAL`, `FILTER`, `ID`, `CHROM`,
//! `POS`, `TYPE`, `N_ALT`, `N_SAMPLES`, the built-ins `AC`, `AN`, `AF`,
//! `MAC`, `MAF`, `N_MISSING`, `F_MISSING` (INFO value when present, else
//! genotype-derived), `INFO/TAG` and bare tags, comparisons
//! (`== = != < <= > >=`), the match operators `~`/`!~` over literal
//! patterns with optional `^`/`$` anchors, logic (`&& || & |`), `!` and
//! parentheses.
//! `FMT/TAG` reads one sample's FORMAT value and is only meaningful through
//! `--set-gts .` (02-cli.md §3.3), where the caller evaluates the
//! expression per sample with the plan's scanned FORMAT table.
//! Not yet supported: subscripts, general regular expressions (any pattern
//! metacharacter beyond the anchors is reported, never reinterpreted),
//! arithmetic.

pub mod lexer;
pub mod parser;

use reed_io::vcf::info::info_fragments;

pub use parser::{CmpOp, Expr, Field, LogicOp, Subscript};

/// Genotype-derived counts for a site (available after the mask stage).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GtCounts {
    /// ALT allele count.
    pub ac: u64,
    /// Called allele count.
    pub an: u64,
    /// Fully missing genotypes.
    pub n_missing: u64,
    /// Samples.
    pub n_samples: u64,
}

/// Everything an expression can read about one site.
#[derive(Clone, Copy, Debug)]
pub struct SiteValues<'a> {
    /// QUAL (`None` when the record has `.`).
    pub qual: Option<f64>,
    /// FILTER column.
    pub filter: &'a [u8],
    /// ID column.
    pub id: &'a [u8],
    /// CHROM column.
    pub chrom: &'a [u8],
    /// 1-based POS.
    pub pos: u64,
    /// Variant class name (`snp`, `mnp`, `indel`, `other`, `ref`).
    pub variant_type: &'a str,
    /// ALT allele count.
    pub n_alt: u32,
    /// Sample count.
    pub n_samples: u64,
    /// Raw INFO text.
    pub info: &'a [u8],
    /// Whether the record's `INFO/AC,AN` still describe the cohort the
    /// rules decide on. False once the sample set was subset or genotypes
    /// were masked: the INFO values then describe the *input* cohort while
    /// the genotype-derived counts describe the kept one, and bcftools
    /// switches to the genotypes in exactly that situation
    /// (`vcfview.c:341-350` recomputes AC/AN from the subset GT before
    /// applying `-q/-Q`; `filter.c:1471` gates its INFO preference on
    /// `!sample_mask_set`). doc_v12 D12-6.
    pub info_trusted: bool,
    /// This record's INFO fragment index (A6-18): set by callers that built
    /// one for the record — `None` falls back to a per-tag scan.
    pub info_index: Option<&'a reed_io::vcf::info::InfoIndex>,
    /// Genotype-derived counts (`None` before the mask stage).
    pub gt: Option<GtCounts>,
    /// The current sample's `FORMAT` values, present when the caller
    /// evaluates per sample (`--set-gts .`, 02-cli.md §3.3): the scanned
    /// `(tag, parsed number)` table borrowed in place — single-valued,
    /// numeric-typed tags only (multi-valued tags need subscripts and are
    /// rejected at compile time).
    pub fmt: Option<&'a [(Vec<u8>, Option<f64>)]>,
    /// The header's declared tags. INFO lookups need the declaration to be
    /// correct: a Flag reads as 1/0 (presence, bcftools
    /// `filters_set_info_flag`, filter.c:964-971) and a multi-valued tag
    /// compares any-element (bcftools `CMP_VECTORS`, filter.c:2738-2790) —
    /// both impossible from the raw text alone (doc_v12 D12-7).
    pub tags: Option<&'a reed_core::tag_table::TagTable>,
}

impl<'a> SiteValues<'a> {
    /// Builds the context from a record line; `gt` carries the
    /// genotype-derived counts when they are known (after the mask stage).
    #[must_use]
    pub fn from_line(
        line: &'a reed_io::vcf::line::LineRef<'a>,
        pos1: u64,
        gt: Option<GtCounts>,
    ) -> Self {
        use reed_core::site::{VariantClass, classify_variant};
        let reference = line.reference();
        let alts = line.alts();
        let primary = match memchr::memchr(b',', alts) {
            Some(i) => &alts[..i],
            None => alts,
        };
        let variant_type = match classify_variant(reference, primary) {
            VariantClass::Snv => "snp",
            VariantClass::Mnv => "mnp",
            VariantClass::Indel => "indel",
            VariantClass::Symbolic => "symbolic",
            VariantClass::Ref => "ref",
            VariantClass::Other | _ => "other",
        };
        let qual = reed_io::vcf::fixed::parse_qual(line.qual_text(), (line.chrom(), pos1))
            .ok()
            .flatten()
            .map(f64::from);
        Self {
            qual,
            filter: line.filter_text(),
            id: line.id(),
            chrom: line.chrom(),
            pos: pos1,
            variant_type,
            n_alt: reed_io::vcf::site::count_alts(alts),
            n_samples: line.n_sample_cols() as u64,
            info: line.info_text(),
            info_trusted: true,
            info_index: None,
            gt,
            fmt: None,
            tags: None,
        }
    }
}

impl SiteValues<'static> {
    /// Context for the genotype-derived fields only (the sketch path's mid
    /// stage, which has counts but no record bytes). Reading a line field
    /// from this context yields the empty/dummy values below; the plan
    /// forces the single-pass strategy whenever a tail expression exists, so
    /// only genotype-derived fields are ever read here.
    #[must_use]
    pub fn from_counts(gt: GtCounts) -> Self {
        Self {
            qual: None,
            filter: b".",
            id: b".",
            chrom: b"",
            pos: 0,
            variant_type: "other",
            n_alt: 0,
            n_samples: gt.n_samples,
            info: b".",
            info_trusted: false,
            info_index: None,
            gt: Some(gt),
            fmt: None,
            tags: None,
        }
    }
}

/// A resolved value.
#[derive(Clone, Debug, PartialEq)]
enum Val<'a> {
    Num(f64),
    Str(&'a [u8]),
    Missing,
}

impl Val<'_> {
    fn as_num(&self) -> Option<f64> {
        match self {
            Self::Num(n) => Some(*n),
            Self::Str(bytes) => std::str::from_utf8(bytes).ok()?.trim().parse().ok(),
            Self::Missing => None,
        }
    }

    /// String rendering for mixed-type comparisons (no allocation for
    /// string values).
    fn as_text(&self) -> String {
        match self {
            Self::Num(n) => format!("{n}"),
            Self::Str(bytes) => String::from_utf8_lossy(bytes).into_owned(),
            Self::Missing => String::new(),
        }
    }
}

/// Parses `text` without a header (for routing decisions); use
/// [`parse_and_check`] for the header-validated form.
///
/// # Errors
///
/// A human-readable syntax error.
pub fn parse(text: &str) -> Result<Expr, String> {
    parser::parse(text)
}

/// Parses and validates `text` against the header's declared INFO tags.
///
/// # Errors
///
/// A human-readable message (syntax error, unknown tag, or a construct the
/// subset does not implement yet).
pub fn parse_and_check(text: &str, tags: &reed_core::tag_table::TagTable) -> Result<Expr, String> {
    let expr = parser::parse(text)?;
    check(&expr, tags, &[])?;
    Ok(expr)
}

/// Parses and validates `text` against the header's INFO and FORMAT tags.
///
/// # Errors
///
/// A human-readable message (syntax error, unknown tag, or a construct the
/// subset does not implement yet).
pub fn parse_and_check_fmt(
    text: &str,
    tags: &reed_core::tag_table::TagTable,
    format_tags: &[String],
) -> Result<Expr, String> {
    let expr = parser::parse(text)?;
    check(&expr, tags, format_tags)?;
    Ok(expr)
}

/// Validates an expression against the declared INFO tags.
///
/// # Errors
///
/// Unknown tags and unimplemented constructs.
pub fn check(
    expr: &Expr,
    tags: &reed_core::tag_table::TagTable,
    format_tags: &[String],
) -> Result<(), String> {
    check_inner(expr, tags, format_tags, false)
}

/// `in_cmp` is true once the walk is inside a comparison: a bare tag is only
/// legal there is no comparison to give it a value (bcftools evaluates a
/// lone token to false, `filter.c:4352-4380`; reed reports it at compile
/// time instead of silently dropping every record — 01-review §2.5-C1).
fn check_inner(
    expr: &Expr,
    tags: &reed_core::tag_table::TagTable,
    format_tags: &[String],
    in_cmp: bool,
) -> Result<(), String> {
    match expr {
        Expr::Fmt(tag) => {
            if !format_tags.iter().any(|t| t == tag) {
                return Err(format!("FORMAT tag {tag:?} is not declared in the header"));
            }
            Ok(())
        }
        Expr::Subscript(inner, _) => {
            // A subscript selects an element of a multi-valued value; the
            // planned columns carry single values only.
            check_inner(inner, tags, format_tags, in_cmp)?;
            Err(
                "subscripts ([0], [*]) are not supported yet (multi-valued FORMAT tags)"
                    .to_string(),
            )
        }
        Expr::Cmp(op @ (CmpOp::Match | CmpOp::NotMatch), lhs, rhs) => {
            let Expr::Str(pattern) = &**rhs else {
                return Err(format!(
                    "the right-hand side of {} must be a string literal",
                    op.as_str()
                ));
            };
            check_pattern(pattern)?;
            check_inner(lhs, tags, format_tags, true)
        }
        Expr::Info(tag) | Expr::Bare(tag) => {
            if Field::from_name(tag).is_none() && tags.info(tag).is_none() {
                return Err(format!("INFO tag {tag:?} is not declared in the header"));
            }
            // A tag outside a comparison is always false in bcftools:
            // `filter_test` never sets `pass_site` without a comparison node
            // (filter.c:4352-4380). Report it instead of filtering every
            // record away (01-review §2.5-C1). Literals stay legal — they
            // are the right-hand sides reached through `Expr::Cmp`.
            if in_cmp {
                return Ok(());
            }
            Err(format!(
                "tag {tag:?} without a comparison is always false in bcftools; \
                 write e.g. `{tag}=1` (flag presence) or `{tag}>0`"
            ))
        }
        // A comparison's operands are values: a tag there is legal.
        Expr::Cmp(_, a, b) => {
            check_inner(a, tags, format_tags, true)?;
            check_inner(b, tags, format_tags, true)
        }
        // Logic operands are whole predicates: the flag passes through.
        Expr::Logic(_, a, b) => {
            check_inner(a, tags, format_tags, in_cmp)?;
            check_inner(b, tags, format_tags, in_cmp)
        }
        Expr::Not(e) => check_inner(e, tags, format_tags, in_cmp),
        _ => Ok(()),
    }
}

/// Evaluates an expression; missing values make comparisons false
/// (bcftools semantics).
#[must_use]
pub fn eval(expr: &Expr, site: &SiteValues<'_>) -> bool {
    match expr {
        Expr::Cmp(op @ (CmpOp::Match | CmpOp::NotMatch), a, b) => {
            let (Some(left), Some(Expr::Str(pattern))) = (value(a, site), Some(&**b)) else {
                return false;
            };
            if matches!(left, Val::Missing) {
                return false;
            }
            let hit = matches_pattern(pattern, &left.as_text());
            if *op == CmpOp::Match { hit } else { !hit }
        }
        Expr::Cmp(op, a, b) => {
            // A multi-valued INFO tag compares any-element: bcftools walks
            // every pair and passes the site when one satisfies the
            // comparison (`CMP_VECTORS`, filter.c:2738-2790; doc_v12 D12-7).
            let left = operand(a, site);
            let right = operand(b, site);
            let (Some(left), Some(right)) = (left, right) else {
                return false;
            };
            left.any_pair(&right, *op)
        }
        Expr::Logic(op, a, b) => match *op {
            LogicOp::And => eval(a, site) && eval(b, site),
            LogicOp::Or => eval(a, site) || eval(b, site),
        },
        Expr::Not(e) => !eval(e, site),
        // A bare expression (no comparison) is true when it is a non-zero
        // number or a non-missing value, matching bcftools' truthiness.
        other => match value(other, site) {
            Some(Val::Num(n)) => n != 0.0 && !n.is_nan(),
            Some(Val::Str(s)) => !s.is_empty(),
            _ => false,
        },
    }
}

/// Resolves one node to a value.
fn value<'a>(expr: &'a Expr, site: &'a SiteValues<'a>) -> Option<Val<'a>> {
    match expr {
        Expr::Number(n) => Some(Val::Num(*n)),
        Expr::Str(s) => Some(Val::Str(s.as_bytes())),
        Expr::Missing => Some(Val::Missing),
        Expr::Field(field) => field_value(*field, site),
        Expr::Fmt(tag) => site.fmt.and_then(|values| {
            values
                .iter()
                .find(|(t, _)| t.as_slice() == tag.as_bytes())
                .and_then(|(_, num)| num.map_or(Some(Val::Missing), |n| Some(Val::Num(n))))
        }),
        Expr::Info(tag) | Expr::Bare(tag) => {
            if let Some(field) = Field::from_name(tag) {
                return field_value(field, site);
            }
            info_operand(site, tag.as_bytes()).scalar()
        }
        Expr::Subscript(_, _) | Expr::Cmp(..) | Expr::Logic(..) | Expr::Not(_) => None,
    }
}

/// Resolves a built-in field: INFO first for the allele-count family, then
/// genotype-derived counts.
fn field_value<'a>(field: Field, site: &'a SiteValues<'a>) -> Option<Val<'a>> {
    match field {
        Field::Qual => site.qual.map_or(Some(Val::Missing), |q| Some(Val::Num(q))),
        Field::Filter => Some(Val::Str(site.filter)),
        Field::Id => Some(Val::Str(site.id)),
        Field::Chrom => Some(Val::Str(site.chrom)),
        Field::Pos => Some(Val::Num(site.pos as f64)),
        Field::Type => Some(Val::Str(site.variant_type.as_bytes())),
        Field::NAlt => Some(Val::Num(f64::from(site.n_alt))),
        Field::NSamples => Some(Val::Num(site.n_samples as f64)),
        Field::Ac => info_or_gt(site, b"AC", |g| g.ac as f64),
        Field::An => info_or_gt(site, b"AN", |g| g.an as f64),
        Field::Af => info_or_gt(site, b"AF", |g| {
            if g.an > 0 {
                g.ac as f64 / g.an as f64
            } else {
                f64::NAN
            }
        }),
        Field::Mac => info_or_gt(site, b"MAC", |g| {
            let major = g.ac.max(g.an.saturating_sub(g.ac));
            (g.an.saturating_sub(major)) as f64
        }),
        Field::Maf => info_or_gt(site, b"MAF", |g| {
            if g.an > 0 {
                let af = g.ac as f64 / g.an as f64;
                af.min(1.0 - af)
            } else {
                f64::NAN
            }
        }),
        Field::NMissing => site.gt.map(|g| Val::Num(g.n_missing as f64)),
        Field::FMissing => site.gt.map(|g| {
            if g.n_samples > 0 {
                Val::Num(g.n_missing as f64 / g.n_samples as f64)
            } else {
                Val::Missing
            }
        }),
    }
}

/// INFO value when the tag is present, else the genotype-derived fallback
/// (bcftools' rule for AC/AN/AF/MAC/MAF).
fn info_or_gt<'a>(
    site: &'a SiteValues<'a>,
    tag: &[u8],
    from_gt: impl Fn(&GtCounts) -> f64,
) -> Option<Val<'a>> {
    if let Some(Val::Num(_) | Val::Str(_)) = info_value(site.info, tag, site.info_index) {
        return info_value(site.info, tag, site.info_index);
    }
    site.gt
        .map(|g| Val::Num(from_gt(&g)))
        .or(Some(Val::Missing))
}

/// Looks one tag up in the raw INFO text (first comma-separated element);
/// through the record's fragment index when the caller built one (A6-18).
fn info_value<'a>(
    info: &'a [u8],
    tag: &[u8],
    index: Option<&reed_io::vcf::info::InfoIndex>,
) -> Option<Val<'a>> {
    let fragment = match index {
        Some(index) => index.get(info, tag),
        None => info_fragments(info)
            .find(|(name, _)| *name == tag)
            .map(|(_, value)| value),
    };
    let Some(value) = fragment? else {
        return Some(Val::Missing);
    };
    if value == b"." || value.is_empty() {
        return Some(Val::Missing);
    }
    let first = value.split(|&b| b == b',').next().unwrap_or(value);
    Some(
        match std::str::from_utf8(first)
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
        {
            Some(n) => Val::Num(n),
            None => Val::Str(first),
        },
    )
}

/// One side of a comparison: a scalar value or the elements of a
/// multi-valued INFO tag.
#[derive(Clone, Debug)]
enum Operand<'a> {
    /// Zero or one value (the bcftools "scalar" view).
    Scalar(Option<Val<'a>>),
    /// Every element of a multi-valued tag (Number=A/R/G or a count > 1).
    Vector(Vec<Val<'a>>),
}

impl<'a> Operand<'a> {
    /// The single value of a scalar operand (a vector collapses to its
    /// first element, which is what a scalar context means).
    fn scalar(self) -> Option<Val<'a>> {
        match self {
            Self::Scalar(v) => v,
            Self::Vector(mut v) => {
                if v.is_empty() {
                    None
                } else {
                    Some(v.swap_remove(0))
                }
            }
        }
    }

    /// bcftools' comparison rule over two operands: any pair satisfies the
    /// comparison, and a missing element on either side only passes for the
    /// operators that treat missing as a value (`=`/`!=`, the
    /// `missing_logic` table of `CMP_VECTORS`). Reed keeps the simpler
    /// "missing makes comparisons false" rule for missing *elements* and
    /// only implements the any-pair walk here.
    fn any_pair(&self, other: &Self, op: CmpOp) -> bool {
        let scalars_a: Vec<Option<Val<'a>>> = match self {
            Self::Scalar(v) => vec![v.clone()],
            Self::Vector(v) => v.iter().cloned().map(Some).collect(),
        };
        let scalars_b: Vec<Option<Val<'a>>> = match other {
            Self::Scalar(v) => vec![v.clone()],
            Self::Vector(v) => v.iter().cloned().map(Some).collect(),
        };
        for a in &scalars_a {
            for b in &scalars_b {
                let (Some(a), Some(b)) = (a, b) else { continue };
                if matches!(a, Val::Missing) || matches!(b, Val::Missing) {
                    continue;
                }
                if compare(op, a, b) {
                    return true;
                }
            }
        }
        false
    }
}

/// Resolves one comparison side: a multi-valued INFO tag becomes a vector,
/// everything else a scalar.
fn operand<'a>(expr: &'a Expr, site: &'a SiteValues<'a>) -> Option<Operand<'a>> {
    if let Expr::Info(tag) | Expr::Bare(tag) = expr
        && Field::from_name(tag).is_none()
    {
        return Some(info_operand(site, tag.as_bytes()));
    }
    value(expr, site).map(|v| Operand::Scalar(Some(v)))
}

/// One INFO tag as a comparison operand, honouring the header's declaration:
///
/// - an undeclared tag has no shape to honour and reads as a scalar;
/// - a Flag reads as 1/0 by presence (bcftools `filters_set_info_flag`,
///   filter.c:964-971; "1 (or 0) to test the presence (or absence) of a
///   flag", doc/bcftools.txt:4067) — a bare flag used to be `Val::Missing`
///   and matched nothing (doc_v12 D12-7);
/// - a multi-valued tag (`Number=A/R/G/.` or `Count(n)` with n > 1) becomes
///   the element vector, compared any-element;
/// - a `Number=1` tag with several values is malformed input: bcftools'
///   scalar reader happens to return 0 there (htslib leaves `v1` at its zero
///   initialiser when `len != 1`, vcf.c:4286), but that is an implementation
///   detail, so reed reads it as missing (01-review §2.5-B5/C2).
fn info_operand<'a>(site: &'a SiteValues<'a>, tag: &[u8]) -> Operand<'a> {
    use reed_core::tag_table::TagNumber;
    let declared = site.tags.and_then(|t| {
        std::str::from_utf8(tag)
            .ok()
            .and_then(|name| t.info(name))
            .map(|def| (def.ty, def.number))
    });
    let fragment = match site.info_index {
        Some(index) => index.get(site.info, tag),
        None => info_fragments(site.info)
            .find(|(name, _)| *name == tag)
            .map(|(_, value)| value),
    };
    let is_flag = matches!(declared, Some((reed_core::tag_table::TagType::Flag, _)));
    let Some(fragment) = fragment else {
        // The tag is absent: a declared Flag reads 0 ("0 to test the
        // absence of a flag", doc/bcftools.txt:4067); anything else has no
        // value to compare.
        return Operand::Scalar(Some(if is_flag { Val::Num(0.0) } else { Val::Missing }));
    };
    let Some(value) = fragment else {
        // A bare key: a declared Flag is present (1); anything else has no
        // value to compare.
        return Operand::Scalar(Some(if is_flag { Val::Num(1.0) } else { Val::Missing }));
    };
    if value == b"." || value.is_empty() {
        return Operand::Scalar(Some(Val::Missing));
    }
    let multi = matches!(
        declared.map(|(_, n)| n),
        Some(TagNumber::A | TagNumber::R | TagNumber::G | TagNumber::Unknown)
            | Some(TagNumber::Count(2..))
    );
    let elements: Vec<Val<'a>> = value.split(|&b| b == b',').map(parse_element).collect();
    if multi {
        return Operand::Vector(elements);
    }
    if elements.len() > 1 {
        // A Number=1 tag carrying several values is malformed (01 §2.5-B5).
        return Operand::Scalar(Some(Val::Missing));
    }
    Operand::Scalar(elements.into_iter().next())
}

/// One INFO element as a value: `.` is missing, a number parses, everything
/// else stays a string for the byte-wise comparisons.
fn parse_element(part: &[u8]) -> Val<'_> {
    if part == b"." || part.is_empty() {
        return Val::Missing;
    }
    match std::str::from_utf8(part)
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
    {
        Some(n) => Val::Num(n),
        None => Val::Str(part),
    }
}

/// Characters that turn a pattern into a general regular expression, which
/// this subset does not implement (02-cli.md §4: report, never reinterpret).
const REGEX_META: &[char] = &['.', '*', '+', '?', '[', ']', '(', ')', '{', '}', '|', '\\'];

/// Validates a `~`/`!~` pattern: literal text with optional `^`/`$` anchors.
fn check_pattern(pattern: &str) -> Result<(), String> {
    let body = pattern.strip_prefix('^').unwrap_or(pattern);
    let body = body.strip_suffix('$').unwrap_or(body);
    if let Some(bad) = body.chars().find(|c| REGEX_META.contains(c)) {
        return Err(format!(
            "regular-expression metacharacter {bad:?} is not supported yet (only literal patterns with ^/$ anchors)"
        ));
    }
    Ok(())
}

/// Matches a validated pattern against a value: `^` anchors the start, `$`
/// the end, otherwise the literal must occur anywhere (bcftools' `~`).
#[must_use]
pub fn matches_pattern(pattern: &str, text: &str) -> bool {
    let anchored_start = pattern.starts_with('^');
    let anchored_end = pattern.ends_with('$') && pattern.len() > usize::from(anchored_start);
    let body = pattern.strip_prefix('^').unwrap_or(pattern);
    let body = if anchored_end {
        &body[..body.len() - 1]
    } else {
        body
    };
    match (anchored_start, anchored_end) {
        (true, true) => text == body,
        (true, false) => text.starts_with(body),
        (false, true) => text.ends_with(body),
        (false, false) => text.contains(body),
    }
}

/// Comparison with bcftools' rules: numeric when both sides parse as
/// numbers, otherwise byte-wise.
///
/// The float comparisons are exact (`==`, `<`, …), deliberately: bcftools
/// evaluates the same operators on the same f64 values, so an epsilon here
/// would change which records survive.
#[allow(clippy::float_cmp)]
fn compare(op: CmpOp, left: &Val<'_>, right: &Val<'_>) -> bool {
    if let (Some(a), Some(b)) = (left.as_num(), right.as_num()) {
        return match op {
            CmpOp::Eq => a == b,
            CmpOp::Ne => a != b,
            CmpOp::Lt => a < b,
            CmpOp::Le => a <= b,
            CmpOp::Gt => a > b,
            CmpOp::Ge => a >= b,
            CmpOp::Match | CmpOp::NotMatch => false,
        };
    }
    let (a, b) = (left.as_text(), right.as_text());
    let (a, b) = (a.as_bytes(), b.as_bytes());
    match op {
        CmpOp::Eq => a == b,
        CmpOp::Ne => a != b,
        CmpOp::Lt => a < b,
        CmpOp::Le => a <= b,
        CmpOp::Gt => a > b,
        CmpOp::Ge => a >= b,
        CmpOp::Match | CmpOp::NotMatch => false,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(info: &[u8], qual: Option<f64>, gt: Option<GtCounts>) -> SiteValues<'_> {
        SiteValues {
            qual,
            filter: b"PASS",
            id: b".",
            chrom: b"chr1",
            pos: 100,
            variant_type: "snp",
            n_alt: 1,
            n_samples: 4,
            info,
            info_trusted: true,
            info_index: None,
            gt,
            fmt: None,
            tags: None,
        }
    }

    /// A TagTable declaring each name as `Type`, `Number` (AC is Number=A).
    fn table(names: &[&str]) -> reed_core::tag_table::TagTable {
        use reed_core::tag_table::{TagNumber, TagTable, TagType};
        let mut t = TagTable::new();
        for name in names {
            if *name == "AC" {
                t.add_info(name, TagType::Integer, TagNumber::A);
            } else {
                t.add_info(name, TagType::Integer, TagNumber::Count(1));
            }
        }
        t
    }

    fn ev(text: &str, s: &SiteValues<'_>, tags: &[&str]) -> bool {
        let expr =
            parse_and_check(text, &table(tags)).unwrap_or_else(|e| panic!("{text}: {e}"));
        eval(&expr, s)
    }

    #[test]
    fn gold_expressions_match_their_typed_equivalents() {
        let tags = ["QD", "FS", "MQ", "SOR"];
        // QUAL >= 30
        assert!(ev("QUAL>=30", &site(b".", Some(50.0), None), &tags));
        assert!(!ev("QUAL>=30", &site(b".", Some(10.0), None), &tags));
        assert!(
            !ev("QUAL>=30", &site(b".", None, None), &tags),
            "missing QUAL fails"
        );
        // QD<4 || FS>50 || MQ<40 || SOR>3 (the gold exclusion)
        let bad = b"QD=2.0;FS=10.0;MQ=60.0;SOR=1.0";
        assert!(ev(
            "QD<4 || FS>50 || MQ<40 || SOR>3",
            &site(bad, Some(50.0), None),
            &tags
        ));
        let good = b"QD=10.0;FS=10.0;MQ=60.0;SOR=1.0";
        assert!(!ev(
            "QD<4 || FS>50 || MQ<40 || SOR>3",
            &site(good, Some(50.0), None),
            &tags
        ));
        // Missing INFO compares false (bcftools).
        assert!(!ev("QD<4", &site(b".", Some(50.0), None), &tags));
    }

    #[test]
    fn strings_logic_and_genotype_fallbacks() {
        let tags: [&str; 0] = [];
        assert!(ev("FILTER=\"PASS\"", &site(b".", Some(1.0), None), &tags));
        assert!(!ev(
            "FILTER=\"LowQual\"",
            &site(b".", Some(1.0), None),
            &tags
        ));
        assert!(ev(
            "CHROM=\"chr1\" && POS==100",
            &site(b".", Some(1.0), None),
            &tags
        ));
        assert!(ev("TYPE=\"snp\"", &site(b".", Some(1.0), None), &tags));
        assert!(ev("!(QUAL<10)", &site(b".", Some(50.0), None), &tags));

        // AC from INFO when present, else from genotypes.
        assert!(ev("AC==5", &site(b"AC=5", Some(1.0), None), &["AC"]));
        let gt = GtCounts {
            ac: 3,
            an: 8,
            n_missing: 2,
            n_samples: 4,
        };
        assert!(ev("AC==3", &site(b".", Some(1.0), Some(gt)), &[]));
        assert!(ev("AF==0.375", &site(b".", Some(1.0), Some(gt)), &[]));
        assert!(ev("F_MISSING==0.5", &site(b".", Some(1.0), Some(gt)), &[]));
        assert!(
            !ev("F_MISSING==0.5", &site(b".", Some(1.0), None), &[]),
            "no counts yet"
        );
        assert!(ev("MAF==0.375", &site(b".", Some(1.0), Some(gt)), &[]));
    }

    #[test]
    fn unsupported_constructs_are_reported() {
        let empty = table(&[]);
        assert!(parse_and_check_fmt("FMT/DP>10", &empty, &["DP".into()]).is_ok());
        assert!(
            parse_and_check_fmt("FMT/DP>10", &empty, &[])
                .unwrap_err()
                .contains("not declared")
        );
        assert!(
            parse_and_check("INFO/AC[0]>1", &table(&["AC"]))
                .unwrap_err()
                .contains("subscripts")
        );
        assert!(
            parse_and_check("ID~\"rs[0-9]\"", &empty)
                .unwrap_err()
                .contains("metacharacter")
        );
        assert!(
            parse_and_check("ID~QSOMETHING", &empty)
                .unwrap_err()
                .contains("string literal")
        );
        assert!(
            parse_and_check("NOSUCH>1", &empty)
                .unwrap_err()
                .contains("not declared")
        );
        assert!(parse_and_check("QUAL>", &empty).is_err());
        // A tag with no comparison is always false in bcftools; reed
        // reports it (doc_v12 D12-7 / 01-review §2.5-C1).
        assert!(
            parse_and_check("DP", &table(&["DP"]))
                .unwrap_err()
                .contains("always false")
        );
    }

    #[test]
    fn literal_patterns_with_anchors() {
        let tags: [&str; 0] = [];
        let s = site(b".", Some(1.0), None);
        // ID is "." here; use FILTER/CHROM for the string values.
        assert!(ev("CHROM~\"chr\"", &s, &tags));
        assert!(ev("CHROM~\"^chr1$\"", &s, &tags));
        assert!(!ev("CHROM~\"^hr1$\"", &s, &tags));
        assert!(ev("CHROM!~\"^chr2\"", &s, &tags));
        assert!(!ev("CHROM~\"nope\"", &s, &tags));
        assert!(ev("FILTER~\"ASS\"", &s, &tags));
        assert!(!ev("FILTER~\"^ASS\"", &s, &tags));
    }

    #[test]
    fn genotype_dependence_placement() {
        assert!(
            !parser::parse("QUAL>=30 && INFO/DP>10")
                .expect("inv")
                .references_genotypes()
        );
        assert!(
            parser::parse("F_MISSING<0.1")
                .expect("inv")
                .references_genotypes()
        );
        assert!(
            parser::parse("N_MISSING<2")
                .expect("inv")
                .references_genotypes()
        );
        assert!(
            !parser::parse("FILTER=\"PASS\"")
                .expect("inv")
                .references_genotypes()
        );
    }
}
