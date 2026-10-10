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
//! expression per sample with a [`FmtValue`] table.
//! Not yet supported: subscripts, general regular expressions (any pattern
//! metacharacter beyond the anchors is reported, never reinterpreted),
//! arithmetic.

pub mod lexer;
pub mod parser;

use reed_io::vcf::info::info_fragments;

pub use parser::{Expr, Field, Subscript};

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

/// One sample's `FORMAT` value for a tag the plan scanned (single-valued,
/// numeric-typed tags only; multi-valued tags need subscripts and are
/// rejected at compile time).
#[derive(Clone, Copy, Debug)]
pub struct FmtValue<'a> {
    /// The FORMAT key.
    pub tag: &'a [u8],
    /// The parsed number (`None` when missing).
    pub num: Option<f64>,
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
    /// Genotype-derived counts (`None` before the mask stage).
    pub gt: Option<GtCounts>,
    /// The current sample's `FORMAT` values, present when the caller
    /// evaluates per sample (`--set-gts .`, 02-cli.md §3.3).
    pub fmt: Option<&'a [FmtValue<'a>]>,
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
        let qual = reed_io::vcf::fixed::parse_qual(
            line.qual_text(),
            &reed_io::Locus {
                contig_name: String::from_utf8_lossy(line.chrom()).into_owned(),
                pos1,
            },
        )
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
            n_alt: u32::try_from(memchr::memchr_iter(b',', alts).count() + 1).unwrap_or(u32::MAX),
            n_samples: line.n_sample_cols() as u64,
            info: line.info_text(),
            gt,
            fmt: None,
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
            gt: Some(gt),
            fmt: None,
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
pub fn parse_and_check(text: &str, info_tags: &[String]) -> Result<Expr, String> {
    let expr = parser::parse(text)?;
    check(&expr, info_tags, &[])?;
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
    info_tags: &[String],
    format_tags: &[String],
) -> Result<Expr, String> {
    let expr = parser::parse(text)?;
    check(&expr, info_tags, format_tags)?;
    Ok(expr)
}

/// Validates an expression against the declared INFO tags.
///
/// # Errors
///
/// Unknown tags and unimplemented constructs.
pub fn check(expr: &Expr, info_tags: &[String], format_tags: &[String]) -> Result<(), String> {
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
            check(inner, info_tags, format_tags)?;
            Err(
                "subscripts ([0], [*]) are not supported yet (multi-valued FORMAT tags)"
                    .to_string(),
            )
        }
        Expr::Cmp(op @ ("~" | "!~"), lhs, rhs) => {
            let Expr::Str(pattern) = &**rhs else {
                return Err(format!(
                    "the right-hand side of {op} must be a string literal"
                ));
            };
            check_pattern(pattern)?;
            check(lhs, info_tags, format_tags)
        }
        Expr::Info(tag) | Expr::Bare(tag) => {
            if Field::from_name(tag).is_none() && !info_tags.iter().any(|t| t == tag) {
                return Err(format!("INFO tag {tag:?} is not declared in the header"));
            }
            Ok(())
        }
        Expr::Cmp(_, a, b) | Expr::Logic(_, a, b) => {
            check(a, info_tags, format_tags)?;
            check(b, info_tags, format_tags)
        }
        Expr::Not(e) => check(e, info_tags, format_tags),
        _ => Ok(()),
    }
}

/// Evaluates an expression; missing values make comparisons false
/// (bcftools semantics).
#[must_use]
pub fn eval(expr: &Expr, site: &SiteValues<'_>) -> bool {
    match expr {
        Expr::Cmp(op @ ("~" | "!~"), a, b) => {
            let (Some(left), Some(Expr::Str(pattern))) = (value(a, site), Some(&**b)) else {
                return false;
            };
            if matches!(left, Val::Missing) {
                return false;
            }
            let hit = matches_pattern(pattern, &left.as_text());
            if *op == "~" { hit } else { !hit }
        }
        Expr::Cmp(op, a, b) => {
            let (Some(left), Some(right)) = (value(a, site), value(b, site)) else {
                return false;
            };
            if matches!(left, Val::Missing) || matches!(right, Val::Missing) {
                return false;
            }
            compare(op, &left, &right)
        }
        Expr::Logic(op, a, b) => match *op {
            "&&" | "&" => eval(a, site) && eval(b, site),
            "||" | "|" => eval(a, site) || eval(b, site),
            _ => false,
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
                .find(|v| v.tag == tag.as_bytes())
                .map(|v| v.num.map_or(Val::Missing, Val::Num))
        }),
        Expr::Info(tag) | Expr::Bare(tag) => {
            if let Some(field) = Field::from_name(tag) {
                return field_value(field, site);
            }
            info_value(site.info, tag.as_bytes())
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
    if let Some(Val::Num(_) | Val::Str(_)) = info_value(site.info, tag) {
        return info_value(site.info, tag);
    }
    site.gt
        .map(|g| Val::Num(from_gt(&g)))
        .or(Some(Val::Missing))
}

/// Looks one tag up in the raw INFO text (first comma-separated element).
fn info_value<'a>(info: &'a [u8], tag: &[u8]) -> Option<Val<'a>> {
    for (name, value) in info_fragments(info) {
        if name == tag {
            let Some(value) = value else {
                return Some(Val::Missing);
            };
            if value == b"." || value.is_empty() {
                return Some(Val::Missing);
            }
            let first = value.split(|&b| b == b',').next().unwrap_or(value);
            return Some(
                match std::str::from_utf8(first)
                    .ok()
                    .and_then(|s| s.parse::<f64>().ok())
                {
                    Some(n) => Val::Num(n),
                    None => Val::Str(first),
                },
            );
        }
    }
    None
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
fn compare(op: &str, left: &Val<'_>, right: &Val<'_>) -> bool {
    if let (Some(a), Some(b)) = (left.as_num(), right.as_num()) {
        return match op {
            "==" | "=" => a == b,
            "!=" => a != b,
            "<" => a < b,
            "<=" => a <= b,
            ">" => a > b,
            ">=" => a >= b,
            _ => false,
        };
    }
    let (a, b) = (left.as_text(), right.as_text());
    let (a, b) = (a.as_bytes(), b.as_bytes());
    match op {
        "==" | "=" => a == b,
        "!=" => a != b,
        "<" => a < b,
        "<=" => a <= b,
        ">" => a > b,
        ">=" => a >= b,
        _ => false,
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
            gt,
            fmt: None,
        }
    }

    fn ev(text: &str, s: &SiteValues<'_>, tags: &[&str]) -> bool {
        let tags: Vec<String> = tags.iter().map(|t| (*t).to_string()).collect();
        let expr = parse_and_check(text, &tags).unwrap_or_else(|e| panic!("{text}: {e}"));
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
        assert!(parse_and_check_fmt("FMT/DP>10", &[], &["DP".into()]).is_ok());
        assert!(
            parse_and_check_fmt("FMT/DP>10", &[], &[])
                .unwrap_err()
                .contains("not declared")
        );
        assert!(
            parse_and_check("INFO/AC[0]>1", &["AC".into()])
                .unwrap_err()
                .contains("subscripts")
        );
        assert!(
            parse_and_check("ID~\"rs[0-9]\"", &[])
                .unwrap_err()
                .contains("metacharacter")
        );
        assert!(
            parse_and_check("ID~QSOMETHING", &[])
                .unwrap_err()
                .contains("string literal")
        );
        assert!(
            parse_and_check("NOSUCH>1", &[])
                .unwrap_err()
                .contains("not declared")
        );
        assert!(parse_and_check("QUAL>", &[]).is_err());
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
