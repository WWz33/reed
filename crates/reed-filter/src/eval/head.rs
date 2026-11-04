//! Head-rule evaluation on raw record bytes (§11.2: targets → variant
//! class/max-alt → FILTER → QUAL → INFO, first failure short-circuits; the
//! evaluation order is semantics and must not be reordered).

use reed_core::site::{VariantClass, classify_variant};
use reed_io::vcf::fixed::parse_qual;
use reed_io::vcf::info::InfoIndex;
use reed_io::vcf::line::LineRef;

use crate::plan::{HeadRule, InfoRule};

/// Per-worker scratch: the INFO fragment index of the current record, built
/// once per record and shared by the typed INFO rules and the expressions
/// (A6-18 — the expression path used to re-scan the INFO text per tag).
#[derive(Default)]
pub struct HeadScratch {
    /// The current record's INFO ranges (valid exactly when
    /// [`Self::index_current`] holds).
    pub index: InfoIndex,
    /// Whether `index` belongs to the record being evaluated right now
    /// (rebuilt by [`eval_head`] when any head rule reads INFO; consumers
    /// outside the head stage check this before reusing the ranges).
    pub index_current: bool,
}

impl HeadScratch {
    /// Fresh scratch.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
}

/// Evaluates the head rules in order; increments `hits[i]` for every rule
/// the record passes (cumulative funnel semantics) and returns whether the
/// record survived the head stage.
///
/// `info_trusted` is the plan's single field (doc_v20 S58): once the cohort
/// changed, the built-in AC/AN/AF/MAF/MAC readers fall back to the genotype
/// values instead of the stale INFO. It used to be hard-coded `true` here,
/// which read a full-cohort `INFO/AC` on a subset run (doc_v25 V25-04).
pub fn eval_head(
    line: &LineRef<'_>,
    pos1: u64,
    rules: &[HeadRule],
    hits: &mut [u64],
    scratch: &mut HeadScratch,
    tags: &reed_core::tag_table::TagTable,
    info_trusted: bool,
) -> bool {
    // The counters are aligned with the rule list; a short slice would make
    // the zip below silently skip rules (a caller invariant violation).
    assert_eq!(
        rules.len(),
        hits.len(),
        "invariant: head counters align with the rules"
    );
    let chrom = line.chrom();
    let reference = line.reference();
    let alts = line.alts();
    // One fragment walk per record serves every INFO reader (A6-18).
    scratch.index_current = rules
        .iter()
        .any(|r| matches!(r, HeadRule::Info(_) | HeadRule::Expr { .. }));
    if scratch.index_current {
        scratch.index.build(line.info_text());
    }
    for (rule, hit) in rules.iter().zip(hits.iter_mut()) {
        let ok = match rule {
            HeadRule::Regions { include, exclude } => {
                crate::spec::region::region_allows(include, exclude, chrom, pos1)
            }
            HeadRule::VariantClass(classes) => {
                let primary = match memchr::memchr(b',', alts) {
                    Some(i) => &alts[..i],
                    None => alts,
                };
                let class = classify_variant(reference, primary);
                classes.iter().any(|c| class_matches(class, c))
            }
            HeadRule::ExcludeVariantClass(classes) => {
                let primary = match memchr::memchr(b',', alts) {
                    Some(i) => &alts[..i],
                    None => alts,
                };
                let class = classify_variant(reference, primary);
                !classes.iter().any(|c| class_matches(class, c))
            }
            // doc_v21 S29: plink2's `just-acgt` excludes N and any other
            // non-ACGT allele, on REF and every ALT.
            HeadRule::JustAcgt => {
                let acgt = |allele: &[u8]| {
                    !allele.is_empty()
                        && allele
                            .iter()
                            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'))
                };
                acgt(reference) && alts.split(|&b| b == b',').all(acgt)
            }
            HeadRule::MinAltAlleles(min) => reed_io::vcf::site::count_alts(alts) >= *min,
            HeadRule::MaxAltAlleles(max) => reed_io::vcf::site::count_alts(alts) <= *max,
            HeadRule::BiallelicGate => !alts.contains(&b','),
            // bcftools `-f` matches any `;`-separated value of the column,
            // not the whole text (`bcftools filter` routinely writes
            // `LowQual;LowMQ`; doc_v20 R1). A listed `.` means "unfiltered":
            // htslib maps `.` to filter id -1, which only matches a column
            // with no filter ids at all, so `LowQual;.` must NOT pass `-f .`
            // (verified against bcftools 1.24; v20 review).
            HeadRule::RequireFilter(set) => {
                let column = line.filter_text();
                set.iter().any(|s| {
                    if s.as_slice() == b"." {
                        column == b"."
                    } else {
                        column.split(|&b| b == b';').any(|value| value == s.as_slice())
                    }
                })
            }
            // The INFO preference holds only while AC/AN describe the
            // cohort being decided on (doc_v12 D12-6): the compile drops
            // these rules when the sample set changes, so a rule that is
            // present always reads trusted INFO here.
            HeadRule::MinAfInfo(bound) => non_reference_af_info(line.info_text(), n_alt_of(alts))
                .is_none_or(|freq| freq >= *bound),
            HeadRule::MaxAfInfo(bound) => non_reference_af_info(line.info_text(), n_alt_of(alts))
                .is_none_or(|freq| freq <= *bound),
            HeadRule::MinQual(q) => {
                // A missing QUAL fails (bcftools semantics).
                parse_qual(line.qual_text(), locus(line, pos1))
                    .ok()
                    .flatten()
                    .is_some_and(|v| f64::from(v) >= f64::from(*q))
            }
            HeadRule::Info(rule) => eval_info(line, rule, scratch),
            HeadRule::Expr { include, expr, .. } => {
                // The cohort-trust decision is the plan's field, not a
                // constant: `INFO/AC` used to be routed here by spelling and
                // read the full-cohort value on a subset run (V25-04).
                let mut values =
                    crate::expr::SiteValues::from_line(line, pos1, None, info_trusted);
                values.info_index = Some(&scratch.index);
                values.tags = Some(tags);
                let truth = crate::expr::eval(expr, &values);
                if *include { truth } else { !truth }
            }
        };
        if !ok {
            return false;
        }
        *hit += 1;
    }
    true
}

/// The record's ALT count from its ALT column (`.` = zero, the same rule
/// as `SiteView::n_alt`; doc_v12 D12-10).
fn n_alt_of(alts: &[u8]) -> u32 {
    if alts == b"." {
        return 0;
    }
    u32::try_from(memchr::memchr_iter(b',', alts).count() + 1).unwrap_or(u32::MAX)
}

/// `AC` (summed over ALTs) / `AN` from the INFO text; `None` when either is
/// missing or malformed (the record then defers to the tail rule, matching
/// bcftools' per-record fallback to genotype counts).
///
/// The reader is [`crate::semantics::info_af`]: head and tail now share the
/// malformed-record guard (`ΣAC > AN` reads as "no usable INFO" on both).
#[must_use]
pub fn non_reference_af_info(info: &[u8], n_alt: u32) -> Option<f64> {
    let counts = crate::semantics::info_af(info, n_alt)?;
    Some((counts.an - counts.per_allele[0]) as f64 / counts.an as f64)
}

fn class_matches(class: VariantClass, wanted: &str) -> bool {
    matches!(
        (class, wanted),
        (VariantClass::Snv, "snv")
            | (VariantClass::Mnv, "mnv")
            | (VariantClass::Indel, "indel")
            | (VariantClass::Symbolic, "symbolic")
            | (VariantClass::Ref, "ref")
    )
}

fn eval_info(line: &LineRef<'_>, rule: &InfoRule, scratch: &mut HeadScratch) -> bool {
    let info = line.info_text();
    // An absent tag is a missing operand, exactly like a `.` value: bcftools
    // reports nvalues 0 for both and runs them through the same
    // `missing_logic` table (doc_v23 V23-59).
    let Some(value) = scratch.index.get(info, &rule.tag) else {
        return missing_value_passes(rule);
    };
    let Some(value) = value else {
        return missing_value_passes(rule);
    };
    if rule.number_a {
        // Comma-separated per-ALT values with a quantifier: `any` is true
        // when an element satisfies the comparison, `all` when every
        // parsable element does (the definitions in `predicate.rs:168-176`
        // and architecture.md §11.3; doc_v12 D12-8).
        let mut any = false;
        let mut all = true;
        for part in value.split(|&b| b == b',') {
            match parse_number(part) {
                Some(v) => {
                    let ok = rule.cmp.apply_f64(v);
                    any |= ok;
                    all &= ok;
                }
                // A `.` element counts as missing, so `!=` passes on it for
                // the same reason a wholly missing tag does (V23-59). The
                // right-hand side is a concrete literal, so this is the
                // one-side-missing row of the table, never both.
                None => {
                    let ok = missing_value_passes(rule);
                    any |= ok;
                    all &= ok;
                }
            }
        }
        match rule.quantifier {
            crate::spec::Quantifier::Any => any,
            crate::spec::Quantifier::All => all,
        }
    } else {
        match parse_number(value) {
            Some(v) => rule.cmp.apply_f64(v),
            // The tag is present but its value is `.` or unparsable: bcftools
            // treats that exactly like an absent tag for `!=` (one side
            // missing) and as false for every other operator
            // (doc_v23 V23-59/V23-23).
            None => missing_value_passes(rule),
        }
    }
}

/// A missing observation in an INFO rule: [`crate::semantics::compare_missing`]
/// with the [`InfoRule`](crate::semantics::MissingCtx::InfoRule) context.
/// `on_missing` still overrides the decision when the caller asked for it:
/// an explicit `pass`/`fail` in the spec wins over the parity default.
fn missing_value_passes(rule: &InfoRule) -> bool {
    // The shared truth table (doc_v25 V25-04): the explicit `on_missing`
    // overrides bcftools' `!=`-passes row; the operator table itself is not
    // spelled out here any more.
    crate::semantics::compare_missing(
        rule.cmp.op.into(),
        crate::semantics::MissingSides::Lhs,
        crate::semantics::MissingCtx::InfoRule(rule.on_missing),
    )
}

/// Parses an INFO value as f64 (Integer widened; `.`/malformed → None).
fn parse_number(value: &[u8]) -> Option<f64> {
    std::str::from_utf8(value).ok().and_then(|s| s.parse().ok())
}

/// A cheap error locus: the borrowed CHROM bytes plus the position; the
/// `Locus` (and its `String`) is built only when an error carries it (A6-11).
fn locus<'a>(line: &'a LineRef<'_>, pos1: u64) -> (&'a [u8], u64) {
    (line.chrom(), pos1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spec::OnMissing;

    /// One synthetic record with the given FILTER column against a single
    /// `RequireFilter` rule carrying `list` (empty tags: no INFO rule runs,
    /// so the table is never read).
    fn require_filter_passes(filter_column: &str, list: &[&str]) -> bool {
        let text = format!("chr1\t100\t.\tA\tG\t50.0\t{filter_column}\t.\tGT\t0/0\t0/1");
        let line = LineRef::parse(text.as_bytes()).expect("invariant: 10 columns");
        let rules = vec![HeadRule::RequireFilter(
            list.iter().map(|s| s.as_bytes().to_vec()).collect(),
        )];
        let mut hits = [0];
        let mut scratch = HeadScratch::new();
        eval_head(
            &line,
            100,
            &rules,
            &mut hits,
            &mut scratch,
            &reed_core::tag_table::TagTable::new(),
            true,
        )
    }

    /// R1 (doc_v20): bcftools `-f` keeps a record when the FILTER column
    /// *contains* any listed value, not when the whole column equals one
    /// entry — `bcftools filter` routinely writes `LowQual;LowMQ`, and
    /// `-f LowQual` silently dropped such records.
    /// doc_v25 V25-04: a record with `ΣAC > AN` is malformed in the same
    /// way on the head and the tail — the head's `-q/-Q` reader returns
    /// `None` (the bound passes) instead of computing an AF above 1.
    #[test]
    fn malformed_ac_an_is_no_usable_info_for_the_head() {
        assert_eq!(non_reference_af_info(b"AC=200;AN=100", 1), None);
        assert_eq!(non_reference_af_info(b"AC=2;AN=4", 1), Some(0.5));
    }

    #[test]
    fn r1_filter_list_matches_any_semicolon_separated_value() {
        // The reported bug: a two-value column dropped by a listed member.
        assert!(require_filter_passes("LowQual;LowMQ", &["LowQual"]));
        // The second column value and a non-first listed entry both count.
        assert!(require_filter_passes("LowQual;LowMQ", &["Other", "LowMQ"]));
        assert!(require_filter_passes("A;B;C", &["B"]));
        // Behaviors that were right before the fix and must stay so:
        // a single-value column still matches exactly.
        assert!(require_filter_passes("LowQual", &["LowQual"]));
        // a listed "." matches the unfiltered "." column.
        assert!(require_filter_passes(".", &["."]));
        // PASS and "." never match a listed name, nor names PASS, unless
        // literally listed.
        assert!(!require_filter_passes("PASS", &["LowQual"]));
        assert!(!require_filter_passes("LowQual", &["PASS"]));
        assert!(!require_filter_passes("LowQual;LowMQ", &["PASS", "."]));
        assert!(require_filter_passes("PASS", &["PASS", "."]));
        // a listed name the column does not contain still drops the record.
        assert!(!require_filter_passes("LowQual;LowMQ", &["Other"]));
        // doc_v20 review of R1: a listed "." matches only a column that is
        // exactly ".", never a "." element inside a multi-value column —
        // bcftools maps "." to filter id -1, which matches only unfiltered
        // columns (verified with bcftools 1.24).
        assert!(!require_filter_passes("LowQual;.", &["."]));
        assert!(!require_filter_passes("PASS;.", &["."]));
        assert!(require_filter_passes(".", &["LowQual", "."]));
    }

    /// One record carrying `info` against a single INFO rule `tag op value`.
    fn info_rule_passes_with(
        info: &[u8],
        tag: &str,
        cmp: &str,
        on_missing: OnMissing,
        number_a: bool,
        quantifier: crate::spec::Quantifier,
    ) -> bool {
        let text = format!("chr1\t100\t.\tA\tG\t50.0\tPASS\t{}\tGT\t0/0\t0/1", String::from_utf8_lossy(info));
        let line = LineRef::parse(text.as_bytes()).expect("invariant: 10 columns");
        let rules = vec![HeadRule::Info(InfoRule {
            tag: tag.as_bytes().to_vec(),
            cmp: crate::spec::predicate::Cmp::parse(cmp).expect("invariant: operator"),
            on_missing,
            quantifier,
            number_a,
        })];
        let mut hits = [0];
        let mut scratch = HeadScratch::new();
        eval_head(
            &line,
            100,
            &rules,
            &mut hits,
            &mut scratch,
            &reed_core::tag_table::TagTable::new(),
            true,
        )
    }

    fn info_rule_passes(info: &[u8], tag: &str, cmp: &str, on_missing: OnMissing) -> bool {
        info_rule_passes_with(
            info,
            tag,
            cmp,
            on_missing,
            false,
            crate::spec::Quantifier::Any,
        )
    }

    /// doc_v23 V23-59 (decision A): with a missing observation, `on_missing`
    /// "pass" means bcftools' operator table decides — `!=` keeps the record,
    /// every other operator drops it. `"fail"` overrides to always drop.
    /// Measured against bcftools 1.24-2 (`DP=.` and `DP` absent).
    #[test]
    fn missing_info_values_follow_the_bcftools_operator_table() {
        // `.` value: the tag is present but missing (the V23-23 shape).
        assert!(info_rule_passes(b"DP=.", "DP", "!=3", OnMissing::Pass));
        assert!(!info_rule_passes(b"DP=.", "DP", "==3", OnMissing::Pass));
        assert!(!info_rule_passes(b"DP=.", "DP", ">3", OnMissing::Pass));
        assert!(!info_rule_passes(b"DP=.", "DP", "<3", OnMissing::Pass));
        // Absent tag: bcftools also treats it as a missing operand.
        assert!(info_rule_passes(b".", "DP", "!=3", OnMissing::Pass));
        assert!(!info_rule_passes(b".", "DP", "==3", OnMissing::Pass));
        // A called value is unaffected.
        assert!(info_rule_passes(b"DP=5", "DP", "!=3", OnMissing::Pass));
        assert!(!info_rule_passes(b"DP=5", "DP", "==3", OnMissing::Pass));
        // `fail` is an explicit override, not the parity default.
        assert!(!info_rule_passes(b"DP=.", "DP", "!=3", OnMissing::Fail));
        assert!(!info_rule_passes(b".", "DP", "!=3", OnMissing::Fail));
    }

    /// doc_v23 V23-59, the `Number=A` branch: a `.` element is a missing
    /// operand, so `!=` passes on it and every other operator fails. The
    /// first case isolates the missing element (the called `3` fails
    /// `!=3`), so only the operator table can make it pass.
    #[test]
    fn number_a_missing_elements_follow_the_operator_table() {
        use crate::spec::Quantifier::{All, Any};
        // `AD=3,.` with `!=3`: the 3 fails, the `.` passes -> any true.
        assert!(info_rule_passes_with(
            b"AD=3,.",
            "AD",
            "!=3",
            OnMissing::Pass,
            true,
            Any
        ));
        // `==` and ordering fail on the missing element; `AD=5,.` ==3 has
        // no passing element either.
        assert!(!info_rule_passes_with(
            b"AD=5,.",
            "AD",
            "==3",
            OnMissing::Pass,
            true,
            Any
        ));
        assert!(!info_rule_passes_with(
            b"AD=1,.",
            "AD",
            ">3",
            OnMissing::Pass,
            true,
            Any
        ));
        // A wholly missing element list: same table.
        assert!(info_rule_passes_with(
            b"AD=.",
            "AD",
            "!=3",
            OnMissing::Pass,
            true,
            Any
        ));
        // `all` requires every element: the missing one passes `!=` but the
        // called 3 does not, so the rule fails; with 5 it passes.
        assert!(!info_rule_passes_with(
            b"AD=3,.",
            "AD",
            "!=3",
            OnMissing::Pass,
            true,
            All
        ));
        assert!(info_rule_passes_with(
            b"AD=5,.",
            "AD",
            "!=3",
            OnMissing::Pass,
            true,
            All
        ));
        // `fail` overrides the table on the missing element too.
        assert!(!info_rule_passes_with(
            b"AD=3,.",
            "AD",
            "!=3",
            OnMissing::Fail,
            true,
            Any
        ));
    }
}
