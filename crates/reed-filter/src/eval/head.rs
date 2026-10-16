//! Head-rule evaluation on raw record bytes (§11.2: targets → variant
//! class/max-alt → FILTER → QUAL → INFO, first failure short-circuits; the
//! evaluation order is semantics and must not be reordered).

use reed_core::site::{VariantClass, classify_variant};
use reed_io::vcf::fixed::parse_qual;
use reed_io::vcf::info::{InfoIndex, info_fragments};
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
pub fn eval_head(
    line: &LineRef<'_>,
    pos1: u64,
    rules: &[HeadRule],
    hits: &mut [u64],
    scratch: &mut HeadScratch,
    tags: &reed_core::tag_table::TagTable,
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
            HeadRule::MinAltAlleles(min) => reed_io::vcf::site::count_alts(alts) >= *min,
            HeadRule::MaxAltAlleles(max) => reed_io::vcf::site::count_alts(alts) <= *max,
            HeadRule::BiallelicGate => !alts.contains(&b','),
            HeadRule::RequireFilter(set) => set.iter().any(|s| s == line.filter_text()),
            // The INFO preference holds only while AC/AN describe the
            // cohort being decided on (doc_v12 D12-6): the compile drops
            // these rules when the sample set changes, so a rule that is
            // present always reads trusted INFO here.
            HeadRule::MinAfInfo(bound) => {
                non_reference_af_info(line.info_text(), n_alt_of(alts))
                    .is_none_or(|freq| freq >= *bound)
            }
            HeadRule::MaxAfInfo(bound) => {
                non_reference_af_info(line.info_text(), n_alt_of(alts))
                    .is_none_or(|freq| freq <= *bound)
            }
            HeadRule::MinQual(q) => {
                // A missing QUAL fails (bcftools semantics).
                parse_qual(line.qual_text(), locus(line, pos1))
                    .ok()
                    .flatten()
                    .is_some_and(|v| f64::from(v) >= f64::from(*q))
            }
            HeadRule::Info(rule) => eval_info(line, rule, scratch),
            HeadRule::Expr { include, expr, .. } => {
                let mut values = crate::expr::SiteValues::from_line(line, pos1, None);
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

/// Per-allele INFO/AC counts and INFO/AN (bcftools' `bcf_calc_ac` prefers
/// them over the genotype tallies; docs/rounds/v5/02-cli.md §3.2). Returns
/// `(counts, an)` with the counts in allele order (REF first), or `None`
/// when either tag is absent, malformed, or when the AC vector does not
/// have one value per ALT allele.
///
/// The length check is htslib's: `bcf_calc_ac` refuses a record whose
/// `ac_len != n_allele - 1` and lets the caller fall back to the genotype
/// counts (`vcfutils.c:56-64`, with a once-per-run warning; doc_v12 D12-13).
/// `n_alt` is the record's ALT count.
#[must_use]
pub fn allele_counts_info(info: &[u8], n_alt: u32) -> Option<(Vec<u64>, u64)> {
    let mut ac: Option<Vec<u64>> = None;
    let mut an: Option<u64> = None;
    for (tag, value) in info_fragments(info) {
        let Some(value) = value else { continue };
        match tag {
            b"AC" => {
                let text = std::str::from_utf8(value).ok()?;
                let mut counts = Vec::new();
                for part in text.split(',') {
                    if part == "." {
                        return None;
                    }
                    counts.push(part.parse::<u64>().ok()?);
                }
                ac = Some(counts);
            }
            b"AN" => {
                an = std::str::from_utf8(value).ok()?.parse::<u64>().ok();
            }
            _ => {}
        }
    }
    match (ac, an) {
        (Some(ac), Some(an)) if an > 0 => {
            if ac.len() != n_alt as usize {
                warn_ac_length_once(ac.len(), n_alt);
                return None;
            }
            Some((ac, an))
        }
        _ => None,
    }
}

/// The record's ALT count from its ALT column (`.` = zero, the same rule
/// as `SiteView::n_alt`; doc_v12 D12-10).
fn n_alt_of(alts: &[u8]) -> u32 {
    if alts == b"." {
        return 0;
    }
    u32::try_from(memchr::memchr_iter(b',', alts).count() + 1).unwrap_or(u32::MAX)
}

/// htslib warns once per run about a mismatched AC vector and then stays
/// quiet (`bcf_calc_ac`'s `warned` static, vcfutils.c:58-63).
fn warn_ac_length_once(got: usize, want: u32) {
    use std::sync::atomic::{AtomicBool, Ordering};
    static WARNED: AtomicBool = AtomicBool::new(false);
    if !WARNED.swap(true, Ordering::Relaxed) {
        tracing::warn!(
            "Incorrect number of AC fields ({got} values for {want} ALT alleles); \
             falling back to the genotype counts (further occurrences are not reported)"
        );
    }
}

/// `AC` (summed over ALTs) / `AN` from the INFO text; `None` when either is
/// missing or malformed (the record then defers to the tail rule, matching
/// bcftools' per-record fallback to genotype counts).
#[must_use]
pub fn non_reference_af_info(info: &[u8], n_alt: u32) -> Option<f64> {
    let (counts, an) = allele_counts_info(info, n_alt)?;
    let sum: u64 = counts.iter().sum();
    if an == 0 {
        return None;
    }
    Some(sum as f64 / an as f64)
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
    let Some(value) = scratch.index.get(info, &rule.tag) else {
        return rule.on_missing == crate::spec::OnMissing::Pass;
    };
    let Some(value) = value else {
        return rule.on_missing == crate::spec::OnMissing::Pass;
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
                None => all = false,
            }
        }
        match rule.quantifier {
            crate::spec::Quantifier::Any => any,
            crate::spec::Quantifier::All => all,
        }
    } else {
        parse_number(value).is_some_and(|v| rule.cmp.apply_f64(v))
    }
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
