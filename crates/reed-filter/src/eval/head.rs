//! Head-rule evaluation on raw record bytes (§11.2: targets → variant
//! class/max-alt → FILTER → QUAL → INFO, first failure short-circuits; the
//! evaluation order is semantics and must not be reordered).

use reed_core::site::{classify_variant, VariantClass};
use reed_io::vcf::fixed::parse_qual;
use reed_io::vcf::info::info_fragments;
use reed_io::vcf::line::LineRef;

use crate::plan::{HeadRule, InfoRule};

/// Per-worker scratch: the INFO fragment ranges of the current record
/// (offsets into the record's INFO text — no per-record allocation).
#[derive(Default)]
pub struct HeadScratch {
    fragments: Vec<FragRange>,
}

#[derive(Clone)]
struct FragRange {
    tag: std::ops::Range<usize>,
    value: Option<std::ops::Range<usize>>,
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
) -> bool {
    // The counters are aligned with the rule list; a short slice would make
    // the zip below silently skip rules (a caller invariant violation).
    assert_eq!(rules.len(), hits.len(), "invariant: head counters align with the rules");
    let chrom = line.chrom();
    let reference = line.reference();
    let alts = line.alts();
    for (rule, hit) in rules.iter().zip(hits.iter_mut()) {
        let ok = match rule {
            HeadRule::Regions { include, exclude } => {
                let included = include.is_empty()
                    || include.iter().any(|(c, s, e)| {
                        c.as_bytes() == chrom
                            && s.is_none_or(|s| pos1 >= s)
                            && e.is_none_or(|e| pos1 <= e)
                    });
                included
                    && !exclude.iter().any(|(c, s, e)| {
                        c.as_bytes() == chrom
                            && s.is_none_or(|s| pos1 >= s)
                            && e.is_none_or(|e| pos1 <= e)
                    })
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
            HeadRule::MinAltAlleles(min) => {
                let n = u32::try_from(alts.iter().filter(|&&b| b == b',').count() + 1)
                    .unwrap_or(u32::MAX);
                n >= *min
            }
            HeadRule::MaxAltAlleles(max) => {
                let n = u32::try_from(alts.iter().filter(|&&b| b == b',').count() + 1)
                    .unwrap_or(u32::MAX);
                n <= *max
            }
            HeadRule::RequireFilter(set) => set.iter().any(|s| s == line.filter_text()),
            HeadRule::MinQual(q) => {
                // A missing QUAL fails (bcftools semantics).
                parse_qual(line.qual_text(), &locus(line, pos1))
                    .ok()
                    .flatten()
                    .is_some_and(|v| f64::from(v) >= f64::from(*q))
            }
            HeadRule::Info(rule) => eval_info(line, rule, scratch),
            HeadRule::Expr { include, expr, .. } => {
                let values = crate::expr::SiteValues::from_line(line, pos1, None);
                let truth = crate::expr::eval(expr, &values);
                if *include {
                    truth
                } else {
                    !truth
                }
            }
        };
        if !ok {
            return false;
        }
        *hit += 1;
    }
    true
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
    // One fragment walk per record serves every rule; the ranges are
    // rebuilt per record (INFO text changes from record to record).
    scratch.fragments.clear();
    if info != b"." {
        for (tag, value) in info_fragments(info) {
            let tag_start = tag.as_ptr() as usize - info.as_ptr() as usize;
            let value = value.map(|v| {
                let start = v.as_ptr() as usize - info.as_ptr() as usize;
                start..start + v.len()
            });
            scratch.fragments.push(FragRange {
                tag: tag_start..tag_start + tag.len(),
                value,
            });
        }
    }
    let Some(frag) = scratch.fragments.iter().find(|f| info[f.tag.clone()] == *rule.tag) else {
        return rule.on_missing == crate::spec::OnMissing::Pass;
    };
    let Some(value) = frag.value.clone() else {
        return rule.on_missing == crate::spec::OnMissing::Pass;
    };
    let value = &info[value];
    if rule.number_a {
        // Comma-separated per-ALT values with a quantifier.
        let mut any = false;
        let mut all = true;
        for part in value.split(|&b| b == b',') {
            match parse_number(part) {
                Some(v) => {
                    any = true;
                    all &= rule.cmp.apply_f64(v);
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

fn locus(line: &LineRef<'_>, pos1: u64) -> reed_io::Locus {
    reed_io::Locus {
        contig_name: String::from_utf8_lossy(line.chrom()).into_owned(),
        pos1,
    }
}
