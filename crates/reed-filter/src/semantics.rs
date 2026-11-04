//! The filter semantics that used to be transcribed per call site (doc_v25
//! V25-04): the missing-value comparison table and the INFO AC/AN reader.
//!
//! bcftools' behaviour here is a *function*, not prose: the expression
//! evaluator, the typed INFO rules and the columnar FORMAT mask all ask the
//! same questions, and the answers must not drift apart.

use crate::expr::CmpOp;
use crate::spec::OnMissing;
use crate::spec::predicate::Op;

/// Which sides of a comparison are missing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MissingSides {
    /// The observation (left-hand side) is missing, the threshold is a
    /// concrete literal.
    Lhs,
    /// The observation exists but the right-hand side is missing (an
    /// expression comparing two record fields).
    Rhs,
    /// Both sides are missing.
    Both,
}

/// The context a missing comparison runs in. The rule differs by context,
/// and the difference is contract, not accident:
///
/// * [`Expr`](Self::Expr) is bcftools' `missing_logic` (filter.c:4449-4455):
///   `!=` passes with one side missing, `==`/`!=` pass with both missing,
///   the ordering operators always fail;
/// * [`InfoRule`](Self::InfoRule) lets an explicit `on_missing` override the
///   same table (`pass`/`fail` in the spec);
/// * [`ColumnarMask`](Self::ColumnarMask) is the FORMAT mask (`--min-dp`,
///   `mask_if` field rules), whose default `keep` means "a missing value
///   never triggers the mask" — deliberately stricter than `-e`, where a
///   missing `DP` satisfies `!= 3` (doc_v25 03 §2).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MissingCtx {
    /// `-i/-e` expression: bcftools' `missing_logic`.
    Expr,
    /// Typed INFO rule (`-q/-Q` short-circuits, `[[site.info]]`).
    InfoRule(OnMissing),
    /// Columnar FORMAT mask (`--min-dp`, `[[genotype.mask_if]]` field rule).
    ColumnarMask(OnMissing),
}

/// Whether a comparison with missing operands passes, by context.
///
/// This is the single implementation of the truth table; every consumer
/// passes its context instead of spelling the table out again.
#[must_use]
pub fn compare_missing(op: CmpOp, sides: MissingSides, ctx: MissingCtx) -> bool {
    match ctx {
        MissingCtx::Expr => match sides {
            // `missing_logic[2]`: only `==` and `!=` pass when both sides
            // are missing.
            MissingSides::Both => matches!(op, CmpOp::Eq | CmpOp::Ne),
            // `missing_logic[1]`: only `!=` passes with one side missing.
            MissingSides::Lhs | MissingSides::Rhs => op == CmpOp::Ne,
        },
        MissingCtx::InfoRule(on_missing) => match on_missing {
            // The INFO-rule default means "do not drop"; within it the
            // operator table decides. `keep`/`mask` are rejected for INFO
            // rules at compile time, but "do not drop" is their reading.
            OnMissing::Pass | OnMissing::Keep => op == CmpOp::Ne,
            OnMissing::Fail | OnMissing::Mask => false,
        },
        // The mask fires exactly when the caller asked missing values to
        // trigger it.
        MissingCtx::ColumnarMask(on_missing) => on_missing == OnMissing::Mask,
    }
}

impl From<Op> for CmpOp {
    fn from(op: Op) -> Self {
        match op {
            Op::Eq => Self::Eq,
            Op::Ne => Self::Ne,
            Op::Lt => Self::Lt,
            Op::Le => Self::Le,
            Op::Gt => Self::Gt,
            Op::Ge => Self::Ge,
        }
    }
}

/// Per-allele AC counts with REF filled in, plus AN (`None` = no usable
/// INFO).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InfoAfCounts {
    /// REF first, then one count per ALT: `n_alt + 1` entries.
    pub per_allele: Vec<u64>,
    /// INFO/AN (always > 0 when `Some`).
    pub an: u64,
}

/// Reads INFO/AC and INFO/AN under htslib's `bcf_calc_ac` contract
/// (vcfutils.c:56-64), shared by `-q/-Q`, the AF/AC built-ins and `--maf`:
///
/// * AC must carry one value per ALT allele (`AC` length == `n_alt`), else
///   the record has no usable INFO (warn once);
/// * AN must be present and > 0;
/// * the summed ALT counts must not exceed AN — a record like
///   `AC=200;AN=100` is malformed and both the head and the tail treat it as
///   "no usable INFO" (the tail used to guard this alone).
///
/// The REF count is `AN − ΣALT`, exactly how `bcf_calc_ac` fills `ac[0]`.
#[must_use]
pub fn info_af(info: &[u8], n_alt: u32) -> Option<InfoAfCounts> {
    use reed_io::vcf::info::info_fragments;

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
    let (alts, an) = match (ac, an) {
        (Some(ac), Some(an)) if an > 0 => (ac, an),
        _ => return None,
    };
    if alts.len() != n_alt as usize {
        warn_ac_length_once(alts.len(), n_alt);
        return None;
    }
    let alt_sum: u64 = alts.iter().sum();
    if alt_sum > an {
        return None;
    }
    let mut per_allele = Vec::with_capacity(alts.len() + 1);
    per_allele.push(an - alt_sum);
    per_allele.extend(alts);
    Some(InfoAfCounts { per_allele, an })
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

#[cfg(test)]
mod tests {
    use super::*;

    fn flags(op: CmpOp) -> [bool; 4] {
        [
            compare_missing(op, MissingSides::Lhs, MissingCtx::Expr),
            compare_missing(op, MissingSides::Rhs, MissingCtx::Expr),
            compare_missing(op, MissingSides::Both, MissingCtx::Expr),
            compare_missing(op, MissingSides::Lhs, MissingCtx::InfoRule(OnMissing::Pass)),
        ]
    }

    /// bcftools `missing_logic` (filter.c:4449-4455): the expression table
    /// (`!=` passes on one side, `==`/`!=` on both) and the INFO-rule
    /// default (same `!=` row), which the explicit `on_missing` overrides.
    #[test]
    fn the_expression_table_matches_filter_c() {
        assert_eq!(flags(CmpOp::Eq), [false, false, true, false]);
        assert_eq!(flags(CmpOp::Ne), [true, true, true, true]);
        for op in [CmpOp::Lt, CmpOp::Le, CmpOp::Gt, CmpOp::Ge] {
            assert_eq!(flags(op), [false, false, false, false], "{op:?}");
        }
        for on_missing in [OnMissing::Fail, OnMissing::Mask] {
            assert!(
                !compare_missing(CmpOp::Ne, MissingSides::Lhs, MissingCtx::InfoRule(on_missing)),
                "{on_missing:?} must veto the != row"
            );
        }
    }

    /// The mask context: only `on_missing = mask` fires on a missing value;
    /// the default `keep` is the documented difference from `-e`.
    #[test]
    fn the_mask_context_needs_an_explicit_opt_in() {
        for op in [
            CmpOp::Eq,
            CmpOp::Ne,
            CmpOp::Lt,
            CmpOp::Le,
            CmpOp::Gt,
            CmpOp::Ge,
        ] {
            assert!(!compare_missing(
                op,
                MissingSides::Lhs,
                MissingCtx::ColumnarMask(OnMissing::Keep)
            ));
            assert!(compare_missing(
                op,
                MissingSides::Lhs,
                MissingCtx::ColumnarMask(OnMissing::Mask)
            ));
        }
    }

    #[test]
    fn spec_ops_map_to_expression_ops() {
        assert_eq!(CmpOp::from(Op::Eq), CmpOp::Eq);
        assert_eq!(CmpOp::from(Op::Ne), CmpOp::Ne);
        assert_eq!(CmpOp::from(Op::Lt), CmpOp::Lt);
        assert_eq!(CmpOp::from(Op::Le), CmpOp::Le);
        assert_eq!(CmpOp::from(Op::Gt), CmpOp::Gt);
        assert_eq!(CmpOp::from(Op::Ge), CmpOp::Ge);
    }

    /// The shared AC/AN reader: REF is `AN − ΣALT`; a length mismatch, an
    /// absent AN, AN=0 and `ΣALT > AN` all read as "no usable INFO".
    #[test]
    fn info_af_fills_ref_and_rejects_malformed_records() {
        let got = info_af(b"AC=2,1;AN=6", 2).expect("usable");
        assert_eq!(got.per_allele, [3, 2, 1]);
        assert_eq!(got.an, 6);

        assert!(info_af(b"AC=2,1;AN=6", 1).is_none(), "AC length mismatch");
        assert!(info_af(b"AC=2;AN=0", 1).is_none(), "AN=0");
        assert!(info_af(b"AC=2", 1).is_none(), "no AN");
        assert!(info_af(b"AN=4", 1).is_none(), "no AC");
        assert!(info_af(b"AC=5;AN=4", 1).is_none(), "ΣALT > AN");
        assert!(info_af(b"AC=.;AN=4", 1).is_none(), "missing AC element");
    }
}