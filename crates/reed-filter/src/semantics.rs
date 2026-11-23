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

/// Per-allele AC counts with REF filled in, plus AN (doc_v30 R02: the type
/// now lives in `reed-core`, so `stats` and the filter cannot drift).
pub use reed_core::count::InfoCounts as InfoAfCounts;

/// What [`info_af_outcome`] found, keeping the malformed cases apart.
pub use reed_core::count::InfoCountsOutcome as InfoAfOutcome;

/// Reads INFO/AC and INFO/AN under htslib's `bcf_calc_ac` contract
/// (vcfutils.c:56-64), shared by `-q/-Q`, the AF/AC built-ins and `--maf`:
///
/// * AC must carry one value per ALT allele (`AC` length == `n_alt`), else
///   the record has no usable INFO (warn once);
/// * AN must be present and > 0;
/// * the summed ALT counts must not exceed AN — a record like
///   `AC=200;AN=100` is malformed; htslib reports `Incorrect AN/AC counts`
///   and aborts the run.
///
/// The REF count is `AN − ΣALT`, exactly how `bcf_calc_ac` fills `ac[0]`.
///
/// This lenient entry point (used by the tail's genotype-fallback path)
/// folds every malformed case into `None`. [`info_af_outcome`] keeps the
/// distinction the typed `-q/-Q` head needs.
#[must_use]
pub fn info_af(info: &[u8], n_alt: u32) -> Option<InfoAfCounts> {
    match info_af_outcome(info, n_alt) {
        InfoAfOutcome::Usable(counts) => Some(counts),
        // bcftools' expression/fallback path recomputes from the genotypes
        // when INFO is unusable; the mismatch warning has already fired.
        InfoAfOutcome::Absent
        | InfoAfOutcome::AcLength { .. }
        | InfoAfOutcome::AcExceedsAn { .. }
        | InfoAfOutcome::AcOverflow { .. } => None,
    }
}

/// The distinguishing reader; the warning-once side effect for an AC
/// length mismatch lives here so every caller shares it.
#[must_use]
pub fn info_af_outcome(info: &[u8], n_alt: u32) -> InfoAfOutcome {
    use reed_core::vcf::info::InfoIndex;

    // doc_v26 V26-04: one lookup convention for every INFO reader. The
    // expression/head path reads through `InfoIndex` (first match, like
    // htslib); this reader used to keep the *last* AC/AN, so
    // `AC=1;AC=3;AN=4` filtered differently depending on the rule that ran.
    // `classify_line` now rejects duplicate tags outright, and this reader
    // uses the same index so the convention cannot drift again.
    //
    // doc_v28 V28-13: a caller that already holds the record's index passes
    // it to [`info_af_indexed`] and this walk disappears.
    let mut index = InfoIndex::default();
    if index.build(info).is_err() {
        // Unreachable for records that passed `classify_line`; a directly
        // constructed call still must not pick an element silently.
        return InfoAfOutcome::Absent;
    }
    info_af_indexed(&index, info, n_alt)
}

/// [`info_af_outcome`] over an index the caller already built for `info`
/// (doc_v28 V28-13). The index must have been built from `info`; a caller
/// cannot satisfy that with a stale table, because the ranges are offsets
/// into this column.
///
/// A failed (empty) index yields [`InfoAfOutcome::Absent`], the same
/// conservative answer as a rebuild failure.
///
/// The AC/AN validation itself is `reed_core::count::InfoCounts::from_info`
/// (doc_v30 R02); this function only supplies the two INFO values and the
/// warn-once side effect.
#[must_use]
pub fn info_af_indexed(
    index: &reed_core::vcf::info::InfoIndex,
    info: &[u8],
    n_alt: u32,
) -> InfoAfOutcome {
    let ac = index.get(info, b"AC").flatten();
    let an = index.get(info, b"AN").flatten();
    let outcome = reed_core::count::InfoCounts::from_info(ac, an, n_alt);
    if let InfoAfOutcome::AcLength { got, want } = outcome {
        warn_ac_length_once(got, want);
    }
    outcome
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
                !compare_missing(
                    CmpOp::Ne,
                    MissingSides::Lhs,
                    MissingCtx::InfoRule(on_missing)
                ),
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

    /// doc_v26 V26-02: AC values that sum past `u64::MAX` are a malformed
    /// value, not a wrapped frequency. `18446744073709551615,1` used to panic
    /// in debug and wrap in release.
    #[test]
    fn info_af_rejects_ac_overflow_instead_of_wrapping() {
        let case = b"AC=18446744073709551615,1;AN=100";
        assert_eq!(
            info_af_outcome(case, 2),
            InfoAfOutcome::AcOverflow { got: 2 },
            "a wrapped sum could pass the > AN test and report a wrong AF"
        );
        assert!(info_af(case, 2).is_none(), "no usable counts");
        // The boundary itself is fine: a sum that equals u64::MAX is
        // representable and then simply exceeds AN.
        let edge = b"AC=18446744073709551614,1;AN=100";
        assert!(matches!(
            info_af_outcome(edge, 2),
            InfoAfOutcome::AcExceedsAn { .. }
        ));
    }

    /// doc_v26 V26-04: AC/AN are read through the same first-match index as
    /// every other INFO reader; a duplicated column yields no counts (the
    /// record is rejected earlier by `classify_line`).
    #[test]
    fn info_af_escalates_duplicate_tags_to_absent() {
        assert_eq!(info_af_outcome(b"AC=1;AC=3;AN=4", 1), InfoAfOutcome::Absent);
    }

    /// doc_v26 V26-02: a value list that cannot match the ALT count is
    /// rejected while counting, before any per-element parse or vector.
    #[test]
    fn info_af_bounds_the_ac_list_before_parsing() {
        // 4 values for a 1-ALT record: the length mismatch is detected at the
        // second value and reported with the count seen so far.
        let got = info_af_outcome(b"AC=1,2,3,4;AN=100", 1);
        assert!(
            matches!(got, InfoAfOutcome::AcLength { want: 1, .. }),
            "{got:?}"
        );
    }
}
