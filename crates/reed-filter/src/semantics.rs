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

/// What [`info_af_outcome`] found, keeping the malformed cases apart.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum InfoAfOutcome {
    /// AC and AN are present and consistent.
    Usable(InfoAfCounts),
    /// AC or AN is absent (or AN is 0): the caller falls back to the
    /// genotype counts, which is what bcftools does when no INFO value
    /// carries the cohort.
    Absent,
    /// AC is present with the wrong number of values. htslib warns once and
    /// `-q/-Q` drops the record; reed keeps D12-13's genotype fallback
    /// (recorded divergence — the expression path does fall back in
    /// bcftools too).
    AcLength {
        /// Values AC carried.
        got: usize,
        /// ALT alleles the record declares.
        want: u32,
    },
    /// ΣAC > AN. htslib prints `[E::bcf_calc_ac] Incorrect AN/AC counts`
    /// and aborts the run (exit 1); the typed `-q/-Q` path must too.
    AcExceedsAn {
        /// Sum of the AC values.
        alt_sum: u64,
        /// INFO/AN.
        an: u64,
    },
    /// The AC values sum past `u64::MAX` (doc_v26 V26-02). The old
    /// `alts.iter().sum()` wrapped in release builds and *panicked* in
    /// debug, so a record like `AC=18446744073709551615,1` was either a
    /// crash or a silently wrong allele frequency. It is a malformed INFO
    /// value, never usable counts.
    AcOverflow {
        /// Values AC carried.
        got: usize,
    },
}

/// The distinguishing reader; the warning-once side effect for an AC
/// length mismatch lives here so every caller shares it.
#[must_use]
pub fn info_af_outcome(info: &[u8], n_alt: u32) -> InfoAfOutcome {
    use reed_io::vcf::info::{InfoIndex, info_fragments};

    // doc_v26 V26-04: one lookup convention for every INFO reader. The
    // expression/head path reads through `InfoIndex` (first match, like
    // htslib); this reader used to keep the *last* AC/AN, so
    // `AC=1;AC=3;AN=4` filtered differently depending on the rule that ran.
    // `classify_line` now rejects duplicate tags outright, and this reader
    // uses the same index so the convention cannot drift again.
    let mut index = InfoIndex::default();
    if index.build(info).is_err() {
        // Unreachable for records that passed `classify_line`; a directly
        // constructed call still must not pick an element silently.
        return InfoAfOutcome::Absent;
    }
    let mut ac: Option<Vec<u64>> = None;
    let mut an: Option<u64> = None;
    for tag in [&b"AC"[..], &b"AN"[..]] {
        let Some(Some(value)) = index.get(info, tag) else {
            continue;
        };
        if tag == b"AC" {
            let Ok(text) = std::str::from_utf8(value) else {
                return InfoAfOutcome::Absent;
            };
            // doc_v26 V26-02: a value list that cannot even match the
            // record's ALT count is rejected before it is materialised, so a
            // multi-megabyte AC field cannot force a proportional allocation.
            // Count first, then parse only when the length is plausible.
            let mut count = 0_usize;
            for part in text.split(',') {
                if count > n_alt as usize {
                    warn_ac_length_once(count, n_alt);
                    return InfoAfOutcome::AcLength {
                        got: count,
                        want: n_alt,
                    };
                }
                // `.` and any non-integer read as no INFO, not as a length
                // mismatch (bcftools' missing-value path).
                if part.parse::<u64>().is_err() {
                    return InfoAfOutcome::Absent;
                }
                count += 1;
            }
            let mut counts = Vec::with_capacity(count);
            for part in text.split(',') {
                let Ok(count_value) = part.parse::<u64>() else {
                    return InfoAfOutcome::Absent;
                };
                counts.push(count_value);
            }
            ac = Some(counts);
        } else {
            an = std::str::from_utf8(value)
                .ok()
                .and_then(|t| t.parse::<u64>().ok());
        }
    }
    let _ = info_fragments;
    let (alts, an) = match (ac, an) {
        (Some(ac), Some(an)) if an > 0 => (ac, an),
        _ => return InfoAfOutcome::Absent,
    };
    if alts.len() != n_alt as usize {
        warn_ac_length_once(alts.len(), n_alt);
        return InfoAfOutcome::AcLength {
            got: alts.len(),
            want: n_alt,
        };
    }
    // doc_v26 V26-02: checked accumulation. A wrapped sum could pass the
    // `> an` test below and produce a wrong AF instead of an error.
    let mut alt_sum: u64 = 0;
    for &count in &alts {
        let Some(next) = alt_sum.checked_add(count) else {
            return InfoAfOutcome::AcOverflow { got: alts.len() };
        };
        alt_sum = next;
    }
    if alt_sum > an {
        return InfoAfOutcome::AcExceedsAn { alt_sum, an };
    }
    let mut per_allele = Vec::with_capacity(alts.len() + 1);
    per_allele.push(an - alt_sum);
    per_allele.extend(alts);
    InfoAfOutcome::Usable(InfoAfCounts { per_allele, an })
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
