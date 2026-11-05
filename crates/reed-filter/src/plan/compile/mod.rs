//! `compile()` split by section (doc_v25 08 §3 ④).
//!
//! The entry keeps only the cross-section wiring; each private submodule
//! owns one section in evaluation order, so a new rule lands where it is
//! evaluated instead of growing one 600-line function. Pure refactor: the
//! sections run in the same order and validate in the same order as before.

mod finish;
mod head;
mod mask;
mod tail;

use reed_core::tag_table::TagTable;

use super::{FilterPlan, FmtExprField, HeadRule, IntField, MaskRule};
use crate::FilterError;
use crate::spec::FilterSpec;

/// Compiles the spec against the input header (doc_v25 08 §3 ④: the entry
/// keeps only the cross-section wiring; head/mask/tail/finish are private
/// submodules in evaluation order).
///
/// `biallelic_gate` is the filter driver's default (§11.6, plink2-style):
/// push [`HeadRule::BiallelicGate`] unless the spec allows multiallelic
/// input. The stats driver passes `false` — its multiallelic handling is per
/// reporter, not a record gate.
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
    let head = head::compile(spec, tags, biallelic_gate)?;
    let mut rules = mask::compile(spec, tags)?;
    let tail = tail::compile(spec, tags, &head.format_tags, &mut rules)?;
    finish::compile(spec, tags, command_line, head, rules, tail)
}

/// What the head section leaves for the later stages.
struct HeadOutcome {
    head: Vec<HeadRule>,
    /// Whether the cohort changed (subset/mask/`--mind`): the tail reads the
    /// genotype counts instead of the stale INFO.
    cohort_changed: bool,
    /// FORMAT tag names, for the `--set-gts .` expression parser.
    format_tags: Vec<String>,
}

/// The accumulating rule/column vectors. The tail section appends the
/// `--set-gts .` mask rules and their FORMAT columns, so the mask section
/// cannot own them exclusively.
struct RuleSet {
    mask: Vec<MaskRule>,
    int_fields: Vec<IntField>,
    float_fields: Vec<Vec<u8>>,
    fmt_expr_fields: Vec<FmtExprField>,
}
