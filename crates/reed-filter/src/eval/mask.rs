//! Genotype mask evaluation (§11.2 Mask stage): per kept sample, the mask
//! rules run in declaration order and the first hit attributes the mask.

use crate::plan::MaskRule;
use crate::spec::OnMissing;

/// Evaluates the mask rules for one sample over its scanned field values.
///
/// `values[i]` corresponds to `rules[i]`'s field observation (`None` =
/// missing). Returns `(masked, rule_index)` — the first rule that fires
/// wins the attribution; `None` when the sample is not masked.
#[must_use]
pub fn eval_mask(values: &[Option<f64>], rules: &[MaskRule]) -> Option<usize> {
    debug_assert_eq!(values.len(), rules.len(), "invariant: one value per rule");
    for (i, (rule, value)) in rules.iter().zip(values).enumerate() {
        let hit = match value {
            Some(v) => rule.cmp.apply_f64(*v),
            None => rule.on_missing == OnMissing::Mask,
        };
        if hit {
            return Some(i);
        }
    }
    None
}
