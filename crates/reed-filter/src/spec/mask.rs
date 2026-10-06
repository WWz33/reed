//! Genotype mask rules (§11.1 `[[genotype.mask_if]]`): any hit sets the
//! sample's GT to `./.` while every other FORMAT field is preserved.

use serde::{Deserialize, Serialize};

use super::predicate::{Cmp, OnMissing};

/// One mask rule over a FORMAT field.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct GenotypeMaskRule {
    /// The FORMAT field (validated against the header at compile time).
    pub field: String,
    /// The comparison, e.g. `"< 5"`.
    pub when: String,
    /// Missing-value behaviour (default `keep`: a missing value never
    /// triggers the mask — bcftools semantics; §11.1).
    #[serde(default = "default_mask_on_missing")]
    pub on_missing: OnMissing,
}

/// The mask-context default (`Keep`, not the enum-wide `Pass`).
fn default_mask_on_missing() -> OnMissing {
    OnMissing::Keep
}

impl GenotypeMaskRule {
    /// Parses `when`.
    ///
    /// # Errors
    ///
    /// [`crate::FilterError::Spec`] when the comparison is malformed.
    pub fn parsed(&self) -> Result<Cmp, crate::FilterError> {
        Cmp::parse(&self.when)
    }
}
