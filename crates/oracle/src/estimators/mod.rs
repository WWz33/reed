//! Estimators: sealed definitions of population-genetics statistics.

pub mod allele_counts;
pub mod dxy;
pub mod fst;
pub mod fstats;
pub mod hwe;
pub mod ld;
pub mod pbs;
pub mod pi;
pub mod site_counts;
pub mod spectrum;

use std::fmt;

use crate::manifest::{Ascertainment, SiteCompleteness};
use crate::population::PopulationView;
use crate::window::WindowBuffer;

/// An estimator refused because the input manifest cannot support its
/// denominator.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Unavailable {
    pub reason: &'static str,
}

impl fmt::Display for Unavailable {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "estimator unavailable: {}", self.reason)
    }
}

impl std::error::Error for Unavailable {}

pub(crate) const ERR_POP_MISMATCH: &str =
    "population view does not match the window's sample count \
     (see PopulationView::validate_against)";

pub(crate) const ERR_GROUP_RANGE: &str = "group index out of range for this population view";

/// Resolves a group's samples after validating the view against the window.
///
/// Every estimator entry point calls this instead of indexing directly, so a
/// view built for a different sample set is a boundary error rather than a
/// deep index panic (§12.2-8).
pub(crate) fn samples_of<'a>(
    pops: &'a PopulationView,
    window: &WindowBuffer,
    group: u32,
) -> Result<&'a [u32], Unavailable> {
    pops.validate_against(window).map_err(|_| Unavailable {
        reason: ERR_POP_MISMATCH,
    })?;
    pops.get(group).ok_or(Unavailable {
        reason: ERR_GROUP_RANGE,
    })
}

/// Bias note for statistics that are **invariant to the presence of
/// monomorphic sites** (their valid site/pair sets are identical with or
/// without them): `VariantOnly(Complete)` yields values identical to
/// all-sites input; only `Filtered` — or `Unknown`, which cannot be
/// distinguished from it — can bias them.
pub(crate) fn monomorphic_invariant_note(
    manifest: SiteCompleteness,
    note: &'static str,
) -> Option<&'static str> {
    match manifest {
        SiteCompleteness::WithInvariants
        | SiteCompleteness::VariantOnly(Ascertainment::Complete) => None,
        SiteCompleteness::VariantOnly(_) => Some(note),
    }
}

/// Bias note for statistics whose site mean **includes monomorphic sites**
/// under all-sites input: any variant-only input changes the quantity
/// (site-density dependence), not merely its precision.
pub(crate) fn variant_only_note(
    manifest: SiteCompleteness,
    note: &'static str,
) -> Option<&'static str> {
    match manifest {
        SiteCompleteness::WithInvariants => None,
        SiteCompleteness::VariantOnly(_) => Some(note),
    }
}
