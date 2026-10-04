//! Estimators: sealed definitions of population-genetics statistics.

pub mod allele_counts;
pub mod dxy;
pub mod fst;
pub mod fstats;
pub mod hwe;
pub mod ld;
pub mod pbs;
pub mod pi;
pub mod spectrum;

use std::fmt;

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
