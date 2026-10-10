//! Result types: exact integer counts first, `f64` derived and auditable.

use crate::estimators::dxy::DxyEstimator;
use crate::estimators::fst::FstEstimator;
use crate::estimators::pi::PiEstimator;
use crate::window::ContigId;

/// Pair-based accounting at a set of sites (Korunes & Samuk 2021 style).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PairCounts {
    /// Pairwise differences summed over sites.
    pub diffs: u64,
    /// Comparable pairs summed over sites.
    pub comparable: u64,
    /// Pairs lost to missing alleles summed over sites.
    pub missing: u64,
}

/// Per-group nucleotide diversity.
#[derive(Clone, Debug)]
pub struct PiResult {
    pub estimator: PiEstimator,
    pub contig: ContigId,
    pub start: u64,
    pub end: u64,
    pub group: u32,
    pub counts: PairCounts,
    /// Numerator of the `AccessibleBp` estimator: sum over sites of the
    /// per-site mean pairwise difference (sites with no comparable pair
    /// contribute zero).
    pub sum_site_means: f64,
    pub n_variant_sites: u64,
    pub n_invariant_sites: u64,
    pub bias: Option<&'static str>,
}

impl PiResult {
    /// Window span in base pairs.
    pub fn window_bp(&self) -> u64 {
        self.end - self.start
    }

    /// The estimate, or `None` when the denominator is empty.
    pub fn value(&self) -> Option<f64> {
        match self.estimator {
            PiEstimator::PerPairComparable => (self.counts.comparable > 0)
                .then(|| self.counts.diffs as f64 / self.counts.comparable as f64),
            PiEstimator::AccessibleBp => {
                let has_comparable = self.counts.comparable > 0;
                (has_comparable && self.window_bp() > 0)
                    .then(|| self.sum_site_means / self.window_bp() as f64)
            }
        }
    }
}

/// Between-group absolute divergence.
#[derive(Clone, Debug)]
pub struct DxyResult {
    pub estimator: DxyEstimator,
    pub contig: ContigId,
    pub start: u64,
    pub end: u64,
    pub group_a: u32,
    pub group_b: u32,
    pub counts: PairCounts,
    pub sum_site_means: f64,
    pub n_variant_sites: u64,
    pub n_invariant_sites: u64,
    pub bias: Option<&'static str>,
}

impl DxyResult {
    /// Window span in base pairs.
    pub fn window_bp(&self) -> u64 {
        self.end - self.start
    }

    /// The estimate, or `None` when the denominator is empty.
    pub fn value(&self) -> Option<f64> {
        match self.estimator {
            DxyEstimator::PerPairComparable => (self.counts.comparable > 0)
                .then(|| self.counts.diffs as f64 / self.counts.comparable as f64),
            DxyEstimator::AccessibleBp => {
                let has_comparable = self.counts.comparable > 0;
                (has_comparable && self.window_bp() > 0)
                    .then(|| self.sum_site_means / self.window_bp() as f64)
            }
        }
    }
}

/// Weir & Cockerham (1984) variance components, summed over sites and alleles.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct WcComponents {
    pub a: f64,
    pub b: f64,
    pub c: f64,
    /// Sites contributing (all groups called, `n_bar > 1`, `n_C > 0`).
    pub n_sites: u64,
    /// Sites skipped by the guards above.
    pub n_skipped: u64,
}

/// Nei (1986) gene diversity components, summed over sites.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NeiComponents {
    pub ht: f64,
    pub hs: f64,
    pub n_sites: u64,
    pub n_skipped: u64,
}

/// Between-group Fst (estimator-specific components; the integer `PairCounts`
/// fields are the pooled raw-count audit trail).
#[derive(Clone, Debug)]
pub struct FstResult {
    pub estimator: FstEstimator,
    pub contig: ContigId,
    pub start: u64,
    pub end: u64,
    pub group_a: u32,
    pub group_b: u32,
    pub within_a: PairCounts,
    pub within_b: PairCounts,
    pub between: PairCounts,
    /// Hudson: Σ over sites of `mpd_AB − (mpd_A + mpd_B)/2`.
    pub hudson_num: Option<f64>,
    /// Hudson: Σ over sites of `mpd_AB`.
    pub hudson_den: Option<f64>,
    pub wc: Option<WcComponents>,
    pub nei: Option<NeiComponents>,
    pub bias: Option<&'static str>,
}

impl FstResult {
    /// The estimate, or `None` when a denominator is empty or the between
    /// diversity is zero (Fst undefined without variation).
    pub fn value(&self) -> Option<f64> {
        match self.estimator {
            FstEstimator::Hudson => {
                let num = self.hudson_num?;
                let den = self.hudson_den?;
                (den != 0.0).then(|| num / den)
            }
            FstEstimator::WeirCockerham84 => {
                let wc = self.wc.as_ref()?;
                let total = wc.a + wc.b + wc.c;
                (total != 0.0).then(|| wc.a / total)
            }
            FstEstimator::Nei86 => {
                let nei = self.nei.as_ref()?;
                (nei.ht != 0.0).then(|| (nei.ht - nei.hs) / nei.ht)
            }
        }
    }
}
