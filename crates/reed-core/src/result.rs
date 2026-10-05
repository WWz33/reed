//! Result types: exact integer counts first, `f64` derived and auditable.

use crate::estimators::dxy::DxyEstimator;
use crate::estimators::fst::FstEstimator;
use crate::estimators::pi::PiEstimator;
use crate::window::ContigId;

/// Pair-based accounting at a set of sites (pixy-style).
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
    /// Which π estimator produced the result.
    pub estimator: PiEstimator,
    /// Contig the window is on.
    pub contig: ContigId,
    /// Window start (0-based, half-open).
    pub start: u64,
    /// Window end (exclusive).
    pub end: u64,
    /// Group the statistic is computed for.
    pub group: u32,
    /// Pooled raw pair-count audit trail.
    pub counts: PairCounts,
    /// Numerator of the `AccessibleBp` estimator: sum over sites of the
    /// per-site mean pairwise difference (sites with no comparable pair
    /// contribute zero).
    pub sum_site_means: f64,
    /// Segregating sites in the window.
    pub n_variant_sites: u64,
    /// Non-segregating sites in the window.
    pub n_invariant_sites: u64,
    /// Manifest caveat, if any.
    pub bias: Option<&'static str>,
}

impl PiResult {
    /// Window span in base pairs.
    pub fn window_bp(&self) -> u64 {
        self.end - self.start
    }

    /// The estimate, or `None` when the denominator is empty.
    //
    // Pair counts are bounded by the window's alleles (≪ 2^53); the ratio
    // casts are exact.
    #[allow(clippy::cast_precision_loss)]
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
    /// Which dxy estimator produced the result.
    pub estimator: DxyEstimator,
    /// Contig the window is on.
    pub contig: ContigId,
    /// Window start (0-based, half-open).
    pub start: u64,
    /// Window end (exclusive).
    pub end: u64,
    /// First group of the pair.
    pub group_a: u32,
    /// Second group of the pair.
    pub group_b: u32,
    /// Pooled raw pair-count audit trail.
    pub counts: PairCounts,
    /// Numerator of the `AccessibleBp` estimator: sum over sites of the
    /// per-site mean cross-pairwise difference (sites with no comparable
    /// pair contribute zero).
    pub sum_site_means: f64,
    /// Segregating sites in the window.
    pub n_variant_sites: u64,
    /// Non-segregating sites in the window.
    pub n_invariant_sites: u64,
    /// Manifest caveat, if any.
    pub bias: Option<&'static str>,
}

impl DxyResult {
    /// Window span in base pairs.
    pub fn window_bp(&self) -> u64 {
        self.end - self.start
    }

    /// The estimate, or `None` when the denominator is empty.
    //
    // Pair counts are bounded by the window's alleles (≪ 2^53); the ratio
    // casts are exact.
    #[allow(clippy::cast_precision_loss)]
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
    /// Between-population component summed over sites and alleles.
    pub a: f64,
    /// Within-population component summed over sites and alleles.
    pub b: f64,
    /// Inbreeding component summed over sites and alleles.
    pub c: f64,
    /// Sites contributing (all groups called, `n_bar > 1`, `n_C > 0`).
    pub n_sites: u64,
    /// Sites skipped by the guards above.
    pub n_skipped: u64,
}

/// Nei (1986) gene diversity components, summed over sites.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NeiComponents {
    /// Total gene diversity summed over sites.
    pub ht: f64,
    /// Within-group gene diversity (unweighted group mean) summed over sites.
    pub hs: f64,
    /// Sites contributing (both groups called).
    pub n_sites: u64,
    /// Sites skipped by the guards above.
    pub n_skipped: u64,
}

/// Between-group Fst (estimator-specific components; the integer `PairCounts`
/// fields are the pooled raw-count audit trail).
#[derive(Clone, Debug)]
pub struct FstResult {
    /// Which Fst estimator produced the result.
    pub estimator: FstEstimator,
    /// Contig the window is on.
    pub contig: ContigId,
    /// Window start (0-based, half-open).
    pub start: u64,
    /// Window end (exclusive).
    pub end: u64,
    /// First group of the pair.
    pub group_a: u32,
    /// Second group of the pair.
    pub group_b: u32,
    /// Within-group-a raw pair-count audit trail.
    pub within_a: PairCounts,
    /// Within-group-b raw pair-count audit trail.
    pub within_b: PairCounts,
    /// Between-group raw pair-count audit trail.
    pub between: PairCounts,
    /// Hudson: Σ over sites of `mpd_AB − (mpd_A + mpd_B)/2`.
    pub hudson_num: Option<f64>,
    /// Hudson: Σ over sites of `mpd_AB`.
    pub hudson_den: Option<f64>,
    /// Weir & Cockerham components (`WeirCockerham84` only).
    pub wc: Option<WcComponents>,
    /// Nei 1986 components (`Nei86` only).
    pub nei: Option<NeiComponents>,
    /// Manifest caveat, if any.
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
