//! Tajima's D in the reference missing-data form (Tajima 1989; the
//! `ThetaRecov` denominator evaluated once at the window-average allele count,
//! Korunes & Samuk 2021 `calc_tajima_d`).
//!
//! `D = (raw_pi − raw_theta) / sqrt(e1·S + e2·S(S−1))` where `raw_pi` is the
//! sum of per-site unbiased π, `raw_theta` the per-allele-class Watterson
//! numerator, `S = Σ(k_i − 1)` and the denominator is evaluated once at
//! `n = round(Σn_i / sites)`.

// Casts in this file are bounded by construction: The mean sample size is bounded by the sample count.
#![allow(clippy::cast_possible_truncation)]

use reed_core::count::SiteGroupCounts;
use reed_core::math::harmonic::{a1, a2};
use reed_core::math::sum::NeumaierSum;

use crate::theta::ThetaAcc;

/// Tajima's D accumulator for one group.
#[derive(Clone, Debug)]
pub struct TajimaAcc {
    group: usize,
    n_chrom: u64,
    /// Per-allele-class mutation/site counts (the θ arm).
    theta: ThetaAcc,
    /// Σ per-site unbiased π over sites with ≥ 2 observed alleles.
    raw_pi: NeumaierSum,
}

impl TajimaAcc {
    /// Accumulator for `group`; `None` when the group has < 2 chromosomes.
    #[must_use]
    pub fn new(group: usize, n_group: u32) -> Option<Self> {
        let theta = ThetaAcc::new(group, n_group)?;
        let n_chrom = u64::from(n_group) * 2;
        Some(Self {
            group,
            n_chrom,
            theta,
            raw_pi: NeumaierSum::new(),
        })
    }

    /// Absorbs one site.
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        let c = &row[self.group];
        let n = u64::from(c.n_ref_alleles) + u64::from(c.n_alt_alleles);
        if n >= 2 {
            let comps = n * (n - 1) / 2;
            let diffs = u64::from(c.n_ref_alleles) * u64::from(c.n_alt_alleles);
            self.raw_pi.add(diffs as f64 / comps as f64);
        }
        self.theta.push(row);
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        self.theta.merge(&other.theta);
        self.raw_pi.merge(&other.raw_pi);
    }

    /// The window's raw π sum (`raw_pi`).
    #[must_use]
    pub fn raw_pi(&self) -> f64 {
        self.raw_pi.value()
    }

    /// The window's raw Watterson θ (`raw_watterson_theta`).
    #[must_use]
    pub fn raw_theta(&self) -> f64 {
        self.theta.raw_theta()
    }

    /// Tajima's D denominator: evaluated once at the rounded window-average
    /// observed allele count (Korunes & Samuk 2021); `None` where
    /// undefined (no sites, no mutations, mean n < 2).
    #[must_use]
    pub fn d_stdev(&self) -> Option<f64> {
        let sites = self.theta.n_sites();
        let mutations = self.theta.n_mutations();
        if sites == 0 || mutations == 0 {
            return None;
        }
        // np.rint semantics: round half to even, matching the reference
        // `int(np.rint(total_allele_count / num_sites))`.
        let mean_n = (self.theta.total_allele_count() as f64 / sites as f64).round_ties_even();
        if mean_n < 2.0 {
            return None;
        }
        let n = usize::try_from(mean_n as u64).unwrap_or(usize::MAX);
        let (e1, e2) = tajima_constants(n)?;
        let s = mutations as f64;
        let variance = e1 * s + e2 * s * (s - 1.0);
        (variance > 0.0).then(|| variance.sqrt())
    }

    /// Window-level D; `None` when `S < min_mutations` or the denominator is
    /// undefined.
    #[must_use]
    pub fn value(&self, min_mutations: u64) -> Option<f64> {
        if self.theta.n_mutations() < min_mutations {
            return None;
        }
        let stdev = self.d_stdev()?;
        (stdev > 0.0).then(|| (self.raw_pi() - self.raw_theta()) / stdev)
    }

    /// Window-total π sum (sum of per-site means).
    #[must_use]
    pub const fn pi_sum(&self) -> &NeumaierSum {
        &self.raw_pi
    }

    /// Contributing site tally.
    #[must_use]
    pub const fn site_tally(&self) -> (u64, u64) {
        self.theta.site_tally()
    }

    /// Segregating sites (any ALT observed).
    #[must_use]
    pub const fn segregating(&self) -> u64 {
        self.theta.segregating()
    }

    /// Chromosome count.
    #[must_use]
    pub const fn n_chrom(&self) -> u64 {
        self.n_chrom
    }
}

/// Tajima 1989's `e1`, `e2` at `n` sampled alleles (`None` for n < 2).
fn tajima_constants(n: usize) -> Option<(f64, f64)> {
    if n < 2 {
        return None;
    }
    let a1n = a1(n);
    let a2n = a2(n);
    if a1n <= 0.0 {
        return None;
    }
    let nf = n as f64;
    let b1 = (nf + 1.0) / (3.0 * (nf - 1.0));
    let b2 = 2.0 * (nf * nf + nf + 3.0) / (9.0 * nf * (nf - 1.0));
    let c1 = b1 - 1.0 / a1n;
    let c2 = b2 - (nf + 2.0) / (a1n * nf) + a2n / (a1n * a1n);
    let e1 = c1 / a1n;
    let e2 = c2 / (a1n * a1n + a2n);
    Some((e1, e2))
}
