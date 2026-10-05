//! Tajima's D over complete-case sites (`docs/architecture-spec.md` §9.4,
//! Tajima 1989 moment variance; `min_sites` guard as in scikit-allel).

use reed_core::count::SiteGroupCounts;
use reed_core::math::harmonic::{a1, a2};
use reed_core::math::sum::NeumaierSum;

/// Tajima's D accumulator for one group.
#[derive(Clone, Debug)]
pub struct TajimaAcc {
    group: usize,
    n_chrom: u64,
    segregating: u64,
    pi: NeumaierSum,
    n_sites: u64,
    n_skipped_missing: u64,
}

impl TajimaAcc {
    /// Accumulator for `group`; `None` when the group has < 2 chromosomes.
    #[must_use]
    pub fn new(group: usize, n_group: u32) -> Option<Self> {
        let n_chrom = u64::from(n_group) * 2;
        (n_chrom >= 2).then(|| Self {
            group,
            n_chrom,
            segregating: 0,
            pi: NeumaierSum::new(),
            n_sites: 0,
            n_skipped_missing: 0,
        })
    }

    /// Absorbs one site.
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        let c = &row[self.group];
        let n_group = (self.n_chrom / 2) as u32;
        if c.n_missing_alleles(n_group) != 0 {
            self.n_skipped_missing += 1;
            return;
        }
        self.n_sites += 1;
        let alt = u64::from(c.n_alt_alleles);
        if alt > 0 && alt < self.n_chrom {
            self.segregating += 1;
        }
        let ref_alleles = u64::from(c.n_ref_alleles);
        let pairs = self.n_chrom * (self.n_chrom - 1) / 2;
        self.pi.add((ref_alleles * alt) as f64 / pairs as f64);
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        self.segregating += other.segregating;
        self.pi.merge(&other.pi);
        self.n_sites += other.n_sites;
        self.n_skipped_missing += other.n_skipped_missing;
    }

    /// Window-level D; `None` when `S < min_sites` or the variance is zero.
    #[must_use]
    pub fn value(&self, min_sites: u64) -> Option<f64> {
        let n = self.n_chrom;
        if n < 2 {
            return None;
        }
        let n = usize::try_from(n).unwrap_or(usize::MAX);
        let a1n = a1(n);
        let a2n = a2(n);
        let theta_w = if a1n > 0.0 { self.segregating as f64 / a1n } else { 0.0 };
        let nf = n as f64;
        let b1 = (nf + 1.0) / (3.0 * (nf - 1.0));
        let b2 = 2.0 * (nf * nf + nf + 3.0) / (9.0 * nf * (nf - 1.0));
        let c1 = b1 - 1.0 / a1n;
        let c2 = b2 - (nf + 2.0) / (a1n * nf) + a2n / (a1n * a1n);
        let e1 = c1 / a1n;
        let e2 = c2 / (a1n * a1n + a2n);
        let s = self.segregating as f64;
        let variance = e1 * s + e2 * s * (s - 1.0);
        (self.segregating >= min_sites && variance > 0.0)
            .then(|| (self.pi.value() - theta_w) / variance.sqrt())
    }

    /// Window-total π (sum over complete-case sites of the per-site mean).
    #[must_use]
    pub const fn pi_sum(&self) -> &NeumaierSum {
        &self.pi
    }

    /// Contributing / skipped site tally.
    #[must_use]
    pub const fn site_tally(&self) -> (u64, u64) {
        (self.n_sites, self.n_skipped_missing)
    }

    /// Segregating-site count.
    #[must_use]
    pub const fn segregating(&self) -> u64 {
        self.segregating
    }

    /// Chromosome count.
    #[must_use]
    pub const fn n_chrom(&self) -> u64 {
        self.n_chrom
    }
}
