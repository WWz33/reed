//! Watterson's θ in the reference missing-data form (Watterson 1975;
//! mutation-count per observed-allele class, Korunes & Samuk 2021,
//! `calc_watterson_theta`).
//!
//! Sites are grouped by their observed allele count `n_i`; each variant
//! site contributes `k_i − 1` mutations (k = distinct observed alleles, so
//! 1 for a biallelic segregating site) and each class contributes
//! `mutations(n) / a(n)` with `a(n) = Σ_{j=1}^{n−1} 1/j`.

use reed_core::count::SiteGroupCounts;
use reed_core::math::harmonic::a1;

/// θ accumulator for one group.
#[derive(Clone, Debug)]
pub struct ThetaAcc {
    group: usize,
    n_chrom: u64,
    /// `mutations[n]` = Σ (k_i − 1) over variant sites with `n` observed alleles.
    mutations: Vec<u64>,
    /// `sites[n]` = number of sites with `n` observed alleles.
    sites: Vec<u64>,
    /// Σ n_i over sites with n_i > 0 (for the window-average n).
    total_allele_count: u64,
    /// Sites with at least one observed allele.
    n_sites: u64,
    /// Variant sites (any ALT observed).
    n_var_sites: u64,
}

impl ThetaAcc {
    /// Accumulator for `group`; `None` when the group has < 2 chromosomes.
    #[must_use]
    pub fn new(group: usize, n_group: u32) -> Option<Self> {
        let n_chrom = u64::from(n_group) * 2;
        if n_chrom < 2 {
            return None;
        }
        let width = usize::try_from(n_chrom)
            .unwrap_or(usize::MAX)
            .saturating_add(1);
        Some(Self {
            group,
            n_chrom,
            mutations: vec![0; width],
            sites: vec![0; width],
            total_allele_count: 0,
            n_sites: 0,
            n_var_sites: 0,
        })
    }

    /// Absorbs one site (all observed-allele counts; missing data included).
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        let c = &row[self.group];
        let n = u64::from(c.n_ref_alleles) + u64::from(c.n_alt_alleles);
        if n == 0 {
            return;
        }
        let n_idx = usize::try_from(n).unwrap_or(usize::MAX);
        if n_idx >= self.sites.len() {
            self.sites.resize(n_idx + 1, 0);
            self.mutations.resize(n_idx + 1, 0);
        }
        self.sites[n_idx] += 1;
        self.total_allele_count += n;
        self.n_sites += 1;
        let alt = u64::from(c.n_alt_alleles);
        if alt > 0 {
            self.n_var_sites += 1;
        }
        // the reference counts mutations as k − 1 (k = distinct observed alleles), so a
        // site fixed for ALT has k = 1 and contributes 0 — it is not
        // segregating. Biallelic sites with both alleles observed contribute 1.
        // Multi-allelic sites need per-allele counts (04-algorithms.md §3),
        // which the count layer does not carry yet; the driver excludes such
        // records for now.
        if alt > 0 && alt < n {
            self.mutations[n_idx] += 1;
        }
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        if other.sites.len() > self.sites.len() {
            self.sites.resize(other.sites.len(), 0);
            self.mutations.resize(other.mutations.len(), 0);
        }
        for (dst, src) in self.sites.iter_mut().zip(&other.sites) {
            *dst += src;
        }
        for (dst, src) in self.mutations.iter_mut().zip(&other.mutations) {
            *dst += src;
        }
        self.total_allele_count += other.total_allele_count;
        self.n_sites += other.n_sites;
        self.n_var_sites += other.n_var_sites;
    }

    /// `Σ_n mutations(n) / a(n)` — the window's raw Watterson θ.
    #[must_use]
    pub fn raw_theta(&self) -> f64 {
        let mut raw = 0.0;
        for (n, &mutations) in self.mutations.iter().enumerate() {
            if mutations == 0 || n < 2 {
                continue;
            }
            let a = a1(n);
            if a > 0.0 {
                raw += mutations as f64 / a;
            }
        }
        raw
    }

    /// Window-average θ = `raw / n_sites` (`avg_watterson_theta`).
    #[must_use]
    pub fn value_total(&self) -> Option<f64> {
        (self.n_sites > 0).then(|| self.raw_theta() / self.n_sites as f64)
    }

    /// `--variant-only` form: `raw / L` (04-algorithms.md §3).
    #[must_use]
    pub fn value_variant_only(&self, window_bp: u64) -> Option<f64> {
        (window_bp > 0).then(|| self.raw_theta() / window_bp as f64)
    }

    /// θ per bp (legacy wide-TSV convention).
    #[must_use]
    pub fn value_per_bp(&self, window_bp: u64) -> Option<f64> {
        (self.n_sites > 0 && window_bp > 0)
            .then(|| (self.raw_theta() / self.n_sites as f64) / window_bp as f64)
    }

    /// The reference `weighted_no_sites`: `Σ_n sites(n) · n / max(n)`.
    #[must_use]
    pub fn weighted_no_sites(&self) -> f64 {
        let Some(max_n) = self.sites.iter().rposition(|&count| count > 0) else {
            return 0.0;
        };
        if max_n == 0 {
            return 0.0;
        }
        let mut weighted = 0.0;
        for (n, &count) in self.sites.iter().enumerate() {
            if count > 0 {
                weighted += count as f64 * n as f64 / max_n as f64;
            }
        }
        weighted
    }

    /// Sites with ≥ 1 observed allele.
    #[must_use]
    pub const fn n_sites(&self) -> u64 {
        self.n_sites
    }

    /// Variant sites.
    #[must_use]
    pub const fn n_var_sites(&self) -> u64 {
        self.n_var_sites
    }

    /// Σ observed alleles over sites with calls (Tajima's mean-n input).
    #[must_use]
    pub const fn total_allele_count(&self) -> u64 {
        self.total_allele_count
    }

    /// `Σ (k_i − 1)` over variant sites — Tajima's `S`.
    #[must_use]
    pub fn n_mutations(&self) -> u64 {
        self.mutations.iter().sum()
    }

    /// Contributing / skipped site tally (skipped is always 0 under the reference
    /// semantics: missing data is incorporated, not dropped).
    #[must_use]
    pub const fn site_tally(&self) -> (u64, u64) {
        (self.n_sites, 0)
    }

    /// Segregating-site count.
    #[must_use]
    pub const fn segregating(&self) -> u64 {
        self.n_var_sites
    }

    /// Chromosome count.
    #[must_use]
    pub const fn n_chrom(&self) -> u64 {
        self.n_chrom
    }
}
