//! Watterson's θ over complete-case sites (`docs/architecture-spec.md` §9.4,
//! oracle semantics: `S / a₁(n)`, per-bp variant divides by window bp).

use reed_core::count::SiteGroupCounts;
use reed_core::math::harmonic::a1;

/// θ accumulator for one group.
#[derive(Clone, Debug)]
pub struct ThetaAcc {
    group: usize,
    n_chrom: u64,
    segregating: u64,
    n_sites: u64,
    n_skipped_missing: u64,
}

impl ThetaAcc {
    /// Accumulator for `group`; `None` when the group has < 2 chromosomes.
    #[must_use]
    pub fn new(group: usize, n_group: u32) -> Option<Self> {
        let n_chrom = u64::from(n_group) * 2;
        if n_chrom >= 2 {
            Some(Self {
                group,
                n_chrom,
                segregating: 0,
                n_sites: 0,
                n_skipped_missing: 0,
            })
        } else {
            None
        }
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
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        self.segregating += other.segregating;
        self.n_sites += other.n_sites;
        self.n_skipped_missing += other.n_skipped_missing;
    }

    /// Window-total θ_W = `S / a₁(n)`; `None` without complete-case sites.
    #[must_use]
    pub fn value_total(&self) -> Option<f64> {
        (self.n_sites > 0).then(|| {
            let n = usize::try_from(self.n_chrom).unwrap_or(usize::MAX);
            self.segregating as f64 / a1(n)
        })
    }

    /// θ_W per bp (scikit-allel windowed `value()`).
    #[must_use]
    pub fn value_per_bp(&self, window_bp: u64) -> Option<f64> {
        (self.n_sites > 0 && window_bp > 0).then(|| {
            let n = usize::try_from(self.n_chrom).unwrap_or(usize::MAX);
            (self.segregating as f64 / a1(n)) / window_bp as f64
        })
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
