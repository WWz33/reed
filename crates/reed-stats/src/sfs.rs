//! Site frequency spectrum over complete-case sites
//! (`docs/architecture-spec.md` §9.4, the oracle's `CompleteCase` semantics;
//! the projected estimator is M7).

use reed_core::count::SiteGroupCounts;

/// Unfolded (or folded) SFS accumulator for one group. Complete-case: any
/// missing allele in the group skips the site (counted separately).
#[derive(Clone, Debug)]
pub struct SfsAcc {
    group: usize,
    n_chrom: u64,
    folded: bool,
    bins: Vec<u64>,
    n_sites: u64,
    n_skipped_missing: u64,
}

impl SfsAcc {
    /// Accumulator for `group`. Errors when `n_chrom < 2` (no spectrum).
    ///
    /// # Errors
    ///
    /// [`reed_core::CoreError::InvalidWindow`] when the group has fewer than
    /// two chromosomes.
    pub fn new(group: usize, n_group: u32, folded: bool) -> Result<Self, reed_core::CoreError> {
        let n_chrom = u64::from(n_group) * 2;
        if n_chrom < 2 {
            return Err(reed_core::CoreError::InvalidWindow);
        }
        let bins = if folded {
            (n_chrom / 2) as usize
        } else {
            n_chrom.saturating_sub(1) as usize
        };
        Ok(Self {
            group,
            n_chrom,
            folded,
            bins: vec![0; bins],
            n_sites: 0,
            n_skipped_missing: 0,
        })
    }

    /// Absorbs one site.
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        let c = &row[self.group];
        let n_group = self.n_chrom / 2;
        if c.n_missing_alleles(n_group as u32) != 0 {
            self.n_skipped_missing += 1;
            return;
        }
        self.n_sites += 1;
        let alt = u64::from(c.n_alt_alleles);
        if alt > 0 && alt < self.n_chrom {
            let index = if self.folded {
                alt.min(self.n_chrom - alt)
            } else {
                alt
            };
            self.bins[(index - 1) as usize] += 1;
        }
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        for (dst, src) in self.bins.iter_mut().zip(&other.bins) {
            *dst += src;
        }
        self.n_sites += other.n_sites;
        self.n_skipped_missing += other.n_skipped_missing;
    }

    /// Bin counts (index `i` = allele count `i+1`, or its minor count when
    /// folded).
    #[must_use]
    pub fn bins(&self) -> &[u64] {
        &self.bins
    }

    /// Contributing / skipped site tally.
    #[must_use]
    pub const fn site_tally(&self) -> (u64, u64) {
        (self.n_sites, self.n_skipped_missing)
    }

    /// Chromosome count the spectrum is defined over.
    #[must_use]
    pub const fn n_chrom(&self) -> u64 {
        self.n_chrom
    }
}
