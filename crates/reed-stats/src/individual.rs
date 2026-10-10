//! Individual-F and the joint (2D) SFS (§9.8–§9.9): per-sample additive
//! vectors over genotype planes (plink `--het` moment estimator) and the
//! pairwise joint spectrum with hypergeometric projection support.

// Casts in this file are bounded by construction: Sample counts and indices fit the accumulator dimensions.
#![allow(clippy::cast_possible_truncation)]

use reed_core::count::SiteGroupCounts;
use reed_core::genotype::PlaneBlock;
use reed_core::math::sum::NeumaierSum;
use reed_core::sample::SampleMask;

/// Per-sample additive counters for plink's `--het` inbreeding
/// coefficient: F = (O_hom − E_hom) / (N − E_hom) with
/// E_hom = N − Σ_k p_k(1−p_k)·… — plink's per-sample moment over
/// *biallelic* sites: O_hom = hom-ref + hom-alt calls, N = called
/// genotypes at biallelic sites, E_hom = N·(1 − 2p(1−p)) where p is the
/// alt frequency of the site. All three are additive per sample.
#[derive(Clone, Debug, Default)]
pub struct IndividualFAcc {
    /// Per sample: homozygous (ref-ref or alt-alt) genotype count at
    /// biallelic sites.
    pub o_hom: Vec<u64>,
    /// Per sample: called genotypes at biallelic sites.
    pub n_called: Vec<u64>,
    /// Per sample: Σ 2p(1−p) over the biallelic sites where the sample
    /// was called (E_hom = n_called − this sum).
    pub hetero_expected: Vec<NeumaierSum>,
    n_samples: usize,
}

impl IndividualFAcc {
    /// Vectors for `n_samples`.
    #[must_use]
    pub fn new(n_samples: usize) -> Self {
        Self {
            o_hom: vec![0; n_samples],
            n_called: vec![0; n_samples],
            hetero_expected: vec![NeumaierSum::new(); n_samples],
            n_samples,
        }
    }

    /// Absorbs one biallelic site's per-sample genotypes from the plane
    /// block. The loop is per selected sample (F6-11 in
    /// docs/rounds/v6/01-review.md: a word-parallel popcount rewrite is
    /// queued for the Y3 columnar-mask work).
    ///
    /// `p_alt` is the site's alt frequency (denominator: called alleles).
    /// Samples outside `mask` contribute nothing.
    pub fn push_site(&mut self, block: &PlaneBlock, site: usize, mask: &SampleMask, p_alt: f64) {
        let planes = block.site_planes(site);
        let w = planes.n_words();
        let het_exp = 2.0 * p_alt * (1.0 - p_alt);
        for i in 0..self.n_samples {
            let word = i / 64;
            let bit = 1_u64 << (i % 64);
            if word >= w || mask.words().get(word).is_none_or(|m| m & bit == 0) {
                continue;
            }
            let s0a = (planes.slot0_alt[word] & bit) != 0;
            let s0m = (planes.slot0_miss[word] & bit) != 0;
            let s1a = (planes.slot1_alt[word] & bit) != 0;
            let s1m = (planes.slot1_miss[word] & bit) != 0;
            if s0m && s1m {
                continue; // missing: not a called genotype
            }
            self.n_called[i] += 1;
            self.hetero_expected[i].add(het_exp);
            // Homozygous: both slots agree (0/0 or 1/1).
            if s0a == s1a {
                self.o_hom[i] += 1;
            }
        }
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        for (dst, src) in self.o_hom.iter_mut().zip(&other.o_hom) {
            *dst += src;
        }
        for (dst, src) in self.n_called.iter_mut().zip(&other.n_called) {
            *dst += src;
        }
        for (dst, src) in self.hetero_expected.iter_mut().zip(&other.hetero_expected) {
            dst.merge(src);
        }
    }

    /// F = (O_hom − E_hom) / (N − E_hom) per sample (plink `--het`).
    ///
    /// `None` when N == E_hom (no information) for that sample.
    #[must_use]
    pub fn f_values(&self) -> Vec<Option<f64>> {
        (0..self.n_samples)
            .map(|i| {
                let n = self.n_called[i] as f64;
                let e_hom = n - self.hetero_expected[i].value();
                let denom = n - e_hom;
                (n > 0.0 && denom != 0.0).then(|| (self.o_hom[i] as f64 - e_hom) / denom)
            })
            .collect()
    }
}

/// The joint (2D) site-frequency spectrum of two groups (§9.9, allel
/// `sf.js`/dadi semantics): bin (i, j) counts sites with i alt alleles in
/// group 1 and j in group 2 over complete-case sites. Folding projects
/// both axes to the minor side. Projection support arrives through the
/// same `reed_stats::sfs::projection_weights` product when needed; this
/// accumulator is the complete-case core (the M7 deliverable — the dadi
/// export is the writer's job).
#[derive(Clone, Debug)]
pub struct JointSfsAcc {
    g1: usize,
    g2: usize,
    n1: u32,
    n2: u32,
    folded: bool,
    /// `bins[i * (n2 + 1) + j]` for i in 0..=n1, j in 0..=n2 (Neumaier for
    /// projected support later; complete-case counts are u64 — the sums
    /// stay Neumaier so a projected variant can share the layout).
    bins: Vec<NeumaierSum>,
    n_sites: u64,
    n_skipped: u64,
}

impl JointSfsAcc {
    /// Accumulator over the `g1`/`g2` pair.
    #[must_use]
    pub fn new(g1: usize, g2: usize, n_chrom1: u32, n_chrom2: u32, folded: bool) -> Self {
        Self {
            g1,
            g2,
            n1: n_chrom1,
            n2: n_chrom2,
            folded,
            bins: vec![NeumaierSum::new(); (n_chrom1 as usize + 1) * (n_chrom2 as usize + 1)],
            n_sites: 0,
            n_skipped: 0,
        }
    }

    /// Absorbs one site (complete-case: any missing allele in either group
    /// skips the site).
    pub fn push(&mut self, row: &[SiteGroupCounts]) {
        let (c1, c2) = (&row[self.g1], &row[self.g2]);
        let n_group1 = self.n1 / 2;
        let n_group2 = self.n2 / 2;
        if c1.n_missing_alleles(n_group1) != 0 || c2.n_missing_alleles(n_group2) != 0 {
            self.n_skipped += 1;
            return;
        }
        self.n_sites += 1;
        let (mut i, mut j) = (u64::from(c1.n_alt_alleles), u64::from(c2.n_alt_alleles));
        if self.folded {
            i = i.min(self.n1 as u64 - i);
            j = j.min(self.n2 as u64 - j);
        }
        let idx = i as usize * (self.n2 as usize + 1) + j as usize;
        self.bins[idx].add(1.0);
    }

    /// Component-wise merge.
    pub fn merge(&mut self, other: &Self) {
        debug_assert_eq!(self.bins.len(), other.bins.len(), "invariant: same shape");
        for (dst, src) in self.bins.iter_mut().zip(&other.bins) {
            dst.merge(src);
        }
        self.n_sites += other.n_sites;
        self.n_skipped += other.n_skipped;
    }

    /// Bin (i, j).
    #[must_use]
    pub fn at(&self, i: u32, j: u32) -> f64 {
        self.bins[i as usize * (self.n2 as usize + 1) + j as usize].value()
    }

    /// Projected / skipped site tally.
    #[must_use]
    pub const fn site_tally(&self) -> (u64, u64) {
        (self.n_sites, self.n_skipped)
    }
}
