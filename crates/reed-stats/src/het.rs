//! Individual-F (§9.8): per-sample additive vectors over genotype planes
//! (plink `--het` moment estimator).

// Casts in this file are bounded by construction: Sample counts and indices fit the accumulator dimensions.
#![allow(clippy::cast_possible_truncation)]

use reed_core::genotype::PlaneBlock;
use reed_core::math::sum::NeumaierSum;
use reed_core::sample::SampleMask;

/// Per-sample additive counters for plink's `--het` inbreeding
/// coefficient: F = (O_hom − E_hom) / (N − E_hom) with
/// E_hom = N − Σ_k 2·p_k(1−p_k) — plink's per-sample moment over
/// *biallelic* sites: O_hom = hom-ref + hom-alt full calls, N = full
/// diploid calls at biallelic sites (a half call is not a full call,
/// doc_v7 D7-14), E_hom = N − Σ 2p(1−p) where p is the alt frequency of
/// the site. All three are additive per sample.
#[derive(Clone, Debug, Default)]
pub struct IndividualFAcc {
    /// Per sample: homozygous (ref-ref or alt-alt) full-call count at
    /// biallelic sites.
    pub o_hom: Vec<u64>,
    /// Per sample: full diploid calls at biallelic sites.
    pub n_called: Vec<u64>,
    /// Per sample: Σ 2p(1−p) over the biallelic sites where the sample
    /// was fully called (E_hom = n_called − this sum).
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
            // Only full diploid calls enter O(HOM)/N: a half call (`./0`,
            // `./1`) is not a called genotype under the allele-level policy,
            // matching `count/mod.rs`'s n_half_called and plink2's
            // `--vcf-half-call m` (doc_v7 D7-14 — the missing slot's alt bit
            // is 0, so a half call used to be miscounted as hom/het).
            if s0m || s1m {
                continue;
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
