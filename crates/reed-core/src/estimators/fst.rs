//! Fst estimators (sealed definitions).
//!
//! - [`FstEstimator::Hudson`] — Hudson 1992 / Bhatia et al. 2013, exactly as
//!   scikit-allel `hudson_fst`: per site, `num_s = mpd_AB − (mpd_A + mpd_B)/2`
//!   and `den_s = mpd_AB` (mean pairwise differences, site skipped when either
//!   group has no called allele or no comparable pair); `Fst = Σnum/Σden`.
//!   The integer `PairCounts` fields are retained as an audit trail of the
//!   pooled raw counts (not the estimator's arithmetic).
//! - [`FstEstimator::WeirCockerham84`] — Weir & Cockerham 1984 variance
//!   components `a/(a+b+c)` ported from scikit-allel (MIT): per-site components
//!   summed over sites and alleles; guards: both groups called, `n_bar > 1`,
//!   `n_C > 0`.
//! - [`FstEstimator::Nei86`] — Nei 1986 GST: `(HT − HS)/HT` with the unweighted
//!   mean over groups, summed over sites where both groups are called.
//!
//! All estimators use the same site set for numerator and denominator, so
//! variant-only input is not structurally broken (unlike π/dxy whose
//! denominator is a site census); the manifest gate is therefore a bias note,
//! not a refusal.

use crate::manifest::{Availability, SiteCompleteness};
use crate::population::PopulationView;
use crate::result::{FstResult, NeiComponents, PairCounts, WcComponents};
use crate::window::WindowBuffer;

use super::Unavailable;

/// Sealed Fst estimator definitions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FstEstimator {
    /// Hudson (1992), elaborated by Bhatia et al. (2013) — the field standard,
    /// same estimator as scikit-allel's `hudson_fst`.
    Hudson,
    /// Weir & Cockerham (1984), components ported from scikit-allel (MIT).
    WeirCockerham84,
    /// Nei (1986) GST.
    Nei86,
}

impl FstEstimator {
    pub fn id(self) -> &'static str {
        match self {
            FstEstimator::Hudson => "fst.hudson",
            FstEstimator::WeirCockerham84 => "fst.weir_cockerham_84",
            FstEstimator::Nei86 => "fst.nei_86",
        }
    }

    pub fn citation(self) -> &'static str {
        match self {
            FstEstimator::Hudson => {
                "Hudson 1992; Bhatia et al. 2013; per-site normalisation as scikit-allel \
                 hudson_fst (mean pairwise differences)"
            }
            FstEstimator::WeirCockerham84 => {
                "Weir & Cockerham 1984 (a/(a+b+c), summed over sites and alleles)"
            }
            FstEstimator::Nei86 => "Nei 1986 (GST, unweighted mean over groups)",
        }
    }

    pub fn availability(self, manifest: SiteCompleteness) -> Availability {
        match manifest {
            SiteCompleteness::WithInvariants => Availability::Available,
            SiteCompleteness::VariantOnly(_) => Availability::Biased(
                "variant-only input: Fst is computed on variant sites only; ascertainment \
                 bias is possible",
            ),
        }
    }
}

/// Computes Fst between two groups.
pub fn fst(
    window: &WindowBuffer,
    pops: &PopulationView,
    group_a: u32,
    group_b: u32,
    estimator: FstEstimator,
) -> Result<FstResult, Unavailable> {
    if group_a == group_b {
        return Err(Unavailable {
            reason: "Fst requires two distinct groups",
        });
    }
    let bias = match estimator.availability(window.manifest) {
        Availability::Available => None,
        Availability::Biased(reason) => Some(reason),
        Availability::Unavailable(reason) => return Err(Unavailable { reason }),
    };

    let samples_a = pops.samples(group_a);
    let samples_b = pops.samples(group_b);
    let n_chrom_a = 2 * samples_a.len() as u64;
    let n_chrom_b = 2 * samples_b.len() as u64;
    let within_total_a = n_chrom_a * n_chrom_a.saturating_sub(1) / 2;
    let within_total_b = n_chrom_b * n_chrom_b.saturating_sub(1) / 2;
    let between_total = n_chrom_a * n_chrom_b;

    let mut within_a = PairCounts::default();
    let mut within_b = PairCounts::default();
    let mut between = PairCounts::default();
    let mut wc = WcComponents::default();
    let mut nei = NeiComponents::default();
    // Hudson: per-site normalised sums (allel semantics).
    let (mut hudson_num, mut hudson_den) = (0.0f64, 0.0f64);

    for site in 0..window.n_sites() {
        let a = window.allele_counts(site, samples_a);
        let b = window.allele_counts(site, samples_b);

        let n_a = a.ref_count + a.alt_count;
        let n_b = b.ref_count + b.alt_count;

        // Audit trail: pooled integer counts (missing pairs included).
        let comps_a = n_a * n_a.saturating_sub(1) / 2;
        let comps_b = n_b * n_b.saturating_sub(1) / 2;
        within_a.comparable += comps_a;
        within_a.diffs += a.ref_count * a.alt_count;
        within_a.missing += within_total_a - comps_a;
        within_b.comparable += comps_b;
        within_b.diffs += b.ref_count * b.alt_count;
        within_b.missing += within_total_b - comps_b;
        between.comparable += n_a * n_b;
        between.diffs += a.ref_count * b.alt_count + a.alt_count * b.ref_count;
        between.missing += between_total - n_a * n_b;

        if estimator == FstEstimator::Hudson {
            // allel inclusion rule: `den` (between) accumulates whenever both
            // groups have ≥1 called allele; `num` also needs both within terms
            // (n ≥ 2 per group), otherwise it is NaN and allel's `nansum` skips
            // it while keeping the site in the denominator. (pixy instead drops
            // the site from both — we follow the cited allel semantics.)
            let mpd_a = (n_a >= 2).then(|| {
                a.ref_count as f64 * a.alt_count as f64 / (n_a * (n_a - 1) / 2) as f64
            });
            let mpd_b = (n_b >= 2).then(|| {
                b.ref_count as f64 * b.alt_count as f64 / (n_b * (n_b - 1) / 2) as f64
            });
            let mpd_ab = (n_a >= 1 && n_b >= 1)
                .then(|| (a.ref_count * b.alt_count + a.alt_count * b.ref_count) as f64
                    / (n_a * n_b) as f64);
            if let Some(pab) = mpd_ab {
                hudson_den += pab;
                if let (Some(pa), Some(pb)) = (mpd_a, mpd_b) {
                    hudson_num += pab - (pa + pb) / 2.0;
                }
            }
        }

        match estimator {
            FstEstimator::Hudson => {}
            FstEstimator::WeirCockerham84 => {
                if n_a == 0 || n_b == 0 {
                    wc.n_skipped += 1;
                    continue;
                }
                // allel semantics (fst.py:163): n = an // 2 — integer floor,
                // losing a lone called allele of a half-missing genotype.
                let (n1, n2) = ((n_a / 2) as f64, (n_b / 2) as f64);
                if n1 == 0.0 || n2 == 0.0 {
                    wc.n_skipped += 1;
                    continue;
                }
                let n_total = n1 + n2;
                let n_bar = n_total / 2.0;
                // r = 2 populations: n_C = (n_total − Σ n_i²/n_total) / (r − 1).
                let n_c = (n_total - (n1 * n1 + n2 * n2) / n_total) / (2.0 - 1.0);
                if n_bar <= 1.0 || n_c <= 0.0 {
                    wc.n_skipped += 1;
                    continue;
                }
                let p_ref_1 = a.ref_count as f64 / n_a as f64;
                let p_ref_2 = b.ref_count as f64 / n_b as f64;
                let p_bar = (a.ref_count + b.ref_count) as f64 / (n_a + n_b) as f64;
                // s² = Σ nᵢ(pᵢ−p̄)² / (n_bar·(r−1)), r = 2.
                let s2 = (n1 * (p_ref_1 - p_bar).powi(2) + n2 * (p_ref_2 - p_bar).powi(2))
                    / (n_bar * (2.0 - 1.0));
                let het_1 = window.genotype_counts(site, samples_a).het as f64;
                let het_2 = window.genotype_counts(site, samples_b).het as f64;
                let h_bar = (het_1 + het_2) / n_total;
                let term = p_bar * (1.0 - p_bar) - (2.0 - 1.0) * s2 / 2.0 - h_bar / 4.0;
                let a_allele = (n_bar / n_c) * (s2 - term / (n_bar - 1.0));
                let b_allele = (n_bar / (n_bar - 1.0))
                    * (p_bar * (1.0 - p_bar) - (2.0 - 1.0) * s2 / 2.0 - (2.0 * n_bar - 1.0)
                        * h_bar
                        / (4.0 * n_bar));
                let c_allele = h_bar / 2.0;
                // Biallelic: the alt-allele components are identical (p̄ ↦ 1−p̄).
                wc.a += 2.0 * a_allele;
                wc.b += 2.0 * b_allele;
                wc.c += 2.0 * c_allele;
                wc.n_sites += 1;
            }
            FstEstimator::Nei86 => {
                if n_a == 0 || n_b == 0 {
                    nei.n_skipped += 1;
                    continue;
                }
                let p1 = a.alt_count as f64 / n_a as f64;
                let p2 = b.alt_count as f64 / n_b as f64;
                let hs = (2.0 * p1 * (1.0 - p1) + 2.0 * p2 * (1.0 - p2)) / 2.0;
                let p_bar = (p1 + p2) / 2.0;
                let ht = 2.0 * p_bar * (1.0 - p_bar);
                nei.hs += hs;
                nei.ht += ht;
                nei.n_sites += 1;
            }
        }
    }

    Ok(FstResult {
        estimator,
        contig: window.contig,
        start: window.start,
        end: window.end,
        group_a,
        group_b,
        within_a,
        within_b,
        between,
        hudson_num: (estimator == FstEstimator::Hudson).then_some(hudson_num),
        hudson_den: (estimator == FstEstimator::Hudson).then_some(hudson_den),
        wc: (estimator == FstEstimator::WeirCockerham84).then_some(wc),
        nei: (estimator == FstEstimator::Nei86).then_some(nei),
        bias,
    })
}
