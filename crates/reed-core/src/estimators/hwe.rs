//! Hardy-Weinberg equilibrium and heterozygosity family.
//!
//! These are **site-local** statistics: every returned value depends only on
//! that site's genotype counts, so results are exact wherever the site is
//! present in the input. There is no cross-site denominator, hence no
//! manifest gate (§4.2 位点内统计豁免) — sites absent from the input are
//! simply not reported.

use crate::estimators::Unavailable;
use crate::population::PopulationView;
use crate::window::{ContigId, WindowBuffer};

use super::samples_of;

/// Exact test for Hardy-Weinberg equilibrium (Wigginton, Cutler & Abecasis
/// 2005), ported from sgkit (Apache-2.0).
///
/// `obs_het` heterozygotes, `obs_hom1` / `obs_hom2` the two homozygote counts
/// (order irrelevant). Returns the two-sided p-value (sum of the probabilities
/// of all genotype configurations no more likely than the observed one), or
/// NaN for an empty sample.
//
// Faithful port of the Wigginton/Cutler/Abecasis reference implementation:
// i64 intermediates with `as` casts, exactly as published. Counts are bounded
// by the site's sample count (≪ 2^53 and non-negative where indexed), so no
// cast loses information in practice.
#[allow(
    clippy::cast_possible_wrap,
    clippy::cast_sign_loss,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss
)]
pub fn hardy_weinberg_p_value(obs_het: u64, obs_hom1: u64, obs_hom2: u64) -> f64 {
    let hom_major = obs_hom1.max(obs_hom2) as i64;
    let hom_minor = obs_hom1.min(obs_hom2) as i64;
    let obs_mac = 2 * hom_minor + obs_het as i64;
    let obs_n = obs_het as i64 + hom_major + hom_minor;
    if obs_n == 0 {
        return f64::NAN;
    }

    let mut probs = vec![0.0f64; obs_mac as usize + 1];
    let mut mid = ((obs_mac * (2 * obs_n - obs_mac)) as f64 / (2 * obs_n) as f64) as i64;
    if (obs_mac & 1) ^ (mid & 1) != 0 {
        mid += 1;
    }
    probs[mid as usize] = 1.0;
    let mut prob_sum = 1.0f64;

    let mut curr_hets = mid;
    let mut curr_hom_minor = (obs_mac - mid) / 2;
    let mut curr_hom_major = obs_n - curr_hets - curr_hom_minor;
    while curr_hets > 1 {
        let next = probs[curr_hets as usize] * curr_hets as f64 * (curr_hets as f64 - 1.0)
            / (4.0 * (curr_hom_minor as f64 + 1.0) * (curr_hom_major as f64 + 1.0));
        probs[(curr_hets - 2) as usize] = next;
        prob_sum += next;
        curr_hom_minor += 1;
        curr_hom_major += 1;
        curr_hets -= 2;
    }

    let mut curr_hets = mid;
    let mut curr_hom_minor = (obs_mac - mid) / 2;
    let mut curr_hom_major = obs_n - curr_hets - curr_hom_minor;
    while curr_hets <= obs_mac - 2 {
        let next = probs[curr_hets as usize] * 4.0 * curr_hom_minor as f64
            * curr_hom_major as f64
            / ((curr_hets as f64 + 2.0) * (curr_hets as f64 + 1.0));
        probs[(curr_hets + 2) as usize] = next;
        prob_sum += next;
        curr_hom_minor -= 1;
        curr_hom_major -= 1;
        curr_hets += 2;
    }

    if prob_sum <= 0.0 {
        return f64::NAN;
    }
    for p in &mut probs {
        *p /= prob_sum;
    }
    let obs_prob = probs[obs_het as usize];
    let p: f64 = probs.iter().filter(|&&x| x <= obs_prob).sum();
    p.clamp(0.0, 1.0)
}

/// Per-site exact HWE test result for one group.
#[derive(Clone, Copy, Debug)]
pub struct HweSiteResult {
    /// Contig the site is on.
    pub contig: ContigId,
    /// Group the genotype counts come from.
    pub group: u32,
    /// 0-based position of the site.
    pub pos: u32,
    /// Individuals homozygous for the REF allele.
    pub hom_ref: u64,
    /// Heterozygous individuals.
    pub het: u64,
    /// Individuals homozygous for the ALT allele.
    pub hom_alt: u64,
    /// Two-sided exact-test p-value.
    pub p_value: f64,
}

/// Per-site exact HWE test over a group; individuals with any missing allele
/// are excluded, sites with no called individual are skipped.
pub fn hwe_exact(
    window: &WindowBuffer,
    pops: &PopulationView,
    group: u32,
) -> Result<Vec<HweSiteResult>, Unavailable> {
    let samples = samples_of(pops, window, group)?;
    let mut out = Vec::new();
    for site in 0..window.n_sites() {
        let counts = window.genotype_counts(site, samples);
        if counts.n_called == 0 {
            continue;
        }
        out.push(HweSiteResult {
            contig: window.contig,
            group,
            pos: window.pos[site],
            hom_ref: counts.hom_ref,
            het: counts.het,
            hom_alt: counts.hom_alt,
            p_value: hardy_weinberg_p_value(counts.het, counts.hom_ref, counts.hom_alt),
        });
    }
    Ok(out)
}

/// Per-site observed/expected heterozygosity and Fis for one group.
#[derive(Clone, Copy, Debug)]
pub struct HetSiteResult {
    /// Contig the site is on.
    pub contig: ContigId,
    /// Group the counts come from.
    pub group: u32,
    /// 0-based position of the site.
    pub pos: u32,
    /// Individuals with both alleles called.
    pub n_called: u64,
    /// Called heterozygous individuals.
    pub n_het: u64,
    /// Observed heterozygosity, `n_het / n_called`.
    pub hobs: f64,
    /// `2p(1-p)` with p the ALT frequency among **called alleles**
    /// (per-allele deletion, matching `allele_counts`).
    pub hexp: f64,
    /// `1 - hobs/hexp`; `None` when `hexp == 0` (monomorphic).
    pub fis: Option<f64>,
}

/// Per-site heterozygosity over a group; sites with no called individual are
/// skipped.
//
// Counts are bounded by the site's called alleles (≪ 2^53), so the ratio
// casts cannot lose precision in practice.
#[allow(clippy::cast_precision_loss)]
pub fn heterozygosity(
    window: &WindowBuffer,
    pops: &PopulationView,
    group: u32,
) -> Result<Vec<HetSiteResult>, Unavailable> {
    let samples = samples_of(pops, window, group)?;
    let mut out = Vec::new();
    for site in 0..window.n_sites() {
        let counts = window.genotype_counts(site, samples);
        if counts.n_called == 0 {
            continue;
        }
        let alleles = window.allele_counts(site, samples);
        let n_alleles_called = alleles.ref_count + alleles.alt_count;
        let p = alleles.alt_count as f64 / n_alleles_called as f64;
        let hobs = counts.het as f64 / counts.n_called as f64;
        let hexp = 2.0 * p * (1.0 - p);
        out.push(HetSiteResult {
            contig: window.contig,
            group,
            pos: window.pos[site],
            n_called: counts.n_called,
            n_het: counts.het,
            hobs,
            hexp,
            fis: (hexp > 0.0).then(|| 1.0 - hobs / hexp),
        });
    }
    Ok(out)
}
