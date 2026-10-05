//! Linkage disequilibrium: Rogers & Huff (2009) genotype r².
//!
//! Manifest gate (§4.2, judged per estimator): a pair's r² depends only on
//! the two sites' dosages, and a monomorphic site can never form a valid pair
//! (zero variance) — the valid-pair set is **invariant to the presence of
//! monomorphic sites**, so `VariantOnly(Complete)` input yields results
//! identical to all-sites input. `Filtered` — and `Unknown`, which cannot be
//! distinguished from it — loses pairs to removed true variants, so the
//! window summary carries a bias note. Per-pair values are pair-local and
//! exact wherever both sites are present.

use crate::estimators::Unavailable;
use crate::population::PopulationView;
use crate::window::{ContigId, WindowBuffer};

use super::samples_of;

/// Bias note for the window summary under `VariantOnly(Filtered/Unknown)`.
const LD_BIAS: &str =
    "variant-only input: pairs involving removed true variants are missing; \
     r2 is exact for pairs fully present, and Unknown cannot be distinguished \
     from Filtered";

/// One valid site pair.
#[derive(Clone, Copy, Debug)]
pub struct LdPair {
    /// Contig both sites are on.
    pub contig: ContigId,
    /// 0-based index of the first site in the window.
    pub site_i: usize,
    /// 0-based index of the second site in the window.
    pub site_j: usize,
    /// 0-based position of the first site.
    pub pos_i: u32,
    /// 0-based position of the second site.
    pub pos_j: u32,
    /// Individuals with both genotypes called.
    pub n_called: u64,
    /// Dosage correlation over pairwise-complete individuals.
    pub r: f64,
    /// Squared dosage correlation.
    pub r2: f64,
}

/// Window-level LD summary over all site pairs.
#[derive(Clone, Copy, Debug)]
pub struct LdSummary {
    /// Contig the window is on.
    pub contig: ContigId,
    /// Window start (0-based, half-open).
    pub start: u64,
    /// Window end (exclusive).
    pub end: u64,
    /// Group the dosages come from.
    pub group: u32,
    /// Sites in the window.
    pub n_sites: usize,
    /// Total pairs attempted: `k * (k - 1) / 2`.
    pub n_pairs: u64,
    /// Pairs with `n_called >= 2` and non-zero variance in both genotypes.
    pub n_valid: u64,
    /// Mean r² over valid pairs, `None` when no pair is valid.
    pub mean_r2: Option<f64>,
    /// Manifest caveat, if any (`VariantOnly(Filtered/Unknown)`).
    pub bias: Option<&'static str>,
}

impl LdSummary {
    /// The window statistic: mean r² over valid pairs.
    pub fn value(&self) -> Option<f64> {
        self.mean_r2
    }
}

/// ALT-allele dosage (0/1/2) of one sample at one site; `None` if missing.
fn dosage(window: &WindowBuffer, site: usize, sample: u32) -> Option<f64> {
    let (a, b) = window.genotype(site, sample as usize);
    (a >= 0 && b >= 0).then(|| f64::from(a + b))
}

/// Correlation of dosages over pairwise-complete individuals.
//
// `n` accumulates one per called sample (a non-negative integer count, at
// most `samples.len()`), so the f64→u64 conversion is exact.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn pair_r(
    window: &WindowBuffer,
    samples: &[u32],
    i: usize,
    j: usize,
) -> Option<(u64, f64, f64)> {
    let (mut n, mut sx, mut sy, mut sxx, mut syy, mut sxy) = (0f64, 0f64, 0f64, 0f64, 0f64, 0f64);
    for &s in samples {
        if let (Some(x), Some(y)) = (dosage(window, i, s), dosage(window, j, s)) {
            n += 1.0;
            sx += x;
            sy += y;
            sxx += x * x;
            syy += y * y;
            sxy += x * y;
        }
    }
    if n < 2.0 {
        return None;
    }
    let cov = sxy - sx * sy / n;
    let vx = sxx - sx * sx / n;
    let vy = syy - sy * sy / n;
    if vx <= 0.0 || vy <= 0.0 {
        return None;
    }
    let r = cov / (vx.sqrt() * vy.sqrt());
    let r = r.clamp(-1.0, 1.0);
    Some((n as u64, r, r * r))
}

/// All valid site pairs with their r and r². Cost is O(k²) in the number of
/// sites; for large windows prefer [`ld_r2_summary`].
pub fn ld_r2_pairs(
    window: &WindowBuffer,
    pops: &PopulationView,
    group: u32,
) -> Result<Vec<LdPair>, Unavailable> {
    let samples = samples_of(pops, window, group)?;
    let k = window.n_sites();
    let mut out = Vec::new();
    for i in 0..k {
        for j in (i + 1)..k {
            if let Some((n_called, r, r2)) = pair_r(window, samples, i, j) {
                out.push(LdPair {
                    contig: window.contig,
                    site_i: i,
                    site_j: j,
                    pos_i: window.pos[i],
                    pos_j: window.pos[j],
                    n_called,
                    r,
                    r2,
                });
            }
        }
    }
    Ok(out)
}

/// Window-level r² summary without materialising every pair.
pub fn ld_r2_summary(
    window: &WindowBuffer,
    pops: &PopulationView,
    group: u32,
) -> Result<LdSummary, Unavailable> {
    let samples = samples_of(pops, window, group)?;
    let bias = super::monomorphic_invariant_note(window.manifest, LD_BIAS);
    let k = window.n_sites();
    let mut n_valid = 0u64;
    let mut sum_r2 = 0.0f64;
    for i in 0..k {
        for j in (i + 1)..k {
            if let Some((_n, _r, r2)) = pair_r(window, samples, i, j) {
                n_valid += 1;
                sum_r2 += r2;
            }
        }
    }
    // `n_valid` counts pairs (≪ 2^53); the ratio cast is exact.
    #[allow(clippy::cast_precision_loss)]
    Ok(LdSummary {
        contig: window.contig,
        start: window.start,
        end: window.end,
        group,
        n_sites: k,
        // k < 2 (e.g. the empty grid windows the reader emits) means no pairs.
        n_pairs: if k >= 2 { (k as u64) * (k as u64 - 1) / 2 } else { 0 },
        n_valid,
        mean_r2: (n_valid > 0).then(|| sum_r2 / n_valid as f64),
        bias,
    })
}
