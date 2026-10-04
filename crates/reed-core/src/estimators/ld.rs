//! Linkage disequilibrium: Rogers & Huff (2009) genotype r².

use crate::population::PopulationView;
use crate::window::WindowBuffer;

/// One valid site pair.
#[derive(Clone, Copy, Debug)]
pub struct LdPair {
    pub site_i: usize,
    pub site_j: usize,
    pub pos_i: u32,
    pub pos_j: u32,
    /// Individuals with both genotypes called.
    pub n_called: u64,
    pub r: f64,
    pub r2: f64,
}

/// Window-level LD summary over all site pairs.
#[derive(Clone, Copy, Debug)]
pub struct LdSummary {
    pub group: u32,
    pub n_sites: usize,
    /// Total pairs attempted: `k * (k - 1) / 2`.
    pub n_pairs: u64,
    /// Pairs with `n_called >= 2` and non-zero variance in both genotypes.
    pub n_valid: u64,
    pub mean_r2: Option<f64>,
}

/// ALT-allele dosage (0/1/2) of one sample at one site; `None` if missing.
fn dosage(window: &WindowBuffer, site: usize, sample: u32) -> Option<f64> {
    let (a, b) = window.genotype(site, sample as usize);
    (a >= 0 && b >= 0).then(|| (a + b) as f64)
}

/// Correlation of dosages over pairwise-complete individuals.
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
pub fn ld_r2_pairs(window: &WindowBuffer, pops: &PopulationView, group: u32) -> Vec<LdPair> {
    let samples = pops.samples(group);
    let k = window.n_sites();
    let mut out = Vec::new();
    for i in 0..k {
        for j in (i + 1)..k {
            if let Some((n_called, r, r2)) = pair_r(window, samples, i, j) {
                out.push(LdPair {
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
    out
}

/// Window-level r² summary without materialising every pair.
pub fn ld_r2_summary(window: &WindowBuffer, pops: &PopulationView, group: u32) -> LdSummary {
    let samples = pops.samples(group);
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
    LdSummary {
        group,
        n_sites: k,
        // k < 2 (e.g. the empty grid windows the reader emits) means no pairs.
        n_pairs: if k >= 2 { (k as u64) * (k as u64 - 1) / 2 } else { 0 },
        n_valid,
        mean_r2: (n_valid > 0).then(|| sum_r2 / n_valid as f64),
    }
}
