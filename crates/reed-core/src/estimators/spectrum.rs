//! Frequency-spectrum family: SFS, Watterson's theta, Tajima's D.
//!
//! All three operate on **complete-case sites only** for a fixed chromosome
//! count `n` (sites where any sample of the group has a missing allele are
//! skipped and counted). This keeps the spectrum definition unambiguous;
//! missing-aware variants are a v0.2 concern.
//!
//! "ALT" is treated as derived (no outgroup); the unfolded spectrum is the ALT
//! count spectrum under that assumption.

use crate::population::PopulationView;
use crate::window::WindowBuffer;

/// Sum of 1/i for i in 1..n (n = chromosome count); `a1` of Tajima's D.
fn a1(n_chroms: u64) -> f64 {
    (1..n_chroms).map(|i| 1.0 / i as f64).sum()
}

/// Sum of 1/i² for i in 1..n; `a2` of Tajima's D.
fn a2(n_chroms: u64) -> f64 {
    (1..n_chroms).map(|i| 1.0 / (i * i) as f64).sum()
}

/// Complete-case allele counts for a group at a site; `None` if any allele of
/// any sample of the group is missing.
fn complete_counts(window: &WindowBuffer, samples: &[u32], site: usize) -> Option<(u64, u64)> {
    let counts = window.allele_counts(site, samples);
    (counts.missing == 0).then_some((counts.ref_count, counts.alt_count))
}

/// Site frequency spectrum over a group.
#[derive(Clone, Debug)]
pub struct SfsResult {
    pub group: u32,
    pub n_chroms: u64,
    pub folded: bool,
    /// `counts[i]` = number of sites with minor (folded) or ALT (unfolded)
    /// allele count `i + 1`; monomorphic complete sites are not counted.
    pub counts: Vec<u64>,
    /// Complete-case sites examined.
    pub n_sites: u64,
    /// Sites skipped because a group sample had a missing allele.
    pub n_skipped_missing: u64,
}

/// Builds the (folded or unfolded) site frequency spectrum over segregating
/// complete-case sites.
pub fn sfs(
    window: &WindowBuffer,
    pops: &PopulationView,
    group: u32,
    folded: bool,
) -> SfsResult {
    let samples = pops.samples(group);
    let n = 2 * samples.len() as u64;
    let bins = if folded { (n / 2) as usize } else { n.saturating_sub(1) as usize };
    let mut counts = vec![0u64; bins];
    let (mut n_sites, mut n_skipped_missing) = (0u64, 0u64);
    if n >= 2 {
        for site in 0..window.n_sites() {
            match complete_counts(window, samples, site) {
                None => n_skipped_missing += 1,
                Some((_ref_count, alt)) => {
                    n_sites += 1;
                    if alt > 0 && alt < n {
                        let index = if folded { alt.min(n - alt) } else { alt };
                        counts[(index - 1) as usize] += 1;
                    }
                }
            }
        }
    }
    SfsResult {
        group,
        n_chroms: n,
        folded,
        counts,
        n_sites,
        n_skipped_missing,
    }
}

/// Watterson's theta over a group.
#[derive(Clone, Debug)]
pub struct WattersonResult {
    pub group: u32,
    pub n_chroms: u64,
    pub segregating: u64,
    pub a_n: f64,
    pub window_bp: u64,
    pub n_sites: u64,
    pub n_skipped_missing: u64,
}

impl WattersonResult {
    /// Theta per base pair: `(S / a_n) / window_bp`.
    pub fn value(&self) -> Option<f64> {
        (self.window_bp > 0 && self.a_n > 0.0)
            .then(|| (self.segregating as f64 / self.a_n) / self.window_bp as f64)
    }
}

/// Computes Watterson's theta from segregating complete-case sites.
pub fn watterson_theta(
    window: &WindowBuffer,
    pops: &PopulationView,
    group: u32,
) -> WattersonResult {
    let samples = pops.samples(group);
    let n = 2 * samples.len() as u64;
    let (mut segregating, mut n_sites, mut n_skipped_missing) = (0u64, 0u64, 0u64);
    if n >= 2 {
        for site in 0..window.n_sites() {
            match complete_counts(window, samples, site) {
                None => n_skipped_missing += 1,
                Some((_ref_count, alt)) => {
                    n_sites += 1;
                    if alt > 0 && alt < n {
                        segregating += 1;
                    }
                }
            }
        }
    }
    WattersonResult {
        group,
        n_chroms: n,
        segregating,
        a_n: a1(n),
        window_bp: window.window_bp(),
        n_sites,
        n_skipped_missing,
    }
}

/// Tajima's D over a group (Tajima 1989; standard moment variance).
#[derive(Clone, Debug)]
pub struct TajimaDResult {
    pub group: u32,
    pub n_chroms: u64,
    pub segregating: u64,
    /// Sum over complete-case sites of the per-site mean pairwise difference.
    pub pi: f64,
    pub theta_w: f64,
    pub variance: f64,
    /// `None` when `n < 2`, `S < min_sites`, or the variance is zero.
    pub d: Option<f64>,
    pub min_sites: u64,
    pub n_sites: u64,
    pub n_skipped_missing: u64,
}

impl TajimaDResult {
    pub fn value(&self) -> Option<f64> {
        self.d
    }
}

/// Computes Tajima's D; `min_sites` is the minimum number of segregating sites
/// required (scikit-allel uses 3 by default).
pub fn tajima_d(
    window: &WindowBuffer,
    pops: &PopulationView,
    group: u32,
    min_sites: u64,
) -> TajimaDResult {
    let samples = pops.samples(group);
    let n = 2 * samples.len() as u64;
    let (mut segregating, mut pi, mut n_sites, mut n_skipped_missing) = (0u64, 0.0f64, 0u64, 0u64);

    if n >= 2 {
        let pairs = n * (n - 1) / 2;
        for site in 0..window.n_sites() {
            match complete_counts(window, samples, site) {
                None => n_skipped_missing += 1,
                Some((ref_count, alt)) => {
                    n_sites += 1;
                    if alt > 0 && alt < n {
                        segregating += 1;
                    }
                    pi += (ref_count * alt) as f64 / pairs as f64;
                }
            }
        }
    }

    let a1n = a1(n);
    let a2n = a2(n);
    let theta_w = if a1n > 0.0 { segregating as f64 / a1n } else { 0.0 };

    let d = if n >= 2 && a1n > 0.0 {
        let nf = n as f64;
        let b1 = (nf + 1.0) / (3.0 * (nf - 1.0));
        let b2 = 2.0 * (nf * nf + nf + 3.0) / (9.0 * nf * (nf - 1.0));
        let c1 = b1 - 1.0 / a1n;
        let c2 = b2 - (nf + 2.0) / (a1n * nf) + a2n / (a1n * a1n);
        let e1 = c1 / a1n;
        let e2 = c2 / (a1n * a1n + a2n);
        let s = segregating as f64;
        let variance = e1 * s + e2 * s * (s - 1.0);
        (segregating >= min_sites && variance > 0.0)
            .then(|| (pi - theta_w) / variance.sqrt())
    } else {
        None
    };

    TajimaDResult {
        group,
        n_chroms: n,
        segregating,
        pi,
        theta_w,
        variance: {
            let nf = n as f64;
            if n >= 2 && a1n > 0.0 {
                let b1 = (nf + 1.0) / (3.0 * (nf - 1.0));
                let b2 = 2.0 * (nf * nf + nf + 3.0) / (9.0 * nf * (nf - 1.0));
                let c1 = b1 - 1.0 / a1n;
                let c2 = b2 - (nf + 2.0) / (a1n * nf) + a2n / (a1n * a1n);
                let e1 = c1 / a1n;
                let e2 = c2 / (a1n * a1n + a2n);
                let s = segregating as f64;
                e1 * s + e2 * s * (s - 1.0)
            } else {
                0.0
            }
        },
        d,
        min_sites,
        n_sites,
        n_skipped_missing,
    }
}
