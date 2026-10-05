//! Frequency-spectrum family: SFS, Watterson's theta, Tajima's D.
//!
//! All three operate on **complete-case sites only** for a fixed chromosome
//! count `n` (sites where any sample of the group has a missing allele are
//! skipped and counted). This keeps the spectrum definition unambiguous;
//! missing-aware variants are a v0.2 concern.
//!
//! "ALT" is treated as derived (no outgroup); the unfolded spectrum is the ALT
//! count spectrum under that assumption.
//!
//! Manifest gates (§4.2, judged per estimator): `θ_W` shares the
//! `AccessibleBp` denominator convention (window bp) and is **refused** under
//! `VariantOnly(Unknown)`. The SFS and Tajima's D are *invariant to the
//! presence of monomorphic sites* — monomorphic sites contribute nothing to
//! any spectrum bin, to S, or to π — so `VariantOnly(Complete)` input yields
//! values identical to all-sites input; `Filtered` (and `Unknown`, which
//! cannot be distinguished from it) carries a bias note.

use crate::estimators::Unavailable;
use crate::manifest::{Ascertainment, Availability, SiteCompleteness};
use crate::population::PopulationView;
use crate::window::{ContigId, WindowBuffer};

use super::samples_of;

/// Bias note for SFS / Tajima's D under `VariantOnly(Filtered/Unknown)`.
const SPECTRUM_BIAS: &str =
    "variant-only input: removed true variants are missing from the site set, \
     so bins (and S) undercount; VariantOnly(Complete) would be exact, and \
     Unknown cannot be distinguished from Filtered";

/// Sum of 1/i for i in 1..n (n = chromosome count); `a1` of Tajima's D.
// n is a chromosome count (2× samples), far below f64's exact-integer range.
#[allow(clippy::cast_precision_loss)]
fn a1(n_chroms: u64) -> f64 {
    (1..n_chroms).map(|i| 1.0 / i as f64).sum()
}

/// Sum of 1/i² for i in 1..n; `a2` of Tajima's D.
#[allow(clippy::cast_precision_loss)]
fn a2(n_chroms: u64) -> f64 {
    (1..n_chroms).map(|i| 1.0 / (i * i) as f64).sum()
}

/// Complete-case allele counts for a group at a site; `None` if any allele of
/// any sample of the group is missing.
fn complete_counts(window: &WindowBuffer, samples: &[u32], site: usize) -> Option<(u64, u64)> {
    let counts = window.allele_counts(site, samples);
    (counts.missing == 0).then_some((counts.ref_count, counts.alt_count))
}

/// `θ_W` shares the `AccessibleBp` denominator convention (window bp), so the
/// same four-way manifest matrix applies (§4.2).
fn theta_availability(manifest: SiteCompleteness) -> Availability {
    use Ascertainment::{Complete, Filtered, Unknown};
    use Availability::{Available, Biased, Unavailable};
    match manifest {
        SiteCompleteness::WithInvariants | SiteCompleteness::VariantOnly(Complete) => Available,
        SiteCompleteness::VariantOnly(Filtered) => Biased(
            "filtered variant-only input: removed true variants are indistinguishable \
             from invariant sites; the segregating-site numerator is undercounted",
        ),
        SiteCompleteness::VariantOnly(Unknown) => Unavailable(
            "variant-only input with unknown ascertainment: theta_W's denominator is \
             the window's base pairs; declare Complete/Filtered or provide an \
             all-sites VCF",
        ),
    }
}

/// Site frequency spectrum over a group.
#[derive(Clone, Debug)]
pub struct SfsResult {
    /// Contig the window is on.
    pub contig: ContigId,
    /// Window start (0-based, half-open).
    pub start: u64,
    /// Window end (exclusive).
    pub end: u64,
    /// Group the spectrum is built for.
    pub group: u32,
    /// Chromosome count `n` (2 per diploid sample).
    pub n_chroms: u64,
    /// Whether the spectrum is folded to the minor allele.
    pub folded: bool,
    /// `counts[i]` = number of sites with minor (folded) or ALT (unfolded)
    /// allele count `i + 1`; monomorphic complete sites are not counted.
    pub counts: Vec<u64>,
    /// Complete-case sites examined.
    pub n_sites: u64,
    /// Sites skipped because a group sample had a missing allele.
    pub n_skipped_missing: u64,
    /// Manifest caveat, if any (`VariantOnly(Filtered/Unknown)`).
    pub bias: Option<&'static str>,
}

/// Builds the (folded or unfolded) site frequency spectrum over segregating
/// complete-case sites.
//
// Indices derive from `n` (chromosome count, ≤ 2× samples); bins never
// approach `usize` limits on any supported target.
#[allow(clippy::cast_possible_truncation)]
pub fn sfs(
    window: &WindowBuffer,
    pops: &PopulationView,
    group: u32,
    folded: bool,
) -> Result<SfsResult, Unavailable> {
    let samples = samples_of(pops, window, group)?;
    let bias = super::monomorphic_invariant_note(window.manifest, SPECTRUM_BIAS);
    let n = 2 * samples.len() as u64;
    let n_bins = if folded { (n / 2) as usize } else { n.saturating_sub(1) as usize };
    let mut counts = vec![0u64; n_bins];
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
    Ok(SfsResult {
        contig: window.contig,
        start: window.start,
        end: window.end,
        group,
        n_chroms: n,
        folded,
        counts,
        n_sites,
        n_skipped_missing,
        bias,
    })
}

/// Watterson's theta over a group.
#[derive(Clone, Debug)]
pub struct WattersonResult {
    /// Contig the window is on.
    pub contig: ContigId,
    /// Window start (0-based, half-open).
    pub start: u64,
    /// Window end (exclusive).
    pub end: u64,
    /// Group the statistic is computed for.
    pub group: u32,
    /// Chromosome count `n` (2 per diploid sample).
    pub n_chroms: u64,
    /// Segregating (complete-case) sites `S`.
    pub segregating: u64,
    /// `a1 = Σ 1/i` for `i` in `1..n`, the Watterson normaliser.
    pub a_n: f64,
    /// Complete-case sites examined.
    pub n_sites: u64,
    /// Sites skipped because a group sample had a missing allele.
    pub n_skipped_missing: u64,
    /// Manifest caveat, if any (`VariantOnly(Filtered)`).
    pub bias: Option<&'static str>,
}

impl WattersonResult {
    /// Window span in base pairs.
    pub fn window_bp(&self) -> u64 {
        self.end - self.start
    }

    /// Window-total theta `S / a_n` — the exact quantity Tajima's D compares
    /// against π. `None` when the window has no complete-case sites (empty or
    /// all-missing input): an absent estimate, not a zero one (§4.4).
    //
    // `segregating` counts sites (≪ 2^53); the ratio cast is exact.
    #[allow(clippy::cast_precision_loss)]
    pub fn total_theta(&self) -> Option<f64> {
        (self.n_sites > 0 && self.a_n > 0.0).then(|| self.segregating as f64 / self.a_n)
    }

    /// Theta per base pair: `(S / a_n) / window_bp` (scikit-allel windowed
    /// semantics). Same `None` conditions as [`Self::total_theta`].
    #[allow(clippy::cast_precision_loss)]
    pub fn theta_per_bp(&self) -> Option<f64> {
        (self.n_sites > 0 && self.window_bp() > 0 && self.a_n > 0.0)
            .then(|| (self.segregating as f64 / self.a_n) / self.window_bp() as f64)
    }

    /// [`Self::theta_per_bp`] — the `value()` every window-level result
    /// exposes.
    pub fn value(&self) -> Option<f64> {
        self.theta_per_bp()
    }
}

/// Computes Watterson's theta from segregating complete-case sites.
///
/// Refused under `VariantOnly(Unknown)`: the per-bp denominator is a site
/// census the input cannot support (§4.2, same row as π `AccessibleBp`).
pub fn watterson_theta(
    window: &WindowBuffer,
    pops: &PopulationView,
    group: u32,
) -> Result<WattersonResult, Unavailable> {
    let samples = samples_of(pops, window, group)?;
    let bias = match theta_availability(window.manifest) {
        Availability::Available => None,
        Availability::Biased(reason) => Some(reason),
        Availability::Unavailable(reason) => return Err(Unavailable { reason }),
    };
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
    Ok(WattersonResult {
        contig: window.contig,
        start: window.start,
        end: window.end,
        group,
        n_chroms: n,
        segregating,
        a_n: a1(n),
        n_sites,
        n_skipped_missing,
        bias,
    })
}

/// Tajima's D over a group (Tajima 1989; standard moment variance).
#[derive(Clone, Debug)]
pub struct TajimaDResult {
    /// Contig the window is on.
    pub contig: ContigId,
    /// Window start (0-based, half-open).
    pub start: u64,
    /// Window end (exclusive).
    pub end: u64,
    /// Group the statistic is computed for.
    pub group: u32,
    /// Chromosome count `n` (2 per diploid sample).
    pub n_chroms: u64,
    /// Segregating (complete-case) sites `S`.
    pub segregating: u64,
    /// Sum over complete-case sites of the per-site mean pairwise difference
    /// (window-total pairs, **not** per bp).
    pub pi: f64,
    /// Window-total `S / a_n` (not per bp) — see
    /// [`WattersonResult::total_theta`].
    pub theta_w: f64,
    /// Tajima 1989 moment variance of `pi − theta_w`.
    pub variance: f64,
    /// `None` when `n < 2`, `S < min_sites`, or the variance is zero.
    pub d: Option<f64>,
    /// Minimum segregating sites required for a value.
    pub min_sites: u64,
    /// Complete-case sites examined.
    pub n_sites: u64,
    /// Sites skipped because a group sample had a missing allele.
    pub n_skipped_missing: u64,
    /// Manifest caveat, if any (`VariantOnly(Filtered/Unknown)`).
    pub bias: Option<&'static str>,
}

impl TajimaDResult {
    /// Tajima's D, or `None` per the documented guards.
    pub fn value(&self) -> Option<f64> {
        self.d
    }
}

/// Computes Tajima's D; `min_sites` is the minimum number of segregating sites
/// required (scikit-allel uses 3 by default).
//
// Counts are bounded by the window's sites and chromosome count (≪ 2^53),
// so the u64→f64 conversions in the moment-variance formulas are exact.
#[allow(clippy::cast_precision_loss)]
pub fn tajima_d(
    window: &WindowBuffer,
    pops: &PopulationView,
    group: u32,
    min_sites: u64,
) -> Result<TajimaDResult, Unavailable> {
    let samples = samples_of(pops, window, group)?;
    let bias = super::monomorphic_invariant_note(window.manifest, SPECTRUM_BIAS);
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

    // Standard moment variance (Tajima 1989), computed once and shared by the
    // guard and the struct field (§12.1 降级项: v0.2 去重).
    let (variance, d) = if n >= 2 && a1n > 0.0 {
        let nf = n as f64;
        let b1 = (nf + 1.0) / (3.0 * (nf - 1.0));
        let b2 = 2.0 * (nf * nf + nf + 3.0) / (9.0 * nf * (nf - 1.0));
        let c1 = b1 - 1.0 / a1n;
        let c2 = b2 - (nf + 2.0) / (a1n * nf) + a2n / (a1n * a1n);
        let e1 = c1 / a1n;
        let e2 = c2 / (a1n * a1n + a2n);
        let s = segregating as f64;
        let variance = e1 * s + e2 * s * (s - 1.0);
        let d = (segregating >= min_sites && variance > 0.0)
            .then(|| (pi - theta_w) / variance.sqrt());
        (variance, d)
    } else {
        (0.0, None)
    };

    Ok(TajimaDResult {
        contig: window.contig,
        start: window.start,
        end: window.end,
        group,
        n_chroms: n,
        segregating,
        pi,
        theta_w,
        variance,
        d,
        min_sites,
        n_sites,
        n_skipped_missing,
        bias,
    })
}
