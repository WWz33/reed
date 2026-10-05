//! Population branching statistic (PBS).
//!
//! Window-level definition, matching scikit-allel `pbs` (MIT): each pair's
//! Hudson Fst uses **per-site-normalised sums** (`num = Σ(mpd_AB −
//! (mpd_A+mpd_B)/2)`, `den = Σ mpd_AB`, exactly `allel.hudson_fst`), is clipped
//! to `[0, 0.99999]` (T capped at ≈ 11.5), and converted via
//! `T = −ln(1 − Fst)`; then
//!
//! - [`PbsEstimator::Yi2010`]: `PBS = (T_ab + T_ac − T_bc) / 2`;
//! - [`PbsEstimator::PbsN1`]: `PBS / (1 + (T_ab + T_ac + T_bc) / 2)`
//!   (Malaspinas et al. 2016), the scikit-allel default.
//!
//! A window yields a value only when all three pairwise Fst denominators are
//! non-zero.
//!
//! Site set: a site enters the estimator only when **all three** groups have
//! at least one called allele (complete intersection — the three branch site
//! sets stay symmetric over the same sites). This differs from scikit-allel's
//! `pbs()`, which runs each pairwise `hudson_fst` independently (per-pair site
//! sets): under missing data, `fst_ab` here can differ from a standalone
//! `fst(A, B)` call on the same window (documented divergence; an
//! Intersection/PerPair switch is a v0.3 question — DESIGN §11).
//!
//! The pooled `PairCounts` audit trail covers **every** site, skipped ones
//! included, with missing pairs fully accounted (fst.rs-aligned; §12.3-1).

use crate::manifest::{Availability, SiteCompleteness};
use crate::population::PopulationView;
use crate::result::PairCounts;
use crate::window::{ContigId, WindowBuffer};

use super::{samples_of, Unavailable};

/// Sealed PBS definitions.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PbsEstimator {
    /// Original PBS of Yi et al. (2010).
    Yi2010,
    /// Normalised `PBSn1` of Malaspinas et al. (2016); scikit-allel's default.
    PbsN1,
}

impl PbsEstimator {
    /// Stable identifier for reports and golden tests.
    pub fn id(self) -> &'static str {
        match self {
            PbsEstimator::Yi2010 => "pbs.yi_2010",
            PbsEstimator::PbsN1 => "pbs.pbsn1",
        }
    }

    /// Literature the definition follows.
    pub fn citation(self) -> &'static str {
        match self {
            PbsEstimator::Yi2010 => {
                "Yi et al. 2010, Science 329:75 (original PBS; Fst via Hudson estimator)"
            }
            PbsEstimator::PbsN1 => {
                "Malaspinas et al. 2016, Nature 538:207 (normalised PBSn1; scikit-allel default)"
            }
        }
    }

    /// Whether this estimator is computable under a manifest.
    pub fn availability(self, manifest: SiteCompleteness) -> Availability {
        match manifest {
            SiteCompleteness::WithInvariants => Availability::Available,
            SiteCompleteness::VariantOnly(_) => Availability::Biased(
                "variant-only input: Fst (and therefore PBS) is computed on variant sites \
                 only; ascertainment bias is possible",
            ),
        }
    }
}

/// PBS result for three groups over one window.
#[derive(Clone, Copy, Debug)]
pub struct PbsResult {
    /// Which PBS definition produced the statistic.
    pub estimator: PbsEstimator,
    /// First group of the triple.
    pub group_a: u32,
    /// Second group of the triple.
    pub group_b: u32,
    /// Third group of the triple.
    pub group_c: u32,
    /// Contig the window is on.
    pub contig: ContigId,
    /// Window start (0-based, half-open).
    pub start: u64,
    /// Window end (exclusive).
    pub end: u64,
    /// Sites where all three groups have at least one called allele.
    pub n_sites: u64,
    /// Sites skipped for lack of coverage in at least one group.
    pub n_skipped: u64,
    /// Window Hudson Fst per pair, clipped to `[0, 0.99999]`.
    pub fst_ab: Option<f64>,
    /// Window Hudson Fst for the (A, C) pair, clipped to `[0, 0.99999]`.
    pub fst_ac: Option<f64>,
    /// Window Hudson Fst for the (B, C) pair, clipped to `[0, 0.99999]`.
    pub fst_bc: Option<f64>,
    /// `T = −ln(1 − Fst)` per pair.
    pub t_ab: Option<f64>,
    /// `T = −ln(1 − Fst)` for the (A, C) pair.
    pub t_ac: Option<f64>,
    /// `T = −ln(1 − Fst)` for the (B, C) pair.
    pub t_bc: Option<f64>,
    /// The statistic; `None` when any pair-level Fst is undefined.
    pub pbs: Option<f64>,
    /// Pooled raw pair counts within group A (audit trail; see [`PairCounts`]).
    pub within_a: PairCounts,
    /// Pooled raw pair counts within group B (audit trail).
    pub within_b: PairCounts,
    /// Pooled raw pair counts within group C (audit trail).
    pub within_c: PairCounts,
    /// Pooled raw cross-pair counts between groups A and B (audit trail).
    pub between_ab: PairCounts,
    /// Pooled raw cross-pair counts between groups A and C (audit trail).
    pub between_ac: PairCounts,
    /// Pooled raw cross-pair counts between groups B and C (audit trail).
    pub between_bc: PairCounts,
    /// Manifest caveat, if any.
    pub bias: Option<&'static str>,
}

impl PbsResult {
    /// The statistic; `None` when any pair-level Fst denominator was empty.
    pub fn value(&self) -> Option<f64> {
        self.pbs
    }
}

/// allel clips window Fst to [0, 0.99999] before the log transform
/// (selection.py:1322-1323), capping T at ≈ 11.5.
const FST_CLIP_MAX: f64 = 0.99999;

fn t_of(fst: f64) -> f64 {
    -(1.0 - fst).ln()
}

/// Per-site mean pairwise differences for one group pair.
struct PairMpd {
    within_a: Option<f64>,
    within_b: Option<f64>,
    between: Option<f64>,
}

/// Per-site mean pairwise differences for one group pair.
//
// Counts are bounded by the window's allele totals (≪ 2^53), so the ratio
// casts are exact in practice.
#[allow(clippy::cast_precision_loss)]
fn pair_mpd(a: crate::window::AlleleCounts, b: crate::window::AlleleCounts) -> PairMpd {
    let n_a = a.ref_count + a.alt_count;
    let n_b = b.ref_count + b.alt_count;
    PairMpd {
        within_a: (n_a >= 2).then(|| a.ref_count as f64 * a.alt_count as f64
            / (n_a * (n_a - 1) / 2) as f64),
        within_b: (n_b >= 2).then(|| b.ref_count as f64 * b.alt_count as f64
            / (n_b * (n_b - 1) / 2) as f64),
        between: (n_a >= 1 && n_b >= 1).then(|| {
            (a.ref_count * b.alt_count + a.alt_count * b.ref_count) as f64 / (n_a * n_b) as f64
        }),
    }
}

/// Computes PBS among three distinct groups over a window.
//
// - Casts: counts are bounded by the window's allele totals (≪ 2^53), so the
//   u64→f64 conversions in the ratio formulas are exact in practice.
// - Length: the pooled audit trail and the three pairwise sums share one
//   loop; splitting it would duplicate the per-site bookkeeping.
// - similar_names: the systematic A/B/C group notation (`fst_ab`, `t_ac`,
//   `between_bc`, …) mirrors the published formulas; any "distinct" rename
//   would obscure the math.
#[allow(
    clippy::cast_precision_loss,
    clippy::too_many_lines,
    clippy::similar_names
)]
pub fn pbs(
    window: &WindowBuffer,
    pops: &PopulationView,
    group_a: u32,
    group_b: u32,
    group_c: u32,
    estimator: PbsEstimator,
) -> Result<PbsResult, Unavailable> {
    if group_a == group_b || group_a == group_c || group_b == group_c {
        return Err(Unavailable {
            reason: "PBS requires three distinct groups",
        });
    }
    let sa = samples_of(pops, window, group_a)?;
    let sb = samples_of(pops, window, group_b)?;
    let sc = samples_of(pops, window, group_c)?;
    let bias = match estimator.availability(window.manifest) {
        Availability::Available => None,
        Availability::Biased(reason) => Some(reason),
        Availability::Unavailable(reason) => return Err(Unavailable { reason }),
    };

    let chrom_a = 2 * sa.len() as u64;
    let chrom_b = 2 * sb.len() as u64;
    let chrom_c = 2 * sc.len() as u64;
    let (within_total_a, within_total_b, within_total_c) = (
        chrom_a * chrom_a.saturating_sub(1) / 2,
        chrom_b * chrom_b.saturating_sub(1) / 2,
        chrom_c * chrom_c.saturating_sub(1) / 2,
    );
    let (between_total_ab, between_total_ac, between_total_bc) =
        (chrom_a * chrom_b, chrom_a * chrom_c, chrom_b * chrom_c);

    let (mut within_a, mut within_b, mut within_c) =
        (PairCounts::default(), PairCounts::default(), PairCounts::default());
    let (mut between_ab, mut between_ac, mut between_bc) =
        (PairCounts::default(), PairCounts::default(), PairCounts::default());
    let (mut n_sites, mut n_skipped) = (0u64, 0u64);
    // Per-site-normalised Hudson sums per pair (allel hudson_fst semantics).
    let (mut num_ab, mut den_ab) = (0.0f64, 0.0f64);
    let (mut num_ac, mut den_ac) = (0.0f64, 0.0f64);
    let (mut num_bc, mut den_bc) = (0.0f64, 0.0f64);

    for site in 0..window.n_sites() {
        let a = window.allele_counts(site, sa);
        let b = window.allele_counts(site, sb);
        let c = window.allele_counts(site, sc);
        let n_a = a.ref_count + a.alt_count;
        let n_b = b.ref_count + b.alt_count;
        let n_c = c.ref_count + c.alt_count;

        // Audit trail for every site, skipped or not: pooled integer counts
        // with missing pairs fully accounted (fst.rs-aligned).
        let comps_a = n_a * n_a.saturating_sub(1) / 2;
        let comps_b = n_b * n_b.saturating_sub(1) / 2;
        let comps_c = n_c * n_c.saturating_sub(1) / 2;
        within_a.comparable += comps_a;
        within_a.diffs += a.ref_count * a.alt_count;
        within_a.missing += within_total_a - comps_a;
        within_b.comparable += comps_b;
        within_b.diffs += b.ref_count * b.alt_count;
        within_b.missing += within_total_b - comps_b;
        within_c.comparable += comps_c;
        within_c.diffs += c.ref_count * c.alt_count;
        within_c.missing += within_total_c - comps_c;

        between_ab.comparable += n_a * n_b;
        between_ab.diffs += a.ref_count * b.alt_count + a.alt_count * b.ref_count;
        between_ab.missing += between_total_ab - n_a * n_b;
        between_ac.comparable += n_a * n_c;
        between_ac.diffs += a.ref_count * c.alt_count + a.alt_count * c.ref_count;
        between_ac.missing += between_total_ac - n_a * n_c;
        between_bc.comparable += n_b * n_c;
        between_bc.diffs += b.ref_count * c.alt_count + b.alt_count * c.ref_count;
        between_bc.missing += between_total_bc - n_b * n_c;

        if n_a == 0 || n_b == 0 || n_c == 0 {
            n_skipped += 1;
            continue;
        }
        n_sites += 1;

        // allel inclusion rule (see fst.rs): `den` whenever both groups are
        // called; `num` additionally requires both within terms.
        let mpd_ab = pair_mpd(a, b);
        if let Some(bet) = mpd_ab.between {
            den_ab += bet;
            if let (Some(wa), Some(wb)) = (mpd_ab.within_a, mpd_ab.within_b) {
                num_ab += bet - f64::midpoint(wa, wb);
            }
        }
        let mpd_ac = pair_mpd(a, c);
        if let Some(bet) = mpd_ac.between {
            den_ac += bet;
            if let (Some(wa), Some(wb)) = (mpd_ac.within_a, mpd_ac.within_b) {
                num_ac += bet - f64::midpoint(wa, wb);
            }
        }
        let mpd_bc = pair_mpd(b, c);
        if let Some(bet) = mpd_bc.between {
            den_bc += bet;
            if let (Some(wa), Some(wb)) = (mpd_bc.within_a, mpd_bc.within_b) {
                num_bc += bet - f64::midpoint(wa, wb);
            }
        }
    }

    let fst_ab = (den_ab != 0.0).then(|| (num_ab / den_ab).clamp(0.0, FST_CLIP_MAX));
    let fst_ac = (den_ac != 0.0).then(|| (num_ac / den_ac).clamp(0.0, FST_CLIP_MAX));
    let fst_bc = (den_bc != 0.0).then(|| (num_bc / den_bc).clamp(0.0, FST_CLIP_MAX));

    let (t_ab, t_ac, t_bc) = (fst_ab.map(t_of), fst_ac.map(t_of), fst_bc.map(t_of));

    let pbs = match (t_ab, t_ac, t_bc) {
        (Some(t_ab), Some(t_ac), Some(t_bc)) => {
            let raw = (t_ab + t_ac - t_bc) / 2.0;
            Some(match estimator {
                PbsEstimator::Yi2010 => raw,
                PbsEstimator::PbsN1 => raw / (1.0 + (t_ab + t_ac + t_bc) / 2.0),
            })
        }
        _ => None,
    };

    Ok(PbsResult {
        estimator,
        group_a,
        group_b,
        group_c,
        contig: window.contig,
        start: window.start,
        end: window.end,
        n_sites,
        n_skipped,
        fst_ab,
        fst_ac,
        fst_bc,
        t_ab,
        t_ac,
        t_bc,
        pbs,
        within_a,
        within_b,
        within_c,
        between_ab,
        between_ac,
        between_bc,
        bias,
    })
}
