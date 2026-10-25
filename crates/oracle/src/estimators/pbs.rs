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
//! Site set: each pair's Hudson sums accumulate **independently per pair**
//! (scikit-allel `pbs` semantics): a site enters the A/B branch whenever A
//! and B each have ≥1 called allele, regardless of group C — `fst_ab` equals
//! a standalone `fst(A, B)` call on the same window, and a fully-uncalled
//! group only flags the site `n_skipped` while the two pairs not involving
//! it keep accumulating (04-algorithms.md §1; doc_v15 D15-12). An earlier
//! revision of this doc claimed a three-way intersection over the sites;
//! the code has been per-pair since it was written — doc corrected, no
//! behavior change (doc_v20 S79; the contract test
//! `s79_pbs_accumulates_per_pair_not_intersection` pins the rule).
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
    /// Normalised PBSn1 of Malaspinas et al. (2016); scikit-allel's default.
    PbsN1,
}

impl PbsEstimator {
    pub fn id(self) -> &'static str {
        match self {
            PbsEstimator::Yi2010 => "pbs.yi_2010",
            PbsEstimator::PbsN1 => "pbs.pbsn1",
        }
    }

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
    pub estimator: PbsEstimator,
    pub group_a: u32,
    pub group_b: u32,
    pub group_c: u32,
    pub contig: ContigId,
    pub start: u64,
    pub end: u64,
    /// Every site in the window (the census; per-pair site sets, so a site
    /// where one group is uncalled still feeds the other two pairs).
    pub n_sites: u64,
    /// Sites where at least one group is fully uncalled (flag only — the
    /// pairs not involving that group still accumulate; doc_v20 S79).
    pub n_skipped: u64,
    /// Window Hudson Fst per pair, clipped to `[0, 0.99999]`.
    pub fst_ab: Option<f64>,
    pub fst_ac: Option<f64>,
    pub fst_bc: Option<f64>,
    /// `T = −ln(1 − Fst)` per pair.
    pub t_ab: Option<f64>,
    pub t_ac: Option<f64>,
    pub t_bc: Option<f64>,
    /// The statistic; `None` when any pair-level Fst is undefined.
    pub pbs: Option<f64>,
    /// Pooled raw pair counts over **all** sites — skipped ones included, with
    /// missing pairs fully accounted (audit trail; fst.rs-aligned).
    pub within_a: PairCounts,
    pub within_b: PairCounts,
    pub within_c: PairCounts,
    pub between_ab: PairCounts,
    pub between_ac: PairCounts,
    pub between_bc: PairCounts,
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

fn pair_mpd(a: crate::window::AlleleCounts, b: crate::window::AlleleCounts) -> PairMpd {
    let n_a = a.ref_count + a.alt_count;
    let n_b = b.ref_count + b.alt_count;
    PairMpd {
        within_a: (n_a >= 2)
            .then(|| a.ref_count as f64 * a.alt_count as f64 / (n_a * (n_a - 1) / 2) as f64),
        within_b: (n_b >= 2)
            .then(|| b.ref_count as f64 * b.alt_count as f64 / (n_b * (n_b - 1) / 2) as f64),
        between: (n_a >= 1 && n_b >= 1).then(|| {
            (a.ref_count * b.alt_count + a.alt_count * b.ref_count) as f64 / (n_a * n_b) as f64
        }),
    }
}

/// Computes PBS among three distinct groups over a window.
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

    let (mut within_a, mut within_b, mut within_c) = (
        PairCounts::default(),
        PairCounts::default(),
        PairCounts::default(),
    );
    let (mut between_ab, mut between_ac, mut between_bc) = (
        PairCounts::default(),
        PairCounts::default(),
        PairCounts::default(),
    );
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

        // Every site counts; a fully-uncalled group flags the site skipped
        // but the two called pairs still accumulate independently — the
        // shipped PbsAcc semantics (04-algorithms.md §1; doc_v15 D15-12).
        n_sites += 1;
        if n_a == 0 || n_b == 0 || n_c == 0 {
            n_skipped += 1;
        }

        // allel inclusion rule (see fst.rs): `den` whenever both groups are
        // called; `num` additionally requires both within terms.
        let mpd_ab = pair_mpd(a, b);
        if let Some(bet) = mpd_ab.between {
            den_ab += bet;
            if let (Some(wa), Some(wb)) = (mpd_ab.within_a, mpd_ab.within_b) {
                num_ab += bet - (wa + wb) / 2.0;
            }
        }
        let mpd_ac = pair_mpd(a, c);
        if let Some(bet) = mpd_ac.between {
            den_ac += bet;
            if let (Some(wa), Some(wb)) = (mpd_ac.within_a, mpd_ac.within_b) {
                num_ac += bet - (wa + wb) / 2.0;
            }
        }
        let mpd_bc = pair_mpd(b, c);
        if let Some(bet) = mpd_bc.between {
            den_bc += bet;
            if let (Some(wa), Some(wb)) = (mpd_bc.within_a, mpd_bc.within_b) {
                num_bc += bet - (wa + wb) / 2.0;
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
