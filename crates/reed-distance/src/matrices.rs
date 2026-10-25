//! The individual-by-individual estimators.
//!
//! Every formula here was read out of a reference implementation; the
//! citation in each doc comment points at `docs/rounds/v17/01-research.md`
//! and, through it, at the exact source line. Do not "simplify" a formula
//! without re-checking it against the reference: the whole point of this
//! module is bit-comparable output.
//!
//! Parallelism (doc_v20 R13): each estimator fans its outer sample loop out
//! over rayon — one task per matrix row, rows collected in index order and
//! scattered serially. No float is ever reduced across threads and every
//! `(i,j)` entry keeps exactly its serial arithmetic, so the output is
//! bit-identical at any pool size (locked by
//! `r13_estimator_results_thread_invariant`). The estimators run on
//! whatever pool is ambient: the rayon global pool by default, or the
//! caller's scoped pool when invoked inside one (e.g. a budget-sized
//! `rayon::ThreadPool::install`).

use rayon::prelude::*;

use crate::standard::Standardized;

/// Which estimator produced a matrix.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MatrixKind {
    /// plink `--distance` (allele-count units).
    AlleleCountDistance,
    /// plink `--distance ibs` (identity by state, proportion).
    Ibs,
    /// plink `--distance 1-ibs` (genomic distance, proportion).
    OneMinusIbs,
    /// plink2 `--make-rel` / GCTA (variance-standardized relationship).
    Grm,
    /// plink2 `--make-rel cov` (mean-centred covariance).
    Covariance,
    /// plink2 `--make-king` (KING-robust kinship).
    King,
    /// plink `--genome` (PI_HAT).
    PiHat,
    /// vcftools `--relatedness2` (Manichaikul phi).
    Phi,
    /// vcftools `--relatedness` (Yang Ajk).
    Ajk,
    /// bcftools `gtcheck` (discordance).
    Discordance,
}

impl MatrixKind {
    /// The output suffix used by the reference tool.
    pub fn suffix(self) -> &'static str {
        match self {
            MatrixKind::AlleleCountDistance => "dist",
            MatrixKind::Ibs => "mibs",
            MatrixKind::OneMinusIbs => "mdist",
            MatrixKind::Grm | MatrixKind::Covariance => "rel",
            MatrixKind::King => "king",
            MatrixKind::PiHat => "genome",
            MatrixKind::Phi => "relatedness2",
            MatrixKind::Ajk => "relatedness",
            MatrixKind::Discordance => "gtcheck",
        }
    }
}

/// A complete symmetric n x n matrix (row-major) plus its provenance.
#[derive(Clone, Debug)]
pub struct MatrixResult {
    /// `n * n` values, row-major, symmetric.
    pub values: Vec<f64>,
    /// Number of individuals.
    pub n: usize,
    /// Which estimator produced it.
    pub kind: MatrixKind,
    /// Sample names in row/column order.
    pub ids: Vec<String>,
    /// Per-pair auxiliary counts, for the estimators that report them
    /// (bcftools' `ncnt`/`nmatch`, vcftools' `N_AaAa`/`N_AAaa`).
    pub aux: Option<AuxCounts>,
}

/// Auxiliary per-pair counts some formats print.
#[derive(Clone, Debug, Default)]
pub struct AuxCounts {
    /// Number of sites compared for the pair.
    pub n_compared: Vec<u64>,
    /// Number of matching genotypes for the pair.
    pub n_matching: Vec<u64>,
    /// Sites where both are heterozygous.
    pub n_hethet: Vec<u64>,
    /// Sites where both are homozygous and differ.
    pub n_opposite_hom: Vec<u64>,
    /// Marginal heterozygous count per individual (vcftools `N_Aa`).
    pub n_het: Vec<u64>,
    /// plink `--genome` Z0/Z1/Z2 per pair (row-major n*n).
    pub z0: Vec<f64>,
    pub z1: Vec<f64>,
    pub z2: Vec<f64>,
}

impl MatrixResult {
    fn new(n: usize, kind: MatrixKind, ids: Vec<String>) -> Self {
        Self { values: vec![0.0; n * n], n, kind, ids, aux: None }
    }

    #[inline]
    #[allow(dead_code)]
    fn get(&self, i: usize, j: usize) -> f64 {
        self.values[i * self.n + j]
    }

    #[inline]
    fn set(&mut self, i: usize, j: usize, v: f64) {
        self.values[i * self.n + j] = v;
        self.values[j * self.n + i] = v;
    }
}

/// Options shared by the distance-family estimators.
#[derive(Clone, Copy, Debug)]
pub struct DistanceOptions {
    /// Apply plink's HWE-weighted missing correction (default `true`).
    ///
    /// plink's `flat-missing` turns this off.
    pub hwe_missing_correction: bool,
}

impl Default for DistanceOptions {
    fn default() -> Self {
        Self { hwe_missing_correction: true }
    }
}

/// plink `--distance`: `A_ij = Σ |dosage_i − dosage_j|` over jointly called
/// sites, with plink's HWE-weighted missing correction (01 §1.1-1.2).
pub fn allele_ct_distance(
    st: &Standardized,
    opts: &DistanceOptions,
) -> MatrixResult {
    let n = st.n_samples();
    let m = st.n_sites;
    let mut out = MatrixResult::new(n, MatrixKind::AlleleCountDistance, Vec::new());

    // Raw sums and per-sample missing weights.
    let (weights, weight_sum) = if opts.hwe_missing_correction {
        hwe_weights(st)
    } else {
        (vec![1.0; m], m as f64)
    };
    let mut miss = vec![0.0f64; n];
    for (i, slot) in miss.iter_mut().enumerate() {
        for j in 0..m {
            if !st.is_called(i, j) {
                *slot += weights[j];
            }
        }
    }

    // Per-sample dosages, hoisted out of the pair loop (doc_v20 R13): the
    // loop used to recompute `w*sd + mean` for both members of every pair
    // at every site — O(n²·m) extra multiplications. The hoisted value is
    // that same expression on the same operands, so each entry below stays
    // bit-identical to the serial recomputation. (Entries for missing
    // calls are computed but never read, exactly as before.)
    // doc_v20 review of R13: stored as i8 (the value is 0/1/2) instead of
    // f64. At 2,000 samples x 1M sites the f64 copy alone was 16 GB on top
    // of `w` and `called`; widening back to f64 in the loop is exact, so
    // every entry stays bit-identical.
    let dose: Vec<i8> = (0..n * m)
        .map(|k| {
            let (i, s) = (k / m, k % m);
            let v = st.w.get(i, s) * st.sd[s] + st.mean[s];
            if v.is_nan() {
                0
            } else {
                v.round().clamp(-128.0, 127.0) as i8
            }
        })
        .collect();

    // One strict-lower-triangle row per task (doc_v20 R13): task `i`
    // produces row `i` only, rows are collected in index order and
    // scattered serially, so the result never depends on the pool size.
    let rows: Vec<Vec<f64>> = (0..n)
        .into_par_iter()
        .map(|i| {
            (0..i)
                .map(|j| {
                    let mut raw = 0.0f64;
                    let mut miss_both = 0.0f64;
                    for s in 0..m {
                        let ci = st.is_called(i, s);
                        let cj = st.is_called(j, s);
                        match (ci, cj) {
                            (true, true) => {
                                // The numerator is the *unweighted* allele-count
                                // difference; plink's HWE weights enter only through
                                // the missing-correction denominator
                                // (1.9/plink_calc.c:8169-8199 vs incr_dists_i at
                                // :973-1024, which popcounts without weighting).
                                raw += f64::from(dose[i * m + s] - dose[j * m + s]).abs();
                            }
                            (false, false) => miss_both += weights[s],
                            _ => {}
                        }
                    }
                    let denom = weight_sum - miss[i] - miss[j] + miss_both;
                    if denom > 0.0 {
                        raw * weight_sum / denom
                    } else {
                        0.0
                    }
                })
                .collect()
        })
        .collect();
    for (i, row) in rows.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            out.set(i, j, v);
        }
    }
    out
}

/// plink's per-variant missing weights: `p(1-p)(p²-p+1)`, rescaled so the
/// total is just under 2^32 (01 §1.2).
fn hwe_weights(st: &Standardized) -> (Vec<f64>, f64) {
    let m = st.n_sites;
    let mut w: Vec<f64> = st
        .freq
        .iter()
        .map(|&p| {
            let q = 1.0 - p;
            (p * q * (p * p - p + 1.0)).max(0.0)
        })
        .collect();
    let total: f64 = w.iter().sum();
    if total <= 0.0 {
        return (vec![1.0; m], m as f64);
    }
    // plink: dyy = (2^32 - marker_ct)/dyy, uii = round(weight*dyy).
    let scale = (4_294_967_296.0 - m as f64) / total;
    let mut sum = 0.0;
    for v in &mut w {
        *v = (*v * scale).round();
        sum += *v;
    }
    (w, sum)
}

/// plink `--distance ibs`: `1 − d·0.5/m` (01 §1.3).
pub fn ibs(st: &Standardized, opts: &DistanceOptions) -> MatrixResult {
    let d = allele_ct_distance(st, opts);
    let m = st.n_sites as f64;
    let mut out = MatrixResult::new(d.n, MatrixKind::Ibs, Vec::new());
    out.values = d
        .values
        .iter()
        .map(|v| 1.0 - v * 0.5 / m)
        .collect();
    out
}

/// plink `--distance 1-ibs`: `d·0.5/m` (01 §1.3).
pub fn one_minus_ibs(st: &Standardized, opts: &DistanceOptions) -> MatrixResult {
    let d = allele_ct_distance(st, opts);
    let m = st.n_sites as f64;
    let mut out = MatrixResult::new(d.n, MatrixKind::OneMinusIbs, Vec::new());
    out.values = d.values.iter().map(|v| v * 0.5 / m).collect();
    out
}

/// Missing-data policy for the relationship matrices.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MissingPolicy {
    /// Each pair divided by the number of jointly called sites
    /// (plink/GCTA default).
    Pairwise,
    /// Every entry divided by the total site count (plink `meanimpute`).
    MeanImpute,
}

/// plink2 `--make-rel` / GCTA relationship matrix (01 §2).
///
/// `variance_standardize = false` gives `--make-rel cov` (mean-centred
/// covariance).
pub fn grm(
    st: &Standardized,
    missing: MissingPolicy,
    variance_standardize: bool,
) -> MatrixResult {
    let n = st.n_samples();
    let m = st.n_sites;
    let kind = if variance_standardize {
        MatrixKind::Grm
    } else {
        MatrixKind::Covariance
    };
    let mut out = MatrixResult::new(n, kind, Vec::new());

    // `w` already holds (dosage - mean)/sd, so the variance-standardized
    // matrix uses it as-is. The covariance variant undoes the division by sd
    // to get (dosage - mean).
    let mut x = vec![0.0f64; n * m];
    for i in 0..n {
        for j in 0..m {
            if !st.is_called(i, j) {
                continue;
            }
            x[i * m + j] = if variance_standardize {
                st.w.get(i, j)
            } else {
                st.w.get(i, j) * st.sd[j]
            };
        }
    }

    // One triangle row (diagonal included) per task (doc_v20 R13); the
    // O(m) `called_both` second pass per pair keeps its math intact. The
    // per-sample `miss` counts the serial loop kept only to silence
    // unused warnings are gone with it: both denominators below read
    // `called_both`/`m` directly.
    let rows: Vec<Vec<f64>> = (0..n)
        .into_par_iter()
        .map(|i| {
            (0..=i)
                .map(|j| {
                    let mut acc = 0.0f64;
                    for s in 0..m {
                        acc += x[i * m + s] * x[j * m + s];
                    }
                    let denom = match missing {
                        MissingPolicy::MeanImpute => m as f64,
                        MissingPolicy::Pairwise => st.called_both(i, j) as f64,
                    };
                    if denom > 0.0 { acc / denom } else { 0.0 }
                })
                .collect()
        })
        .collect();
    for (i, row) in rows.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            out.set(i, j, v);
        }
    }
    out
}

/// Counts behind the KING estimator (plink2 `--make-king`, 01 §3).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct KingCounts {
    /// Opposite homozygotes.
    pub ibs0: u64,
    /// Both heterozygous.
    pub hethet: u64,
    /// i het, j hom-alt.
    pub het1hom2: u64,
    /// i hom-alt, j het.
    pub het2hom1: u64,
}

impl KingCounts {
    /// `0.5 − (4·IBS0 + H1H2 + H2H1) / (4·(HETHET + min(H1H2, H2H1)))`.
    pub fn kinship(&self) -> f64 {
        let smaller_het = self.hethet + self.het1hom2.min(self.het2hom1);
        if smaller_het == 0 {
            return f64::NAN;
        }
        0.5 - (4.0 * self.ibs0 as f64 + self.het1hom2 as f64 + self.het2hom1 as f64)
            / (4.0 * smaller_het as f64)
    }
}

/// plink2 `--make-king`: KING-robust kinship (01 §3.1).
///
/// The caller is responsible for restricting to autosomes, as plink2 does.
pub fn king(st: &Standardized) -> (MatrixResult, Vec<KingCounts>) {
    let n = st.n_samples();
    let m = st.n_sites;
    let mut out = MatrixResult::new(n, MatrixKind::King, Vec::new());
    let mut counts = vec![KingCounts::default(); n * n];

    // Dosages 0/1/2, -1 for missing.
    let mut dose = vec![0i8; n * m];
    for i in 0..n {
        for j in 0..m {
            if st.is_called(i, j) {
                dose[i * m + j] =
                    (st.w.get(i, j) * st.sd[j] + st.mean[j]).round() as i8;
            } else {
                dose[i * m + j] = -1;
            }
        }
    }

    // One triangle row of values plus its counts per task (doc_v20 R13).
    let rows: Vec<(Vec<f64>, Vec<KingCounts>)> = (0..n)
        .into_par_iter()
        .map(|i| {
            let mut vals = Vec::with_capacity(i + 1);
            let mut cs = Vec::with_capacity(i + 1);
            for j in 0..=i {
                let mut c = KingCounts::default();
                for s in 0..m {
                    let a = dose[i * m + s];
                    let b = dose[j * m + s];
                    if a < 0 || b < 0 {
                        continue;
                    }
                    match (a, b) {
                        (1, 1) => c.hethet += 1,
                        (0, 2) | (2, 0) => c.ibs0 += 1,
                        (1, 0) | (1, 2) => c.het1hom2 += 1,
                        (0, 1) | (2, 1) => c.het2hom1 += 1,
                        _ => {}
                    }
                }
                vals.push(c.kinship());
                cs.push(c);
            }
            (vals, cs)
        })
        .collect();
    for (i, (vals, cs)) in rows.iter().enumerate() {
        for (j, (&v, &c)) in vals.iter().zip(cs.iter()).enumerate() {
            out.set(i, j, v);
            counts[i * n + j] = c;
            counts[j * n + i] = c;
        }
    }
    (out, counts)
}

/// Options for the PI_HAT estimator.
#[derive(Clone, Copy, Debug, Default)]
pub struct PiHatOptions {
    /// Skip the [0,1] clipping (plink `unbounded`).
    pub unbounded: bool,
    /// Snap to a theoretically possible configuration when
    /// `PI_HAT² < Z2` (plink `nudge`).
    pub nudge: bool,
}

/// The plink `--genome` quantities for one pair.
#[derive(Clone, Copy, Debug, Default)]
pub struct PiHatResult {
    /// P(IBD=0).
    pub z0: f64,
    /// P(IBD=1).
    pub z1: f64,
    /// P(IBD=2).
    pub z2: f64,
    /// `Z1/2 + Z2`.
    pub pi_hat: f64,
    /// Number of jointly called sites.
    pub n_valid: u64,
    /// Sites with one IBS allele (het/hom or hom/het).
    pub ibs1: u64,
    /// Sites with zero IBS alleles (opposite homozygotes).
    pub ibs0: u64,
}

/// plink `--genome` PI_HAT (01 §4).
///
/// The expected-IBS coefficients `e00..e12` are computed from the allele
/// frequencies as plink does, then averaged over the sites in use.
pub fn pi_hat(st: &Standardized, opts: &PiHatOptions) -> (MatrixResult, Vec<PiHatResult>) {
    let n = st.n_samples();
    let m = st.n_sites;
    let mut out = MatrixResult::new(n, MatrixKind::PiHat, Vec::new());
    let mut results = vec![PiHatResult::default(); n * n];

    let mut dose = vec![0i8; n * m];
    for i in 0..n {
        for j in 0..m {
            dose[i * m + j] = if st.is_called(i, j) {
                (st.w.get(i, j) * st.sd[j] + st.mean[j]).round() as i8
            } else {
                -1
            };
        }
    }

    // Expected IBS counts given IBD state, averaged over sites (plink
    // 1.9/plink_calc.c:4848-4868). For a site with allele frequency p:
    //   E[IBS=0 | IBD=0] = 2 p²(1-p)²
    //   E[IBS=0 | IBD=2] = 0
    // and the coefficients are the sums over sites, later normalized.
    let (e00, e01, e02, e11, e12) = expected_ibs_coefficients(st);

    // One triangle row of (PI_HAT, triplet) per task (doc_v20 R13).
    let rows: Vec<Vec<(f64, PiHatResult)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            (0..=i)
                .map(|j| {
                    let mut n_valid = 0u64;
                    let mut ibs0 = 0u64;
                    let mut ibs1 = 0u64;
                    for s in 0..m {
                        let a = dose[i * m + s];
                        let b = dose[j * m + s];
                        if a < 0 || b < 0 {
                            continue;
                        }
                        n_valid += 1;
                        let d = (a - b).unsigned_abs();
                        match d {
                            0 => {}
                            1 => ibs1 += 1,
                            _ => ibs0 += 1,
                        }
                    }
                    let nn = n_valid as f64;
                    let oo = nn - ibs1 as f64 - ibs0 as f64;
                    let (mut z0, mut z1, mut z2) = if nn > 0.0 {
                        let z0 = ibs0 as f64 / (e00 * nn);
                        let z1 = (ibs1 as f64 - z0 * e01 * nn) / (e11 * nn);
                        let z2 = (oo - nn * (z0 * e02 + z1 * e12)) / nn;
                        (z0, z1, z2)
                    } else {
                        (0.0, 0.0, 0.0)
                    };
                    if !opts.unbounded {
                        (z0, z1, z2) = clip_ibd(z0, z1, z2);
                    }
                    let mut pi = z1 * 0.5 + z2;
                    if opts.nudge && pi * pi < z2 {
                        z0 = (1.0 - pi) * (1.0 - pi);
                        z1 = 2.0 * pi * (1.0 - pi);
                        z2 = pi * pi;
                        pi = z1 * 0.5 + z2;
                    }
                    let r = PiHatResult { z0, z1, z2, pi_hat: pi, n_valid, ibs1, ibs0 };
                    (pi, r)
                })
                .collect()
        })
        .collect();
    for (i, row) in rows.iter().enumerate() {
        for (j, &(v, r)) in row.iter().enumerate() {
            out.set(i, j, v);
            results[i * n + j] = r;
            results[j * n + i] = r;
        }
    }
    // Keep the triplet alongside the PI_HAT value for the table writer.
    let aux = AuxCounts {
        z0: results.iter().map(|r| r.z0).collect(),
        z1: results.iter().map(|r| r.z1).collect(),
        z2: results.iter().map(|r| r.z2).collect(),
        ..AuxCounts::default()
    };
    out.aux = Some(aux);
    (out, results)
}

/// plink's `e00,e01,e02,e11,e12`: expected IBS counts per IBD state,
/// accumulated per site and normalized by the site count.
fn expected_ibs_coefficients(st: &Standardized) -> (f64, f64, f64, f64, f64) {
    let m = st.n_sites as f64;
    let mut e00 = 0.0;
    let mut e01 = 0.0;
    let mut e02 = 0.0;
    let mut e11 = 0.0;
    let mut e12 = 0.0;
    for &p in &st.freq {
        let q = 1.0 - p;
        // Expected number of zero-IBS alleles for each IBD state.
        e00 += 2.0 * p * p * q * q;
        e01 += 0.0;
        e02 += 0.0;
        // Expected one-IBS alleles.
        e11 += 2.0 * p * q * (p * p + q * q);
        e12 += 2.0 * p * q;
    }
    if m > 0.0 {
        (e00 / m, e01 / m, e02 / m, e11 / m, e12 / m)
    } else {
        (1.0, 1.0, 1.0, 1.0, 1.0)
    }
}

/// plink's clipping of the IBD triplet (`1.9/plink_calc.c:4386-4416`).
fn clip_ibd(z0: f64, z1: f64, z2: f64) -> (f64, f64, f64) {
    let mut z = [z0, z1, z2];
    if z.iter().any(|v| *v > 1.0) {
        let k = z.iter().position(|v| *v > 1.0).unwrap_or(0);
        z = [0.0; 3];
        z[k] = 1.0;
        return (z[0], z[1], z[2]);
    }
    if z.iter().any(|v| *v < 0.0) {
        let k = z.iter().position(|v| *v < 0.0).unwrap_or(0);
        z[k] = 0.0;
        let others: Vec<usize> = (0..3).filter(|&i| i != k).collect();
        let sum = z[others[0]] + z[others[1]];
        if sum > 0.0 {
            for &i in &others {
                z[i] /= sum;
            }
        }
    }
    (z[0], z[1], z[2])
}

/// vcftools `--relatedness2`: Manichaikul phi (01 §5).
///
/// Reproduces vcftools exactly, including its **marginal** denominator
/// `N_Aa[i] + N_Aa[j]` (not pairwise-complete) and the full symmetric matrix
/// with a 0.5 diagonal.
pub fn phi(st: &Standardized) -> (MatrixResult, AuxCounts) {
    let n = st.n_samples();
    let m = st.n_sites;
    let mut out = MatrixResult::new(n, MatrixKind::Phi, Vec::new());
    let mut aux = AuxCounts {
        n_hethet: vec![0; n * n],
        n_opposite_hom: vec![0; n * n],
        ..AuxCounts::default()
    };

    let mut dose = vec![0i8; n * m];
    for i in 0..n {
        for j in 0..m {
            dose[i * m + j] = if st.is_called(i, j) {
                (st.w.get(i, j) * st.sd[j] + st.mean[j]).round() as i8
            } else {
                -1
            };
        }
    }
    // Marginal heterozygote counts.
    let n_het: Vec<u64> = (0..n)
        .map(|i| (0..m).filter(|&s| dose[i * m + s] == 1).count() as u64)
        .collect();
    aux.n_het = n_het.clone();

    // One full row of (phi, hethet, opposite) per task (doc_v20 R13); phi
    // walks the complete ordered matrix, so no mirror writes are needed.
    let rows: Vec<Vec<(f64, u64, u64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            (0..n)
                .map(|j| {
                    let mut hethet = 0u64;
                    let mut opposite = 0u64;
                    for s in 0..m {
                        let a = dose[i * m + s];
                        let b = dose[j * m + s];
                        if a < 0 || b < 0 {
                            continue;
                        }
                        if a == 1 && b == 1 {
                            hethet += 1;
                        } else if a == b {
                            // Both homozygous and equal: nothing.
                        } else if a != 1 && b != 1 {
                            opposite += 1;
                        }
                    }
                    let denom = (n_het[i] + n_het[j]) as f64;
                    let val = if denom > 0.0 {
                        (hethet as f64 - 2.0 * opposite as f64) / denom
                    } else {
                        f64::NAN
                    };
                    (val, hethet, opposite)
                })
                .collect()
        })
        .collect();
    for (i, row) in rows.iter().enumerate() {
        for (j, &(val, hethet, opposite)) in row.iter().enumerate() {
            out.set(i, j, val);
            aux.n_hethet[i * n + j] = hethet;
            aux.n_opposite_hom[i * n + j] = opposite;
        }
    }
    (out, aux)
}

/// Options for the Yang Ajk estimator.
#[derive(Clone, Copy, Debug, Default)]
pub struct AjkOptions {
    /// Write the diagonal as `1 + ...` (vcftools does).
    pub diagonal: bool,
}

/// vcftools `--relatedness`: Yang's Ajk (01 §6).
///
/// Pairwise-complete denominator; requires `0 < p < 1` per site.
pub fn ajk(st: &Standardized, opts: &AjkOptions) -> MatrixResult {
    let n = st.n_samples();
    let m = st.n_sites;
    let mut out = MatrixResult::new(n, MatrixKind::Ajk, Vec::new());

    let mut dose = vec![0i8; n * m];
    for i in 0..n {
        for j in 0..m {
            dose[i * m + j] = if st.is_called(i, j) {
                (st.w.get(i, j) * st.sd[j] + st.mean[j]).round() as i8
            } else {
                -1
            };
        }
    }

    // vcftools skips sites with p at the boundary.
    let usable: Vec<bool> = st
        .freq
        .iter()
        .map(|&p| p > f64::EPSILON && p < 1.0 - f64::EPSILON)
        .collect();

    // One triangle row per task (doc_v20 R13).
    let rows: Vec<Vec<f64>> = (0..n)
        .into_par_iter()
        .map(|i| {
            (0..=i)
                .map(|j| {
                    let mut acc = 0.0f64;
                    let mut sites = 0u64;
                    for s in 0..m {
                        if !usable[s] {
                            continue;
                        }
                        let a = dose[i * m + s];
                        let b = dose[j * m + s];
                        if a < 0 || b < 0 {
                            continue;
                        }
                        let p = st.freq[s];
                        let div = 1.0 / (2.0 * p * (1.0 - p));
                        if i == j {
                            if opts.diagonal {
                                let x = a as f64;
                                acc += (x * x - (1.0 + 2.0 * p) * x + 2.0 * p * p) * div;
                            }
                        } else {
                            acc += (a as f64 - 2.0 * p) * (b as f64 - 2.0 * p) * div;
                        }
                        sites += 1;
                    }
                    if sites == 0 {
                        f64::NAN
                    } else if i == j && opts.diagonal {
                        1.0 + acc / sites as f64
                    } else {
                        acc / sites as f64
                    }
                })
                .collect()
        })
        .collect();
    for (i, row) in rows.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            out.set(i, j, v);
        }
    }
    out
}

/// How bcftools `gtcheck` scores a pair.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DiscordanceMode {
    /// Count mismatching genotypes (`gtcheck -E 0`).
    Count,
    /// Sum of `min_G[-log P_i(G) - log P_j(G)]` (gtcheck's default).
    Probability,
}

/// Options for the discordance estimator.
#[derive(Clone, Copy, Debug)]
pub struct DiscordanceOptions {
    pub mode: DiscordanceMode,
    /// Phred-scaled genotyping error probability (bcftools `-E`, default 40).
    pub error_probability: u32,
}

impl Default for DiscordanceOptions {
    fn default() -> Self {
        Self { mode: DiscordanceMode::Probability, error_probability: 40 }
    }
}

/// bcftools `gtcheck` discordance (01 §7).
///
/// Returns the score matrix plus the per-pair `n_compared` and `n_matching`
/// counts bcftools prints alongside it.
pub fn discordance(
    st: &Standardized,
    opts: &DiscordanceOptions,
) -> (MatrixResult, AuxCounts) {
    let n = st.n_samples();
    let m = st.n_sites;
    let mut out = MatrixResult::new(n, MatrixKind::Discordance, Vec::new());
    let mut aux = AuxCounts {
        n_compared: vec![0; n * n],
        n_matching: vec![0; n * n],
        ..AuxCounts::default()
    };

    let mut dose = vec![0i8; n * m];
    for i in 0..n {
        for j in 0..m {
            dose[i * m + j] = if st.is_called(i, j) {
                (st.w.get(i, j) * st.sd[j] + st.mean[j]).round() as i8
            } else {
                -1
            };
        }
    }

    // bcftools' GT error model: e = 10^(-E/10);
    // P(observed | true) has a 1 on the diagonal and e^|delta| off it.
    let e = 10f64.powf(-0.1 * f64::from(opts.error_probability));
    let neg_log_e = -e.ln();
    let gt_cost = |observed: i8, truth: i8| -> f64 {
        let delta = (i16::from(observed) - i16::from(truth)).unsigned_abs();
        match delta {
            0 => 0.0,
            1 => neg_log_e,
            _ => 2.0 * neg_log_e,
        }
    };

    // One triangle row of (score, compared, matching) per task (doc_v20
    // R13).
    let rows: Vec<Vec<(f64, u64, u64)>> = (0..n)
        .into_par_iter()
        .map(|i| {
            (0..=i)
                .map(|j| {
                    let mut count_mismatch = 0u64;
                    let mut compared = 0u64;
                    let mut matching = 0u64;
                    let mut prob_sum = 0.0f64;
                    for s in 0..m {
                        let a = dose[i * m + s];
                        let b = dose[j * m + s];
                        if a < 0 || b < 0 {
                            continue;
                        }
                        compared += 1;
                        if a == b {
                            matching += 1;
                        } else {
                            count_mismatch += 1;
                        }
                        if opts.mode == DiscordanceMode::Probability {
                            // min over the true genotype G of the summed costs.
                            let mut best = f64::INFINITY;
                            for g in 0..3i8 {
                                let cost = gt_cost(a, g) + gt_cost(b, g);
                                if cost < best {
                                    best = cost;
                                }
                            }
                            prob_sum += best;
                        }
                    }
                    let val = match opts.mode {
                        DiscordanceMode::Count => count_mismatch as f64,
                        DiscordanceMode::Probability => prob_sum,
                    };
                    (val, compared, matching)
                })
                .collect()
        })
        .collect();
    for (i, row) in rows.iter().enumerate() {
        for (j, &(val, compared, matching)) in row.iter().enumerate() {
            out.set(i, j, val);
            aux.n_compared[i * n + j] = compared;
            aux.n_matching[i * n + j] = matching;
            aux.n_compared[j * n + i] = compared;
            aux.n_matching[j * n + i] = matching;
        }
    }
    (out, aux)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::standard::standardize_binom2;
    use reed_core::linalg::Mat;

    /// A fixture with a known pair structure: samples 0-2 are similar to
    /// each other, 3-5 to each other.
    fn fixture() -> Mat {
        Mat::from_fn(6, 8, |i, j| {
            let table = [
                [0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 0.0, 1.0],
                [0.0, 1.0, 2.0, 0.0, 2.0, 2.0, 1.0, 1.0],
                [1.0, 0.0, 0.0, 1.0, 2.0, 1.0, 0.0, 2.0],
                [2.0, 2.0, 1.0, 1.0, 0.0, 0.0, 2.0, 0.0],
                [2.0, 1.0, 2.0, 2.0, 0.0, 0.0, 1.0, 1.0],
                [1.0, 2.0, 0.0, 1.0, 0.0, 1.0, 2.0, 2.0],
            ];
            table[i][j]
        })
    }

    #[test]
    fn allele_ct_distance_matches_manual_sum() {
        let x = fixture();
        let st = standardize_binom2(&x);
        let d = allele_ct_distance(&st, &DistanceOptions { hwe_missing_correction: false });
        // No missing data, so the value is the plain sum of dosage differences.
        for i in 0..6 {
            for j in 0..6 {
                let mut want = 0.0;
                for s in 0..8 {
                    want += (x.get(i, s) - x.get(j, s)).abs();
                }
                assert!(
                    (d.get(i, j) - want).abs() < 1e-9,
                    "({i},{j}): {} vs {want}",
                    d.get(i, j)
                );
            }
        }
    }

    #[test]
    fn ibs_is_one_minus_scaled_distance() {
        let st = standardize_binom2(&fixture());
        let opts = DistanceOptions::default();
        let d = allele_ct_distance(&st, &opts);
        let b = ibs(&st, &opts);
        let o = one_minus_ibs(&st, &opts);
        for k in 0..d.values.len() {
            assert!((b.values[k] + o.values[k] - 1.0).abs() < 1e-12);
            assert!((o.values[k] - d.values[k] * 0.5 / 8.0).abs() < 1e-12);
        }
    }

    #[test]
    fn grm_diagonal_reflects_the_theoretical_variance() {
        // binom2 divides by the *theoretical* sd sqrt(2p(1-p)), not the
        // empirical one, so the diagonal is sum(w^2)/m which equals the ratio
        // of empirical to theoretical variance. It is 1 only when every site
        // is exactly in HWE; otherwise it deviates. plink2 behaves the same
        // way (its GRM diagonal is not forced to 1).
        let st = standardize_binom2(&fixture());
        let g = grm(&st, MissingPolicy::Pairwise, true);
        for i in 0..6 {
            let mut want = 0.0;
            for s in 0..st.n_sites {
                want += st.w.get(i, s) * st.w.get(i, s);
            }
            want /= st.n_sites as f64;
            let got = g.get(i, i);
            assert!((got - want).abs() < 1e-12, "diagonal {i}: {got} vs {want}");
            assert!(got > 0.0 && got < 5.0, "diagonal {i} implausible: {got}");
        }
    }

    #[test]
    fn grm_pairwise_and_meanimpute_agree_without_missing() {
        let st = standardize_binom2(&fixture());
        let a = grm(&st, MissingPolicy::Pairwise, true);
        let b = grm(&st, MissingPolicy::MeanImpute, true);
        for k in 0..a.values.len() {
            assert!((a.values[k] - b.values[k]).abs() < 1e-12);
        }
    }

    #[test]
    fn r9_fixed_alt_column_standardizes_to_zero_not_the_mean() {
        // doc_v20 R9: an all-ALT column (every sample 1/1, dosage 2, p = 1)
        // has theoretical sd 0. The old code wrote the raw mean dosage (2.0)
        // into `w` there — 2 in every cell — which added `4k/m` to every
        // GRM entry (off-diagonals included) and a rank-1 constant direction
        // to PCA. GCTA and plink2 write the centred value 0. An all-REF
        // column cannot tell the two conventions apart (its mean is 0), which
        // is why the all-REF pinning test below never caught this.
        let with_alt = Mat::from_fn(6, 5, |i, j| {
            let table = [
                [0.0, 2.0, 1.0, 2.0, 0.0],
                [1.0, 2.0, 0.0, 2.0, 1.0],
                [2.0, 2.0, 1.0, 0.0, 2.0],
                [0.0, 2.0, 2.0, 1.0, 0.0],
                [2.0, 2.0, 0.0, 0.0, 1.0],
                [1.0, 2.0, 2.0, 2.0, 2.0],
            ];
            table[i][j]
        });
        let st = standardize_binom2(&with_alt);
        assert_eq!(st.sd[1], 0.0, "column 1 is all-ALT, so sd must be 0");
        for i in 0..6 {
            assert_eq!(
                st.w.get(i, 1),
                0.0,
                "fixed-ALT column must hold the centred 0, got {}",
                st.w.get(i, 1)
            );
        }
        // A fixed-ALT column must hit the GRM exactly like a fixed-REF
        // column (both contribute 0 to every numerator) instead of
        // inflating it; the two matrices are bit-identical.
        let with_ref = Mat::from_fn(6, 5, |i, j| {
            let v = with_alt.get(i, j);
            if j == 1 { 0.0 } else { v }
        });
        let a = grm(&standardize_binom2(&with_alt), MissingPolicy::Pairwise, true);
        let b = grm(&standardize_binom2(&with_ref), MissingPolicy::Pairwise, true);
        assert_eq!(a.values, b.values, "fixed-ALT column must not inflate the GRM");
    }

    #[test]
    fn monomorphic_all_ref_column_standardizes_to_zero() {
        // A dedicated matrix with one all-0 column: p = 0 so sd = 0. Since
        // doc_v20 R9 the degenerate branch writes the centred 0 (the
        // GCTA/plink2 convention); for an all-REF column that coincides with
        // the column mean, which is why the pre-R9 mean-fill never showed
        // here. The fixed-ALT case, where the two conventions differ, is
        // pinned by `r9_fixed_alt_column_standardizes_to_zero_not_the_mean`.
        let mono = Mat::from_fn(6, 5, |i, j| {
            // Column 1 is monomorphic at REF (all 0): p = 0, so the
            // theoretical sd sqrt(2p(1-p)) is 0. (An all-1 column has p = 0.5
            // and therefore sd = 0.707 — monomorphic at the *dosage* level is
            // not the same as zero theoretical variance.)
            let table = [
                [0.0, 0.0, 1.0, 2.0, 0.0],
                [1.0, 0.0, 0.0, 2.0, 1.0],
                [2.0, 0.0, 1.0, 0.0, 2.0],
                [0.0, 0.0, 2.0, 1.0, 0.0],
                [2.0, 0.0, 0.0, 0.0, 1.0],
                [1.0, 0.0, 2.0, 2.0, 2.0],
            ];
            table[i][j]
        });
        let st = standardize_binom2(&mono);
        assert_eq!(st.sd[1], 0.0, "column 1 should be monomorphic");
        assert_eq!(st.w.get(0, 1), 0.0, "monomorphic column is centred to 0");
        let g = grm(&st, MissingPolicy::Pairwise, true);
        // The monomorphic column contributes mean^2/m to the diagonal.
        let expected = (0..st.n_sites)
            .map(|s| st.w.get(0, s) * st.w.get(0, s))
            .sum::<f64>()
            / st.n_sites as f64;
        assert!((g.get(0, 0) - expected).abs() < 1e-12);
    }

    #[test]
    fn covariance_differs_from_grm_by_the_standardization() {
        let st = standardize_binom2(&fixture());
        let g = grm(&st, MissingPolicy::Pairwise, true);
        let c = grm(&st, MissingPolicy::Pairwise, false);
        let mut differs = false;
        for k in 0..g.values.len() {
            if (g.values[k] - c.values[k]).abs() > 1e-6 {
                differs = true;
            }
        }
        assert!(differs, "cov must not equal the variance-standardized GRM");
    }

    #[test]
    fn king_formula_matches_the_reference_expression() {
        let c = KingCounts { ibs0: 2, hethet: 5, het1hom2: 3, het2hom1: 3 };
        let want = 0.5 - (4.0 * 2.0 + 3.0 + 3.0) / (4.0 * (5.0 + 3.0));
        assert!((c.kinship() - want).abs() < 1e-12);
    }

    #[test]
    fn king_separates_similar_from_different_pairs() {
        let st = standardize_binom2(&fixture());
        let (k, _) = king(&st);
        let within = k.get(0, 1);
        let across = k.get(0, 3);
        assert!(
            within > across,
            "same-group kinship {within} should exceed cross-group {across}"
        );
    }

    #[test]
    fn pi_hat_is_bounded_and_sums_to_one() {
        let st = standardize_binom2(&fixture());
        let (p, res) = pi_hat(&st, &PiHatOptions::default());
        for i in 0..6 {
            for j in 0..6 {
                let r = &res[i * 6 + j];
                let sum = r.z0 + r.z1 + r.z2;
                assert!(
                    (sum - 1.0).abs() < 1e-6,
                    "Z0+Z1+Z2 = {sum} at ({i},{j})"
                );
                let v = p.get(i, j);
                assert!((0.0..=1.0).contains(&v), "PI_HAT {v} out of range");
            }
        }
    }

    #[test]
    fn phi_diagonal_is_half_as_in_vcftools() {
        let st = standardize_binom2(&fixture());
        let (m, _) = phi(&st);
        for i in 0..6 {
            assert!(
                (m.get(i, i) - 0.5).abs() < 1e-12,
                "diagonal {i} = {}",
                m.get(i, i)
            );
        }
    }

    #[test]
    fn phi_is_symmetric_and_matches_manual_counts() {
        let st = standardize_binom2(&fixture());
        let (m, aux) = phi(&st);
        for i in 0..6 {
            for j in 0..6 {
                assert!((m.get(i, j) - m.get(j, i)).abs() < 1e-12);
            }
        }
        // Recount pair (0,1) directly from the fixture.
        let x = fixture();
        let mut hethet = 0u64;
        let mut opposite = 0u64;
        for s in 0..8 {
            let a = x.get(0, s);
            let b = x.get(1, s);
            if a == 1.0 && b == 1.0 {
                hethet += 1;
            } else if a != 1.0 && b != 1.0 && a != b {
                opposite += 1;
            }
        }
        assert_eq!(aux.n_hethet[1], hethet);
        assert_eq!(aux.n_opposite_hom[1], opposite);
        let n_het0 = (0..8).filter(|&s| x.get(0, s) == 1.0).count() as f64;
        let n_het1 = (0..8).filter(|&s| x.get(1, s) == 1.0).count() as f64;
        let want = (hethet as f64 - 2.0 * opposite as f64) / (n_het0 + n_het1);
        assert!((m.get(0, 1) - want).abs() < 1e-12);
    }

    #[test]
    fn ajk_matches_manual_computation() {
        let x = fixture();
        let st = standardize_binom2(&x);
        let a = ajk(&st, &AjkOptions::default());
        // Recompute pair (0,3) by hand from the definition.
        let mut acc = 0.0;
        let mut sites = 0;
        for s in 0..8 {
            let p = st.freq[s];
            if p <= f64::EPSILON || p >= 1.0 - f64::EPSILON {
                continue;
            }
            let xi = x.get(0, s);
            let xj = x.get(3, s);
            acc += (xi - 2.0 * p) * (xj - 2.0 * p) / (2.0 * p * (1.0 - p));
            sites += 1;
        }
        let want = acc / sites as f64;
        assert!((a.get(0, 3) - want).abs() < 1e-9, "{} vs {want}", a.get(0, 3));
    }

    #[test]
    fn discordance_count_matches_genotype_differences() {
        let x = fixture();
        let st = standardize_binom2(&x);
        let opts = DiscordanceOptions { mode: DiscordanceMode::Count, error_probability: 0 };
        let (d, aux) = discordance(&st, &opts);
        for i in 0..6 {
            for j in 0..6 {
                let mut want = 0.0;
                for s in 0..8 {
                    if x.get(i, s) != x.get(j, s) {
                        want += 1.0;
                    }
                }
                assert!((d.get(i, j) - want).abs() < 1e-12);
                assert_eq!(aux.n_compared[i * 6 + j], 8);
            }
        }
    }

    #[test]
    fn discordance_probability_is_zero_for_identical_samples() {
        let x = fixture();
        let st = standardize_binom2(&x);
        let opts = DiscordanceOptions::default();
        let (d, _) = discordance(&st, &opts);
        for i in 0..6 {
            assert!(d.get(i, i).abs() < 1e-12, "self discordance {}", d.get(i, i));
        }
    }

    #[test]
    fn r13_estimator_results_thread_invariant() {
        // doc_v20 R13: the estimator pair loops are row-parallel. Every
        // (i,j) entry keeps exactly its serial arithmetic and no float is
        // ever summed across threads, so the matrices must be bit-identical
        // under any pool size. Exercise the hard inputs together: missing
        // data, a degenerate all-ALT column (R9) and plink's HWE missing
        // correction.
        let mut x = fixture();
        x.set(0, 0, f64::NAN);
        x.set(4, 6, f64::NAN);
        for i in 0..6 {
            x.set(i, 5, 2.0); // all-ALT: degenerate column
        }
        let st = standardize_binom2(&x);
        let one = rayon::ThreadPoolBuilder::new().num_threads(1).build().unwrap();
        let four = rayon::ThreadPoolBuilder::new().num_threads(4).build().unwrap();

        let bits = |a: &[f64], b: &[f64], what: &str| {
            assert_eq!(a.len(), b.len(), "{what} length");
            for (x, y) in a.iter().zip(b) {
                assert_eq!(x.to_bits(), y.to_bits(), "{what}: {x} vs {y}");
            }
        };

        let g1 = one.install(|| grm(&st, MissingPolicy::Pairwise, true));
        let g4 = four.install(|| grm(&st, MissingPolicy::Pairwise, true));
        bits(&g1.values, &g4.values, "grm");

        let d1 = one.install(|| allele_ct_distance(&st, &DistanceOptions::default()));
        let d4 = four.install(|| allele_ct_distance(&st, &DistanceOptions::default()));
        bits(&d1.values, &d4.values, "allele_ct_distance");

        let (k1, c1) = one.install(|| king(&st));
        let (k4, c4) = four.install(|| king(&st));
        bits(&k1.values, &k4.values, "king");
        assert_eq!(c1, c4, "king counts");
    }

    #[test]
    fn missing_data_is_skipped_everywhere() {
        let mut x = fixture();
        x.set(0, 0, f64::NAN);
        let st = standardize_binom2(&x);
        assert!(!st.is_called(0, 0));
        assert_eq!(st.called_both(0, 1), 7);
        // The pairwise denominator shrinks with the missing site, so compare
        // the diagonal against the fully-called case rather than against 1.
        let full = standardize_binom2(&fixture());
        let g_missing = grm(&st, MissingPolicy::Pairwise, true);
        let g_full = grm(&full, MissingPolicy::Pairwise, true);
        assert!(
            (g_missing.get(0, 0) - g_full.get(0, 0)).abs() > 1e-9,
            "a missing site must change the diagonal"
        );
    }
}
