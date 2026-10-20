//! Standardization of genotype dosages, shared by every matrix estimator.
//!
//! This module is the migration target of `reed-reduce::pca`'s
//! `Standardized`/`Stand`/`Grm` (doc_v17): the PCA path and the distance
//! matrices must standardize identically or their numbers diverge.

use reed_core::linalg::Mat;

/// Standardization convention.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Stand {
    /// `sd = sqrt(2 p (1-p))`. GCTA, plink2 and flashpca's default.
    Binom2,
    /// `sd = sqrt(p (1-p))`. The older EIGENSOFT (Price 2006) convention.
    /// Differs from `Binom2` by a constant factor sqrt(2) per SNP, so the
    /// eigenvectors are identical and the eigenvalues are exactly halved.
    Binom,
}

impl Stand {
    /// Parses a CLI spelling.
    ///
    /// # Errors
    ///
    /// An unknown name.
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "binom2" => Ok(Stand::Binom2),
            "binom" => Ok(Stand::Binom),
            other => Err(format!(
                "unknown standardization {other:?} (binom2|binom)"
            )),
        }
    }

    /// The CLI spelling.
    pub fn name(self) -> &'static str {
        match self {
            Stand::Binom2 => "binom2",
            Stand::Binom => "binom",
        }
    }

    /// The multiplier inside the sd: 2 for binom2, 1 for binom.
    fn mult(self) -> f64 {
        match self {
            Stand::Binom2 => 2.0,
            Stand::Binom => 1.0,
        }
    }
}

/// Standardized genotype matrix plus the per-SNP statistics it used.
#[derive(Clone, Debug)]
pub struct Standardized {
    /// n samples x m sites; missing entries are 0 (the column mean after
    /// centring), i.e. mean imputation.
    pub w: Mat,
    /// Per-site mean dosage (`2p`).
    pub mean: Vec<f64>,
    /// Per-site sd used for the division.
    pub sd: Vec<f64>,
    /// Per-site allele frequency `p = mean/2`.
    pub freq: Vec<f64>,
    /// `n_samples * n_sites`, row-major: whether the call was present.
    ///
    /// The distance and count estimators need "both called" masks; deriving
    /// them from `w` is impossible because a missing call and a call at the
    /// mean both map to 0.
    pub called: Vec<bool>,
    /// Number of sites.
    pub n_sites: usize,
}

impl Standardized {
    /// Number of samples.
    pub fn n_samples(&self) -> usize {
        self.w.rows
    }

    /// Whether sample `i` was called at site `j`.
    #[inline]
    pub fn is_called(&self, i: usize, j: usize) -> bool {
        self.called[i * self.n_sites + j]
    }

    /// Number of called sites for sample `i`.
    pub fn called_count(&self, i: usize) -> usize {
        (0..self.n_sites).filter(|&j| self.is_called(i, j)).count()
    }

    /// Number of sites where both `i` and `j` are called.
    pub fn called_both(&self, i: usize, j: usize) -> usize {
        (0..self.n_sites)
            .filter(|&s| self.is_called(i, s) && self.is_called(j, s))
            .count()
    }
}

/// Standardizes `x` (samples x sites, dosages in {0,1,2}, `NaN` missing).
///
/// Degenerate columns (`sd <= 1e-9`, i.e. monomorphic) get the column mean as
/// their value, matching flashpca's matrix path. GCTA and plink2 write 0
/// there instead; doc_v16 §3.1 records that divergence.
pub fn standardize(x: &Mat, stand: Stand) -> Standardized {
    let n = x.rows;
    let m = x.cols;
    let mut w = Mat::zeros(n, m);
    let mut mean = vec![0.0; m];
    let mut sd = vec![0.0; m];
    let mut freq = vec![0.0; m];
    let mut called = vec![false; n * m];
    let mult = stand.mult();

    for j in 0..m {
        let mut sum = 0.0;
        let mut cnt = 0usize;
        for i in 0..n {
            let v = x.get(i, j);
            if !v.is_nan() {
                sum += v;
                cnt += 1;
                called[i * m + j] = true;
            }
        }
        if cnt == 0 {
            continue;
        }
        let mu = sum / cnt as f64;
        let r = mu / 2.0;
        let s = (mult * r * (1.0 - r)).sqrt();
        mean[j] = mu;
        sd[j] = s;
        freq[j] = r;
        for i in 0..n {
            let v = x.get(i, j);
            let z = if v.is_nan() {
                0.0
            } else if s > 1e-9 {
                (v - mu) / s
            } else {
                mu
            };
            w.set(i, j, z);
        }
    }
    Standardized { w, mean, sd, freq, called, n_sites: m }
}

/// binom2 standardization, the common case.
pub fn standardize_binom2(x: &Mat) -> Standardized {
    standardize(x, Stand::Binom2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toy() -> Mat {
        Mat::from_fn(5, 4, |i, j| {
            let table = [
                [0.0, 1.0, 2.0, 0.0],
                [1.0, 2.0, 1.0, 0.0],
                [2.0, 2.0, 0.0, 1.0],
                [0.0, 0.0, 1.0, 2.0],
                [1.0, 1.0, f64::NAN, 2.0],
            ];
            table[i][j]
        })
    }

    #[test]
    fn columns_are_centred_and_called_mask_tracks_missing() {
        let st = standardize_binom2(&toy());
        for j in 0..st.n_sites {
            let mean: f64 =
                (0..st.n_samples()).map(|i| st.w.get(i, j)).sum::<f64>() / st.n_samples() as f64;
            assert!(mean.abs() < 1e-12, "column {j} mean {mean}");
        }
        assert!(!st.is_called(4, 2), "the NaN must be recorded as missing");
        assert!(st.is_called(0, 2));
        assert_eq!(st.called_count(4), 3);
        assert_eq!(st.called_both(4, 0), 3);
    }

    #[test]
    fn binom_is_half_of_binom2() {
        let x = toy();
        let a = standardize(&x, Stand::Binom2);
        let b = standardize(&x, Stand::Binom);
        for (sa, sb) in a.sd.iter().zip(&b.sd) {
            if *sa > 0.0 && *sb > 0.0 {
                let ratio = sa / sb;
                assert!((ratio - 2f64.sqrt()).abs() < 1e-12, "ratio {ratio}");
            }
        }
        for (fa, fb) in a.freq.iter().zip(&b.freq) {
            assert!((fa - fb).abs() < 1e-15);
        }
    }
}
