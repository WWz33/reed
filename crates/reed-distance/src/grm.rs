//! The relationship matrix (GRM), migrated from `reed-reduce::pca` (doc_v17).
//!
//! PCA and the distance matrices must build the relationship matrix the same
//! way; keeping two implementations is exactly the "multiple implementations
//! drift apart" failure mode doc_v15 documented. This module is the single
//! implementation; `reed-reduce::pca` re-exports it.

use reed_core::linalg::Mat;

use crate::standard::Stand;

/// Missing-data policy for the relationship matrix.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MissingPolicy {
    /// Each pair divided by the number of jointly called sites
    /// (plink/GCTA default).
    Pairwise,
    /// Every entry divided by the total site count (plink `meanimpute`).
    MeanImpute,
}

/// The relationship matrix and, for the pairwise policy, the per-pair counts
/// that make it expressible.
#[derive(Clone, Debug)]
pub struct Grm {
    /// n x n relationship matrix.
    pub a: Mat,
    /// n x n counts of jointly non-missing sites (only with [`Grm::pairwise`]).
    pub n_pairs: Option<Mat>,
    /// Number of sites used.
    pub m: usize,
}

impl Grm {
    /// Mean-imputed GRM: `W W' / m`.
    pub fn mean_impute(w: &Mat) -> Self {
        let m = w.cols;
        let n = w.rows;
        let mut a = Mat::zeros(n, n);
        for i in 0..n {
            for j in 0..=i {
                let mut acc = 0.0;
                for s in 0..m {
                    acc += w.get(i, s) * w.get(j, s);
                }
                let val = if m > 0 { acc / m as f64 } else { 0.0 };
                a.set(i, j, val);
                a.set(j, i, val);
            }
        }
        Self { a, n_pairs: None, m }
    }

    /// Pairwise-complete GRM with the default (binom2) standardization.
    pub fn pairwise(x: &Mat) -> Self {
        Self::pairwise_with(x, Stand::Binom2)
    }

    /// Pairwise-complete GRM: each entry divided by the number of sites where
    /// both individuals are called (GCTA's `.grm.N` semantics).
    pub fn pairwise_with(x: &Mat, stand: Stand) -> Self {
        let n = x.rows;
        let m = x.cols;
        let st = crate::standard::standardize(x, stand);
        let mut a = Mat::zeros(n, n);
        let mut counts = Mat::zeros(n, n);
        for i in 0..n {
            for j in 0..=i {
                let mut acc = 0.0;
                let mut cnt = 0.0;
                for s in 0..m {
                    if st.is_called(i, s) && st.is_called(j, s) {
                        acc += st.w.get(i, s) * st.w.get(j, s);
                        cnt += 1.0;
                    }
                }
                let val = if cnt > 0.0 { acc / cnt } else { 0.0 };
                a.set(i, j, val);
                a.set(j, i, val);
                counts.set(i, j, cnt);
                counts.set(j, i, cnt);
            }
        }
        Self { a, n_pairs: Some(counts), m }
    }

    /// Builds from an already-standardized matrix.
    pub fn from_standardized(st: &crate::standard::Standardized, policy: MissingPolicy) -> Self {
        let n = st.n_samples();
        let m = st.n_sites;
        let mut a = Mat::zeros(n, n);
        let mut counts = Mat::zeros(n, n);
        for i in 0..n {
            for j in 0..=i {
                let mut acc = 0.0;
                let mut cnt = 0.0;
                for s in 0..m {
                    match policy {
                        MissingPolicy::MeanImpute => {
                            acc += st.w.get(i, s) * st.w.get(j, s);
                            cnt = m as f64;
                        }
                        MissingPolicy::Pairwise => {
                            if st.is_called(i, s) && st.is_called(j, s) {
                                acc += st.w.get(i, s) * st.w.get(j, s);
                                cnt += 1.0;
                            }
                        }
                    }
                }
                let val = if cnt > 0.0 { acc / cnt } else { 0.0 };
                a.set(i, j, val);
                a.set(j, i, val);
                counts.set(i, j, cnt);
                counts.set(j, i, cnt);
            }
        }
        Self {
            a,
            n_pairs: (policy == MissingPolicy::Pairwise).then_some(counts),
            m,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toy() -> Mat {
        Mat::from_fn(5, 6, |i, j| {
            let table = [
                [0.0, 1.0, 2.0, 0.0, 1.0, 2.0],
                [1.0, 2.0, 1.0, 0.0, 0.0, 1.0],
                [2.0, 2.0, 0.0, 1.0, 2.0, 0.0],
                [0.0, 0.0, 1.0, 2.0, 1.0, 1.0],
                [1.0, 1.0, f64::NAN, 2.0, 2.0, 0.0],
            ];
            table[i][j]
        })
    }

    #[test]
    fn pairwise_counts_match_calls() {
        let g = Grm::pairwise(&toy());
        let counts = g.n_pairs.expect("pairwise keeps the counts");
        assert_eq!(counts.get(0, 0), 6.0);
        assert_eq!(counts.get(4, 4), 5.0, "sample 4 has one missing call");
        assert_eq!(counts.get(0, 4), 5.0);
    }

    #[test]
    fn mean_impute_and_pairwise_agree_without_missing() {
        let x = Mat::from_fn(4, 5, |i, j| ((i * 3 + j) % 3) as f64);
        let st = crate::standard::standardize_binom2(&x);
        let a = Grm::mean_impute(&st.w);
        let b = Grm::from_standardized(&st, MissingPolicy::Pairwise);
        for i in 0..4 {
            for j in 0..4 {
                assert!((a.a.get(i, j) - b.a.get(i, j)).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn grm_is_symmetric() {
        let g = Grm::pairwise(&toy());
        for i in 0..g.a.rows {
            for j in 0..g.a.cols {
                assert!((g.a.get(i, j) - g.a.get(j, i)).abs() < 1e-15);
            }
        }
    }
}
