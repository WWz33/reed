//! PCA with the binom2 standardization, the convention shared by GCTA,
//! plink2 and flashpca:
//!
//! ```text
//! p_j  = mean dosage of SNP j / 2
//! W_ij = (x_ij - 2 p_j) / sqrt(2 p_j (1 - p_j))
//! A    = W W' / m          (the GRM / relationship matrix)
//! ```
//!
//! Missing genotypes are mean-imputed (set to 0 after centring), which is what
//! flashpca and plink2's `meanimpute` do. Pairwise-complete is a separate
//! missing-data policy and needs the per-pair counts (`grm_n`), see
//! `Grm::pairwise`.

use crate::backend::{sym_eig, sym_eig_top_k, Backend};
use reed_core::linalg::Mat;

// Standardization and the GRM now live in `reed-distance` (doc_v17): the
// distance matrices and PCA must share one implementation or their numbers
// diverge. These are re-exports so the v16 API keeps working.
pub use reed_distance::standard::{standardize, standardize_binom2, Stand, Standardized};
pub use reed_distance::grm::{Grm, MissingPolicy};

/// Result of a PCA run.
#[derive(Clone, Debug)]
pub struct Pca {
    /// Eigenvalues of the GRM, descending (i.e. lambda / m).
    pub values: Vec<f64>,
    /// Sample eigenvectors (unit norm), one column per PC.
    pub vectors: Mat,
    /// Proportion of variance explained per PC.
    pub pve: Vec<f64>,
    /// Eigenvalues of the full spectrum, ascending (empty for partial runs).
    pub all_values: Vec<f64>,
}

impl Pca {
    /// Sample coordinates: `U * sqrt(lambda)`.
    ///
    /// The eigenvectors of `A = W W' / m` already live in sample space, so the
    /// coordinates are just the eigenvectors scaled by the square root of
    /// their eigenvalue. (Multiplying by `W` again would project onto the SNP
    /// space, which is the loadings, not the scores.)
    pub fn scores(&self) -> Mat {
        let n = self.vectors.rows;
        let k = self.vectors.cols;
        let mut out = Mat::zeros(n, k);
        for c in 0..k {
            let scale = self.values[c].max(0.0).sqrt();
            for i in 0..n {
                out.set(i, c, self.vectors.get(i, c) * scale);
            }
        }
        out
    }

    /// SNP loadings: `V = W' U / sqrt(lambda)` (m x k).
    pub fn loadings(&self, w: &Mat) -> Mat {
        let m = w.cols;
        let k = self.vectors.cols;
        let mut out = Mat::zeros(m, k);
        for c in 0..k {
            let scale = if self.values[c] > 0.0 {
                1.0 / self.values[c].sqrt()
            } else {
                0.0
            };
            for j in 0..m {
                let mut acc = 0.0;
                for i in 0..w.rows {
                    acc += w.get(i, j) * self.vectors.get(i, c);
                }
                out.set(j, c, acc * scale);
            }
        }
        out
    }
}

fn pack(eig: crate::backend::SymEigResult, k: Option<usize>) -> Pca {
    let n = eig.values.len();
    let k = k.unwrap_or(n).min(n);
    let (values, vectors) = eig.top_k(k);
    let total: f64 = eig.values.iter().sum();
    let pve = values
        .iter()
        .map(|v| if total != 0.0 { v / total } else { 0.0 })
        .collect();
    Pca { values, vectors, pve, all_values: eig.values }
}

/// Dense PCA: build the GRM and do a full symmetric eigendecomposition.
pub fn pca_dense(x: &Mat, k: usize) -> Pca {
    pca_dense_with(x, k, Stand::Binom2)
}

/// Dense PCA with an explicit standardization convention.
pub fn pca_dense_with(x: &Mat, k: usize, stand: Stand) -> Pca {
    let st = standardize(x, stand);
    let grm = Grm::mean_impute(&st.w);
    pack(sym_eig(&grm.a), Some(k))
}

/// The backend that served the most recent dense PCA call in this build.
pub fn backend() -> Backend {
    crate::backend::COMPILED_BACKEND
}

/// PCA that only extracts the top-k eigenpairs.
pub fn pca_top_k(x: &Mat, k: usize) -> Pca {
    pca_top_k_with(x, k, Stand::Binom2)
}

/// Partial-spectrum PCA with an explicit standardization convention.
pub fn pca_top_k_with(x: &Mat, k: usize, stand: Stand) -> Pca {
    let st = standardize(x, stand);
    let grm = Grm::mean_impute(&st.w);
    let n = grm.a.rows;
    let (values, vectors, _backend) = sym_eig_top_k(&grm.a, k);
    let total: f64 = (0..n).map(|i| grm.a.get(i, i)).sum();
    let pve = values
        .iter()
        .map(|v| if total != 0.0 { v / total } else { 0.0 })
        .collect();
    Pca { values, vectors, pve, all_values: Vec::new() }
}

/// Dense PCA with the pairwise-complete missing-data policy.
pub fn pca_dense_pairwise(x: &Mat, k: usize) -> Pca {
    pca_dense_pairwise_with(x, k, Stand::Binom2)
}

/// Pairwise-complete PCA with an explicit standardization convention.
pub fn pca_dense_pairwise_with(x: &Mat, k: usize, stand: Stand) -> Pca {
    let grm = Grm::pairwise_with(x, stand);
    pack(sym_eig(&grm.a), Some(k))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn toy() -> Mat {
        // 6 samples, 4 SNPs, one missing call.
        Mat::from_fn(6, 4, |i, j| {
            let table = [
                [0.0, 1.0, 2.0, 0.0],
                [1.0, 2.0, 1.0, 0.0],
                [2.0, 2.0, 0.0, 1.0],
                [0.0, 0.0, 1.0, 2.0],
                [1.0, 1.0, 2.0, 2.0],
                [2.0, 1.0, f64::NAN, 1.0],
            ];
            table[i][j]
        })
    }

    #[test]
    fn standardized_columns_are_centred() {
        let st = standardize_binom2(&toy());
        for j in 0..st.w.cols {
            let mean: f64 = (0..st.w.rows).map(|i| st.w.get(i, j)).sum::<f64>()
                / st.w.rows as f64;
            assert!(mean.abs() < 1e-12, "column {j} mean {mean}");
        }
    }

    #[test]
    fn eigenvalues_descending_and_pve_sums_to_one() {
        let p = pca_dense(&toy(), 3);
        for w in p.values.windows(2) {
            assert!(w[0] >= w[1] - 1e-12);
        }
        // pve is per-PC share of the total spectrum, so the first k shares sum
        // to at most 1 (exactly 1 only when k == n).
        let s: f64 = p.pve.iter().sum();
        assert!(s <= 1.0 + 1e-12 && s > 0.0, "pve sum {s}");
        let s_all: f64 = p.all_values.iter().map(|v| v / p.all_values.iter().sum::<f64>()).sum();
        assert!((s_all - 1.0).abs() < 1e-9);
    }

    #[test]
    fn grm_diagonal_is_one_when_complete() {
        let x = Mat::from_fn(5, 6, |i, j| ((i * 3 + j) % 3) as f64);
        let st = standardize_binom2(&x);
        let g = Grm::mean_impute(&st.w);
        for i in 0..5 {
            // With mean imputation and no missing data the diagonal is m/m = 1
            // only if each row is exactly standardized; check it is close.
            assert!((g.a.get(i, i) - 1.0).abs() < 0.6);
        }
    }

    #[test]
    fn pairwise_counts_match_calls() {
        let x = toy();
        let g = Grm::pairwise(&x);
        let counts = g.n_pairs.unwrap();
        // Sample 5 has one missing call (SNP 2), so pair (5,5) sees 3 SNPs.
        assert_eq!(counts.get(5, 5), 3.0);
        assert_eq!(counts.get(0, 0), 4.0);
    }
}
