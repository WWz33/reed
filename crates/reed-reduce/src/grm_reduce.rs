//! GRM-driven dimension reduction.
//!
//! The statistical point: for genotype data the meaningful individual-level
//! geometry is the relationship matrix `A = W W' / m` (the GRM), not the raw
//! Euclidean distance between dosage rows. Feeding `A` (or a distance derived
//! from it) into the graph/kernel methods gives embeddings that respect
//! relatedness and population structure.
//!
//! Two derived operators are provided:
//!
//! * **GRM distance**: `d_ij = sqrt(max(0, a_ii + a_jj - 2 a_ij))`, i.e. the
//!   Euclidean distance in the standardized-genotype space, which is what
//!   `A` encodes. This is the correct input for PCoA/Isomap/UMAP.
//! * **GRM similarity**: `s_ij = a_ij / sqrt(a_ii a_jj)` (a correlation-like
//!   normalization), used directly as a kernel.
//!
//! Both are computed from a `Grm`, so the missing-data policy chosen when the
//! GRM was built (mean impute vs pairwise) propagates to every method here.

use crate::graph::{diffusion_map, laplacian_eigenmaps, pcoa, Embedding, Kernel};
use reed_core::linalg::Mat;
use crate::pca::Grm;

/// Distance matrix induced by a GRM: `d_ij^2 = a_ii + a_jj - 2 a_ij`.
pub fn grm_distances(grm: &Grm) -> Mat {
    let n = grm.a.rows;
    let mut d = Mat::zeros(n, n);
    for i in 0..n {
        let aii = grm.a.get(i, i);
        for j in 0..i {
            let a_jj = grm.a.get(j, j);
            let v = (aii + a_jj - 2.0 * grm.a.get(i, j)).max(0.0).sqrt();
            d.set(i, j, v);
            d.set(j, i, v);
        }
    }
    d
}

/// Similarity matrix induced by a GRM: `s_ij = a_ij / sqrt(a_ii a_jj)`.
pub fn grm_similarity(grm: &Grm) -> Mat {
    let n = grm.a.rows;
    let mut s = Mat::zeros(n, n);
    for i in 0..n {
        let di = grm.a.get(i, i).max(1e-12);
        for j in 0..n {
            let dj = grm.a.get(j, j).max(1e-12);
            s.set(i, j, grm.a.get(i, j) / (di * dj).sqrt());
        }
    }
    s
}

/// PCoA on GRM-induced distances (the "PCA of the GRM" view, with the
/// advantage that any distance-shaped operator can be swapped in).
pub fn grm_pcoa(grm: &Grm, k: usize) -> Embedding {
    pcoa(&grm_distances(grm), k)
}

/// Kernel PCA using the normalized GRM as the kernel.
pub fn grm_kernel_pca(grm: &Grm, k: usize) -> Embedding {
    // kernel_pca expects data, so use the similarity matrix directly through
    // the same centring + eigendecomposition path.
    let mut s = grm_similarity(grm);
    crate::graph::centre_kernel(&mut s);
    let e = crate::backend::sym_eig(&s);
    let (vals, vecs) = e.top_k(k);
    let n = s.rows;
    let mut coords = Mat::zeros(n, k);
    for c in 0..k {
        let scale = vals[c].abs().sqrt();
        for i in 0..n {
            coords.set(i, c, vecs.get(i, c) * scale);
        }
    }
    Embedding { coords, values: vals }
}

/// Laplacian eigenmaps on a kNN graph built from GRM distances.
pub fn grm_laplacian(grm: &Grm, k: usize, n_neighbors: usize) -> Embedding {
    let d = grm_distances(grm);
    let w = knn_from_distances(&d, n_neighbors);
    laplacian_eigenmaps(&w, k)
}

/// Diffusion map on a kNN graph built from GRM distances.
pub fn grm_diffusion(grm: &Grm, k: usize, n_neighbors: usize, alpha: f64) -> Embedding {
    let d = grm_distances(grm);
    let w = knn_from_distances(&d, n_neighbors);
    diffusion_map(&w, k, alpha)
}

/// Isomap on GRM distances.
pub fn grm_isomap(grm: &Grm, k: usize, n_neighbors: usize) -> Embedding {
    let d = grm_distances(grm);
    isomap_from_distances(&d, n_neighbors, k)
}

/// UMAP-style graph on GRM distances: returns the fuzzy weights so the caller
/// can run the SGD stage.
pub fn grm_knn_weights(grm: &Grm, n_neighbors: usize) -> Mat {
    let d = grm_distances(grm);
    knn_from_distances(&d, n_neighbors)
}

/// t-SNE on GRM distances.
///
/// # Errors
///
/// [`crate::extra::TsneError`]; see [`crate::extra::tsne`].
pub fn grm_tsne(
    grm: &Grm,
    cfg: &crate::extra::TsneConfig,
) -> Result<reed_core::linalg::Mat, crate::extra::TsneError> {
    let d = grm_distances(grm);
    let mut d2 = reed_core::linalg::Mat::zeros(d.rows, d.cols);
    for i in 0..d.rows {
        for j in 0..d.cols {
            let v = d.get(i, j);
            d2.set(i, j, v * v);
        }
    }
    crate::extra::tsne_from_sq_distances(&d2, cfg)
}

/// UMAP on GRM distances.
pub fn grm_umap(grm: &Grm, cfg: &crate::umap::UmapConfig) -> reed_core::linalg::Mat {
    let d = grm_distances(grm);
    crate::umap::umap_from_distances(&d, cfg)
}

/// LLE with GRM-distance neighbourhoods; `features` supplies the local Gram
/// matrix (pass the standardized genotype matrix).
pub fn grm_lle(
    grm: &Grm,
    features: &reed_core::linalg::Mat,
    k: usize,
    n_components: usize,
) -> reed_core::linalg::Mat {
    let d = grm_distances(grm);
    let mut d2 = reed_core::linalg::Mat::zeros(d.rows, d.cols);
    for i in 0..d.rows {
        for j in 0..d.cols {
            let v = d.get(i, j);
            d2.set(i, j, v * v);
        }
    }
    crate::extra::lle_from_sq_distances(features, &d2, k, n_components)
}

/// Build a Gaussian-weighted kNN graph from a distance matrix.
fn knn_from_distances(d: &Mat, k: usize) -> Mat {
    let n = d.rows;
    let mut w = Mat::zeros(n, n);
    // Per-point kernel width: the distance to the k-th neighbour. A single
    // global width fails when the GRM distances have a narrow range (then the
    // Gaussian collapses to zero for every edge and the graph becomes a set of
    // isolated points, which is exactly what a global median produced here).
    for i in 0..n {
        let mut order: Vec<(f64, usize)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (d.get(i, j), j))
            .collect();
        order.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let sigma = order
            .get(k.min(order.len().saturating_sub(1)))
            .map(|&(v, _)| v)
            .unwrap_or(1.0)
            .max(1e-12);
        let gamma = 1.0 / (sigma * sigma);
        for &(dist, j) in order.iter().take(k) {
            let val = (-gamma * dist * dist).exp();
            if val > w.get(i, j) {
                w.set(i, j, val);
            }
        }
    }
    for i in 0..n {
        for j in 0..i {
            let v = w.get(i, j).max(w.get(j, i));
            w.set(i, j, v);
            w.set(j, i, v);
        }
    }
    w
}

/// Isomap from a precomputed distance matrix (geodesics on the kNN graph).
fn isomap_from_distances(d_in: &Mat, k: usize, n_components: usize) -> Embedding {
    let n = d_in.rows;
    let mut d = d_in.clone();
    let mut inf = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            if i != j {
                inf = inf.max(d.get(i, j));
            }
        }
    }
    let inf = inf * (n as f64 + 1.0);
    let mut keep = vec![vec![false; n]; n];
    for i in 0..n {
        let mut order: Vec<(f64, usize)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (d.get(i, j), j))
            .collect();
        order.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        for &(_, j) in order.iter().take(k) {
            keep[i][j] = true;
            keep[j][i] = true;
        }
    }
    for i in 0..n {
        for j in 0..n {
            if i != j && !keep[i][j] {
                d.set(i, j, inf);
            }
        }
    }
    for kk in 0..n {
        for i in 0..n {
            let dik = d.get(i, kk);
            if dik >= inf {
                continue;
            }
            for j in 0..n {
                let cand = dik + d.get(kk, j);
                if cand < d.get(i, j) {
                    d.set(i, j, cand);
                }
            }
        }
    }
    let mut max_fin = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            let v = d.get(i, j);
            if v < inf && v > max_fin {
                max_fin = v;
            }
        }
    }
    for i in 0..n {
        for j in 0..n {
            if d.get(i, j) >= inf {
                d.set(i, j, max_fin);
            }
        }
    }
    pcoa(&d, n_components)
}

/// Convenience: kernel PCA with a linear kernel equals PCoA on the same
/// distances; kept explicit so callers can pick the RBF/IBS variants too.
pub fn grm_kernel(grm: &Grm, k: usize, kind: Kernel) -> Embedding {
    match kind {
        Kernel::Linear => grm_kernel_pca(grm, k),
        _ => {
            // For non-linear kernels we need features; the GRM distances are
            // the feature-space metric, so use them directly.
            let d = grm_distances(grm);
            let n = d.rows;
            let mut kk = Mat::zeros(n, n);
            for i in 0..n {
                for j in 0..n {
                    let v = match kind {
                        Kernel::Rbf => {
                            let sigma = 1.0f64;
                            (-(d.get(i, j) * d.get(i, j)) / (2.0 * sigma * sigma)).exp()
                        }
                        Kernel::Ibs => 1.0 - d.get(i, j),
                        Kernel::Linear => unreachable!(),
                    };
                    kk.set(i, j, v);
                }
            }
            crate::graph::centre_kernel(&mut kk);
            let e = crate::backend::sym_eig(&kk);
            let (vals, vecs) = e.top_k(k);
            let mut coords = Mat::zeros(n, k);
            for c in 0..k {
                let scale = vals[c].abs().sqrt();
                for i in 0..n {
                    coords.set(i, c, vecs.get(i, c) * scale);
                }
            }
            Embedding { coords, values: vals }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pca::{standardize_binom2, Grm};

    fn geno() -> Mat {
        let mut state = 0x9E3779B97F4A7C15u64;
        let mut rnd = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let n_snp = 200;
        let base: Vec<f64> = (0..n_snp).map(|_| 0.05 + 0.9 * rnd()).collect();
        Mat::from_fn(45, n_snp, move |i, j| {
            let group = i / 15;
            let shift = [0.0, 0.18, -0.15][group];
            let p = (base[j] + shift).clamp(0.001, 0.999);
            let u = rnd();
            if u < p * p {
                2.0
            } else if u < p * p + 2.0 * p * (1.0 - p) {
                1.0
            } else {
                0.0
            }
        })
    }

    fn ratio(m: &Mat, n_per: usize) -> f64 {
        let dim = m.cols;
        let mut cent = vec![vec![0.0f64; dim]; 3];
        for g in 0..3 {
            for i in g * n_per..(g + 1) * n_per {
                for d in 0..dim {
                    cent[g][d] += m.get(i, d) / n_per as f64;
                }
            }
        }
        let mut within = 0.0;
        for g in 0..3 {
            for i in g * n_per..(g + 1) * n_per {
                let mut s = 0.0;
                for d in 0..dim {
                    s += (m.get(i, d) - cent[g][d]).powi(2);
                }
                within += s;
            }
        }
        let within = (within / (3 * n_per) as f64).sqrt().max(1e-12);
        let mut sep = f64::MAX;
        for a in 0..3 {
            for b in (a + 1)..3 {
                let mut s = 0.0;
                for d in 0..dim {
                    s += (cent[a][d] - cent[b][d]).powi(2);
                }
                sep = sep.min(s.sqrt());
            }
        }
        sep / within
    }

    #[test]
    fn grm_distances_are_metric_like() {
        let x = geno();
        let st = standardize_binom2(&x);
        let g = Grm::mean_impute(&st.w);
        let d = grm_distances(&g);
        for i in 0..d.rows {
            assert!(d.get(i, i).abs() < 1e-9, "diagonal not zero");
            for j in 0..d.cols {
                assert!((d.get(i, j) - d.get(j, i)).abs() < 1e-12);
                assert!(d.get(i, j) >= 0.0);
            }
        }
        // Triangle inequality on a few triples.
        for i in 0..5 {
            for j in 0..5 {
                for k in 0..5 {
                    assert!(d.get(i, k) <= d.get(i, j) + d.get(j, k) + 1e-9);
                }
            }
        }
    }

    #[test]
    fn grm_pcoa_separates_populations() {
        let x = geno();
        let st = standardize_binom2(&x);
        let g = Grm::mean_impute(&st.w);
        let e = grm_pcoa(&g, 2);
        assert!(ratio(&e.coords, 15) > 2.0, "grm pcoa ratio {}", ratio(&e.coords, 15));
    }

    #[test]
    fn grm_kernel_and_graph_methods_separate_populations() {
        let x = geno();
        let st = standardize_binom2(&x);
        let g = Grm::mean_impute(&st.w);
        for (name, e) in [
            ("kernel_pca", grm_kernel_pca(&g, 2)),
            ("laplacian", grm_laplacian(&g, 2, 8)),
            ("diffusion", grm_diffusion(&g, 2, 8, 1.0)),
            ("isomap", grm_isomap(&g, 2, 8)),
            ("kernel_rbf", grm_kernel(&g, 2, Kernel::Rbf)),
        ] {
            let r = ratio(&e.coords, 15);
            assert!(r > 1.5, "{name} ratio {r}");
        }
    }

    #[test]
    fn all_distance_based_families_accept_grm_geometry() {
        let x = geno();
        let st = standardize_binom2(&x);
        let g = Grm::mean_impute(&st.w);
        let tsne_cfg = crate::extra::TsneConfig {
            perplexity: 8.0,
            n_iter: 300,
            seed: 3,
            ..Default::default()
        };
        let umap_cfg = crate::umap::UmapConfig {
            n_neighbors: 8,
            n_epochs: 200,
            seed: 7,
            ..Default::default()
        };
        // Deterministic methods: assert real separation.
        let det: Vec<(&str, Mat)> = vec![
            ("lle", grm_lle(&g, &st.w, 8, 2)),
        ];
        for (name, m) in det {
            let r = ratio(&m, 15);
            assert!(r > 1.2, "{name} on GRM geometry: ratio {r}");
        }
        // Stochastic methods: assert they run, produce finite output, and are
        // not degenerate. The exact ratio depends on the SGD trajectory, which
        // legitimately differs between backends (different summation order).
        let stoch: Vec<(&str, Mat)> = vec![
            ("tsne", grm_tsne(&g, &tsne_cfg).expect("grm tsne fixture")),
            ("umap", grm_umap(&g, &umap_cfg)),
        ];
        for (name, m) in stoch {
            assert_eq!(m.rows, st.w.rows, "{name}: wrong sample count");
            assert!(m.cols >= 2, "{name}: need at least 2 components");
            let spread = m
                .data
                .iter()
                .cloned()
                .fold(0.0f64, |a, v| a.max(v.abs()));
            assert!(spread.is_finite() && spread > 1e-6, "{name}: degenerate output");
        }
    }

    #[test]
    fn distance_entry_points_match_feature_entry_points() {
        // Passing the same geometry explicitly must reproduce the
        // feature-based call bit for bit.
        let x = geno();
        let n = x.rows;
        let mut d = Mat::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                d.set(i, j, reed_core::linalg::dist2(x.row(i), x.row(j)).sqrt());
            }
        }
        let cfg = crate::umap::UmapConfig {
            n_neighbors: 6,
            n_epochs: 60,
            seed: 11,
            ..Default::default()
        };
        let a = crate::umap::umap(&x, &cfg);
        let b = crate::umap::umap_from_distances(&d, &cfg);
        let mut worst: f64 = 0.0;
        for i in 0..a.rows {
            for j in 0..a.cols {
                worst = worst.max((a.get(i, j) - b.get(i, j)).abs());
            }
        }
        assert!(worst < 1e-12, "umap feature vs distance paths differ by {worst}");
    }

    #[test]
    fn pairwise_and_mean_impute_agree_without_missing_data() {
        let x = geno();
        let st = standardize_binom2(&x);
        let g1 = Grm::mean_impute(&st.w);
        let g2 = Grm::pairwise(&x);
        for i in 0..g1.a.rows {
            for j in 0..g1.a.cols {
                let diff = (g1.a.get(i, j) - g2.a.get(i, j)).abs();
                assert!(diff < 1e-10, "entry ({i},{j}) differs by {diff}");
            }
        }
    }
}
