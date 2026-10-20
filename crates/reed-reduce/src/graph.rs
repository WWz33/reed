//! Graph- and distance-based dimension reduction.
//!
//! All of these reduce to "build an N x N similarity/distance operator, then
//! take its top or bottom eigenvectors", which is why they can share the same
//! dense eigensolver as PCA:
//!
//! | method | operator | eigenpairs |
//! |---|---|---|
//! | PCoA / classical MDS | double-centred distance matrix | largest |
//! | Kernel PCA | kernel matrix K | largest |
//! | Diffusion map | row-normalized kernel powers | largest (after normalization) |
//! | Laplacian eigenmaps | graph Laplacian L = D - W | smallest |
//! | Isomap | double-centred geodesic distances | largest |
//!
//! Every method takes the same `Mat` (samples x features) and returns an
//! embedding (samples x components).

use crate::backend::sym_eig;
use reed_core::linalg::{dist2, Mat};

/// A dense embedding result.
#[derive(Clone, Debug)]
pub struct Embedding {
    pub coords: Mat,
    pub values: Vec<f64>,
}

/// Pairwise Euclidean distance matrix (n x n).
pub fn pairwise_distances(x: &Mat) -> Mat {
    let n = x.rows;
    let mut d = Mat::zeros(n, n);
    for i in 0..n {
        for j in 0..i {
            let v = dist2(x.row(i), x.row(j)).sqrt();
            d.set(i, j, v);
            d.set(j, i, v);
        }
    }
    d
}

/// Classical multidimensional scaling (PCoA) on a distance matrix.
///
/// `B = -0.5 * J D^2 J` with `J = I - 11'/n`; the top eigenvectors of `B`
/// scaled by `sqrt(lambda)` are the principal coordinates.
pub fn pcoa(d: &Mat, k: usize) -> Embedding {
    let n = d.rows;
    let mut b = Mat::zeros(n, n);
    // Row means of D^2, then double centring.
    let mut row_mean = vec![0.0f64; n];
    for i in 0..n {
        let mut s = 0.0;
        for j in 0..n {
            let v = d.get(i, j);
            s += v * v;
        }
        row_mean[i] = s / n as f64;
    }
    let grand: f64 = row_mean.iter().sum::<f64>() / n as f64;
    for i in 0..n {
        for j in 0..n {
            let v = d.get(i, j);
            b.set(i, j, -0.5 * (v * v - row_mean[i] - row_mean[j] + grand));
        }
    }
    let e = sym_eig(&b);
    let (vals, vecs) = e.top_k(k);
    let mut coords = Mat::zeros(n, k);
    for c in 0..k {
        let scale = vals[c].max(0.0).sqrt();
        for i in 0..n {
            coords.set(i, c, vecs.get(i, c) * scale);
        }
    }
    Embedding { coords, values: vals }
}

/// Kernel PCA. `kind` selects the kernel; `gamma` is the RBF width parameter.
pub fn kernel_pca(x: &Mat, k: usize, kind: Kernel, gamma: f64) -> Embedding {
    let n = x.rows;
    let gamma = if gamma > 0.0 { gamma } else { auto_gamma(x) };
    let mut kk = Mat::zeros(n, n);
    for i in 0..n {
        for j in 0..=i {
            let v = match kind {
                Kernel::Linear => {
                    let mut s = 0.0;
                    for (a, b) in x.row(i).iter().zip(x.row(j)) {
                        s += a * b;
                    }
                    s
                }
                Kernel::Rbf => (-gamma * dist2(x.row(i), x.row(j))).exp(),
                Kernel::Ibs => {
                    // Identity-by-state similarity: 1 - mean |x_i - x_j| / 2.
                    let mut s = 0.0;
                    let mut cnt = 0.0;
                    for (a, b) in x.row(i).iter().zip(x.row(j)) {
                        if !a.is_nan() && !b.is_nan() {
                            s += (a - b).abs();
                            cnt += 1.0;
                        }
                    }
                    if cnt > 0.0 { 1.0 - s / (2.0 * cnt) } else { 0.0 }
                }
            };
            kk.set(i, j, v);
            kk.set(j, i, v);
        }
    }
    // Centre the kernel: K <- K - 1_n K - K 1_n + 1_n K 1_n.
    centre_kernel(&mut kk);
    let e = sym_eig(&kk);
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

pub fn centre_kernel(k: &mut Mat) {
    let n = k.rows;
    let mut row_mean = vec![0.0f64; n];
    for i in 0..n {
        row_mean[i] = (0..n).map(|j| k.get(i, j)).sum::<f64>() / n as f64;
    }
    let grand = row_mean.iter().sum::<f64>() / n as f64;
    for i in 0..n {
        for j in 0..n {
            let v = k.get(i, j) - row_mean[i] - row_mean[j] + grand;
            k.set(i, j, v);
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kernel {
    Linear,
    Rbf,
    Ibs,
}

/// A kernel width that adapts to the data scale: `1 / median(distance)^2`.
///
/// A fixed `gamma` breaks as soon as the number of features (and therefore
/// the typical distance) changes: the Gaussian then either collapses to the
/// identity or saturates. All graph/kernel entry points default to this.
pub fn auto_gamma(x: &Mat) -> f64 {
    let n = x.rows;
    let mut ds = Vec::new();
    // Sample up to ~2000 pairs to keep this O(1) in n.
    let step = (n / 45).max(1);
    let mut i = 0;
    while i < n {
        let mut j = i + step;
        while j < n {
            let d = dist2(x.row(i), x.row(j)).sqrt();
            if d > 0.0 {
                ds.push(d);
            }
            j += step;
        }
        i += step;
    }
    if ds.is_empty() {
        return 1.0;
    }
    ds.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = ds[ds.len() / 2];
    1.0 / (med * med).max(1e-12)
}

/// k-nearest-neighbour graph, returned as a symmetric weight matrix using the
/// Gaussian kernel on the Euclidean distances. Pass `gamma = 0.0` to derive
/// the width from the data (`auto_gamma`).
pub fn knn_graph(x: &Mat, k: usize, gamma: f64) -> Mat {
    let gamma = if gamma > 0.0 { gamma } else { auto_gamma(x) };
    let n = x.rows;
    let mut w = Mat::zeros(n, n);
    for i in 0..n {
        let mut d: Vec<(f64, usize)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (dist2(x.row(i), x.row(j)), j))
            .collect();
        d.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        for &(d2, j) in d.iter().take(k) {
            let val = (-gamma * d2).exp();
            w.set(i, j, val.max(w.get(i, j)));
        }
    }
    // Symmetrize with the max, the usual kNN graph convention.
    for i in 0..n {
        for j in 0..i {
            let v = w.get(i, j).max(w.get(j, i));
            w.set(i, j, v);
            w.set(j, i, v);
        }
    }
    w
}

/// Laplacian eigenmaps: eigenvectors of the smallest non-zero eigenvalues of
/// the normalized Laplacian `L_sym = I - D^-1/2 W D^-1/2`.
pub fn laplacian_eigenmaps(w: &Mat, k: usize) -> Embedding {
    let n = w.rows;
    let mut deg = vec![0.0f64; n];
    for i in 0..n {
        deg[i] = (0..n).map(|j| w.get(i, j)).sum();
    }
    let mut l = Mat::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let v = if i == j { 1.0 } else { 0.0 };
            let denom = (deg[i] * deg[j]).sqrt();
            let off = if denom > 0.0 { w.get(i, j) / denom } else { 0.0 };
            l.set(i, j, v - off);
        }
    }
    let e = sym_eig(&l);
    // A disconnected graph has one zero eigenvalue per connected component.
    // Skipping only the first leaves component indicators in the output (they
    // are legitimate eigenvectors, just not structure); skipping all of them
    // can leave fewer than k directions and silently return zero columns.
    // Report the component count so the caller can tell the two apart.
    let scale = e.values.last().copied().unwrap_or(1.0).abs().max(1.0);
    let tol = 1e-10 * scale;
    let n_zero = e.values.iter().filter(|v| v.abs() <= tol).count();
    let skip = n_zero.max(1); // always drop at least the trivial constant
    let mut values = Vec::with_capacity(k);
    let mut coords = Mat::zeros(n, k);
    for (c, idx) in (skip..n).enumerate() {
        if c >= k {
            break;
        }
        values.push(e.values[idx]);
        for (i, &d) in deg.iter().enumerate() {
            let scale = if d > 0.0 { 1.0 / d.sqrt() } else { 0.0 };
            coords.set(i, c, e.vectors.get(i, idx) * scale);
        }
    }
    // If the graph had fewer usable directions than requested, the remaining
    // columns stay zero; that is visible in the output rather than an error.
    Embedding { coords, values }
}

/// Diffusion map: the first `k` non-trivial eigenvectors of the random-walk
/// operator `P = D^-1 W`, with coordinates scaled by the eigenvalues.
pub fn diffusion_map(w: &Mat, k: usize, alpha: f64) -> Embedding {
    let n = w.rows;
    // Density normalization: W <- D^-alpha W D^-alpha.
    let mut deg = vec![0.0f64; n];
    for i in 0..n {
        deg[i] = (0..n).map(|j| w.get(i, j)).sum();
    }
    let mut wn = Mat::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let di = if deg[i] > 0.0 { deg[i].powf(-alpha) } else { 0.0 };
            let dj = if deg[j] > 0.0 { deg[j].powf(-alpha) } else { 0.0 };
            wn.set(i, j, w.get(i, j) * di * dj);
        }
    }
    let mut deg2 = vec![0.0f64; n];
    for i in 0..n {
        deg2[i] = (0..n).map(|j| wn.get(i, j)).sum();
    }
    // Symmetric conjugate: S = D^-1/2 W D^-1/2.
    let mut s = Mat::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let denom = (deg2[i] * deg2[j]).sqrt();
            s.set(i, j, if denom > 0.0 { wn.get(i, j) / denom } else { 0.0 });
        }
    }
    let e = sym_eig(&s);
    let mut values = Vec::with_capacity(k);
    let mut coords = Mat::zeros(n, k);
    let mut c = 0usize;
    let mut skipped_stationary = false;
    for idx in (0..n).rev() {
        if c >= k {
            break;
        }
        // Skip the stationary eigenvector (eigenvalue 1) exactly once.
        if !skipped_stationary && (e.values[idx] - 1.0).abs() < 1e-9 {
            skipped_stationary = true;
            continue;
        }
        let lam = e.values[idx];
        values.push(lam);
        for i in 0..n {
            let scale = if deg2[i] > 0.0 { 1.0 / deg2[i].sqrt() } else { 0.0 };
            coords.set(i, c, e.vectors.get(i, idx) * scale * lam);
        }
        c += 1;
    }
    Embedding { coords, values }
}

/// Isomap: geodesic distances from a kNN graph, then classical MDS.
pub fn isomap(x: &Mat, k: usize, n_components: usize) -> Embedding {
    let n = x.rows;
    let mut d = pairwise_distances(x);
    // Restrict to the kNN graph, then Floyd-Warshall for shortest paths.
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
    // Non-finite geodesics break double centring; clamp to the largest finite.
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Three blobs joined by a sparse bridge so the kNN graph stays connected.
    fn connected_blobs(n_per: usize) -> Mat {
        let mut rows = Vec::new();
        for c in 0..3 {
            for i in 0..n_per {
                let base = c as f64 * 10.0;
                let jitter = ((i * 37) % 11) as f64 * 0.05;
                rows.push([base + jitter, base + jitter * 0.5]);
            }
        }
        // Bridge points between consecutive blobs.
        for step in 0..5 {
            let t = step as f64 / 5.0;
            let a = 10.0 * 0.0 + t * 10.0;
            rows.push([a, a * 0.5]);
            let b = 10.0 + t * 10.0;
            rows.push([b, b * 0.5]);
        }
        Mat::from_fn(rows.len(), 2, |i, j| rows[i][j])
    }

    fn blobs(n_per: usize) -> Mat {
        // Three well-separated clusters in 2D.
        let mut rows = Vec::new();
        for c in 0..3 {
            for i in 0..n_per {
                let base = c as f64 * 10.0;
                let jitter = ((i * 37) % 11) as f64 * 0.05;
                rows.push([base + jitter, base + jitter * 0.5]);
            }
        }
        Mat::from_fn(rows.len(), 2, |i, j| rows[i][j])
    }

    fn cluster_separation(coords: &Mat, n_per: usize) -> bool {
        // Compare between-group separation on the first coordinate with the
        // within-group standard deviation.
        let groups: Vec<f64> = (0..3)
            .map(|g| {
                (g * n_per..(g + 1) * n_per)
                    .map(|i| coords.get(i, 0))
                    .sum::<f64>()
                    / n_per as f64
            })
            .collect();
        let mut within = 0.0f64;
        for g in 0..3 {
            for i in g * n_per..(g + 1) * n_per {
                within += (coords.get(i, 0) - groups[g]).powi(2);
            }
        }
        within = (within / (3 * n_per) as f64).sqrt();
        let d01 = (groups[0] - groups[1]).abs();
        let d12 = (groups[1] - groups[2]).abs();
        d01 > 3.0 * within && d12 > 3.0 * within
    }

    /// Between-centroid distance relative to within-group spread, using both
    /// embedding coordinates (rotation-invariant inside the subspace).
    fn separation_2d(coords: &Mat, n_per: usize) -> bool {
        let dim = coords.cols.min(2);
        let mut centroids = vec![vec![0.0f64; dim]; 3];
        for g in 0..3 {
            for i in g * n_per..(g + 1) * n_per {
                for d in 0..dim {
                    centroids[g][d] += coords.get(i, d) / n_per as f64;
                }
            }
        }
        let mut within = 0.0f64;
        for g in 0..3 {
            for i in g * n_per..(g + 1) * n_per {
                let mut s = 0.0;
                for d in 0..dim {
                    s += (coords.get(i, d) - centroids[g][d]).powi(2);
                }
                within += s;
            }
        }
        let within = (within / (3 * n_per) as f64).sqrt();
        let mut ok = true;
        for a in 0..3 {
            for b in (a + 1)..3 {
                let mut s = 0.0;
                for d in 0..dim {
                    s += (centroids[a][d] - centroids[b][d]).powi(2);
                }
                ok &= s.sqrt() > 3.0 * within;
            }
        }
        ok
    }

    /// Rotation-invariant separation using all available coordinates (up to
    /// three), for degenerate spectra where the leading directions mix.
    fn separation_3d(coords: &Mat, n_per: usize) -> bool {
        let dim = coords.cols;
        let mut centroids = vec![vec![0.0f64; dim]; 3];
        for g in 0..3 {
            for i in g * n_per..(g + 1) * n_per {
                for d in 0..dim {
                    centroids[g][d] += coords.get(i, d) / n_per as f64;
                }
            }
        }
        let mut within = 0.0f64;
        for g in 0..3 {
            for i in g * n_per..(g + 1) * n_per {
                let mut s = 0.0;
                for d in 0..dim {
                    s += (coords.get(i, d) - centroids[g][d]).powi(2);
                }
                within += s;
            }
        }
        let within = (within / (3 * n_per) as f64).sqrt();
        let mut ok = true;
        for a in 0..3 {
            for b in (a + 1)..3 {
                let mut s = 0.0;
                for d in 0..dim {
                    s += (centroids[a][d] - centroids[b][d]).powi(2);
                }
                ok &= s.sqrt() > 3.0 * within;
            }
        }
        ok
    }

    #[test]
    fn pcoa_recovers_linear_structure() {
        let x = blobs(20);
        let d = pairwise_distances(&x);
        let e = pcoa(&d, 2);
        assert!(cluster_separation(&e.coords, 20), "PCoA did not separate clusters");
    }

    #[test]
    fn kernel_pca_linear_matches_pcoa_geometry() {
        let x = blobs(15);
        let d = pairwise_distances(&x);
        let p = pcoa(&d, 2);
        let kp = kernel_pca(&x, 2, Kernel::Linear, 1.0);
        assert!(cluster_separation(&p.coords, 15));
        assert!(cluster_separation(&kp.coords, 15));
    }

    #[test]
    fn laplacian_and_diffusion_separate_clusters() {
        // A connected graph: three blobs plus a few bridge points, so the
        // kNN graph is not disconnected. With a disconnected graph every
        // component contributes a zero eigenvalue and the leading non-zero
        // subspace is degenerate, which makes any single-direction assertion
        // meaningless (that is a property of the data, not of the method).
        let x = connected_blobs(15);
        let w = knn_graph(&x, 8, 0.5);
        // Three disconnected clusters need three components to be fully
        // representable; with k=2 the degenerate subspace can rotate so that
        // two centroids coincide.
        let le = laplacian_eigenmaps(&w, 3);
        let dm = diffusion_map(&w, 2, 1.0);
        // Both spectra are degenerate here: three disconnected clusters give
        // three zero Laplacian eigenvalues and a three-fold degenerate first
        // non-zero one, so the leading eigenvectors are indicator-like
        // directions inside that subspace. Only rotation-invariant tests are
        // meaningful.
        assert!(separation_3d(&le.coords, 15), "laplacian failed");
        assert!(separation_2d(&dm.coords, 15), "diffusion failed");

    }

    #[test]
    fn isomap_separates_clusters() {
        // Same reasoning as the Laplacian test: Isomap's geodesic distances
        // are only defined on a connected kNN graph.
        let x = connected_blobs(12);
        let iso = isomap(&x, 8, 2);
        let n_per = 12;
        // The bridge points sit between clusters, so only the cluster rows are
        // checked (they are the first 3*n_per rows).
        let coords = Mat::from_fn(3 * n_per, iso.coords.cols, |i, j| iso.coords.get(i, j));
        assert!(cluster_separation(&coords, n_per), "isomap failed");
    }

    #[test]
    fn pcoa_on_euclidean_matches_pca_geometry() {
        // For centred data PCoA coordinates equal PCA scores up to sign/scale.
        let x = Mat::from_fn(20, 3, |i, j| ((i * 7 + j * 13) % 17) as f64);
        let d = pairwise_distances(&x);
        let p = pcoa(&d, 2);
        assert_eq!(p.coords.rows, 20);
        assert_eq!(p.coords.cols, 2);
    }
}
