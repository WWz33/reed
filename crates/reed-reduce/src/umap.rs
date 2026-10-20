//! UMAP prototype, following `umap-learn` 0.5.x line by line:
//!
//! 1. kNN graph on the input (`n_neighbors`, default 15).
//! 2. `smooth_knn_dist`: per-point binary search for sigma such that the
//!    fuzzy cardinality equals `log2(k)`; rho = distance to the nearest
//!    non-zero neighbour (`local_connectivity = 1`).
//! 3. Membership strengths `exp(-(d - rho)/sigma)`.
//! 4. Fuzzy union of the directed graph with its transpose:
//!    `A + A' - A * A'`.
//! 5. Spectral initialization (largest eigenvectors of the normalized
//!    Laplacian), then SGD on the cross-entropy with negative sampling.
//!
//! Source references (umap-learn 0.5.12):
//! `umap/umap_.py:153` smooth_knn_dist, `:374` compute_membership_strengths,
//! `:455` fuzzy_simplicial_set, `umap/layouts.py:63` the SGD epoch,
//! `umap/umap_.py:1478` find_ab_params.

use crate::backend::sym_eig;
use reed_core::linalg::{dist2, Mat};

const SMOOTH_K_TOLERANCE: f64 = 1e-5;
const MIN_K_DIST_SCALE: f64 = 1e-3;

/// UMAP hyper-parameters with the umap-learn defaults.
#[derive(Clone, Debug)]
pub struct UmapConfig {
    pub n_neighbors: usize,
    pub n_components: usize,
    pub n_epochs: usize,
    pub min_dist: f64,
    pub spread: f64,
    pub negative_sample_rate: usize,
    pub learning_rate: f64,
    pub local_connectivity: f64,
    pub set_op_mix_ratio: f64,
    pub gamma: f64,
    pub seed: u64,
}

impl Default for UmapConfig {
    fn default() -> Self {
        Self {
            n_neighbors: 15,
            n_components: 2,
            n_epochs: 200,
            min_dist: 0.1,
            spread: 1.0,
            negative_sample_rate: 5,
            learning_rate: 1.0,
            local_connectivity: 1.0,
            set_op_mix_ratio: 1.0,
            gamma: 1.0,
            seed: 42,
        }
    }
}

struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Self(seed.wrapping_mul(2685821657736338717).wrapping_add(1))
    }
    fn next_u64(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }
}

/// `smooth_knn_dist` (umap_.py:153): returns (sigmas, rhos).
fn smooth_knn_dist(knn_dists: &Mat, k: f64, local_connectivity: f64) -> (Vec<f64>, Vec<f64>) {
    let n = knn_dists.rows;
    let n_neighbors = knn_dists.cols;
    let target = k.log2();
    let mut rho = vec![0.0f64; n];
    let mut result = vec![0.0f64; n];

    // Mean of all finite distances.
    let mut sum = 0.0;
    let mut cnt = 0usize;
    for i in 0..n {
        for j in 0..n_neighbors {
            let v = knn_dists.get(i, j);
            if v.is_finite() {
                sum += v;
                cnt += 1;
            }
        }
    }
    let mean_distances = if cnt > 0 { sum / cnt as f64 } else { 0.0 };

    for i in 0..n {
        let mut lo = 0.0f64;
        let mut hi = f64::MAX;
        let mut mid = 1.0f64;

        // rho from the nearest non-zero distances.
        let mut non_zero: Vec<f64> = (0..n_neighbors)
            .map(|j| knn_dists.get(i, j))
            .filter(|v| v.is_finite() && *v > 0.0)
            .collect();
        non_zero.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if non_zero.len() as f64 >= local_connectivity {
            let index = local_connectivity.floor() as usize;
            let interpolation = local_connectivity - index as f64;
            if index > 0 {
                rho[i] = non_zero[index - 1];
                if interpolation > SMOOTH_K_TOLERANCE {
                    rho[i] += interpolation * (non_zero[index] - non_zero[index - 1]);
                }
            } else {
                rho[i] = interpolation * non_zero[0];
            }
        } else if !non_zero.is_empty() {
            rho[i] = *non_zero.last().unwrap();
        }

        for _ in 0..64 {
            let mut psum = 0.0;
            for j in 1..n_neighbors {
                let d = knn_dists.get(i, j) - rho[i];
                if d > 0.0 {
                    psum += (-(d / mid)).exp();
                } else {
                    psum += 1.0;
                }
            }
            if (psum - target).abs() < SMOOTH_K_TOLERANCE {
                break;
            }
            if psum > target {
                hi = mid;
                mid = (lo + hi) / 2.0;
            } else {
                lo = mid;
                if hi >= f64::MAX {
                    mid *= 2.0;
                } else {
                    mid = (lo + hi) / 2.0;
                }
            }
        }
        result[i] = mid;

        // Floor on sigma (MIN_K_DIST_SCALE).
        if rho[i] > 0.0 {
            let mean_i: f64 = (0..n_neighbors)
                .map(|j| knn_dists.get(i, j))
                .filter(|v| v.is_finite())
                .sum::<f64>()
                / n_neighbors as f64;
            if result[i] < MIN_K_DIST_SCALE * mean_i {
                result[i] = MIN_K_DIST_SCALE * mean_i;
            }
        } else if result[i] < MIN_K_DIST_SCALE * mean_distances {
            result[i] = MIN_K_DIST_SCALE * mean_distances;
        }
    }
    (result, rho)
}

/// The fuzzy simplicial set: a dense symmetric membership matrix.
pub struct FuzzyGraph {
    pub weights: Mat,
    pub knn_indices: Vec<Vec<usize>>,
    pub knn_dists: Mat,
}

/// Build the fuzzy graph from data.
pub fn fuzzy_simplicial_set(x: &Mat, cfg: &UmapConfig) -> FuzzyGraph {
    let n = x.rows;
    let mut d = Mat::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            d.set(i, j, dist2(x.row(i), x.row(j)).sqrt());
        }
    }
    fuzzy_simplicial_set_from_distances(&d, cfg)
}

/// Build the fuzzy graph from a precomputed distance matrix.
///
/// This is the GRM entry point: `grm_distances(&grm)` feeds relatedness-aware
/// neighbourhoods into the same fuzzy-set + SGD machinery.
pub fn fuzzy_simplicial_set_from_distances(d: &Mat, cfg: &UmapConfig) -> FuzzyGraph {
    let n = d.rows;
    let k = cfg.n_neighbors.min(n - 1);

    // Exact kNN (the prototype does not implement NN-descent).
    let mut knn_indices = vec![vec![0usize; k + 1]; n];
    let mut knn_dists = Mat::zeros(n, k + 1);
    for i in 0..n {
        let mut order: Vec<(f64, usize)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (d.get(i, j), j))
            .collect();
        order.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        knn_indices[i][0] = i;
        knn_dists.set(i, 0, 0.0);
        for (slot, &(d, j)) in order.iter().take(k).enumerate() {
            knn_indices[i][slot + 1] = j;
            knn_dists.set(i, slot + 1, d);
        }
    }

    let (sigmas, rhos) = smooth_knn_dist(&knn_dists, k as f64, cfg.local_connectivity);

    // Membership strengths (umap_.py:374).
    let mut directed = Mat::zeros(n, n);
    for i in 0..n {
        for slot in 0..=k {
            let j = knn_indices[i][slot];
            let d = knn_dists.get(i, slot);
            let val = if j == i {
                0.0
            } else if d - rhos[i] <= 0.0 || sigmas[i] == 0.0 {
                1.0
            } else {
                (-((d - rhos[i]) / sigmas[i])).exp()
            };
            if val > directed.get(i, j) {
                directed.set(i, j, val);
            }
        }
    }

    // Fuzzy union: A + A' - A * A' (set_op_mix_ratio = 1).
    let mut weights = Mat::zeros(n, n);
    let r = cfg.set_op_mix_ratio;
    for i in 0..n {
        for j in 0..n {
            let a = directed.get(i, j);
            let b = directed.get(j, i);
            let prod = a * b;
            let union = a + b - prod;
            let v = if r >= 1.0 {
                union
            } else if r <= 0.0 {
                prod
            } else {
                r * union + (1.0 - r) * prod
            };
            weights.set(i, j, v);
        }
    }

    FuzzyGraph { weights, knn_indices, knn_dists }
}

/// `find_ab_params` (umap_.py:1478): fit `1/(1 + a x^{2b})` to the target
/// curve by a small Levenberg-Marquardt style least-squares solve.
pub fn find_ab_params(spread: f64, min_dist: f64) -> (f64, f64) {
    // The original uses scipy curve_fit. Reproduce the fit with a coarse
    // grid + Gauss-Newton refinement over the same objective.
    let xs: Vec<f64> = (0..300).map(|i| spread * 3.0 * i as f64 / 299.0).collect();
    let ys: Vec<f64> = xs
        .iter()
        .map(|&x| if x < min_dist { 1.0 } else { (-(x - min_dist) / spread).exp() })
        .collect();

    let model = |a: f64, b: f64, x: f64| 1.0 / (1.0 + a * x.powf(2.0 * b));

    // Grid search for a starting point (scipy's curve_fit starts at 1,1).
    let mut best = (1.0f64, 1.0f64);
    let mut best_err = f64::MAX;
    for ai in 0..40 {
        let a = 0.05 * 1.18f64.powi(ai);
        for bi in 0..40 {
            let b = 0.2 + 0.05 * bi as f64;
            let err: f64 = xs
                .iter()
                .zip(&ys)
                .map(|(&x, &y)| (model(a, b, x) - y).powi(2))
                .sum();
            if err < best_err {
                best_err = err;
                best = (a, b);
            }
        }
    }

    // Gauss-Newton refinement with numerical derivatives.
    let (mut a, mut b) = best;
    for _ in 0..50 {
        let mut jtj = [[0.0f64; 2]; 2];
        let mut jtr = [0.0f64; 2];
        for (&x, &y) in xs.iter().zip(&ys) {
            let f = model(a, b, x);
            let r = f - y;
            let da = {
                let h = 1e-6 * a.abs().max(1.0);
                (model(a + h, b, x) - f) / h
            };
            let db = {
                let h = 1e-6;
                (model(a, b + h, x) - f) / h
            };
            jtj[0][0] += da * da;
            jtj[0][1] += da * db;
            jtj[1][0] += da * db;
            jtj[1][1] += db * db;
            jtr[0] += da * r;
            jtr[1] += db * r;
        }
        let det = jtj[0][0] * jtj[1][1] - jtj[0][1] * jtj[1][0];
        if det.abs() < 1e-18 {
            break;
        }
        let d_a = (jtj[1][1] * jtr[0] - jtj[0][1] * jtr[1]) / det;
        let d_b = (-jtj[1][0] * jtr[0] + jtj[0][0] * jtr[1]) / det;
        a -= d_a;
        b -= d_b;
        if d_a.abs() < 1e-12 && d_b.abs() < 1e-12 {
            break;
        }
    }
    (a.max(1e-9), b.max(1e-9))
}

/// Spectral initialization: the top eigenvectors of the normalized Laplacian
/// of the fuzzy graph (`umap/spectral.py`).
fn spectral_init(g: &Mat, dim: usize, rng: &mut Rng) -> Mat {
    let n = g.rows;
    let mut deg = vec![0.0f64; n];
    for i in 0..n {
        deg[i] = (0..n).map(|j| g.get(i, j)).sum();
    }
    let mut l = Mat::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let diag = if i == j { 1.0 } else { 0.0 };
            let denom = (deg[i] * deg[j]).sqrt();
            let off = if denom > 0.0 { g.get(i, j) / denom } else { 0.0 };
            l.set(i, j, diag - off);
        }
    }
    let e = sym_eig(&l);
    let mut coords = Mat::zeros(n, dim);
    // Skip the trivial smallest eigenvalue; take the next `dim`.
    for c in 0..dim {
        let idx = 1 + c;
        if idx >= n {
            break;
        }
        for i in 0..n {
            coords.set(i, c, e.vectors.get(i, idx));
        }
    }
    // Scale to a spread comparable to the data, as umap does.
    let mut max_abs: f64 = 0.0;
    for v in coords.data.iter() {
        max_abs = max_abs.max(v.abs());
    }
    if max_abs > 0.0 {
        let scale = 10.0 / max_abs;
        for v in coords.data.iter_mut() {
            *v *= scale;
        }
    }
    // Tiny jitter breaks the symmetry of degenerate spectra.
    for v in coords.data.iter_mut() {
        *v += 1e-4 * (rng.next_f64() - 0.5);
    }
    coords
}

/// UMAP embedding.
pub fn umap(x: &Mat, cfg: &UmapConfig) -> Mat {
    let n = x.rows;
    let mut d = Mat::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            d.set(i, j, dist2(x.row(i), x.row(j)).sqrt());
        }
    }
    umap_from_distances(&d, cfg)
}

/// UMAP on a precomputed distance matrix (the GRM entry point).
pub fn umap_from_distances(d: &Mat, cfg: &UmapConfig) -> Mat {
    let n = d.rows;
    let g = fuzzy_simplicial_set_from_distances(d, cfg);
    let (a, b) = find_ab_params(cfg.spread, cfg.min_dist);
    let mut rng = Rng::new(cfg.seed);
    let mut emb = spectral_init(&g.weights, cfg.n_components, &mut rng);

    // Edge list from the fuzzy graph (both directions).
    let mut head = Vec::new();
    let mut tail = Vec::new();
    for i in 0..n {
        for j in 0..n {
            if g.weights.get(i, j) > 0.0 {
                head.push(i);
                tail.push(j);
            }
        }
    }
    let n_edges = head.len();

    // Epoch scheduling (umap's epochs_per_sample / epochs_per_negative_sample).
    let mut epochs_per_sample = vec![0.0f64; n_edges];
    let max_weight = g
        .weights
        .data
        .iter()
        .fold(0.0f64, |acc, &v| acc.max(v));
    for e in 0..n_edges {
        let w = g.weights.get(head[e], tail[e]);
        epochs_per_sample[e] = if max_weight > 0.0 { max_weight / w } else { 1.0 };
    }
    let epochs_per_negative_sample: Vec<f64> = epochs_per_sample
        .iter()
        .map(|&e| e / cfg.negative_sample_rate as f64)
        .collect();

    let dim = cfg.n_components;
    let mut next_sample = epochs_per_sample.clone();
    let mut next_neg = epochs_per_negative_sample.clone();

    for epoch in 0..cfg.n_epochs {
        // umap decays alpha linearly to 1/4 of its initial value.
        let progress = epoch as f64 / cfg.n_epochs as f64;
        let alpha = cfg.learning_rate * (1.0 - progress).max(0.25);

        for e in 0..n_edges {
            if next_sample[e] > epoch as f64 {
                continue;
            }
            let j = head[e];
            let k = tail[e];

            // Attractive gradient.
            let mut d2 = 0.0;
            for d in 0..dim {
                let diff = emb.get(j, d) - emb.get(k, d);
                d2 += diff * diff;
            }
            let grad_coeff = if d2 > 0.0 {
                -2.0 * a * b * d2.powf(b - 1.0) / (a * d2.powf(b) + 1.0)
            } else {
                0.0
            };
            for d in 0..dim {
                let gd = (grad_coeff * (emb.get(j, d) - emb.get(k, d))).clamp(-4.0, 4.0);
                emb.set(j, d, emb.get(j, d) + gd * alpha);
                emb.set(k, d, emb.get(k, d) - gd * alpha);
            }
            next_sample[e] += epochs_per_sample[e];

            // Negative sampling.
            let mut n_neg = ((epoch as f64 - next_neg[e]) / epochs_per_negative_sample[e]) as usize;
            while n_neg > 0 {
                let kk = rng.below(n);
                let mut d2 = 0.0;
                for d in 0..dim {
                    let diff = emb.get(j, d) - emb.get(kk, d);
                    d2 += diff * diff;
                }
                let grad_coeff = if d2 > 0.0 && kk != j {
                    2.0 * cfg.gamma * b / ((0.001 + d2) * (a * d2.powf(b) + 1.0))
                } else {
                    0.0
                };
                for d in 0..dim {
                    let gd = (grad_coeff * (emb.get(j, d) - emb.get(kk, d))).clamp(-4.0, 4.0);
                    emb.set(j, d, emb.get(j, d) + gd * alpha);
                }
                n_neg -= 1;
                next_neg[e] += epochs_per_negative_sample[e];
            }
        }
    }
    emb
}

#[cfg(test)]
mod tests {
    use super::*;

    fn blobs(n_per: usize) -> Mat {
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

    #[test]
    fn ab_params_are_in_the_expected_range() {
        let (a, b) = find_ab_params(1.0, 0.1);
        // umap-learn gives a ~= 1.5769, b ~= 0.8951 for (spread, min_dist) =
        // (1.0, 0.1).
        assert!((a - 1.577).abs() < 0.05, "a = {a}");
        assert!((b - 0.895).abs() < 0.02, "b = {b}");
    }

    #[test]
    fn fuzzy_graph_is_symmetric_and_bounded() {
        let x = blobs(10);
        let cfg = UmapConfig { n_neighbors: 5, ..Default::default() };
        let g = fuzzy_simplicial_set(&x, &cfg);
        for i in 0..g.weights.rows {
            for j in 0..g.weights.cols {
                let v = g.weights.get(i, j);
                assert!((0.0..=1.0).contains(&v), "weight {v} out of range");
                assert!((v - g.weights.get(j, i)).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn umap_separates_three_blobs() {
        let x = blobs(20);
        let cfg = UmapConfig {
            n_neighbors: 10,
            n_epochs: 200,
            seed: 7,
            ..Default::default()
        };
        let emb = umap(&x, &cfg);
        // Cluster separation: compare centroid distances to the within spread.
        let n_per = 20;
        let dim = 2;
        let mut centroids = vec![vec![0.0f64; dim]; 3];
        for g in 0..3 {
            for i in g * n_per..(g + 1) * n_per {
                for d in 0..dim {
                    centroids[g][d] += emb.get(i, d) / n_per as f64;
                }
            }
        }
        let mut within = 0.0f64;
        for g in 0..3 {
            for i in g * n_per..(g + 1) * n_per {
                let mut s = 0.0;
                for d in 0..dim {
                    s += (emb.get(i, d) - centroids[g][d]).powi(2);
                }
                within += s;
            }
        }
        let within = (within / (3 * n_per) as f64).sqrt();
        let mut min_sep = f64::MAX;
        for a in 0..3 {
            for b in (a + 1)..3 {
                let mut s = 0.0;
                for d in 0..dim {
                    s += (centroids[a][d] - centroids[b][d]).powi(2);
                }
                min_sep = min_sep.min(s.sqrt());
            }
        }
        assert!(min_sep > 2.0 * within, "min separation {min_sep} vs spread {within}");
    }

    #[test]
    fn smooth_knn_dist_cardinality_is_log2_k() {
        let x = blobs(10);
        let cfg = UmapConfig { n_neighbors: 8, ..Default::default() };
        let g = fuzzy_simplicial_set(&x, &cfg);
        // For each point, sum of exp(-(d-rho)/sigma) over neighbours should be
        // close to log2(k).
        let n = x.rows;
        let k = cfg.n_neighbors;
        let (sigmas, rhos) = smooth_knn_dist(&g.knn_dists, k as f64, 1.0);
        for i in 0..n {
            let mut psum = 0.0;
            for j in 1..=k {
                let d = g.knn_dists.get(i, j) - rhos[i];
                psum += if d > 0.0 { (-(d / sigmas[i])).exp() } else { 1.0 };
            }
            let target = (k as f64).log2();
            assert!(
                (psum - target).abs() < 0.05,
                "point {i}: cardinality {psum} vs target {target}"
            );
        }
    }
}
