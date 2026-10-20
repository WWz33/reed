//! Additional dimension-reduction families: t-SNE, NMF, LLE and PPCA.
//!
//! t-SNE follows van der Maaten & Hinton (2008): Gaussian affinities with a
//! per-point binary search on the perplexity, symmetrization, then the
//! Student-t heavy-tailed kernel optimized by gradient descent with momentum.
//!
//! NMF is Lee & Seung (1999) multiplicative updates.
//!
//! LLE is Roweis & Saul (2000): reconstruct each point from its neighbours,
//! then take the bottom eigenvectors of `(I - W)'(I - W)`.
//!
//! PPCA is Tipping & Bishop (1999) solved by EM, which handles missing values
//! natively (the E step simply skips unobserved entries).

use crate::backend::sym_eig;
use reed_core::linalg::{dist2, Mat};

fn rng_new(seed: u64) -> impl FnMut() -> f64 {
    let mut s = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    move || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// t-SNE embedding.
pub fn tsne(x: &Mat, cfg: &TsneConfig) -> Mat {
    let n = x.rows;
    let mut d2 = Mat::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            d2.set(i, j, dist2(x.row(i), x.row(j)));
        }
    }
    tsne_from_sq_distances(&d2, cfg)
}

/// t-SNE on precomputed **squared** distances.
///
/// This is the entry point for GRM-derived geometry: pass
/// `grm_distances(&grm)` squared and the embedding respects relatedness
/// instead of raw dosage differences.
pub fn tsne_from_sq_distances(d2: &Mat, cfg: &TsneConfig) -> Mat {
    let n = d2.rows;

    // Per-point sigma from the perplexity target.
    let target = cfg.perplexity.ln();
    let mut p = Mat::zeros(n, n);
    for i in 0..n {
        let mut lo = 1e-20f64;
        let mut hi = 1e20f64;
        let mut beta = 1.0f64;
        for _ in 0..64 {
            let mut sum = 0.0;
            let mut probs = vec![0.0f64; n];
            for j in 0..n {
                if i == j {
                    continue;
                }
                let v = (-beta * d2.get(i, j)).exp();
                probs[j] = v;
                sum += v;
            }
            if sum == 0.0 {
                break;
            }
            let mut h = 0.0;
            for j in 0..n {
                if i == j || probs[j] == 0.0 {
                    continue;
                }
                let pj = probs[j] / sum;
                h -= pj * pj.ln();
            }
            if (h - target).abs() < 1e-5 {
                for j in 0..n {
                    if i != j {
                        p.set(i, j, probs[j] / sum);
                    }
                }
                break;
            }
            if h > target {
                lo = beta;
                beta = if hi >= 1e20 { beta * 2.0 } else { (beta + hi) / 2.0 };
            } else {
                hi = beta;
                beta = (beta + lo) / 2.0;
            }
        }
    }

    // Symmetrize and normalize.
    let mut pj = Mat::zeros(n, n);
    let mut total = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            if i != j {
                let v = (p.get(i, j) + p.get(j, i)) / (2.0 * n as f64);
                pj.set(i, j, v);
                total += v;
            }
        }
    }
    if total > 0.0 {
        for v in pj.data.iter_mut() {
            *v = (*v / total).max(1e-12);
        }
    }

    // Initialize from a small Gaussian.
    let mut rnd = rng_new(cfg.seed);
    let dim = cfg.n_components;
    let mut y = Mat::from_fn(n, dim, |_, _| (rnd() - 0.5) * 1e-4);
    let mut vel = Mat::zeros(n, dim);
    let mut gains = Mat::from_fn(n, dim, |_, _| 1.0);

    for iter in 0..cfg.n_iter {
        let mut q = Mat::zeros(n, n);
        let mut sum_q = 0.0f64;
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let mut d2 = 0.0;
                for d in 0..dim {
                    let diff = y.get(i, d) - y.get(j, d);
                    d2 += diff * diff;
                }
                let v = 1.0 / (1.0 + d2);
                q.set(i, j, v);
                sum_q += v;
            }
        }
        if sum_q == 0.0 {
            break;
        }

        // Gradient: 4 * sum_j (p_ij - q_ij) * q_ij * (y_i - y_j)
        let mut grad = Mat::zeros(n, dim);
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let qij = q.get(i, j) / sum_q;
                let mult = 4.0 * (pj.get(i, j) - qij) * q.get(i, j);
                for d in 0..dim {
                    grad.set(i, d, grad.get(i, d) + mult * (y.get(i, d) - y.get(j, d)));
                }
            }
        }

        // Momentum with per-parameter gains (the original schedule).
        let momentum = if iter < 250 { 0.5 } else { 0.8 };
        for i in 0..n {
            for d in 0..dim {
                let g = grad.get(i, d);
                let same_sign = (g > 0.0) == (vel.get(i, d) > 0.0);
                let gain = if same_sign {
                    (gains.get(i, d) * 0.8).max(0.01)
                } else {
                    gains.get(i, d) + 0.2
                };
                gains.set(i, d, gain);
                let v = momentum * vel.get(i, d) - cfg.learning_rate * gain * g;
                vel.set(i, d, v);
                y.set(i, d, y.get(i, d) + v);
            }
        }
        // Zero-centre every iteration.
        for d in 0..dim {
            let mean: f64 = (0..n).map(|i| y.get(i, d)).sum::<f64>() / n as f64;
            for i in 0..n {
                y.set(i, d, y.get(i, d) - mean);
            }
        }
    }
    y
}

#[derive(Clone, Debug)]
pub struct TsneConfig {
    pub perplexity: f64,
    pub n_components: usize,
    pub n_iter: usize,
    pub learning_rate: f64,
    pub seed: u64,
}

impl Default for TsneConfig {
    fn default() -> Self {
        Self { perplexity: 30.0, n_components: 2, n_iter: 1000, learning_rate: 200.0, seed: 42 }
    }
}

/// Non-negative matrix factorization, `X ~ W H`, multiplicative updates.
pub fn nmf(x: &Mat, rank: usize, n_iter: usize, seed: u64) -> (Mat, Mat) {
    let n = x.rows;
    let m = x.cols;
    let mut rnd = rng_new(seed);
    let mut w = Mat::from_fn(n, rank, |_, _| rnd() * 0.1 + 1e-6);
    let mut h = Mat::from_fn(rank, m, |_, _| rnd() * 0.1 + 1e-6);

    for _ in 0..n_iter {
        // H <- H * (W'X) / (W'WH)
        let wt = w.transpose();
        let wtx = wt.matmul(x);
        let wtw = wt.matmul(&w);
        let wtw_h = wtw.matmul(&h);
        for i in 0..rank {
            for j in 0..m {
                let denom = wtw_h.get(i, j).max(1e-12);
                h.set(i, j, h.get(i, j) * wtx.get(i, j) / denom);
            }
        }
        // W <- W * (XH') / (WHH')
        let ht = h.transpose();
        let xht = x.matmul(&ht);
        let wh = w.matmul(&h);
        let whh = wh.matmul(&ht);
        for i in 0..n {
            for j in 0..rank {
                let denom = whh.get(i, j).max(1e-12);
                w.set(i, j, w.get(i, j) * xht.get(i, j) / denom);
            }
        }
    }
    (w, h)
}

/// Locally linear embedding.
pub fn lle(x: &Mat, k: usize, n_components: usize) -> Mat {
    let n = x.rows;
    let mut d2 = Mat::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            d2.set(i, j, dist2(x.row(i), x.row(j)));
        }
    }
    lle_from_sq_distances(x, &d2, k, n_components)
}

/// LLE with neighbourhoods taken from precomputed **squared** distances.
///
/// `features` is still needed for the local Gram matrix (the reconstruction
/// weights are defined in the feature space). Pass the standardized genotype
/// matrix with GRM distances to get relatedness-aware neighbourhoods.
pub fn lle_from_sq_distances(features: &Mat, d2: &Mat, k: usize, n_components: usize) -> Mat {
    let x = features;
    let n = x.rows;
    // Neighbours from the supplied distances.
    let mut nbrs = vec![vec![0usize; k]; n];
    for i in 0..n {
        let mut order: Vec<(f64, usize)> = (0..n)
            .filter(|&j| j != i)
            .map(|j| (d2.get(i, j), j))
            .collect();
        order.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        for (slot, &(_, j)) in order.iter().take(k).enumerate() {
            nbrs[i][slot] = j;
        }
    }

    // Reconstruction weights from the local Gram matrix, solved with a small
    // regularized linear system (G + reg*I) w = 1, then normalized.
    let mut wmat = Mat::zeros(n, n);
    for i in 0..n {
        let mut g = Mat::zeros(k, k);
        for a in 0..k {
            for b in 0..k {
                let mut s = 0.0;
                for d in 0..x.cols {
                    s += (x.get(i, d) - x.get(nbrs[i][a], d))
                        * (x.get(i, d) - x.get(nbrs[i][b], d));
                }
                g.set(a, b, s);
            }
        }
        let trace: f64 = (0..k).map(|a| g.get(a, a)).sum();
        let reg = 1e-3 * trace.max(1e-12) / k as f64;
        for a in 0..k {
            g.set(a, a, g.get(a, a) + reg);
        }
        let weights = solve_sym(&g, &vec![1.0; k]);
        let sum: f64 = weights.iter().sum();
        for a in 0..k {
            wmat.set(i, nbrs[i][a], if sum != 0.0 { weights[a] / sum } else { 0.0 });
        }
    }

    // M = (I - W)'(I - W); take the bottom eigenvectors (skip the trivial one).
    let mut m = Mat::zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            let wi = if i == j { 1.0 } else { 0.0 } - wmat.get(i, j);
            let mut acc = 0.0;
            for kk in 0..n {
                let a = if kk == i { 1.0 } else { 0.0 } - wmat.get(kk, i);
                let b = if kk == j { 1.0 } else { 0.0 } - wmat.get(kk, j);
                acc += a * b;
            }
            m.set(i, j, acc);
            let _ = wi;
        }
    }
    for i in 0..n {
        for j in 0..i {
            let s = 0.5 * (m.get(i, j) + m.get(j, i));
            m.set(i, j, s);
            m.set(j, i, s);
        }
    }
    let e = sym_eig(&m);
    let mut coords = Mat::zeros(n, n_components);
    for c in 0..n_components {
        let idx = 1 + c;
        if idx >= n {
            break;
        }
        for i in 0..n {
            coords.set(i, c, e.vectors.get(i, idx));
        }
    }
    coords
}

/// Solve a symmetric positive-definite system with Gaussian elimination.
fn solve_sym(a: &Mat, b: &[f64]) -> Vec<f64> {
    let n = a.rows;
    let mut m = a.clone();
    let mut x = b.to_vec();
    for col in 0..n {
        // Pivot.
        let mut piv = col;
        for r in (col + 1)..n {
            if m.get(r, col).abs() > m.get(piv, col).abs() {
                piv = r;
            }
        }
        if piv != col {
            for j in 0..n {
                let tmp = m.get(col, j);
                m.set(col, j, m.get(piv, j));
                m.set(piv, j, tmp);
            }
            x.swap(col, piv);
        }
        let d = m.get(col, col);
        if d.abs() < 1e-14 {
            continue;
        }
        for r in (col + 1)..n {
            let factor = m.get(r, col) / d;
            for j in col..n {
                let v = m.get(r, j) - factor * m.get(col, j);
                m.set(r, j, v);
            }
            x[r] -= factor * x[col];
        }
    }
    let mut out = vec![0.0f64; n];
    for col in (0..n).rev() {
        let d = m.get(col, col);
        if d.abs() < 1e-14 {
            continue;
        }
        let mut s = x[col];
        for j in (col + 1)..n {
            s -= m.get(col, j) * out[j];
        }
        out[col] = s / d;
    }
    out
}

/// Probabilistic PCA by EM; missing entries (NaN) are skipped in the E step.
pub fn ppca(x: &Mat, n_components: usize, n_iter: usize, seed: u64) -> Mat {
    let n = x.rows;
    let m = x.cols;

    // Column means over observed entries.
    let mut mean = vec![0.0f64; m];
    for j in 0..m {
        let mut s = 0.0;
        let mut c = 0usize;
        for i in 0..n {
            let v = x.get(i, j);
            if !v.is_nan() {
                s += v;
                c += 1;
            }
        }
        mean[j] = if c > 0 { s / c as f64 } else { 0.0 };
    }

    let mut rnd = rng_new(seed);
    let q = n_components;
    let mut w = Mat::from_fn(m, q, |_, _| (rnd() - 0.5) * 0.1);
    let mut sigma2 = 1.0f64;

    for _ in 0..n_iter {
        // E step: per-sample posterior mean and covariance.
        let mut mz = vec![vec![0.0f64; q]; n];
        for i in 0..n {
            // M = sigma2*I + W_obs' W_obs  (q x q), using only observed cols.
            let mut mm = Mat::zeros(q, q);
            for a in 0..q {
                mm.set(a, a, sigma2);
            }
            let mut rhs = vec![0.0f64; q];
            for j in 0..m {
                let v = x.get(i, j);
                if v.is_nan() {
                    continue;
                }
                let centered = v - mean[j];
                for a in 0..q {
                    rhs[a] += w.get(j, a) * centered;
                    for b in 0..q {
                        mm.set(a, b, mm.get(a, b) + w.get(j, a) * w.get(j, b));
                    }
                }
            }
            let sol = solve_sym(&mm, &rhs);
            for a in 0..q {
                mz[i][a] = sol[a];
            }
        }

        // M step: W and sigma2.
        let mut wnum = Mat::zeros(m, q);
        let mut ezzt = Mat::zeros(q, q);
        for i in 0..n {
            for a in 0..q {
                for b in 0..q {
                    ezzt.set(a, b, ezzt.get(a, b) + mz[i][a] * mz[i][b]);
                }
            }
        }
        for j in 0..m {
            for a in 0..q {
                let mut s = 0.0;
                for i in 0..n {
                    let v = x.get(i, j);
                    if v.is_nan() {
                        continue;
                    }
                    s += (v - mean[j]) * mz[i][a];
                }
                wnum.set(j, a, s);
            }
        }
        // W = Wnum * inv(EZZ' + sigma2*I). Build the inverse by solving
        // against each unit vector, then multiply as a matrix.
        let mut reg = ezzt.clone();
        for a in 0..q {
            reg.set(a, a, reg.get(a, a) + sigma2 + 1e-9);
        }
        let mut inv = Mat::zeros(q, q);
        for a in 0..q {
            let mut unit = vec![0.0f64; q];
            unit[a] = 1.0;
            let col = solve_sym(&reg, &unit);
            for r in 0..q {
                inv.set(r, a, col[r]);
            }
        }
        let new_w = wnum.matmul(&inv);
        w = new_w;
        // sigma2 = mean squared residual over observed entries.
        let mut resid = 0.0f64;
        let mut cnt = 0usize;
        for i in 0..n {
            for j in 0..m {
                let v = x.get(i, j);
                if v.is_nan() {
                    continue;
                }
                let mut recon = mean[j];
                for a in 0..q {
                    recon += w.get(j, a) * mz[i][a];
                }
                resid += (v - recon).powi(2);
                cnt += 1;
            }
        }
        sigma2 = if cnt > 0 { resid / cnt as f64 } else { sigma2 };
        let _ = &mut wnum;
    }

    // Latent coordinates from the final E step.
    let mut out = Mat::zeros(n, q);
    for i in 0..n {
        let mut mm = Mat::zeros(q, q);
        for a in 0..q {
            mm.set(a, a, sigma2);
        }
        let mut rhs = vec![0.0f64; q];
        for j in 0..m {
            let v = x.get(i, j);
            if v.is_nan() {
                continue;
            }
            let centered = v - mean[j];
            for a in 0..q {
                rhs[a] += w.get(j, a) * centered;
                for b in 0..q {
                    mm.set(a, b, mm.get(a, b) + w.get(j, a) * w.get(j, b));
                }
            }
        }
        let sol = solve_sym(&mm, &rhs);
        for a in 0..q {
            out.set(i, a, sol[a]);
        }
    }
    out
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

    fn separation(coords: &Mat, n_per: usize) -> f64 {
        let dim = coords.cols;
        let mut cent = vec![vec![0.0f64; dim]; 3];
        for g in 0..3 {
            for i in g * n_per..(g + 1) * n_per {
                for d in 0..dim {
                    cent[g][d] += coords.get(i, d) / n_per as f64;
                }
            }
        }
        let mut within = 0.0f64;
        for g in 0..3 {
            for i in g * n_per..(g + 1) * n_per {
                let mut s = 0.0;
                for d in 0..dim {
                    s += (coords.get(i, d) - cent[g][d]).powi(2);
                }
                within += s;
            }
        }
        let within = (within / (3 * n_per) as f64).sqrt().max(1e-12);
        let mut min_sep = f64::MAX;
        for a in 0..3 {
            for b in (a + 1)..3 {
                let mut s = 0.0;
                for d in 0..dim {
                    s += (cent[a][d] - cent[b][d]).powi(2);
                }
                min_sep = min_sep.min(s.sqrt());
            }
        }
        min_sep / within
    }

    #[test]
    fn tsne_separates_blobs() {
        let x = blobs(20);
        let cfg = TsneConfig { perplexity: 10.0, n_iter: 300, seed: 3, ..Default::default() };
        let y = tsne(&x, &cfg);
        let sep = separation(&y, 20);
        assert!(sep > 2.0, "t-SNE separation ratio {sep}");
    }

    #[test]
    fn nmf_reconstructs_nonnegative_data() {
        let x = Mat::from_fn(30, 8, |i, j| (((i * 7 + j * 11) % 13) as f64) * 0.5);
        let (w, h) = nmf(&x, 3, 200, 1);
        let recon = w.matmul(&h);
        let mut err = 0.0;
        let mut base = 0.0;
        for i in 0..30 {
            for j in 0..8 {
                err += (x.get(i, j) - recon.get(i, j)).powi(2);
                base += x.get(i, j).powi(2);
            }
        }
        assert!(err < 0.25 * base, "relative reconstruction error {}", err / base);
        for v in w.data.iter().chain(h.data.iter()) {
            assert!(*v >= 0.0);
        }
    }

    #[test]
    fn lle_separates_blobs() {
        let x = blobs(15);
        let y = lle(&x, 8, 2);
        let sep = separation(&y, 15);
        assert!(sep > 1.5, "LLE separation ratio {sep}");
    }

    #[test]
    fn ppca_handles_missing_values() {
        let mut x = blobs(10);
        // Punch holes in 20% of the entries.
        for i in 0..x.rows {
            for j in 0..x.cols {
                if (i * 7 + j * 13) % 5 == 0 {
                    x.set(i, j, f64::NAN);
                }
            }
        }
        let y = ppca(&x, 2, 60, 5);
        let sep = separation(&y, 10);
        assert!(sep > 1.5, "PPCA separation ratio {sep}");
    }
}
