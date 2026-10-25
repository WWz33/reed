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
use reed_core::rng::Rng;

/// Errors from t-SNE.
#[derive(Clone, Debug, PartialEq)]
pub enum TsneError {
    /// The perplexity was below 1.0: a probability distribution's Shannon
    /// entropy is non-negative, so `perplexity = exp(H) >= 1` and no beta can
    /// reach a smaller target.
    PerplexityBelowOne {
        /// The requested perplexity.
        perplexity: f64,
    },
    /// The perplexity was at least `n - 1`: the per-point binary search
    /// targets the entropy `ln(perplexity)`, but a distribution over the
    /// `n - 1` other points tops out at the uniform `ln(n - 1)`.
    PerplexityTooLarge {
        /// The requested perplexity.
        perplexity: f64,
        /// The number of samples.
        n: usize,
    },
}

impl std::fmt::Display for TsneError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::PerplexityBelowOne { perplexity } => {
                write!(f, "perplexity {perplexity} is below 1, the smallest any probability distribution can realize")
            }
            Self::PerplexityTooLarge { perplexity, n } => write!(
                f,
                "perplexity {perplexity} is at least n - 1 = {}, the largest a distribution over the other samples can realize",
                n.saturating_sub(1)
            ),
        }
    }
}

impl std::error::Error for TsneError {}

/// t-SNE embedding.
///
/// # Errors
///
/// [`TsneError`] when the perplexity cannot be realized for this sample
/// count; before doc_v20 S29 the per-point binary search then failed to
/// converge and silently zeroed every affinity row, so the default
/// perplexity of 30 on `n <= 31` samples returned a pure-noise embedding
/// instead of an error.
pub fn tsne(x: &Mat, cfg: &TsneConfig) -> Result<Mat, TsneError> {
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
///
/// # Errors
///
/// [`TsneError`]; see [`tsne`].
pub fn tsne_from_sq_distances(d2: &Mat, cfg: &TsneConfig) -> Result<Mat, TsneError> {
    let n = d2.rows;

    // doc_v20 S29: reject perplexities no beta can realize before running
    // the search. The binary search targets the entropy ln(perplexity), but
    // a distribution over the n - 1 other points spans only [0, ln(n - 1)],
    // so perplexities below 1 or at/above n - 1 never converge — and an
    // unconverged row used to be left all-zero, which with the default 30
    // on n <= 31 samples silently turned the whole embedding into noise.
    if cfg.perplexity < 1.0 {
        return Err(TsneError::PerplexityBelowOne { perplexity: cfg.perplexity });
    }
    if cfg.perplexity >= n.saturating_sub(1) as f64 {
        return Err(TsneError::PerplexityTooLarge { perplexity: cfg.perplexity, n });
    }

    // Per-point sigma from the perplexity target.
    let target = cfg.perplexity.ln();
    let mut p = Mat::zeros(n, n);
    for i in 0..n {
        let mut lo = 1e-20f64;
        let mut hi = 1e20f64;
        // doc_v20 review of S29: start the search at the scale the data
        // actually lives on. Dosage distances are O(number of sites), so
        // `exp(-1 * d^2)` underflows for any realistic genotype matrix and
        // a fixed `beta = 1` start collapsed every row to zero (pure
        // repulsion, no error). `1 / mean(d^2)` is the scale at which the
        // weights are O(1) before the search starts.
        let mean_d2 = (0..n)
            .filter(|&j| j != i)
            .map(|j| d2.get(i, j))
            .sum::<f64>()
            / (n.saturating_sub(1)).max(1) as f64;
        let mut beta = if mean_d2 > 0.0 && mean_d2.is_finite() { 1.0 / mean_d2 } else { 1.0 };
        let mut probs = vec![0.0f64; n];
        let mut sum = 0.0f64;
        for _ in 0..64 {
            sum = 0.0;
            for j in 0..n {
                if i == j {
                    continue;
                }
                let v = (-beta * d2.get(i, j)).exp();
                probs[j] = v;
                sum += v;
            }
            if sum == 0.0 {
                // Every weight underflowed: beta is too cold. Move the upper
                // bound down and keep searching (used to `break`, which left
                // the row unwritten — doc_v20 review of S29).
                hi = beta;
                beta = (beta + lo) / 2.0;
                continue;
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
        // doc_v20 S29: always write the Shannon-normalized row computed from
        // the final beta, converged or not. Exhausting the 64 iterations
        // used to leave the row all-zero so the sample exerted no force at
        // all; with the final weights an unconverged row degrades gracefully
        // instead. (sum == 0 means every exp underflowed and there is no
        // distribution to normalize.)
        if sum > 0.0 {
            for j in 0..n {
                if i != j {
                    p.set(i, j, probs[j] / sum);
                }
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
    let mut rng = Rng::from_seed(cfg.seed);
    let dim = cfg.n_components;
    let mut y = Mat::from_fn(n, dim, |_, _| (rng.next_f64() - 0.5) * 1e-4);
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

        let grad = tsne_grad(&pj, &q, sum_q, &y);

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
    Ok(y)
}

/// Gradient of the t-SNE KL cost w.r.t. the embedding points.
///
/// With `w_ij = 1 / (1 + ||y_i - y_j||^2)` and `q_ij = w_ij / sum_q`, the
/// cost `C = sum_{i!=j} p_ij * ln(p_ij / q_ij)` differentiates to
/// `dC/dy_i = 4 * sum_j (p_ij - q_ij) * w_ij * (y_i - y_j)`: the residual
/// uses the NORMALIZED `q_ij`, the last factor is the RAW kernel `w_ij`
/// (the `1/Z` in `q` cancels against the `Z` that differentiating the
/// softmax denominator produces). Verified against a central finite
/// difference of the cost in `r8_tsne_gradient_matches_the_kl_cost`.
fn tsne_grad(pj: &Mat, q: &Mat, sum_q: f64, y: &Mat) -> Mat {
    let n = y.rows;
    let dim = y.cols;
    let mut grad = Mat::zeros(n, dim);
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let qij = q.get(i, j) / sum_q;
            // Normalized residual, raw kernel: doc_v20 review of R8.
            let mult = 4.0 * (pj.get(i, j) - qij) * q.get(i, j);
            for d in 0..dim {
                grad.set(i, d, grad.get(i, d) + mult * (y.get(i, d) - y.get(j, d)));
            }
        }
    }
    grad
}

/// The t-SNE KL cost `C = sum_{i!=j} p_ij * ln(p_ij / q_ij)` with
/// `q_ij = w_ij / sum_q`. Test-only: the gradient above is checked against
/// a central finite difference of this function.
#[cfg(test)]
fn tsne_kl_cost(pj: &Mat, y: &Mat) -> f64 {
    let n = y.rows;
    let dim = y.cols;
    let mut w = Mat::zeros(n, n);
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
            w.set(i, j, v);
            sum_q += v;
        }
    }
    let mut cost = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let qij = w.get(i, j) / sum_q;
            let pij = pj.get(i, j);
            cost += pij * (pij / qij).ln();
        }
    }
    cost
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
    let mut rng = Rng::from_seed(seed);
    let mut w = Mat::from_fn(n, rank, |_, _| rng.next_f64() * 0.1 + 1e-6);
    let mut h = Mat::from_fn(rank, m, |_, _| rng.next_f64() * 0.1 + 1e-6);

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
        // total_cmp: NaN squared distances sort last instead of panicking
        // (doc_v20 S31).
        order.sort_by(|a, b| a.0.total_cmp(&b.0));
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

    let mut rng = Rng::from_seed(seed);
    let q = n_components;
    let mut w = Mat::from_fn(m, q, |_, _| (rng.next_f64() - 0.5) * 0.1);
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
    fn s29_impossible_perplexity_is_rejected() {
        // doc_v20 S29: a distribution over the n - 1 other points has at
        // most ln(n - 1) entropy, so perplexity >= n - 1 can never be
        // reached; with the default 30 and n <= 31 the binary search then
        // failed for every row, which silently zeroed all affinities and
        // returned a pure-noise embedding (pre-fix: Ok, not Err).
        let x = Mat::from_fn(4, 2, |i, j| (i * 3 + j) as f64 * 0.7);
        for bad in [3.0, 30.0] {
            let cfg = TsneConfig { perplexity: bad, n_iter: 20, ..Default::default() };
            assert_eq!(
                tsne(&x, &cfg),
                Err(TsneError::PerplexityTooLarge { perplexity: bad, n: 4 }),
                "perplexity {bad} with n = 4 must be rejected"
            );
        }
        // Below 1 is equally impossible: exp(H) >= 1 for any distribution.
        let cfg = TsneConfig { perplexity: 0.5, n_iter: 20, ..Default::default() };
        assert_eq!(
            tsne(&x, &cfg),
            Err(TsneError::PerplexityBelowOne { perplexity: 0.5 })
        );
        // Just below the ceiling stays accepted (max entropy ln 3).
        let cfg = TsneConfig { perplexity: 2.9, n_iter: 20, ..Default::default() };
        assert!(tsne(&x, &cfg).is_ok());
    }

    #[test]
    fn s29_nonconvergent_rows_keep_their_final_beta_weights() {
        // doc_v20 S29: distance scales spanning 1..4e36 with perplexity 3.5
        // (of a maximum 4) need beta ~ 1e-36, which 64 bisection steps from
        // 1.0 cannot reach: no row converges. Pre-fix every row was left
        // all-zero, the affinity matrix vanished, and the dynamics became
        // pure repulsion, which pushes every pair apart at the same rate.
        // Post-fix each row keeps its final-beta weights: the pair with the
        // dominant affinity (3, 4 — each other's only non-underflowed
        // neighbour, p_34 = 0.2 vs 1/15 within {0,1,2} and ~0 across) stays
        // relatively close, while the zero-affinity pair (1, 3) is held apart
        // by repulsion alone. Vanilla t-SNE has no step-size control, so the
        // absolute scale grows without bound on this fixture either way; the
        // regression signal is the distance RATIO, which pure repulsion
        // cannot produce (it equalizes pair distances instead).
        let x = Mat::from_fn(5, 1, |i, _| [-2.0, 0.0, 2.0, 2e18, 2e18 + 2.0][i]);
        let cfg = TsneConfig {
            perplexity: 3.5,
            n_iter: 300,
            seed: 3,
            learning_rate: 20.0,
            ..Default::default()
        };
        let y = tsne(&x, &cfg).expect("perplexity below the ceiling");
        for v in y.data.iter() {
            assert!(v.is_finite(), "coordinate {v} not finite");
        }
        let d01 = (y.get(0, 0) - y.get(1, 0)).abs();
        let d34 = (y.get(3, 0) - y.get(4, 0)).abs();
        let d13 = (y.get(1, 0) - y.get(3, 0)).abs();
        assert!(d13 > 0.0, "degenerate embedding: all distances collapsed");
        // Both attracted pairs stay at least 4x closer than the repulsion-only
        // pair. Measured post-fix: d34/d13 ~= 0.06, d01/d13 ~= 0.06.
        assert!(
            d34 < d13 / 4.0,
            "affinity weights were dropped: d34 = {d34:e} not << d13 = {d13:e}"
        );
        assert!(
            d01 < d13 / 4.0,
            "affinity weights were dropped: d01 = {d01:e} not << d13 = {d13:e}"
        );
    }

    #[test]
    fn s29_converging_case_embedding_is_unchanged() {
        // doc_v20 S29: the always-write change must only add rows for
        // non-converged points; a fully converging run keeps a deterministic
        // embedding. Pinned after the R8 gradient correction and the S29
        // data-scale beta start (doc_v20 review). Both earlier pins
        // (73.6595911798..., 2.2772981587e3) were artifacts of the bugs:
        // the first of the over-normalized gradient, the second of the
        // fixed beta = 1 search start. This test's job is regression
        // detection, not agreement with a buggy baseline.
        let x = blobs(5);
        let cfg = TsneConfig { perplexity: 5.0, n_iter: 300, seed: 3, ..Default::default() };
        let y = tsne(&x, &cfg).expect("converging perplexity");
        let sum_abs: f64 = y.data.iter().map(|v| v.abs()).sum();
        assert!(
            (sum_abs - 2.264348394311807e3).abs() < 1e-6,
            "converged embedding changed: sum|y| = {sum_abs:e}"
        );
        assert!(separation(&y, 5) > 1.5, "t-SNE lost the blob structure");
    }

    #[test]
    fn s31_lle_knn_sort_tolerates_nan_distances() {
        // doc_v20 S31: the neighbourhood sort in lle_from_sq_distances used
        // partial_cmp().unwrap(); one NaN squared distance panicked it.
        let pts = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [3.0, 0.0]];
        let x = Mat::from_fn(4, 2, |i, j| pts[i][j]);
        let mut d2 = Mat::zeros(4, 4);
        for i in 0..4 {
            for j in 0..4 {
                d2.set(i, j, dist2(x.row(i), x.row(j)));
            }
        }
        d2.set(1, 2, f64::NAN);
        d2.set(2, 1, f64::NAN);
        let y = lle_from_sq_distances(&x, &d2, 2, 2);
        for v in y.data.iter() {
            assert!(v.is_finite(), "coordinate {v} not finite");
        }
    }

    #[test]
    fn s29_genotype_scale_distances_do_not_zero_every_row() {
        // doc_v20 review of S29: dosage distances are O(number of sites), so
        // d^2 ~ 1e2..1e6 and `exp(-1 * d^2)` underflows immediately. With
        // the old fixed `beta = 1` start every row broke on `sum == 0` and
        // the embedding was pure repulsion. The search must find a usable
        // beta from the data scale instead.
        let n = 12;
        let x = Mat::from_fn(n, 3, |i, j| {
            // Three groups of four, spread over a realistic dosage range.
            (i / 4) as f64 * 40.0 + (j as f64) * 7.0
        });
        let y = tsne(&x, &TsneConfig { perplexity: 4.0, n_iter: 120, seed: 7, ..Default::default() })
            .expect("genotype-scale fixture");
        for v in y.data.iter() {
            assert!(v.is_finite(), "coordinate {v} not finite");
        }
        let sum_abs: f64 = y.data.iter().map(|v| v.abs()).sum();
        assert!(
            sum_abs > 1e-2,
            "embedding collapsed to the initialization scale (pure repulsion): sum|y| = {sum_abs:e}"
        );
        assert!(separation(&y, 4) > 1.5, "no group structure recovered");
    }

    #[test]
    fn tsne_separates_blobs() {
        let x = blobs(20);
        let cfg = TsneConfig { perplexity: 10.0, n_iter: 300, seed: 3, ..Default::default() };
        let y = tsne(&x, &cfg).expect("converging t-SNE fixture");
        let sep = separation(&y, 20);
        assert!(sep > 2.0, "t-SNE separation ratio {sep}");
    }

    #[test]
    fn r8_tsne_gradient_matches_the_kl_cost() {
        // doc_v20 review of R8: the earlier version of this test compared
        // `tsne_grad` against a line-for-line copy of itself, so it could
        // not fail for any input. The gradient is now checked against a
        // central finite difference of the cost it claims to minimize —
        // an independent oracle. With the wrong (fully normalized) factor
        // the error is O(1e-1); with the right one it is O(1e-10).
        let y = Mat::from_fn(3, 2, |i, d| {
            match (i, d) {
                (0, 0) => 0.0,
                (0, 1) => 1.0,
                (1, 0) => 2.0,
                (1, 1) => -1.0,
                (2, 0) => 0.5,
                _ => 0.25,
            }
        });
        // A joint distribution over ordered pairs, symmetrized and
        // normalized exactly as `tsne_from_sq_distances` builds `pj`.
        let mut raw = Mat::zeros(3, 3);
        let pairs = [(0, 1, 0.20), (0, 2, 0.15), (1, 0, 0.18), (1, 2, 0.10), (2, 0, 0.05), (2, 1, 0.22)];
        for (i, j, p) in pairs {
            raw.set(i, j, p);
        }
        let mut pj = Mat::zeros(3, 3);
        let mut total = 0.0f64;
        for i in 0..3 {
            for j in 0..3 {
                if i != j {
                    let v = (raw.get(i, j) + raw.get(j, i)) / 2.0;
                    pj.set(i, j, v);
                    total += v;
                }
            }
        }
        for v in pj.data.iter_mut() {
            *v /= total;
        }
        // Kernel and normalizer, exactly as the optimizer builds them.
        let mut q = Mat::zeros(3, 3);
        let mut sum_q = 0.0f64;
        for i in 0..3 {
            for j in 0..3 {
                if i == j {
                    continue;
                }
                let dd: f64 = (0..2).map(|d| (y.get(i, d) - y.get(j, d)).powi(2)).sum();
                let v = 1.0 / (1.0 + dd);
                q.set(i, j, v);
                sum_q += v;
            }
        }
        let got = tsne_grad(&pj, &q, sum_q, &y);

        // Central finite difference of the cost.
        let h = 1e-6;
        let mut worst = 0.0f64;
        let mut scale = 0.0f64;
        for i in 0..3 {
            for d in 0..2 {
                let mut yp = y.clone();
                let mut ym = y.clone();
                yp.set(i, d, y.get(i, d) + h);
                ym.set(i, d, y.get(i, d) - h);
                let num = (tsne_kl_cost(&pj, &yp) - tsne_kl_cost(&pj, &ym)) / (2.0 * h);
                worst = worst.max((got.get(i, d) - num).abs());
                scale = scale.max(num.abs());
            }
        }
        assert!(
            worst < 1e-6 * scale.max(1.0),
            "t-SNE gradient is not the derivative of the KL cost: max abs error {worst:.3e}, gradient scale {scale:.3e}"
        );
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
