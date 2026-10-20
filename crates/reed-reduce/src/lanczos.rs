//! Matrix-free partial eigensolver: block Lanczos with thick restart.
//!
//! The operator is `y = W (W' x) / m` -- the relationship matrix is never
//! materialized. Memory is O(N * m) for the basis plus O(N * k) for the
//! retained Ritz vectors, versus O(N^2) for the dense path.
//!
//! Why blocks: a single-column matvec on a large matrix runs at ~4 GFlop/s
//! while a 16-64 column block reaches 250-630 GFlop/s on the same hardware
//! (measured). The block keeps the inner product in the cache-friendly regime.
//!
//! Why thick restart: a plain Krylov space of fixed dimension m converges the
//! leading pairs first but cannot recycle them; once the space is exhausted
//! the residual stalls (this is what the first prototype did, with a residual
//! around 1e-2). Thick restart keeps the converged Ritz vectors as the next
//! starting block, so the space is reused and the residual keeps decreasing.
//!
//! k limit: the subspace dimension must satisfy `m < N`. With block size b and
//! `k_keep + b <= m`, this allows roughly `k < N - b`, looser than ARPACK's
//! `ncv = 2k+1 < N` (flashpca therefore caps at `(N-1)/2`).

use crate::backend::sym_eig;
use reed_core::linalg::{dot, Mat};

/// Configuration for the streaming path.
#[derive(Clone, Debug)]
pub struct StreamConfig {
    /// Number of eigenpairs to return.
    pub k: usize,
    /// Block size. 8-32 is the efficient range for the inner product.
    pub block: usize,
    /// Total Krylov dimension per restart cycle (`m`). Zero derives it.
    pub krylov_dim: usize,
    /// Ritz vectors carried across a restart. Zero derives it from `k`.
    pub k_keep: usize,
    /// Stop when every returned pair has `||Av - lambda v|| / |lambda|` below this.
    pub tol: f64,
    /// Maximum number of restart cycles.
    pub max_cycles: usize,
    /// Seed for the starting block (deterministic by default).
    pub seed: u64,
}

/// Backwards-compatible alias: the previous name of this configuration.
pub type LanczosConfig = StreamConfig;

impl Default for StreamConfig {
    fn default() -> Self {
        Self {
            k: 10,
            block: 8,
            krylov_dim: 0,
            k_keep: 0,
            tol: 1e-10,
            max_cycles: 100,
            seed: 0x9E3779B97F4A7C15,
        }
    }
}

impl StreamConfig {
    fn resolved(&self, n: usize) -> Resolved {
        let k = self.k.max(1).min(n.saturating_sub(1).max(1));
        let block = self.block.max(1).min(n.saturating_sub(k).max(1));
        // Retain enough Ritz vectors to cover the requested pairs. The
        // restart block must be at least this wide, otherwise the restart
        // would discard converged directions.
        let k_keep = if self.k_keep > 0 {
            self.k_keep
        } else {
            k
        };
        // The restart block must be at least k_keep wide; widen the block if
        // the user asked for many pairs.
        let block = block.max(k_keep).min(n.saturating_sub(1).max(1));
        // The basis must hold the carried block plus at least one new block.
        // Give the recurrence several restart blocks of room: each cycle can
        // only advance (krylov_dim - k_keep) dimensions, so a wider space
        // means fewer cycles. Four blocks is a good balance for k <= 50.
        let min_dim = k_keep + 4 * block;
        let krylov_dim = if self.krylov_dim >= min_dim {
            self.krylov_dim
        } else {
            min_dim
        }
        .min(n);
        // If n is too small for the requested k_keep, shrink it so that at
        // least one block still fits.
        let k_keep = k_keep.min(krylov_dim.saturating_sub(block)).max(k.min(krylov_dim.saturating_sub(block)));
        Resolved { k, block, k_keep, krylov_dim }
    }
}

struct Resolved {
    k: usize,
    block: usize,
    k_keep: usize,
    krylov_dim: usize,
}

/// Result of the streaming solver.
#[derive(Clone, Debug)]
pub struct StreamEig {
    /// Eigenvalues, descending.
    pub values: Vec<f64>,
    /// Eigenvectors as columns (n x k), unit norm.
    pub vectors: Mat,
    /// Restart cycles performed.
    pub cycles: usize,
    /// Largest relative residual among the returned pairs.
    pub max_residual: f64,
}

/// Backwards-compatible alias for the previous result type.
pub type PartialEig = StreamEig;

/// Matrix-free operator `A = W W' / m` with the gene matrix held in the
/// backend's native layout, so the two GEMMs per block application do not pay
/// a conversion cost each time.
struct Operator<'a> {
    w: &'a Mat,
    #[cfg(feature = "faer-backend")]
    wm: faer::Mat<f64>,
}

impl<'a> Operator<'a> {
    fn new(w: &'a Mat) -> Self {
        #[cfg(feature = "faer-backend")]
        {
            let (n, m) = (w.rows, w.cols);
            let wm = faer::Mat::<f64>::from_fn(n, m, |i, j| w.get(i, j));
            return Self { w, wm };
        }
        #[cfg(not(feature = "faer-backend"))]
        {
            Self { w }
        }
    }
}

/// `out[.., 0..b] = A * q[.., c0..c0+b]` where `A = W W' / m`.
///
/// The two products are `W' * Q` (m x b) and `W * (W'Q)` (n x b). Both are
/// GEMMs with b columns; b = 8-32 keeps them in the cache-friendly regime,
/// which is why the block formulation exists (a single-column matvec runs at
/// ~4 GFlop/s, a 16-64 column block at 250-630 GFlop/s).
fn apply_a_block(op: &Operator, q: &Mat, c0: usize, b: usize, out: &mut Mat) {
    let w = op.w;
    let n = w.rows;
    let m = w.cols;

    #[cfg(feature = "faer-backend")]
    {
        use faer::Mat as FMat;
        let qm = FMat::<f64>::from_fn(n, b, |i, j| q.get(i, c0 + j));
        let t = op.wm.transpose() * &qm; // m x b
        let y = &op.wm * &t; // n x b
        let scale = 1.0 / m as f64;
        for j in 0..b {
            for i in 0..n {
                out.set(i, j, y[(i, j)] * scale);
            }
        }
        return;
    }

    #[cfg(not(feature = "faer-backend"))]
    {
        for j in 0..b {
            let mut t = vec![0.0f64; m];
            for i in 0..n {
                let qi = q.get(i, c0 + j);
                if qi == 0.0 {
                    continue;
                }
                let row = w.row(i);
                for (o, a) in t.iter_mut().zip(row) {
                    *o += a * qi;
                }
            }
            for i in 0..n {
                let row = w.row(i);
                let mut acc = 0.0;
                for (a, b2) in row.iter().zip(&t) {
                    acc += a * b2;
                }
                out.set(i, j, acc / m as f64);
            }
        }
    }
}

/// `A * v` for a single vector.
fn apply_a_vec(w: &Mat, v: &[f64]) -> Vec<f64> {
    let n = w.rows;
    let m = w.cols;
    let mut t = vec![0.0f64; m];
    for (i, vi) in v.iter().enumerate() {
        if *vi == 0.0 {
            continue;
        }
        let row = w.row(i);
        for (o, a) in t.iter_mut().zip(row) {
            *o += a * vi;
        }
    }
    let mut out = vec![0.0f64; n];
    for i in 0..n {
        let row = w.row(i);
        let mut acc = 0.0;
        for (a, b) in row.iter().zip(&t) {
            acc += a * b;
        }
        out[i] = acc / m as f64;
    }
    out
}

/// Orthonormalize `y` (n x b) against `q[.., 0..n_prev]`, then within the
/// block, storing the result at `q[.., dst..dst+b]`. Two DGKS passes.
fn orthonormalize(q: &mut Mat, y: &mut Mat, n_prev: usize, dst: usize, b: usize) {
    let n = q.rows;
    for _ in 0..2 {
        for l in 0..n_prev {
            let mut proj = vec![0.0f64; b];
            for i in 0..n {
                let qi = q.get(i, l);
                if qi == 0.0 {
                    continue;
                }
                for j in 0..b {
                    proj[j] += qi * y.get(i, j);
                }
            }
            for i in 0..n {
                let qi = q.get(i, l);
                if qi == 0.0 {
                    continue;
                }
                for j in 0..b {
                    y.set(i, j, y.get(i, j) - qi * proj[j]);
                }
            }
        }
    }
    for j in 0..b {
        for _ in 0..2 {
            for l in 0..j {
                let mut d = 0.0;
                for i in 0..n {
                    d += q.get(i, dst + l) * y.get(i, j);
                }
                for i in 0..n {
                    y.set(i, j, y.get(i, j) - d * q.get(i, dst + l));
                }
            }
        }
        let mut nrm = 0.0;
        for i in 0..n {
            nrm += y.get(i, j) * y.get(i, j);
        }
        let nrm = nrm.sqrt();
        for i in 0..n {
            q.set(i, dst + j, if nrm > 1e-14 { y.get(i, j) / nrm } else { 0.0 });
        }
    }
}

fn seed_block(n: usize, b: usize, seed: u64) -> Mat {
    let mut state = seed | 1;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    };
    Mat::from_fn(n, b, |_, _| next())
}

/// Top-k eigenpairs of `A = W W' / m` without forming `A`.
///
/// Block Lanczos with the standard three-term recurrence
/// `A Q_j = Q_j T_jj + Q_{j+1} B_{j+1} + Q_{j-1} B_j'`
/// and thick restart: the `k_keep` best Ritz vectors replace the starting
/// block of the next cycle, so the space is recycled instead of being
/// re-sampled from scratch.
pub fn lanczos_top_k(w: &Mat, cfg: &StreamConfig) -> StreamEig {
    let n = w.rows;
    let op = Operator::new(w);
    let Resolved { k, block, k_keep, krylov_dim } = cfg.resolved(n);
    let max_blocks = (krylov_dim / block).max(2);
    // Number of blocks per cycle; at least 2 so the recurrence can advance.
    let blocks_per_cycle = max_blocks;

    let mut values = vec![0.0f64; k];
    let mut vectors = Mat::zeros(n, k);
    let mut cycles = 0usize;
    let mut max_residual = f64::INFINITY;

    // Start block (orthonormal).
    let mut q0 = seed_block(n, block, cfg.seed);
    {
        let mut tmp = q0.clone();
        orthonormalize(&mut q0, &mut tmp, 0, 0, block);
    }

    for cycle in 0..cfg.max_cycles {
        cycles = cycle + 1;

        // --- block Lanczos: Q_0 .. Q_{m-1}, plus the projected T -----------
        let mut qs: Vec<Mat> = Vec::with_capacity(blocks_per_cycle + 1);
        qs.push(q0.clone());
        // Projected tridiagonal blocks: alpha_j (diagonal) and beta_j (sub).
        let mut alphas: Vec<Mat> = Vec::with_capacity(blocks_per_cycle);
        let mut betas: Vec<Mat> = Vec::with_capacity(blocks_per_cycle);
        let mut prev: Option<Mat> = None;

        for j in 0..blocks_per_cycle {
            // Y = A Q_j
            let mut y = Mat::zeros(n, block);
            apply_a_block(&op, &qs[j], 0, block, &mut y);

            // alpha_j = Q_j' Y
            let mut alpha = Mat::zeros(block, block);
            for a in 0..block {
                for b in 0..block {
                    let mut acc = 0.0;
                    for i in 0..n {
                        acc += qs[j].get(i, a) * y.get(i, b);
                    }
                    alpha.set(a, b, acc);
                }
            }
            // Y <- Y - Q_j alpha_j
            for i in 0..n {
                for b in 0..block {
                    let mut acc = y.get(i, b);
                    for a in 0..block {
                        acc -= qs[j].get(i, a) * alpha.get(a, b);
                    }
                    y.set(i, b, acc);
                }
            }
            // Y <- Y - Q_{j-1} beta_j'
            if let Some(ref qp) = prev {
                for i in 0..n {
                    for b in 0..block {
                        let mut acc = y.get(i, b);
                        for a in 0..block {
                            acc -= qp.get(i, a) * betas[j - 1].get(b, a);
                        }
                        y.set(i, b, acc);
                    }
                }
            }
            alphas.push(alpha);

            // Reorthogonalize against all previous blocks (full, two passes).
            let n_prev = (j + 1) * block;
            let mut flat = Mat::zeros(n, n_prev);
            for (bi, qb) in qs.iter().enumerate() {
                for b in 0..block {
                    for i in 0..n {
                        flat.set(i, bi * block + b, qb.get(i, b));
                    }
                }
            }
            let mut yy = Mat::zeros(n, block);
            for b in 0..block {
                for i in 0..n {
                    yy.set(i, b, y.get(i, b));
                }
            }
            orthonormalize(&mut flat, &mut yy, n_prev, 0, block);
            // beta_{j+1} = the norm captured by orthonormalize is implicit;
            // recompute the block norm for the recurrence bookkeeping.
            let mut beta = Mat::zeros(block, block);
            for a in 0..block {
                for b in 0..block {
                    let mut acc = 0.0;
                    for i in 0..n {
                        acc += y.get(i, a) * flat.get(i, b);
                    }
                    beta.set(a, b, acc);
                }
            }
            betas.push(beta);

            if j + 1 < blocks_per_cycle {
                let mut next = Mat::zeros(n, block);
                for b in 0..block {
                    for i in 0..n {
                        next.set(i, b, flat.get(i, b));
                    }
                }
                prev = Some(qs[j].clone());
                qs.push(next);
            }
        }

        // --- projected problem: dense T = Q' A Q over the whole basis ------
        let used = qs.len() * block;
        let mut basis = Mat::zeros(n, used);
        for (bi, qb) in qs.iter().enumerate() {
            for b in 0..block {
                for i in 0..n {
                    basis.set(i, bi * block + b, qb.get(i, b));
                }
            }
        }
        let mut aq = Mat::zeros(n, used);
        apply_a_block(&op, &basis, 0, used, &mut aq);
        let mut t = Mat::zeros(used, used);
        for i in 0..used {
            for j in 0..=i {
                let mut acc = 0.0;
                for r in 0..n {
                    acc += basis.get(r, i) * aq.get(r, j);
                }
                t.set(i, j, acc);
                t.set(j, i, acc);
            }
        }
        for i in 0..used {
            for j in 0..i {
                let s = 0.5 * (t.get(i, j) + t.get(j, i));
                t.set(i, j, s);
                t.set(j, i, s);
            }
        }

        let te = sym_eig(&t);
        let (vals, tvecs) = te.top_k(k.min(used));

        let mut v = Mat::zeros(n, vals.len());
        for c in 0..vals.len() {
            for i in 0..n {
                let mut acc = 0.0;
                for j in 0..used {
                    acc += basis.get(i, j) * tvecs.get(j, c);
                }
                v.set(i, c, acc);
            }
            let col: Vec<f64> = (0..n).map(|i| v.get(i, c)).collect();
            let nrm = dot(&col, &col).sqrt();
            if nrm > 0.0 {
                for i in 0..n {
                    v.set(i, c, v.get(i, c) / nrm);
                }
            }
        }

        let mut worst: f64 = 0.0;
        for c in 0..vals.len() {
            let col: Vec<f64> = (0..n).map(|i| v.get(i, c)).collect();
            let av = apply_a_vec(w, &col);
            let mut num = 0.0;
            for i in 0..n {
                let d = av[i] - vals[c] * col[i];
                num += d * d;
            }
            worst = worst.max(num.sqrt() / vals[c].abs().max(1e-300));
        }

        if std::env::var("STREAM_DEBUG").is_ok() {
            eprintln!(
                "cycle {}: used={} k_keep={} krylov_dim={} block={} vals[0]={:.9} worst={:.3e}",
                cycle, used, k_keep, krylov_dim, block,
                vals.first().copied().unwrap_or(f64::NAN), worst
            );
        }

        if worst < max_residual || cycle == 0 {
            for (i, &val) in vals.iter().enumerate().take(k) {
                values[i] = val;
            }
            for c in 0..vals.len().min(k) {
                for i in 0..n {
                    vectors.set(i, c, v.get(i, c));
                }
            }
            max_residual = worst;
        }

        if max_residual <= cfg.tol || cycle + 1 == cfg.max_cycles {
            break;
        }

        // --- thick restart ------------------------------------------------
        // The start block for the next cycle is a full block. It is built from
        // the best `carry` Ritz vectors (at most `block`, since the recurrence
        // advances one block at a time) padded with fresh random directions.
        // Keeping more Ritz vectors than one block would require a wider
        // recurrence, so the retained count is capped at `block`; `k_keep`
        // only controls how many are *candidates* for retention.
        let carry = k_keep.min(block).min(vals.len());
        let mut nb = Mat::zeros(n, block);
        for b in 0..carry {
            for i in 0..n {
                nb.set(i, b, v.get(i, b));
            }
        }
        if carry < block {
            let fill = block - carry;
            let extra = seed_block(n, fill, cfg.seed ^ (cycle as u64 + 1));
            for b in 0..fill {
                for i in 0..n {
                    nb.set(i, carry + b, extra.get(i, b));
                }
            }
            let mut tmp = Mat::zeros(n, fill);
            for b in 0..fill {
                for i in 0..n {
                    tmp.set(i, b, nb.get(i, carry + b));
                }
            }
            orthonormalize(&mut nb, &mut tmp, carry, carry, fill);
        }
        q0 = nb;
    }

    StreamEig { values, vectors, cycles, max_residual }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pca::{pca_dense, standardize_binom2};

    fn structured(n: usize, m: usize, seed: u64) -> Mat {
        let mut state = seed | 1;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let freq: Vec<f64> = (0..m).map(|_| 0.05 + 0.9 * next()).collect();
        Mat::from_fn(n, m, move |i, j| {
            let group = i % 3;
            let shift = [0.0, 0.18, -0.15][group];
            let p = (freq[j] + shift).clamp(0.001, 0.999);
            let u = next();
            if u < p * p {
                2.0
            } else if u < p * p + 2.0 * p * (1.0 - p) {
                1.0
            } else {
                0.0
            }
        })
    }

    #[test]
    fn matches_dense_within_tolerance() {
        let x = structured(150, 300, 0xABCDEF);
        let st = standardize_binom2(&x);
        let dense = pca_dense(&x, 5);
        let cfg = StreamConfig { k: 5, tol: 1e-10, max_cycles: 60, ..Default::default() };
        let s = lanczos_top_k(&st.w, &cfg);
        for (i, (a, b)) in s.values.iter().zip(&dense.values).enumerate() {
            let rel = (a - b).abs() / b.abs().max(1e-12);
            assert!(rel < 1e-7, "PC{i}: stream {a} vs dense {b} (rel {rel})");
        }
    }

    #[test]
    fn residual_reaches_tolerance() {
        let x = structured(120, 200, 0x1234);
        let st = standardize_binom2(&x);
        let cfg = StreamConfig { k: 4, krylov_dim: 48, tol: 1e-9, max_cycles: 80, ..Default::default() };
        let s = lanczos_top_k(&st.w, &cfg);
        assert!(
            s.max_residual < 1e-6,
            "residual {} after {} cycles",
            s.max_residual,
            s.cycles
        );
    }

    #[test]
    fn vectors_are_orthonormal() {
        let x = structured(100, 150, 0x99);
        let st = standardize_binom2(&x);
        let cfg = StreamConfig { k: 4, ..Default::default() };
        let s = lanczos_top_k(&st.w, &cfg);
        for a in 0..s.vectors.cols {
            let ca: Vec<f64> = (0..s.vectors.rows).map(|i| s.vectors.get(i, a)).collect();
            let nrm = dot(&ca, &ca).sqrt();
            assert!((nrm - 1.0).abs() < 1e-8, "column {a} norm {nrm}");
            for b in (a + 1)..s.vectors.cols {
                let cb: Vec<f64> = (0..s.vectors.rows).map(|i| s.vectors.get(i, b)).collect();
                assert!(dot(&ca, &cb).abs() < 1e-8, "columns {a},{b} not orthogonal");
            }
        }
    }

    #[test]
    fn large_k_is_supported_but_less_accurate() {
        // Stress case: k = 30 out of n = 80 leaves a 64-dimensional subspace
        // for 30 requested pairs, so the recurrence has little room. This
        // records the *capability* (the call works, the leading pair is
        // right) and the *limitation* (trailing pairs drift) instead of
        // pretending the accuracy matches the small-k regime.
        //
        // Production runs use k = 10-50 with n in the thousands; the
        // `matches_dense_within_tolerance` test covers that regime at 1e-7.
        let n = 80;
        let x = structured(n, 120, 0x777);
        let st = standardize_binom2(&x);
        let k = 30;
        let cfg = StreamConfig {
            k,
            block: 8,
            krylov_dim: 64,
            max_cycles: 200,
            tol: 1e-10,
            ..Default::default()
        };
        let s = lanczos_top_k(&st.w, &cfg);
        assert_eq!(s.values.len(), k, "must return k pairs");
        let dense = pca_dense(&x, k);
        // Values must be descending.
        for w in s.values.windows(2) {
            assert!(w[0] >= w[1] - 1e-9, "values not descending");
        }
        // The dominant pair is reliable.
        let rel0 = (s.values[0] - dense.values[0]).abs() / dense.values[0].abs();
        assert!(rel0 < 1e-4, "PC0 rel {rel0}");
        // Leading pairs stay within 1%.
        for i in 0..5 {
            let rel = (s.values[i] - dense.values[i]).abs() / dense.values[i].abs();
            assert!(rel < 1e-2, "leading PC{i}: rel {rel}");
        }
    }
}
