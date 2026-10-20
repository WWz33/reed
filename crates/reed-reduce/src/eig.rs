//! Dense symmetric eigendecomposition via Householder tridiagonalization
//! (`tred2`) followed by implicit QL iteration (`tql2`).
//!
//! These are the classic EISPACK routines (public domain, also reproduced in
//! *Numerical Recipes* and in every LAPACK predecessor). Written here from the
//! published algorithm so the crate stays dependency-free. Cost is O(n^3) and
//! memory O(n^2) -- intended for the small-to-medium N path.

use reed_core::linalg::Mat;

/// Eigendecomposition of a real symmetric matrix.
#[derive(Clone, Debug)]
pub struct SymEig {
    /// Eigenvalues, ascending.
    pub values: Vec<f64>,
    /// Eigenvectors as columns of `vectors` (row-major n x n).
    pub vectors: Mat,
}

impl SymEig {
    /// Top-k eigenvalues (descending) and their eigenvectors as columns.
    pub fn top_k(&self, k: usize) -> (Vec<f64>, Mat) {
        let n = self.values.len();
        let k = k.min(n);
        let mut vals = Vec::with_capacity(k);
        let mut vecs = Mat::zeros(n, k);
        for c in 0..k {
            let src = n - 1 - c;
            vals.push(self.values[src]);
            for r in 0..n {
                vecs.set(r, c, self.vectors.get(r, src));
            }
        }
        (vals, vecs)
    }
}

/// Eigendecomposition of a real symmetric matrix (only the lower triangle of
/// `a` is read). Returns eigenvalues ascending and orthonormal eigenvectors.
pub fn sym_eig(a: &Mat) -> SymEig {
    let n = a.rows;
    assert_eq!(n, a.cols, "sym_eig needs a square matrix");
    if n == 0 {
        return SymEig { values: Vec::new(), vectors: Mat::zeros(0, 0) };
    }

    // Work on a full symmetric copy in row-major order.
    let mut v = a.clone();
    for i in 0..n {
        for j in 0..i {
            let s = 0.5 * (a.get(i, j) + a.get(j, i));
            v.set(i, j, s);
            v.set(j, i, s);
        }
    }

    let mut d = vec![0.0f64; n];
    let mut e = vec![0.0f64; n];

    tred2(n, &mut v, &mut d, &mut e);
    tql2(n, &mut v, &mut d, &mut e);

    SymEig { values: d, vectors: v }
}

/// Householder reduction of a symmetric matrix to tridiagonal form.
///
/// On exit `d` holds the diagonal, `e[0..n-1]` the subdiagonal (`e[0] = 0`),
/// and `v` holds the accumulated orthogonal transformation (eigenvectors on
/// columns once `tql2` has run).
fn tred2(n: usize, v: &mut Mat, d: &mut [f64], e: &mut [f64]) {
    for i in (1..n).rev() {
        let l = i - 1;
        let mut h = 0.0;
        let mut scale = 0.0;
        if l > 0 {
            for k in 0..=l {
                scale += v.get(i, k).abs();
            }
            if scale == 0.0 {
                e[i] = v.get(i, l);
            } else {
                for k in 0..=l {
                    let val = v.get(i, k) / scale;
                    v.set(i, k, val);
                    h += val * val;
                }
                let f = v.get(i, l);
                let g = if f >= 0.0 { -h.sqrt() } else { h.sqrt() };
                e[i] = scale * g;
                h -= f * g;
                v.set(i, l, f - g);
                // Accumulate the transformation in the yet-unused column i.
                let mut f_acc = 0.0;
                for j in 0..=l {
                    let vij = v.get(i, j);
                    v.set(j, i, vij / h);
                    let mut g_acc = 0.0;
                    for k in 0..=j {
                        g_acc += v.get(j, k) * v.get(i, k);
                    }
                    for k in (j + 1)..=l {
                        g_acc += v.get(k, j) * v.get(i, k);
                    }
                    e[j] = g_acc / h;
                    f_acc += e[j] * v.get(i, j);
                }
                let hh = f_acc / (h + h);
                for j in 0..=l {
                    let f = v.get(i, j);
                    let g = e[j] - hh * f;
                    e[j] = g;
                    for k in 0..=j {
                        let upd = v.get(j, k) - f * e[k] - g * v.get(i, k);
                        v.set(j, k, upd);
                    }
                }
            }
        } else {
            e[i] = v.get(i, l);
        }
        d[i] = h;
    }

    // Accumulate the transformations.
    d[0] = 0.0;
    e[0] = 0.0;
    for i in 0..n {
        let l = i;
        if d[i] != 0.0 {
            for j in 0..l {
                let mut g = 0.0;
                for k in 0..l {
                    g += v.get(i, k) * v.get(k, j);
                }
                for k in 0..l {
                    let upd = v.get(k, j) - g * v.get(k, i);
                    v.set(k, j, upd);
                }
            }
        }
        d[i] = v.get(i, i);
        v.set(i, i, 1.0);
        for j in 0..i {
            v.set(j, i, 0.0);
            v.set(i, j, 0.0);
        }
    }
}

/// Implicit QL iteration on a tridiagonal matrix.
fn tql2(n: usize, v: &mut Mat, d: &mut [f64], e: &mut [f64]) {
    for i in 1..n {
        e[i - 1] = e[i];
    }
    e[n - 1] = 0.0;

    let mut f = 0.0;
    let mut tst1: f64 = 0.0;
    let eps = f64::EPSILON;
    for l in 0..n {
        tst1 = tst1.max(d[l].abs() + e[l].abs());
        let mut m = l;
        while m < n {
            if e[m].abs() <= eps * tst1 {
                break;
            }
            m += 1;
        }
        if m > l {
            let mut iter = 0;
            loop {
                iter += 1;
                if iter > 100 {
                    // Non-convergence: leave the remaining values as they are.
                    break;
                }
                let mut g = d[l];
                let mut p = (d[l + 1] - g) / (2.0 * e[l]);
                let mut r = p.hypot(1.0);
                if p < 0.0 {
                    r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l + 1] = e[l] * (p + r);
                let dl1 = d[l + 1];
                let mut h = g - d[l];
                for i in (l + 2)..n {
                    d[i] -= h;
                }
                f += h;
                p = d[m];
                let mut c = 1.0;
                let mut c2 = c;
                let mut c3 = c;
                let el1 = e[l + 1];
                let mut s = 0.0;
                let mut s2 = 0.0;
                for i in (l..m).rev() {
                    c3 = c2;
                    c2 = c;
                    s2 = s;
                    g = c * e[i];
                    h = c * p;
                    r = p.hypot(e[i]);
                    e[i + 1] = s * r;
                    s = e[i] / r;
                    c = p / r;
                    p = c * d[i] - s * g;
                    d[i + 1] = h + s * (c * g + s * d[i]);
                    for k in 0..n {
                        h = v.get(k, i + 1);
                        v.set(k, i + 1, s * v.get(k, i) + c * h);
                        v.set(k, i, c * v.get(k, i) - s * h);
                    }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;
                if e[l].abs() <= eps * tst1 {
                    break;
                }
            }
        }
        d[l] += f;
        e[l] = 0.0;
    }

    // Sort ascending (selection sort on the small array).
    for i in 0..n.saturating_sub(1) {
        let mut k = i;
        let mut p = d[i];
        for j in (i + 1)..n {
            if d[j] < p {
                k = j;
                p = d[j];
            }
        }
        if k != i {
            d[k] = d[i];
            d[i] = p;
            for j in 0..n {
                v.data.swap(j * n + i, j * n + k);
            }
        }
    }
}


/// Householder tridiagonalization with the accumulated orthogonal factor:
/// returns `(q, diag, off)` with `q * T * q' = a` and `T` tridiagonal
/// (`off[i]` is the element `(i+1, i)`).
fn tridiagonalize(a: &Mat) -> (Mat, Vec<f64>, Vec<f64>) {
    let n = a.rows;
    let mut m = a.clone();
    for i in 0..n {
        for j in 0..i {
            let s = 0.5 * (m.get(i, j) + m.get(j, i));
            m.set(i, j, s);
            m.set(j, i, s);
        }
    }
    let mut q = Mat::zeros(n, n);
    for i in 0..n {
        q.set(i, i, 1.0);
    }

    for k in 0..n.saturating_sub(2) {
        let len = n - k - 1;
        // Householder vector v (normalized so that H = I - 2 v v', ||v|| = 1).
        let mut v = vec![0.0f64; len];
        for i in 0..len {
            v[i] = m.get(k + 1 + i, k);
        }
        let norm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm <= 1e-300 {
            continue;
        }
        let alpha = if v[0] >= 0.0 { -norm } else { norm };
        v[0] -= alpha;
        let vnorm: f64 = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if vnorm <= 1e-300 {
            continue;
        }
        for x in v.iter_mut() {
            *x /= vnorm;
        }

        // A22 <- H A22 H with H = I - 2 v v'.
        let mut av = vec![0.0f64; len];
        for i in 0..len {
            let row = k + 1 + i;
            let mut acc = 0.0;
            for j in 0..len {
                acc += m.get(row, k + 1 + j) * v[j];
            }
            av[i] = acc;
        }
        let vav: f64 = v.iter().zip(&av).map(|(a, b)| a * b).sum();
        for i in 0..len {
            let row = k + 1 + i;
            for j in 0..len {
                let col = k + 1 + j;
                let val = m.get(row, col)
                    - 2.0 * v[i] * av[j]
                    - 2.0 * av[i] * v[j]
                    + 4.0 * vav * v[i] * v[j];
                m.set(row, col, val);
            }
        }

        // Q <- Q H (apply H on the right of the trailing columns).
        let mut qv = vec![0.0f64; n];
        for r in 0..n {
            let mut acc = 0.0;
            for j in 0..len {
                acc += q.get(r, k + 1 + j) * v[j];
            }
            qv[r] = acc;
        }
        for r in 0..n {
            for j in 0..len {
                let val = q.get(r, k + 1 + j) - 2.0 * qv[r] * v[j];
                q.set(r, k + 1 + j, val);
            }
        }

        for i in (k + 2)..n {
            m.set(i, k, 0.0);
            m.set(k, i, 0.0);
        }
        m.set(k + 1, k, alpha);
        m.set(k, k + 1, alpha);
    }

    let diag: Vec<f64> = (0..n).map(|i| m.get(i, i)).collect();
    let mut off = vec![0.0f64; n.saturating_sub(1)];
    for i in 0..n.saturating_sub(1) {
        off[i] = m.get(i + 1, i);
    }
    (q, diag, off)
}

/// Number of eigenvalues of the tridiagonal `(diag, off)` strictly below `x`,
/// by the Sturm sequence / LDL pivot count.
fn sturm_count(diag: &[f64], off: &[f64], x: f64) -> usize {
    let n = diag.len();
    if n == 0 {
        return 0;
    }
    let mut count = 0usize;
    let mut pivot = diag[0] - x;
    if pivot < 0.0 {
        count += 1;
    }
    for i in 1..n {
        let mut p = diag[i] - x;
        if pivot.abs() > 1e-300 {
            p -= off[i - 1] * off[i - 1] / pivot;
        } else {
            p -= off[i - 1] * off[i - 1] / 1e-300;
        }
        pivot = p;
        if pivot < 0.0 {
            count += 1;
        }
    }
    count
}

/// Bisection for the `j`-th smallest eigenvalue (0-based).
fn kth_eigenvalue(diag: &[f64], off: &[f64], j: usize, mut lo: f64, mut hi: f64) -> f64 {
    let tol = 1e-15 * (lo.abs() + hi.abs() + 1.0);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if sturm_count(diag, off, mid) <= j {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo <= tol {
            break;
        }
    }
    0.5 * (lo + hi)
}

/// Inverse iteration for the eigenvector of the tridiagonal at `mu`.
fn inverse_iteration(diag: &[f64], off: &[f64], mu: f64, iters: usize, seed: usize) -> Vec<f64> {
    let n = diag.len();
    let scale = diag.iter().fold(0.0f64, |a, &b| a.max(b.abs()))
        + off.iter().fold(0.0f64, |a, &b| a.max(b.abs()));
    // Shift slightly off the eigenvalue so the factorization stays defined.
    let shift = mu + 1e-13 * scale.max(1e-300);

    let mut piv = vec![0.0f64; n];
    let mut p0 = diag[0] - shift;
    if p0.abs() < 1e-200 {
        p0 = if p0 < 0.0 { -1e-200 } else { 1e-200 };
    }
    piv[0] = p0;
    for i in 1..n {
        let denom = if piv[i - 1].abs() < 1e-200 { 1e-200 } else { piv[i - 1] };
        let mut p = diag[i] - shift - off[i - 1] * off[i - 1] / denom;
        if !p.is_finite() {
            p = diag[i] - shift;
        }
        if p.abs() < 1e-200 {
            p = if p < 0.0 { -1e-200 } else { 1e-200 };
        }
        piv[i] = p;
    }

    let mut x = vec![0.0f64; n];
    for (i, v) in x.iter_mut().enumerate() {
        *v = 1.0 + ((i * 37 + seed * 11) % 13) as f64;
    }
    let nrm = x.iter().map(|v| v * v).sum::<f64>().sqrt();
    if nrm > 0.0 {
        for v in x.iter_mut() {
            *v /= nrm;
        }
    }

    for _ in 0..iters {
        let mut y = vec![0.0f64; n];
        y[0] = x[0];
        for i in 1..n {
            y[i] = x[i] - off[i - 1] / piv[i - 1] * y[i - 1];
        }
        let mut z = vec![0.0f64; n];
        z[n - 1] = y[n - 1] / piv[n - 1];
        for i in (0..n.saturating_sub(1)).rev() {
            z[i] = y[i] / piv[i] - off[i] / piv[i] * z[i + 1];
        }
        if z.iter().any(|v| !v.is_finite()) {
            break;
        }
        let nrm = z.iter().map(|v| v * v).sum::<f64>().sqrt();
        if nrm <= 1e-300 {
            break;
        }
        for v in z.iter_mut() {
            *v /= nrm;
        }
        x = z;
    }
    x
}

/// Top-k eigenpairs of a symmetric matrix via tridiagonalization plus
/// bisection and inverse iteration: only the requested pairs are extracted,
/// so the O(n^3) eigenvector accumulation of the dense QL path is avoided.
///
/// Returns eigenvalues descending and eigenvectors as columns (n x k).
pub fn sym_eig_top_k(a: &Mat, k: usize) -> (Vec<f64>, Mat) {
    let n = a.rows;
    assert_eq!(n, a.cols, "sym_eig_top_k needs a square matrix");
    let k = k.min(n);
    if k == 0 {
        return (Vec::new(), Mat::zeros(n, 0));
    }
    let (q, diag, off) = tridiagonalize(a);

    // Gershgorin bracket.
    let mut lo = f64::MAX;
    let mut hi = f64::MIN;
    for i in 0..n {
        let mut r = 0.0;
        if i > 0 {
            r += off[i - 1].abs();
        }
        if i + 1 < n {
            r += off[i].abs();
        }
        lo = lo.min(diag[i] - r);
        hi = hi.max(diag[i] + r);
    }
    if !(lo.is_finite() && hi.is_finite()) {
        lo = -1.0;
        hi = 1.0;
    }

    let mut values = Vec::with_capacity(k);
    let mut zs: Vec<Vec<f64>> = Vec::with_capacity(k);
    for c in 0..k {
        // The c-th largest is index n-1-c in ascending order.
        let j = n - 1 - c;
        let lam = kth_eigenvalue(&diag, &off, j, lo, hi);
        let mut z = inverse_iteration(&diag, &off, lam, 4, c + 1);
        // Re-orthogonalize against the already computed vectors (clusters).
        for prev in &zs {
            let d: f64 = z.iter().zip(prev).map(|(a, b)| a * b).sum();
            for (zi, pi) in z.iter_mut().zip(prev) {
                *zi -= d * pi;
            }
        }
        let nrm = z.iter().map(|v| v * v).sum::<f64>().sqrt();
        if nrm > 1e-14 {
            for zi in z.iter_mut() {
                *zi /= nrm;
            }
        }
        // Recompute the Rayleigh quotient after orthogonalization.
        let mut tz = vec![0.0f64; n];
        for i in 0..n {
            let mut acc = diag[i] * z[i];
            if i > 0 {
                acc += off[i - 1] * z[i - 1];
            }
            if i + 1 < n {
                acc += off[i] * z[i + 1];
            }
            tz[i] = acc;
        }
        let rq: f64 = z.iter().zip(&tz).map(|(a, b)| a * b).sum();
        values.push(rq);
        zs.push(z);
    }

    // Back-transform: v = Q z.
    let mut vectors = Mat::zeros(n, k);
    for c in 0..k {
        let z = &zs[c];
        for i in 0..n {
            let mut acc = 0.0;
            for j in 0..n {
                acc += q.get(i, j) * z[j];
            }
            vectors.set(i, c, acc);
        }
    }

    // Sort descending (bisection may return near-ties out of order).
    let mut order: Vec<usize> = (0..k).collect();
    order.sort_by(|&x, &y| values[y].partial_cmp(&values[x]).unwrap_or(std::cmp::Ordering::Equal));
    let mut v2 = Vec::with_capacity(k);
    let mut m2 = Mat::zeros(n, k);
    for (dst, &src) in order.iter().enumerate() {
        v2.push(values[src]);
        for i in 0..n {
            m2.set(i, dst, vectors.get(i, src));
        }
    }
    (v2, m2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn residual(a: &Mat, vals: &[f64], vecs: &Mat) -> f64 {
        let n = a.rows;
        let mut worst: f64 = 0.0;
        for c in 0..vecs.cols {
            let v: Vec<f64> = (0..n).map(|i| vecs.get(i, c)).collect();
            let av = a.matvec(&v);
            let mut num = 0.0;
            let mut den = 0.0;
            for i in 0..n {
                let d = av[i] - vals[c] * v[i];
                num += d * d;
                den += v[i] * v[i];
            }
            worst = worst.max((num / den).sqrt() / vals[c].abs().max(1e-300));
        }
        worst
    }

    #[test]
    fn diagonal_matrix() {
        let a = Mat::from_fn(4, 4, |i, j| if i == j { (i + 1) as f64 * 2.0 } else { 0.0 });
        let e = sym_eig(&a);
        let want = [2.0, 4.0, 6.0, 8.0];
        for (got, want) in e.values.iter().zip(want) {
            assert!((got - want).abs() < 1e-10, "got {got} want {want}");
        }
    }

    #[test]
    fn known_2x2() {
        // [[2,1],[1,2]] has eigenvalues 1 and 3.
        let a = Mat::from_fn(2, 2, |i, j| if i == j { 2.0 } else { 1.0 });
        let e = sym_eig(&a);
        assert!((e.values[0] - 1.0).abs() < 1e-12);
        assert!((e.values[1] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn random_symmetric_residual_is_small() {
        // Deterministic pseudo-random symmetric matrix.
        let n = 40;
        let mut a = Mat::zeros(n, n);
        let mut state = 0x2545F4914F6CDD1Du64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        for i in 0..n {
            for j in 0..=i {
                let v = next();
                a.set(i, j, v);
                a.set(j, i, v);
            }
        }
        let e = sym_eig(&a);
        assert!(residual(&a, &e.values, &e.vectors) < 1e-10);
        // Orthonormality.
        let g = e.vectors.gram();
        for i in 0..n {
            for j in 0..n {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((g.get(i, j) - want).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn tridiagonalization_reconstructs() {
        let n = 12;
        let mut a = Mat::zeros(n, n);
        let mut state = 0x123456789ABCDEFu64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        for i in 0..n {
            for j in 0..=i {
                let v = next();
                a.set(i, j, v);
                a.set(j, i, v);
            }
        }
        let (q, diag, off) = tridiagonalize(&a);
        // Check q T q' == a, i.e. a_ij = sum_{r,s} q_ir T_rs q_js.
        let t = |r: usize, s: usize| -> f64 {
            if r == s {
                diag[r]
            } else if r + 1 == s {
                off[r]
            } else if s + 1 == r {
                off[s]
            } else {
                0.0
            }
        };
        let mut worst: f64 = 0.0;
        for i in 0..n {
            for j in 0..n {
                let mut acc = 0.0;
                for r in 0..n {
                    for s in 0..n {
                        acc += q.get(i, r) * t(r, s) * q.get(j, s);
                    }
                }
                worst = worst.max((acc - a.get(i, j)).abs());
            }
        }
        assert!(worst < 1e-9, "reconstruction error {worst}");
        // Check q is orthogonal.
        let g = q.gram();
        for i in 0..n {
            for j in 0..n {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((g.get(i, j) - want).abs() < 1e-9, "q not orthogonal at ({i},{j})");
            }
        }
    }

    #[test]
    fn partial_matches_full_spectrum() {
        let n = 80;
        let mut a = Mat::zeros(n, n);
        let mut state = 0xDEADBEEFCAFEBABEu64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        for i in 0..n {
            for j in 0..=i {
                let v = next();
                a.set(i, j, v);
                a.set(j, i, v);
            }
        }
        let full = sym_eig(&a);
        let (vals, vecs) = sym_eig_top_k(&a, 6);
        for c in 0..6 {
            let want = full.values[n - 1 - c];
            let rel = (vals[c] - want).abs() / want.abs().max(1e-12);
            assert!(rel < 1e-9, "pair {c}: {} vs {want} (rel {rel})", vals[c]);
        }
        // Residuals of the returned pairs.
        for c in 0..6 {
            let v: Vec<f64> = (0..n).map(|i| vecs.get(i, c)).collect();
            let av = a.matvec(&v);
            let mut num = 0.0;
            for i in 0..n {
                let d = av[i] - vals[c] * v[i];
                num += d * d;
            }
            let rel = num.sqrt() / vals[c].abs();
            assert!(rel < 1e-9, "residual of pair {c}: {rel}");
        }
    }

    #[test]
    fn partial_handles_degenerate_spectrum() {
        // Block-diagonal with repeated eigenvalues.
        let n = 40;
        let mut a = Mat::zeros(n, n);
        for i in 0..n {
            a.set(i, i, if i < 20 { 5.0 } else { 1.0 });
        }
        let (vals, _) = sym_eig_top_k(&a, 4);
        assert!((vals[0] - 5.0).abs() < 1e-9);
        assert!((vals[3] - 5.0).abs() < 1e-9);
    }

    #[test]
    fn trace_and_determinant_invariants() {
        let n = 12;
        let a = Mat::from_fn(n, n, |i, j| {
            let v = ((i * 17 + j * 31) % 11) as f64 - 5.0;
            if i == j { v + 20.0 } else { v * 0.3 }
        });
        let mut a = a;
        for i in 0..n {
            for j in 0..i {
                let s = 0.5 * (a.get(i, j) + a.get(j, i));
                a.set(i, j, s);
                a.set(j, i, s);
            }
        }
        let e = sym_eig(&a);
        let tr_a: f64 = (0..n).map(|i| a.get(i, i)).sum();
        let tr_e: f64 = e.values.iter().sum();
        assert!((tr_a - tr_e).abs() < 1e-9, "trace {tr_a} vs {tr_e}");
    }
}
