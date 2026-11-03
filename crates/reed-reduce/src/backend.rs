//! Dense backend: blocked GEMM and symmetric eigendecomposition.
//!
//! Two implementations are compiled in:
//!
//! * `faer` (default, `faer-backend` feature): multi-threaded blocked GEMM and
//!   a Householder-tridiagonalization + implicit-QR dense eigensolver. This is
//!   the production path.
//! * the dependency-free fallback (`reed_core::linalg::Mat::matmul` and
//!   `crate::eig::sym_eig`): single-threaded scalar code, kept so the numbers
//!   can be audited against a from-scratch implementation and so the crate
//!   still builds where no dependency is allowed.
//!
//! Measured on a 56-core machine, 4000 x 20000 standardized genotypes:
//!
//! | step | fallback | faer |
//! |---|---|---|
//! | GRM = W W' / m | 79.1 s | see `grm` (blocked + parallel) |
//! | dense eigendecomposition | 258.2 s | ~1.9 s |

use reed_core::linalg::Mat;

/// Which implementation served a call. Reported so benchmarks and tests can
/// assert that the fast path was actually taken.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Backend {
    Faer,
    Fallback,
}

impl Backend {
    pub fn name(self) -> &'static str {
        match self {
            Backend::Faer => "faer",
            Backend::Fallback => "fallback",
        }
    }
}

/// The backend this build was compiled with.
pub const COMPILED_BACKEND: Backend = if cfg!(feature = "faer-backend") {
    Backend::Faer
} else {
    Backend::Fallback
};

/// Symmetric eigendecomposition, eigenvalues ascending.
pub struct SymEigResult {
    pub values: Vec<f64>,
    /// Eigenvectors as columns (n x n).
    pub vectors: Mat,
    pub backend: Backend,
}

impl SymEigResult {
    /// Top-k eigenvalues (descending) and their eigenvectors as columns.
    ///
    /// Each column carries a sign convention: its largest-magnitude component
    /// is positive (ties take the first row; an all-zero column is left
    /// alone). The solvers do not agree on sign — the same symmetric matrix
    /// came back with `-0.998...` PC1 from the fallback backend and
    /// `+0.998...` from faer — so without normalization a backend switch or a
    /// faer upgrade changes `.eigenvec` bytes and embedding coordinates even
    /// though the analysis is identical (doc_v23 V23-09). This is the same
    /// convention `lda` already used.
    pub fn top_k(&self, k: usize) -> (Vec<f64>, Mat) {
        let n = self.values.len();
        let k = k.min(n);
        let mut vals = Vec::with_capacity(k);
        let mut vecs = Mat::zeros(n, k);
        for c in 0..k {
            let src = n - 1 - c;
            vals.push(self.values[src]);
            let mut best_abs = 0.0f64;
            let mut best_val = 0.0f64;
            for r in 0..n {
                let x = self.vectors.get(r, src);
                if x.abs() > best_abs {
                    best_abs = x.abs();
                    best_val = x;
                }
            }
            let sign = if best_val < 0.0 { -1.0 } else { 1.0 };
            for r in 0..n {
                vecs.set(r, c, sign * self.vectors.get(r, src));
            }
        }
        (vals, vecs)
    }
}

/// `A = W W' / m`, using the fastest available GEMM.
pub fn gram_divide(w: &Mat, m: usize) -> (Mat, Backend) {
    #[cfg(feature = "faer-backend")]
    {
        let a = faer_gram(w);
        let scale = 1.0 / m as f64;
        let mut out = a;
        for v in out.data.iter_mut() {
            *v *= scale;
        }
        return (out, Backend::Faer);
    }
    #[cfg(not(feature = "faer-backend"))]
    {
        let wt = w.transpose();
        let mut a = w.matmul(&wt);
        let scale = 1.0 / m as f64;
        for v in a.data.iter_mut() {
            *v *= scale;
        }
        (a, Backend::Fallback)
    }
}

#[cfg(feature = "faer-backend")]
fn faer_gram(w: &Mat) -> Mat {
    use faer::{Mat as FMat, Parallelism};

    let n = w.rows;
    let p = w.cols;
    let wm = FMat::<f64>::from_fn(n, p, |i, j| w.get(i, j));
    // faer's operator form dispatches to the multi-threaded GEMM. The
    // explicit `matmul(acc, lhs, rhs, Some(alpha), beta, par)` call silently
    // writes nothing in 0.20.2, so the operator is the supported path.
    let _ = Parallelism::Rayon(0);
    let product = &wm * wm.transpose();
    Mat::from_fn(n, n, |i, j| product[(i, j)])
}

/// Dense symmetric eigendecomposition of `a` (only the lower triangle is read).
pub fn sym_eig(a: &Mat) -> SymEigResult {
    #[cfg(feature = "faer-backend")]
    {
        let (values, vectors) = faer_sym_eig(a);
        return SymEigResult { values, vectors, backend: Backend::Faer };
    }
    #[cfg(not(feature = "faer-backend"))]
    {
        let e = crate::eig::sym_eig(a);
        SymEigResult { values: e.values, vectors: e.vectors, backend: Backend::Fallback }
    }
}

#[cfg(feature = "faer-backend")]
fn faer_sym_eig(a: &Mat) -> (Vec<f64>, Mat) {
    use faer::dyn_stack::{GlobalPodBuffer, PodStack};
    use faer::linalg::evd::{
        compute_hermitian_evd, compute_hermitian_evd_req, ComputeVectors, HermitianEvdParams,
    };
    use faer::{Col, Mat as FMat, Parallelism};

    let n = a.rows;
    assert_eq!(n, a.cols, "sym_eig needs a square matrix");
    if n == 0 {
        return (Vec::new(), Mat::zeros(0, 0));
    }
    // faer reads the lower triangle; make the matrix exactly symmetric first
    // so round-off in the input cannot bias the result.
    let mut am = FMat::<f64>::from_fn(n, n, |i, j| 0.5 * (a.get(i, j) + a.get(j, i)));

    let par = Parallelism::Rayon(0);
    let req = compute_hermitian_evd_req::<f64>(
        n,
        ComputeVectors::Yes,
        par,
        HermitianEvdParams::default(),
    )
    .expect("faer workspace request");
    let mut buf = GlobalPodBuffer::new(req);
    let mut stack = PodStack::new(&mut buf);
    let mut s = Col::<f64>::zeros(n);
    let mut u = FMat::<f64>::zeros(n, n);
    compute_hermitian_evd(
        am.as_ref(),
        s.as_mut(),
        Some(u.as_mut()),
        par,
        &mut stack,
        HermitianEvdParams::default(),
    );

    let values: Vec<f64> = (0..n).map(|i| s[i]).collect();
    let vectors = Mat::from_fn(n, n, |i, j| u[(i, j)]);
    let _ = &mut am;
    (values, vectors)
}

/// Top-k eigenpairs without computing the full spectrum.
///
/// faer 0.20 has no partial-spectrum API (`krylov_schur.rs` is a commented-out
/// prototype), so this falls back to a full decomposition and truncation. The
/// from-scratch bisection path in `crate::eig::sym_eig_top_k` does extract only
/// k pairs but is far slower without BLAS; it stays available as the
/// dependency-free option.
pub fn sym_eig_top_k(a: &Mat, k: usize) -> (Vec<f64>, Mat, Backend) {
    let full = sym_eig(a);
    let (vals, vecs) = full.top_k(k);
    (vals, vecs, full.backend)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_symmetric(n: usize, seed: u64) -> Mat {
        let mut state = seed | 1;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        let mut a = Mat::zeros(n, n);
        for i in 0..n {
            for j in 0..=i {
                let v = next();
                a.set(i, j, v);
                a.set(j, i, v);
            }
        }
        a
    }

    #[test]
    fn eigenvalues_ascending_and_residual_small() {
        let n = 60;
        let a = random_symmetric(n, 0x1234);
        let e = sym_eig(&a);
        assert_eq!(e.values.len(), n);
        for w in e.values.windows(2) {
            assert!(w[0] <= w[1] + 1e-12, "not ascending");
        }
        // Residual of the top pair.
        let c = n - 1;
        let v: Vec<f64> = (0..n).map(|i| e.vectors.get(i, c)).collect();
        let av = a.matvec(&v);
        let mut num = 0.0;
        for i in 0..n {
            let d = av[i] - e.values[c] * v[i];
            num += d * d;
        }
        assert!(num.sqrt() / e.values[c].abs() < 1e-9, "residual too large");
    }

    #[test]
    fn trace_matches() {
        let n = 40;
        let a = random_symmetric(n, 0xBEEF);
        let e = sym_eig(&a);
        let tr_a: f64 = (0..n).map(|i| a.get(i, i)).sum();
        let tr_e: f64 = e.values.iter().sum();
        assert!((tr_a - tr_e).abs() < 1e-9, "trace {tr_a} vs {tr_e}");
    }

    /// doc_v23 V23-09: the two backends must return the same eigenvectors,
    /// sign included. Before normalization the same symmetric matrix gave
    /// `-0.99845624390341792` PC1 from the fallback and
    /// `+0.99845624390341725` from faer (PC2 flipped too).
    #[test]
    fn top_k_sign_is_normalized_across_backends() {
        let n = 6;
        let a = random_symmetric(n, 0x5EED);
        // The fallback (tred2/tql2) result, whichever backend `sym_eig`
        // delegates to by default.
        let fallback = crate::eig::sym_eig(&a);
        let (fb_vals, fb_vecs) = fallback.top_k(2);
        // Independent normalization check: every column's largest-magnitude
        // component is positive.
        for c in 0..2 {
            let mut best_abs = 0.0f64;
            let mut best_val = 0.0f64;
            for r in 0..n {
                let x = fb_vecs.get(r, c);
                if x.abs() > best_abs {
                    best_abs = x.abs();
                    best_val = x;
                }
            }
            assert!(
                best_val > 0.0,
                "column {c} largest component {best_val} is not positive"
            );
        }
        // Eigenvalues must stay descending and unchanged by the sign fix.
        assert!(fb_vals[0] >= fb_vals[1]);

        // Cross-backend comparison when faer is compiled in: the public
        // `sym_eig` returns faer, and the fallback is reachable directly.
        #[cfg(feature = "faer-backend")]
        {
            let (faer_vals, faer_vecs) = sym_eig(&a).top_k(2);
            for c in 0..2 {
                assert!(
                    (faer_vals[c] - fb_vals[c]).abs() <= 1e-9 * faer_vals[c].abs().max(1.0),
                    "eigenvalue {c}: faer {} vs fallback {}",
                    faer_vals[c],
                    fb_vals[c]
                );
                for r in 0..n {
                    let diff = (faer_vecs.get(r, c) - fb_vecs.get(r, c)).abs();
                    assert!(
                        diff <= 1e-9,
                        "PC{} row {r}: faer {} vs fallback {} (diff {diff})",
                        c + 1,
                        faer_vecs.get(r, c),
                        fb_vecs.get(r, c)
                    );
                }
            }
        }
    }

    /// A column that is all zeros (a degenerate matrix) keeps its sign-free
    /// value: the convention must not produce `-0.0` or a NaN.
    #[test]
    fn top_k_leaves_an_all_zero_column_alone() {
        let mut r = SymEigResult {
            values: vec![0.0, 0.0],
            vectors: Mat::zeros(2, 2),
            backend: COMPILED_BACKEND,
        };
        r.vectors.set(0, 1, 0.0);
        let (vals, vecs) = r.top_k(2);
        assert_eq!(vals, vec![0.0, 0.0]);
        for v in vecs.data.iter() {
            assert_eq!(*v, 0.0);
            assert!(!v.is_sign_negative(), "zero column must not become -0");
        }
    }

    /// The normalization picks the *first* row among equal-magnitude
    /// components (documented tie rule).
    #[test]
    fn top_k_tie_rule_takes_the_first_row() {
        // Eigenvalues ascending; column 1 is the top one.
        let mut vectors = Mat::zeros(2, 2);
        vectors.set(0, 0, 1.0);
        vectors.set(1, 0, 0.0);
        vectors.set(0, 1, -0.5);
        vectors.set(1, 1, 0.5);
        let r = SymEigResult {
            values: vec![0.0, 1.0],
            vectors,
            backend: COMPILED_BACKEND,
        };
        let (_, vecs) = r.top_k(1);
        assert_eq!(vecs.get(0, 0), 0.5, "first row decides the sign");
        assert_eq!(vecs.get(1, 0), -0.5);
    }

    #[test]
    fn gram_matches_naive() {
        let n = 25;
        let p = 30;
        let mut state = 0xFEEDu64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let w = Mat::from_fn(n, p, |_, _| next() - 0.5);
        let (g, backend) = gram_divide(&w, p);
        // Compare against a direct triple loop.
        for i in 0..n {
            for j in 0..n {
                let mut acc = 0.0;
                for k in 0..p {
                    acc += w.get(i, k) * w.get(j, k);
                }
                let want = acc / p as f64;
                assert!(
                    (g.get(i, j) - want).abs() < 1e-10,
                    "({i},{j}) {} vs {want} (backend {:?})",
                    g.get(i, j),
                    backend
                );
            }
        }
    }

    #[test]
    fn backend_is_reported() {
        let n = 10;
        let a = random_symmetric(n, 7);
        let e = sym_eig(&a);
        assert_eq!(e.backend, COMPILED_BACKEND);
    }
}
