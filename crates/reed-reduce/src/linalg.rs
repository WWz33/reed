//! Minimal dense linear algebra on `Vec<f64>` with no external crates.
//!
//! Everything here is deliberately small and explicit: the point is to have a
//! reference implementation that can be audited against LAPACK/faer results,
//! not to replace them.

/// A dense row-major matrix.
#[derive(Clone, Debug, PartialEq)]
pub struct Mat {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<f64>,
}

impl Mat {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self { rows, cols, data: vec![0.0; rows * cols] }
    }

    pub fn from_fn(rows: usize, cols: usize, mut f: impl FnMut(usize, usize) -> f64) -> Self {
        let mut m = Self::zeros(rows, cols);
        for i in 0..rows {
            for j in 0..cols {
                m.data[i * cols + j] = f(i, j);
            }
        }
        m
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        self.data[i * self.cols + j]
    }

    #[inline]
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        self.data[i * self.cols + j] = v;
    }

    pub fn row(&self, i: usize) -> &[f64] {
        &self.data[i * self.cols..(i + 1) * self.cols]
    }

    pub fn row_mut(&mut self, i: usize) -> &mut [f64] {
        &mut self.data[i * self.cols..(i + 1) * self.cols]
    }

    /// `self * other` (row-major, blocked, single-threaded).
    pub fn matmul(&self, other: &Mat) -> Mat {
        assert_eq!(self.cols, other.rows, "matmul shape mismatch");
        let mut out = Mat::zeros(self.rows, other.cols);
        const BLOCK: usize = 32;
        for ii in (0..self.rows).step_by(BLOCK) {
            let i_end = (ii + BLOCK).min(self.rows);
            for kk in (0..self.cols).step_by(BLOCK) {
                let k_end = (kk + BLOCK).min(self.cols);
                for jj in (0..other.cols).step_by(BLOCK) {
                    let j_end = (jj + BLOCK).min(other.cols);
                    for i in ii..i_end {
                        for k in kk..k_end {
                            let a = self.data[i * self.cols + k];
                            if a == 0.0 {
                                continue;
                            }
                            let dst = &mut out.data[i * other.cols + jj..i * other.cols + j_end];
                            let src = &other.data[k * other.cols + jj..k * other.cols + j_end];
                            for (d, s) in dst.iter_mut().zip(src) {
                                *d += a * s;
                            }
                        }
                    }
                }
            }
        }
        out
    }

    pub fn transpose(&self) -> Mat {
        Mat::from_fn(self.cols, self.rows, |i, j| self.get(j, i))
    }

    /// `self * v` where `v` is a single column.
    pub fn matvec(&self, v: &[f64]) -> Vec<f64> {
        assert_eq!(self.cols, v.len(), "matvec shape mismatch");
        let mut out = vec![0.0; self.rows];
        for (i, o) in out.iter_mut().enumerate() {
            let row = self.row(i);
            let mut acc = 0.0;
            for (a, b) in row.iter().zip(v) {
                acc += a * b;
            }
            *o = acc;
        }
        out
    }

    /// `self' * v`.
    pub fn tmatvec(&self, v: &[f64]) -> Vec<f64> {
        assert_eq!(self.rows, v.len(), "tmatvec shape mismatch");
        let mut out = vec![0.0; self.cols];
        for (i, vi) in v.iter().enumerate() {
            if *vi == 0.0 {
                continue;
            }
            let row = self.row(i);
            for (o, a) in out.iter_mut().zip(row) {
                *o += a * vi;
            }
        }
        out
    }

    /// `self' * self`.
    pub fn gram(&self) -> Mat {
        let mut out = Mat::zeros(self.cols, self.cols);
        for i in 0..self.cols {
            for j in 0..=i {
                let mut acc = 0.0;
                for r in 0..self.rows {
                    acc += self.get(r, i) * self.get(r, j);
                }
                out.set(i, j, acc);
                out.set(j, i, acc);
            }
        }
        out
    }
}

pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

pub fn norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

/// Squared Euclidean distance.
pub fn dist2(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y) * (x - y)).sum()
}

/// In-place normalization to unit length; returns the original norm.
pub fn normalize(v: &mut [f64]) -> f64 {
    let n = norm(v);
    if n > 0.0 {
        for x in v.iter_mut() {
            *x /= n;
        }
    }
    n
}

/// Stable sort of indices by a key, ascending.
pub fn argsort_by(v: &[f64]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap_or(std::cmp::Ordering::Equal));
    idx
}

/// Stable sort of indices by a key, descending.
pub fn argsort_desc(v: &[f64]) -> Vec<usize> {
    let mut idx = argsort_by(v);
    idx.reverse();
    idx
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matmul_matches_naive() {
        let a = Mat::from_fn(7, 5, |i, j| (i * 3 + j) as f64 * 0.5);
        let b = Mat::from_fn(5, 9, |i, j| (i + j * 2) as f64 * 0.25);
        let c = a.matmul(&b);
        for i in 0..7 {
            for j in 0..9 {
                let mut want = 0.0;
                for k in 0..5 {
                    want += a.get(i, k) * b.get(k, j);
                }
                assert!((c.get(i, j) - want).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn gram_is_symmetric() {
        let a = Mat::from_fn(11, 4, |i, j| ((i * 7 + j * 3) % 5) as f64);
        let g = a.gram();
        for i in 0..4 {
            for j in 0..4 {
                assert_eq!(g.get(i, j), g.get(j, i));
            }
        }
    }
}
