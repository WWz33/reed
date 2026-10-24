//! Fisher linear discriminant analysis.
//!
//! The algorithm is `MASS:::lda.default`'s `method = "moment"` branch,
//! which is what `adegenet::dapc` calls (with `tol = 1e-30` so that no
//! singular value is discarded for being "small"). It is written here from
//! the published discriminant-analysis formulation (Fisher 1936; Rao 1948)
//! rather than translated: the steps below are the two-stage whitening that
//! `lda` performs, spelled out.
//!
//! Stage 1 divides every column by its **within-group** residual standard
//! deviation and whitens the within-group scatter. Stage 2 projects the
//! class means onto that whitened space and takes the leading directions.
//! The result is the same Fisher discriminant basis, in the same units.
//!
//! # Sign convention
//!
//! Singular vectors are defined up to a sign. This implementation fixes the
//! sign of every returned direction so that its largest-magnitude component
//! is positive, which makes the output reproducible. `MASS`/LAPACK leave the
//! sign to the SVD routine, so a parity check must compare `scaling` columns
//! by `|dot|`, not entry-wise. Posterior probabilities and assignments are
//! invariant to the sign.

use reed_core::linalg::Mat;

/// Errors from fitting or applying an LDA.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LdaError {
    /// The feature matrix had no rows or no columns.
    NoObservations,
    /// `groups` length did not match the number of rows.
    GroupCountMismatch {
        /// Rows in the feature matrix.
        rows: usize,
        /// Entries in `groups`.
        groups: usize,
    },
    /// A group index was `>= ng`.
    GroupIndexOutOfRange {
        /// The offending index.
        index: u32,
        /// Number of groups.
        ng: usize,
    },
    /// A group had no observations.
    EmptyGroup {
        /// Index of the empty group.
        group: usize,
    },
    /// Fewer than two groups were supplied.
    TooFewGroups(usize),
    /// A feature was constant within every group.
    ConstantWithinGroups {
        /// Column index of the offending feature.
        column: usize,
    },
    /// The between-group scatter had no positive singular value.
    Degenerate,
    /// `dimen` was zero.
    ZeroDimen,
    /// `project` received a different feature count than `fit` saw.
    FeatureCountMismatch {
        /// Features the model was fitted on.
        expected: usize,
        /// Features supplied to `project`.
        got: usize,
    },
}

impl std::fmt::Display for LdaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NoObservations => write!(f, "the feature matrix is empty"),
            Self::GroupCountMismatch { rows, groups } => {
                write!(f, "group vector has {groups} entries for {rows} rows")
            }
            Self::GroupIndexOutOfRange { index, ng } => {
                write!(f, "group index {index} is out of range for {ng} groups")
            }
            Self::EmptyGroup { group } => {
                write!(f, "group {group} has no observations")
            }
            Self::TooFewGroups(ng) => {
                write!(f, "discriminant analysis needs at least 2 groups, got {ng}")
            }
            Self::ConstantWithinGroups { column } => {
                write!(f, "feature {column} is constant within every group")
            }
            Self::Degenerate => write!(f, "group means are numerically identical"),
            Self::ZeroDimen => write!(f, "the number of discriminant axes must be at least 1"),
            Self::FeatureCountMismatch { expected, got } => {
                write!(f, "the model expects {expected} features, got {got}")
            }
        }
    }
}

impl std::error::Error for LdaError {}

/// A fitted linear discriminant model.
#[derive(Clone, Debug)]
pub struct Lda {
    /// Class priors (proportions, summing to 1).
    pub prior: Vec<f64>,
    /// Class counts.
    pub counts: Vec<usize>,
    /// Class means in the input feature space (`ng` rows, `p` columns).
    pub means: Mat,
    /// Class labels, in the order every other field uses.
    pub levels: Vec<String>,
    /// Discriminant directions in the input feature space (`p` rows,
    /// `rank` columns).
    pub scaling: Mat,
    /// Singular values of the between-class scatter, descending.
    pub svd: Vec<f64>,
    /// Number of observations the model was fitted on.
    pub n: usize,
    /// Group index per observation, in the order `fit` received them.
    pub groups: Vec<u32>,
}

/// A fitted model applied to data.
#[derive(Clone, Debug)]
pub struct LdaProjection {
    /// Discriminant coordinates (`n` rows, `dimen` columns).
    pub coords: Mat,
    /// Posterior class probabilities (`n` rows, `ng` columns).
    pub posterior: Mat,
    /// Assigned class index per observation.
    pub assign: Vec<u32>,
}

/// Fits an LDA on `x` (`n` rows, `p` features) with the given group labels.
///
/// `levels` names the groups; `groups[i]` indexes into it.
///
/// # Errors
///
/// See [`LdaError`]. A feature that is constant within every group is an
/// error (as in `MASS::lda`), not a silently dropped column.
pub fn fit(x: &Mat, groups: &[u32], levels: &[String], tol: f64) -> Result<Lda, LdaError> {
    let n = x.rows;
    let p = x.cols;
    let ng = levels.len();
    if n == 0 || p == 0 {
        return Err(LdaError::NoObservations);
    }
    if groups.len() != n {
        return Err(LdaError::GroupCountMismatch {
            rows: n,
            groups: groups.len(),
        });
    }
    if ng < 2 {
        return Err(LdaError::TooFewGroups(ng));
    }
    for &g in groups {
        if g as usize >= ng {
            return Err(LdaError::GroupIndexOutOfRange { index: g, ng });
        }
    }

    // Class counts, priors and means.
    let mut counts = vec![0usize; ng];
    for &g in groups {
        counts[g as usize] += 1;
    }
    if let Some(empty) = counts.iter().position(|&c| c == 0) {
        return Err(LdaError::EmptyGroup { group: empty });
    }
    let prior: Vec<f64> = counts.iter().map(|&c| c as f64 / n as f64).collect();
    let mut means = Mat::zeros(ng, p);
    for (i, &g) in groups.iter().enumerate() {
        let row = means.row_mut(g as usize);
        for j in 0..p {
            row[j] += x.get(i, j);
        }
    }
    for c in 0..ng {
        let row = means.row_mut(c);
        for v in row.iter_mut() {
            *v /= counts[c] as f64;
        }
    }

    // Stage 1: divide by the within-group residual standard deviation
    // (denominator n-1, as `var()` in MASS) and whiten.
    let mut f1 = vec![0.0f64; p];
    for (i, &g) in groups.iter().enumerate() {
        for j in 0..p {
            let r = x.get(i, j) - means.get(g as usize, j);
            f1[j] += r * r;
        }
    }
    for (j, v) in f1.iter_mut().enumerate() {
        *v = (*v / (n as f64 - 1.0)).sqrt();
        if *v <= tol || v.is_nan() {
            return Err(LdaError::ConstantWithinGroups { column: j });
        }
    }

    let fac_w = 1.0 / (n as f64 - ng as f64);
    let mut xw = Mat::zeros(n, p);
    for (i, &g) in groups.iter().enumerate() {
        for j in 0..p {
            let r = x.get(i, j) - means.get(g as usize, j);
            xw.set(i, j, fac_w.sqrt() * r / f1[j]);
        }
    }
    let (d_w, v_w) = svd_right(&xw, tol);
    let rank_w = d_w.len();
    if rank_w == 0 {
        return Err(LdaError::Degenerate);
    }

    // scaling1 = diag(1/f1) * V_w * diag(1/d_w)   (p x rank_w)
    let mut scaling1 = Mat::zeros(p, rank_w);
    for j in 0..p {
        for c in 0..rank_w {
            scaling1.set(j, c, v_w.get(j, c) / (f1[j] * d_w[c]));
        }
    }

    // Stage 2: project the class means onto the whitened space and take the
    // leading directions.
    let mut xbar = vec![0.0f64; p];
    for c in 0..ng {
        for j in 0..p {
            xbar[j] += prior[c] * means.get(c, j);
        }
    }
    let fac_b = 1.0 / (ng as f64 - 1.0);
    let mut xb = Mat::zeros(ng, rank_w);
    for c in 0..ng {
        let w = (n as f64 * prior[c] * fac_b).sqrt();
        for k in 0..rank_w {
            let mut acc = 0.0;
            for j in 0..p {
                acc += (means.get(c, j) - xbar[j]) * scaling1.get(j, k);
            }
            xb.set(c, k, w * acc);
        }
    }
    let (d_b, v_b) = svd_right(&xb, tol * d_b_first(&xb));
    let rank_b = d_b.len();
    if rank_b == 0 {
        return Err(LdaError::Degenerate);
    }

    // scaling = scaling1 * V_b[:, :rank_b]   (p x rank_b)
    let mut scaling = Mat::zeros(p, rank_b);
    for j in 0..p {
        for k in 0..rank_b {
            let mut acc = 0.0;
            for c in 0..rank_w {
                acc += scaling1.get(j, c) * v_b.get(c, k);
            }
            scaling.set(j, k, acc);
        }
    }

    Ok(Lda {
        prior,
        counts,
        means,
        levels: levels.to_vec(),
        scaling,
        svd: d_b,
        n,
        groups: groups.to_vec(),
    })
}

impl Lda {
    /// Number of discriminant axes available.
    #[must_use]
    pub fn rank(&self) -> usize {
        self.svd.len()
    }

    /// Projects `x` onto the first `dimen` discriminant axes.
    ///
    /// # Errors
    ///
    /// `ZeroDimen` for `dimen == 0`; `FeatureCountMismatch` when `x` has the
    /// wrong number of columns.
    pub fn project(&self, x: &Mat, dimen: usize) -> Result<LdaProjection, LdaError> {
        let p = self.means.cols;
        if x.cols != p {
            return Err(LdaError::FeatureCountMismatch {
                expected: p,
                got: x.cols,
            });
        }
        if dimen == 0 {
            return Err(LdaError::ZeroDimen);
        }
        let dimen = dimen.min(self.rank());
        if dimen == 0 {
            return Err(LdaError::ZeroDimen);
        }
        let ng = self.prior.len();
        let n = x.rows;

        // Centring point: the prior-weighted mean of the class means.
        let mut centre = vec![0.0f64; p];
        for c in 0..ng {
            for j in 0..p {
                centre[j] += self.prior[c] * self.means.get(c, j);
            }
        }

        // xc = (x - centre) * scaling, truncated to `dimen`.
        let mut xc = Mat::zeros(n, dimen);
        for i in 0..n {
            for k in 0..dimen {
                let mut acc = 0.0;
                for j in 0..p {
                    acc += (x.get(i, j) - centre[j]) * self.scaling.get(j, k);
                }
                xc.set(i, k, acc);
            }
        }

        // dm = (means - centre) * scaling, truncated to `dimen`.
        let mut dm = Mat::zeros(ng, dimen);
        for c in 0..ng {
            for k in 0..dimen {
                let mut acc = 0.0;
                for j in 0..p {
                    acc += (self.means.get(c, j) - centre[j]) * self.scaling.get(j, k);
                }
                dm.set(c, k, acc);
            }
        }

        // dist = 0.5*|dm_c|^2 - log(prior_c) - xc . dm_c
        //
        // This is `0.5*|xc - dm_c|^2 - 0.5*|xc|^2 - log(prior_c)`: a squared
        // distance to the class centroid minus a term shared by every class,
        // so the *smallest* `dist` is the best class and `posterior` is
        // proportional to `exp(-dist)`.
        let mut posterior = Mat::zeros(n, ng);
        let mut assign = Vec::with_capacity(n);
        let mut dist = vec![0.0f64; ng];
        for i in 0..n {
            let mut best = 0usize;
            let mut best_v = f64::INFINITY;
            for c in 0..ng {
                let mut half_norm = 0.0;
                let mut dot = 0.0;
                for k in 0..dimen {
                    let m = dm.get(c, k);
                    half_norm += m * m;
                    dot += xc.get(i, k) * m;
                }
                let v = 0.5 * half_norm - self.prior[c].ln() - dot;
                dist[c] = v;
                if v < best_v {
                    best_v = v;
                    best = c;
                }
            }
            // exp(-(dist - min)) keeps the best class at 1.
            let mut sum = 0.0;
            for c in 0..ng {
                let e = (-(dist[c] - best_v)).exp();
                posterior.set(i, c, e);
                sum += e;
            }
            for c in 0..ng {
                posterior.set(i, c, posterior.get(i, c) / sum);
            }
            assign.push(best as u32);
        }

        Ok(LdaProjection {
            coords: xc,
            posterior,
            assign,
        })
    }
}

/// First singular value of `a`, used for the relative tolerance of the
/// second stage (`tol * d[0]` in `MASS::lda`).
fn d_b_first(a: &Mat) -> f64 {
    let mut best = 0.0f64;
    for j in 0..a.cols {
        let mut s = 0.0;
        for i in 0..a.rows {
            let v = a.get(i, j);
            s += v * v;
        }
        if s > best {
            best = s;
        }
    }
    best.sqrt()
}

/// Right singular vectors and singular values of `a`, descending, keeping
/// only those strictly greater than `tol`.
///
/// One-sided Jacobi: repeatedly rotate pairs of columns until they are
/// orthogonal. This is the numerically accurate variant (it never forms
/// `A^T A`), and for the small feature counts of an LDA it is also fast.
///
/// The sign of each returned vector is fixed so that its largest-magnitude
/// component is positive; see the module docs.
fn svd_right(a: &Mat, tol: f64) -> (Vec<f64>, Mat) {
    let n = a.rows;
    let p = a.cols;
    let mut b = a.clone();
    let mut v = Mat::zeros(p, p);
    for i in 0..p {
        v.set(i, i, 1.0);
    }

    // A sweep count cap keeps a pathological input from looping forever;
    // 60 sweeps is far past convergence for these sizes.
    for _sweep in 0..60 {
        let mut rotated = false;
        for i in 0..p {
            for j in (i + 1)..p {
                let mut alpha = 0.0;
                let mut beta = 0.0;
                let mut gamma = 0.0;
                for r in 0..n {
                    let x = b.get(r, i);
                    let y = b.get(r, j);
                    alpha += x * x;
                    beta += y * y;
                    gamma += x * y;
                }
                if gamma == 0.0 {
                    continue;
                }
                let scale = (alpha * beta).sqrt();
                if gamma.abs() <= 1e-15 * scale {
                    continue;
                }
                let zeta = (beta - alpha) / (2.0 * gamma);
                let t = if zeta >= 0.0 {
                    1.0 / (zeta + (1.0 + zeta * zeta).sqrt())
                } else {
                    -1.0 / (-zeta + (1.0 + zeta * zeta).sqrt())
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = c * t;
                for r in 0..n {
                    let x = b.get(r, i);
                    let y = b.get(r, j);
                    b.set(r, i, c * x - s * y);
                    b.set(r, j, s * x + c * y);
                }
                for r in 0..p {
                    let x = v.get(r, i);
                    let y = v.get(r, j);
                    v.set(r, i, c * x - s * y);
                    v.set(r, j, s * x + c * y);
                }
                rotated = true;
            }
        }
        if !rotated {
            break;
        }
    }

    let mut pairs: Vec<(f64, usize)> = (0..p)
        .map(|j| {
            let mut s = 0.0;
            for r in 0..n {
                let x = b.get(r, j);
                s += x * x;
            }
            (s.sqrt(), j)
        })
        .filter(|(d, _)| *d > tol)
        .collect();
    pairs.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap_or(std::cmp::Ordering::Equal));

    let rank = pairs.len();
    let mut d = Vec::with_capacity(rank);
    let mut out = Mat::zeros(p, rank);
    for (col, &(val, src)) in pairs.iter().enumerate() {
        d.push(val);
        // Sign convention: the largest-magnitude component is positive.
        let mut best_abs = 0.0f64;
        let mut best_val = 0.0f64;
        for r in 0..p {
            let x = v.get(r, src);
            if x.abs() > best_abs {
                best_abs = x.abs();
                best_val = x;
            }
        }
        let sign = if best_val < 0.0 { -1.0 } else { 1.0 };
        for r in 0..p {
            out.set(r, col, sign * v.get(r, src));
        }
    }
    (d, out)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn levels(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("g{i}")).collect()
    }

    /// Two groups, two features: the discriminating direction is x.
    fn two_group_fixture() -> (Mat, Vec<u32>) {
        let mut x = Mat::zeros(8, 2);
        for i in 0..4 {
            x.set(i, 0, -1.0 - i as f64 * 0.1);
            x.set(i, 1, (i as f64) * 0.05);
        }
        for i in 4..8 {
            x.set(i, 0, 1.0 + (i - 4) as f64 * 0.1);
            x.set(i, 1, (i as f64) * 0.05);
        }
        (x, vec![0, 0, 0, 0, 1, 1, 1, 1])
    }

    #[test]
    fn priors_and_means_match_the_groups() {
        let (x, g) = two_group_fixture();
        let m = fit(&x, &g, &levels(2), 1e-30).unwrap();
        assert_eq!(m.counts, vec![4, 4]);
        assert!((m.prior[0] - 0.5).abs() < 1e-15);
        assert!((m.means.get(0, 0) + 1.15).abs() < 1e-12);
        assert!((m.means.get(1, 0) - 1.15).abs() < 1e-12);
    }

    #[test]
    fn two_groups_yield_one_discriminant_axis() {
        let (x, g) = two_group_fixture();
        let m = fit(&x, &g, &levels(2), 1e-30).unwrap();
        assert_eq!(m.rank(), 1);
        assert_eq!(m.scaling.cols, 1);
        // The first feature separates the groups; the second does not.
        assert!(m.scaling.get(0, 0).abs() > m.scaling.get(1, 0).abs());
    }

    #[test]
    fn projections_are_perfectly_separated() {
        let (x, g) = two_group_fixture();
        let m = fit(&x, &g, &levels(2), 1e-30).unwrap();
        let p = m.project(&x, 1).unwrap();
        for i in 0..4 {
            assert_eq!(p.assign[i], 0);
            assert!(p.posterior.get(i, 0) > 0.5);
        }
        for i in 4..8 {
            assert_eq!(p.assign[i], 1);
            assert!(p.posterior.get(i, 1) > 0.5);
        }
    }

    #[test]
    fn posteriors_sum_to_one() {
        let (x, g) = two_group_fixture();
        let m = fit(&x, &g, &levels(2), 1e-30).unwrap();
        let p = m.project(&x, 1).unwrap();
        for i in 0..x.rows {
            let s: f64 = (0..2).map(|c| p.posterior.get(i, c)).sum();
            assert!((s - 1.0).abs() < 1e-12, "row {i} sums to {s}");
        }
    }

    #[test]
    fn constant_within_groups_is_rejected() {
        let x = Mat::from_fn(6, 2, |i, j| if j == 1 { 7.0 } else { i as f64 });
        let g = vec![0, 0, 0, 1, 1, 1];
        let err = fit(&x, &g, &levels(2), 1e-30).unwrap_err();
        assert_eq!(err, LdaError::ConstantWithinGroups { column: 1 });
    }

    #[test]
    fn empty_group_is_rejected() {
        let x = Mat::from_fn(4, 1, |i, _| i as f64);
        let g = vec![0, 0, 0, 0];
        let err = fit(&x, &g, &levels(2), 1e-30).unwrap_err();
        assert_eq!(err, LdaError::EmptyGroup { group: 1 });
    }

    #[test]
    fn single_group_is_rejected() {
        let x = Mat::from_fn(4, 1, |i, _| i as f64);
        let g = vec![0, 0, 0, 0];
        let err = fit(&x, &g, &levels(1), 1e-30).unwrap_err();
        assert_eq!(err, LdaError::TooFewGroups(1));
    }

    #[test]
    fn group_index_out_of_range_is_rejected() {
        let x = Mat::from_fn(4, 1, |i, _| i as f64);
        let g = vec![0, 0, 0, 2];
        let err = fit(&x, &g, &levels(2), 1e-30).unwrap_err();
        assert_eq!(err, LdaError::GroupIndexOutOfRange { index: 2, ng: 2 });
    }

    #[test]
    fn three_groups_yield_two_axes() {
        let mut x = Mat::zeros(9, 2);
        for i in 0..3 {
            x.set(i, 0, -3.0 + i as f64 * 0.1);
            x.set(i, 1, 0.0);
        }
        for i in 3..6 {
            x.set(i, 0, 3.0 + (i - 3) as f64 * 0.1);
            x.set(i, 1, 0.0);
        }
        for i in 6..9 {
            x.set(i, 0, 0.0);
            x.set(i, 1, 3.0 + (i - 6) as f64 * 0.1);
        }
        let g = vec![0, 0, 0, 1, 1, 1, 2, 2, 2];
        let m = fit(&x, &g, &levels(3), 1e-30).unwrap();
        assert_eq!(m.rank(), 2);
        let p = m.project(&x, 2).unwrap();
        for i in 0..9 {
            assert_eq!(p.assign[i], g[i] as u32, "sample {i}");
        }
    }

    #[test]
    fn dimen_is_clamped_to_rank() {
        let (x, g) = two_group_fixture();
        let m = fit(&x, &g, &levels(2), 1e-30).unwrap();
        let p = m.project(&x, 5).unwrap();
        assert_eq!(p.coords.cols, 1);
    }

    #[test]
    fn zero_dimen_is_rejected() {
        let (x, g) = two_group_fixture();
        let m = fit(&x, &g, &levels(2), 1e-30).unwrap();
        assert_eq!(m.project(&x, 0).unwrap_err(), LdaError::ZeroDimen);
    }

    #[test]
    fn feature_count_mismatch_is_rejected() {
        let (x, g) = two_group_fixture();
        let m = fit(&x, &g, &levels(2), 1e-30).unwrap();
        let bad = Mat::zeros(8, 3);
        let err = m.project(&bad, 1).unwrap_err();
        assert_eq!(
            err,
            LdaError::FeatureCountMismatch {
                expected: 2,
                got: 3
            }
        );
    }

    #[test]
    fn svd_right_recovers_known_singular_values() {
        // A diagonal matrix: singular values are the absolute diagonal.
        let a = Mat::from_fn(3, 3, |i, j| if i == j { [3.0, 2.0, 1.0][i] } else { 0.0 });
        let (d, v) = svd_right(&a, 1e-30);
        assert_eq!(d.len(), 3);
        assert!((d[0] - 3.0).abs() < 1e-12);
        assert!((d[1] - 2.0).abs() < 1e-12);
        assert!((d[2] - 1.0).abs() < 1e-12);
        assert_eq!(v.rows, 3);
    }
}
