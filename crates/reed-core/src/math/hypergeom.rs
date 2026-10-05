//! Hypergeometric projection weights for the SFS (`docs/architecture-spec.md`
//! §9.4).
//!
//! Projects the allele count of a site with `m` called alleles onto a fixed
//! sample size `n ≤ m`:
//!
//! ```text
//! P(k | n, m, alt) = C(alt, k) · C(m − alt, n − k) / C(m, n)
//! ```
//!
//! (Marth et al. 2004; dadi `project`; easySFS.) Weights are produced by the
//! ratio recurrence `P(k+1)/P(k) = (alt−k)(n−k) / ((k+1)(m−alt−n+k+1))`,
//! started at the lower support bound in log space, truncated below a
//! relative weight of 1e−15 and renormalized to sum 1.

/// Hypergeometric projection onto a fixed sample size.
///
/// `weights[i]` is the probability mass of bin `offset_k + i`; bins outside
/// the stored range have probability 0.
#[derive(Clone, Debug, PartialEq)]
pub struct ProjectionWeights {
    /// First bin `k` with non-zero weight.
    pub offset_k: u32,
    /// Non-zero weights, contiguous from `offset_k`.
    pub weights: Vec<f64>,
}

impl ProjectionWeights {
    /// Weight of bin `k` (0 outside the stored support).
    #[must_use]
    pub fn get(&self, k: u32) -> f64 {
        let idx = (k - self.offset_k) as usize; // u32 subtraction wraps; guarded below
        if k < self.offset_k || idx >= self.weights.len() {
            0.0
        } else {
            self.weights[idx]
        }
    }

    /// Sum of the stored weights (1 after renormalization, up to rounding).
    #[must_use]
    pub fn sum(&self) -> f64 {
        self.weights.iter().sum()
    }
}

/// Projects a site with `n_called` called alleles and `n_alt` ALT alleles
/// onto `n_target`.
///
/// Returns `None` when the site has fewer called alleles than the target
/// (`m < n`) — the caller skips it and counts the skip (§9.4).
#[must_use]
pub fn project(n_called: u32, n_alt: u32, n_target: u32) -> Option<ProjectionWeights> {
    let (m, n) = (u64::from(n_called), u64::from(n_target));
    let alt = u64::from(n_alt);
    if m < n || alt > m {
        return None;
    }
    // Identity cases: the projection of a fully-called site is a point mass,
    // and so are the all-REF / all-ALT extremes.
    if m == n {
        return Some(point_mass(n_alt));
    }
    if alt == 0 {
        return Some(point_mass(0));
    }
    if alt == m {
        return Some(point_mass(n_target));
    }

    let k_min = n.saturating_sub(m - alt);
    let k_max = n.min(alt);
    // Log-weight of the first supported bin.
    let mut log_w = ln_choose(alt, k_min) + ln_choose(m - alt, n - k_min) - ln_choose(m, n);
    let mut log_weights = Vec::new();
    let mut k = k_min;
    while k <= k_max {
        log_weights.push(log_w);
        // P(k+1)/P(k) = (alt−k)(n−k) / ((k+1)(m−alt−n+k+1)); the
        // denominators/numerators are strictly positive inside the support.
        // Counts are ≤ m ≤ u32::MAX, far below f64's exact-integer range.
        #[allow(clippy::cast_precision_loss)]
        let num = (alt - k) as f64 * (n - k) as f64;
        #[allow(clippy::cast_precision_loss)]
        let den = (k + 1) as f64 * (m - alt - n + k + 1) as f64;
        if num <= 0.0 || den <= 0.0 {
            break;
        }
        log_w += num.ln() - den.ln();
        k += 1;
    }

    let max_log = log_weights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mut weights: Vec<f64> = log_weights
        .iter()
        .map(|&lw| (lw - max_log).exp())
        .collect();
    // Truncate below 1e−15 of the maximum, then renormalize.
    for w in &mut weights {
        if *w < 1e-15 {
            *w = 0.0;
        }
    }
    let total: f64 = weights.iter().sum();
    for w in &mut weights {
        *w /= total;
    }
    // Trim trailing zeros so the stored support is tight.
    while weights.last() == Some(&0.0) {
        weights.pop();
    }
    Some(ProjectionWeights {
        offset_k: u32::try_from(k_min).ok()?,
        weights,
    })
}

fn point_mass(k: u32) -> ProjectionWeights {
    let mut weights = vec![0.0; k as usize + 1];
    weights[k as usize] = 1.0;
    ProjectionWeights {
        offset_k: 0,
        weights,
    }
}

/// `ln Γ(x)` via the Lanczos approximation (g = 7, n = 9 coefficients).
//
// The only cast is the loop's coefficient index (i ≤ 8); exact in f64.
#[allow(clippy::cast_precision_loss)]
fn ln_gamma(x: f64) -> f64 {
    const LANCSOS_G: f64 = 7.0;
    // Lanczos g=7, n=9 coefficients; digits trimmed to f64 round-trip.
    const COEFFS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    debug_assert!(x > 0.0, "invariant: ln_gamma called on positive argument");
    if x < 0.5 {
        // Reflection formula.
        let s = std::f64::consts::PI / (std::f64::consts::PI * x).sin();
        (std::f64::consts::PI * x).ln() + s.ln() - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let mut a = COEFFS[0];
        let t = x + LANCSOS_G + 0.5;
        for (i, &c) in COEFFS.iter().enumerate().skip(1) {
            a += c / (x + i as f64);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

/// `ln C(n, k)` for `k ≤ n`.
fn ln_choose(n: u64, k: u64) -> f64 {
    if k > n {
        return f64::NEG_INFINITY;
    }
    // Arguments are u64 counts bounded by the site's called-allele count
    // (≤ u32::MAX in every caller), far below f64's exact-integer range.
    #[allow(clippy::cast_precision_loss)]
    let (nf, kf) = (n as f64, k as f64);
    ln_gamma(nf + 1.0) - ln_gamma(kf + 1.0) - ln_gamma(nf - kf + 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_when_target_equals_called() {
        let w = project(10, 3, 10).unwrap();
        assert_eq!(w.offset_k, 0);
        assert_eq!(w.weights.iter().position(|&p| p == 1.0), Some(3));
    }

    #[test]
    fn skips_sites_below_target() {
        assert!(project(5, 2, 10).is_none());
    }

    #[test]
    fn weights_sum_to_one_and_match_direct_hypergeometric() {
        // m=20 called, alt=6, project onto n=10.
        let m = 20_u64;
        let alt = 6_u64;
        let n = 10_u64;
        let w = project(20, 6, 10).unwrap();
        assert!((w.sum() - 1.0).abs() < 1e-12);
        // Direct evaluation of the pmf at every k.
        let ln_norm = ln_choose(m, n);
        for k in 0..=n {
            let expected = (ln_choose(alt, k) + ln_choose(m - alt, n - k) - ln_norm).exp();
            assert!(
                (w.get(k as u32) - expected).abs() <= 1e-12 + 1e-12 * expected,
                "k={k} got {} expected {expected}",
                w.get(k as u32)
            );
        }
    }

    #[test]
    fn extreme_skew_still_normalizes() {
        // m=1000, alt=1: the mass concentrates on k=0..2 but must sum to 1.
        let w = project(1000, 1, 500).unwrap();
        assert!((w.sum() - 1.0).abs() < 1e-12);
        assert!(w.get(0) > 0.4);
        assert!(w.get(1) > 0.4);
    }

    #[test]
    fn symmetry_of_balanced_sites() {
        // alt = m − alt makes the projected SFS symmetric: P(k) = P(n−k).
        let w = project(40, 20, 10).unwrap();
        for k in 0..=10_u32 {
            assert!((w.get(k) - w.get(10 - k)).abs() < 1e-12);
        }
    }
}
