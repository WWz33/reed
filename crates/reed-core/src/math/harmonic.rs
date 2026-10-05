//! Harmonic number tables for Watterson's θ and Tajima's D.
//!
//! `a1(n) = H(n−1) = Σ_{k=1}^{n−1} 1/k` and
//! `a2(n) = Σ_{k=1}^{n−1} 1/k²` (Tajima 1989, equations A1–A3).
//!
//! Values are accumulated with [`NeumaierSum`] once into a grow-on-demand
//! table, so every caller sees bit-identical constants.

use std::sync::OnceLock;

use super::sum::NeumaierSum;

/// Table covers `1..=MAX_N`; larger arguments fall back to direct (still
/// compensated) summation.
const MAX_N: usize = 1 << 16;

struct HarmonicTable {
    h1: Vec<f64>,
    h2: Vec<f64>,
}

fn table() -> &'static HarmonicTable {
    static TABLE: OnceLock<HarmonicTable> = OnceLock::new();
    TABLE.get_or_init(|| {
        let mut h1 = Vec::with_capacity(MAX_N + 1);
        let mut h2 = Vec::with_capacity(MAX_N + 1);
        h1.push(0.0);
        h2.push(0.0);
        let mut s1 = NeumaierSum::new();
        let mut s2 = NeumaierSum::new();
        for k in 1..=MAX_N {
            // k ≤ MAX_N = 2^16; exact in f64.
            #[allow(clippy::cast_precision_loss)]
            let kf = k as f64;
            s1.add(1.0 / kf);
            s2.add(1.0 / (kf * kf));
            h1.push(s1.value());
            h2.push(s2.value());
        }
        HarmonicTable { h1, h2 }
    })
}

/// `H(n)` = `1 + 1/2 + … + 1/n` (`H(0) = 0`).
// Harmonic sums are only meaningful for `n` far below f64's exact-integer
// range (2^53); larger `usize` arguments are unreachable in practice.
#[allow(clippy::cast_precision_loss)]
#[must_use]
pub fn harmonic1(n: usize) -> f64 {
    if n <= MAX_N {
        table().h1[n]
    } else {
        let mut s = NeumaierSum::new();
        for k in 1..=n {
            s.add(1.0 / k as f64);
        }
        s.value()
    }
}

/// `Σ_{k=1}^{n} 1/k²` (`0` for `n = 0`).
#[allow(clippy::cast_precision_loss)]
#[must_use]
pub fn harmonic2(n: usize) -> f64 {
    if n <= MAX_N {
        table().h2[n]
    } else {
        let mut s = NeumaierSum::new();
        for k in 1..=n {
            s.add(1.0 / (k as f64 * k as f64));
        }
        s.value()
    }
}

/// Watterson's `a1(n) = H(n−1)` (`0` for `n < 2`).
#[must_use]
pub fn a1(n: usize) -> f64 {
    if n < 2 {
        0.0
    } else {
        harmonic1(n - 1)
    }
}

/// Tajima's `a2(n) = Σ_{k=1}^{n−1} 1/k²` (`0` for `n < 2`).
#[must_use]
pub fn a2(n: usize) -> f64 {
    if n < 2 {
        0.0
    } else {
        harmonic2(n - 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_values() {
        assert_eq!(harmonic1(0), 0.0);
        assert_eq!(harmonic1(1), 1.0);
        assert!((harmonic1(2) - 1.5).abs() < 1e-15);
        assert!((harmonic1(4) - 25.0 / 12.0).abs() < 1e-15);
        assert!((harmonic2(2) - 1.25).abs() < 1e-15);
        assert!((a1(3) - 1.5).abs() < 1e-15);
        assert_eq!(a1(0), 0.0);
        assert_eq!(a1(1), 0.0);
        assert_eq!(a2(1), 0.0);
    }

    #[test]
    fn beyond_table_matches_direct() {
        // MAX_N is 65536; compare a table value against a fresh direct sum.
        let n = 65_537;
        let mut s = NeumaierSum::new();
        for k in 1..=n {
            s.add(1.0 / k as f64);
        }
        assert!((harmonic1(n) - s.value()).abs() <= 4.0 * f64::EPSILON * s.value().abs());
    }
}
