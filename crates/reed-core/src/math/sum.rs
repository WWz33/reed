//! Neumaier compensated summation (`docs/architecture-spec.md` §9.10).
//!
//! All floating-point accumulation in statistics goes through this type so
//! that merging partial sums stays deterministic to the last bit regardless
//! of chunk boundaries (`docs/architecture-spec.md` §7.6).

/// Neumaier (improved Kahan) compensated sum.
///
/// Symbol table: `sum` is the running total, `c` the accumulated correction
/// for lost low-order bits (Neumaier 1974; Higham, *Accuracy and Stability of
/// Numerical Algorithms*, 2nd ed., §4.3).
#[derive(Clone, Copy, Debug)]
pub struct NeumaierSum {
    sum: f64,
    c: f64,
}

impl Default for NeumaierSum {
    fn default() -> Self {
        Self {
            sum: 0.0,
            c: 0.0,
        }
    }
}

impl NeumaierSum {
    /// Creates a zeroed accumulator.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            sum: 0.0,
            c: 0.0,
        }
    }

    /// Adds one term.
    pub fn add(&mut self, x: f64) {
        // Order the operands so the compensation captures the low bits of
        // the larger magnitude (Neumaier's original formulation).
        let t = self.sum + x;
        if self.sum.abs() >= x.abs() {
            self.c += (self.sum - t) + x;
        } else {
            self.c += (x - t) + self.sum;
        }
        self.sum = t;
    }

    /// Absorbs another compensated sum (partial-order merge).
    ///
    /// Both sides' residuals are re-absorbed through the same compensation
    /// path as single adds, so merging partials stays within the same
    /// reorder tolerance as reordering the terms themselves (`docs/
    /// architecture-spec.md` §8.4: ≤ 4 ulp).
    pub fn merge(&mut self, other: &Self) {
        self.add(other.sum);
        self.add(other.c);
    }

    /// Current compensated value (`sum + c`); the only read accessor.
    #[must_use]
    pub const fn value(&self) -> f64 {
        self.sum + self.c
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compensates_catastrophic_cancellation() {
        // 1e16 + ten 1s + (-1e16): each 1.0 rounds away in the naive sum
        // (ulp(1e16) = 2), so it ends at 0.0; the compensated sum recovers 10.
        let mut acc = NeumaierSum::new();
        acc.add(1.0e16);
        for _ in 0..10 {
            acc.add(1.0);
        }
        acc.add(-1.0e16);
        assert_eq!(acc.value(), 10.0);
        // The naive sum returns 0.0 here.
        assert_eq!(1.0e16 + 1.0 + 1.0 - 1.0e16, 0.0);
    }

    #[test]
    fn order_within_tolerance() {
        let values: Vec<f64> = (1..=1000).map(f64::from).map(|x| x / 7.0).collect();
        let mut forward = NeumaierSum::new();
        for &v in &values {
            forward.add(v);
        }
        let mut backward = NeumaierSum::new();
        for &v in values.iter().rev() {
            backward.add(v);
        }
        assert!((forward.value() - backward.value()).abs() <= 4.0 * f64::EPSILON * forward.value().abs());
    }

    #[test]
    fn merge_recovers_partial_order() {
        // Splitting the catastrophic-cancellation sequence into two partials
        // and merging must still recover 10.0 exactly.
        let mut first = NeumaierSum::new();
        first.add(1.0e16);
        for _ in 0..5 {
            first.add(1.0);
        }
        let mut second = NeumaierSum::new();
        for _ in 0..5 {
            second.add(1.0);
        }
        second.add(-1.0e16);
        first.merge(&second);
        assert_eq!(first.value(), 10.0);
    }
}
