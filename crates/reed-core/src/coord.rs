//! Coordinates (`docs/architecture-spec.md` §5.2).
//!
//! Internal computation uses 0-based half-open intervals ([`Pos0`],
//! [`Interval0`]); VCF text, user regions, error messages and TSV output use
//! 1-based closed coordinates ([`Pos1`], [`Interval1`]). All unit conversion
//! must go through the [`From`] impls in this module — bare `pos - 1`
//! arithmetic at call sites is a review reject.

/// A 0-based position.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Pos0(pub u64);

/// A 1-based position (= VCF `POS`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Pos1(pub u64);

/// A 0-based half-open interval `[start, end)`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Interval0 {
    /// First position in the interval (0-based).
    pub start: Pos0,
    /// Position one past the interval end (0-based).
    pub end: Pos0,
}

/// A 1-based closed interval `[start, end]` (VCF/TSV convention).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Interval1 {
    /// First position in the interval (1-based).
    pub start: Pos1,
    /// Last position in the interval (1-based).
    pub end: Pos1,
}

impl From<Pos1> for Pos0 {
    fn from(pos1: Pos1) -> Self {
        Self(pos1.0 - 1)
    }
}

impl From<Pos0> for Pos1 {
    fn from(pos0: Pos0) -> Self {
        Self(pos0.0 + 1)
    }
}

impl Pos0 {
    /// Converts to the 1-based VCF convention.
    #[must_use]
    pub const fn to_pos1(self) -> Pos1 {
        Pos1(self.0 + 1)
    }
}

impl Pos1 {
    /// Converts to the 0-based internal convention.
    #[must_use]
    pub const fn to_pos0(self) -> Pos0 {
        Pos0(self.0 - 1)
    }
}

impl Interval0 {
    /// Builds an interval; the caller guarantees `start <= end`.
    #[must_use]
    pub const fn new(start: Pos0, end: Pos0) -> Self {
        Self { start, end }
    }

    /// Span in base pairs (`end - start`).
    #[must_use]
    pub const fn span_bp(&self) -> u64 {
        self.end.0.saturating_sub(self.start.0)
    }

    /// Whether the interval contains the position.
    #[must_use]
    pub const fn contains(&self, pos: Pos0) -> bool {
        pos.0 >= self.start.0 && pos.0 < self.end.0
    }

    /// Converts to the 1-based closed convention.
    #[must_use]
    pub const fn to_interval1(&self) -> Interval1 {
        Interval1 {
            start: Pos1(self.start.0 + 1),
            end: Pos1(self.end.0),
        }
    }
}

impl Interval1 {
    /// Converts to the 0-based half-open convention.
    #[must_use]
    pub const fn to_interval0(&self) -> Interval0 {
        Interval0 {
            start: Pos0(self.start.0 - 1),
            end: Pos0(self.end.0),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pos_roundtrip() {
        assert_eq!(Pos1(42).to_pos0(), Pos0(41));
        assert_eq!(Pos0(41).to_pos1(), Pos1(42));
    }

    #[test]
    fn interval_conventions() {
        // 1-based closed [10, 20] == 0-based half-open [9, 20).
        let i1 = Interval1 { start: Pos1(10), end: Pos1(20) };
        let i0 = i1.to_interval0();
        assert_eq!(i0.start, Pos0(9));
        assert_eq!(i0.end, Pos0(20));
        assert_eq!(i0.span_bp(), 11);
        assert!(i0.contains(Pos0(9)));
        assert!(!i0.contains(Pos0(20)));
        assert_eq!(i0.to_interval1(), i1);
    }
}
