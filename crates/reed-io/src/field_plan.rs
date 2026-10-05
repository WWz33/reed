//! Field plans: which columns a pass needs decoded
//! (`docs/architecture-spec.md` §6.3). Columns outside the plan are never
//! parsed — the head predicates of a filter never touch the sample region.

use bitflags::bitflags;
use reed_core::TagId;

bitflags! {
    /// The fixed columns a pass reads.
    #[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
    pub struct FixedColumns: u16 {
        /// CHROM.
        const CHROM = 1;
        /// POS.
        const POS = 1 << 1;
        /// ID.
        const ID = 1 << 2;
        /// REF.
        const REF = 1 << 3;
        /// ALT.
        const ALT = 1 << 4;
        /// QUAL.
        const QUAL = 1 << 5;
        /// FILTER.
        const FILTER = 1 << 6;
        /// INFO.
        const INFO = 1 << 7;
    }
}

/// How much of the GT column a pass needs.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum GtNeed {
    /// GT is not read at all.
    #[default]
    None,
    /// Genotype bit planes (4 per sample per site).
    Planes,
    /// Planes plus the phased plane.
    PlanesWithPhase,
}

impl GtNeed {
    /// The stricter of two needs.
    #[must_use]
    pub fn max(self, other: Self) -> Self {
        match (self, other) {
            (Self::PlanesWithPhase, _) | (_, Self::PlanesWithPhase) => Self::PlanesWithPhase,
            (Self::Planes, _) | (_, Self::Planes) => Self::Planes,
            (Self::None, Self::None) => Self::None,
        }
    }

    /// Whether planes are needed at all.
    #[must_use]
    pub const fn needs_planes(&self) -> bool {
        !matches!(self, Self::None)
    }
}

/// The set of columns one pass (or the union of passes) needs decoded.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct FieldPlan {
    /// Fixed columns required.
    pub fixed: FixedColumns,
    /// Required INFO tags (header-definition order ids).
    pub info: Vec<TagId>,
    /// Required FORMAT keys other than GT (header-definition order ids).
    pub format: Vec<TagId>,
    /// GT requirement.
    pub gt: GtNeed,
}

impl FieldPlan {
    /// The union of two plans (a pass over the merged stream needs every
    /// column either consumer needs).
    #[must_use]
    pub fn union(&self, other: &Self) -> Self {
        let mut info = self.info.clone();
        for tag in &other.info {
            if !info.contains(tag) {
                info.push(*tag);
            }
        }
        let mut format = self.format.clone();
        for tag in &other.format {
            if !format.contains(tag) {
                format.push(*tag);
            }
        }
        Self {
            fixed: self.fixed.union(other.fixed),
            info,
            format,
            gt: self.gt.max(other.gt),
        }
    }

    /// Whether the plan needs the sample region at all (GT or FORMAT keys).
    #[must_use]
    pub fn needs_samples(&self) -> bool {
        self.gt.needs_planes() || !self.format.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn union_takes_the_stricter_need() {
        let a = FieldPlan {
            fixed: FixedColumns::POS | FixedColumns::QUAL,
            info: vec![TagId(1)],
            format: vec![TagId(4)],
            gt: GtNeed::Planes,
        };
        let b = FieldPlan {
            fixed: FixedColumns::POS | FixedColumns::FILTER,
            info: vec![TagId(2)],
            format: vec![TagId(4)],
            gt: GtNeed::PlanesWithPhase,
        };
        let u = a.union(&b);
        assert_eq!(u.fixed, FixedColumns::POS | FixedColumns::QUAL | FixedColumns::FILTER);
        assert_eq!(u.info, vec![TagId(1), TagId(2)]);
        assert_eq!(u.format, vec![TagId(4)]);
        assert_eq!(u.gt, GtNeed::PlanesWithPhase);
        assert!(u.needs_samples());
    }

    #[test]
    fn head_only_plan_skips_samples() {
        let plan = FieldPlan {
            fixed: FixedColumns::POS,
            ..FieldPlan::default()
        };
        assert!(!plan.needs_samples());
    }
}
