//! Population assignment: which samples belong to which group.

use std::fmt;

/// Errors from constructing a [`PopulationView`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PopError {
    /// `n_groups` was zero.
    NoGroups,
    /// Assignment length does not match the sample count.
    AssignmentLength { expected: usize, got: usize },
    /// A sample was assigned to a group index outside `0..n_groups`.
    GroupOutOfRange { sample: usize, group: u32, n_groups: u32 },
    /// A declared group has no samples.
    EmptyGroup(u32),
}

impl fmt::Display for PopError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            PopError::NoGroups => write!(f, "n_groups must be at least 1"),
            PopError::AssignmentLength { expected, got } => {
                write!(f, "assignment length {got} does not match {expected} samples")
            }
            PopError::GroupOutOfRange { sample, group, n_groups } => write!(
                f,
                "sample {sample} assigned to group {group}, but only {n_groups} groups exist"
            ),
            PopError::EmptyGroup(g) => write!(f, "group {g} has no samples"),
        }
    }
}

impl std::error::Error for PopError {}

/// A partition of sample indices into named-by-index groups.
///
/// Sample order matches the VCF header order; the reader is responsible for
/// mapping sample names to indices once, at construction time.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PopulationView {
    groups: Vec<Vec<u32>>,
}

impl PopulationView {
    /// Builds a view from a per-sample group assignment.
    pub fn from_assignment(
        n_samples: usize,
        assignment: &[u32],
        n_groups: u32,
    ) -> Result<Self, PopError> {
        if n_groups == 0 {
            return Err(PopError::NoGroups);
        }
        if assignment.len() != n_samples {
            return Err(PopError::AssignmentLength {
                expected: n_samples,
                got: assignment.len(),
            });
        }
        let mut groups = vec![Vec::new(); n_groups as usize];
        for (sample, &group) in assignment.iter().enumerate() {
            if group >= n_groups {
                return Err(PopError::GroupOutOfRange {
                    sample,
                    group,
                    n_groups,
                });
            }
            groups[group as usize].push(sample as u32);
        }
        for (g, samples) in groups.iter().enumerate() {
            if samples.is_empty() {
                return Err(PopError::EmptyGroup(g as u32));
            }
        }
        Ok(Self { groups })
    }

    /// Number of groups.
    pub fn n_groups(&self) -> usize {
        self.groups.len()
    }

    /// Sample indices of a group.
    ///
    /// # Panics
    ///
    /// Panics if `group` is out of range.
    pub fn samples(&self, group: u32) -> &[u32] {
        &self.groups[group as usize]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_partition() {
        let view = PopulationView::from_assignment(4, &[0, 0, 1, 1], 2).unwrap();
        assert_eq!(view.n_groups(), 2);
        assert_eq!(view.samples(0), &[0, 1]);
        assert_eq!(view.samples(1), &[2, 3]);
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(
            PopulationView::from_assignment(2, &[0, 1], 0),
            Err(PopError::NoGroups)
        );
        assert_eq!(
            PopulationView::from_assignment(3, &[0, 1], 2),
            Err(PopError::AssignmentLength { expected: 3, got: 2 })
        );
        assert_eq!(
            PopulationView::from_assignment(2, &[0, 2], 2),
            Err(PopError::GroupOutOfRange { sample: 1, group: 2, n_groups: 2 })
        );
        assert_eq!(
            PopulationView::from_assignment(2, &[0, 0], 2),
            Err(PopError::EmptyGroup(1))
        );
    }
}
