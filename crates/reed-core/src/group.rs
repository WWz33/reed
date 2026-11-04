//! Named population groups over sample masks (`docs/spec/architecture.md` §5.5).
//!
//! Supersedes the legacy `PopulationView` (which required a complete,
//! non-overlapping partition): groups here are named, may overlap (e.g. an
//! auto-generated `all` group) and samples may belong to no group at all.
//! "Population" is the documentation word; the type-level word is "group".

use std::io::BufRead;
use std::path::Path;

use crate::error::CoreError;
use crate::id::{GroupId, SampleId};
use crate::sample::{SampleMask, SampleTable};

/// Options for reading a group file.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct GroupFileOptions {
    /// When set, group-file rows naming samples absent from the VCF header are
    /// skipped with a warning instead of failing; the caller reports the count.
    pub ignore_missing_samples: bool,
}

/// One named group: its id, name and the mask of its member samples.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Group {
    /// Dense id (index into the group table).
    pub id: GroupId,
    /// Group name as given by the assignment.
    pub name: String,
    /// Member samples.
    pub mask: SampleMask,
}

impl Group {
    /// Number of member samples.
    #[must_use]
    pub const fn n_samples(&self) -> u32 {
        self.mask.n_selected()
    }
}

/// A table of named groups over a sample set.
#[derive(Debug, Clone, Default)]
pub struct GroupTable {
    groups: Vec<Group>,
    n_samples: usize,
}

impl GroupTable {
    /// Builds a table from an ordered list of `(name, member sample ids)`
    /// pairs. Duplicate group names are rejected; empty groups are rejected
    /// (a group with no samples cannot support any denominator).
    pub fn from_members<I, S>(n_samples: usize, groups: I) -> Result<Self, CoreError>
    where
        I: IntoIterator<Item = (S, Vec<SampleId>)>,
        S: Into<String>,
    {
        let mut table = Self {
            groups: Vec::new(),
            n_samples,
        };
        let mut seen = std::collections::HashSet::new();
        for (name, member_ids) in groups {
            let name = name.into();
            if !seen.insert(name.clone()) {
                return Err(CoreError::DuplicateGroup { group_name: name });
            }
            let mask = SampleMask::from_ids(n_samples, member_ids)?;
            if mask.n_selected() == 0 {
                return Err(CoreError::EmptyGroup { group_name: name });
            }
            let group_id = GroupId(
                u16::try_from(table.groups.len())
                    .map_err(|_| CoreError::PolicyRejected("group count exceeds u16::MAX"))?,
            );
            table.groups.push(Group {
                id: group_id,
                name,
                mask,
            });
        }
        Ok(table)
    }

    /// Reads a `sample<TAB>population` TSV (the conventional used
    /// by `clean.3pop.txt`): one assignment per line, `#` comments and blank
    /// lines skipped. Group ids follow first-appearance order of the group
    /// name in the file.
    ///
    /// This is the one sanctioned I/O entry point in reed-core; it exists so
    /// every consumer shares one group-file parser.
    pub fn from_tsv(
        path: &Path,
        samples: &SampleTable,
        opts: GroupFileOptions,
    ) -> Result<Self, CoreError> {
        let file = std::fs::File::open(path).map_err(|err| CoreError::GroupFile {
            message: format!("cannot open {}: {err}", path.display()),
        })?;
        let mut names_in_order: Vec<String> = Vec::new();
        let mut members: Vec<Vec<SampleId>> = Vec::new();
        let mut n_ignored = 0_u64;
        for line in std::io::BufReader::new(file).lines() {
            let line = line.map_err(|err| CoreError::GroupFile {
                message: format!("read error in {}: {err}", path.display()),
            })?;
            let line = line.trim_end_matches(['\r', '\n']);
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let (sample_name, group_name) =
                line.split_once('\t').ok_or_else(|| CoreError::GroupFile {
                    message: format!(
                        "{}: expected `sample<TAB>population`, got {line:?}",
                        path.display()
                    ),
                })?;
            let Some(sample_id) = samples.get(sample_name.trim()) else {
                if opts.ignore_missing_samples {
                    n_ignored += 1;
                    continue;
                }
                return Err(CoreError::UnknownSample {
                    sample_name: sample_name.trim().to_string(),
                });
            };
            let group_name = group_name.trim().to_string();
            if let Some(group_idx) = names_in_order.iter().position(|n| *n == group_name) {
                members[group_idx].push(sample_id);
            } else {
                names_in_order.push(group_name);
                members.push(vec![sample_id]);
            }
        }
        let table =
            Self::from_members(samples.n_samples(), names_in_order.into_iter().zip(members))?;
        if n_ignored > 0 {
            tracing::warn!(
                n_ignored,
                "group file referenced samples absent from the header"
            );
        }
        Ok(table)
    }

    /// Builds a table from `(sample name, group name)` pairs.
    pub fn from_assignment(
        samples: &SampleTable,
        pairs: &[(&str, &str)],
    ) -> Result<Self, CoreError> {
        let mut names_in_order: Vec<String> = Vec::new();
        let mut members: Vec<Vec<SampleId>> = Vec::new();
        for &(sample_name, group_name) in pairs {
            let sample_id = samples.get(sample_name).ok_or(CoreError::UnknownSample {
                sample_name: sample_name.to_string(),
            })?;
            let group_name = group_name.to_string();
            if let Some(group_idx) = names_in_order.iter().position(|n| *n == group_name) {
                members[group_idx].push(sample_id);
            } else {
                names_in_order.push(group_name);
                members.push(vec![sample_id]);
            }
        }
        Self::from_members(samples.n_samples(), names_in_order.into_iter().zip(members))
    }

    /// Restricts every group mask to the kept samples (cohort filtering);
    /// group ids are unchanged and emptied groups are kept so existing
    /// `GroupId`s stay valid (statistics over them will yield `None`).
    ///
    /// # Errors
    ///
    /// [`CoreError::LengthMismatch`] when the table and the keep mask were
    /// built for different sample counts. The old fallback kept the original
    /// (unfiltered) group masks and reported success, so every statistic ran
    /// on the wrong cohort without a diagnostic (doc_v25 V25-07b).
    pub fn restrict(&self, keep: &SampleMask) -> Result<GroupTable, CoreError> {
        let groups = self
            .groups
            .iter()
            .map(|group| {
                let mask = group.mask.and(keep)?;
                Ok(Group {
                    id: group.id,
                    name: group.name.clone(),
                    mask,
                })
            })
            .collect::<Result<Vec<_>, CoreError>>()?;
        Ok(Self {
            groups,
            n_samples: self.n_samples,
        })
    }

    /// Number of groups.
    #[must_use]
    pub fn n_groups(&self) -> usize {
        self.groups.len()
    }

    /// Number of samples the table was built for.
    #[must_use]
    pub const fn n_samples(&self) -> usize {
        self.n_samples
    }

    /// Group by id.
    #[must_use]
    pub fn get(&self, group_id: GroupId) -> Option<&Group> {
        self.groups.get(group_id.get() as usize)
    }

    /// Group id by name.
    #[must_use]
    pub fn by_name(&self, name: &str) -> Option<GroupId> {
        self.groups.iter().find(|g| g.name == name).map(|g| g.id)
    }

    /// All groups in id order.
    //
    // No `IntoIterator for &GroupTable`: adding one would extend the public
    // API surface this iteration cycle; `iter()` stays the single accessor.
    #[allow(clippy::iter_without_into_iter)]
    pub fn iter(&self) -> std::slice::Iter<'_, Group> {
        self.groups.iter()
    }

    /// Word masks of every group, for the counting kernels.
    #[must_use]
    pub fn masks(&self) -> Vec<SampleMask> {
        self.groups.iter().map(|g| g.mask.clone()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table() -> (SampleTable, GroupTable) {
        let samples = SampleTable::from_names(["s1", "s2", "s3", "s4"]).unwrap();
        let groups =
            GroupTable::from_assignment(&samples, &[("s1", "A"), ("s2", "A"), ("s3", "B")])
                .unwrap();
        (samples, groups)
    }

    #[test]
    fn assignment_and_lookup() {
        let (_, groups) = table();
        assert_eq!(groups.n_groups(), 2);
        let a = groups.by_name("A").unwrap();
        let b = groups.by_name("B").unwrap();
        assert_eq!(groups.get(a).unwrap().n_samples(), 2);
        assert_eq!(groups.get(b).unwrap().n_samples(), 1);
        assert!(groups.by_name("C").is_none());
        // s4 belongs to no group.
        assert!(!groups.get(a).unwrap().mask.contains(SampleId(3)));
        assert!(!groups.get(b).unwrap().mask.contains(SampleId(3)));
    }

    #[test]
    fn overlapping_and_unassigned_samples_allowed() {
        let _samples = SampleTable::from_names(["s1", "s2"]).unwrap();
        let groups = GroupTable::from_members(
            2,
            [
                ("all", vec![SampleId(0), SampleId(1)]),
                ("sub", vec![SampleId(0)]),
            ],
        )
        .unwrap();
        assert_eq!(groups.n_groups(), 2);
        assert_eq!(groups.by_name("sub").unwrap().get(), 1);
    }

    #[test]
    fn restrict_keeps_group_ids() {
        let (_, groups) = table();
        let keep = SampleMask::from_ids(4, [SampleId(0), SampleId(2)]).unwrap();
        let reduced = groups.restrict(&keep).expect("same sample count");
        assert_eq!(reduced.n_groups(), groups.n_groups());
        let a = reduced.by_name("A").unwrap();
        let b = reduced.by_name("B").unwrap();
        assert_eq!(reduced.get(a).unwrap().n_samples(), 1);
        assert_eq!(reduced.get(b).unwrap().n_samples(), 1);
    }

    /// doc_v25 V25-07b: a mask built for a different sample count is an
    /// error, not a silently unfiltered table.
    #[test]
    fn restrict_rejects_a_mask_of_a_different_length() {
        let (_, groups) = table();
        let keep = SampleMask::all(groups.n_samples() + 1);
        let err = groups.restrict(&keep).expect_err("length mismatch");
        assert!(matches!(
            err,
            CoreError::LengthMismatch { expected, got }
                if expected == groups.n_samples() && got == groups.n_samples() + 1
        ));
    }

    #[test]
    fn rejects_duplicates_and_empties() {
        let _samples = SampleTable::from_names(["s1", "s2"]).unwrap();
        assert!(
            GroupTable::from_members(2, [("A", vec![SampleId(0)]), ("A", vec![SampleId(1)]),],)
                .is_err()
        );
        assert!(GroupTable::from_members(2, [("A", vec![SampleId(0)]), ("B", vec![])]).is_err());
    }

    #[test]
    fn unknown_sample_rejected() {
        let samples = SampleTable::from_names(["s1"]).unwrap();
        assert!(GroupTable::from_assignment(&samples, &[("ghost", "A")]).is_err());
    }
}
