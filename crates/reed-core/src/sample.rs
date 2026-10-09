//! Sample table and sample masks (`docs/spec/architecture.md` §5.3).

use crate::bitset::Bitset;
use crate::error::CoreError;
use crate::id::SampleId;

/// Bidirectional lookup between sample names and header-order sample ids.
#[derive(Debug, Clone, Default)]
pub struct SampleTable {
    names: Vec<String>,
    index: std::collections::HashMap<String, SampleId>,
}

impl SampleTable {
    /// Builds the table from header-order sample names.
    ///
    /// Duplicate names are a header error and are rejected here.
    pub fn from_names<I, S>(names: I) -> Result<Self, CoreError>
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        let names: Vec<String> = names.into_iter().map(Into::into).collect();
        let mut index = std::collections::HashMap::with_capacity(names.len());
        for (sample_idx, name) in names.iter().enumerate() {
            let sample_id = SampleId(
                u32::try_from(sample_idx)
                    .map_err(|_| CoreError::PolicyRejected("sample count exceeds u32::MAX"))?,
            );
            if index.insert(name.clone(), sample_id).is_some() {
                return Err(CoreError::DuplicateSample {
                    sample_name: name.clone(),
                });
            }
        }
        Ok(Self { names, index })
    }

    /// Number of samples.
    #[must_use]
    pub fn n_samples(&self) -> usize {
        self.names.len()
    }

    /// Name of a sample id.
    #[must_use]
    pub fn name(&self, sample_id: SampleId) -> Option<&str> {
        self.names.get(sample_id.get() as usize).map(String::as_str)
    }

    /// Id of a sample name.
    #[must_use]
    pub fn get(&self, name: &str) -> Option<SampleId> {
        self.index.get(name).copied()
    }
}

/// Which samples are selected, indexed by header sample order.
///
/// The length always equals the sample count: there is no implicit
/// "shorter mask keeps the rest" default (`docs/spec/architecture.md` §5.3,
/// fixing legacy issue A13). The word "keep" is reserved for sample-level
/// retention; genotype-level missingness uses "mask".
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SampleMask {
    bits: Bitset,
    n_selected: u32,
}

impl SampleMask {
    /// A mask selecting all `n_samples` samples.
    #[must_use]
    pub fn all(n_samples: usize) -> Self {
        let mut bits = Bitset::new(n_samples);
        for w in bits.words_mut() {
            *w = u64::MAX;
        }
        bits.clear_padding();
        Self {
            n_selected: u32::try_from(n_samples).unwrap_or(u32::MAX),
            bits,
        }
    }

    /// A mask selecting no samples.
    #[must_use]
    pub fn none(n_samples: usize) -> Self {
        Self {
            bits: Bitset::new(n_samples),
            n_selected: 0,
        }
    }

    /// A mask selecting exactly the given sample ids.
    ///
    /// Errors when a sample id is out of range.
    pub fn from_ids(
        n_samples: usize,
        ids: impl IntoIterator<Item = SampleId>,
    ) -> Result<Self, CoreError> {
        let mut mask = Self::none(n_samples);
        for sample_id in ids {
            let idx = sample_id.get() as usize;
            if idx >= n_samples {
                return Err(CoreError::LengthMismatch {
                    expected: n_samples,
                    got: idx + 1,
                });
            }
            mask.bits.set(idx, true);
            mask.n_selected += 1;
        }
        Ok(mask)
    }

    /// Whether the mask selects the sample.
    #[must_use]
    pub fn contains(&self, sample_id: SampleId) -> bool {
        self.bits.get(sample_id.get() as usize)
    }

    /// Number of selected samples.
    #[must_use]
    pub const fn n_selected(&self) -> u32 {
        self.n_selected
    }

    /// Total length (the sample count this mask is bound to).
    #[must_use]
    pub const fn n_samples(&self) -> usize {
        self.bits.len()
    }

    /// Word-wise intersection; both masks must cover the same sample count.
    pub fn and(&self, other: &SampleMask) -> Result<SampleMask, CoreError> {
        let bits = self.bits.and(&other.bits)?;
        Ok(Self {
            n_selected: u32::try_from(bits.count_ones())
                .map_err(|_| CoreError::PolicyRejected("sample count exceeds u32::MAX"))?,
            bits,
        })
    }

    /// Backing words for the counting kernels (padding bits are zero).
    #[must_use]
    pub fn words(&self) -> &[u64] {
        self.bits.words()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_none_and_ids() {
        assert_eq!(SampleMask::all(5).n_selected(), 5);
        assert_eq!(SampleMask::none(5).n_selected(), 0);
        let m = SampleMask::from_ids(5, [SampleId(1), SampleId(4)]).unwrap();
        assert_eq!(m.n_selected(), 2);
        assert!(!m.contains(SampleId(0)));
        assert!(m.contains(SampleId(1)));
        assert!(m.contains(SampleId(4)));
        assert!(SampleMask::from_ids(5, [SampleId(5)]).is_err());
    }

    #[test]
    fn all_mask_word_padding() {
        let m = SampleMask::all(70);
        assert_eq!(m.words().len(), 2);
        // The second word must not set padding bits beyond bit 70.
        assert_eq!(m.words()[1], (1 << 6) - 1);
    }

    #[test]
    fn and_requires_same_length() {
        let a = SampleMask::all(4);
        let b = SampleMask::all(5);
        assert!(a.and(&b).is_err());
        let c = SampleMask::from_ids(4, [SampleId(0), SampleId(3)]).unwrap();
        assert_eq!(a.and(&c).unwrap().n_selected(), 2);
    }

    #[test]
    fn table_lookup() {
        let t = SampleTable::from_names(["a", "b", "c"]).unwrap();
        assert_eq!(t.n_samples(), 3);
        assert_eq!(t.get("b"), Some(SampleId(1)));
        assert_eq!(t.name(SampleId(2)), Some("c"));
        assert!(SampleTable::from_names(["a", "a"]).is_err());
    }
}
