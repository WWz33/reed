//! FORMAT-string layout cache (`docs/rounds/v2/03-architecture.md` §3.2): maps a
//! FORMAT string to the subfield index of GT and of every planned field,
//! plus the last subfield any consumer reads.
//!
//! Layouts are interned once per distinct FORMAT string and handed out as an
//! index, so the scanner can hold a `&Layout` view of the cache while writing
//! its own output buffers (no per-record clone, fixing 01-review D10).

use std::borrow::Cow;
use std::collections::HashMap;

/// One planned integer FORMAT column. The policy is **explicit**: the
/// fallback is chosen by the consumer (a plan's depth mask, a stats report),
/// never inferred from the key's name (doc_v22 D22-1).
///
/// The primary key is always one integer per sample; a fallback is a
/// comma-separated vector summed per sample (`DP` → `AD`). A record whose
/// FORMAT carries `key` never falls back, even when a sample's value is
/// missing.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IntSpec {
    /// The FORMAT key this column reads.
    pub key: Cow<'static, [u8]>,
    /// Key to read when `key` is absent from a record's FORMAT; its values
    /// are summed.
    pub fallback: Option<Cow<'static, [u8]>>,
}

impl IntSpec {
    /// A plain scalar column.
    #[must_use]
    pub fn scalar(key: &'static [u8]) -> Self {
        Self {
            key: Cow::Borrowed(key),
            fallback: None,
        }
    }

    /// A scalar column with a summed fallback key.
    #[must_use]
    pub fn scalar_with_fallback(key: &'static [u8], fallback: &'static [u8]) -> Self {
        Self {
            key: Cow::Borrowed(key),
            fallback: Some(Cow::Borrowed(fallback)),
        }
    }

    /// The sequencing-depth column used by `--min-dp`/`--max-dp` and the
    /// depth reports: `DP`, or the sum of the sample's `AD` values when the
    /// record's FORMAT has no DP key (bcftools `calc_sample_depth`).
    #[must_use]
    pub fn depth() -> Self {
        Self::scalar_with_fallback(b"DP", b"AD")
    }

    /// A column built from owned names (a compiled plan's fields).
    #[must_use]
    pub fn owned(key: Vec<u8>, fallback: Option<Vec<u8>>) -> Self {
        Self {
            key: Cow::Owned(key),
            fallback: fallback.map(Cow::Owned),
        }
    }
}

/// One integer column's resolved subfield positions in a record's FORMAT.
///
/// The index is a `usize`: FORMAT keys are counted, not stored, and an
/// artificial `u8` ceiling silently turned every key past the 256th into
/// "absent" (V22-R2-1). htslib rejects FORMAT columns with more than 255
/// keys at the text parser; reed accepts them and reads the real position.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IntLayout {
    /// Index of the primary key (`None` = absent from this FORMAT).
    pub index: Option<usize>,
    /// Index of the fallback key (`None` = absent or no fallback declared).
    pub fallback_index: Option<usize>,
}

/// Subfield layout of one FORMAT string.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Layout {
    /// Index of `GT` among the colon-separated keys (`None` when absent).
    pub gt_index: Option<usize>,
    /// One entry per planned integer column.
    pub ints: Vec<IntLayout>,
    /// Index of each requested float key.
    pub float_indices: Vec<Option<usize>>,
    /// One past the highest subfield index a consumer reads (0 = none): the
    /// scanner stops splitting subfields there and skips to the next tab.
    pub last_needed: usize,
}

/// Per-worker cache of [`Layout`]s, keyed by the raw FORMAT text.
///
/// The cache clears itself at [`CACHE_LIMIT`] distinct FORMAT strings so an
/// adversarial input cannot grow it without bound; a returned index is valid
/// until the next [`FormatLayoutCache::intern`].
pub struct FormatLayoutCache {
    int_specs: Vec<IntSpec>,
    float_keys: Vec<Box<[u8]>>,
    layouts: Vec<Layout>,
    index: HashMap<Box<[u8]>, usize>,
    /// Last-hit cache: consecutive records share their FORMAT string.
    last: (Vec<u8>, usize),
}

/// Distinct FORMAT strings kept before the cache clears itself.
///
/// 1024 (was 64): a real callset commonly carries more than 64 distinct
/// FORMAT strings once per-sample tags vary, and the old limit cleared the
/// whole table on every crossing, so the hit rate collapsed to zero
/// (doc_v23 V23-34). This is a size cap, not an LRU; an adversarial input
/// still cannot grow the cache without bound.
pub const CACHE_LIMIT: usize = 1024;

impl FormatLayoutCache {
    /// Creates a cache for GT plus the planned integer columns (each with
    /// its optional fallback key) and the named float keys.
    #[must_use]
    pub fn new_with_specs(int_specs: &[IntSpec], float_keys: &[&[u8]]) -> Self {
        Self {
            int_specs: int_specs.to_vec(),
            float_keys: float_keys
                .iter()
                .map(|k| k.to_vec().into_boxed_slice())
                .collect(),
            layouts: Vec::new(),
            index: HashMap::new(),
            last: (Vec::new(), usize::MAX),
        }
    }

    /// Number of integer key columns.
    #[must_use]
    pub fn n_int_keys(&self) -> usize {
        self.int_specs.len()
    }

    /// Number of float key columns.
    #[must_use]
    pub fn n_float_keys(&self) -> usize {
        self.float_keys.len()
    }

    /// Interns `format` and returns the layout index.
    ///
    /// The index is valid until the next call that interns a *new* FORMAT
    /// string at [`CACHE_LIMIT`] (which clears the table).
    pub fn intern(&mut self, format: &[u8]) -> usize {
        // The seed sentinel is `(empty, usize::MAX)`: an empty FORMAT column
        // must not match it (doc_v20 R2) — it interns as a GT-less layout.
        if self.last.1 != usize::MAX && self.last.0 == format {
            return self.last.1;
        }
        if let Some(&ix) = self.index.get(format) {
            self.last = (format.to_vec(), ix);
            return ix;
        }
        if self.layouts.len() >= CACHE_LIMIT {
            self.layouts.clear();
            self.index.clear();
        }
        let keys: Vec<&[u8]> = format.split(|&b| b == b':').collect();
        let index_of = |key: &[u8]| -> Option<usize> { keys.iter().position(|k| *k == key) };
        let gt_index = index_of(b"GT");
        let ints: Vec<IntLayout> = self
            .int_specs
            .iter()
            .map(|s| IntLayout {
                index: index_of(&s.key),
                fallback_index: s.fallback.as_deref().and_then(index_of),
            })
            .collect();
        let float_indices: Vec<Option<usize>> =
            self.float_keys.iter().map(|k| index_of(k)).collect();
        let last_needed = gt_index
            .into_iter()
            .chain(
                ints.iter()
                    .flat_map(|i| [i.index, i.fallback_index])
                    .flatten(),
            )
            .chain(float_indices.iter().copied().flatten())
            .map(|i| i + 1)
            .max()
            .unwrap_or(0);
        let ix = self.layouts.len();
        self.layouts.push(Layout {
            gt_index,
            ints,
            float_indices,
            last_needed,
        });
        self.index.insert(format.to_vec().into_boxed_slice(), ix);
        self.last = (format.to_vec(), ix);
        ix
    }

    /// The interned layout.
    ///
    /// # Panics
    ///
    /// Panics for an index this cache never returned.
    #[must_use]
    pub fn get(&self, ix: usize) -> &Layout {
        &self.layouts[ix]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intern_reports_indices_and_last_needed() {
        let mut cache = FormatLayoutCache::new_with_specs(
            &[IntSpec::scalar(b"DP"), IntSpec::scalar(b"GQ")],
            &[],
        );
        let ix = cache.intern(b"GT:DP:GQ:PL");
        let layout = cache.get(ix);
        assert_eq!(layout.gt_index, Some(0));
        assert_eq!(
            layout.ints,
            vec![
                IntLayout {
                    index: Some(1),
                    fallback_index: None,
                },
                IntLayout {
                    index: Some(2),
                    fallback_index: None,
                },
            ]
        );
        assert_eq!(layout.last_needed, 3);
        // Cache hit on the same bytes: same index.
        assert_eq!(cache.intern(b"GT:DP:GQ:PL"), ix);
    }

    #[test]
    fn absent_keys_are_none_and_do_not_extend_last_needed() {
        let mut cache = FormatLayoutCache::new_with_specs(&[IntSpec::scalar(b"DP")], &[b"AB"]);
        let ix = cache.intern(b"GT:AB");
        let layout = cache.get(ix);
        assert_eq!(layout.gt_index, Some(0));
        assert_eq!(layout.ints[0].index, None);
        assert_eq!(layout.float_indices, vec![Some(1)]);
        assert_eq!(layout.last_needed, 2);
    }

    #[test]
    fn fallback_key_extends_last_needed_and_resolves_per_record() {
        let mut cache =
            FormatLayoutCache::new_with_specs(&[IntSpec::depth()], &[]);
        // The primary key is absent: the AD subfield (index 2) is still
        // split because the fallback is read.
        let ix = cache.intern(b"GT:PL:AD");
        let layout = cache.get(ix);
        assert_eq!(layout.ints[0].index, None);
        assert_eq!(layout.ints[0].fallback_index, Some(2));
        assert_eq!(layout.last_needed, 3);
        // A record carrying DP resolves the primary only; the fallback stays
        // visible for the (never taken) branch but is not needed.
        let ix = cache.intern(b"GT:DP:AD");
        let layout = cache.get(ix);
        assert_eq!(layout.ints[0].index, Some(1));
        assert_eq!(layout.ints[0].fallback_index, Some(2));
        assert_eq!(layout.last_needed, 3);
    }

    #[test]
    fn layout_without_plan_keys_stops_early() {
        let mut cache = FormatLayoutCache::new_with_specs(&[], &[]);
        let ix = cache.intern(b"GT:DP:GQ");
        assert_eq!(cache.get(ix).last_needed, 1);
        let ix = cache.intern(b"DP:GQ");
        assert_eq!(cache.get(ix).gt_index, None);
        assert_eq!(cache.get(ix).last_needed, 0);
    }

    /// V22-R2-1: key positions are unbounded. A `u8` conversion used to
    /// return `None` for the 256th key, making GT/GQ/DP/AD silently absent.
    #[test]
    fn keys_past_index_255_keep_their_position() {
        use std::fmt::Write as _;
        let mut filler = String::new();
        for i in 0..256 {
            let _ = write!(filler, "K{i}:");
        }
        // GQ at position 256 (0-based), GT at 257.
        let format = format!("{filler}GQ:GT");
        let mut cache = FormatLayoutCache::new_with_specs(
            &[IntSpec::scalar(b"DP"), IntSpec::depth(), IntSpec::scalar(b"GQ")],
            &[b"AB"],
        );
        let ix = cache.intern(format.as_bytes());
        let layout = cache.get(ix);
        assert_eq!(layout.ints[0].index, None, "DP is absent");
        assert_eq!(layout.ints[1].index, None, "primary is absent");
        assert_eq!(
            layout.ints[1].fallback_index, None,
            "AD is absent, not the 256th key"
        );
        assert_eq!(layout.ints[2].index, Some(256));
        assert_eq!(layout.gt_index, Some(257));
        assert_eq!(layout.float_indices, vec![None]);
        assert_eq!(layout.last_needed, 258);
    }

    /// doc_v23 V23-34: more than the old 64-entry limit is kept. Interning
    /// an early key after filling the cache must not clear the table.
    #[test]
    fn cache_keeps_more_than_the_old_limit() {
        let mut cache = FormatLayoutCache::new_with_specs(&[], &[]);
        let keys: Vec<Vec<u8>> = (0..=64).map(|i| format!("GT:XX{i}").into_bytes()).collect();
        for key in &keys {
            cache.intern(key);
        }
        assert_eq!(cache.layouts.len(), 65, "65 distinct layouts stay resident");
        let first = cache.intern(&keys[0]);
        assert_eq!(first, 0, "the first layout was not evicted");
        assert_eq!(cache.layouts.len(), 65, "a hit must not clear the table");
    }
}
