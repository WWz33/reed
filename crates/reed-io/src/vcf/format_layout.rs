//! FORMAT-string layout cache (`doc_v2/03-architecture.md` §3.2): maps a
//! FORMAT string to the subfield index of GT and of every planned field,
//! plus the last subfield any consumer reads.
//!
//! Layouts are interned once per distinct FORMAT string and handed out as an
//! index, so the scanner can hold a `&Layout` view of the cache while writing
//! its own output buffers (no per-record clone, fixing 01-review D10).

use std::collections::HashMap;

/// Subfield layout of one FORMAT string.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Layout {
    /// Index of `GT` among the colon-separated keys (`None` when absent).
    pub gt_index: Option<u8>,
    /// Index of each requested integer key (aligned with the cache's key
    /// list; `None` = the key is absent from this FORMAT).
    pub int_indices: Vec<Option<u8>>,
    /// Index of each requested float key.
    pub float_indices: Vec<Option<u8>>,
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
    int_keys: Vec<Box<[u8]>>,
    float_keys: Vec<Box<[u8]>>,
    layouts: Vec<Layout>,
    index: HashMap<Box<[u8]>, usize>,
    /// Last-hit cache: consecutive records share their FORMAT string.
    last: (Vec<u8>, usize),
}

/// Distinct FORMAT strings kept before the cache clears itself.
pub const CACHE_LIMIT: usize = 64;

impl FormatLayoutCache {
    /// Creates a cache for GT plus the named integer and float keys, in
    /// output-column order.
    #[must_use]
    pub fn new(int_keys: &[&[u8]], float_keys: &[&[u8]]) -> Self {
        Self {
            int_keys: int_keys.iter().map(|k| k.to_vec().into_boxed_slice()).collect(),
            float_keys: float_keys.iter().map(|k| k.to_vec().into_boxed_slice()).collect(),
            layouts: Vec::new(),
            index: HashMap::new(),
            last: (Vec::new(), usize::MAX),
        }
    }

    /// Number of integer key columns.
    #[must_use]
    pub fn n_int_keys(&self) -> usize {
        self.int_keys.len()
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
        if self.last.0 == format {
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
        let index_of = |key: &[u8]| -> Option<u8> {
            keys.iter()
                .position(|k| *k == key)
                .and_then(|i| u8::try_from(i).ok())
        };
        let gt_index = index_of(b"GT");
        let int_indices: Vec<Option<u8>> = self.int_keys.iter().map(|k| index_of(k)).collect();
        let float_indices: Vec<Option<u8>> = self.float_keys.iter().map(|k| index_of(k)).collect();
        let last_needed = gt_index
            .into_iter()
            .chain(int_indices.iter().copied().flatten())
            .chain(float_indices.iter().copied().flatten())
            .map(|i| usize::from(i) + 1)
            .max()
            .unwrap_or(0);
        let ix = self.layouts.len();
        self.layouts.push(Layout { gt_index, int_indices, float_indices, last_needed });
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
        let mut cache = FormatLayoutCache::new(&[b"DP", b"GQ"], &[]);
        let ix = cache.intern(b"GT:DP:GQ:PL");
        let layout = cache.get(ix);
        assert_eq!(layout.gt_index, Some(0));
        assert_eq!(layout.int_indices, vec![Some(1), Some(2)]);
        assert_eq!(layout.last_needed, 3);
        // Cache hit on the same bytes: same index.
        assert_eq!(cache.intern(b"GT:DP:GQ:PL"), ix);
    }

    #[test]
    fn absent_keys_are_none_and_do_not_extend_last_needed() {
        let mut cache = FormatLayoutCache::new(&[b"DP"], &[b"AB"]);
        let ix = cache.intern(b"GT:AB");
        let layout = cache.get(ix);
        assert_eq!(layout.gt_index, Some(0));
        assert_eq!(layout.int_indices, vec![None]);
        assert_eq!(layout.float_indices, vec![Some(1)]);
        assert_eq!(layout.last_needed, 2);
    }

    #[test]
    fn layout_without_plan_keys_stops_early() {
        let mut cache = FormatLayoutCache::new(&[], &[]);
        let ix = cache.intern(b"GT:DP:GQ");
        assert_eq!(cache.get(ix).last_needed, 1);
        let ix = cache.intern(b"DP:GQ");
        assert_eq!(cache.get(ix).gt_index, None);
        assert_eq!(cache.get(ix).last_needed, 0);
    }
}
