//! FORMAT-string layout cache (`docs/rounds/v2/03-architecture.md` §3.2): maps a
//! FORMAT string to the subfield index of GT and of every planned field,
//! plus the last subfield any consumer reads.
//!
//! Layouts are interned once per distinct FORMAT string and handed out as an
//! index, so the scanner can hold a `&Layout` view of the cache while writing
//! its own output buffers (no per-record clone, fixing 01-review D10).

use std::borrow::Cow;
use std::collections::HashMap;

/// Index of `key` among `format_text`'s colon-separated keys.
///
/// The one key-position primitive (doc_v25 V25-02): the layout cache and the
/// non-cache consumers (encoder, convert, pca) all resolve a FORMAT key
/// through it instead of writing their own `split`/`position` loop.
#[must_use]
pub fn key_index(format_text: &[u8], key: &[u8]) -> Option<usize> {
    format_text.split(|&b| b == b':').position(|k| k == key)
}

/// Index of `GT` among the FORMAT keys; `None` when the record has no GT.
///
/// This is the single GT-position lookup (doc_v25 V25-02): the scanner, the
/// filter's edit path, the `-Ob` encoder, `convert` and `pca` all call it.
#[must_use]
pub fn gt_index(format_text: &[u8]) -> Option<usize> {
    key_index(format_text, b"GT")
}

/// The GT subfield of one sample column, resolved by the record's
/// [`gt_index`].
///
/// `None` when the column has fewer subfields than the GT position — the
/// genotype is absent, which consumers treat as missing (never as the first
/// subfield: `FORMAT=DP:GT` with a `30`-only column has no GT text).
#[must_use]
pub fn gt_subfield(sample_col: &[u8], gt_index: Option<usize>) -> Option<&[u8]> {
    gt_index.and_then(|i| subfield(sample_col, i))
}

/// Subfield `key` of one sample column, or `None` when the column has fewer
/// subfields than `key` — the one colon-scan primitive (doc_v28 V28-24).
///
/// Callers that splice the column rather than read it use [`subfield_span`].
#[must_use]
pub fn subfield(sample_col: &[u8], key: usize) -> Option<&[u8]> {
    sample_col.split(|&b| b == b':').nth(key)
}

/// The byte range of subfield `key` within `sample_col`, or `None` when the
/// column has fewer subfields.
///
/// The range is relative to `sample_col`; a caller working in full-line
/// offsets adds the column's own start.
#[must_use]
pub fn subfield_span(sample_col: &[u8], key: usize) -> Option<std::ops::Range<usize>> {
    let mut from = 0_usize;
    for (index, part) in sample_col.split(|&b| b == b':').enumerate() {
        if index == key {
            return Some(from..from + part.len());
        }
        from += part.len() + 1;
    }
    None
}

/// The ranges of the first `needed` subfields of `sample_col` (relative to
/// the column), written into `out` (cleared first).
///
/// The scanner stops splitting at its plan's `last_needed` and the encoder
/// walks every subfield; both want this same colon walk (doc_v28 V28-24).
pub fn subfield_spans(sample_col: &[u8], needed: usize, out: &mut Vec<std::ops::Range<usize>>) {
    out.clear();
    if needed == 0 {
        return;
    }
    // One pass over the (short) field instead of one `memchr` per subfield
    // (doc_v3 03 O6, level 1: the scan is call-bound on ~30-byte fields).
    let mut from = 0_usize;
    for (i, &b) in sample_col.iter().enumerate() {
        if b == b':' {
            out.push(from..i);
            from = i + 1;
            if out.len() == needed {
                return;
            }
        }
    }
    out.push(from..sample_col.len());
}

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

/// A handle to one interned layout. The generation makes ids from a replaced
/// generation detectably stale instead of aliasing a newer entry.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LayoutId {
    epoch: u32,
    index: u32,
}

/// One cache generation: its own key map and layout storage.
#[derive(Default)]
struct Generation {
    epoch: u32,
    layouts: Vec<Layout>,
    index: HashMap<Box<[u8]>, usize>,
}

/// Per-worker, two-generation cache of [`Layout`]s keyed by the raw FORMAT
/// text.
///
/// Each generation holds [`CACHE_LIMIT`] distinct FORMAT strings. When the
/// current generation is full it becomes the previous one and a fresh
/// current starts (`doc_v24 V24-06`); the old previous is dropped. The
/// cache therefore holds at most `2 * CACHE_LIMIT` layouts and never
/// "clears everything" on crossing a threshold — a one-shot scan of
/// 65+ distinct FORMAT strings keeps its recent working set instead of
/// losing the hit rate entirely.
///
/// A lookup checks current, then previous. Previous hits are served but not
/// promoted (no per-hit allocation); an id stays valid until the generation
/// holding it is replaced, i.e. at most two rotations.
pub struct FormatLayoutCache {
    int_specs: Vec<IntSpec>,
    float_keys: Vec<Box<[u8]>>,
    current: Generation,
    previous: Generation,
    /// Last-hit cache: consecutive records share their FORMAT string.
    last: (Vec<u8>, Option<LayoutId>),
    /// Generation number for the next rotation (starts at 1).
    next_generation: u32,
}

/// Distinct FORMAT strings per generation (two generations are resident).
pub const CACHE_LIMIT: usize = 256;

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
            current: Generation::default(),
            previous: Generation::default(),
            last: (Vec::new(), None),
            next_generation: 1,
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

    /// Interns `format` and returns the layout handle.
    ///
    /// The handle is valid while its generation is resident (current or
    /// previous); after two further rotations [`Self::get`] panics rather
    /// than returning another entry.
    pub fn intern(&mut self, format: &[u8]) -> LayoutId {
        // The seed sentinel is `(empty, None)`: an empty FORMAT column must
        // not match it (doc_v20 R2) — it interns as a GT-less layout.
        if let Some(id) = self.last.1
            && self.last.0 == format
        {
            return id;
        }
        if let Some(&ix) = self.current.index.get(format) {
            let id = LayoutId {
                epoch: self.current.epoch,
                index: u32::try_from(ix).expect("invariant: bounded by CACHE_LIMIT"),
            };
            self.last = (format.to_vec(), Some(id));
            return id;
        }
        if let Some(&ix) = self.previous.index.get(format) {
            let id = LayoutId {
                epoch: self.previous.epoch,
                index: u32::try_from(ix).expect("invariant: bounded by CACHE_LIMIT"),
            };
            self.last = (format.to_vec(), Some(id));
            return id;
        }
        if self.current.layouts.len() >= CACHE_LIMIT {
            self.previous = std::mem::take(&mut self.current);
            self.current = Generation {
                epoch: self.next_generation,
                ..Generation::default()
            };
            self.next_generation = self.next_generation.wrapping_add(1);
        }
        let index_of = |key: &[u8]| -> Option<usize> { key_index(format, key) };
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
        let ix = self.current.layouts.len();
        self.current.layouts.push(Layout {
            gt_index,
            ints,
            float_indices,
            last_needed,
        });
        self.current
            .index
            .insert(format.to_vec().into_boxed_slice(), ix);
        let id = LayoutId {
            epoch: self.current.epoch,
            index: u32::try_from(ix).expect("invariant: bounded by CACHE_LIMIT"),
        };
        self.last = (format.to_vec(), Some(id));
        id
    }

    /// The interned layout.
    ///
    /// # Panics
    ///
    /// Panics for a handle whose generation has been replaced (stale id), or
    /// one this cache never returned.
    #[must_use]
    pub fn get(&self, id: LayoutId) -> &Layout {
        let layouts = if id.epoch == self.current.epoch {
            &self.current.layouts
        } else if id.epoch == self.previous.epoch {
            &self.previous.layouts
        } else {
            unreachable!("invariant: stale layout id from a replaced generation");
        };
        &layouts[id.index as usize]
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
        let mut cache = FormatLayoutCache::new_with_specs(&[IntSpec::depth()], &[]);
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
            &[
                IntSpec::scalar(b"DP"),
                IntSpec::depth(),
                IntSpec::scalar(b"GQ"),
            ],
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

    /// doc_v24 V24-06: crossing the old 64-entry cliff no longer clears the
    /// table. A 65th distinct FORMAT rotates; the first 64 stay resident in
    /// the previous generation and are still served.
    #[test]
    fn rotation_keeps_the_previous_generation() {
        let mut cache = FormatLayoutCache::new_with_specs(&[], &[]);
        let keys: Vec<Vec<u8>> = (0..CACHE_LIMIT)
            .map(|i| format!("GT:XX{i}").into_bytes())
            .collect();
        let ids: Vec<LayoutId> = keys.iter().map(|k| cache.intern(k)).collect();
        // One more distinct FORMAT rotates current -> previous.
        let extra = cache.intern(b"GT:NEW");
        assert_eq!(cache.current.layouts.len(), 1);
        assert_eq!(cache.previous.layouts.len(), CACHE_LIMIT);
        // Every old id is still readable, and the byte lookup hits too.
        for (key, id) in keys.iter().zip(&ids) {
            assert_eq!(cache.intern(key), *id, "previous-generation hit");
            assert_eq!(cache.get(*id).gt_index, Some(0));
        }
        assert_eq!(cache.get(extra).gt_index, Some(0));
    }

    /// The structure is bounded: many more distinct FORMAT strings than two
    /// generations keep at most `2 * CACHE_LIMIT` layouts resident.
    #[test]
    fn cache_is_bounded_under_distinct_formats() {
        let mut cache = FormatLayoutCache::new_with_specs(&[], &[]);
        for i in 0..(CACHE_LIMIT * 8) {
            let key = format!("GT:MANY{i}");
            cache.intern(key.as_bytes());
        }
        assert!(cache.current.layouts.len() <= CACHE_LIMIT);
        assert!(cache.previous.layouts.len() <= CACHE_LIMIT);
        assert!(cache.current.index.len() <= CACHE_LIMIT);
        assert!(cache.previous.index.len() <= CACHE_LIMIT);
    }

    /// doc_v25 V25-02: one GT-position lookup and one subfield resolver.
    /// `DP:GT` must find the second subfield, and a short column must not
    /// fall back to the first subfield (the `30:0/1` bug).
    #[test]
    fn gt_index_and_subfield_follow_the_format_order() {
        assert_eq!(gt_index(b"GT:DP:GQ"), Some(0));
        assert_eq!(gt_index(b"DP:GT"), Some(1));
        assert_eq!(gt_index(b"AD:DP:GQ:GT"), Some(3));
        assert_eq!(gt_index(b"DP:GQ"), None);
        assert_eq!(gt_index(b""), None);

        assert_eq!(gt_subfield(b"0/1:30", Some(0)), Some(&b"0/1"[..]));
        assert_eq!(gt_subfield(b"30:0/1", Some(1)), Some(&b"0/1"[..]));
        assert_eq!(gt_subfield(b"5,6:10:40:1|2", Some(3)), Some(&b"1|2"[..]));
        assert_eq!(gt_subfield(b".", Some(0)), Some(&b"."[..]));
        // A short column has no GT — never fall back to the first subfield.
        assert_eq!(gt_subfield(b"30", Some(1)), None);
        assert_eq!(gt_subfield(b"30", None), None);
    }

    /// doc_v24 V24-06: ids from a dropped generation are loudly stale, not
    /// silently aliased to a newer layout at the same index.
    #[test]
    #[should_panic(expected = "stale layout id")]
    fn stale_ids_do_not_alias_a_new_generation() {
        let mut cache = FormatLayoutCache::new_with_specs(&[], &[]);
        let old = cache.intern(b"GT:FIRST");
        // Rotation 1: generation 0 becomes previous.
        for i in 0..CACHE_LIMIT {
            let key = format!("GT:A{i}");
            cache.intern(key.as_bytes());
        }
        // Rotation 2: generation 0 is dropped.
        for i in 0..CACHE_LIMIT {
            let key = format!("GT:B{i}");
            cache.intern(key.as_bytes());
        }
        let _ = cache.get(old);
    }
}
