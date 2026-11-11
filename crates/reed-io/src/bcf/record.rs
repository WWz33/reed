//! One BCF record's byte layout (hts-specs §6.3): `l_shared`, `l_indiv`,
//! then the shared and individual blocks.
//!
//! Decoding is lazy: the field walkers below hand out borrowed slices, so
//! rendering and consumption read each byte once and allocate nothing.

// Casts in this file are bounded by construction: Record fields are bounded by the 1 GiB record ceiling.
#![allow(clippy::cast_possible_truncation)]

use crate::error::{IoError, Locus};

/// BCF value types (the low nibble of a typed descriptor).
pub mod value_type {
    /// Sentinel placeholder (no payload; a INFO flag).
    pub const NULL: u8 = 0;
    /// 8-bit signed integer.
    pub const INT8: u8 = 1;
    /// 16-bit signed integer.
    pub const INT16: u8 = 2;
    /// 32-bit signed integer.
    pub const INT32: u8 = 3;
    /// 32-bit IEEE float.
    pub const FLOAT: u8 = 5;
    /// NUL-padded character vector.
    pub const CHAR: u8 = 7;
}

/// Byte width of one element of `ty`.
#[must_use]
pub const fn type_width(ty: u8) -> usize {
    match ty {
        value_type::INT8 | value_type::CHAR => 1,
        value_type::INT16 => 2,
        value_type::INT32 | value_type::FLOAT => 4,
        _ => 0,
    }
}

/// The missing and vector-end sentinels of one integer type.
#[must_use]
pub const fn int_sentinels(ty: u8) -> (i64, i64) {
    match ty {
        value_type::INT8 => (i8::MIN as i64, i8::MIN as i64 + 1),
        value_type::INT16 => (i16::MIN as i64, i16::MIN as i64 + 1),
        _ => (i32::MIN as i64, i32::MIN as i64 + 1),
    }
}

/// One integer element of `ty` (both renderers and the column decoder read
/// the same little-endian bytes).
pub(super) fn int_element(ty: u8, elem: &[u8]) -> i64 {
    match ty {
        value_type::INT8 => i64::from(elem[0] as i8),
        value_type::INT16 => i64::from(i16::from_le_bytes([elem[0], elem[1]])),
        _ => i64::from(i32::from_le_bytes([elem[0], elem[1], elem[2], elem[3]])),
    }
}

/// One GT element: `None` at vector-end, `Some(code)` otherwise.
///
/// BCF GT codes are non-negative (`(allele + 1) << 1 | phased`, with 0/1
/// for the missing allele). The only other legal negative value is the
/// type's vector-end sentinel. A negative element anywhere else is corrupt
/// input; the old renderer formatted it as allele 0 and every consumer saw
/// a homozygous REF call (doc_v28 V28-04).
pub(super) fn gt_element(ty: u8, elem: &[u8]) -> Result<Option<i64>, &'static str> {
    if !matches!(ty, value_type::INT8 | value_type::INT16 | value_type::INT32) {
        return Err("BCF GT field is not integer-typed");
    }
    let (_, vector_end) = int_sentinels(ty);
    let value = int_element(ty, elem);
    if value == vector_end {
        Ok(None)
    } else if value < 0 {
        Err("BCF GT element is negative")
    } else {
        Ok(Some(value))
    }
}

/// htslib's legal value ranges (`BCF_MIN_BT_INT*`/`BCF_MAX_BT_INT*`); the
/// reserved sentinel values sit below `INT*_MIN_VALUE`.
pub const INT8_MIN_VALUE: i64 = -120;
/// See [`INT8_MIN_VALUE`].
pub const INT8_MAX_VALUE: i64 = 127;
/// See [`INT8_MIN_VALUE`].
pub const INT16_MIN_VALUE: i64 = -32760;
/// See [`INT8_MIN_VALUE`].
pub const INT16_MAX_VALUE: i64 = 32767;
/// See [`INT8_MIN_VALUE`].
pub const INT32_MIN_VALUE: i64 = -2_147_483_640;
/// See [`INT8_MIN_VALUE`].
pub const INT32_MAX_VALUE: i64 = 2_147_483_647;

/// The float bit patterns for missing and vector-end.
pub const FLOAT_MISSING: u32 = 0x7F80_0001;
/// See [`FLOAT_MISSING`].
pub const FLOAT_VECTOR_END: u32 = 0x7F80_0002;

/// One BCF record's blocks.
#[derive(Clone, Copy, Debug)]
pub struct BcfRecordRef<'a> {
    /// The shared block (site-level fields).
    pub shared: &'a [u8],
    /// The individual block (FORMAT fields × samples).
    pub indiv: &'a [u8],
}

/// A malformed-BCF error with a placeholder locus (the record position is
/// unknown until the shared block is decoded).
#[must_use]
pub fn malformed(reason: &'static str) -> IoError {
    IoError::MalformedRecord {
        locus: Locus {
            contig_name: "<bcf>".into(),
            pos1: 0,
        },
        reason: reason.to_string(),
    }
}

/// One chunk's framing decision against a record prefix carried from the
/// previous chunk (doc_v28 V28-25).
///
/// The filter's BCF input and the index job each stitched chunk seams with
/// their own copy of this arithmetic; they differ in what they retain and
/// decode, but the question "does the length prefix complete, how long is the
/// record, where does the first whole record end" is the same.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FrameStep {
    /// Fewer than the eight length-prefix bytes are available across the
    /// carried prefix and this chunk.
    PrefixIncomplete,
    /// The crossing record declares `total` bytes and does not end inside
    /// this chunk.
    Spans {
        /// The record's declared total length.
        total: usize,
    },
    /// The crossing record ends `body_start` bytes into this chunk; whole
    /// records run `[body_start, body_end)`; the next crossing record starts
    /// at `body_end`.
    Complete {
        /// Offset where this chunk's own whole records begin.
        body_start: usize,
        /// Offset after the last whole record.
        body_end: usize,
    },
}

/// Frames one chunk whose start continues a record carried from the previous
/// chunk.
///
/// `prefix` holds the carried bytes the caller kept (it may be shorter than
/// `prefix_len` when the caller caps its buffer), `prefix_len` is how many
/// bytes of the crossing record precede `text`, and `start` is the offset of
/// the first record when nothing is carried. The function allocates nothing,
/// stores nothing, and knows no virtual offsets or ceilings — those stay with
/// the caller (doc_v28 V28-25).
#[must_use]
pub fn frame_chunk(prefix: &[u8], prefix_len: usize, text: &[u8], start: usize) -> FrameStep {
    if prefix_len == 0 {
        return FrameStep::Complete {
            body_start: start,
            body_end: record_boundary(text, start),
        };
    }
    let Some(total) = declared_len_across(prefix, text) else {
        return FrameStep::PrefixIncomplete;
    };
    let need = total.saturating_sub(prefix_len);
    if need > text.len() {
        return FrameStep::Spans { total };
    }
    FrameStep::Complete {
        body_start: need,
        body_end: record_boundary(text, need),
    }
}

/// The offset after the last whole record starting at `from`; `text.len()`
/// when the tail crosses the seam.
#[must_use]
pub fn record_boundary(text: &[u8], from: usize) -> usize {
    let mut p = from;
    while let Some(len) = record_len(&text[p..]) {
        p += len;
    }
    p
}

/// The record's declared total length, taking the missing prefix bytes from
/// `text` (a chunk may end inside the eight-byte prefix).
fn declared_len_across(prefix: &[u8], text: &[u8]) -> Option<usize> {
    if prefix.len() >= 8 {
        return declared_len(prefix);
    }
    let need = 8 - prefix.len();
    let available = text.len().min(need);
    let mut head = [0_u8; 8];
    head[..prefix.len()].copy_from_slice(prefix);
    head[prefix.len()..prefix.len() + available].copy_from_slice(&text[..available]);
    // A caller that capped its carried buffer may hold fewer than eight
    // bytes while `prefix_len` says more were seen; the length prefix is
    // then simply not recoverable from what is in hand.
    if prefix.len() + available < 8 {
        return None;
    }
    declared_len(&head)
}

/// The total byte length of the record starting at `bytes`, when the whole
/// record is present. This is the length-prefix that frames a BCF stream.
#[must_use]
pub fn record_len(bytes: &[u8]) -> Option<usize> {
    let head = bytes.get(..8)?;
    let l_shared = u32::from_le_bytes([head[0], head[1], head[2], head[3]]) as usize;
    let l_indiv = u32::from_le_bytes([head[4], head[5], head[6], head[7]]) as usize;
    let total = 8_usize.checked_add(l_shared)?.checked_add(l_indiv)?;
    (bytes.len() >= total).then_some(total)
}

/// The total length a record declares in its first eight bytes, whether or
/// not the whole record is present yet.
#[must_use]
pub fn declared_len(head: &[u8]) -> Option<usize> {
    let head = head.get(..8)?;
    let l_shared = u32::from_le_bytes([head[0], head[1], head[2], head[3]]) as usize;
    let l_indiv = u32::from_le_bytes([head[4], head[5], head[6], head[7]]) as usize;
    8_usize.checked_add(l_shared)?.checked_add(l_indiv)
}

/// Splits one complete record's bytes into its blocks.
///
/// # Errors
///
/// [`IoError::MalformedRecord`] when the length prefixes overrun `bytes`.
pub fn record_ref(bytes: &[u8]) -> Result<BcfRecordRef<'_>, IoError> {
    // `record_len` already includes the "is the whole record present" test,
    // so a second `total > bytes.len()` check here could never fire
    // (doc_v28 V28-18).
    record_len(bytes).ok_or_else(|| malformed("BCF record is truncated"))?;
    let head = &bytes[..8];
    let l_shared = u32::from_le_bytes([head[0], head[1], head[2], head[3]]) as usize;
    let l_indiv = u32::from_le_bytes([head[4], head[5], head[6], head[7]]) as usize;
    Ok(BcfRecordRef {
        shared: &bytes[8..8 + l_shared],
        indiv: &bytes[8 + l_shared..8 + l_shared + l_indiv],
    })
}

/// Maps an inflated-text offset within one chunk back to a virtual offset
/// through the chunk's block table (BGZF's `(compressed << 16) | offset`).
///
/// A thin alias over [`crate::bgzf::address::BlockTable`] (doc_v28 V28-27):
/// the BCF record path and the index job used to carry identical copies of
/// this walk.
#[derive(Clone, Copy, Debug)]
pub struct VirtualMap<'a> {
    /// Compressed offset of the chunk's first block.
    pub chunk_start: u64,
    /// Per-block `(compressed length, inflated length)` pairs.
    pub blocks: &'a [(u32, u32)],
}

impl VirtualMap<'_> {
    /// The virtual offset of `text_off` in the chunk's inflated text.
    #[must_use]
    pub fn voff(&self, text_off: usize) -> u64 {
        crate::bgzf::address::BlockTable::new(self.blocks)
            .voff(self.chunk_start, text_off)
            .raw()
    }
}

/// The `(contig index, 0-based position, reference length)` of a record's
/// shared header (the first 12 bytes), for index building.
///
/// # Errors
///
/// [`IoError::MalformedRecord`] when fewer than 12 bytes are available.
pub fn record_site(bytes: &[u8]) -> Result<(i32, i32, i32), IoError> {
    let head = bytes
        .get(..12)
        .ok_or_else(|| malformed("BCF record header is truncated"))?;
    let int = |i: usize| i32::from_le_bytes([head[i], head[i + 1], head[i + 2], head[i + 3]]);
    Ok((int(0), int(4), int(8)))
}

/// The byte length of a typed field with `count` elements per sample across
/// `n_sample` samples, checked against overflow.
///
/// On 64-bit targets the product cannot overflow (count ≤ i32::MAX, n_sample
/// ≤ 2^24−1, width ≤ 4), but the multiplication must not *depend* on that;
/// an overflow is a malformed record, never a wrap (doc_v10 N8, D10-8).
///
/// # Errors
///
/// [`IoError::MalformedRecord`] when the product overflows.
pub fn field_bytes(ty: u8, count: usize, n_sample: usize) -> Result<usize, IoError> {
    type_width(ty)
        .checked_mul(count)
        .and_then(|v| v.checked_mul(n_sample))
        .ok_or_else(|| malformed("BCF field size overflows"))
}

/// Walks the individual block's FORMAT descriptors, returning
/// `(tag id, type, count, payload offset)` in field order (V22-R11).
///
/// The column decoder and the text renderer share this walk; `bad` lets each
/// caller attach its own locus to the two malformed-record cases (a negative
/// key or a truncated payload).
///
/// # Errors
///
/// [`IoError::MalformedRecord`] for a negative key, a truncated descriptor or
/// an overflowing field size.
pub(crate) fn format_fields(
    indiv: &[u8],
    n_fmt: usize,
    n_sample: usize,
    bad: impl Fn(&'static str) -> IoError,
) -> Result<Vec<(u32, u8, usize, usize)>, IoError> {
    let mut fields = Vec::with_capacity(n_fmt);
    walk_format_fields(indiv, n_fmt, n_sample, bad, &mut fields)?;
    Ok(fields)
}

/// The one FORMAT descriptor walk (doc_v28 V28-21): both the allocating
/// [`format_fields`] and the caller-owned [`FormatTable`] are views over this,
/// so a change to the framing rule cannot drift between them.
fn walk_format_fields(
    indiv: &[u8],
    n_fmt: usize,
    n_sample: usize,
    bad: impl Fn(&'static str) -> IoError,
    out: &mut Vec<(u32, u8, usize, usize)>,
) -> Result<(), IoError> {
    let mut ci = Cursor::new(indiv);
    for _ in 0..n_fmt {
        let key = ci.typed_int()?;
        let key = u32::try_from(key).map_err(|_| bad("BCF FORMAT key is negative"))?;
        let (ty, count) = ci.typed()?;
        let start = indiv.len() - ci.remaining();
        ci.take(field_bytes(ty, count, n_sample)?)?;
        out.push((key, ty, count, start));
    }
    Ok(())
}

/// A reusable FORMAT descriptor table for one record (doc_v28 V28-16).
///
/// The GT decoder used to find its own `GT` entry and then the column
/// decoder walked the same descriptors again, allocating a fresh `Vec` per
/// record. One caller-owned table serves both: the walk happens once, and
/// its capacity survives across records.
#[derive(Default, Debug)]
pub struct FormatTable {
    fields: Vec<(u32, u8, usize, usize)>,
}

impl FormatTable {
    /// A fresh table.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Rebuilds the table for `record`'s individual block.
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] as [`format_fields`].
    pub fn rebuild(
        &mut self,
        record: &BcfRecordRef<'_>,
        n_fmt: usize,
        n_sample: usize,
    ) -> Result<(), IoError> {
        self.fields.clear();
        walk_format_fields(record.indiv, n_fmt, n_sample, malformed, &mut self.fields)
    }

    /// The entries in field order: `(tag id, type, count, payload offset)`.
    #[must_use]
    pub fn fields(&self) -> &[(u32, u8, usize, usize)] {
        &self.fields
    }
}

/// The fixed 24-byte shared-block header (hts-specs §6.3.1), decoded once
/// (doc_v28 V28-21).
///
/// Three callers used to unpack these six words themselves — the GT site
/// reader, the text renderer and the native encoder — with the same shifts and
/// the same zero-allele rule written out three times. This type owns the
/// framing and that one validation; everything semantic (contig lookup, the
/// telomere sentinel, what to render) stays with the callers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SharedHeader {
    /// Contig dictionary index.
    pub chrom: i32,
    /// 0-based position.
    pub pos0: i32,
    /// Reference length in bases.
    pub rlen: i32,
    /// Raw QUAL bits (`0x7F80_0001` is missing).
    pub qual: u32,
    /// Allele count, REF included.
    pub n_allele: usize,
    /// INFO field count.
    pub n_info: usize,
    /// FORMAT field count.
    pub n_fmt: usize,
    /// Samples declared by the record.
    pub n_sample: usize,
}

impl SharedHeader {
    /// Decodes the fixed header and leaves the cursor on the ID string.
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] when the shared block is shorter than its
    /// header, and when the record declares no allele (doc_v26 V26-05).
    pub fn parse(shared: &[u8]) -> Result<(Self, Cursor<'_>), IoError> {
        let mut c = Cursor::new(shared);
        let chrom = c.i32()?;
        let pos0 = c.i32()?;
        let rlen = c.i32()?;
        let qual = c.u32()?;
        let na_ni = c.u32()?;
        let nf_ns = c.u32()?;
        let n_allele = (na_ni >> 16) as usize;
        // doc_v26 V26-05: REF is required. A zero allele count means the
        // shared block has no allele string at all; htslib rejects such a
        // record, so every reed path must too.
        if n_allele == 0 {
            return Err(malformed(
                "BCF record declares zero alleles (n_allele = 0); REF is required",
            ));
        }
        Ok((
            Self {
                chrom,
                pos0,
                rlen,
                qual,
                n_allele,
                n_info: (na_ni & 0xFFFF) as usize,
                n_fmt: (nf_ns >> 24) as usize,
                n_sample: (nf_ns & 0x00FF_FFFF) as usize,
            },
            c,
        ))
    }
}

/// A cursor over one block.
pub struct Cursor<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    /// A cursor at the block's start.
    #[must_use]
    pub const fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }

    /// Bytes not yet consumed.
    #[must_use]
    pub const fn remaining(&self) -> usize {
        self.data.len() - self.pos
    }

    /// Takes `n` bytes.
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] at the block's end.
    pub fn take(&mut self, n: usize) -> Result<&'a [u8], IoError> {
        let end = self
            .pos
            .checked_add(n)
            .ok_or_else(|| malformed("BCF block overrun"))?;
        let slice = self
            .data
            .get(self.pos..end)
            .ok_or_else(|| malformed("BCF block ends mid-field"))?;
        self.pos = end;
        Ok(slice)
    }

    /// One unsigned byte.
    ///
    /// # Errors
    ///
    /// As [`Cursor::take`].
    pub fn u8(&mut self) -> Result<u8, IoError> {
        Ok(self.take(1)?[0])
    }

    /// One little-endian `u32`.
    ///
    /// # Errors
    ///
    /// As [`Cursor::take`].
    pub fn u32(&mut self) -> Result<u32, IoError> {
        let b = self.take(4)?;
        Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    /// One little-endian `i32`.
    ///
    /// # Errors
    ///
    /// As [`Cursor::take`].
    pub fn i32(&mut self) -> Result<i32, IoError> {
        Ok(self.u32()? as i32)
    }

    /// One little-endian `f32`.
    ///
    /// # Errors
    ///
    /// As [`Cursor::take`].
    pub fn f32(&mut self) -> Result<f32, IoError> {
        Ok(f32::from_bits(self.u32()?))
    }

    /// A typed descriptor: `(type, count)`; counts of 15 read the extended
    /// typed integer that follows (§6.3.3).
    ///
    /// Only one extension level is legal. htslib's `bcf_dec_size_safe`
    /// (vcf.c:2013-2026) reads the extension with `bcf_dec_typed_int1_safe`,
    /// which decodes the *type* nibble and a value and never recurses; a
    /// second 15 is rejected there because 15 is not an integer type. Reed
    /// previously routed the extension back through `typed`, so a run of
    /// `0xF1` bytes recursed `typed → typed_int → typed` until the stack
    /// overflowed (doc_v23 V23-01).
    ///
    /// A NULL type with a nonzero count is malformed — NULL carries no
    /// payload, so its count must be zero (an INFO flag). htslib rejects it
    /// in `bcf_record_check` (vcf.c:2231-2236); without the check the empty
    /// payload reaches the renderers, which index past it (doc_v10 D10-2).
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] for an invalid descriptor.
    pub fn typed(&mut self) -> Result<(u8, usize), IoError> {
        let desc = self.u8()?;
        let ty = desc & 0x0F;
        let mut count = usize::from(desc >> 4);
        if count == 15 {
            let ext = self.extended_count()?;
            count = ext;
        }
        if type_width(ty) == 0 && ty != value_type::NULL {
            return Err(malformed("BCF descriptor has an unknown type"));
        }
        if ty == value_type::NULL && count > 0 {
            return Err(malformed("BCF null-typed descriptor has a nonzero count"));
        }
        Ok((ty, count))
    }

    /// The extended count after a `15` nibble: exactly one descriptor byte,
    /// whose own count nibble must not be another extension.
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] for a second-level extension, a
    /// non-integer type, a count other than one or a negative value.
    fn extended_count(&mut self) -> Result<usize, IoError> {
        let desc = self.u8()?;
        let ty = desc & 0x0F;
        let count = desc >> 4;
        if count == 15 {
            return Err(malformed("BCF descriptor has a second-level extension"));
        }
        if count != 1 {
            return Err(malformed("BCF typed integer has a count other than 1"));
        }
        let ext = self.int_value(ty)?;
        usize::try_from(ext).map_err(|_| malformed("BCF descriptor count is negative"))
    }

    /// A typed integer (the key/index encoding).
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] for an invalid descriptor.
    pub fn typed_int(&mut self) -> Result<i64, IoError> {
        let desc = self.u8()?;
        let ty = desc & 0x0F;
        if desc >> 4 == 15 {
            return Err(malformed("BCF descriptor has a second-level extension"));
        }
        if desc >> 4 != 1 {
            return Err(malformed("BCF typed integer has a count other than 1"));
        }
        self.int_value(ty)
    }

    /// One integer element of type `ty` (no sentinel interpretation).
    ///
    /// # Errors
    ///
    /// As [`Cursor::take`].
    pub fn int_value(&mut self, ty: u8) -> Result<i64, IoError> {
        match ty {
            value_type::INT8 => Ok(i64::from(self.take(1)?[0] as i8)),
            value_type::INT16 => {
                let b = self.take(2)?;
                Ok(i64::from(i16::from_le_bytes([b[0], b[1]])))
            }
            value_type::INT32 => Ok(i64::from(self.i32()?)),
            _ => Err(malformed("BCF integer field has a non-integer type")),
        }
    }

    /// A typed string: `(bytes)`; the bytes are NUL-padded but the stored
    /// length is the descriptor count.
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] when the descriptor is not CHAR.
    pub fn typed_string(&mut self) -> Result<&'a [u8], IoError> {
        let (ty, count) = self.typed()?;
        if ty != value_type::CHAR {
            return Err(malformed("BCF string field is not CHAR-typed"));
        }
        self.take(count)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn record_length_and_split() {
        // l_shared = 4, l_indiv = 2: one 14-byte record.
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&4_u32.to_le_bytes());
        bytes.extend_from_slice(&2_u32.to_le_bytes());
        bytes.extend_from_slice(b"abcd");
        bytes.extend_from_slice(b"ef");
        assert_eq!(record_len(&bytes), Some(14));
        let r = record_ref(&bytes).expect("invariant: complete");
        assert_eq!(r.shared, b"abcd");
        assert_eq!(r.indiv, b"ef");
        // Trailing data (the next record) does not leak in.
        bytes.push(0);
        assert_eq!(record_len(&bytes), Some(14));
        assert_eq!(
            record_ref(&buf(&bytes)).expect("invariant: parses").indiv,
            b"ef"
        );
    }

    fn buf(b: &[u8]) -> Vec<u8> {
        b.to_vec()
    }

    #[test]
    fn typed_descriptors_read_the_extended_count() {
        // CHAR with count 15 -> extended typed int 20.
        let data = [0xF7, 0x11, 20];
        let mut c = Cursor::new(&data);
        assert_eq!(
            c.typed().expect("invariant: parses"),
            (value_type::CHAR, 20)
        );
    }

    #[test]
    fn second_level_extension_is_malformed_not_a_stack_overflow() {
        // A run of 0xF1 (count nibble 15, INT8) used to recurse
        // typed -> typed_int -> typed until the stack overflowed (V23-01).
        for repeat in [2_usize, 8, 64, 4096] {
            let data = vec![0xF1_u8; repeat];
            let mut c = Cursor::new(&data);
            let err = c.typed().expect_err("invariant: second extension rejected");
            assert!(err.to_string().contains("second-level extension"), "{err}");
            // The key path has the same guard.
            let mut c = Cursor::new(&data);
            assert!(c.typed_int().is_err());
        }
    }

    #[test]
    fn extended_count_rejects_non_integer_and_non_unit() {
        // Extension must be one integer of count 1: a CHAR type right after
        // the 15 nibble is malformed, as is a count nibble other than 1.
        for desc in [0x17_u8, 0x27, 0x25, 0x65] {
            let data = [0xF7, desc, 0, 0, 0, 0];
            let mut c = Cursor::new(&data);
            assert!(
                c.typed().is_err(),
                "descriptor {desc:#04x} must be rejected"
            );
        }
    }

    #[test]
    fn typed_int_reads_a_single_descriptor_without_extending() {
        // htslib bcf_enc_int1: one count-1 descriptor byte then the value.
        let data = [0x11, 0x07];
        let mut c = Cursor::new(&data);
        assert_eq!(c.typed_int().expect("invariant: parses"), 7);
        // Count 2 is not a typed integer.
        let mut c = Cursor::new(&[0x21_u8, 0x07, 0x08]);
        assert!(c.typed_int().is_err());
        // A 15 nibble is a size extension, never a typed integer.
        let mut c = Cursor::new(&[0xF1_u8, 0x11, 0x07]);
        assert!(c.typed_int().is_err());
    }

    #[test]
    fn truncated_record_is_an_error() {
        let bytes = [4_u8, 0, 0, 0, 2, 0, 0, 0, b'a'];
        assert_eq!(record_len(&bytes), None);
        assert!(record_ref(&bytes).is_err());
    }
}
