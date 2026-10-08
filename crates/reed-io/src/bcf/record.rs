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
        locus: Locus { contig_name: "<bcf>".into(), pos1: 0 },
        reason,
    }
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
    let total = record_len(bytes).ok_or_else(|| malformed("BCF record is truncated"))?;
    let head = &bytes[..8];
    let l_shared = u32::from_le_bytes([head[0], head[1], head[2], head[3]]) as usize;
    let l_indiv = u32::from_le_bytes([head[4], head[5], head[6], head[7]]) as usize;
    if total > bytes.len() {
        return Err(malformed("BCF record is truncated"));
    }
    Ok(BcfRecordRef {
        shared: &bytes[8..8 + l_shared],
        indiv: &bytes[8 + l_shared..8 + l_shared + l_indiv],
    })
}

/// Maps an inflated-text offset within one chunk back to a virtual offset
/// through the chunk's block table (BGZF's `(compressed << 16) | offset`).
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
        let mut remaining = text_off as u64;
        let mut coff = self.chunk_start;
        for &(clen, ilen) in self.blocks {
            if remaining < u64::from(ilen) {
                return (coff << 16) | remaining;
            }
            remaining -= u64::from(ilen);
            coff += u64::from(clen);
        }
        coff << 16
    }
}

/// The `(contig index, 0-based position, reference length)` of a record's
/// shared header (the first 12 bytes), for index building.
///
/// # Errors
///
/// [`IoError::MalformedRecord`] when fewer than 12 bytes are available.
pub fn record_site(bytes: &[u8]) -> Result<(i32, i32, i32), IoError> {
    let head = bytes.get(..12).ok_or_else(|| malformed("BCF record header is truncated"))?;
    let int = |i: usize| i32::from_le_bytes([head[i], head[i + 1], head[i + 2], head[i + 3]]);
    Ok((int(0), int(4), int(8)))
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
        let end = self.pos.checked_add(n).ok_or_else(|| malformed("BCF block overrun"))?;
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
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] for an invalid descriptor.
    pub fn typed(&mut self) -> Result<(u8, usize), IoError> {
        let desc = self.u8()?;
        let ty = desc & 0x0F;
        let mut count = usize::from(desc >> 4);
        if count == 15 {
            let ext = self.typed_int()?;
            count = usize::try_from(ext)
                .map_err(|_| malformed("BCF descriptor count is negative"))?;
        }
        if type_width(ty) == 0 && ty != value_type::NULL {
            return Err(malformed("BCF descriptor has an unknown type"));
        }
        Ok((ty, count))
    }

    /// A typed integer (the key/index encoding).
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] for an invalid descriptor.
    pub fn typed_int(&mut self) -> Result<i64, IoError> {
        let (ty, count) = self.typed()?;
        if count != 1 {
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

    /// Skips a typed value of `(ty, count)`.
    ///
    /// # Errors
    ///
    /// As [`Cursor::take`].
    pub fn skip(&mut self, ty: u8, count: usize) -> Result<(), IoError> {
        self.take(type_width(ty) * count)?;
        Ok(())
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
        assert_eq!(record_ref(&buf(&bytes)).expect("invariant: parses").indiv, b"ef");
    }

    fn buf(b: &[u8]) -> Vec<u8> {
        b.to_vec()
    }

    #[test]
    fn typed_descriptors_read_the_extended_count() {
        // CHAR with count 15 -> extended typed int 20.
        let data = [0xF7, 0x11, 20];
        let mut c = Cursor::new(&data);
        assert_eq!(c.typed().expect("invariant: parses"), (value_type::CHAR, 20));
    }

    #[test]
    fn truncated_record_is_an_error() {
        let bytes = [4_u8, 0, 0, 0, 2, 0, 0, 0, b'a'];
        assert_eq!(record_len(&bytes), None);
        assert!(record_ref(&bytes).is_err());
    }
}
