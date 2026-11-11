//! Native FORMAT column extraction (doc_v4 O-4, doc_v22 R22-2).
//!
//! The filter's mask stage reads a few typed columns without rendering the
//! record to text. This is the BCF analogue of the text scanner's
//! `GenotypeScanner`; `tests/format_columns_vs_bcf.rs` pins both
//! implementations to the same values (doc_v22 R22-4).
//!
//! The descriptor walk collects the field table because an `IntSpec` may
//! reference any key (and a fallback); the renderer collects the same table.
//! Both use [`FormatTable`], so the field-major layout rule lives in one
//! place (V22-R11, doc_v28 V28-16).

// Casts in this file are bounded by construction: Element widths and counts come from the record descriptor.
#![allow(clippy::cast_possible_truncation)]

use super::dict::BcfDict;
use super::gt::SiteHeader;
use super::record::{
    BcfRecordRef, FormatTable, int_element, int_sentinels, type_width, value_type,
};
use crate::error::IoError;
use crate::vcf::format_layout::IntSpec;
use crate::vcf::genotype_scan::{INT_MISSING, clamp_int_sum};

/// Decodes the requested FORMAT columns for every sample without rendering
/// them (doc_v4 O-4): the typed values feed the filter's mask stage directly.
///
/// Each [`IntSpec`] names its primary key and optional fallback; the policy
/// is explicit at the call site (doc_v22 D22-1). The primary key reads one
/// integer per sample; the fallback sums every value of the sample's vector
/// (missing values skipped, vector-end truncates), the text scanner's
/// [`crate::vcf::genotype_scan::parse_int_sum`] rule.
///
/// Missing values use the text scanner's sentinels — [`INT_MISSING`] and
/// `NaN` — and a negative or non-integer value maps to `INT_MISSING`,
/// exactly as the text scanner's `parse_int` does, so the native and text
/// paths agree bit for bit.
///
/// This builds its own descriptor table; a caller that already walked the
/// record's FORMAT (the GT decoder does) passes that table to
/// [`decode_format_columns_with_table`] and skips the second walk (doc_v28
/// V28-16).
///
/// # Errors
///
/// [`IoError::MalformedRecord`] for a truncated individual block.
pub fn decode_format_columns(
    rec: &BcfRecordRef<'_>,
    header: &SiteHeader,
    dict: &BcfDict,
    int_specs: &[IntSpec],
    float_keys: &[&[u8]],
    int_columns: &mut [Vec<u32>],
    float_columns: &mut [Vec<f32>],
) -> Result<(), IoError> {
    let mut table = FormatTable::new();
    if header.n_fmt > 0 {
        table.rebuild(rec, header.n_fmt, header.n_sample)?;
    }
    decode_format_columns_with_table(
        rec,
        header,
        dict,
        &table,
        int_specs,
        float_keys,
        int_columns,
        float_columns,
    )
}

/// [`decode_format_columns`] over a table the caller already built for
/// `rec` (doc_v28 V28-16). The GT decoder builds it once per record, so the
/// descriptor walk and its allocation happen once instead of twice.
///
/// # Errors
///
/// [`IoError::MalformedRecord`] for a truncated individual block.
#[allow(clippy::too_many_arguments)]
pub fn decode_format_columns_with_table(
    rec: &BcfRecordRef<'_>,
    header: &SiteHeader,
    dict: &BcfDict,
    table: &FormatTable,
    int_specs: &[IntSpec],
    float_keys: &[&[u8]],
    int_columns: &mut [Vec<u32>],
    float_columns: &mut [Vec<f32>],
) -> Result<(), IoError> {
    let n_sample = header.n_sample;
    for column in int_columns.iter_mut() {
        column.clear();
        column.resize(n_sample, INT_MISSING);
    }
    for column in float_columns.iter_mut() {
        column.clear();
        column.resize(n_sample, f32::NAN);
    }
    if header.n_fmt == 0 || n_sample == 0 || (int_specs.is_empty() && float_keys.is_empty()) {
        return Ok(());
    }
    let indiv = rec.indiv;
    let fields = table.fields();
    // Integer columns: each spec resolves its primary field (or the fallback
    // when this record's FORMAT lacks the primary) once, then the column is
    // filled sample by sample.
    for (spec, column) in int_specs.iter().zip(int_columns.iter_mut()) {
        let find = |key: &[u8]| {
            fields
                .iter()
                .find(|(id, ..)| dict.tag(*id) == Some(key))
                .copied()
        };
        let (ty, count, start, sum) = if let Some((_, ty, count, start)) = find(&spec.key) {
            (ty, count, start, false)
        } else if let Some((_, ty, count, start)) = spec.fallback.as_deref().and_then(find) {
            (ty, count, start, true)
        } else {
            continue;
        };
        let width = type_width(ty);
        let is_int = matches!(ty, value_type::INT8 | value_type::INT16 | value_type::INT32);
        if width == 0 || !is_int || (!sum && count != 1) {
            continue;
        }
        for (s, slot) in column.iter_mut().enumerate() {
            let data = &indiv[start + s * count * width..start + (s + 1) * count * width];
            *slot = if sum {
                int_sum_column(ty, count, data)
            } else {
                int_element_column(ty, data)
            };
        }
    }
    for &(key, ty, count, start) in fields {
        let Some(name) = dict.tag(key) else { continue };
        let width = type_width(ty);
        if width == 0 || count != 1 {
            continue;
        }
        if let Some(i) = float_keys.iter().position(|k| *k == name) {
            let Some(column) = float_columns.get_mut(i) else {
                continue;
            };
            for (s, slot) in column.iter_mut().enumerate() {
                let elem = &indiv[start + s * width..start + (s + 1) * width];
                *slot = float_element_column(ty, elem);
            }
        }
    }
    Ok(())
}

/// One sample's integer vector summed the way the text scanner's
/// `parse_int_sum` reads it: missing elements are skipped, vector-end
/// truncates, and a sample with no value at all is [`INT_MISSING`].
///
/// The accumulator is `u64` and clamped like the text path, so a vector of
/// large elements cannot overflow into [`INT_MISSING`] (V22-R2-2).
fn int_sum_column(ty: u8, count: usize, data: &[u8]) -> u32 {
    let width = type_width(ty);
    let (missing, vector_end) = int_sentinels(ty);
    let mut sum = 0_u64;
    let mut any = false;
    for j in 0..count {
        let v = int_element(ty, &data[j * width..(j + 1) * width]);
        if v == vector_end {
            break;
        }
        if v == missing || v < 0 {
            continue;
        }
        sum = sum.saturating_add(u64::try_from(v).unwrap_or(u64::MAX));
        any = true;
    }
    if any { clamp_int_sum(sum) } else { INT_MISSING }
}

/// One integer FORMAT element as the text scanner would parse it.
fn int_element_column(ty: u8, elem: &[u8]) -> u32 {
    let (missing, _) = int_sentinels(ty);
    let value = int_element(ty, elem);
    if value == missing || value < 0 {
        INT_MISSING
    } else {
        u32::try_from(value).unwrap_or(INT_MISSING)
    }
}

/// One float FORMAT element as the text scanner would parse it.
///
/// A descriptor that is neither FLOAT nor an integer type cannot be read as
/// a number: a CHAR element is one byte wide, so `int_element` would index
/// past the slice (doc_v23 V23-02). Such an element is missing.
fn float_element_column(ty: u8, elem: &[u8]) -> f32 {
    use super::record::{FLOAT_MISSING, FLOAT_VECTOR_END};
    if ty == value_type::FLOAT {
        let bits = u32::from_le_bytes([elem[0], elem[1], elem[2], elem[3]]);
        if bits == FLOAT_MISSING || bits == FLOAT_VECTOR_END {
            f32::NAN
        } else {
            f32::from_bits(bits)
        }
    } else if matches!(ty, value_type::INT8 | value_type::INT16 | value_type::INT32) {
        let (missing, vector_end) = int_sentinels(ty);
        let value = int_element(ty, elem);
        if value == missing || value == vector_end {
            f32::NAN
        } else {
            value as f32
        }
    } else {
        f32::NAN
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// doc_v23 V23-02: a CHAR-typed element is one byte wide; reading it as
    /// int32 must not index past the slice. It decodes as missing.
    #[test]
    fn char_typed_float_column_element_is_missing_not_a_panic() {
        // One-byte slice: the old code read elem[1..4] and panicked.
        assert!(float_element_column(value_type::CHAR, b"0").is_nan());
        assert!(float_element_column(value_type::NULL, b"").is_nan());
        // Integer and FLOAT types keep their values.
        assert_eq!(float_element_column(value_type::INT8, &[7]), 7.0);
        assert!(float_element_column(value_type::INT8, &[i8::MIN as u8]).is_nan());
        // doc_v28 V28-05: the integer vector-end sentinel is missing too;
        // leaking it made an integer-backed float column read -127.0.
        assert!(float_element_column(value_type::INT8, &[(i8::MIN + 1) as u8]).is_nan());
        assert!(float_element_column(value_type::INT16, &(i16::MIN + 1).to_le_bytes()).is_nan());
        assert!(float_element_column(value_type::INT32, &(i32::MIN + 1).to_le_bytes()).is_nan());
        assert_eq!(
            float_element_column(value_type::FLOAT, &1.5_f32.to_le_bytes()),
            1.5
        );
    }

    /// The integer column helper's sentinel rules stay unchanged: a
    /// negative value is missing, a positive one is kept. The CHAR case is
    /// gated in `decode_format_columns` (`is_int`), so it never reaches
    /// `int_element` through the public path.
    #[test]
    fn integer_column_element_keeps_its_sentinel_rules() {
        assert_eq!(int_element_column(value_type::INT8, &[5]), 5);
        // A CHAR path never reaches `int_element_column` from
        // `decode_format_columns` (its `is_int` gate), but the helper stays
        // total for direct callers.
        assert_eq!(
            int_element_column(value_type::INT32, &(-3_i32).to_le_bytes()),
            INT_MISSING
        );
    }
}
