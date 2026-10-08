//! BCF → VCF text rendering, byte-compatible with `bcftools view`
//! (doc_v3 02-bcf.md B3a).
//!
//! Rules come from htslib `vcf.c`: `bcf_fmt_array` / `bcf_fmt_array1` for
//! numeric and string values, `bcf_format_gt_v2` for genotypes (the
//! separator before an allele uses *that* allele's phased bit), `kputd` for
//! floats (6 significant digits, `%g` beyond `[1e-4, 1e6)`), missing →
//! `.`, vector-end → truncation, `;` between INFO fields, `:` between
//! FORMAT fields and `\t` between samples.

// Casts in this file are bounded by construction: Element widths and counts come from the record descriptor.
#![allow(clippy::cast_possible_truncation)]


use super::dict::BcfDict;
use super::record::{
    int_sentinels, malformed, type_width, value_type, BcfRecordRef, Cursor, FLOAT_MISSING,
    FLOAT_VECTOR_END,
};
use crate::error::{IoError, Locus};
use crate::vcf::line_edit::{format_g6_into};

/// htslib's `bcf_str_missing`: a single char value equal to this prints `.`
/// in the `n == 1` FORMAT path.
const STR_MISSING: u8 = 0x07;

/// Renders one BCF record as a VCF line (newline included).
///
/// # Errors
///
/// [`IoError::MalformedRecord`] for a truncated or invalid record, and
/// [`IoError::Header`]-level dictionary misses (reported as malformed, since
/// they can only come from a header/record mismatch).
pub fn render_record(
    rec: &BcfRecordRef<'_>,
    dict: &BcfDict,
    out: &mut Vec<u8>,
) -> Result<(), IoError> {
    let mut c = Cursor::new(rec.shared);
    let chrom = c.i32()?;
    let pos0 = c.i32()?;
    let _rlen = c.i32()?;
    let qual_bits = c.u32()?;
    let na_ni = c.u32()?;
    let nf_ns = c.u32()?;
    let n_allele = (na_ni >> 16) as usize;
    let n_info = (na_ni & 0xFFFF) as usize;
    let n_fmt = (nf_ns >> 24) as usize;
    let n_sample = (nf_ns & 0x00FF_FFFF) as usize;

    let locus = Locus {
        contig_name: String::from_utf8_lossy(
            dict.contig(chrom)
                .ok_or_else(|| malformed("BCF contig index is not in the header dictionary"))?,
        )
        .into_owned(),
        pos1: u64::try_from(pos0).unwrap_or(0) + 1,
    };
    out.extend_from_slice(locus.contig_name.as_bytes());
    out.push(b'\t');
    push_u64(out, u64::try_from(pos0).unwrap_or(0) + 1);

    // ID (empty prints `.`).
    let id = c.typed_string()?;
    out.push(b'\t');
    if id.is_empty() {
        out.push(b'.');
    } else {
        out.extend_from_slice(id);
    }

    // Alleles: REF then ALT (a lone allele prints `.`).
    out.push(b'\t');
    if n_allele == 0 {
        out.push(b'.');
    } else {
        out.extend_from_slice(c.typed_string()?);
    }
    out.push(b'\t');
    if n_allele <= 1 {
        out.push(b'.');
    } else {
        for i in 1..n_allele {
            if i > 1 {
                out.push(b',');
            }
            out.extend_from_slice(c.typed_string()?);
        }
    }

    // QUAL.
    out.push(b'\t');
    if qual_bits == FLOAT_MISSING || qual_bits == FLOAT_VECTOR_END {
        out.push(b'.');
    } else {
        float_into(f32::from_bits(qual_bits), out);
    }

    // FILTER (dictionary indices; empty prints `.`).
    out.push(b'\t');
    let (f_ty, f_count) = c.typed()?;
    if f_count == 0 {
        out.push(b'.');
    } else {
        for i in 0..f_count {
            if i > 0 {
                out.push(b';');
            }
            let idx = c.int_value(f_ty)?;
            let idx = u32::try_from(idx).map_err(|_| bad(&locus, "BCF FILTER index is negative"))?;
            out.extend_from_slice(
                dict.tag(idx)
                    .ok_or_else(|| bad(&locus, "BCF FILTER index is not in the header dictionary"))?,
            );
        }
    }

    // INFO.
    out.push(b'\t');
    if n_info == 0 {
        out.push(b'.');
    } else {
        for i in 0..n_info {
            if i > 0 {
                out.push(b';');
            }
            let key = c.typed_int()?;
            let key = u32::try_from(key).map_err(|_| bad(&locus, "BCF INFO key is negative"))?;
            out.extend_from_slice(
                dict.tag(key)
                    .ok_or_else(|| bad(&locus, "BCF INFO key is not in the header dictionary"))?,
            );
            let (ty, count) = c.typed()?;
            if count == 0 {
                continue; // a flag: key only
            }
            out.push(b'=');
            let data = c.take(type_width(ty) * count)?;
            fmt_array(out, ty, count, data);
        }
    }

    // FORMAT + samples.
    out.push(b'\t');
    if n_fmt == 0 {
        for _ in 0..=n_sample {
            out.extend_from_slice(b"\t.");
        }
    } else {
        // Field descriptors live in the individual block; values are
        // field-major: each field's descriptor is followed by its
        // `n_sample × count` elements (htslib's packed layout).
        let indiv = rec.indiv;
        let mut ci = Cursor::new(indiv);
        let mut fields: Vec<(u32, u8, usize, usize)> = Vec::with_capacity(n_fmt);
        for _ in 0..n_fmt {
            let key = ci.typed_int()?;
            let key = u32::try_from(key).map_err(|_| bad(&locus, "BCF FORMAT key is negative"))?;
            let (ty, count) = ci.typed()?;
            let start = indiv.len() - ci.remaining();
            ci.take(type_width(ty) * count * n_sample)?;
            fields.push((key, ty, count, start));
        }

        let mut first = true;
        for &(key, ..) in &fields {
            if !first {
                out.push(b':');
            }
            first = false;
            out.extend_from_slice(
                dict.tag(key)
                    .ok_or_else(|| bad(&locus, "BCF FORMAT key is not in the header dictionary"))?,
            );
        }
        // Hoist everything the per-sample loop would otherwise recompute.
        let mut layout: Vec<(bool, u8, usize, usize, usize)> = Vec::with_capacity(fields.len());
        for &(key, ty, count, field_start) in &fields {
            let is_gt = dict.tag(key) == Some(b"GT".as_slice());
            layout.push((is_gt, ty, count, field_start, type_width(ty) * count));
        }
        for s in 0..n_sample {
            out.push(b'\t');
            for (i, &(is_gt, ty, count, base, stride)) in layout.iter().enumerate() {
                if i > 0 {
                    out.push(b':');
                }
                let start = base + s * stride;
                let data = &indiv[start..start + stride];
                if is_gt {
                    fmt_gt(out, ty, count, data);
                } else if count == 1 {
                    fmt_array1(out, ty, data);
                } else {
                    fmt_array(out, ty, count, data);
                }
            }
        }
    }
    out.push(b'\n');
    Ok(())
}

/// A malformed record at a known locus.
fn bad(locus: &Locus, reason: &'static str) -> IoError {
    IoError::MalformedRecord { locus: locus.clone(), reason }
}

/// `bcf_fmt_array`: elements joined by `,`; vector-end truncates the vector;
/// missing prints `.`; char vectors print up to the first NUL (empty → `.`).
fn fmt_array(out: &mut Vec<u8>, ty: u8, count: usize, data: &[u8]) {
    if count == 0 {
        out.push(b'.');
        return;
    }
    if ty == value_type::CHAR {
        let end = data.iter().position(|&b| b == 0).unwrap_or(data.len());
        out.extend_from_slice(&data[..end]);
        return;
    }
    let width = type_width(ty);
    for j in 0..count {
        let elem = &data[j * width..(j + 1) * width];
        match ty {
            value_type::FLOAT => {
                let bits = u32::from_le_bytes([elem[0], elem[1], elem[2], elem[3]]);
                if bits == FLOAT_VECTOR_END {
                    break;
                }
                if j > 0 {
                    out.push(b',');
                }
                if bits == FLOAT_MISSING {
                    out.push(b'.');
                } else {
                    float_into(f32::from_bits(bits), out);
                }
            }
            _ => {
                let v = int_element(ty, elem);
                let (missing, vector_end) = int_sentinels(ty);
                if v == vector_end {
                    break;
                }
                if j > 0 {
                    out.push(b',');
                }
                if v == missing {
                    out.push(b'.');
                } else {
                    push_i64(out, v);
                }
            }
        }
    }
}

/// `bcf_fmt_array1`: the `n == 1` fast path — a char equal to
/// `bcf_str_missing` prints `.`, a vector-end prints nothing.
fn fmt_array1(out: &mut Vec<u8>, ty: u8, data: &[u8]) {
    match ty {
        value_type::CHAR => {
            let b = data[0];
            out.push(if b == STR_MISSING { b'.' } else { b });
        }
        value_type::FLOAT => {
            let bits = u32::from_le_bytes([data[0], data[1], data[2], data[3]]);
            if bits == FLOAT_VECTOR_END {
                return;
            }
            if bits == FLOAT_MISSING {
                out.push(b'.');
            } else {
                float_into(f32::from_bits(bits), out);
            }
        }
        _ => {
            let v = int_element(ty, data);
            let (missing, vector_end) = int_sentinels(ty);
            if v == vector_end {
                return;
            }
            if v == missing {
                out.push(b'.');
            } else {
                push_i64(out, v);
            }
        }
    }
}

/// `bcf_format_gt_v2`: `(allele+1)<<1 | phased`; the separator before an
/// allele uses that allele's phased bit; a zero allele prints `.`; an empty
/// genotype prints `.`.
fn fmt_gt(out: &mut Vec<u8>, ty: u8, count: usize, data: &[u8]) {
    let width = type_width(ty);
    let (_, vector_end) = int_sentinels(ty);
    let mut printed = false;
    for i in 0..count {
        let v = int_element(ty, &data[i * width..(i + 1) * width]);
        if v == vector_end {
            break;
        }
        if i > 0 {
            out.push(if v & 1 != 0 { b'|' } else { b'/' });
        }
        let allele = v >> 1;
        if allele == 0 {
            out.push(b'.');
        } else {
            push_u64(out, u64::try_from(allele - 1).unwrap_or(0));
        }
        printed = true;
    }
    if !printed {
        out.push(b'.');
    }
}

/// One integer element of `ty`.
fn int_element(ty: u8, elem: &[u8]) -> i64 {
    match ty {
        value_type::INT8 => i64::from(elem[0] as i8),
        value_type::INT16 => i64::from(i16::from_le_bytes([elem[0], elem[1]])),
        _ => i64::from(i32::from_le_bytes([elem[0], elem[1], elem[2], elem[3]])),
    }
}

/// htslib `kputd`: 6 significant digits with `%g` beyond `[1e-4, 1e6)`;
/// `-0` keeps its sign.
fn float_into(v: f32, out: &mut Vec<u8>) {
    if v == 0.0 && v.is_sign_negative() {
        out.extend_from_slice(b"-0");
        return;
    }
    format_g6_into(out, f64::from(v));
}

/// Decimal rendering without the `core::fmt` machinery (this is the hot
/// path: ~3500 numbers per record on a 348-sample file).
fn push_u64(out: &mut Vec<u8>, mut v: u64) {
    if v == 0 {
        out.push(b'0');
        return;
    }
    let mut buf = [0_u8; 20];
    let mut i = buf.len();
    while v > 0 {
        i -= 1;
        buf[i] = b'0' + (v % 10) as u8;
        v /= 10;
    }
    out.extend_from_slice(&buf[i..]);
}

/// Signed decimal, see [`push_u64`].
fn push_i64(out: &mut Vec<u8>, v: i64) {
    if v < 0 {
        out.push(b'-');
        push_u64(out, v.unsigned_abs());
    } else {
        push_u64(out, v as u64);
    }
}
