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
    BcfRecordRef, Cursor, FLOAT_MISSING, FLOAT_VECTOR_END, field_bytes, int_sentinels, malformed,
    type_width, value_type,
};
use crate::error::IoError;
use crate::vcf::genotype_scan::LocusSource;
use crate::vcf::line_edit::format_g6_into;

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
    let (n_sample, n_fmt) = render_site_part(rec, dict, out)?;
    render_sample_part(rec, dict, n_sample, n_fmt, out)?;
    out.push(b'\n');
    Ok(())
}

/// Renders a record's site columns — CHROM, POS, ID, REF, ALT, QUAL, FILTER
/// and the whole INFO column — into `out`, returning the record's
/// `(n_sample, n_fmt)`.
///
/// The native filter path uses this for the line it evaluates head rules and
/// expressions on (doc_v4 O-4): the sample columns are decoded as typed
/// values instead of being rendered.
///
/// # Errors
///
/// As [`render_record`].
pub fn render_site_part(
    rec: &BcfRecordRef<'_>,
    dict: &BcfDict,
    out: &mut Vec<u8>,
) -> Result<(usize, usize), IoError> {
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

    let contig_name = dict
        .contig(chrom)
        .ok_or_else(|| malformed("BCF contig index is not in the header dictionary"))?;
    // (doc_v20 S43, corrected in review) pos0 == -1 is htslib's telomere
    // sentinel and renders as POS 0 (vcf.c:1649-1650); below -1 is malformed.
    // The first version turned both into an error, rejecting BCFs bcftools
    // writes.
    let pos1 = if pos0 == -1 {
        0
    } else {
        u64::try_from(pos0)
            .map_err(|_| bad((contig_name, 0), "BCF position is below the telomere sentinel (-1)"))?
            + 1
    };
    // The error locus stays borrowed (A6-11): the tuple's `LocusSource` impl
    // builds a `Locus` only when an error actually carries it.
    let locus = (contig_name, pos1);
    out.extend_from_slice(contig_name);
    out.push(b'\t');
    push_u64(out, pos1);

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
            let idx = u32::try_from(idx).map_err(|_| bad(locus, "BCF FILTER index is negative"))?;
            out.extend_from_slice(
                dict.tag(idx).ok_or_else(|| {
                    bad(locus, "BCF FILTER index is not in the header dictionary")
                })?,
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
            let key = u32::try_from(key).map_err(|_| bad(locus, "BCF INFO key is negative"))?;
            out.extend_from_slice(
                dict.tag(key)
                    .ok_or_else(|| bad(locus, "BCF INFO key is not in the header dictionary"))?,
            );
            let (ty, count) = c.typed()?;
            if count == 0 {
                continue; // a flag: key only
            }
            out.push(b'=');
            let data = c.take(field_bytes(ty, count, 1)?)?;
            fmt_array(out, ty, count, data);
        }
    }

    Ok((n_sample, n_fmt))
}

/// Renders the FORMAT keys and every sample's values (the `render_record`
/// tail).
///
/// # Errors
///
/// As [`render_record`].
pub fn render_sample_part(
    rec: &BcfRecordRef<'_>,
    dict: &BcfDict,
    n_sample: usize,
    n_fmt: usize,
    out: &mut Vec<u8>,
) -> Result<(), IoError> {
    // The locus is only needed for error messages; keep it borrowed (A6-11).
    let mut hc = Cursor::new(rec.shared);
    let chrom = hc.i32().unwrap_or(-1);
    let pos0 = hc.i32().unwrap_or(-1);
    let locus = (
        dict.contig(chrom).unwrap_or(b"<bcf>"),
        u64::try_from(pos0).unwrap_or(0) + 1,
    );
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
            let key = u32::try_from(key).map_err(|_| bad(locus, "BCF FORMAT key is negative"))?;
            let (ty, count) = ci.typed()?;
            let start = indiv.len() - ci.remaining();
            ci.take(field_bytes(ty, count, n_sample)?)?;
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
                    .ok_or_else(|| bad(locus, "BCF FORMAT key is not in the header dictionary"))?,
            );
        }
        // Hoist everything the per-sample loop would otherwise recompute.
        let mut layout: Vec<(bool, u8, usize, usize, usize)> = Vec::with_capacity(fields.len());
        for &(key, ty, count, field_start) in &fields {
            let is_gt = dict.tag(key) == Some(b"GT".as_slice());
            layout.push((is_gt, ty, count, field_start, field_bytes(ty, count, 1)?));
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
                    fmt_gt(out, ty, count, data, dict.vcf44());
                } else if count == 1 {
                    fmt_array1(out, ty, data);
                } else {
                    fmt_array(out, ty, count, data);
                }
            }
        }
    }
    Ok(())
}

/// Renders the line the native filter path decides on (doc_v4 O-4): the
/// site columns plus a placeholder FORMAT column and one `.` per sample, so
/// the sample *count* stays visible to expressions (`N_SAMPLES`) while the
/// values come from [`decode_format_columns`]. No trailing newline.
///
/// # Errors
///
/// As [`render_record`].
pub fn render_site_line(
    rec: &BcfRecordRef<'_>,
    dict: &BcfDict,
    out: &mut Vec<u8>,
) -> Result<(), IoError> {
    let (n_sample, _) = render_site_part(rec, dict, out)?;
    out.extend_from_slice(b"\t.");
    for _ in 0..n_sample {
        out.extend_from_slice(b"\t.");
    }
    Ok(())
}

/// Decodes the requested `Number=1` FORMAT columns for every sample without
/// rendering them (doc_v4 O-4): the typed values feed the filter's mask
/// stage directly.
///
/// Missing values use the text scanner's sentinels — `INT_MISSING`
/// ([`crate::vcf::genotype_scan::INT_MISSING`]) and `NaN` — and a negative or
/// non-integer value maps to `INT_MISSING`, exactly as the text scanner's
/// `parse_int` does, so the native and text paths agree bit for bit.
///
/// # Errors
///
/// [`IoError::MalformedRecord`] for a truncated individual block.
pub fn decode_format_columns(
    rec: &BcfRecordRef<'_>,
    header: &super::gt::SiteHeader,
    dict: &BcfDict,
    int_keys: &[&[u8]],
    float_keys: &[&[u8]],
    int_columns: &mut [Vec<u32>],
    float_columns: &mut [Vec<f32>],
) -> Result<(), IoError> {
    let n_sample = header.n_sample;
    for column in int_columns.iter_mut() {
        column.clear();
        column.resize(n_sample, crate::vcf::genotype_scan::INT_MISSING);
    }
    for column in float_columns.iter_mut() {
        column.clear();
        column.resize(n_sample, f32::NAN);
    }
    if header.n_fmt == 0 || n_sample == 0 || (int_keys.is_empty() && float_keys.is_empty()) {
        return Ok(());
    }
    // The individual block is field-major: descriptors first, then each
    // field's `n_sample × count` elements.
    let indiv = rec.indiv;
    let mut ci = Cursor::new(indiv);
    let mut fields: Vec<(u32, u8, usize, usize)> = Vec::with_capacity(header.n_fmt);
    for _ in 0..header.n_fmt {
        let key = ci.typed_int()?;
        let key = u32::try_from(key).map_err(|_| malformed("BCF FORMAT key is negative"))?;
        let (ty, count) = ci.typed()?;
        let start = indiv.len() - ci.remaining();
        ci.take(field_bytes(ty, count, n_sample)?)?;
        fields.push((key, ty, count, start));
    }
    for &(key, ty, count, start) in &fields {
        let Some(name) = dict.tag(key) else { continue };
        let width = type_width(ty);
        if width == 0 || count != 1 {
            continue;
        }
        if let Some(i) = int_keys.iter().position(|k| *k == name) {
            let Some(column) = int_columns.get_mut(i) else {
                continue;
            };
            for (s, slot) in column.iter_mut().enumerate() {
                let elem = &indiv[start + s * width..start + (s + 1) * width];
                *slot = int_element_column(ty, elem);
            }
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

/// One integer FORMAT element as the text scanner would parse it.
fn int_element_column(ty: u8, elem: &[u8]) -> u32 {
    let (missing, _) = int_sentinels(ty);
    let value = int_element(ty, elem);
    if value == missing || value < 0 {
        crate::vcf::genotype_scan::INT_MISSING
    } else {
        u32::try_from(value).unwrap_or(crate::vcf::genotype_scan::INT_MISSING)
    }
}

/// One float FORMAT element as the text scanner would parse it.
fn float_element_column(ty: u8, elem: &[u8]) -> f32 {
    if ty == value_type::FLOAT {
        let bits = u32::from_le_bytes([elem[0], elem[1], elem[2], elem[3]]);
        if bits == FLOAT_MISSING || bits == FLOAT_VECTOR_END {
            f32::NAN
        } else {
            f32::from_bits(bits)
        }
    } else {
        let (missing, _) = int_sentinels(ty);
        let value = int_element(ty, elem);
        if value == missing {
            f32::NAN
        } else {
            value as f32
        }
    }
}

/// A malformed record at a known locus (borrowed: the `Locus` is built only
/// when the error actually carries it).
fn bad(locus: impl LocusSource, reason: &'static str) -> IoError {
    IoError::MalformedRecord {
        locus: locus.error_locus(),
        reason,
    }
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
/// genotype prints `.`. With `vcf44` (the header declares VCF ≥ 4.4) the
/// first allele's phasing is rendered as a `|`/`/` prefix exactly when
/// dropping it would change its inferred value (htslib vcf.c:6495-6517,
/// A6-14).
///
/// The BCF GT fast path ([`super::gt`]) builds its byte-pair table with this
/// function (at `vcf44 = false`, since the table is header-independent), so
/// a decoded code is exactly what the text pipeline would see.
pub(crate) fn fmt_gt(out: &mut Vec<u8>, ty: u8, count: usize, data: &[u8], vcf44: bool) {
    let width = type_width(ty);
    let (_, vector_end) = int_sentinels(ty);
    let mark = out.len();
    let mut val0 = 0_i64;
    let mut ploidy = 0_usize;
    let mut anyunphased = false;
    for i in 0..count {
        let v = int_element(ty, &data[i * width..(i + 1) * width]);
        if v == vector_end {
            break;
        }
        if i == 0 {
            val0 = v;
        } else {
            out.push(if v & 1 != 0 { b'|' } else { b'/' });
            anyunphased |= v & 1 == 0;
        }
        let allele = v >> 1;
        if allele == 0 {
            out.push(b'.');
        } else {
            push_u64(out, u64::try_from(allele - 1).unwrap_or(0));
        }
        ploidy += 1;
    }
    if ploidy == 0 {
        out.push(b'.');
    }
    if vcf44 {
        if val0 & 1 != 0 {
            // A phased first allele needs an explicit prefix only when the
            // inference rules would otherwise read it as unphased.
            if (ploidy > 1 && anyunphased) || (ploidy <= 1 && val0 >> 1 == 0) {
                out.insert(mark, b'|');
            }
        } else if (ploidy <= 1 && val0 != 0) || (ploidy > 1 && !anyunphased) {
            // An unphased first allele the inference would read as phased.
            out.insert(mark, b'/');
        }
    }
}

/// One integer element of `ty`.
pub(super) fn int_element(ty: u8, elem: &[u8]) -> i64 {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bcf::BcfDict;

    /// doc_v20 S43: BCF POS is a signed 32-bit 0-based field, so a negative
    /// pos0 in malformed file bytes must error instead of rendering POS=1
    /// (the `unwrap_or(0) + 1` used to swallow the sign; htslib rejects the
    /// record).
    #[test]
    fn s43_render_rejects_negative_pos0() {
        let header_text = "##fileformat=VCFv4.2\n\
            ##contig=<ID=Chr01,length=1000>\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let header: noodles::vcf::Header = header_text.parse().expect("invariant: header");
        let dict_text = crate::bcf::bcf_header_text(&header).expect("invariant: bcf header text");
        let dict = BcfDict::from_header(&dict_text).expect("invariant: dict");
        let mut shared = Vec::new();
        shared.extend_from_slice(&0_i32.to_le_bytes()); // chrom = Chr01 (tid 0)
        shared.extend_from_slice(&(-5_i32).to_le_bytes()); // below the telomere sentinel
        shared.extend_from_slice(&1_i32.to_le_bytes()); // rlen
        shared.extend_from_slice(&FLOAT_MISSING.to_le_bytes()); // qual = .
        shared.extend_from_slice(&(1_u32 << 16).to_le_bytes()); // 1 allele, 0 INFO
        shared.extend_from_slice(&0_u32.to_le_bytes()); // 0 FORMAT, 0 samples
        shared.extend_from_slice(&[0x07]); // ID: empty char vector
        shared.extend_from_slice(&[0x17, b'A']); // REF "A"
        shared.extend_from_slice(&[0x01]); // FILTER: count 0 (".")
        let rec = BcfRecordRef {
            shared: &shared,
            indiv: &[],
        };
        let mut out = Vec::new();
        let err = render_site_part(&rec, &dict, &mut out)
            .expect_err("a pos0 below -1 must be an error");
        assert!(
            err.to_string().contains("telomere sentinel"),
            "{err}"
        );
    }

    /// htslib `bcf_format_gt_v2` (vcf.c:6457-6518): the first allele's phasing
    /// renders as a prefix only when inference would otherwise get it wrong
    /// (A6-14).
    #[test]
    fn fmt_gt_prefixed_phasing_matrix() {
        let cases: &[(&[u8], bool, &[u8])] = &[
            // [allele0 phased, allele1 unphased]: the prefix is required.
            (&[3, 4], true, b"|0/1"),
            (&[3, 4], false, b"0/1"),
            // [unphased 0, phased 1]: explicit `/` for the first allele.
            (&[2, 5], true, b"/0|1"),
            (&[2, 5], false, b"0|1"),
            // All phased: no prefix (inference recovers it).
            (&[3, 5], true, b"0|1"),
            // Haploid unphased: explicit `/`; phased: bare allele.
            (&[4], true, b"/1"),
            (&[5], true, b"1"),
            // Phased missing haploid: `|.` (inference would drop the bit).
            (&[1], true, b"|."),
            (&[1], false, b"."),
            (&[0], true, b"."),
            (&[0, 0], true, b"./."),
        ];
        for (codes, vcf44, want) in cases {
            let mut out = Vec::new();
            fmt_gt(&mut out, value_type::INT8, codes.len(), codes, *vcf44);
            assert_eq!(out, want.to_vec(), "codes {codes:02x?} vcf44={vcf44}");
        }
    }
}
