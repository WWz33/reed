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
    BcfRecordRef, Cursor, FLOAT_MISSING, FLOAT_VECTOR_END, field_bytes, format_fields, gt_element,
    int_element, int_sentinels, malformed, type_width, value_type,
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
    // doc_v28 V28-21: the fixed header decode (including doc_v26 V26-05's
    // zero-allele rejection) is shared, so the renderer cannot drift from the
    // GT path and the encoder.
    let (header, mut c) = super::record::SharedHeader::parse(rec.shared)?;
    let chrom = header.chrom;
    let pos0 = header.pos0;
    let qual_bits = header.qual;
    let n_allele = header.n_allele;
    let n_info = header.n_info;
    let n_fmt = header.n_fmt;
    let n_sample = header.n_sample;

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
        u64::try_from(pos0).map_err(|_| {
            bad(
                (contig_name, 0),
                "BCF position is below the telomere sentinel (-1)",
            )
        })? + 1
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

    // Alleles: REF then ALT (a lone allele prints `.` for ALT).
    out.push(b'\t');
    out.extend_from_slice(c.typed_string()?);
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
    // A sites-only BCF (no FORMAT keys, no samples) has no FORMAT/sample
    // region at all: htslib renders 8 columns, so opening the FORMAT field
    // here would add a phantom ninth column (doc_v22 V22-R1).
    if n_fmt == 0 && n_sample == 0 {
        return Ok(());
    }
    // FORMAT + samples.
    out.push(b'\t');
    if n_fmt == 0 {
        // htslib vcf_format: FORMAT `.` plus one `.` per sample. The
        // leading `\t` opened the FORMAT field, so each item here must
        // carry its own separator — one `\t.` per column, not one extra
        // `\t.` past the (empty) FORMAT column.
        out.push(b'.');
        for _ in 0..n_sample {
            out.extend_from_slice(b"\t.");
        }
    } else {
        // Field descriptors live in the individual block; values are
        // field-major: each field's descriptor is followed by its
        // `n_sample × count` elements (htslib's packed layout). The walk is
        // shared with the column decoder (V22-R11).
        let indiv = rec.indiv;
        let fields = format_fields(indiv, n_fmt, n_sample, |reason| bad(locus, reason))?;

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
                    // htslib `bcf_format_gt_v2` rejects a non-integer GT
                    // type (`default: return -2`); a CHAR-typed GT has a
                    // 1-byte element, so reading it as int32 panics
                    // (doc_v23 V23-02).
                    fmt_gt(out, ty, count, data, dict.vcf44())
                        .map_err(|reason| bad(locus, reason))?;
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
/// values come from [`super::columns::decode_format_columns`]. No trailing newline.
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
    // Sites-only BCF: no FORMAT placeholder, the site line stays 8 columns
    // (doc_v22 V22-R1).
    if n_sample == 0 {
        return Ok(());
    }
    out.extend_from_slice(b"\t.");
    for _ in 0..n_sample {
        out.extend_from_slice(b"\t.");
    }
    Ok(())
}

/// A malformed record at a known locus (borrowed: the `Locus` is built only
/// when the error actually carries it).
fn bad(locus: impl LocusSource, reason: &'static str) -> IoError {
    IoError::MalformedRecord {
        locus: locus.error_locus(),
        reason: reason.to_string(),
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
///
/// # Errors
///
/// The descriptor's type is not an integer type. htslib returns an error for
/// that (`bcf_format_gt_v2`'s `default: return -2`); without the guard a
/// CHAR element (1 byte wide) was read as int32 and panicked on the slice
/// (doc_v23 V23-02).
pub(crate) fn fmt_gt(
    out: &mut Vec<u8>,
    ty: u8,
    count: usize,
    data: &[u8],
    vcf44: bool,
) -> Result<(), &'static str> {
    // htslib represents an all-missing GT field as BCF_BT_NULL with count 0;
    // every sample renders as `.` (doc_v28 V28-08).
    if ty == value_type::NULL && count == 0 {
        out.push(b'.');
        return Ok(());
    }
    if !matches!(ty, value_type::INT8 | value_type::INT16 | value_type::INT32) {
        return Err("BCF GT field is not integer-typed");
    }
    let width = type_width(ty);
    let mark = out.len();
    let mut val0 = 0_i64;
    let mut ploidy = 0_usize;
    let mut anyunphased = false;
    for i in 0..count {
        let Some(v) = gt_element(ty, &data[i * width..(i + 1) * width])? else {
            break;
        };
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
    Ok(())
}

/// One integer element of `ty`.
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
        let err =
            render_site_part(&rec, &dict, &mut out).expect_err("a pos0 below -1 must be an error");
        assert!(err.to_string().contains("telomere sentinel"), "{err}");
    }

    /// doc_v26 V26-05: a record with `n_allele == 0` has no REF at all.
    /// The renderer used to print `.` for REF and ALT and hand on a "valid"
    /// line; it must report `MalformedRecord` instead.
    #[test]
    fn v26_render_rejects_zero_alleles() {
        let header_text = "##fileformat=VCFv4.2\n\
            ##contig=<ID=Chr01,length=1000>\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let header: noodles::vcf::Header = header_text.parse().expect("invariant: header");
        let dict_text = crate::bcf::bcf_header_text(&header).expect("invariant: bcf header text");
        let dict = BcfDict::from_header(&dict_text).expect("invariant: dict");
        let mut shared = Vec::new();
        shared.extend_from_slice(&0_i32.to_le_bytes()); // chrom = Chr01
        shared.extend_from_slice(&9_i32.to_le_bytes()); // 0-based pos 9
        shared.extend_from_slice(&1_i32.to_le_bytes()); // rlen
        shared.extend_from_slice(&FLOAT_MISSING.to_le_bytes()); // qual = .
        shared.extend_from_slice(&0_u32.to_le_bytes()); // n_allele = 0, n_info = 0
        shared.extend_from_slice(&0_u32.to_le_bytes()); // 0 FORMAT, 0 samples
        shared.extend_from_slice(&[0x07]); // ID: empty
        shared.extend_from_slice(&[0x01]); // FILTER: count 0
        let rec = BcfRecordRef {
            shared: &shared,
            indiv: &[],
        };
        let mut out = Vec::new();
        let err = render_site_part(&rec, &dict, &mut out)
            .expect_err("n_allele = 0 must be MalformedRecord");
        assert!(matches!(err, IoError::MalformedRecord { .. }), "{err:?}");
        assert!(err.to_string().contains("zero alleles"), "{err}");
        assert!(
            out.is_empty(),
            "nothing may be rendered for a rejected record"
        );
    }

    /// A BCF record with no FORMAT keys (`n_fmt == 0`, what htslib writes
    /// for a text FORMAT-`.` record) renders FORMAT as `.` and one `.` per
    /// sample — htslib vcf_format's output, byte-checked against
    /// `bcftools view` (found while cross-checking noodles 0.117→0.119
    /// `3b5b1a8c0`: the render used to emit an empty FORMAT column and an
    /// extra phantom sample field, which htslib reads back with a
    /// `FORMAT '' is not defined in the header` warning).
    #[test]
    fn fmt_no_format_keys_render_as_missing() {
        let header_text = "##fileformat=VCFv4.2\n\
            ##contig=<ID=Chr01,length=1000>\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let header: noodles::vcf::Header = header_text.parse().expect("invariant: header");
        let dict_text = crate::bcf::bcf_header_text(&header).expect("invariant: bcf header text");
        let dict = BcfDict::from_header(&dict_text).expect("invariant: dict");
        let mut shared = Vec::new();
        shared.extend_from_slice(&0_i32.to_le_bytes()); // chrom = Chr01 (tid 0)
        shared.extend_from_slice(&19_i32.to_le_bytes()); // pos0 19 → POS 20
        shared.extend_from_slice(&1_i32.to_le_bytes()); // rlen
        shared.extend_from_slice(&FLOAT_MISSING.to_le_bytes()); // qual = .
        shared.extend_from_slice(&(1_u32 << 16).to_le_bytes()); // 1 allele, 0 INFO
        shared.extend_from_slice(&2_u32.to_le_bytes()); // 0 FORMAT, 2 samples
        shared.extend_from_slice(&[0x07]); // ID: empty char vector
        shared.extend_from_slice(&[0x17, b'A']); // REF "A"
        shared.extend_from_slice(&[0x01]); // FILTER: count 0 (".")
        let rec = BcfRecordRef {
            shared: &shared,
            indiv: &[],
        };
        let mut out = Vec::new();
        render_record(&rec, &dict, &mut out).expect("invariant: render");
        assert_eq!(
            out, b"Chr01\t20\t.\tA\t.\t.\t.\t.\t.\t.\t.\n",
            "n_fmt == 0 renders FORMAT `.` plus one `.` per sample (11 columns)"
        );
    }

    /// A sites-only BCF (`n_fmt == 0`, `n_sample == 0`) renders 8 columns:
    /// no FORMAT column and no phantom sample field (doc_v22 V22-R1).
    #[test]
    fn sites_only_bcf_stays_eight_columns() {
        let header_text = "##fileformat=VCFv4.2\n\
            ##contig=<ID=Chr01,length=1000>\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
        let header: noodles::vcf::Header = header_text.parse().expect("invariant: header");
        let dict_text = crate::bcf::bcf_header_text(&header).expect("invariant: bcf header text");
        let dict = BcfDict::from_header(&dict_text).expect("invariant: dict");
        let mut shared = Vec::new();
        shared.extend_from_slice(&0_i32.to_le_bytes()); // chrom = Chr01 (tid 0)
        shared.extend_from_slice(&19_i32.to_le_bytes()); // pos0 19 → POS 20
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
        render_record(&rec, &dict, &mut out).expect("invariant: render");
        assert_eq!(
            out, b"Chr01\t20\t.\tA\t.\t.\t.\t.\n",
            "sites-only BCF renders exactly 8 columns"
        );
        // The native filter path's decision line is 8 columns too.
        let mut site_line = Vec::new();
        render_site_line(&rec, &dict, &mut site_line).expect("invariant: site line");
        assert_eq!(site_line, b"Chr01\t20\t.\tA\t.\t.\t.\t.");
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
            fmt_gt(&mut out, value_type::INT8, codes.len(), codes, *vcf44)
                .expect("invariant: INT8 is integer-typed");
            assert_eq!(out, want.to_vec(), "codes {codes:02x?} vcf44={vcf44}");
        }
    }

    /// doc_v28 V28-04: a negative GT code is malformed, not REF allele 0.
    #[test]
    fn negative_gt_element_is_rejected_not_rendered_as_ref() {
        let mut out = Vec::new();
        let err = fmt_gt(
            &mut out,
            value_type::INT8,
            2,
            &[i8::MIN as u8, i8::MIN as u8],
            false,
        )
        .expect_err("negative GT codes are corrupt");
        assert!(err.contains("negative"), "{err}");
        // Vector-end still truncates the allele vector.
        let mut out = Vec::new();
        fmt_gt(
            &mut out,
            value_type::INT8,
            2,
            &[2, (i8::MIN + 1) as u8],
            false,
        )
        .expect("vector-end is legal");
        assert_eq!(out, b"0");

        // V28-08: NULL/count=0 is an all-missing GT, not a type error.
        let mut out = Vec::new();
        fmt_gt(&mut out, value_type::NULL, 0, &[], true).expect("NULL GT is all-missing");
        assert_eq!(out, b".");
    }

    /// doc_v23 V23-02: a GT descriptor whose type nibble is CHAR has a 1-byte
    /// element, so reading it as int32 must be an error, not a panic.
    /// htslib `bcf_format_gt_v2` rejects the type (`default: return -2`).
    #[test]
    fn char_typed_gt_is_malformed_not_an_out_of_bounds_panic() {
        let mut out = Vec::new();
        let err = fmt_gt(&mut out, value_type::CHAR, 2, b"01", true)
            .expect_err("CHAR is not an integer GT type");
        assert_eq!(err, "BCF GT field is not integer-typed");
        assert!(out.is_empty(), "nothing is rendered for a rejected field");
        // FLOAT is rejected the same way (its element would be readable but
        // semantically invalid).
        assert!(fmt_gt(&mut out, value_type::FLOAT, 1, &[0, 0, 0, 0], false).is_err());
    }

    /// doc_v23 V23-02: the render path (not just `fmt_gt`) reports a
    /// CHAR-typed GT as a malformed record at the record's locus.
    #[test]
    fn render_record_rejects_char_typed_gt() {
        let header_text = "##fileformat=VCFv4.2\n\
            ##contig=<ID=Chr01,length=1000>\n\
            ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";
        let header: noodles::vcf::Header = header_text.parse().expect("invariant: header");
        let dict_text = crate::bcf::bcf_header_text(&header).expect("invariant: bcf header text");
        let dict = BcfDict::from_header(&dict_text).expect("invariant: dict");
        let mut shared = Vec::new();
        shared.extend_from_slice(&0_i32.to_le_bytes()); // chrom = Chr01 (tid 0)
        shared.extend_from_slice(&19_i32.to_le_bytes()); // pos0 19 -> POS 20
        shared.extend_from_slice(&1_i32.to_le_bytes()); // rlen
        shared.extend_from_slice(&FLOAT_MISSING.to_le_bytes()); // qual = .
        shared.extend_from_slice(&(1_u32 << 16).to_le_bytes()); // 1 allele, 0 INFO
        shared.extend_from_slice(&(1_u32 << 24 | 1).to_le_bytes()); // 1 FORMAT, 1 sample
        shared.extend_from_slice(&[0x07]); // ID: empty char vector
        shared.extend_from_slice(&[0x17, b'A']); // REF "A"
        shared.extend_from_slice(&[0x01]); // FILTER: count 0 (".")
        // FORMAT: GT key, then a descriptor with count 2 and type CHAR.
        let gt = dict
            .key_of(b"GT")
            .expect("invariant: GT is in the dictionary");
        let mut indiv = Vec::new();
        indiv.push(0x11); // typed int: count 1, INT8
        indiv.push(u8::try_from(gt).expect("invariant: small dictionary"));
        indiv.extend_from_slice(&[0x27, b'0', b'1']); // count 2, CHAR, "01"
        let rec = BcfRecordRef {
            shared: &shared,
            indiv: &indiv,
        };
        let mut out = Vec::new();
        let err = render_record(&rec, &dict, &mut out)
            .expect_err("a CHAR-typed GT is a malformed record");
        let text = err.to_string();
        assert!(text.contains("Chr01:20"), "{text}");
        assert!(text.contains("not integer-typed"), "{text}");
    }
}
