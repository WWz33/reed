//! Native record → BCF encoder (doc_v3 02-bcf.md B4, doc_v6 03 §8).
//!
//! The noodles encoder costs ~800 µs per 348-sample record (it re-parses
//! every formatted value), which breaks the `-Ob` budget. This encoder
//! applies the filter's edits while it walks the *source* record — VCF text
//! ([`NativeEncoder::encode_text`]) or a framed BCF record
//! ([`NativeEncoder::encode_bcf`]) — so the edited record never round-trips
//! through text. Unmodified BCF FORMAT columns copy their raw element bytes;
//! only INFO and the GT column are rewritten.
//!
//! Widths follow htslib's `bcf_enc_vint` (the smallest of int8/int16/int32
//! covering every element, legal ranges −120…127 / −32760…32767), and
//! dictionary indices come from the output header (see
//! [`super::encode::bcf_header_text`]).

// Casts in this file are bounded by construction: Integer widths are chosen from the parsed value ranges.
#![allow(clippy::cast_possible_truncation)]

use std::collections::HashMap;

use reed_core::SampleId;

use crate::error::{IoError, Locus};
use crate::vcf::info::info_fragments;
use crate::vcf::line::LineRef;
use crate::vcf::line_edit::{EditPlan, InfoValueOut, render_info_value};
use crate::vcf::site::SiteView;
use noodles::vcf;

use super::dict::BcfDict;
use super::record::{
    BcfRecordRef, Cursor, FLOAT_MISSING, FLOAT_VECTOR_END, INT8_MAX_VALUE, INT8_MIN_VALUE,
    INT16_MAX_VALUE, INT16_MIN_VALUE, INT32_MAX_VALUE, INT32_MIN_VALUE, field_bytes, malformed,
    type_width, value_type,
};
use super::render::int_element;

/// A per-worker native encoder.
pub struct NativeEncoder {
    contig_index: HashMap<Vec<u8>, i32>,
    /// Tag name → dictionary index (PASS included).
    tags: HashMap<Vec<u8>, u32>,
    /// INFO tag → declared value shape.
    info_types: HashMap<Vec<u8>, InfoType>,
    /// FORMAT tag → declared BCF type.
    format_types: HashMap<Vec<u8>, u8>,
    /// Samples in the *output* header (the sample subset after pruning).
    n_samples: usize,
    /// Scratch buffers, reused across records.
    shared: Vec<u8>,
    indiv: Vec<u8>,
    /// Byte ranges of the current record's sample subfields, flattened as
    /// `sample * n_keys + key` (one pass over the sample text per record).
    cells: Vec<(u32, u32)>,
    ints: Vec<i64>,
    /// One sample's decoded GT codes (reused; the text GT walk allocates
    /// nothing per record).
    codes: Vec<i32>,
    /// Field headers of the current BCF source record `(key, ty, count,
    /// payload offset)`, reused across records.
    bcf_fields: Vec<(u32, u8, usize, usize)>,
}

/// The value shape an INFO tag declares.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum InfoType {
    Flag,
    Integer,
    Float,
    String,
}

impl NativeEncoder {
    /// Builds an encoder for `header` (the output header).
    ///
    /// # Errors
    ///
    /// [`IoError::Io`] when the header's maps cannot be built.
    pub fn new(header: &vcf::Header) -> Result<Self, IoError> {
        let maps = vcf::header::StringMaps::try_from(header)
            .map_err(|e| IoError::Io(format!("BCF header maps: {e}")))?;
        use noodles::vcf::header::record::value::map::{format, info};
        let mut info_types = HashMap::new();
        for (key, map) in header.infos() {
            let ty = match map.ty() {
                info::Type::Integer => InfoType::Integer,
                info::Type::Float => InfoType::Float,
                info::Type::Flag => InfoType::Flag,
                info::Type::String | info::Type::Character => InfoType::String,
            };
            info_types.insert(key.as_bytes().to_vec(), ty);
        }
        let mut format_types = HashMap::new();
        for (key, map) in header.formats() {
            let ty = match map.ty() {
                format::Type::Integer => value_type::INT32,
                format::Type::Float => value_type::FLOAT,
                _ => value_type::CHAR,
            };
            format_types.insert(key.as_bytes().to_vec(), ty);
        }
        let mut tags = HashMap::new();
        for key in header
            .filters()
            .keys()
            .chain(header.infos().keys())
            .chain(header.formats().keys())
        {
            if let Some(idx) = maps.strings().get_index_of(key) {
                tags.insert(
                    key.as_bytes().to_vec(),
                    u32::try_from(idx).unwrap_or(u32::MAX),
                );
            }
        }
        let mut contig_index = HashMap::new();
        for name in header.contigs().keys() {
            if let Some(idx) = maps.contigs().get_index_of(name) {
                contig_index.insert(
                    name.as_bytes().to_vec(),
                    i32::try_from(idx).unwrap_or(i32::MAX),
                );
            }
        }
        Ok(Self {
            contig_index,
            tags,
            info_types,
            format_types,
            n_samples: header.sample_names().len(),
            shared: Vec::with_capacity(4096),
            indiv: Vec::with_capacity(1 << 16),
            cells: Vec::new(),
            ints: Vec::with_capacity(64),
            codes: Vec::with_capacity(8),
            bcf_fields: Vec::new(),
        })
    }

    /// Number of sample columns the encoder expects.
    #[must_use]
    pub const fn n_samples(&self) -> usize {
        self.n_samples
    }

    /// Encodes one VCF text record into `out` (the BCF record bytes,
    /// including the two length prefixes), applying `edit` in the same walk
    /// — with a default edit the output equals the pre-Y4 encode of the
    /// verbatim line.
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] for malformed input and
    /// [`IoError::Header`] for tags the output header does not declare.
    pub fn encode_text(
        &mut self,
        line: &LineRef<'_>,
        site: &SiteView<'_>,
        n_samples_in: usize,
        edit: &EditPlan<'_>,
        out: &mut Vec<u8>,
    ) -> Result<(), IoError> {
        debug_assert_eq!(
            n_samples_out(edit, n_samples_in),
            self.n_samples,
            "encoder built from the output header"
        );
        self.shared.clear();
        self.indiv.clear();
        self.encode_shared(line, site, edit)?;
        self.encode_indiv(line, n_samples_in, edit)?;
        self.finish(out)?;
        Ok(())
    }

    /// Encodes one framed BCF record into `out`, applying `edit` in the walk
    /// (doc_v6 03 §8): dictionary indices are translated through the input
    /// dictionary, unmodified FORMAT columns copy their raw element bytes,
    /// and only INFO and the GT column are rewritten.
    ///
    /// # Errors
    ///
    /// As [`Self::encode_text`].
    pub fn encode_bcf(
        &mut self,
        rec: &BcfRecordRef<'_>,
        dict: &BcfDict,
        edit: &EditPlan<'_>,
        out: &mut Vec<u8>,
    ) -> Result<(), IoError> {
        let mut head = Cursor::new(rec.shared);
        let in_chrom = head.i32()?;
        let pos0 = head.i32()?;
        let rlen = head.i32()?;
        let qual = head.u32()?;
        let na_ni = head.u32()?;
        let nf_ns = head.u32()?;
        let n_allele = (na_ni >> 16) as usize;
        let n_info = (na_ni & 0xFFFF) as usize;
        let n_fmt = (nf_ns >> 24) as usize;
        let n_sample_in = (nf_ns & 0x00FF_FFFF) as usize;
        debug_assert_eq!(
            n_samples_out(edit, n_sample_in),
            self.n_samples,
            "encoder built from the output header"
        );
        self.shared.clear();
        self.indiv.clear();
        self.encode_shared_bcf(
            rec, dict, edit, in_chrom, pos0, rlen, qual, n_allele, n_info,
        )?;
        self.encode_indiv_bcf(rec, dict, edit, n_fmt, n_sample_in)?;
        self.finish(out)?;
        Ok(())
    }

    /// Appends the length prefixes and the two blocks to `out`.
    fn finish(&mut self, out: &mut Vec<u8>) -> Result<(), IoError> {
        let l_shared = u32::try_from(self.shared.len())
            .map_err(|_| bad("BCF shared block exceeds u32::MAX bytes"))?;
        let l_indiv = u32::try_from(self.indiv.len())
            .map_err(|_| bad("BCF individual block exceeds u32::MAX bytes"))?;
        out.extend_from_slice(&l_shared.to_le_bytes());
        out.extend_from_slice(&l_indiv.to_le_bytes());
        out.extend_from_slice(&self.shared);
        out.extend_from_slice(&self.indiv);
        Ok(())
    }

    fn undeclared(what: &str, name: &[u8]) -> IoError {
        IoError::Header(format!(
            "BCF output header does not declare {what} {:?}",
            String::from_utf8_lossy(name)
        ))
    }

    fn encode_shared(
        &mut self,
        lr: &LineRef<'_>,
        site: &SiteView<'_>,
        edit: &EditPlan<'_>,
    ) -> Result<(), IoError> {
        let chrom = lr.chrom();
        let contig = self
            .contig_index
            .get(chrom)
            .copied()
            .ok_or_else(|| Self::undeclared("contig", chrom))?;
        let pos0 = i32::try_from(site.pos1.saturating_sub(1)).unwrap_or(i32::MAX);
        let rlen = i32::try_from(site.reference.len()).unwrap_or(1).max(1);
        // The site header, ID, alleles and FILTER: a scoped borrow of the
        // shared block (the INFO section below calls `&mut self` helpers).
        {
            let s = &mut self.shared;
            s.extend_from_slice(&contig.to_le_bytes());
            s.extend_from_slice(&pos0.to_le_bytes());
            s.extend_from_slice(&rlen.to_le_bytes());
            s.extend_from_slice(&site.qual.map_or(FLOAT_MISSING, f32::to_bits).to_le_bytes());
            s.extend_from_slice(&[0; 8]); // n_allele|n_info, n_fmt|n_sample

            // ID: an empty ID (".") is stored as a zero-length string.
            let id = lr.id();
            encode_typed_string(s, if id == b"." { &[] } else { id });

            // Alleles: REF then the ALT list ("." means no ALT).
            let alts = lr.alts();
            if alts == b"." {
                encode_typed_string(s, lr.reference());
            } else {
                encode_typed_string(s, lr.reference());
                for alt in alts.split(|&b| b == b',') {
                    encode_typed_string(s, alt);
                }
            }
        }
        let n_allele = if lr.alts() == b"." {
            1_u32
        } else {
            // REF + the comma-separated ALT list.
            2 + u32::try_from(memchr::memchr_iter(b',', lr.alts()).count()).unwrap_or(u32::MAX)
        };

        // FILTER.
        let filter = lr.filter_text();
        if filter == b"." || filter.is_empty() {
            encode_size(&mut self.shared, 0, value_type::INT8);
        } else {
            let mut ints = std::mem::take(&mut self.ints);
            ints.clear();
            for name in filter.split(|&b| b == b';') {
                let idx = self
                    .tags
                    .get(name)
                    .copied()
                    .ok_or_else(|| Self::undeclared("FILTER", name))?;
                ints.push(i64::from(idx));
            }
            encode_int_vector(&mut self.shared, &ints);
            self.ints = ints;
        }

        // INFO: one walk of the source fragments with the edit applied —
        // stale tags drop out, recomputed tags replace in place (typed, no
        // text round-trip), the rest copy their parsed source values.
        let mut n_info = 0_u32;
        let mut seen = SeenMarks::new(edit.info_set.len());
        let info = lr.info_text();
        if info != b"." && !info.is_empty() {
            for (key, value) in info_fragments(info) {
                if edit.info_remove.contains(&key) {
                    continue;
                }
                let idx = self
                    .tags
                    .get(key)
                    .copied()
                    .ok_or_else(|| Self::undeclared("INFO", key))?;
                let declared = self
                    .info_types
                    .get(key)
                    .copied()
                    .unwrap_or(InfoType::String);
                encode_typed_int(&mut self.shared, i64::from(idx));
                if let Some(i) = edit.info_set.iter().position(|(t, _)| *t == key) {
                    seen.mark(i);
                    encode_info_typed(&mut self.shared, declared, &edit.info_set[i].1);
                } else {
                    self.encode_info_text(declared, value);
                }
                n_info += 1;
            }
        }
        self.encode_info_appends(edit, &seen, &mut n_info)?;
        // Patch the packed counters (n_allele<<16|n_info; n_fmt|n_sample is
        // filled in by `encode_indiv`).
        let packed = pack_allele_info(n_allele, n_info)?;
        self.shared[16..20].copy_from_slice(&packed.to_le_bytes());
        Ok(())
    }

    /// One INFO value parsed from the source text (the pre-Y4 behavior for
    /// values the edit does not touch).
    ///
    /// The value's shape follows htslib's `vcf_parse_info`:
    ///
    /// - a Flag **with** a value is stored as a string — its comment reads
    ///   "if Flag has a value, treat it as a string" (vcf.c:3949-3951). The
    ///   VCF spec forbids the combination (a Flag "does not contain a Value
    ///   entry", v4.5 §1.4.2), so this is the recovery path for foreign
    ///   input; dropping the value (or flipping `DB=0` to a set flag) is
    ///   worse than keeping it (doc_v12 D12-4);
    /// - a String value of `.` is a **literal** dot. Only the whole INFO
    ///   column being `.` means "no INFO" (v4.5 §1.6.1); inside a value the
    ///   dot has no special meaning, so it must round-trip as `TAG=.`
    ///   (doc_v12 D12-5). The empty-string special case stays where htslib
    ///   really has one: the ID column.
    fn encode_info_text(&mut self, declared: InfoType, value: Option<&[u8]>) {
        let s = &mut self.shared;
        match (declared, value) {
            // A bare flag carries no value at all.
            (InfoType::Flag, None) => encode_size(s, 0, value_type::NULL),
            (InfoType::Integer, Some(v)) => {
                let mut ints = std::mem::take(&mut self.ints);
                ints.clear();
                parse_int_vector(v, &mut ints);
                encode_int_vector(s, &ints);
                self.ints = ints;
            }
            (InfoType::Float, Some(v)) => {
                let mut floats = Vec::new();
                parse_float_vector(v, &mut floats);
                encode_size(s, floats.len(), value_type::FLOAT);
                for f in floats {
                    s.extend_from_slice(&f.to_le_bytes());
                }
            }
            (_, Some(v)) => encode_typed_string(s, v),
            (InfoType::Integer, None) => {
                encode_size(s, 1, value_type::INT8);
                s.push(i8::MIN as u8);
            }
            (InfoType::Float, None) => {
                encode_size(s, 1, value_type::FLOAT);
                s.extend_from_slice(&FLOAT_MISSING.to_le_bytes());
            }
            (InfoType::String, None) => encode_typed_string(s, b""),
        }
    }

    /// The `info_set` entries the source walk did not replace in place,
    /// appended in `info_set` order (the `edit_info` semantics).
    fn encode_info_appends(
        &mut self,
        edit: &EditPlan<'_>,
        seen: &SeenMarks,
        n_info: &mut u32,
    ) -> Result<(), IoError> {
        for (i, (tag, value)) in edit.info_set.iter().enumerate() {
            if seen.get(i) {
                continue;
            }
            let idx = self
                .tags
                .get(*tag)
                .copied()
                .ok_or_else(|| Self::undeclared("INFO", tag))?;
            let declared = self
                .info_types
                .get(*tag)
                .copied()
                .unwrap_or(InfoType::String);
            encode_typed_int(&mut self.shared, i64::from(idx));
            encode_info_typed(&mut self.shared, declared, value);
            *n_info += 1;
        }
        Ok(())
    }

    fn encode_indiv(
        &mut self,
        lr: &LineRef<'_>,
        n_samples_in: usize,
        edit: &EditPlan<'_>,
    ) -> Result<(), IoError> {
        let format = lr.format_text();
        let n_fmt = if format == b"." || format.is_empty() {
            0_u32
        } else {
            u32::try_from(memchr::memchr_iter(b':', format).count() + 1).unwrap_or(u32::MAX)
        };
        if n_fmt > 0 {
            let samples = lr.samples_text();
            let keys: Vec<&[u8]> = format.split(|&b| b == b':').collect();
            let n_keys = keys.len();
            let gt_index = keys.iter().position(|k| *k == b"GT");
            // One pass over the sample text: every subfield's byte range.
            self.cells.clear();
            self.cells.resize(n_samples_in * n_keys, ABSENT);
            let mut column_start = 0_usize;
            for sample in 0..n_samples_in {
                let rest = &samples[column_start..];
                let column_len = rest.iter().position(|&b| b == b'\t').unwrap_or(rest.len());
                let column = &rest[..column_len];
                let mut sub_start = 0_usize;
                for k in 0..n_keys {
                    let (value, next) = match column[sub_start..].iter().position(|&b| b == b':') {
                        Some(i) => (&column[sub_start..sub_start + i], sub_start + i + 1),
                        None => (&column[sub_start..], usize::MAX),
                    };
                    let base = u32::try_from(column_start + sub_start).unwrap_or(u32::MAX);
                    self.cells[sample * n_keys + k] =
                        (base, base + u32::try_from(value.len()).unwrap_or(u32::MAX));
                    if next == usize::MAX {
                        break;
                    }
                    sub_start = next;
                }
                column_start += column_len + 1;
                if column_start > samples.len() {
                    break;
                }
            }
            for (k, key) in keys.iter().enumerate() {
                let idx = self
                    .tags
                    .get(*key)
                    .copied()
                    .ok_or_else(|| Self::undeclared("FORMAT", key))?;
                let ty = self
                    .format_types
                    .get(*key)
                    .copied()
                    .unwrap_or(value_type::CHAR);
                encode_typed_int(&mut self.indiv, i64::from(idx));
                if Some(k) == gt_index {
                    self.encode_gt_column(samples, k, n_keys, n_samples_in, edit);
                } else {
                    match ty {
                        value_type::FLOAT => {
                            self.encode_float_column(samples, k, n_keys, n_samples_in, edit);
                        }
                        value_type::CHAR => {
                            self.encode_string_column(samples, k, n_keys, n_samples_in, edit);
                        }
                        _ => self.encode_int_column(samples, k, n_keys, n_samples_in, edit)?,
                    }
                }
            }
        }
        let s = &mut self.shared;
        s[20..24].copy_from_slice(
            &pack_fmt_sample(n_fmt, n_samples_out(edit, n_samples_in) as u32)?.to_le_bytes(),
        );
        Ok(())
    }

    /// One sample's subfield `k` (`ABSENT` when the column is shorter).
    fn cell<'a>(
        cells: &[(u32, u32)],
        samples: &'a [u8],
        sample: usize,
        k: usize,
        n_keys: usize,
    ) -> &'a [u8] {
        match cells.get(sample * n_keys + k).copied() {
            Some((start, end)) if start != u32::MAX => &samples[start as usize..end as usize],
            _ => b".",
        }
    }

    fn encode_gt_column(
        &mut self,
        samples: &[u8],
        k: usize,
        n_keys: usize,
        n_samples_in: usize,
        edit: &EditPlan<'_>,
    ) {
        let mut max_ploidy = 0_usize;
        let mut max_allele = 0_i64;
        for s in 0..n_samples_in {
            if !kept(edit, s) {
                continue;
            }
            if masked(edit, s) {
                // `./.`: two missing alleles.
                max_ploidy = max_ploidy.max(2);
                continue;
            }
            let (ploidy, allele) = gt_shape(Self::cell(&self.cells, samples, s, k, n_keys));
            max_ploidy = max_ploidy.max(ploidy);
            max_allele = max_allele.max(allele);
        }
        let max_ploidy = max_ploidy.max(1);
        let (ty, width) = gt_width(max_allele);
        encode_size(&mut self.indiv, max_ploidy, ty);
        let (_, vector_end) = super::record::int_sentinels(ty);
        let vector_end = vector_end as i32;
        for s in 0..n_samples_in {
            if !kept(edit, s) {
                continue;
            }
            if masked(edit, s) {
                push_int(&mut self.indiv, 0, width);
                push_int(&mut self.indiv, 0, width);
                for _ in 2..max_ploidy {
                    push_int(&mut self.indiv, vector_end, width);
                }
                continue;
            }
            gt_codes_into(
                Self::cell(&self.cells, samples, s, k, n_keys),
                &mut self.codes,
            );
            for i in 0..max_ploidy {
                push_int(
                    &mut self.indiv,
                    self.codes.get(i).copied().unwrap_or(vector_end),
                    width,
                );
            }
        }
    }

    fn encode_int_column(
        &mut self,
        samples: &[u8],
        k: usize,
        n_keys: usize,
        n_samples_in: usize,
        edit: &EditPlan<'_>,
    ) -> Result<(), IoError> {
        let mut max_len = 0_usize;
        let (mut min, mut max) = (i64::MAX, i64::MIN);
        for s in 0..n_samples_in {
            if !kept(edit, s) {
                continue;
            }
            let sub = Self::cell(&self.cells, samples, s, k, n_keys);
            if sub == b"." || sub.is_empty() {
                max_len = max_len.max(1);
                continue;
            }
            let mut len = 0_usize;
            for part in sub.split(|&b| b == b',') {
                len += 1;
                if part == b"." || part.is_empty() {
                    continue;
                }
                let v = parse_i64(part).ok_or_else(|| bad("non-integer FORMAT value"))?;
                min = min.min(v);
                max = max.max(v);
            }
            max_len = max_len.max(len);
        }
        let max_len = max_len.max(1);
        let (ty, width) = int_width(min, max);
        encode_size(&mut self.indiv, max_len, ty);
        let (missing, vector_end) = super::record::int_sentinels(ty);
        let (missing, vector_end) = (missing as i32, vector_end as i32);
        for s in 0..n_samples_in {
            if !kept(edit, s) {
                continue;
            }
            let sub = Self::cell(&self.cells, samples, s, k, n_keys);
            let mut written = 0_usize;
            if sub != b"." && !sub.is_empty() {
                for part in sub.split(|&b| b == b',') {
                    let v = if part == b"." || part.is_empty() {
                        missing
                    } else {
                        let v = parse_i64(part).ok_or_else(|| bad("non-integer FORMAT value"))?;
                        i32::try_from(v.clamp(INT32_MIN_VALUE, INT32_MAX_VALUE)).unwrap_or(missing)
                    };
                    push_int(&mut self.indiv, v, width);
                    written += 1;
                }
            }
            while written < max_len {
                let v = if written == 0 { missing } else { vector_end };
                push_int(&mut self.indiv, v, width);
                written += 1;
            }
        }
        Ok(())
    }

    fn encode_float_column(
        &mut self,
        samples: &[u8],
        k: usize,
        n_keys: usize,
        n_samples_in: usize,
        edit: &EditPlan<'_>,
    ) {
        let mut max_len = 0_usize;
        for s in 0..n_samples_in {
            if !kept(edit, s) {
                continue;
            }
            let sub = Self::cell(&self.cells, samples, s, k, n_keys);
            if sub == b"." || sub.is_empty() {
                max_len = max_len.max(1);
            } else {
                max_len = max_len.max(sub.split(|&b| b == b',').count());
            }
        }
        let max_len = max_len.max(1);
        encode_size(&mut self.indiv, max_len, value_type::FLOAT);
        for s in 0..n_samples_in {
            if !kept(edit, s) {
                continue;
            }
            let sub = Self::cell(&self.cells, samples, s, k, n_keys);
            let mut written = 0_usize;
            if sub != b"." && !sub.is_empty() {
                for part in sub.split(|&b| b == b',') {
                    let bits = if part == b"." || part.is_empty() {
                        FLOAT_MISSING
                    } else {
                        std::str::from_utf8(part)
                            .ok()
                            .and_then(|t| t.parse::<f32>().ok())
                            .map_or(FLOAT_MISSING, f32::to_bits)
                    };
                    self.indiv.extend_from_slice(&bits.to_le_bytes());
                    written += 1;
                }
            }
            while written < max_len {
                let bits = if written == 0 {
                    FLOAT_MISSING
                } else {
                    FLOAT_VECTOR_END
                };
                self.indiv.extend_from_slice(&bits.to_le_bytes());
                written += 1;
            }
        }
    }

    fn encode_string_column(
        &mut self,
        samples: &[u8],
        k: usize,
        n_keys: usize,
        n_samples_in: usize,
        edit: &EditPlan<'_>,
    ) {
        // The width is the longest value; a missing value is the literal `.`.
        let mut width = 1_usize;
        for s in 0..n_samples_in {
            if !kept(edit, s) {
                continue;
            }
            let sub = Self::cell(&self.cells, samples, s, k, n_keys);
            width = width.max(if sub.is_empty() { 1 } else { sub.len() });
        }
        encode_size(&mut self.indiv, width, value_type::CHAR);
        for s in 0..n_samples_in {
            if !kept(edit, s) {
                continue;
            }
            let sub = Self::cell(&self.cells, samples, s, k, n_keys);
            let value: &[u8] = if sub.is_empty() { b"." } else { sub };
            self.indiv.extend_from_slice(value);
            for _ in value.len()..width {
                self.indiv.push(0);
            }
        }
    }

    /// The shared block of a BCF source record: dictionary indices
    /// (contig/FILTER/INFO keys) are translated to the output header; ID,
    /// alleles and untouched INFO values copy their raw bytes.
    #[allow(clippy::too_many_arguments)]
    fn encode_shared_bcf(
        &mut self,
        rec: &BcfRecordRef<'_>,
        dict: &BcfDict,
        edit: &EditPlan<'_>,
        in_chrom: i32,
        pos0: i32,
        rlen: i32,
        qual: u32,
        n_allele: usize,
        n_info_in: usize,
    ) -> Result<(), IoError> {
        let contig = dict
            .contig(in_chrom)
            .ok_or_else(|| malformed("BCF contig index is not in the dictionary"))?;
        let out_contig = self
            .contig_index
            .get(contig)
            .copied()
            .ok_or_else(|| Self::undeclared("contig", contig))?;
        {
            let s = &mut self.shared;
            s.extend_from_slice(&out_contig.to_le_bytes());
            s.extend_from_slice(&pos0.to_le_bytes());
            s.extend_from_slice(&rlen.to_le_bytes());
            s.extend_from_slice(&qual.to_le_bytes());
            s.extend_from_slice(&[0; 8]); // n_allele|n_info, n_fmt|n_sample
        }
        let mut c = Cursor::new(&rec.shared[24..]);
        // ID and alleles are format-neutral typed strings: raw copy.
        let id = c.typed_string()?;
        encode_typed_string(&mut self.shared, id);
        for _ in 0..n_allele {
            let allele = c.typed_string()?;
            encode_typed_string(&mut self.shared, allele);
        }
        // FILTER: translate the indices.
        let (f_ty, f_count) = c.typed()?;
        let mut ints = std::mem::take(&mut self.ints);
        ints.clear();
        for _ in 0..f_count {
            let key = c.int_value(f_ty)?;
            let name = dict
                .tag(u32::try_from(key).map_err(|_| malformed("BCF FILTER key is negative"))?)
                .ok_or_else(|| malformed("BCF FILTER index is not in the dictionary"))?;
            let idx = self
                .tags
                .get(name)
                .copied()
                .ok_or_else(|| Self::undeclared("FILTER", name))?;
            ints.push(i64::from(idx));
        }
        encode_int_vector(&mut self.shared, &ints);
        self.ints = ints;
        // INFO: the edit walk of encode_shared, with raw value bytes.
        let mut n_info = 0_u32;
        let mut seen = SeenMarks::new(edit.info_set.len());
        for _ in 0..n_info_in {
            let key = c.typed_int()?;
            let name = dict
                .tag(u32::try_from(key).map_err(|_| malformed("BCF INFO key is negative"))?)
                .ok_or_else(|| malformed("BCF INFO index is not in the dictionary"))?;
            let (ty, count) = c.typed()?;
            let payload = c.take(field_bytes(ty, count, 1)?)?;
            if edit.info_remove.contains(&name) {
                continue;
            }
            let idx = self
                .tags
                .get(name)
                .copied()
                .ok_or_else(|| Self::undeclared("INFO", name))?;
            encode_typed_int(&mut self.shared, i64::from(idx));
            if let Some(i) = edit.info_set.iter().position(|(t, _)| *t == name) {
                seen.mark(i);
                let declared = self
                    .info_types
                    .get(name)
                    .copied()
                    .unwrap_or(InfoType::String);
                encode_info_typed(&mut self.shared, declared, &edit.info_set[i].1);
            } else {
                encode_size(&mut self.shared, count, ty);
                self.shared.extend_from_slice(payload);
            }
            n_info += 1;
        }
        self.encode_info_appends(edit, &seen, &mut n_info)?;
        let n_allele = u32::try_from(n_allele).unwrap_or(u32::MAX);
        let packed = pack_allele_info(n_allele, n_info)?;
        self.shared[16..20].copy_from_slice(&packed.to_le_bytes());
        Ok(())
    }

    /// The individual block of a BCF source record: unmodified FORMAT
    /// columns copy their raw element bytes (dropped samples skip their
    /// slice); the GT column is re-encoded (mask/drop) from the source
    /// codes.
    fn encode_indiv_bcf(
        &mut self,
        rec: &BcfRecordRef<'_>,
        dict: &BcfDict,
        edit: &EditPlan<'_>,
        n_fmt: usize,
        n_sample_in: usize,
    ) -> Result<(), IoError> {
        let indiv = rec.indiv;
        let mut ci = Cursor::new(indiv);
        self.bcf_fields.clear();
        for _ in 0..n_fmt {
            let key = ci.typed_int()?;
            let key = u32::try_from(key).map_err(|_| malformed("BCF FORMAT key is negative"))?;
            let (ty, count) = ci.typed()?;
            let start = indiv.len() - ci.remaining();
            ci.take(field_bytes(ty, count, n_sample_in)?)?;
            self.bcf_fields.push((key, ty, count, start));
        }
        for f in 0..self.bcf_fields.len() {
            let (key, ty, count, start) = self.bcf_fields[f];
            let name = dict
                .tag(key)
                .ok_or_else(|| malformed("BCF FORMAT index is not in the dictionary"))?;
            let out_idx = self
                .tags
                .get(name)
                .copied()
                .ok_or_else(|| Self::undeclared("FORMAT", name))?;
            encode_typed_int(&mut self.indiv, i64::from(out_idx));
            let width = type_width(ty);
            let payload = &indiv[start..start + width * count * n_sample_in];
            if name == b"GT" {
                self.encode_gt_column_bcf(payload, ty, count, n_sample_in, edit)?;
            } else {
                encode_size(&mut self.indiv, count, ty);
                if edit.keep_samples.is_none() {
                    self.indiv.extend_from_slice(payload);
                } else {
                    let stride = width * count;
                    for s in 0..n_sample_in {
                        if kept(edit, s) {
                            self.indiv
                                .extend_from_slice(&payload[s * stride..(s + 1) * stride]);
                        }
                    }
                }
            }
        }
        let n_fmt = u32::try_from(n_fmt).unwrap_or(u32::MAX);
        self.shared[20..24].copy_from_slice(
            &pack_fmt_sample(n_fmt, n_samples_out(edit, n_sample_in) as u32)?.to_le_bytes(),
        );
        Ok(())
    }

    /// The GT column of a BCF source record: codes copy through (the width
    /// follows the largest allele, as the text path's does), masked samples
    /// become two missing alleles, dropped samples are skipped.
    fn encode_gt_column_bcf(
        &mut self,
        payload: &[u8],
        ty: u8,
        count: usize,
        n_sample_in: usize,
        edit: &EditPlan<'_>,
    ) -> Result<(), IoError> {
        let width = type_width(ty);
        if width == 0 {
            return Err(malformed("BCF GT field is not integer-typed"));
        }
        let (_, vector_end) = super::record::int_sentinels(ty);
        let stride = width * count;
        let mut max_ploidy = 0_usize;
        let mut max_allele = 0_i64;
        for s in 0..n_sample_in {
            if !kept(edit, s) {
                continue;
            }
            if masked(edit, s) {
                max_ploidy = max_ploidy.max(2);
                continue;
            }
            let sample = &payload[s * stride..(s + 1) * stride];
            let mut ploidy = 0_usize;
            for j in 0..count {
                let v = int_element(ty, &sample[j * width..(j + 1) * width]);
                if v == vector_end {
                    break;
                }
                ploidy += 1;
                if v > 0 {
                    max_allele = max_allele.max((v >> 1) - 1);
                }
            }
            max_ploidy = max_ploidy.max(ploidy);
        }
        let count_out = max_ploidy.max(1);
        let (ty_out, width_out) = gt_width(max_allele);
        encode_size(&mut self.indiv, count_out, ty_out);
        let (_, vector_end_out) = super::record::int_sentinels(ty_out);
        let vector_end_out = vector_end_out as i32;
        for s in 0..n_sample_in {
            if !kept(edit, s) {
                continue;
            }
            if masked(edit, s) {
                push_int(&mut self.indiv, 0, width_out);
                push_int(&mut self.indiv, 0, width_out);
                for _ in 2..count_out {
                    push_int(&mut self.indiv, vector_end_out, width_out);
                }
                continue;
            }
            let sample = &payload[s * stride..(s + 1) * stride];
            let mut written = 0_usize;
            for j in 0..count {
                let v = int_element(ty, &sample[j * width..(j + 1) * width]);
                if v == vector_end {
                    break;
                }
                push_int(
                    &mut self.indiv,
                    i32::try_from(v).unwrap_or(vector_end_out),
                    width_out,
                );
                written += 1;
            }
            while written < count_out {
                push_int(&mut self.indiv, vector_end_out, width_out);
                written += 1;
            }
        }
        Ok(())
    }
}

/// The output sample count under `edit`.
fn n_samples_out(edit: &EditPlan<'_>, n_in: usize) -> usize {
    edit.keep_samples.map_or(n_in, |k| k.n_selected() as usize)
}

/// Whether sample `s` survives the sample pruning.
fn kept(edit: &EditPlan<'_>, s: usize) -> bool {
    edit.keep_samples
        .is_none_or(|k| k.contains(SampleId(s as u32)))
}

/// Whether a surviving sample's GT becomes `./.`.
fn masked(edit: &EditPlan<'_>, s: usize) -> bool {
    edit.mask_gt.is_some_and(|m| m.get(s))
}

/// The `info_set` replacement marks: on the stack for the common handful of
/// tags, on the heap beyond eight (the `edit_info` pattern).
enum SeenMarks {
    Stack([bool; 8]),
    Heap(Vec<bool>),
}

impl SeenMarks {
    fn new(n: usize) -> Self {
        if n <= 8 {
            Self::Stack([false; 8])
        } else {
            Self::Heap(vec![false; n])
        }
    }

    fn mark(&mut self, i: usize) {
        match self {
            Self::Stack(marks) => marks[i] = true,
            Self::Heap(marks) => marks[i] = true,
        }
    }

    fn get(&self, i: usize) -> bool {
        match self {
            Self::Stack(marks) => marks[i],
            Self::Heap(marks) => marks[i],
        }
    }
}

/// Marker for a sample column shorter than the FORMAT key list.
const ABSENT: (u32, u32) = (u32::MAX, u32::MAX);

/// One recomputed INFO value, encoded per the tag's declared type. The
/// outcomes mirror the text round-trip exactly: an integer value on a Float
/// tag widens, a float on an Integer tag becomes the missing sentinel (the
/// text "0.25" never parsed as an integer), and a Float keeps its original
/// f32 bits — `g6(parse(g6(v))) == g6(v)`, so the read-back text is what the
/// `-Oz` path prints.
fn encode_info_typed(s: &mut Vec<u8>, declared: InfoType, value: &InfoValueOut) {
    match (declared, value) {
        (InfoType::Flag, _) => encode_size(s, 0, value_type::NULL),
        (InfoType::Integer, InfoValueOut::Int(v)) => encode_int_vector(s, &[*v]),
        (InfoType::Integer, InfoValueOut::Float(_)) => {
            encode_size(s, 1, value_type::INT8);
            s.push(i8::MIN as u8);
        }
        (InfoType::Float, InfoValueOut::Int(v)) => {
            encode_size(s, 1, value_type::FLOAT);
            s.extend_from_slice(&(*v as f32).to_le_bytes());
        }
        (InfoType::Float, InfoValueOut::Float(v)) => {
            encode_size(s, 1, value_type::FLOAT);
            s.extend_from_slice(&v.to_le_bytes());
        }
        (InfoType::String, value) => {
            let mut text = Vec::new();
            render_info_value(value, &mut text);
            encode_typed_string(s, &text);
        }
    }
}

// ---- typed value encoders (the mirror of render.rs) ----

/// BCF record-header field widths (hts-specs §6.3; the packed words are
/// `n_allele<<16|n_info` and `n_fmt<<24|n_sample`). htslib refuses an
/// over-large count at parse time (`vcf.c:3912` for INFO, `:4174` for
/// alleles, `vcf_parse_format_dict2` for FORMAT keys); reed writes records
/// it did not parse from text (BCF in, BCF out), so it must check before
/// packing instead of letting the high bits leak into the neighbouring
/// field (doc_v12 D12-9).
const MAX_N_INFO: u32 = 0xFFFF;
/// See [`MAX_N_INFO`].
const MAX_N_ALLELE: u32 = 0xFFFF;
/// See [`MAX_N_INFO`].
const MAX_N_FMT: u32 = 0xFF;
/// See [`MAX_N_INFO`].
const MAX_N_SAMPLE: u32 = 0xFF_FFFF;

/// Packs `n_allele<<16 | n_info`, rejecting counts that do not fit the
/// fields (doc_v12 D12-9).
fn pack_allele_info(n_allele: u32, n_info: u32) -> Result<u32, IoError> {
    if n_allele > MAX_N_ALLELE {
        return Err(bad("BCF record has too many alleles (max 65535)"));
    }
    if n_info > MAX_N_INFO {
        return Err(bad("BCF record has too many INFO fields (max 65535)"));
    }
    Ok(n_allele << 16 | n_info)
}

/// Packs `n_fmt<<24 | n_sample`, rejecting counts that do not fit the
/// fields (doc_v12 D12-9).
fn pack_fmt_sample(n_fmt: u32, n_sample: u32) -> Result<u32, IoError> {
    if n_fmt > MAX_N_FMT {
        return Err(bad("BCF record has too many FORMAT fields (max 255)"));
    }
    if n_sample > MAX_N_SAMPLE {
        return Err(bad("BCF record has too many samples (max 16777215)"));
    }
    Ok(n_fmt << 24 | n_sample)
}

/// A typed string: descriptor + bytes (ID, alleles, INFO strings).
fn encode_typed_string(out: &mut Vec<u8>, value: &[u8]) {
    encode_size(out, value.len(), value_type::CHAR);
    out.extend_from_slice(value);
}

/// A typed integer (dictionary keys, extended counts).
fn encode_typed_int(out: &mut Vec<u8>, value: i64) {
    if (INT8_MIN_VALUE..=INT8_MAX_VALUE).contains(&value) {
        out.push(0x10 | value_type::INT8);
        out.push(value as i8 as u8);
    } else if (INT16_MIN_VALUE..=INT16_MAX_VALUE).contains(&value) {
        out.push(0x10 | value_type::INT16);
        out.extend_from_slice(&(value as i16).to_le_bytes());
    } else {
        out.push(0x10 | value_type::INT32);
        out.extend_from_slice(&(value as i32).to_le_bytes());
    }
}

/// A descriptor: type in the low nibble, count in the high nibble (counts
/// of 15 or more continue in a typed integer).
fn encode_size(out: &mut Vec<u8>, count: usize, ty: u8) {
    if count < 15 {
        out.push(((count as u8) << 4) | ty);
    } else {
        out.push(0xF0 | ty);
        encode_typed_int(out, count as i64);
    }
}

/// The smallest integer type covering `[min, max]` (htslib `bcf_enc_vint`;
/// `i64::MIN` marks a value set that was empty — all missing).
fn int_width(min: i64, max: i64) -> (u8, usize) {
    if min > max || (min >= INT8_MIN_VALUE && max <= INT8_MAX_VALUE) {
        (value_type::INT8, 1)
    } else if min >= INT16_MIN_VALUE && max <= INT16_MAX_VALUE {
        (value_type::INT16, 2)
    } else {
        (value_type::INT32, 4)
    }
}

/// The GT integer width for the largest allele index, judged on the
/// *encoded value* domain: allele `a` encodes to `(a+1)<<1 | phased`, whose
/// maximum is `2a+3`. INT8 ends at allele 62 (`2*62+3 = 127`) and INT16 at
/// 16382 (`2*16382+3 = 32767`); larger alleles need INT32 (htslib judges
/// the encoded values in `bcf_enc_vint`, vcf.c:2955-2993).
///
/// Both GT encoders share this: the thresholds are load-bearing, and one
/// stale copy silently truncates codes with `as i16` (doc_v10 D10-1).
fn gt_width(max_allele: i64) -> (u8, usize) {
    if max_allele <= 62 {
        (value_type::INT8, 1)
    } else if max_allele <= 16382 {
        (value_type::INT16, 2)
    } else {
        (value_type::INT32, 4)
    }
}

/// A typed int vector with htslib's width choice.
///
/// Values outside `[BCF_MIN_BT_INT32, BCF_MAX_BT_INT32]` are written as the
/// missing sentinel instead of being truncated by the final `as i32`
/// (doc_v12 D12-3; [`parse_int_vector`] already filters the text path, this
/// covers the recomputed-value path).
fn encode_int_vector(out: &mut Vec<u8>, values: &[i64]) {
    let (mut min, mut max) = (i64::MAX, i64::MIN);
    for &v in values {
        if v == i64::MIN {
            continue; // the missing sentinel never widens the type
        }
        min = min.min(v);
        max = max.max(v);
    }
    let (ty, width) = int_width(min, max);
    encode_size(out, values.len(), ty);
    let (missing, _) = super::record::int_sentinels(ty);
    for &v in values {
        let v = if v == i64::MIN || !(INT32_MIN_VALUE..=INT32_MAX_VALUE).contains(&v) {
            missing as i32
        } else {
            v as i32
        };
        push_int(out, v, width);
    }
}

/// Appends one integer element of `width` bytes (little-endian, signed).
fn push_int(out: &mut Vec<u8>, value: i32, width: usize) {
    match width {
        1 => out.push(value as i8 as u8),
        2 => out.extend_from_slice(&(value as i16).to_le_bytes()),
        _ => out.extend_from_slice(&value.to_le_bytes()),
    }
}

/// Parses a comma-separated integer vector; `.` and empty parts become the
/// missing sentinel (`i64::MIN`).
///
/// A value outside the BCF-legal range also becomes missing rather than
/// being truncated to 32 bits: htslib's parser checks
/// `overflow || val < BCF_MIN_BT_INT32 || val > BCF_MAX_BT_INT32` and
/// substitutes the missing sentinel (vcf.c:4004-4013), and VCFv4.5 §1.3
/// forbids the reserved band `[-2^31, -2^31+7]` in both VCF and BCF
/// (doc_v12 D12-3). Without the check `3000000000` wrapped to
/// `-1294967296` through the `as i32` in [`encode_int_vector`].
fn parse_int_vector(value: &[u8], out: &mut Vec<i64>) {
    for part in value.split(|&b| b == b',') {
        if part == b"." || part.is_empty() {
            out.push(i64::MIN);
            continue;
        }
        let parsed = parse_i64(part).filter(|v| {
            let ok = (INT32_MIN_VALUE..=INT32_MAX_VALUE).contains(v);
            if !ok {
                warn_once_out_of_range(part);
            }
            ok
        });
        out.push(parsed.unwrap_or(i64::MIN));
    }
}

/// htslib warns once per run about an out-of-range integer and then stays
/// quiet (`extreme_int_warned`, vcf.c:4008-4012).
fn warn_once_out_of_range(value: &[u8]) {
    use std::sync::atomic::{AtomicBool, Ordering};
    static WARNED: AtomicBool = AtomicBool::new(false);
    if !WARNED.swap(true, Ordering::Relaxed) {
        tracing::warn!(
            "Integer INFO value {:?} is outside the BCF range and was set to missing \
             (further occurrences are not reported)",
            String::from_utf8_lossy(value)
        );
    }
}

/// Parses a comma-separated float vector (`.` → the float missing pattern).
fn parse_float_vector(value: &[u8], out: &mut Vec<f32>) {
    for part in value.split(|&b| b == b',') {
        if part == b"." || part.is_empty() {
            out.push(f32::from_bits(FLOAT_MISSING));
        } else {
            out.push(
                std::str::from_utf8(part)
                    .ok()
                    .and_then(|t| t.parse::<f32>().ok())
                    .unwrap_or(f32::from_bits(FLOAT_MISSING)),
            );
        }
    }
}

/// `(ploidy, largest allele index)` of one text genotype.
fn gt_shape(sub: &[u8]) -> (usize, i64) {
    let mut ploidy = 0_usize;
    let mut max_allele = 0_i64;
    for allele in strip_gt_prefix(sub).split(|&b| b == b'/' || b == b'|') {
        ploidy += 1;
        if allele != b"."
            && !allele.is_empty()
            && let Some(v) = parse_i64(allele)
        {
            max_allele = max_allele.max(v);
        }
    }
    (ploidy, max_allele)
}

/// Drops a VCF 4.4 phasing prefix (htslib consumes it only for ≥ 4.4
/// headers; accepted unconditionally here, as in the scanner's decoder).
fn strip_gt_prefix(sub: &[u8]) -> &[u8] {
    match sub.first() {
        Some(b'|' | b'/') => &sub[1..],
        _ => sub,
    }
}

/// One text genotype's BCF codes in `out` (`(allele+1)<<1 | phased`),
/// exactly htslib's parser (vcf.c:3496-3558): allele `i > 0` carries the bit
/// of the separator before it; allele 0's bit is the explicit prefix when
/// present, else inferred (a haploid is phased, else phased iff every later
/// allele is); a missing allele stores the bare phasing bit.
fn gt_codes_into(sub: &[u8], out: &mut Vec<i32>) {
    out.clear();
    let (prefix, mut rest) = match sub.first() {
        Some(b'|') => (Some(true), &sub[1..]),
        Some(b'/') => (Some(false), &sub[1..]),
        _ => (None, sub),
    };
    let mut ploidy = 0_usize;
    let mut anyunphased = false;
    let mut bit = prefix.unwrap_or(false);
    loop {
        match rest.iter().position(|&b| b == b'/' || b == b'|') {
            Some(i) => {
                push_gt_code(out, &rest[..i], bit);
                ploidy += 1;
                bit = rest[i] == b'|';
                anyunphased |= !bit;
                rest = &rest[i + 1..];
            }
            None => {
                push_gt_code(out, rest, bit);
                ploidy += 1;
                break;
            }
        }
    }
    if prefix.is_none() {
        // No explicit phasing for allele 0: infer it the v44 way (a missing
        // haploid allele — code 0 — stays unphased, htslib's `unknown1`).
        let set = if ploidy == 1 {
            out.first().is_some_and(|&c| c != 0)
        } else {
            !anyunphased
        };
        if set && let Some(first) = out.first_mut() {
            *first |= 1;
        }
    }
}

/// One allele token's code: the index plus the phasing bit, or the bare bit
/// for a missing allele (htslib stores `.` as `is_phased`, vcf.c:3519).
fn push_gt_code(out: &mut Vec<i32>, allele: &[u8], phased: bool) {
    let bit = i32::from(phased);
    let code = match parse_i64(allele) {
        // An allele index beyond i32 cannot be represented in any BCF GT
        // encoding; treat it as missing rather than truncating. The bound is
        // the last allele whose `(a+1)<<1 | bit` still fits i32
        // (a = 1073741822; 1073741823 would wrap, doc_v10 D10-1).
        Some(a) if (0..=(i64::from(i32::MAX >> 1) - 1)).contains(&a) => {
            (((a + 1) << 1) as i32) | bit
        }
        _ => bit,
    };
    out.push(code);
}

fn bad(reason: &'static str) -> IoError {
    IoError::MalformedRecord {
        locus: Locus {
            contig_name: "<bcf-encode>".into(),
            pos1: 0,
        },
        reason,
    }
}

/// Strict integer parse (optional leading '-', digits only).
fn parse_i64(text: &[u8]) -> Option<i64> {
    if text.is_empty() {
        return None;
    }
    let (neg, digits) = match text[0] {
        b'-' => (true, &text[1..]),
        _ => (false, text),
    };
    if digits.is_empty() {
        return None;
    }
    let mut value: i64 = 0;
    for &b in digits {
        if !b.is_ascii_digit() {
            return None;
        }
        value = value.checked_mul(10)?.checked_add(i64::from(b - b'0'))?;
    }
    Some(if neg { -value } else { value })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bcf::{BcfEncoder, bcf_header_text, record_ref};
    use crate::vcf::line_edit::{InfoValueOut, edit_vcf_line};
    use reed_core::Bitset;
    use reed_core::sample::SampleMask;

    const HEADER: &[u8] = b"##fileformat=VCFv4.2\n\
        ##FILTER=<ID=PASS,Description=\"All\">\n\
        ##contig=<ID=Chr01,length=1000>\n\
        ##INFO=<ID=AC,Number=A,Type=Integer,Description=\"ac\">\n\
        ##INFO=<ID=AN,Number=1,Type=Integer,Description=\"an\">\n\
        ##INFO=<ID=NS,Number=1,Type=Integer,Description=\"ns\">\n\
        ##INFO=<ID=MAF,Number=1,Type=Float,Description=\"maf\">\n\
        ##INFO=<ID=AF,Number=A,Type=Float,Description=\"af\">\n\
        ##INFO=<ID=QD,Number=1,Type=Float,Description=\"qd\">\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"dp\">\n\
        ##FORMAT=<ID=GP,Number=3,Type=Float,Description=\"gp\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n";

    const LINE: &[u8] = b"Chr01\t100\t.\tA\tG\t50\tPASS\tAC=9;AF=0.5;QD=4.5\tGT:DP:GP\t0/1:15:0.5,0.25,0.25\t1|1:12:0.25,0.5,0.25\t./.:0:.,.,.\t0/0:8:0.25,0.25,0.5\n";

    fn header() -> vcf::Header {
        std::str::from_utf8(HEADER).unwrap().parse().unwrap()
    }

    /// The reference path: edit the text, then noodles encodes the result.
    fn reference_encode(
        line: &[u8],
        plan: &EditPlan<'_>,
        gt_idx: Option<u8>,
        n_out: usize,
    ) -> Vec<u8> {
        let lr = LineRef::parse(line).expect("invariant: line");
        let mut edited = Vec::new();
        edit_vcf_line(&lr, plan, gt_idx, &mut edited).expect("invariant: edit");
        let header_text = std::str::from_utf8(HEADER).unwrap();
        let out_text = match plan.keep_samples {
            None => header_text.to_string(),
            Some(mask) => {
                let mut names = header_text.lines().last().unwrap().split('\t');
                let mut kept: Vec<&str> = (0..9).map(|_| names.next().unwrap()).collect();
                for (i, name) in names.enumerate() {
                    if mask.contains(SampleId(i as u32)) {
                        kept.push(name);
                    }
                }
                let mut h: Vec<&str> = header_text.lines().collect();
                let last = kept.join("\t");
                let chrom_line = h.len() - 1;
                h[chrom_line] = &last;
                h.join("\n") + "\n"
            }
        };
        let out_header: vcf::Header = out_text.parse().expect("invariant: header");
        assert_eq!(out_header.sample_names().len(), n_out);
        let mut enc = BcfEncoder::new().expect("invariant: noodles");
        enc.write_header(&out_header).expect("invariant: header");
        let mut out = Vec::new();
        enc.encode_record(&out_header, &edited, &mut out)
            .expect("invariant: noodles encode");
        out
    }

    fn mask_of(n: usize, keeps: &[usize]) -> SampleMask {
        SampleMask::from_ids(n, keeps.iter().map(|&i| SampleId(i as u32))).unwrap()
    }

    #[test]
    fn encode_text_applies_the_edit_like_the_text_path() {
        let mut bits = Bitset::new(4);
        bits.set(1, true);
        let keep = mask_of(4, &[0, 2, 3]);
        let plan = EditPlan {
            keep_samples: Some(&keep),
            mask_gt: Some(&bits),
            info_remove: &[b"AF"],
            info_set: &[
                (b"AC".as_slice(), InfoValueOut::Int(3)),
                (b"NS".as_slice(), InfoValueOut::Int(3)),
                (b"MAF".as_slice(), InfoValueOut::Float(0.25)),
            ],
        };
        let lr = LineRef::parse(LINE).unwrap();
        let site = SiteView::parse(&lr).unwrap();
        // The encoder expects the *output* sample count: build it from the
        // pruned header (the production path does the same).
        let pruned: vcf::Header = std::str::from_utf8(HEADER)
            .unwrap()
            .replace("\ts4\n", "\n")
            .parse()
            .expect("invariant: pruned header");
        let mut enc = NativeEncoder::new(&pruned).unwrap();
        let mut got = Vec::new();
        enc.encode_text(&lr, &site, 4, &plan, &mut got)
            .expect("invariant: encode");
        let want = reference_encode(LINE, &plan, Some(0), 3);
        assert_eq!(got, want, "encode_text vs edit+noodles");
    }

    #[test]
    fn encode_bcf_matches_the_text_reference() {
        let pruned: vcf::Header = std::str::from_utf8(HEADER)
            .unwrap()
            .replace("\ts4\n", "\n")
            .parse()
            .expect("invariant: pruned header");
        let header = header();
        let dict_text = bcf_header_text(&header).expect("invariant: bcf header");
        let dict = BcfDict::from_header(&dict_text).expect("invariant: dict");
        // The framed source record (verbatim encode of LINE).
        let mut framed = Vec::new();
        {
            let mut enc = NativeEncoder::new(&header).unwrap();
            let lr = LineRef::parse(LINE).unwrap();
            let site = SiteView::parse(&lr).unwrap();
            enc.encode_text(&lr, &site, 4, &EditPlan::default(), &mut framed)
                .unwrap();
        }
        let mut bits = Bitset::new(4);
        bits.set(1, true);
        let keep = mask_of(4, &[0, 2, 3]);
        let plan = EditPlan {
            keep_samples: Some(&keep),
            mask_gt: Some(&bits),
            info_remove: &[b"AF"],
            info_set: &[
                (b"AC".as_slice(), InfoValueOut::Int(3)),
                (b"NS".as_slice(), InfoValueOut::Int(3)),
                (b"MAF".as_slice(), InfoValueOut::Float(0.25)),
            ],
        };
        let rec = record_ref(&framed).expect("invariant: framed");
        let mut enc = NativeEncoder::new(&pruned).unwrap();
        let mut got = Vec::new();
        enc.encode_bcf(&rec, &dict, &plan, &mut got)
            .expect("invariant: encode_bcf");
        let want = reference_encode(LINE, &plan, Some(0), 3);
        assert_eq!(
            got, want,
            "encode_bcf vs edit+noodles\n got: {got:02x?}\nwant: {want:02x?}"
        );
    }

    #[test]
    fn gt_codes_follow_the_htslib_parser() {
        let cases: &[(&[u8], &[i32])] = &[
            (b"0/1", &[2, 4]),
            (b"0|1", &[3, 5]),
            (b"|0/1", &[3, 4]), // 4.4 prefix: allele 0 explicitly phased
            (b"/0|1", &[2, 5]), // explicit unphased first allele
            (b"1", &[5]),       // haploid: implicitly phased
            (b".", &[0]),       // missing haploid stays unphased
            (b"|.", &[1]),      // prefixed missing: the bare phasing bit
            (b"0|.", &[3, 1]),  // all phased: allele 0 inferred; missing carries the bit
            (b"0/.", &[2, 0]),
            (b"./.", &[0, 0]),
            (b".|.", &[1, 1]), // all phased -> first allele inferred phased
        ];
        let mut codes = Vec::new();
        for (text, want) in cases {
            gt_codes_into(text, &mut codes);
            assert_eq!(codes, *want, "GT {}", String::from_utf8_lossy(text));
        }
    }

    #[test]
    fn gt_shape_counts_alleles_after_the_prefix() {
        assert_eq!(gt_shape(b"|0/1"), (2, 1));
        assert_eq!(gt_shape(b"0/1"), (2, 1));
        assert_eq!(gt_shape(b"|1"), (1, 1));
        assert_eq!(gt_shape(b"."), (1, 0));
    }

    /// doc_v10 D10-1: the GT integer width follows the *encoded value* domain
    /// `(allele+1)<<1 | phased` — whose maximum for allele `a` is `2a+3` —
    /// not the allele domain. So INT8 ends at allele 62 (`2*62+3 = 127`) and
    /// INT16 at 16382 (`2*16382+3 = 32767`); allele 16383 encodes to 32768
    /// and must take INT32. The old allele-domain threshold (`<= 32766`
    /// alleles → INT16) silently truncated those codes with `as i16`
    /// (htslib judges the encoded values: `bcf_enc_vint`, vcf.c:2955-2993,
    /// `BCF_MAX_BT_INT16 = 0x7fff`, vcf.h:1539).
    #[test]
    fn gt_width_follows_the_encoded_value_domain() {
        const GT_HEADER: &[u8] = b"##fileformat=VCFv4.2\n\
            ##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
            ##contig=<ID=Chr01,length=1000>\n\
            ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n";
        let header: vcf::Header = std::str::from_utf8(GT_HEADER)
            .expect("invariant: utf-8")
            .parse()
            .expect("invariant: header");
        let mut enc = NativeEncoder::new(&header).expect("invariant: encoder");
        // (allele, expected descriptor type, expected first code)
        let cases: &[(u64, u8, i64)] = &[
            (62, value_type::INT8, 126),
            (63, value_type::INT16, 128),
            (16382, value_type::INT16, 32766),
            (16383, value_type::INT32, 32768),
            (32766, value_type::INT32, 65534),
            (32767, value_type::INT32, 65536),
        ];
        for &(allele, want_ty, want_code) in cases {
            let line = format!("Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT\t{allele}/{allele}\n");
            let lr = LineRef::parse(line.as_bytes()).expect("invariant: 10 columns");
            let site = SiteView::parse(&lr).expect("invariant: site");
            let mut out = Vec::new();
            enc.encode_text(&lr, &site, 1, &EditPlan::default(), &mut out)
                .expect("invariant: encode");
            let rec = record_ref(&out).expect("invariant: framed record");
            let mut c = Cursor::new(rec.indiv);
            let _key = c.typed_int().expect("invariant: GT key");
            let (ty, count) = c.typed().expect("invariant: GT descriptor");
            assert_eq!(count, 2, "allele {allele}");
            assert_eq!(ty, want_ty, "allele {allele}: descriptor type");
            assert_eq!(
                c.int_value(ty).expect("invariant: code"),
                want_code,
                "allele {allele}: encoded code"
            );
        }
    }
}
