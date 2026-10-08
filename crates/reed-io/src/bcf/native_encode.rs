//! Native text → BCF encoder (doc_v3 02-bcf.md B4).
//!
//! The noodles encoder costs ~800 µs per 348-sample record (it re-parses
//! every formatted value), which breaks the `-Ob` budget. This encoder
//! walks the record's text once — the mirror image of [`super::render`] —
//! and writes the BCF bytes directly, with no per-value allocation.
//!
//! Widths follow htslib's `bcf_enc_vint` (the smallest of int8/int16/int32
//! covering every element, legal ranges −120…127 / −32760…32767), and
//! dictionary indices come from the output header (see
//! [`super::encode::bcf_header_text`]).

use std::collections::HashMap;

use crate::error::{IoError, Locus};
use crate::noodles_vcf as vcf;
use crate::vcf::line::LineRef;
use crate::vcf::site::SiteView;

/// Marker for a sample column shorter than the FORMAT key list.
const ABSENT: (u32, u32) = (u32::MAX, u32::MAX);

use super::record::{
    value_type, FLOAT_MISSING, FLOAT_VECTOR_END, INT16_MAX_VALUE, INT16_MIN_VALUE,
    INT32_MAX_VALUE, INT32_MIN_VALUE, INT8_MAX_VALUE, INT8_MIN_VALUE,
};

/// A per-worker native encoder.
pub struct NativeEncoder {
    contig_index: HashMap<Vec<u8>, i32>,
    /// Tag name → dictionary index (PASS included).
    tags: HashMap<Vec<u8>, u32>,
    /// INFO tag → declared value shape.
    info_types: HashMap<Vec<u8>, InfoType>,
    /// FORMAT tag → declared BCF type.
    format_types: HashMap<Vec<u8>, u8>,
    n_samples: usize,
    /// Scratch buffers, reused across records.
    shared: Vec<u8>,
    indiv: Vec<u8>,
    /// Byte ranges of the current record's sample subfields, flattened as
    /// `sample * n_keys + key` (one pass over the sample text per record).
    cells: Vec<(u32, u32)>,
    ints: Vec<i64>,
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
                tags.insert(key.as_bytes().to_vec(), u32::try_from(idx).unwrap_or(u32::MAX));
            }
        }
        let mut contig_index = HashMap::new();
        for name in header.contigs().keys() {
            if let Some(idx) = maps.contigs().get_index_of(name) {
                contig_index.insert(name.as_bytes().to_vec(), i32::try_from(idx).unwrap_or(i32::MAX));
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
        })
    }

    /// Number of sample columns the encoder expects.
    #[must_use]
    pub const fn n_samples(&self) -> usize {
        self.n_samples
    }

    /// Encodes one VCF text record into `out` (the BCF record bytes,
    /// including the two length prefixes).
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] for malformed input and
    /// [`IoError::Header`] for tags the output header does not declare.
    pub fn encode(&mut self, line: &[u8], out: &mut Vec<u8>) -> Result<(), IoError> {
        let Some(lr) = LineRef::parse(line) else {
            return Err(bad("record line has fewer than 9 columns"));
        };
        let site = SiteView::parse(&lr)?;
        self.shared.clear();
        self.indiv.clear();
        self.encode_shared(&lr, &site)?;
        self.encode_indiv(&lr)?;
        out.extend_from_slice(&(self.shared.len() as u32).to_le_bytes());
        out.extend_from_slice(&(self.indiv.len() as u32).to_le_bytes());
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

    fn encode_shared(&mut self, lr: &LineRef<'_>, site: &SiteView<'_>) -> Result<(), IoError> {
        let chrom = lr.chrom();
        let contig = self
            .contig_index
            .get(chrom)
            .copied()
            .ok_or_else(|| Self::undeclared("contig", chrom))?;
        let pos0 = i32::try_from(site.pos1.saturating_sub(1)).unwrap_or(i32::MAX);
        let rlen = i32::try_from(site.reference.len()).unwrap_or(1).max(1);
        let s = &mut self.shared;
        s.extend_from_slice(&contig.to_le_bytes());
        s.extend_from_slice(&pos0.to_le_bytes());
        s.extend_from_slice(&rlen.to_le_bytes());
        s.extend_from_slice(
            &site.qual.map_or(FLOAT_MISSING, f32::to_bits).to_le_bytes(),
        );
        s.extend_from_slice(&[0; 8]); // n_allele|n_info, n_fmt|n_sample

        // ID: an empty ID (".") is stored as a zero-length string.
        let id = lr.id();
        encode_typed_string(s, if id == b"." { &[] } else { id });

        // Alleles: REF then the ALT list ("." means no ALT).
        let alts = lr.alts();
        let n_allele = if alts == b"." {
            encode_typed_string(s, lr.reference());
            1_u32
        } else {
            encode_typed_string(s, lr.reference());
            let mut n = 1_u32;
            for alt in alts.split(|&b| b == b',') {
                encode_typed_string(s, alt);
                n += 1;
            }
            n
        };

        // FILTER.
        let filter = lr.filter_text();
        let n_flt = if filter == b"." || filter.is_empty() {
            encode_size(s, 0, value_type::INT8);
            0_u32
        } else {
            let mut indices = Vec::new();
            for name in filter.split(|&b| b == b';') {
                let idx = self
                    .tags
                    .get(name)
                    .copied()
                    .ok_or_else(|| Self::undeclared("FILTER", name))?;
                indices.push(i64::from(idx));
            }
            encode_int_vector(s, &indices);
            u32::try_from(indices.len()).unwrap_or(u32::MAX)
        };

        // INFO: resolve the dictionary keys and declared types first (the
        // writes below hold a mutable borrow of the shared buffer).
        let info = lr.info_text();
        let mut info_items: Vec<(u32, InfoType, Option<&[u8]>)> = Vec::new();
        if info != b"." && !info.is_empty() {
            for entry in info.split(|&b| b == b';') {
                if entry.is_empty() {
                    continue;
                }
                let (key, value) = match entry.iter().position(|&b| b == b'=') {
                    Some(i) => (&entry[..i], Some(&entry[i + 1..])),
                    None => (entry, None),
                };
                let idx = self
                    .tags
                    .get(key)
                    .copied()
                    .ok_or_else(|| Self::undeclared("INFO", key))?;
                let declared = self.info_types.get(key).copied().unwrap_or(InfoType::String);
                info_items.push((idx, declared, value));
            }
        }
        let n_info = if info_items.is_empty() {
            0_u32
        } else {
            let mut n = 0_u32;
            for (idx, declared, value) in info_items {
                encode_typed_int(s, i64::from(idx));
                match (declared, value) {
                    // A flag carries no value at all.
                    (InfoType::Flag, _) => encode_size(s, 0, value_type::NULL),
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
                    (InfoType::String, Some(v)) | (_, Some(v)) => {
                        encode_typed_string(s, if v == b"." { &[] } else { v });
                    }
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
                n += 1;
            }
            n
        };
        // Patch the packed counters (n_allele<<16|n_info, n_fmt<<24|n_sample;
        // n_fmt is filled in by `encode_indiv`).
        s[16..20].copy_from_slice(&(n_allele << 16 | n_info).to_le_bytes());
        let _ = n_flt;
        Ok(())
    }

    fn encode_indiv(&mut self, lr: &LineRef<'_>) -> Result<(), IoError> {
        let format = lr.format_text();
        let n_fmt = if format == b"." || format.is_empty() {
            0_u32
        } else {
            u32::try_from(format.iter().filter(|&&b| b == b':').count() + 1).unwrap_or(u32::MAX)
        };
        if n_fmt > 0 {
            let samples = lr.samples_text();
            let keys: Vec<&[u8]> = format.split(|&b| b == b':').collect();
            let n_keys = keys.len();
            // One pass over the sample text: every subfield's byte range.
            self.cells.clear();
            self.cells.resize(self.n_samples * n_keys, ABSENT);
            let mut column_start = 0_usize;
            for sample in 0..self.n_samples {
                let rest = &samples[column_start..];
                let column_len = rest
                    .iter()
                    .position(|&b| b == b'\t')
                    .unwrap_or(rest.len());
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
                let ty = self.format_types.get(*key).copied().unwrap_or(value_type::CHAR);
                encode_typed_int(&mut self.indiv, i64::from(idx));
                if *key == b"GT" {
                    self.encode_gt_column(samples, k, n_keys);
                } else {
                    match ty {
                        value_type::FLOAT => self.encode_float_column(samples, k, n_keys),
                        value_type::CHAR => self.encode_string_column(samples, k, n_keys),
                        _ => self.encode_int_column(samples, k, n_keys)?,
                    }
                }
            }
        }
        let s = &mut self.shared;
        s[20..24].copy_from_slice(&(n_fmt << 24 | self.n_samples as u32).to_le_bytes());
        Ok(())
    }

    /// One sample's subfield `k` (`ABSENT` when the column is shorter).
    fn cell<'a>(cells: &[(u32, u32)], samples: &'a [u8], sample: usize, k: usize, n_keys: usize) -> &'a [u8] {
        match cells.get(sample * n_keys + k).copied() {
            Some((start, end)) if start != u32::MAX => &samples[start as usize..end as usize],
            _ => b".",
        }
    }

    fn encode_gt_column(&mut self, samples: &[u8], k: usize, n_keys: usize) {
        let mut max_ploidy = 0_usize;
        let mut max_allele = 0_i64;
        for s in 0..self.n_samples {
            let (ploidy, allele) = gt_shape(Self::cell(&self.cells, samples, s, k, n_keys));
            max_ploidy = max_ploidy.max(ploidy);
            max_allele = max_allele.max(allele);
        }
        let max_ploidy = max_ploidy.max(1);
        let (ty, width) = if max_allele <= 62 {
            (value_type::INT8, 1_usize)
        } else if max_allele <= 32766 {
            (value_type::INT16, 2)
        } else {
            (value_type::INT32, 4)
        };
        encode_size(&mut self.indiv, max_ploidy, ty);
        let (_, vector_end) = super::record::int_sentinels(ty);
        let vector_end = vector_end as i32;
        for s in 0..self.n_samples {
            let codes = gt_codes(Self::cell(&self.cells, samples, s, k, n_keys));
            for i in 0..max_ploidy {
                push_int(
                    &mut self.indiv,
                    codes.get(i).copied().unwrap_or(vector_end),
                    width,
                );
            }
        }
    }

    fn encode_int_column(&mut self, samples: &[u8], k: usize, n_keys: usize) -> Result<(), IoError> {
        let mut max_len = 0_usize;
        let (mut min, mut max) = (i64::MAX, i64::MIN);
        for s in 0..self.n_samples {
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
        for s in 0..self.n_samples {
            let sub = Self::cell(&self.cells, samples, s, k, n_keys);
            let mut written = 0_usize;
            if sub != b"." && !sub.is_empty() {
                for part in sub.split(|&b| b == b',') {
                    let v = if part == b"." || part.is_empty() {
                        missing
                    } else {
                        let v = parse_i64(part).ok_or_else(|| bad("non-integer FORMAT value"))?;
                        i32::try_from(v.clamp(INT32_MIN_VALUE, INT32_MAX_VALUE))
                            .unwrap_or(missing)
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

    fn encode_float_column(&mut self, samples: &[u8], k: usize, n_keys: usize) {
        let mut max_len = 0_usize;
        for s in 0..self.n_samples {
            let sub = Self::cell(&self.cells, samples, s, k, n_keys);
            if sub == b"." || sub.is_empty() {
                max_len = max_len.max(1);
            } else {
                max_len = max_len.max(sub.split(|&b| b == b',').count());
            }
        }
        let max_len = max_len.max(1);
        encode_size(&mut self.indiv, max_len, value_type::FLOAT);
        for s in 0..self.n_samples {
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
                let bits = if written == 0 { FLOAT_MISSING } else { FLOAT_VECTOR_END };
                self.indiv.extend_from_slice(&bits.to_le_bytes());
                written += 1;
            }
        }
    }

    fn encode_string_column(&mut self, samples: &[u8], k: usize, n_keys: usize) {
        // The width is the longest value; a missing value is the literal `.`.
        let mut width = 1_usize;
        for s in 0..self.n_samples {
            let sub = Self::cell(&self.cells, samples, s, k, n_keys);
            width = width.max(if sub.is_empty() { 1 } else { sub.len() });
        }
        encode_size(&mut self.indiv, width, value_type::CHAR);
        for s in 0..self.n_samples {
            let sub = Self::cell(&self.cells, samples, s, k, n_keys);
            let value: &[u8] = if sub.is_empty() { b"." } else { sub };
            self.indiv.extend_from_slice(value);
            for _ in value.len()..width {
                self.indiv.push(0);
            }
        }
    }
}

// ---- typed value encoders (the mirror of render.rs) ----

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

/// A typed int vector with htslib's width choice.
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
        let v = if v == i64::MIN { missing as i32 } else { v as i32 };
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
fn parse_int_vector(value: &[u8], out: &mut Vec<i64>) {
    for part in value.split(|&b| b == b',') {
        if part == b"." || part.is_empty() {
            out.push(i64::MIN);
        } else {
            out.push(parse_i64(part).unwrap_or(i64::MIN));
        }
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
    for allele in sub.split(|&b| b == b'/' || b == b'|') {
        ploidy += 1;
        if allele != b"." && !allele.is_empty()
            && let Some(v) = parse_i64(allele)
        {
            max_allele = max_allele.max(v);
        }
    }
    (ploidy, max_allele)
}

/// One text genotype's BCF codes (`(allele+1)<<1 | phased`), following
/// htslib's parser: allele `i > 0` carries the bit of the separator before
/// it, and allele 0 is phased when a haploid or when every later allele is
/// (htslib's VCF 4.4 "prefixed phasing" normalisation, `vcf.c`).
fn gt_codes(sub: &[u8]) -> Vec<i32> {
    let mut alleles: Vec<&[u8]> = Vec::with_capacity(2);
    let mut phased_seps: Vec<bool> = Vec::with_capacity(2);
    let mut rest = sub;
    loop {
        match rest.iter().position(|&b| b == b'/' || b == b'|') {
            Some(i) => {
                alleles.push(&rest[..i]);
                phased_seps.push(rest[i] == b'|');
                rest = &rest[i + 1..];
            }
            None => {
                alleles.push(rest);
                break;
            }
        }
    }
    let ploidy = alleles.len();
    let all_phased = phased_seps.iter().all(|&s| s);
    let mut codes = Vec::with_capacity(ploidy);
    for (i, allele) in alleles.iter().enumerate() {
        let phased = if i == 0 {
            ploidy == 1 || all_phased
        } else {
            phased_seps[i - 1]
        };
        let code = match parse_i64(allele) {
            Some(a) => (((a + 1) << 1) as i32) | i32::from(phased),
            // htslib stores a missing allele as the bare `bcf_gt_missing`
            // sentinel 0 (it cannot carry a phasing bit).
            None => 0,
        };
        codes.push(code);
    }
    codes
}

fn bad(reason: &'static str) -> IoError {
    IoError::MalformedRecord {
        locus: Locus { contig_name: "<bcf-encode>".into(), pos1: 0 },
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
