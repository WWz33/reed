//! BCF GT decoding without the text round-trip (doc_v3 02-bcf.md §2.4).
//!
//! The int8 two-element case — what every biallelic diploid file uses — goes
//! through one 128 KiB table indexed by the two raw bytes. Each entry is
//! built by rendering the pair with the renderer's `fmt_gt` and decoding the
//! result with the text scanner's decoder, so the fast path *is* the text
//! path by construction: the same half-call, ploidy and multiallelic
//! policies, the same counters, the same results as `bcftools view` +
//! rescan. Unusual encodings (int16/int32 GT, other element counts, empty
//! vectors) take the same render-then-decode route per record.
//!
//! A record's shared block past the header is *skipped*, not decoded:
//! `stats` needs only contig/position/allele count, so ID, alleles, FILTER
//! and INFO cost nothing (the §2.5 advantage over rendering them to text).

use reed_core::genotype::{GenotypePolicies, GtCode, PolicyCounters};
use reed_core::sample::KeepMask;

use super::dict::BcfDict;
use super::record::{BcfRecordRef, Cursor, FormatTable, malformed, type_width, value_type};
use super::render::fmt_gt;
use crate::error::{IoError, Locus};
use crate::vcf::genotype_scan::{LocusSource, decode_gt_text};

/// The shared-block header fields a numeric consumer needs (hts-specs
/// §6.3.1): everything before the ID string.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SiteHeader {
    /// Contig dictionary index.
    pub chrom: i32,
    /// 0-based position.
    pub pos0: i32,
    /// Allele count, REF included.
    pub n_allele: usize,
    /// Samples declared by the record.
    pub n_sample: usize,
    /// FORMAT fields declared by the record.
    pub n_fmt: usize,
}

/// Reads a record's header fields; the rest of the shared block is skipped.
///
/// # Errors
///
/// [`IoError::MalformedRecord`] when the shared block is shorter than its
/// header, and when the record declares no allele (doc_v26 V26-05):
/// `n_allele` counts REF plus the ALTs, so zero means the shared block has no
/// allele string at all. reed used to accept it — the renderer printed `.`
/// for REF and ALT and the GT path treated every code as "0 alleles" — while
/// htslib rejects such a record; a BCF that no other tool accepts must not
/// flow through filter/stats looking valid.
pub fn site_header(record: &BcfRecordRef<'_>) -> Result<SiteHeader, IoError> {
    let (header, _) = super::record::SharedHeader::parse(record.shared)?;
    Ok(SiteHeader {
        chrom: header.chrom,
        pos0: header.pos0,
        n_allele: header.n_allele,
        n_sample: header.n_sample,
        n_fmt: header.n_fmt,
    })
}

/// The entry's policies reject this pair; the decoder re-runs the text path
/// to recover the exact error. Trigger bits 0..`FIELDS.len()` come from
/// [`PolicyCounters::FIELDS`] (doc_v25 V25-07c), so this is the only bit
/// defined here.
const TRIG_ERROR: u8 = 1 << 7;

/// One table entry.
#[derive(Clone, Copy)]
struct GtEntry {
    code: GtCode,
    triggers: u8,
}

/// The GT fast path: the byte-pair table plus the policies it was built
/// under. Build once per run and share it across workers.
pub struct GtTable {
    entries: Box<[GtEntry]>,
    /// The same table rendered under the VCF 4.4 prefixed-phasing rules.
    /// The first allele's phasing bit is otherwise inferred away, so a
    /// header-aware decode needs its own entries (doc_v28 V28-07).
    entries_vcf44: Box<[GtEntry]>,
    policies: GenotypePolicies,
}

impl GtTable {
    /// Builds the table for `policies` by rendering each int8 byte pair and
    /// running the text decoder on it.
    #[must_use]
    pub fn new(policies: GenotypePolicies) -> Self {
        Self {
            entries: Self::build_entries(&policies, false),
            entries_vcf44: Self::build_entries(&policies, true),
            policies,
        }
    }

    fn build_entries(policies: &GenotypePolicies, vcf44: bool) -> Box<[GtEntry]> {
        let locus = Locus {
            contig_name: "<bcf-gt-table>".into(),
            pos1: 0,
        };
        let mut entries = vec![
            GtEntry {
                code: GtCode::MISSING,
                triggers: TRIG_ERROR
            };
            65_536
        ];
        let mut text = Vec::with_capacity(8);
        for pair in 0..65_536_u32 {
            // The table index is the little-endian byte pair as it appears in
            // the record, so `decode_pair` can index with the two raw bytes.
            let raw = u16::try_from(pair)
                .expect("invariant: the loop is bounded by 65,536")
                .to_le_bytes();
            text.clear();
            // INT8 is an integer type, so rendering succeeds except for a
            // negative GT element. Those pairs stay at `TRIG_ERROR` and the
            // lookup re-runs the slow path for the exact error (V28-04).
            let Ok(()) = fmt_gt(&mut text, value_type::INT8, raw.len(), &raw, vcf44) else {
                continue;
            };
            let mut counters = PolicyCounters::default();
            let mut token = [None; 2];
            // The table spans every allele an int8 pair can name (0..=127);
            // each record's own bound is applied at lookup (doc_v25 V25-01).
            if let Ok(code) =
                decode_gt_text(&text, &locus, 127, policies, &mut counters, &mut token)
            {
                entries[pair as usize] = GtEntry {
                    code,
                    triggers: pack_triggers(&counters),
                };
            }
        }
        entries.into_boxed_slice()
    }

    /// The policies the table (and its slow path) decode under.
    #[must_use]
    pub const fn policies(&self) -> &GenotypePolicies {
        &self.policies
    }

    /// Decodes one int8 two-element genotype.
    ///
    /// `n_allele` is the record's allele count (REF included, from the
    /// shared-block header). The stored values are the BCF GT codes: `0`/`1`
    /// are the missing allele (bare phasing bit), `>= 2` is
    /// `(allele + 1) << 1 | phased` (htslib's `bcf_format_gt_v2`). An encoded
    /// allele at or above `n_allele` is malformed input, refused before the
    /// table lookup (doc_v25 V25-01).
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] for an allele index beyond the record's
    /// alleles; the entries' policies reject this genotype (`Reject`-style
    /// ploidy or multiallelic policies).
    pub fn decode_pair<S: LocusSource>(
        &self,
        b0: u8,
        b1: u8,
        locus: &S,
        n_allele: usize,
        counters: &mut PolicyCounters,
    ) -> Result<GtCode, IoError> {
        self.decode_pair_with(b0, b1, locus, n_allele, counters, false)
    }

    fn decode_pair_with<S: LocusSource>(
        &self,
        b0: u8,
        b1: u8,
        locus: &S,
        n_allele: usize,
        counters: &mut PolicyCounters,
        vcf44: bool,
    ) -> Result<GtCode, IoError> {
        for b in [b0, b1] {
            let v = i64::from(b as i8);
            if v < 0 {
                // The vector-end sentinel is legal and handled by the table;
                // every other negative code is corrupt (V28-04).
                let (_, vector_end) = super::record::int_sentinels(value_type::INT8);
                if v == vector_end {
                    // Vector-end ends the genotype: later elements are
                    // padding and must be ignored, not validated.
                    break;
                }
                return Err(IoError::MalformedRecord {
                    locus: locus.error_locus(),
                    reason: "BCF GT element is negative".to_string(),
                });
            }
            if v < 2 {
                // 0/1: missing allele (bare phasing bit).
                continue;
            }
            let allele = usize::from((b >> 1) - 1);
            if allele >= n_allele {
                return Err(IoError::MalformedRecord {
                    locus: locus.error_locus(),
                    reason: format!(
                        "GT allele index {allele} exceeds ALT count {}",
                        n_allele.saturating_sub(1)
                    ),
                });
            }
        }
        let entries = if vcf44 {
            &self.entries_vcf44
        } else {
            &self.entries
        };
        let entry = entries[usize::from(b0) | (usize::from(b1) << 8)];
        if entry.triggers == TRIG_ERROR {
            // Recover the exact error (and locus) through the text path.
            // The pair's bytes were domain-checked above, so the slow path
            // can run unbounded (its errors here are syntax or policy).
            let mut text = Vec::new();
            return self.decode_slow(
                value_type::INT8,
                &[b0, b1],
                locus,
                u32::MAX,
                counters,
                &mut text,
                vcf44,
            );
        }
        apply_triggers(entry.triggers, counters);
        Ok(entry.code)
    }

    /// Decodes one sample's GT field of any integer encoding through the
    /// text path (render, then parse). `n_alt` bounds the allele domain
    /// (`u32::MAX` = unbounded; doc_v25 V25-01).
    #[allow(clippy::too_many_arguments)]
    fn decode_slow<S: LocusSource>(
        &self,
        ty: u8,
        data: &[u8],
        locus: &S,
        n_alt: u32,
        counters: &mut PolicyCounters,
        text: &mut Vec<u8>,
        vcf44: bool,
    ) -> Result<GtCode, IoError> {
        if !matches!(ty, value_type::INT8 | value_type::INT16 | value_type::INT32) {
            return Err(IoError::MalformedRecord {
                locus: locus.error_locus(),
                reason: "BCF GT field is not integer-typed".to_string(),
            });
        }
        let count = data.len() / type_width(ty);
        text.clear();
        // `ty` was just checked above, so the renderer cannot reject it.
        fmt_gt(text, ty, count, data, vcf44).map_err(|reason| IoError::MalformedRecord {
            locus: locus.error_locus(),
            reason: reason.to_string(),
        })?;
        let mut token = [None; 2];
        decode_gt_text(text, locus, n_alt, &self.policies, counters, &mut token)
    }
}

/// One bit per counter the entry fires (each entry is one genotype, so every
/// field is 0 or 1).
fn pack_triggers(c: &PolicyCounters) -> u8 {
    let mut t = 0_u8;
    for (bit, field) in PolicyCounters::FIELDS.iter().enumerate() {
        if (field.get)(c) > 0 {
            t |= 1 << bit;
        }
    }
    t
}

fn apply_triggers(t: u8, c: &mut PolicyCounters) {
    for (bit, field) in PolicyCounters::FIELDS.iter().enumerate() {
        if t & (1 << bit) != 0 {
            *(field.get_mut)(c) += 1;
        }
    }
}

/// What the fast path needs to decode one record's GT column: the FORMAT
/// index of `GT` (resolved once per run), the expected sample count, and the
/// locus pieces for error reporting.
#[derive(Clone, Copy, Debug)]
pub struct GtQuery<'a> {
    /// The FORMAT dictionary index of `GT`, when the header declares it.
    pub gt_key: Option<u32>,
    /// Samples the header declares.
    pub n_samples: usize,
    /// Contig name for the error locus.
    pub contig_name: &'a [u8],
    /// 1-based position for the error locus.
    pub pos1: u64,
    /// The header is VCF ≥ 4.4: apply the prefixed first-allele phasing
    /// rules when decoding (`BcfDict::vcf44`).
    pub vcf44: bool,
}

/// The FORMAT entry of `GT` for one record, without applying any genotype
/// policy: `(type, count, payload offset)`. `None` when the record carries no
/// GT field.
///
/// This is the raw view a consumer needs when it does not want the policies
/// applied — `pca`/`distance` accept a haploid genotype as missing data
/// rather than rejecting it, which is the text path's rule (doc_v28 V28-15).
///
/// # Errors
///
/// [`IoError::MalformedRecord`] for a truncated or oversized individual
/// block; the `formats` table, when given, is rebuilt for `record`.
pub fn gt_field(
    record: &BcfRecordRef<'_>,
    header: &SiteHeader,
    gt_key: Option<u32>,
    formats: Option<&mut FormatTable>,
) -> Result<Option<(u8, usize, usize)>, IoError> {
    let Some(gt_key) = gt_key.filter(|_| header.n_fmt > 0) else {
        return Ok(None);
    };
    if let Some(table) = formats {
        table.rebuild(record, header.n_fmt, header.n_sample)?;
        return Ok(table
            .fields()
            .iter()
            .find(|(key, ..)| *key == gt_key)
            .map(|&(_, ty, count, start)| (ty, count, start)));
    }
    let indiv = record.indiv;
    let mut ci = Cursor::new(indiv);
    for _ in 0..header.n_fmt {
        let key = ci.typed_int()?;
        let key = u32::try_from(key).map_err(|_| malformed("BCF FORMAT key is negative"))?;
        let (ty, count) = ci.typed()?;
        let start = indiv.len() - ci.remaining();
        let bytes = type_width(ty)
            .checked_mul(count)
            .and_then(|v| v.checked_mul(header.n_sample))
            .ok_or_else(|| malformed("BCF FORMAT field size overflows"))?;
        ci.take(bytes)?;
        if key == gt_key {
            return Ok(Some((ty, count, start)));
        }
    }
    Ok(None)
}

/// One sample's raw GT allele indices from a decoded `GT` field: `None` for a
/// missing element, and the vector-end sentinel truncates the genotype (the
/// remaining elements are padding).
///
/// Values are the BCF codes, not text: `0`/`1` are the missing allele with
/// the phasing bit, and `>= 2` is `(allele + 1) << 1 | phased`. The function
/// returns allele indices only — no policy, no counters — so a consumer with
/// its own genotype rule (`pca`) reads the same bytes as the renderer.
///
/// # Errors
///
/// [`IoError::MalformedRecord`] for a non-integer GT type, an allele index
/// beyond `n_allele`, or a corrupt negative element (doc_v28 V28-04).
pub fn gt_raw_alleles(
    ty: u8,
    count: usize,
    data: &[u8],
    n_allele: usize,
    contig_name: &[u8],
    pos1: u64,
) -> Result<([Option<u16>; 2], usize), IoError> {
    let malformed_at = |reason: &str| IoError::MalformedRecord {
        locus: Locus {
            contig_name: String::from_utf8_lossy(contig_name).into_owned(),
            pos1,
        },
        reason: reason.to_string(),
    };
    if !matches!(ty, value_type::INT8 | value_type::INT16 | value_type::INT32) {
        return Err(malformed_at("BCF GT field is not integer-typed"));
    }
    let width = type_width(ty);
    let mut slots = [None; 2];
    let mut n = 0_usize;
    for i in 0..count {
        let elem = &data[i * width..(i + 1) * width];
        let Some(code) = super::record::gt_element(ty, elem).map_err(malformed_at)? else {
            break;
        };
        if code < 2 {
            // The missing allele carries only the phasing bit.
            if n < 2 {
                slots[n] = None;
                n += 1;
            }
            continue;
        }
        let allele = usize::from(u16::try_from(code).unwrap_or(u16::MAX) >> 1).saturating_sub(1);
        if allele >= n_allele {
            return Err(malformed_at(&format!(
                "GT allele index {allele} exceeds ALT count {}",
                n_allele.saturating_sub(1)
            )));
        }
        if n < 2 {
            slots[n] = Some(u16::try_from(allele).unwrap_or(u16::MAX));
            n += 1;
        }
    }
    Ok((slots, n))
}

/// Decodes a whole record's GT column into `codes` (one per sample),
/// mirroring the text scanner: an absent GT field yields [`GtCode::MISSING`]
/// for every sample, and a sample-count disagreement is the scanner's error.
/// The `bool` result is whether the record carried a GT field at all (the
/// all-`MISSING` outcome of a GT-less record must stay distinguishable from
/// a GT whose values are all `./.`).
///
/// The error locus is built lazily from `query` — the happy path allocates
/// nothing.
///
/// # Errors
///
/// [`IoError::MalformedRecord`] for a malformed record and
/// [`IoError::PolicyRejected`] for genotypes the policies forbid.
pub fn decode_record_gt(
    record: &BcfRecordRef<'_>,
    header: &SiteHeader,
    query: &GtQuery<'_>,
    table: &GtTable,
    codes: &mut Vec<GtCode>,
    counters: &mut PolicyCounters,
) -> Result<bool, IoError> {
    decode_record_gt_impl(record, header, query, table, codes, counters, None, None)
}

/// [`decode_record_gt`] restricted to `active` samples: an inactive column
/// yields [`GtCode::MISSING`] without running the genotype policies.
///
/// This is the BCF-native counterpart of
/// [`GenotypeScanner::with_active_samples`](crate::vcf::genotype_scan::GenotypeScanner::with_active_samples):
/// excluding a sample must be equivalent to removing its column before
/// decoding (v29 V29-02 / D29-03).
///
/// # Errors
///
/// As [`decode_record_gt`].
pub fn decode_record_gt_masked(
    record: &BcfRecordRef<'_>,
    header: &SiteHeader,
    query: &GtQuery<'_>,
    table: &GtTable,
    active: &KeepMask,
    codes: &mut Vec<GtCode>,
    counters: &mut PolicyCounters,
) -> Result<bool, IoError> {
    debug_assert_eq!(active.len(), query.n_samples);
    decode_record_gt_impl(
        record,
        header,
        query,
        table,
        codes,
        counters,
        None,
        Some(active),
    )
}

/// [`decode_record_gt`], also filling the reusable FORMAT table so a caller
/// that needs typed columns afterwards does not walk the descriptors a
/// second time (doc_v28 V28-16).
///
/// `formats` is rebuilt for every record, including a record without a GT
/// key, so a caller may always read it after a successful decode.
///
/// # Errors
///
/// As [`decode_record_gt`].
pub fn decode_record_gt_with_format(
    record: &BcfRecordRef<'_>,
    header: &SiteHeader,
    query: &GtQuery<'_>,
    table: &GtTable,
    codes: &mut Vec<GtCode>,
    counters: &mut PolicyCounters,
    formats: Option<&mut FormatTable>,
) -> Result<bool, IoError> {
    decode_record_gt_impl(record, header, query, table, codes, counters, formats, None)
}

#[allow(clippy::too_many_arguments)]
fn decode_record_gt_impl(
    record: &BcfRecordRef<'_>,
    header: &SiteHeader,
    query: &GtQuery<'_>,
    table: &GtTable,
    codes: &mut Vec<GtCode>,
    counters: &mut PolicyCounters,
    formats: Option<&mut FormatTable>,
    active: Option<&KeepMask>,
) -> Result<bool, IoError> {
    let GtQuery {
        gt_key,
        n_samples,
        contig_name,
        pos1,
        vcf44,
    } = *query;
    codes.clear();
    if header.n_sample != n_samples {
        let reason = if header.n_sample > n_samples {
            "more sample columns than the header declares"
        } else {
            "sample column count disagrees with the header"
        };
        return Err(IoError::MalformedRecord {
            locus: Locus {
                contig_name: String::from_utf8_lossy(contig_name).into_owned(),
                pos1,
            },
            reason: reason.to_string(),
        });
    }
    let mut found: Option<(u8, usize, usize)> = None;
    // doc_v28 V28-16: one descriptor walk serves the GT lookup and, when the
    // caller asked, the FORMAT column table. The old code walked the same
    // descriptors twice and allocated a fresh `Vec` for the second walk.
    match formats {
        Some(table) => {
            if header.n_fmt > 0 {
                table.rebuild(record, header.n_fmt, header.n_sample)?;
            }
            if let Some(gt_key) = gt_key.filter(|_| header.n_fmt > 0) {
                found = table
                    .fields()
                    .iter()
                    .find(|(key, ..)| *key == gt_key)
                    .map(|&(_, ty, count, start)| (ty, count, start));
            }
        }
        None => {
            if let Some(gt_key) = gt_key.filter(|_| header.n_fmt > 0) {
                let indiv = record.indiv;
                let mut ci = Cursor::new(indiv);
                for _ in 0..header.n_fmt {
                    let key = ci.typed_int()?;
                    let key =
                        u32::try_from(key).map_err(|_| malformed("BCF FORMAT key is negative"))?;
                    let (ty, count) = ci.typed()?;
                    let start = indiv.len() - ci.remaining();
                    let bytes = type_width(ty)
                        .checked_mul(count)
                        .and_then(|v| v.checked_mul(header.n_sample))
                        .ok_or_else(|| malformed("BCF FORMAT field size overflows"))?;
                    ci.take(bytes)?;
                    if key == gt_key {
                        found = Some((ty, count, start));
                    }
                }
            }
        }
    }
    let Some((ty, count, start)) = found else {
        codes.resize(n_samples, GtCode::MISSING);
        return Ok(false);
    };
    // `(NULL, count=0)` is how htslib stores an entirely missing GT column
    // (bcftools renders one `.` per sample). It carries no elements to
    // decode, but the GT field still exists (doc_v28 V28-08).
    if ty == value_type::NULL && count == 0 {
        codes.resize(n_samples, GtCode::MISSING);
        return Ok(true);
    }
    // The error locus stays borrowed (A6-11: no per-record String; the
    // `Locus` is built only when an error carries it).
    let locus = (contig_name, pos1);
    let mut text = Vec::new();
    let stride = type_width(ty) * count;
    // The record's allele bound (REF included → ALT count is one less):
    // an allele index at or above it names an allele the record does not
    // declare (doc_v25 V25-01).
    let n_alt = u32::try_from(header.n_allele.saturating_sub(1)).unwrap_or(u32::MAX);
    if ty == value_type::INT8 && count == 2 {
        for s in 0..header.n_sample {
            if active.is_some_and(|mask| !mask[s]) {
                codes.push(GtCode::MISSING);
                continue;
            }
            let b = &record.indiv[start + s * 2..start + s * 2 + 2];
            codes.push(table.decode_pair_with(
                b[0],
                b[1],
                &locus,
                header.n_allele,
                counters,
                vcf44,
            )?);
        }
    } else {
        for s in 0..header.n_sample {
            if active.is_some_and(|mask| !mask[s]) {
                codes.push(GtCode::MISSING);
                continue;
            }
            let off = start + s * stride;
            let data = &record.indiv[off..off + stride];
            codes.push(table.decode_slow(ty, data, &locus, n_alt, counters, &mut text, vcf44)?);
        }
    }
    Ok(true)
}

/// The FORMAT dictionary index of `tag`, resolved once per run (the fast
/// path compares integers, never names).
#[must_use]
pub fn format_key(dict: &BcfDict, tag: &[u8]) -> Option<u32> {
    dict.key_of(tag)
}

#[cfg(test)]
mod tests {
    use super::*;
    use reed_core::genotype::PloidyPolicy;

    fn policies() -> GenotypePolicies {
        GenotypePolicies {
            ploidy: PloidyPolicy::HaploidAsHalfCall,
            ..GenotypePolicies::default()
        }
    }

    /// doc_v26 V26-05: `site_header` is the numeric path's only record
    /// entry; a record declaring zero alleles has no REF and must be
    /// rejected there, before any consumer treats it as a valid site.
    #[test]
    fn v26_site_header_rejects_zero_alleles() {
        let mut shared = Vec::new();
        shared.extend_from_slice(&0_i32.to_le_bytes()); // chrom
        shared.extend_from_slice(&9_i32.to_le_bytes()); // pos0
        shared.extend_from_slice(&1_i32.to_le_bytes()); // rlen
        shared.extend_from_slice(&0_u32.to_le_bytes()); // qual
        shared.extend_from_slice(&0_u32.to_le_bytes()); // n_allele = 0, n_info = 0
        shared.extend_from_slice(&0_u32.to_le_bytes()); // 0 FORMAT, 0 samples
        let rec = BcfRecordRef {
            shared: &shared,
            indiv: &[],
        };
        let err = site_header(&rec).expect_err("n_allele = 0 must be MalformedRecord");
        assert!(matches!(err, IoError::MalformedRecord { .. }), "{err:?}");
        // The control: one allele (REF only) is a valid sites-only record.
        let mut valid = Vec::new();
        valid.extend_from_slice(&0_i32.to_le_bytes());
        valid.extend_from_slice(&9_i32.to_le_bytes());
        valid.extend_from_slice(&1_i32.to_le_bytes());
        valid.extend_from_slice(&0_u32.to_le_bytes());
        valid.extend_from_slice(&(1_u32 << 16).to_le_bytes());
        valid.extend_from_slice(&0_u32.to_le_bytes());
        valid.extend_from_slice(&[0x07]);
        valid.extend_from_slice(&[0x17, b'A']);
        valid.extend_from_slice(&[0x01]);
        let rec = BcfRecordRef {
            shared: &valid,
            indiv: &[],
        };
        assert_eq!(
            site_header(&rec).expect("invariant: one allele").n_allele,
            1
        );
    }

    /// Every table entry equals the text path on the same rendered bytes
    /// (§4's unit gate), and the trigger bits match the text decoder's
    /// counters.
    #[test]
    fn table_matches_the_text_path_for_every_byte_pair() {
        let table = GtTable::new(policies());
        let locus = Locus {
            contig_name: "c".into(),
            pos1: 1,
        };
        for pair in 0..65_536_u32 {
            let b0 = u8::try_from(pair & 0xFF).expect("invariant: low byte");
            let b1 = u8::try_from(pair >> 8).expect("invariant: high byte");
            let mut native_counters = PolicyCounters::default();
            // The int8 domain: every allele a byte pair can name is legal,
            // so the comparison exercises the table, not the record bound.
            let got = table.decode_pair(b0, b1, &locus, 128, &mut native_counters);
            // The text path: render the same two bytes and decode.
            let mut rendered = Vec::new();
            if let Err(_negative) = fmt_gt(&mut rendered, value_type::INT8, 2, &[b0, b1], false) {
                // Negative GT elements are rejected by both paths (V28-04).
                assert!(got.is_err(), "pair {b0:#04x},{b1:#04x}: native {got:?}");
                continue;
            }
            let mut text_counters = PolicyCounters::default();
            let mut token = [None; 2];
            let want = decode_gt_text(
                &rendered,
                &locus,
                127,
                &policies(),
                &mut text_counters,
                &mut token,
            );
            match (got, want) {
                (Ok(a), Ok(b)) => {
                    assert_eq!(a.raw(), b.raw(), "pair {b0:#04x},{b1:#04x}");
                    assert_eq!(
                        pack_triggers(&native_counters),
                        pack_triggers(&text_counters),
                        "counters for pair {b0:#04x},{b1:#04x}"
                    );
                }
                (Err(_), Err(_)) => {}
                (a, b) => panic!("pair {b0:#04x},{b1:#04x}: native {a:?} vs text {b:?}"),
            }
        }
    }

    /// doc_v28 V28-07: the table is built for 4.3 text, but a 4.4 header's
    /// first-allele prefix can carry phasing the inferred form drops.
    #[test]
    fn vcf44_first_allele_phasing_uses_the_header_aware_entries() {
        let table = GtTable::new(policies());
        let locus = Locus {
            contig_name: "c".into(),
            pos1: 1,
        };
        let mut counters = PolicyCounters::default();
        // `|0/1`: first REF allele phased, second unphased.
        let old = table
            .decode_pair_with(3, 4, &locus, 2, &mut counters, false)
            .expect("invariant: int8 pair");
        assert!(!old.phased(), "the inferred 4.3 form is unphased");
        let new = table
            .decode_pair_with(3, 4, &locus, 2, &mut counters, true)
            .expect("invariant: int8 pair");
        assert!(new.phased(), "the explicit 4.4 prefix keeps phasing");

        let mut text = Vec::new();
        let missing = table
            .decode_slow(
                value_type::INT8,
                &[1],
                &locus,
                2,
                &mut counters,
                &mut text,
                true,
            )
            .expect("invariant: phased missing haploid");
        assert!(missing.phased() && missing.is_fully_missing());
    }

    /// doc_v28 V28-08: htslib stores an all-missing GT as NULL/count=0. The
    /// decoder must report `has_gt = true` and one missing code per sample,
    /// not reject the record as non-integer.
    #[test]
    fn null_count_zero_gt_decodes_as_missing() {
        let dict =
            BcfDict::from_header(b"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n")
                .expect("invariant: dict");
        let gt_key = dict.key_of(b"GT").expect("invariant: GT declared");
        let mut shared = Vec::new();
        shared.extend_from_slice(&0_i32.to_le_bytes()); // chrom
        shared.extend_from_slice(&0_i32.to_le_bytes()); // pos0
        shared.extend_from_slice(&1_i32.to_le_bytes()); // rlen
        shared.extend_from_slice(&0_u32.to_le_bytes()); // qual
        shared.extend_from_slice(&(2_u32 << 16).to_le_bytes()); // REF+ALT
        shared.extend_from_slice(&((1_u32 << 24) | 1).to_le_bytes()); // 1 fmt, 1 sample
        // ID, REF, ALT as typed strings, then empty FILTER and INFO.
        for text in [b".".as_slice(), b"A", b"G"] {
            shared.push(u8::try_from((text.len() << 4) | usize::from(value_type::CHAR)).unwrap());
            shared.extend_from_slice(text);
        }
        shared.push(0x01); // FILTER: INT8, count 0
        shared.push(0x01); // INFO: INT8, count 0
        let indiv = vec![
            0x11, // FORMAT key: INT8, count 1
            u8::try_from(gt_key).expect("invariant: small dictionary"),
            0x00, // GT: NULL, count 0
        ];
        let rec = BcfRecordRef {
            shared: &shared,
            indiv: &indiv,
        };
        let header = site_header(&rec).expect("invariant: valid site header");
        let table = GtTable::new(policies());
        let mut codes = Vec::new();
        let mut counters = PolicyCounters::default();
        let has_gt = decode_record_gt(
            &rec,
            &header,
            &GtQuery {
                gt_key: Some(gt_key),
                n_samples: 1,
                contig_name: b"c",
                pos1: 1,
                vcf44: false,
            },
            &table,
            &mut codes,
            &mut counters,
        )
        .expect("NULL GT must decode, not error");
        assert!(has_gt);
        assert_eq!(codes, vec![GtCode::MISSING]);
    }

    /// doc_v25 V25-01 (F1): the BCF decode applies the same allele-domain
    /// rule as the text scanner. A framed record naming allele 9 on a
    /// single-ALT site must be refused at decode, not handed to consumers
    /// (the phantom allele then reaches the `-Ob` re-encode as `(9+1)<<1`,
    /// which htslib refuses to read back).
    #[test]
    fn decode_record_gt_rejects_alleles_beyond_the_record() {
        use crate::bcf::{BcfDict, BcfEncoder, bcf_header_text, record_ref};

        let header_text = b"##fileformat=VCFv4.2\n\
            ##FILTER=<ID=PASS,Description=\"All\">\n\
            ##contig=<ID=Chr01,length=1000>\n\
            ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n";
        let header: noodles::vcf::Header = std::str::from_utf8(header_text)
            .expect("invariant: utf8")
            .parse()
            .expect("invariant: header");
        let dict_text = bcf_header_text(&header).expect("invariant: bcf header");
        let dict = BcfDict::from_header(&dict_text).expect("invariant: dict");

        // Frame the record through the noodles encoder: it does not check
        // the allele domain, which is exactly how such a record reaches reed.
        let mut enc = BcfEncoder::new();
        enc.write_header(&header).expect("invariant: header");
        let mut framed = Vec::new();
        enc.encode_record(
            &header,
            b"Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/9\n",
            &mut framed,
        )
        .expect("invariant: noodles accepts the out-of-domain pair");
        let rec = record_ref(&framed).expect("invariant: framed");
        let site = site_header(&rec).expect("invariant: site header");

        let table = GtTable::new(GenotypePolicies {
            ploidy: PloidyPolicy::HaploidAsHalfCall,
            multiallelic: reed_core::genotype::MultiallelicPolicy::CollapseAlt,
            ..GenotypePolicies::default()
        });
        let mut codes = Vec::new();
        let mut counters = PolicyCounters::default();
        let err = decode_record_gt(
            &rec,
            &site,
            &GtQuery {
                gt_key: format_key(&dict, b"GT"),
                n_samples: 1,
                contig_name: b"Chr01",
                pos1: 100,
                vcf44: false,
            },
            &table,
            &mut codes,
            &mut counters,
        )
        .expect_err("allele 9 on a single-ALT record must be refused at decode");
        assert!(matches!(err, IoError::MalformedRecord { .. }), "{err:?}");
        assert!(
            err.to_string().contains("exceeds ALT count 1"),
            "unexpected error: {err}"
        );
    }

    /// v29 V29-02 / D29-03: the BCF-native path must decode exactly the
    /// active columns, so an excluded sample's out-of-domain GT cannot veto
    /// the record.
    #[test]
    fn decode_record_gt_masked_skips_inactive_samples() {
        use crate::bcf::{BcfDict, BcfEncoder, bcf_header_text, record_ref};

        let header_text = b"##fileformat=VCFv4.2\n\
            ##FILTER=<ID=PASS,Description=\"All\">\n\
            ##contig=<ID=Chr01,length=1000>\n\
            ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"gt\">\n\
            #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n";
        let header: noodles::vcf::Header = std::str::from_utf8(header_text)
            .expect("invariant: utf8")
            .parse()
            .expect("invariant: header");
        let dict_text = bcf_header_text(&header).expect("invariant: bcf header");
        let dict = BcfDict::from_header(&dict_text).expect("invariant: dict");
        let mut enc = BcfEncoder::new();
        enc.write_header(&header).expect("invariant: header");
        let mut framed = Vec::new();
        enc.encode_record(
            &header,
            b"Chr01\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/9\n",
            &mut framed,
        )
        .expect("invariant: noodles accepts the out-of-domain pair");
        let rec = record_ref(&framed).expect("invariant: framed");
        let site = site_header(&rec).expect("invariant: site header");
        let table = GtTable::new(GenotypePolicies {
            ploidy: PloidyPolicy::HaploidAsHalfCall,
            multiallelic: reed_core::genotype::MultiallelicPolicy::CollapseAlt,
            ..GenotypePolicies::default()
        });
        let query = GtQuery {
            gt_key: format_key(&dict, b"GT"),
            n_samples: 2,
            contig_name: b"Chr01",
            pos1: 100,
            vcf44: false,
        };
        let active = KeepMask::try_new(2, vec![true, false]).expect("invariant: length");
        let mut codes = Vec::new();
        let mut counters = PolicyCounters::default();
        decode_record_gt_masked(
            &rec,
            &site,
            &query,
            &table,
            &active,
            &mut codes,
            &mut counters,
        )
        .expect("the excluded sample's GT must not be decoded");
        assert_eq!(codes.len(), 2);
        assert_ne!(codes[0], GtCode::MISSING);
        assert_eq!(codes[1], GtCode::MISSING);
    }
}
