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

use super::dict::BcfDict;
use super::record::{BcfRecordRef, Cursor, malformed, type_width, value_type};
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
/// header.
pub fn site_header(record: &BcfRecordRef<'_>) -> Result<SiteHeader, IoError> {
    let mut c = Cursor::new(record.shared);
    let chrom = c.i32()?;
    let pos0 = c.i32()?;
    let _rlen = c.i32()?;
    let _qual = c.u32()?;
    let na_ni = c.u32()?;
    let nf_ns = c.u32()?;
    Ok(SiteHeader {
        chrom,
        pos0,
        n_allele: (na_ni >> 16) as usize,
        n_sample: (nf_ns & 0x00FF_FFFF) as usize,
        n_fmt: (nf_ns >> 24) as usize,
    })
}

/// A policy trigger fired by a table entry (one bit per counter).
const TRIG_HAPLOID_HALF: u8 = 1 << 0;
const TRIG_HAPLOID_HOMOZYGOUS: u8 = 1 << 1;
const TRIG_HAPLOID_MISSING: u8 = 1 << 2;
const TRIG_HALF_CALL_MISSING: u8 = 1 << 3;
const TRIG_MULTIALLELIC: u8 = 1 << 4;
/// The entry's policies reject this pair; the decoder re-runs the text path
/// to recover the exact error.
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
    policies: GenotypePolicies,
}

impl GtTable {
    /// Builds the table for `policies` by rendering each int8 byte pair and
    /// running the text decoder on it.
    #[must_use]
    pub fn new(policies: GenotypePolicies) -> Self {
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
            fmt_gt(&mut text, value_type::INT8, raw.len(), &raw, false);
            let mut counters = PolicyCounters::default();
            let mut token = [None; 2];
            if let Ok(code) = decode_gt_text(&text, &locus, &policies, &mut counters, &mut token) {
                entries[pair as usize] = GtEntry {
                    code,
                    triggers: pack_triggers(&counters),
                };
            }
        }
        Self {
            entries: entries.into_boxed_slice(),
            policies,
        }
    }

    /// The policies the table (and its slow path) decode under.
    #[must_use]
    pub const fn policies(&self) -> &GenotypePolicies {
        &self.policies
    }

    /// Decodes one int8 two-element genotype.
    ///
    /// # Errors
    ///
    /// The entries' policies reject this genotype (`Reject`-style ploidy or
    /// multiallelic policies).
    pub fn decode_pair<S: LocusSource>(
        &self,
        b0: u8,
        b1: u8,
        locus: &S,
        counters: &mut PolicyCounters,
    ) -> Result<GtCode, IoError> {
        let entry = self.entries[usize::from(b0) | (usize::from(b1) << 8)];
        if entry.triggers == TRIG_ERROR {
            // Recover the exact error (and locus) through the text path.
            let mut text = Vec::new();
            return self.decode_slow(value_type::INT8, &[b0, b1], locus, counters, &mut text);
        }
        apply_triggers(entry.triggers, counters);
        Ok(entry.code)
    }

    /// Decodes one sample's GT field of any integer encoding through the
    /// text path (render, then parse).
    fn decode_slow<S: LocusSource>(
        &self,
        ty: u8,
        data: &[u8],
        locus: &S,
        counters: &mut PolicyCounters,
        text: &mut Vec<u8>,
    ) -> Result<GtCode, IoError> {
        if !matches!(ty, value_type::INT8 | value_type::INT16 | value_type::INT32) {
            return Err(IoError::MalformedRecord {
                locus: locus.error_locus(),
                reason: "BCF GT field is not integer-typed",
            });
        }
        let count = data.len() / type_width(ty);
        text.clear();
        fmt_gt(text, ty, count, data, false);
        let mut token = [None; 2];
        decode_gt_text(text, locus, &self.policies, counters, &mut token)
    }
}

/// One bit per counter the entry fires (each entry is one genotype, so every
/// field is 0 or 1).
fn pack_triggers(c: &PolicyCounters) -> u8 {
    let mut t = 0_u8;
    if c.n_haploid_as_half_call > 0 {
        t |= TRIG_HAPLOID_HALF;
    }
    if c.n_haploid_as_homozygous > 0 {
        t |= TRIG_HAPLOID_HOMOZYGOUS;
    }
    if c.n_haploid_as_missing > 0 {
        t |= TRIG_HAPLOID_MISSING;
    }
    if c.n_half_call_as_missing > 0 {
        t |= TRIG_HALF_CALL_MISSING;
    }
    if c.n_multiallelic_collapsed > 0 {
        t |= TRIG_MULTIALLELIC;
    }
    t
}

fn apply_triggers(t: u8, c: &mut PolicyCounters) {
    if t & TRIG_HAPLOID_HALF != 0 {
        c.n_haploid_as_half_call += 1;
    }
    if t & TRIG_HAPLOID_HOMOZYGOUS != 0 {
        c.n_haploid_as_homozygous += 1;
    }
    if t & TRIG_HAPLOID_MISSING != 0 {
        c.n_haploid_as_missing += 1;
    }
    if t & TRIG_HALF_CALL_MISSING != 0 {
        c.n_half_call_as_missing += 1;
    }
    if t & TRIG_MULTIALLELIC != 0 {
        c.n_multiallelic_collapsed += 1;
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
    let GtQuery {
        gt_key,
        n_samples,
        contig_name,
        pos1,
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
            reason,
        });
    }
    let mut found: Option<(u8, usize, usize)> = None;
    if let Some(gt_key) = gt_key.filter(|_| header.n_fmt > 0) {
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
                found = Some((ty, count, start));
            }
        }
    }
    let Some((ty, count, start)) = found else {
        codes.resize(n_samples, GtCode::MISSING);
        return Ok(false);
    };
    // The error locus stays borrowed (A6-11: no per-record String; the
    // `Locus` is built only when an error carries it).
    let locus = (contig_name, pos1);
    let mut text = Vec::new();
    let stride = type_width(ty) * count;
    if ty == value_type::INT8 && count == 2 {
        for s in 0..header.n_sample {
            let b = &record.indiv[start + s * 2..start + s * 2 + 2];
            codes.push(table.decode_pair(b[0], b[1], &locus, counters)?);
        }
    } else {
        for s in 0..header.n_sample {
            let off = start + s * stride;
            let data = &record.indiv[off..off + stride];
            codes.push(table.decode_slow(ty, data, &locus, counters, &mut text)?);
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
            let got = table.decode_pair(b0, b1, &locus, &mut native_counters);
            // The text path: render the same two bytes and decode.
            let mut rendered = Vec::new();
            fmt_gt(&mut rendered, value_type::INT8, 2, &[b0, b1], false);
            let mut text_counters = PolicyCounters::default();
            let mut token = [None; 2];
            let want = decode_gt_text(
                &rendered,
                &locus,
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
}
