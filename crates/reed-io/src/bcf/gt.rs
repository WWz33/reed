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
/// header, and when the record declares no allele (doc_v26 V26-05):
/// `n_allele` counts REF plus the ALTs, so zero means the shared block has no
/// allele string at all. reed used to accept it — the renderer printed `.`
/// for REF and ALT and the GT path treated every code as "0 alleles" — while
/// htslib rejects such a record; a BCF that no other tool accepts must not
/// flow through filter/stats looking valid.
pub fn site_header(record: &BcfRecordRef<'_>) -> Result<SiteHeader, IoError> {
    let mut c = Cursor::new(record.shared);
    let chrom = c.i32()?;
    let pos0 = c.i32()?;
    let _rlen = c.i32()?;
    let _qual = c.u32()?;
    let na_ni = c.u32()?;
    let nf_ns = c.u32()?;
    let n_allele = (na_ni >> 16) as usize;
    if n_allele == 0 {
        return Err(malformed(
            "BCF record declares zero alleles (n_allele = 0); REF is required",
        ));
    }
    Ok(SiteHeader {
        chrom,
        pos0,
        n_allele,
        n_sample: (nf_ns & 0x00FF_FFFF) as usize,
        n_fmt: (nf_ns >> 24) as usize,
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
            // INT8 is an integer type, so this never fails (the guard only
            // rejects CHAR/FLOAT, doc_v23 V23-02).
            fmt_gt(&mut text, value_type::INT8, raw.len(), &raw, false)
                .expect("invariant: INT8 is integer-typed");
            let mut counters = PolicyCounters::default();
            let mut token = [None; 2];
            // The table spans every allele an int8 pair can name (0..=127);
            // each record's own bound is applied at lookup (doc_v25 V25-01).
            if let Ok(code) = decode_gt_text(&text, &locus, 127, &policies, &mut counters, &mut token) {
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
        for b in [b0, b1] {
            let v = b as i8;
            if v < 2 {
                // 0/1: missing allele (bare phasing bit). Negatives are
                // invalid encodings the table itself rejects (TRIG_ERROR).
                continue;
            }
            let allele = usize::from((v as u8 >> 1) - 1);
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
        let entry = self.entries[usize::from(b0) | (usize::from(b1) << 8)];
        if entry.triggers == TRIG_ERROR {
            // Recover the exact error (and locus) through the text path.
            // The pair's bytes were domain-checked above, so the slow path
            // can run unbounded (its errors here are syntax or policy).
            let mut text = Vec::new();
            return self.decode_slow(value_type::INT8, &[b0, b1], locus, u32::MAX, counters, &mut text);
        }
        apply_triggers(entry.triggers, counters);
        Ok(entry.code)
    }

    /// Decodes one sample's GT field of any integer encoding through the
    /// text path (render, then parse). `n_alt` bounds the allele domain
    /// (`u32::MAX` = unbounded; doc_v25 V25-01).
    fn decode_slow<S: LocusSource>(
        &self,
        ty: u8,
        data: &[u8],
        locus: &S,
        n_alt: u32,
        counters: &mut PolicyCounters,
        text: &mut Vec<u8>,
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
        fmt_gt(text, ty, count, data, false)
            .map_err(|reason| IoError::MalformedRecord {
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
            reason: reason.to_string(),
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
    // The record's allele bound (REF included → ALT count is one less):
    // an allele index at or above it names an allele the record does not
    // declare (doc_v25 V25-01).
    let n_alt = u32::try_from(header.n_allele.saturating_sub(1)).unwrap_or(u32::MAX);
    if ty == value_type::INT8 && count == 2 {
        for s in 0..header.n_sample {
            let b = &record.indiv[start + s * 2..start + s * 2 + 2];
            codes.push(table.decode_pair(b[0], b[1], &locus, header.n_allele, counters)?);
        }
    } else {
        for s in 0..header.n_sample {
            let off = start + s * stride;
            let data = &record.indiv[off..off + stride];
            codes.push(table.decode_slow(ty, data, &locus, n_alt, counters, &mut text)?);
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
            fmt_gt(&mut rendered, value_type::INT8, 2, &[b0, b1], false)
                .expect("invariant: INT8 is integer-typed");
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
        let mut enc = BcfEncoder::new().expect("invariant: noodles");
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
            },
            &table,
            &mut codes,
            &mut counters,
        )
        .expect_err("allele 9 on a single-ALT record must be refused at decode");
        assert!(
            matches!(err, IoError::MalformedRecord { .. }),
            "{err:?}"
        );
        assert!(
            err.to_string().contains("exceeds ALT count 1"),
            "unexpected error: {err}"
        );
    }
}
