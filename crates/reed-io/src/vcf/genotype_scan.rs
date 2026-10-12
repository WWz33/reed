//! [`GenotypeScanner`]: the one record decoder (`docs/rounds/v2/03-architecture.md`
//! §3.4). One pass over the sample region produces a [`GtCode`] per sample
//! plus the planned integer/float FORMAT columns; the genotype policies
//! execute here and only here, and every buffer is reused across records
//! (fixing 01-review D10: no per-sample allocation, no layout clone, no
//! second scan for the sample count).
//!
//! Correctness is pinned by a differential test against noodles
//! (`tests/genotype_scan_vs_noodles.rs`).

use std::ops::Range;

use reed_core::genotype::{GenotypePolicies, GtCode, Ploidy, PolicyCounters};

use crate::error::{IoError, Locus};
use crate::vcf::format_layout::FormatLayoutCache;
use crate::vcf::line::LineRef;
use crate::vcf::site::SiteView;

/// Sentinel for a missing integer sample value: real DP/GQ values never
/// reach `u32::MAX`.
pub const INT_MISSING: u32 = u32::MAX;

/// One record's decoded genotypes and planned FORMAT columns, borrowed from
/// the scanner (valid until its next `scan`).
#[derive(Clone, Copy, Debug)]
pub struct ScanOutput<'s> {
    /// One code per sample, in sample order.
    pub codes: &'s [GtCode],
    /// Planned integer columns, aligned with the scanner's integer keys;
    /// missing values are [`INT_MISSING`].
    pub int_columns: &'s [Vec<u32>],
    /// Planned float columns, aligned with the scanner's float keys; missing
    /// values are `NaN`.
    pub float_columns: &'s [Vec<f32>],
    /// Raw allele indices (2 per sample) when
    /// [`GenotypeScanner::tracking_alleles`] is on: the multi-allelic index
    /// survives here for consumers that count per allele.
    pub alleles: Option<&'s [Option<u16>]>,
    /// The record's FORMAT carried a GT key. `false` decodes every sample as
    /// missing — distinct from a GT whose values are all `./.` (bcftools
    /// skips `-g/--genotype` when the GT key is absent).
    pub has_gt: bool,
}

/// Single-pass sample-region decoder with reusable buffers.
pub struct GenotypeScanner {
    n_samples: usize,
    policies: GenotypePolicies,
    layouts: FormatLayoutCache,
    codes: Vec<GtCode>,
    int_columns: Vec<Vec<u32>>,
    float_columns: Vec<Vec<f32>>,
    alleles: Vec<Option<u16>>,
    track_alleles: bool,
    /// Subfield ranges of the current sample column, up to `last_needed`.
    subfields: Vec<Range<usize>>,
    /// Raw allele indices of the current sample.
    token: [Option<u16>; 2],
}

impl GenotypeScanner {
    /// Creates a scanner for `n_samples` header columns, decoding GT plus the
    /// named integer/float keys, under `policies`.
    #[must_use]
    pub fn new(
        n_samples: usize,
        int_keys: &[&[u8]],
        float_keys: &[&[u8]],
        policies: GenotypePolicies,
    ) -> Self {
        let layouts = FormatLayoutCache::new(int_keys, float_keys);
        let n_int = layouts.n_int_keys();
        let n_float = layouts.n_float_keys();
        Self {
            n_samples,
            policies,
            layouts,
            codes: Vec::new(),
            int_columns: (0..n_int).map(|_| Vec::new()).collect(),
            float_columns: (0..n_float).map(|_| Vec::new()).collect(),
            alleles: Vec::new(),
            track_alleles: false,
            subfields: Vec::new(),
            token: [None, None],
        }
    }

    /// Also records each sample's raw allele indices (`ScanOutput::alleles`).
    #[must_use]
    pub fn tracking_alleles(mut self) -> Self {
        self.track_alleles = true;
        self
    }

    /// Number of integer columns.
    #[must_use]
    pub fn n_int_columns(&self) -> usize {
        self.int_columns.len()
    }

    /// Number of float columns.
    #[must_use]
    pub fn n_float_columns(&self) -> usize {
        self.float_columns.len()
    }

    /// Decodes one record.
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] for a sample-column count that disagrees
    /// with the header, unparsable POS/GT text or unparsable integer values;
    /// [`IoError::PolicyRejected`] for genotypes the declared policies
    /// forbid (called haploid under `RequireDiploid`, allele index above 1
    /// under `Reject`, polyploid genotypes).
    pub fn scan<'s>(
        &'s mut self,
        line: &LineRef<'_>,
        site: &SiteView<'_>,
        counters: &mut PolicyCounters,
    ) -> Result<ScanOutput<'s>, IoError> {
        let Self {
            n_samples,
            policies,
            layouts,
            codes,
            int_columns,
            float_columns,
            alleles,
            track_alleles,
            subfields,
            token,
        } = self;
        codes.clear();
        for column in int_columns.iter_mut() {
            column.clear();
        }
        for column in float_columns.iter_mut() {
            column.clear();
        }
        if *track_alleles {
            alleles.clear();
        }

        let layout_ix = layouts.intern(line.format_text());
        let gt_index = layouts.get(layout_ix).gt_index;
        let last_needed = layouts.get(layout_ix).last_needed;

        let samples = line.samples_text();
        let mut field_start = 0_usize;
        let mut sample_count = 0_usize;
        while field_start <= samples.len() {
            let field_end = memchr::memchr(b'\t', &samples[field_start..])
                .map_or(samples.len(), |rel| field_start + rel);
            if sample_count >= *n_samples {
                return Err(IoError::MalformedRecord {
                    locus: site.locus(),
                    reason: "more sample columns than the header declares",
                });
            }
            let field = &samples[field_start..field_end];
            split_needed(subfields, field, last_needed);

            let code = match gt_index.and_then(|gi| subfields.get(usize::from(gi))) {
                Some(range) => decode_gt(&field[range.clone()], site, policies, counters, token)?,
                None => {
                    *token = [None, None];
                    GtCode::MISSING
                }
            };
            codes.push(code);
            if *track_alleles {
                alleles.extend_from_slice(token);
            }

            for (key, column) in layouts
                .get(layout_ix)
                .int_indices
                .iter()
                .zip(int_columns.iter_mut())
            {
                column.push(
                    key.and_then(|ki| subfields.get(usize::from(ki)))
                        .map_or(INT_MISSING, |range| parse_int(&field[range.clone()])),
                );
            }
            for (key, column) in layouts
                .get(layout_ix)
                .float_indices
                .iter()
                .zip(float_columns.iter_mut())
            {
                column.push(
                    key.and_then(|ki| subfields.get(usize::from(ki)))
                        .map_or(f32::NAN, |range| parse_float(&field[range.clone()])),
                );
            }

            sample_count += 1;
            field_start = field_end + 1;
        }
        if sample_count != *n_samples {
            return Err(IoError::MalformedRecord {
                locus: site.locus(),
                reason: "sample column count disagrees with the header",
            });
        }

        Ok(ScanOutput {
            codes,
            int_columns,
            float_columns,
            alleles: track_alleles.then_some(alleles.as_slice()),
            has_gt: gt_index.is_some(),
        })
    }
}

/// Splits the subfields a record's plan needs: ranges 0..`needed`, then the
/// rest of the field is skipped (trailing subfields are never scanned).
fn split_needed(ranges: &mut Vec<Range<usize>>, field: &[u8], needed: usize) {
    ranges.clear();
    if needed == 0 {
        return;
    }
    // One pass over the (short) field instead of one `memchr` per subfield
    // (doc_v3 03 O6, level 1: the scan is call-bound on ~30-byte fields).
    let mut from = 0;
    for (i, &b) in field.iter().enumerate() {
        if b == b':' {
            ranges.push(from..i);
            from = i + 1;
            if ranges.len() == needed {
                return;
            }
        }
    }
    ranges.push(from..field.len());
}

/// A cheap source of a record's error locus: the text scanner passes the
/// [`SiteView`] (no allocation on the happy path), the BCF fast path passes
/// the pieces of a locus it already has. Implementors build the `Locus`
/// (and its `String`) only when an error actually carries it.
pub trait LocusSource {
    /// The locus for an error message.
    fn error_locus(&self) -> Locus;
}

impl LocusSource for SiteView<'_> {
    fn error_locus(&self) -> Locus {
        self.locus()
    }
}

impl LocusSource for Locus {
    fn error_locus(&self) -> Locus {
        self.clone()
    }
}

impl LocusSource for (&[u8], u64) {
    fn error_locus(&self) -> Locus {
        Locus {
            contig_name: String::from_utf8_lossy(self.0).into_owned(),
            pos1: self.1,
        }
    }
}

/// Decodes one GT *text* subfield exactly as the scanner does, against any
/// cheap locus source (doc_v3 02-bcf.md §2.4: the BCF GT table is built by
/// rendering each byte pair and running this decoder, so the fast path is
/// the text path by construction).
///
/// # Errors
///
/// As [`GenotypeScanner::scan`] for the GT column.
pub(crate) fn decode_gt_text<S: LocusSource>(
    gt: &[u8],
    locus: &S,
    policies: &GenotypePolicies,
    counters: &mut PolicyCounters,
    token: &mut [Option<u16>; 2],
) -> Result<GtCode, IoError> {
    decode_gt(gt, locus, policies, counters, token)
}

/// Decodes one sample's GT subfield, applying the policies and filling
/// `token` with the raw allele indices (2 per sample).
fn decode_gt<S: LocusSource>(
    gt: &[u8],
    site: &S,
    policies: &GenotypePolicies,
    counters: &mut PolicyCounters,
    token: &mut [Option<u16>; 2],
) -> Result<GtCode, IoError> {
    // Fast path: exactly `[0|1|.][/|][0|1|.]`.
    if let Some((s0, s1, phased)) = fast_gt(gt) {
        *token = [
            (s0 != reed_core::genotype::STATE_MISSING)
                .then(|| u16::from(s0 == reed_core::genotype::STATE_ALT)),
            (s1 != reed_core::genotype::STATE_MISSING)
                .then(|| u16::from(s1 == reed_core::genotype::STATE_ALT)),
        ];
        return Ok(GtCode::from_states(s0, s1, phased, false));
    }
    if gt == b"." {
        // Haploid missing: both slots missing, no policy involved.
        *token = [None, None];
        return Ok(GtCode::from_states(
            reed_core::genotype::STATE_MISSING,
            reed_core::genotype::STATE_MISSING,
            false,
            true,
        ));
    }
    // Slow path: split on the separators. A leading `|`/`/` is the VCF 4.4
    // prefixed phasing of the first allele (htslib only honors it for ≥ 4.4
    // headers; accepted unconditionally here — a pre-4.4 file carrying it is
    // invalid input anyway).
    let (prefix_phased, gt) = match gt.first() {
        Some(b'|') => (true, &gt[1..]),
        Some(b'/') => (false, &gt[1..]),
        _ => (false, gt),
    };
    let mut slots = [None::<u16>; 2];
    let mut n_alleles = 0_usize;
    let mut phased = prefix_phased;
    let mut from = 0_usize;
    while from <= gt.len() {
        let end = memchr::memchr2(b'/', b'|', &gt[from..]).map_or(gt.len(), |rel| from + rel);
        if n_alleles >= 2 {
            return Err(IoError::PolicyRejected {
                locus: site.error_locus(),
                reason: "polyploid genotype (more than two alleles)",
            });
        }
        slots[n_alleles] = parse_allele(&gt[from..end], site)?;
        n_alleles += 1;
        if end == gt.len() {
            break;
        }
        phased |= gt[end] == b'|';
        from = end + 1;
    }
    *token = slots;
    let (ploidy, slots) = if n_alleles == 1 {
        (Ploidy::Haploid, slots)
    } else {
        (Ploidy::Diploid, slots)
    };
    GtCode::from_alleles(slots, ploidy, phased, policies, counters).map_err(|err| match err {
        reed_core::CoreError::PolicyRejected(reason) => IoError::PolicyRejected {
            locus: site.error_locus(),
            reason,
        },
        reed_core::CoreError::Unimplemented(reason) => IoError::UnsupportedInput(reason),
        other => IoError::MalformedRecord {
            locus: site.error_locus(),
            reason: reason_str(&other),
        },
    })
}

/// The 3-byte GT fast path: `[0|1|.][/|][0|1|.]`, returned as slot states.
fn fast_gt(gt: &[u8]) -> Option<(u8, u8, bool)> {
    use reed_core::genotype::{STATE_ALT, STATE_MISSING, STATE_REF};
    if gt.len() != 3 {
        return None;
    }
    let slot = |b: u8| match b {
        b'0' => Some(STATE_REF),
        b'1' => Some(STATE_ALT),
        b'.' => Some(STATE_MISSING),
        _ => None,
    };
    let s0 = slot(gt[0])?;
    let s1 = slot(gt[2])?;
    if gt[1] != b'/' && gt[1] != b'|' {
        return None;
    }
    Some((s0, s1, gt[1] == b'|'))
}

/// One allele token: `None` for `.`, else the index (multi-digit capable,
/// checked arithmetic — 01-review D14).
fn parse_allele<S: LocusSource>(token: &[u8], site: &S) -> Result<Option<u16>, IoError> {
    if token == b"." {
        return Ok(None);
    }
    if token.is_empty() || token.len() > 5 {
        return Err(IoError::MalformedRecord {
            locus: site.error_locus(),
            reason: "GT allele index is empty or too long",
        });
    }
    let mut value: u16 = 0;
    for &b in token {
        if !b.is_ascii_digit() {
            return Err(IoError::MalformedRecord {
                locus: site.error_locus(),
                reason: "GT allele contains a non-digit",
            });
        }
        // 5 digits can overflow u16 (99_999 > 65_535): reject the wrap.
        value = value
            .checked_mul(10)
            .and_then(|v| v.checked_add(u16::from(b - b'0')))
            .ok_or(IoError::MalformedRecord {
                locus: site.error_locus(),
                reason: "GT allele index exceeds u16",
            })?;
    }
    Ok(Some(value))
}

/// Integer FORMAT value; `.`, empty or non-digit text is [`INT_MISSING`].
#[must_use]
pub fn parse_int(text: &[u8]) -> u32 {
    if text == b"." || text.is_empty() {
        return INT_MISSING;
    }
    let mut value: u32 = 0;
    for &b in text {
        if !b.is_ascii_digit() {
            return INT_MISSING;
        }
        value = value.saturating_mul(10).saturating_add(u32::from(b - b'0'));
    }
    value
}

/// Float FORMAT value; `.`, empty or unparsable text is `NaN`.
#[must_use]
pub fn parse_float(text: &[u8]) -> f32 {
    if text == b"." || text.is_empty() {
        return f32::NAN;
    }
    std::str::from_utf8(text)
        .ok()
        .and_then(|s| s.parse::<f32>().ok())
        .unwrap_or(f32::NAN)
}

/// Short fallback text for `CoreError` variants without a static reason.
fn reason_str(err: &reed_core::CoreError) -> &'static str {
    match err {
        reed_core::CoreError::PolicyRejected(reason) => reason,
        _ => "genotype rejected by policy",
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use reed_core::genotype::{MultiallelicPolicy, PloidyPolicy};

    /// VCF 4.4 prefixed phasing parses like the equivalent unprefixed form
    /// (the bit of the first allele joins the genotype's phased flag).
    #[test]
    fn decode_gt_accepts_the_phasing_prefix() {
        let locus = ("Chr01".as_bytes(), 1_u64);
        let policies = permissive();
        let decode = |gt: &[u8]| {
            let mut counters = PolicyCounters::default();
            let mut token = [None; 2];
            decode_gt_text(gt, &locus, &policies, &mut counters, &mut token)
        };
        assert_eq!(decode(b"|0/1").unwrap(), decode(b"0|1").unwrap());
        assert_eq!(decode(b"/0|1").unwrap(), decode(b"0|1").unwrap());
        // The prefix is the first allele's explicit phasing bit.
        let prefixed = decode(b"|1").unwrap();
        assert!(prefixed.phased() && prefixed.is_haploid());
        assert_eq!(prefixed.slot_state(0), reed_core::genotype::STATE_ALT);
        let prefixed_missing = decode(b"|.").unwrap();
        assert!(prefixed_missing.is_fully_missing() && prefixed_missing.is_haploid());
        assert!(decode(b"|").is_err(), "a bare prefix has no allele");
    }

    fn line(text: &[u8]) -> LineRef<'_> {
        LineRef::parse(text).expect("invariant: test line has 10 columns")
    }

    fn permissive() -> GenotypePolicies {
        GenotypePolicies {
            ploidy: PloidyPolicy::HaploidAsHalfCall,
            multiallelic: MultiallelicPolicy::CollapseAlt,
            ..GenotypePolicies::default()
        }
    }

    fn scan<'s>(
        scanner: &'s mut GenotypeScanner,
        text: &'s [u8],
    ) -> (ScanOutput<'s>, PolicyCounters) {
        let lr = line(text);
        let site = SiteView::parse(&lr).expect("invariant: well-formed");
        let mut counters = PolicyCounters::default();
        let out = scanner
            .scan(&lr, &site, &mut counters)
            .expect("invariant: test record decodes");
        (out, counters)
    }

    #[test]
    fn decodes_slots_and_phasing() {
        let mut scanner = GenotypeScanner::new(3, &[], &[], GenotypePolicies::default());
        let (out, counters) = scan(
            &mut scanner,
            b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t1|0\t./.\n",
        );
        assert_eq!(counters, PolicyCounters::default());
        let states: Vec<u8> = out.codes.iter().map(|c| c.state()).collect();
        assert_eq!(states, vec![0b0100, 0b0001, 0b1010]);
        assert!(out.codes[1].phased() && !out.codes[0].phased());
        assert!(out.codes[2].is_fully_missing());
        assert!(out.alleles.is_none());
    }

    #[test]
    fn haploid_and_multiallelic_follow_the_policies() {
        let mut scanner = GenotypeScanner::new(3, &[], &[], permissive()).tracking_alleles();
        let (out, counters) = scan(
            &mut scanner,
            b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT\t1\t10/11\t0/.\n",
        );
        assert_eq!(out.codes[0].state(), 0b1001); // haploid ALT → half call
        assert!(out.codes[0].is_haploid());
        assert_eq!(out.codes[1].state(), 0b0101); // multi-digit collapses to ALT
        assert_eq!(out.codes[2].state(), 0b1000); // 0/. → slot0 REF, slot1 missing
        assert_eq!(counters.n_haploid_as_half_call, 1);
        assert_eq!(counters.n_multiallelic_collapsed, 2);
        assert_eq!(
            out.alleles.expect("invariant: tracked"),
            vec![Some(1), None, Some(10), Some(11), Some(0), None]
        );
    }

    #[test]
    fn default_policies_reject_called_haploids() {
        let mut scanner = GenotypeScanner::new(1, &[], &[], GenotypePolicies::default());
        let lr = line(b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT\t1\n");
        let site = SiteView::parse(&lr).expect("invariant: well-formed");
        let mut counters = PolicyCounters::default();
        let err = scanner.scan(&lr, &site, &mut counters).unwrap_err();
        assert!(matches!(err, IoError::PolicyRejected { .. }));
    }

    #[test]
    fn int_and_float_columns_with_missing_variants() {
        let mut scanner = GenotypeScanner::new(3, &[b"DP", b"GQ"], &[b"AB"], permissive());
        let (out, _) = scan(
            &mut scanner,
            b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ:AB\t0/1:9:40:0.5\t1/1:.:.:.\t1/1\n",
        );
        assert_eq!(out.int_columns[0], vec![9, INT_MISSING, INT_MISSING]);
        assert_eq!(out.int_columns[1], vec![40, INT_MISSING, INT_MISSING]);
        assert_eq!(out.float_columns[0][0], 0.5);
        assert!(out.float_columns[0][1].is_nan() && out.float_columns[0][2].is_nan());
    }

    #[test]
    fn gt_absent_from_format_is_all_missing() {
        let mut scanner = GenotypeScanner::new(3, &[b"DP"], &[], permissive());
        let (out, _) = scan(&mut scanner, b"c\t1\t.\tA\tG\t.\tPASS\t.\tDP\t5\t6\t7\n");
        assert!(out.codes.iter().all(|c| c.is_fully_missing()));
        assert_eq!(out.int_columns[0], vec![5, 6, 7]);
    }

    #[test]
    fn column_count_mismatch_is_an_error() {
        let mut scanner = GenotypeScanner::new(3, &[], &[], permissive());
        let lr = line(b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\n");
        let site = SiteView::parse(&lr).expect("invariant: well-formed");
        let mut counters = PolicyCounters::default();
        let err = scanner.scan(&lr, &site, &mut counters).unwrap_err();
        assert!(matches!(err, IoError::MalformedRecord { .. }));
    }

    #[test]
    fn polyploid_and_garbage_are_rejected() {
        let mut scanner = GenotypeScanner::new(3, &[], &[], permissive()).tracking_alleles();
        let lr = line(b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/1/2\t0/0\t0/0\n");
        let site = SiteView::parse(&lr).expect("invariant: well-formed");
        let mut counters = PolicyCounters::default();
        let err = scanner.scan(&lr, &site, &mut counters).unwrap_err();
        assert!(matches!(err, IoError::PolicyRejected { .. }));

        let lr = line(b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT\tx/1\t0/0\t0/0\n");
        let site = SiteView::parse(&lr).expect("invariant: well-formed");
        let err = scanner.scan(&lr, &site, &mut counters).unwrap_err();
        assert!(matches!(err, IoError::MalformedRecord { .. }));
    }

    #[test]
    fn subfield_split_stops_at_last_needed() {
        let mut ranges = Vec::new();
        split_needed(&mut ranges, b"0/1:9:40:60:99", 2);
        assert_eq!(ranges, vec![0..3, 4..5]);
        split_needed(&mut ranges, b"0/1", 3);
        assert_eq!(ranges, vec![0..3]);
        split_needed(&mut ranges, b"0/1:9", 0);
        assert_eq!(ranges.len(), 0);
    }
}
