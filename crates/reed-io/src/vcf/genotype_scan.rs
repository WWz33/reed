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
use crate::vcf::format_layout::{FormatLayoutCache, IntSpec};
use crate::vcf::line::{LineRef, SampleColumnsError};
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
    /// named integer/float keys, under `policies`. Every integer key is a
    /// plain scalar.
    #[must_use]
    pub fn new(
        n_samples: usize,
        int_keys: &[&[u8]],
        float_keys: &[&[u8]],
        policies: GenotypePolicies,
    ) -> Self {
        let specs: Vec<IntSpec> = int_keys
            .iter()
            .map(|k| IntSpec::owned(k.to_vec(), None))
            .collect();
        Self::new_with_specs(n_samples, &specs, float_keys, policies)
    }

    /// Creates a scanner whose integer columns carry a decode mode and an
    /// optional fallback key (see [`IntSpec`]).
    #[must_use]
    pub fn new_with_specs(
        n_samples: usize,
        int_specs: &[IntSpec],
        float_keys: &[&[u8]],
        policies: GenotypePolicies,
    ) -> Self {
        let layouts = FormatLayoutCache::new_with_specs(int_specs, float_keys);
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

        // The structural contract lives in `LineRef` (doc_v22 R22-3): the
        // scanner does not re-derive it from the tab count.
        line.validate_sample_cols(*n_samples)
            .map_err(|err| IoError::MalformedRecord {
                locus: site.locus(),
                // Name the counts: the old message dropped the
                // `found`/`expected` both variants carry, which is the one
                // detail a mismatched record needs (doc_v23 V23-29).
                reason: match err {
                    SampleColumnsError::TooFew { found, expected } => format!(
                        "sample column count disagrees with the header \
                         (found {found}, expected {expected})"
                    ),
                    SampleColumnsError::TooMany { found, expected } => format!(
                        "more sample columns than the header declares \
                         (first non-empty extra column at {found}, expected {expected})"
                    ),
                },
            })?;

        let layout_ix = layouts.intern(line.format_text());
        let gt_index = layouts.get(layout_ix).gt_index;
        let last_needed = layouts.get(layout_ix).last_needed;

        let samples = line.samples_text();
        let mut field_start = 0_usize;
        // `validate_sample_cols` already fixed the structural contract: at
        // least `n_samples` fields, and every field past them empty
        // (V22-R16). Decode exactly the declared samples; the trailing tabs
        // need no second check here.
        for _ in 0..*n_samples {
            let field_end = memchr::memchr(b'\t', &samples[field_start..])
                .map_or(samples.len(), |rel| field_start + rel);
            let field = &samples[field_start..field_end];
            split_needed(subfields, field, last_needed);

            let code = match gt_index.and_then(|gi| subfields.get(gi)) {
                Some(range) => {
                    decode_gt(&field[range.clone()], site, site.n_alt(), policies, counters, token)?
                }
                None => {
                    *token = [None, None];
                    GtCode::MISSING
                }
            };
            codes.push(code);
            if *track_alleles {
                alleles.extend_from_slice(token);
            }

            for (int_layout, column) in layouts
                .get(layout_ix)
                .ints
                .iter()
                .zip(int_columns.iter_mut())
            {
                // A record carrying the primary key never falls back, even
                // when a sample's value is missing (bcftools
                // `calc_sample_depth`).
                let value = if let Some(ki) = int_layout.index {
                    subfields
                        .get(ki)
                        .map_or(INT_MISSING, |range| parse_int(&field[range.clone()]))
                } else {
                    int_layout
                        .fallback_index
                        .and_then(|fi| subfields.get(fi))
                        .map_or(INT_MISSING, |range| parse_int_sum(&field[range.clone()]))
                };
                column.push(value);
            }
            for (key, column) in layouts
                .get(layout_ix)
                .float_indices
                .iter()
                .zip(float_columns.iter_mut())
            {
                column.push(
                    key.and_then(|ki| subfields.get(ki))
                        .map_or(f32::NAN, |range| parse_float(&field[range.clone()])),
                );
            }

            field_start = field_end + 1;
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
    n_alt: u32,
    policies: &GenotypePolicies,
    counters: &mut PolicyCounters,
    token: &mut [Option<u16>; 2],
) -> Result<GtCode, IoError> {
    decode_gt(gt, locus, n_alt, policies, counters, token)
}

/// Decodes one sample's GT subfield, applying the policies and filling
/// `token` with the raw allele indices (2 per sample).
///
/// `n_alt` is the record's ALT count (0 for `ALT=.`): an allele index above
/// it names an allele the record does not declare, which is malformed input
/// — not a policy question (doc_v25 V25-01).
fn decode_gt<S: LocusSource>(
    gt: &[u8],
    site: &S,
    n_alt: u32,
    policies: &GenotypePolicies,
    counters: &mut PolicyCounters,
    token: &mut [Option<u16>; 2],
) -> Result<GtCode, IoError> {
    // Fast path: exactly `[0|1|.][/|][0|1|.]`.
    if let Some((s0, s1, phased)) = fast_gt(gt) {
        // `ALT='.'` leaves only the REF allele: the fast path's `1` (ALT)
        // names an allele that does not exist.
        if n_alt == 0
            && (s0 == reed_core::genotype::STATE_ALT || s1 == reed_core::genotype::STATE_ALT)
        {
            return Err(IoError::MalformedRecord {
                locus: site.error_locus(),
                reason: "GT allele index 1 exceeds ALT count 0".to_string(),
            });
        }
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
        slots[n_alleles] = parse_allele(&gt[from..end], site, n_alt)?;
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
            reason: reason_str(&other).to_string(),
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
/// checked arithmetic — 01-review D14). The index must name a declared
/// allele: `0..=n_alt` (0 is REF; `n_alt == 0` means `ALT=.`), anything
/// above is malformed input (doc_v25 V25-01).
fn parse_allele<S: LocusSource>(
    token: &[u8],
    site: &S,
    n_alt: u32,
) -> Result<Option<u16>, IoError> {
    if token == b"." {
        return Ok(None);
    }
    if token.is_empty() || token.len() > 5 {
        return Err(IoError::MalformedRecord {
            locus: site.error_locus(),
            reason: "GT allele index is empty or too long".to_string(),
        });
    }
    let mut value: u16 = 0;
    for &b in token {
        if !b.is_ascii_digit() {
            return Err(IoError::MalformedRecord {
                locus: site.error_locus(),
                reason: "GT allele contains a non-digit".to_string(),
            });
        }
        // 5 digits can overflow u16 (99_999 > 65_535): reject the wrap.
        value = value
            .checked_mul(10)
            .and_then(|v| v.checked_add(u16::from(b - b'0')))
            .ok_or(IoError::MalformedRecord {
                locus: site.error_locus(),
                reason: "GT allele index exceeds u16".to_string(),
            })?;
    }
    if u32::from(value) > n_alt {
        return Err(IoError::MalformedRecord {
            locus: site.error_locus(),
            reason: format!("GT allele index {value} exceeds ALT count {n_alt}"),
        });
    }
    Ok(Some(value))
}

/// One integer FORMAT element under the strict VCF Integer grammar
/// ([`crate::vcf::info::parse_info_int`]: optional `+`/`-`, digits only),
/// mapped to the column's unsigned domain.
///
/// htslib stores FORMAT integers in `int32`; text outside `i32` (and a
/// negative value, which the `u32` column cannot carry) is missing. The old
/// saturating `u32` parser turned `4294967295` and every overflowing value
/// into [`INT_MISSING`] by coincidence and `5000000000` into a saturated
/// `u32::MAX` *value* — silently disagreeing with htslib (V22-R2-2).
/// Negative FORMAT Integers are a representation limit of the `u32` columns,
/// deferred to the signed-column round (V22-R2-8).
fn parse_field_int(text: &[u8]) -> Option<u32> {
    let value = crate::vcf::info::parse_info_int(text)?;
    let value = u32::try_from(value).ok()?;
    (value <= i32::MAX.unsigned_abs()).then_some(value)
}

/// Integer FORMAT value; `.`, empty, non-numeric, negative or outside
/// `i32` is [`INT_MISSING`].
#[must_use]
pub fn parse_int(text: &[u8]) -> u32 {
    parse_field_int(text).unwrap_or(INT_MISSING)
}

/// Comma-separated integer FORMAT value (`AD`): the sum of its readable
/// values. `INT_MISSING` when no value is present — htslib's text parse maps
/// `.`/unparsable tokens to missing and pads a short vector with vector-end,
/// which the comma split represents by the token simply ending.
///
/// The sum is accumulated in `u64` and clamped to `INT_MISSING - 1`, so a
/// vector of large (but individually valid) elements can never collide with
/// the missing sentinel: three `i32::MAX` elements sum to `INT_MISSING - 1`
/// rather than disappearing as missing (V22-R2-2).
#[must_use]
pub fn parse_int_sum(text: &[u8]) -> u32 {
    let mut sum = 0_u64;
    let mut any = false;
    for value in text.split(|&b| b == b',') {
        if let Some(v) = parse_field_int(value) {
            sum = sum.saturating_add(u64::from(v));
            any = true;
        }
    }
    if any { clamp_int_sum(sum) } else { INT_MISSING }
}

/// Clamps an `AD`-style sum to the largest representable non-missing value.
#[must_use]
pub fn clamp_int_sum(sum: u64) -> u32 {
    u32::try_from(sum.min(u64::from(INT_MISSING - 1))).unwrap_or(INT_MISSING - 1)
}

/// Float FORMAT value; `.`, empty or unparsable text is `NaN`.
///
/// The grammar is Rust's `f32` parser, which matches htslib's `strtod`
/// acceptance for the VCF Float regex: `inf`/`-inf` are values, `nan` is a
/// NaN value (V22-R2-3 — this is a compatibility choice, not a bug).
///
/// `NaN` is also the column's missing sentinel: a literal `nan` cannot be
/// told apart from a missing value downstream (V22-R2-7). Fixing that needs
/// `Option<f32>` columns or the BCF bit patterns and is deferred to its own
/// round; the BCF reader's `float_element_column` collapses its missing and
/// vector-end patterns the same way, so both input paths agree.
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
            decode_gt_text(gt, &locus, u32::MAX, &policies, &mut counters, &mut token)
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
        // Twelve ALTs keep `10/11` inside the allele domain (doc_v25 V25-01):
        // the multi-digit indices still collapse onto the ALT bit.
        let (out, counters) = scan(
            &mut scanner,
            b"c\t1\t.\tA\tG,T,A,C,G,T,A,C,G,T,A,C\t.\tPASS\t.\tGT\t1\t10/11\t0/.\n",
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

    /// doc_v25 V25-01 (F1): a GT allele index beyond the record's ALT count
    /// is malformed input, not a policy question. `0/9` on a single-ALT
    /// record used to decode as a plain ALT call under `CollapseAlt` and
    /// flow into every consumer — the `-Ob` encoder wrote `(9+1)<<1`,
    /// which htslib then refuses to read back.
    /// doc_v25 V25-02: GT is found by FORMAT key order, never assumed to be
    /// the first subfield. `DP:GT` used to hand `30` to the genotype parser
    /// in convert/pca (the scanner's layout was already correct; this pins
    /// the shared lookup).
    #[test]
    fn gt_position_follows_the_format_order() {
        let mut scanner = GenotypeScanner::new(1, &[], &[], permissive());
        let (out, _) = scan(&mut scanner, b"c\t1\t.\tA\tG\t.\tPASS\t.\tDP:GT\t30:0/1\n");
        assert_eq!(out.codes[0].state(), 0b0100, "DP:GT decodes the 0/1");
        let (out, _) = scan(
            &mut scanner,
            b"c\t1\t.\tA\tG,T\t.\tPASS\t.\tAD:DP:GQ:GT\t5,6:10:40:1|2\n",
        );
        assert_eq!(out.codes[0].state(), 0b0101, "GT in the last position");
    }

    #[test]
    fn gt_allele_beyond_the_alt_count_is_malformed() {
        // (ALT, GT) the scanner must refuse.
        let bad: &[(&[u8], &[u8])] = &[
            (b"G", b"0/9"),
            (b"G", b"1|2"),
            (b"G", b"9"),    // haploid out-of-domain call
            (b"G", b"9/9"),
            (b"G,T", b"0/3"),
            // ALT='.' leaves only the REF allele: allele 1 does not exist.
            (b".", b"0/1"),
        ];
        for (alt, gt) in bad {
            let mut text = String::from("c\t1\t.\tA\t");
            text.push_str(String::from_utf8_lossy(alt).as_ref());
            text.push_str("\t.\tPASS\t.\tGT\t");
            text.push_str(String::from_utf8_lossy(gt).as_ref());
            text.push('\n');
            let mut scanner = GenotypeScanner::new(1, &[], &[], permissive());
            let lr = line(text.as_bytes());
            let site = SiteView::parse(&lr).expect("invariant: well-formed");
            let mut counters = PolicyCounters::default();
            let Err(err) = scanner.scan(&lr, &site, &mut counters) else {
                panic!("ALT={alt:?} GT={gt:?} is out of the allele domain and must be refused");
            };
            assert!(
                matches!(err, IoError::MalformedRecord { .. }),
                "ALT={alt:?} GT={gt:?}: {err:?}"
            );
            let msg = err.to_string();
            let n_alt = if *alt == b"." {
                0
            } else {
                memchr::memchr_iter(b',', alt).count() + 1
            };
            assert!(
                msg.contains(&format!("exceeds ALT count {n_alt}")),
                "ALT={alt:?} GT={gt:?}: {msg}"
            );
        }
        // In-domain controls decode normally under the same policies.
        let good: &[(&[u8], &[u8])] = &[(b"G", b"0/1"), (b"G,T", b"1|2"), (b".", b"0/0")];
        for (alt, gt) in good {
            let mut text = String::from("c\t1\t.\tA\t");
            text.push_str(String::from_utf8_lossy(alt).as_ref());
            text.push_str("\t.\tPASS\t.\tGT\t");
            text.push_str(String::from_utf8_lossy(gt).as_ref());
            text.push('\n');
            let mut scanner = GenotypeScanner::new(1, &[], &[], permissive());
            let lr = line(text.as_bytes());
            let site = SiteView::parse(&lr).expect("invariant: well-formed");
            let mut counters = PolicyCounters::default();
            scanner
                .scan(&lr, &site, &mut counters)
                .unwrap_or_else(|e| panic!("ALT={alt:?} GT={gt:?} is in-domain: {e}"));
        }
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
    fn depth_column_falls_back_to_ad_sum() {
        use crate::vcf::format_layout::IntSpec;
        let mut scanner =
            GenotypeScanner::new_with_specs(3, &[IntSpec::depth()], &[], permissive());
        // No DP key: the samples' AD values are summed (Number=R, but every
        // present value counts — htslib does not truncate).
        let (out, _) = scan(
            &mut scanner,
            b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT:PL:AD\t0/1:.:10,20\t1/1:.:.,25\t0/0:.:.\n",
        );
        assert_eq!(out.int_columns[0], vec![30, 25, INT_MISSING]);
        // FORMAT with DP: DP wins and a missing DP stays missing (no AD
        // fallback), even though AD is present.
        let (out, _) = scan(
            &mut scanner,
            b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT:DP:AD\t0/1:.:10,20\t1/1:11:5,6\t0/0:0:1,2\n",
        );
        assert_eq!(out.int_columns[0], vec![INT_MISSING, 11, 0]);
    }

    #[test]
    fn parse_int_sum_matches_htslib_missing_rules() {
        assert_eq!(parse_int_sum(b"10,20"), 30);
        assert_eq!(parse_int_sum(b"10,.,20"), 30);
        assert_eq!(parse_int_sum(b"."), INT_MISSING);
        assert_eq!(parse_int_sum(b""), INT_MISSING);
        assert_eq!(parse_int_sum(b".,."), INT_MISSING);
        assert_eq!(parse_int_sum(b"7"), 7);
    }

    /// V22-R2-2/R2-5: one integer grammar (optional sign, htslib's `i32`
    /// domain). Saturating `u32` used to turn `4294967295` into the missing
    /// sentinel coincidentally and `5000000000` into a valid saturated value.
    #[test]
    fn parse_int_uses_the_htslib_i32_domain() {
        assert_eq!(parse_int(b"0"), 0);
        assert_eq!(parse_int(b"+50"), 50);
        assert_eq!(parse_int(b"0007"), 7);
        assert_eq!(parse_int(b"2147483647"), 2_147_483_647);
        assert_eq!(parse_int(b"2147483648"), INT_MISSING);
        assert_eq!(parse_int(b"4294967295"), INT_MISSING);
        assert_eq!(parse_int(b"5000000000"), INT_MISSING);
        // Negative FORMAT integers need signed columns (V22-R2-8).
        assert_eq!(parse_int(b"-5"), INT_MISSING);
        assert_eq!(parse_int(b"+"), INT_MISSING);
        assert_eq!(parse_int(b"-"), INT_MISSING);
        assert_eq!(parse_int(b"1.5"), INT_MISSING);
        assert_eq!(parse_int(b" 5"), INT_MISSING);
        assert_eq!(parse_int(b"5x"), INT_MISSING);
    }

    /// V22-R2-2: the `AD` sum clamps to the largest non-missing value and
    /// skips out-of-domain elements instead of either wrapping or vanishing.
    #[test]
    fn int_sum_never_collides_with_the_missing_sentinel() {
        // 2 × i32::MAX is exactly INT_MISSING - 1.
        assert_eq!(parse_int_sum(b"2147483647,2147483647"), INT_MISSING - 1);
        assert_eq!(
            parse_int_sum(b"2147483647,2147483647,2147483647"),
            INT_MISSING - 1
        );
        // Out-of-domain and negative elements are missing elements.
        assert_eq!(parse_int_sum(b"2147483648,5"), 5);
        assert_eq!(parse_int_sum(b"-5,7"), 7);
        assert_eq!(parse_int_sum(b"+5,+0"), 5);
        assert_eq!(clamp_int_sum(u64::MAX), INT_MISSING - 1);
        assert_eq!(clamp_int_sum(0), 0);
    }

    /// V22-R2-3: `inf`/`nan` are accepted exactly like htslib's `strtod`
    /// (this pins the decision; the GQ example in the review used the
    /// integer parser, which rejects both). V22-R2-7: `NaN` is also the
    /// missing sentinel, so a literal `nan` is indistinguishable from a
    /// missing value.
    #[test]
    fn parse_float_accepts_htslib_inf_and_nan() {
        assert_eq!(parse_float(b"1.5"), 1.5);
        assert!(parse_float(b"inf").is_infinite() && parse_float(b"inf") > 0.0);
        assert!(parse_float(b"-inf").is_infinite() && parse_float(b"-inf") < 0.0);
        assert!(parse_float(b"nan").is_nan());
        assert!(parse_float(b"NaN").is_nan());
        // Not numbers under the VCF Float grammar.
        assert!(parse_float(b"0x10").is_nan());
        assert!(parse_float(b"1,5").is_nan());
        assert!(parse_float(b"").is_nan());
        assert!(parse_float(b".").is_nan());
        // The GQ path is Integer, so `inf` is missing there, not a value.
        assert_eq!(parse_int(b"inf"), INT_MISSING);
    }

    /// V22-R2-1: positions at and past index 256 are ordinary subfields.
    /// The `u8` layout index used to make GQ/AD/GT silently absent there.
    #[test]
    fn keys_past_index_255_decode() {
        use crate::vcf::format_layout::IntSpec;
        let mut keys: Vec<String> = (0..256).map(|i| format!("K{i}")).collect();
        keys.extend(["GQ".to_string(), "AD".to_string(), "GT".to_string()]);
        let format = keys.join(":");
        let mut values = vec!["."; 256];
        values.extend(["7", "10,20", "0/1"]);
        let sample = values.join(":");
        let text = format!("c\t1\t.\tA\tG\t.\tPASS\t.\t{format}\t{sample}\n");
        let mut scanner = GenotypeScanner::new_with_specs(
            1,
            &[IntSpec::scalar(b"GQ"), IntSpec::depth()],
            &[],
            permissive(),
        );
        let (out, _) = scan(&mut scanner, text.as_bytes());
        assert_eq!(out.int_columns[0], vec![7], "GQ at index 256");
        assert_eq!(out.int_columns[1], vec![30], "AD sum at index 257");
        assert!(out.has_gt, "GT at index 258");
        assert_eq!(out.codes[0], GtCode::from_states(0, 1, false, false));
    }

    #[test]
    fn sites_only_and_trailing_empty_columns_decode() {
        // 8 columns, no samples: nothing to decode, no phantom sample column.
        let mut scanner = GenotypeScanner::new(0, &[], &[], permissive());
        let (out, _) = scan(
            &mut scanner,
            b"c\t1\t.\tA\tG\t.\tPASS\tDP=5\n",
        );
        assert_eq!(out.codes, []);
        assert!(!out.has_gt);
        // A record ending in one (or more) empty sample columns keeps the
        // declared samples and ignores the extras, as htslib does.
        let mut scanner = GenotypeScanner::new(2, &[b"DP"], &[], permissive());
        let (out, _) = scan(
            &mut scanner,
            b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT:DP\t0/1:9\t1/1:8\t\t\n",
        );
        assert_eq!(out.codes.len(), 2);
        assert_eq!(out.int_columns[0], vec![9, 8]);
        // A non-empty extra column is still malformed.
        let lr = line(b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT:DP\t0/1:9\t1/1:8\t0/0:7\n");
        let site = SiteView::parse(&lr).expect("invariant: well-formed");
        let mut counters = PolicyCounters::default();
        let rejected = scanner.scan(&lr, &site, &mut counters);
        assert!(matches!(rejected, Err(IoError::MalformedRecord { .. })));

        // A FORMAT column with no sample columns is malformed when the
        // header declares samples (bcftools: "FORMAT column with no sample
        // columns"), not one phantom missing sample (doc_v22 R22-3).
        let mut scanner = GenotypeScanner::new(1, &[b"DP"], &[], permissive());
        let lr = line(b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT\n");
        let site = SiteView::parse(&lr).expect("invariant: well-formed");
        let mut counters = PolicyCounters::default();
        let rejected = scanner.scan(&lr, &site, &mut counters);
        assert!(matches!(rejected, Err(IoError::MalformedRecord { .. })));
    }

    #[test]
    fn gt_absent_from_format_is_all_missing() {
        let mut scanner = GenotypeScanner::new(3, &[b"DP"], &[], permissive());
        let (out, _) = scan(&mut scanner, b"c\t1\t.\tA\tG\t.\tPASS\t.\tDP\t5\t6\t7\n");
        assert!(out.codes.iter().all(|c| c.is_fully_missing()));
        assert_eq!(out.int_columns[0], vec![5, 6, 7]);
    }

    #[test]
    fn r2_empty_format_column_does_not_panic() {
        // doc_v20 R2: a fresh layout cache seeds its last-hit entry with
        // `(empty, usize::MAX)`, so a record with an *empty* FORMAT column
        // matched the sentinel and `get(usize::MAX)` panicked the scanner.
        // The empty FORMAT interns like any GT-less layout: every sample
        // decodes as missing, `has_gt` is false, and the cached index stays
        // valid on a repeat of the same (empty) FORMAT.
        let mut scanner = GenotypeScanner::new(1, &[b"DP"], &[], permissive());
        let text = b"c\t1\t.\tA\tG\t.\tPASS\t.\t\t0/0\n";
        let (out, _) = scan(&mut scanner, text);
        assert!(!out.has_gt);
        assert!(out.codes.iter().all(|c| c.is_fully_missing()));
        assert_eq!(out.int_columns[0], vec![INT_MISSING]);
        let (out, _) = scan(&mut scanner, text);
        assert!(!out.has_gt && out.codes[0].is_fully_missing());
    }

    #[test]
    fn column_count_mismatch_is_an_error() {
        let mut scanner = GenotypeScanner::new(3, &[], &[], permissive());
        let lr = line(b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\n");
        let site = SiteView::parse(&lr).expect("invariant: well-formed");
        let mut counters = PolicyCounters::default();
        let err = scanner.scan(&lr, &site, &mut counters).unwrap_err();
        // doc_v23 V23-29: the diagnostic carries both counts.
        let text = err.to_string();
        assert!(matches!(err, IoError::MalformedRecord { .. }));
        assert!(text.contains("found 1"), "{text}");
        assert!(text.contains("expected 3"), "{text}");
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
