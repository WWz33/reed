//! Single-pass, multi-field sample-column scanner (`docs/architecture-spec.md`
//! §6.4, decision D1). For K requested FORMAT fields the sample region is
//! walked exactly once; GT has a 3-byte fast path that writes bit-plane bits
//! directly. Correctness is pinned by a differential test against noodles
//! (`tests/sample_scan_vs_noodles.rs`) and gold parity.

use std::collections::HashMap;
use std::ops::Range;

use reed_core::genotype::{
    write_genotype, GenotypePolicies, PlaneBlock, Ploidy, PolicyCounters,
};

use crate::error::IoError;
use crate::vcf::line::LineRef;

/// Sentinel for a missing integer sample value (`docs/architecture-spec.md`
/// §6.4): real DP/GQ values never reach `u32::MAX`.
pub const INT_MISSING: u32 = u32::MAX;

/// Subfield index layout of one FORMAT string, cached by its bytes.
#[derive(Debug, Clone)]
struct Layout {
    /// Index of `GT` among the colon-separated keys (None when absent).
    gt_index: Option<u8>,
    /// Index of each requested integer key, aligned with the scanner's key
    /// list (None = key absent from this FORMAT).
    int_indices: Vec<Option<u8>>,
}

/// Worker-local scanner: format-layout cache plus reusable scratch.
///
/// # ponytail: clear-on-overflow format cache (64 entries); a real LRU only
/// if a corpus with adversarial FORMAT variety ever shows up.
pub struct SampleScanner {
    int_keys: Vec<Box<[u8]>>,
    layouts: HashMap<Box<[u8]>, Layout>,
    subfields: Vec<Range<usize>>,
    alleles: Vec<Option<u16>>,
}

impl SampleScanner {
    /// Creates a scanner for GT plus the named integer keys (e.g.
    /// `&[b"DP", b"GQ"]`), in output order.
    #[must_use]
    pub fn new(int_keys: &[&[u8]]) -> Self {
        Self {
            int_keys: int_keys.iter().map(|k| k.to_vec().into_boxed_slice()).collect(),
            layouts: HashMap::new(),
            subfields: Vec::new(),
            alleles: Vec::new(),
        }
    }

    /// Number of integer columns the scanner fills.
    #[must_use]
    pub fn n_int_columns(&self) -> usize {
        self.int_keys.len()
    }

    /// Scans one record.
    ///
    /// Writes the genotype planes of row `site_idx` (one per sample) and
    /// appends one integer value per sample per requested key to `int_out`
    /// ([`INT_MISSING`] when the key or the value is absent).
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] for column-count mismatches and unparsable
    /// values; [`IoError::PolicyRejected`] for genotypes the declared
    /// policies forbid (polyploid, multi-allelic under `Reject`, haploid
    /// under `RequireDiploid`).
    #[allow(clippy::too_many_arguments)]
    pub fn scan_record(
        &mut self,
        line: &LineRef<'_>,
        site_idx: usize,
        pos1: u64,
        planes: &mut PlaneBlock,
        int_out: &mut [Vec<u32>],
        policies: &GenotypePolicies,
        counters: &mut PolicyCounters,
    ) -> Result<(), IoError> {
        debug_assert_eq!(int_out.len(), self.int_keys.len(), "invariant: one output column per key");
        let n_samples = planes.n_samples();
        if line.n_sample_cols() != n_samples {
            return Err(self.malformed(line, pos1, "sample column count disagrees with the header"));
        }

        let layout = self.layout_for(line.format_text())?;

        // Sample fields, one tab-delimited column each.
        let samples = line.samples_text();
        let mut field_start = 0;
        for sample_idx in 0..n_samples {
            let field_end = memchr::memchr(b'\t', &samples[field_start..])
                .map_or(samples.len(), |rel| field_start + rel);
            let field = &samples[field_start..field_end];
            self.split_subfields(field);

            // Genotype.
            let gt_bytes = layout
                .gt_index
                .and_then(|gi| self.subfields.get(usize::from(gi)))
                .map(|r| &field[r.clone()]);
            self.write_gt(site_idx, sample_idx, gt_bytes, planes, policies, counters, line, pos1)?;

            // Integer columns.
            for (key_index, value_out) in int_out.iter_mut().enumerate() {
                let value = layout.int_indices[key_index]
                    .and_then(|ki| self.subfields.get(usize::from(ki)))
                    .map(|r| &field[r.clone()])
                    .map_or(INT_MISSING, parse_int);
                value_out.push(value);
            }

            field_start = field_end + 1;
        }
        Ok(())
    }

    /// Scans one record with a callback: per sample, the GT subfield bytes
    /// (`None` when the FORMAT has no GT) and one parsed integer per
    /// requested key (`None` when the key or the value is absent) are
    /// handed to `f` in sample order. The filter's per-site pipeline runs
    /// on this — no plane materialization, single pass.
    ///
    /// # Errors
    ///
    /// [`IoError::MalformedRecord`] for column-count mismatches and
    /// unparsable integer values.
    pub fn scan_fields<F>(&mut self, line: &LineRef<'_>, mut f: F) -> Result<(), IoError>
    where
        F: FnMut(usize, Option<&[u8]>, &[Option<&[u8]>]),
    {
        let layout = self.layout_for(line.format_text())?;
        let samples = line.samples_text();
        let mut field_start = 0;
        let mut sample_idx = 0_usize;
        while field_start < samples.len() {
            let field_end = memchr::memchr(b'\t', &samples[field_start..])
                .map_or(samples.len(), |rel| field_start + rel);
            let field = &samples[field_start..field_end];
            self.split_subfields(field);
            let gt = layout
                .gt_index
                .and_then(|gi| self.subfields.get(usize::from(gi)))
                .map(|r| &field[r.clone()]);
            let ints: Vec<Option<&[u8]>> = layout
                .int_indices
                .iter()
                .map(|ki| ki.and_then(|i| self.subfields.get(usize::from(i))).map(|r| &field[r.clone()]))
                .collect();
            f(sample_idx, gt, &ints);
            sample_idx += 1;
            field_start = field_end + 1;
        }
        Ok(())
    }

    /// Cached FORMAT layout (parsed on first sight of each distinct FORMAT).
    fn layout_for(&mut self, format: &[u8]) -> Result<Layout, IoError> {
        if self.layouts.len() > 64 {
            self.layouts.clear();
        }
        if let Some(layout) = self.layouts.get(format) {
            return Ok(layout.clone());
        }
        let keys: Vec<&[u8]> = format.split(|&b| b == b':').collect();
        let index_of = |key: &[u8]| -> Option<u8> {
            keys.iter()
                .position(|k| *k == key)
                .and_then(|i| u8::try_from(i).ok())
        };
        let layout = Layout {
            gt_index: index_of(b"GT"),
            int_indices: self.int_keys.iter().map(|k| index_of(k)).collect(),
        };
        self.layouts.insert(format.to_vec().into_boxed_slice(), layout.clone());
        Ok(layout)
    }

    /// Splits one sample field into colon-separated subfield ranges.
    fn split_subfields(&mut self, field: &[u8]) {
        self.subfields.clear();
        let mut from = 0;
        loop {
            match memchr::memchr(b':', &field[from..]) {
                Some(rel) => {
                    self.subfields.push(from..from + rel);
                    from += rel + 1;
                }
                None => {
                    self.subfields.push(from..field.len());
                    break;
                }
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn write_gt(
        &mut self,
        site_idx: usize,
        sample_idx: usize,
        gt: Option<&[u8]>,
        planes: &mut PlaneBlock,
        policies: &GenotypePolicies,
        counters: &mut PolicyCounters,
        line: &LineRef<'_>,
        pos1: u64,
    ) -> Result<(), IoError> {
        let Some(gt) = gt else {
            // GT absent from this record's FORMAT: every slot missing.
            planes.set_slot(site_idx, sample_idx, 0, false, true);
            planes.set_slot(site_idx, sample_idx, 1, false, true);
            planes.set_phased(site_idx, sample_idx, false);
            return Ok(());
        };
        // Fast path: exactly `[0|1|.][/|][0|1|.]`.
        if let Some((alt0, miss0, alt1, miss1, phased)) = fast_gt(gt) {
            planes.set_slot(site_idx, sample_idx, 0, alt0, miss0);
            planes.set_slot(site_idx, sample_idx, 1, alt1, miss1);
            planes.set_phased(site_idx, sample_idx, phased);
            return Ok(());
        }
        // Slow path.
        if gt == b"." {
            // Haploid missing: both slots missing, no policy involved.
            planes.set_slot(site_idx, sample_idx, 0, false, true);
            planes.set_slot(site_idx, sample_idx, 1, false, true);
            planes.set_phased(site_idx, sample_idx, false);
            return Ok(());
        }
        self.alleles.clear();
        let mut phased = false;
        let mut from = 0;
        while from <= gt.len() {
            let end = memchr::memchr2(b'/', b'|', &gt[from..])
                .map_or(gt.len(), |rel| from + rel);
            self.alleles.push(parse_allele(&gt[from..end], line, pos1)?);
            if end == gt.len() {
                break;
            }
            phased |= gt[end] == b'|';
            from = end + 1;
        }
        let mut write = |planes: &mut PlaneBlock, slots: [Option<u16>; 2], ploidy: Ploidy| {
            write_genotype(planes, site_idx, sample_idx, slots, ploidy, phased, policies, counters)
                .map_err(|err| match err {
                    reed_core::CoreError::PolicyRejected(reason) => {
                        IoError::PolicyRejected { locus: locus_of(line, pos1), reason }
                    }
                    reed_core::CoreError::Unimplemented(reason) => IoError::UnsupportedInput(reason),
                    other => IoError::MalformedRecord {
                        locus: locus_of(line, pos1),
                        reason: reason_str(&other),
                    },
                })
        };
        match self.alleles.as_slice() {
            [allele] => write(planes, [*allele, None], Ploidy::Haploid),
            [a, b] => write(planes, [*a, *b], Ploidy::Diploid),
            _ => Err(self.rejected(line, pos1, "polyploid genotype (more than two alleles)")),
        }
    }

    fn locus(&self, line: &LineRef<'_>, pos1: u64) -> crate::Locus {
        crate::Locus {
            contig_name: String::from_utf8_lossy(line.chrom()).into_owned(),
            pos1,
        }
    }

    fn malformed(&self, line: &LineRef<'_>, pos1: u64, reason: &'static str) -> IoError {
        IoError::MalformedRecord { locus: self.locus(line, pos1), reason }
    }

    fn rejected(&self, line: &LineRef<'_>, pos1: u64, reason: &'static str) -> IoError {
        IoError::PolicyRejected { locus: self.locus(line, pos1), reason }
    }
}

/// The 3-byte GT fast path: returns plane bits when the genotype is exactly
/// one of `[0|1|.] [/|] [0|1|.]`.
fn fast_gt(gt: &[u8]) -> Option<(bool, bool, bool, bool, bool)> {
    if gt.len() != 3 {
        return None;
    }
    let slot = |b: u8| match b {
        b'0' => Some((false, false)),
        b'1' => Some((true, false)),
        b'.' => Some((false, true)),
        _ => None,
    };
    let (alt0, miss0) = slot(gt[0])?;
    let (alt1, miss1) = slot(gt[2])?;
    if gt[1] != b'/' && gt[1] != b'|' {
        return None;
    }
    Some((alt0, miss0, alt1, miss1, gt[1] == b'|'))
}

/// One allele token: `None` for `.`, else the index (multi-digit capable).
fn parse_allele(token: &[u8], line: &LineRef<'_>, pos1: u64) -> Result<Option<u16>, IoError> {
    if token == b"." {
        return Ok(None);
    }
    if token.is_empty() || token.len() > 5 {
        return Err(IoError::MalformedRecord {
            locus: locus_of(line, pos1),
            reason: "GT allele index is empty or too long",
        });
    }
    let mut value: u16 = 0;
    for &b in token {
        if !b.is_ascii_digit() {
            return Err(IoError::MalformedRecord {
                locus: locus_of(line, pos1),
                reason: "GT allele contains a non-digit",
            });
        }
        value = value * 10 + u16::from(b - b'0');
    }
    Ok(Some(value))
}

fn parse_int(text: &[u8]) -> u32 {
    if text == b"." || text.is_empty() {
        return INT_MISSING;
    }
    let mut value: u32 = 0;
    for &b in text {
        if !b.is_ascii_digit() {
            // Negative or malformed: count as missing (matches the legacy
            // decoder, which yields None for non-integer values).
            return INT_MISSING;
        }
        value = value.saturating_mul(10).saturating_add(u32::from(b - b'0'));
    }
    value
}

/// Short fallback text for `CoreError` variants without a static reason.
fn reason_str(err: &reed_core::CoreError) -> &'static str {
    match err {
        reed_core::CoreError::PolicyRejected(reason) => reason,
        _ => "genotype rejected by policy",
    }
}

fn locus_of(line: &LineRef<'_>, pos1: u64) -> crate::Locus {
    crate::Locus {
        contig_name: String::from_utf8_lossy(line.chrom()).into_owned(),
        pos1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use reed_core::genotype::{MultiallelicPolicy, PloidyPolicy};

    fn scan_line(
        scanner: &mut SampleScanner,
        text: &[u8],
        policies: &GenotypePolicies,
    ) -> (PlaneBlock, Vec<Vec<u32>>) {
        let line = LineRef::parse(text).expect("invariant: well-formed test line");
        let mut planes = PlaneBlock::new(1, 3, true);
        let mut int_out = vec![Vec::new(); scanner.n_int_columns()];
        let mut counters = PolicyCounters::default();
        scanner
            .scan_record(&line, 0, 1, &mut planes, &mut int_out, policies, &mut counters)
            .expect("invariant: test line scans");
        (planes, int_out)
    }

    fn policies() -> GenotypePolicies {
        GenotypePolicies {
            ploidy: PloidyPolicy::HaploidAsHalfCall,
            multiallelic: MultiallelicPolicy::CollapseAlt,
            ..GenotypePolicies::default()
        }
    }

    fn bit(planes: &PlaneBlock, plane: usize, sample: usize) -> bool {
        let words = planes.plane(0, plane);
        words[sample / 64] >> (sample % 64) & 1 == 1
    }

    #[test]
    fn fast_path_and_phasing() {
        let mut scanner = SampleScanner::new(&[]);
        let (planes, _) = scan_line(
            &mut scanner,
            b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\t1|0\t./.\n",
            &GenotypePolicies::default(),
        );
        // s0 "0/1": slot0 called REF, slot1 ALT (plane 2).
        assert!(!bit(&planes, 0, 0) && !bit(&planes, 1, 0) && bit(&planes, 2, 0) && !bit(&planes, 3, 0));
        // s1 "1|0": slot0 ALT, slot1 called REF, phased bit set.
        assert!(bit(&planes, 0, 1) && !bit(&planes, 1, 1) && !bit(&planes, 2, 1) && !bit(&planes, 3, 1));
        assert!(bit(&planes, 4, 1));
        // s2 "./.": both slots missing (planes 1 and 3).
        assert!(bit(&planes, 1, 2) && bit(&planes, 3, 2));
    }

    #[test]
    fn haploid_half_call_and_multidigit() {
        let mut scanner = SampleScanner::new(&[]);
        let (planes, _) = scan_line(&mut scanner, b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT\t1\t10/11\t0/.\n", &policies());
        // s0: haploid ALT → slot0 alt, slot1 missing (half call).
        assert!(bit(&planes, 0, 0) && bit(&planes, 3, 0));
        // s1: multi-digit alleles collapse onto ALT.
        assert!(bit(&planes, 0, 1) && bit(&planes, 2, 1));
        // s2: half call 0/. → slot0 ref called, slot1 missing.
        assert!(!bit(&planes, 0, 2) && bit(&planes, 3, 2));
    }

    #[test]
    fn int_columns_with_missing_variants() {
        let mut scanner = SampleScanner::new(&[b"DP", b"GQ"]);
        let (planes, cols) = scan_line(
            &mut scanner,
            b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT:DP:GQ\t0/1:9:40\t1/1:.:.\t1/1\n",
            &policies(),
        );
        assert_eq!(cols[0], vec![9, INT_MISSING, INT_MISSING]);
        assert_eq!(cols[1], vec![40, INT_MISSING, INT_MISSING]);
        let _ = planes;
    }

    #[test]
    fn gt_absent_from_format_yields_all_missing() {
        let mut scanner = SampleScanner::new(&[b"DP"]);
        let (planes, cols) = scan_line(
            &mut scanner,
            b"c\t1\t.\tA\tG\t.\tPASS\t.\tDP\t5\t6\t7\n",
            &GenotypePolicies::default(),
        );
        for sample in 0..3 {
            assert!(bit(&planes, 1, sample) && bit(&planes, 3, sample), "sample {sample} must be all-missing");
        }
        assert_eq!(cols[0], vec![5, 6, 7]);
    }

    #[test]
    fn column_count_mismatch_is_an_error() {
        let mut scanner = SampleScanner::new(&[]);
        let line = LineRef::parse(b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\n").expect("invariant");
        let mut planes = PlaneBlock::new(1, 3, false);
        let mut counters = PolicyCounters::default();
        let err = scanner
            .scan_record(&line, 0, 1, &mut planes, &mut [], &policies(), &mut counters)
            .unwrap_err();
        assert!(matches!(err, IoError::MalformedRecord { .. }));
    }

    #[test]
    fn polyploid_is_rejected() {
        let mut scanner = SampleScanner::new(&[]);
        let line = LineRef::parse(b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/1/2\t0/0\t0/0\n").expect("invariant");
        let mut planes = PlaneBlock::new(1, 3, false);
        let mut counters = PolicyCounters::default();
        let err = scanner
            .scan_record(&line, 0, 1, &mut planes, &mut [], &policies(), &mut counters)
            .unwrap_err();
        assert!(matches!(err, IoError::PolicyRejected { .. }));
    }

    #[test]
    fn require_diploid_rejects_haploid() {
        let mut scanner = SampleScanner::new(&[]);
        let line = LineRef::parse(b"c\t1\t.\tA\tG\t.\tPASS\t.\tGT\t1\t0/0\t0/0\n").expect("invariant");
        let mut planes = PlaneBlock::new(1, 3, false);
        let mut counters = PolicyCounters::default();
        let err = scanner
            .scan_record(&line, 0, 1, &mut planes, &mut [], &GenotypePolicies::default(), &mut counters)
            .unwrap_err();
        assert!(matches!(err, IoError::PolicyRejected { .. }));
    }
}
