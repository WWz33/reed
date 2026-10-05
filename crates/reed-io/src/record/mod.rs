//! Record-level VCF/BCF transport: stream, decode, edit, count.
//!
//! Everything here works on `vcf::variant::RecordBuf` (noodles' in-memory
//! record with full INFO/FORMAT fidelity):
//!
//! - [`for_each_record`]: stream every record of a file (dispatched by
//!   extension: `.bcf`, `.vcf.gz`/`.gz`, `.vcf`) as a `RecordBuf`.
//! - [`record_header`]: read only the header of a file (the full
//!   `vcf::Header`, for output-side work on the consumer side).
//! - [`decode`]: `RecordBuf` -> [`RecordView`], the FieldSet-decoded site
//!   view (GT/DP/GQ plus the four GATK INFO tags).
//! - [`apply`]: the surgical edit applied after a site passed every filter:
//!   prune sample columns, mask GT to `./.`, drop stale INFO tags, write the
//!   recomputed counts. Every other field passes through untouched.
//! - [`filter_header`]: prune the #CHROM sample-name list to the kept
//!   columns.
//! - [`bcf::to_record_buf`]: the `bcf::Record` -> `RecordBuf` conversion
//!   (used by [`for_each_record`]; public for callers that read BCF records
//!   themselves).
//!
//! The mask rule ([`GtMask`]) travels with the record layer: it is expressed
//! purely in DP/GQ values and carries no policy context.

pub mod bcf;
pub mod vcf;

pub use bcf::to_record_buf;
pub use vcf::{apply, count_missing, decode, filter_header, STALE_INFO};

use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use noodles::bgzf;
use noodles::vcf::variant::RecordBuf;

use reed_core::estimators::site_counts::count_site;

use crate::ReadError;

/// Exact per-site allele and genotype counts (the canonical definition lives
/// in reed-core's `estimators::site_counts`; the record layer only hands it
/// genotype vectors).
pub use reed_core::estimators::site_counts::SiteCounts;

/// The FieldSet-decoded view of one site.
///
/// gt length = 2 x n_samples with -1 for a missing allele; dp/gq length =
/// n_samples with None when that column is missing for that sample.
/// `masked` is filled by [`RecordView::mask`]: which samples this rule set
/// to missing.
#[derive(Clone, Debug, Default)]
pub struct RecordView {
    /// Contig (CHROM) name.
    pub contig: String,
    /// 1-based VCF position (POS).
    pub pos: u64,
    /// None = the QUAL column is '.'.
    pub qual: Option<f32>,
    /// Whether the FILTER column is exactly PASS ('.' and any other filter
    /// yield false).
    pub is_pass: bool,
    /// REF bases.
    pub ref_allele: String,
    /// ALT alleles.
    pub alt_alleles: Vec<String>,
    /// The four GATK tags from the raw INFO (None = missing; a sub-rule
    /// whose tag is missing passes, matching bcftools).
    pub info_qd: Option<f64>,
    /// Fisher-strand bias (FS) from the raw INFO (None = missing).
    pub info_fs: Option<f64>,
    /// RMS mapping quality (MQ) from the raw INFO (None = missing).
    pub info_mq: Option<f64>,
    /// Symmetric odds ratio (SOR) from the raw INFO (None = missing).
    pub info_sor: Option<f64>,
    /// 2 x n_samples; -1 = missing allele. A half-call (one side missing)
    /// is stored as (a, -1); a haploid call as (a, -1).
    pub gt: Vec<i8>,
    /// Per-sample read depth (None = the DP column is missing for that
    /// sample).
    pub dp: Vec<Option<u32>>,
    /// Per-sample genotype quality (None = the GQ column is missing for that
    /// sample).
    pub gq: Vec<Option<u32>>,
    /// Number of sample columns.
    pub n_samples: usize,
    /// Filled by [`RecordView::mask`]: the samples this rule set to
    /// missing (including already-missing ones).
    pub masked: Vec<bool>,
}

/// Genotype-level mask rule: kept samples whose DP/GQ fall outside the
/// thresholds become `./.`.
///
/// Samples with missing DP/GQ are not masked (a bcftools numeric comparison
/// is false on missing). The mask changes GT only and keeps the original
/// DP/GQ values.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GtMask {
    /// Mask kept samples whose DP is below this (None = no lower bound).
    pub dp_min: Option<u32>,
    /// Mask kept samples whose DP is above this (None = no upper bound).
    pub dp_max: Option<u32>,
    /// Mask kept samples whose GQ is below this (None = no bound).
    pub gq_min: Option<u32>,
}

impl GtMask {
    /// Whether one sample is masked. A missing value never fails (see the
    /// type doc).
    pub fn fails(&self, dp: Option<u32>, gq: Option<u32>) -> bool {
        let dp_fail = match (self.dp_min, dp) {
            (Some(m), Some(v)) => v < m,
            _ => false,
        };
        let dp_hi_fail = match (self.dp_max, dp) {
            (Some(m), Some(v)) => v > m,
            _ => false,
        };
        let gq_fail = match (self.gq_min, gq) {
            (Some(m), Some(v)) => v < m,
            _ => false,
        };
        dp_fail || dp_hi_fail || gq_fail
    }
}

/// The set of INFO tags recomputed and overwritten after masking.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Recompute {
    /// Recompute the AC (allele count) INFO tag.
    pub ac: bool,
    /// Recompute the AN (allele number) INFO tag.
    pub an: bool,
    /// Recompute the NS (number of samples with data) INFO tag.
    pub ns: bool,
    /// Recompute the MAF (minor allele frequency) INFO tag.
    pub maf: bool,
    /// Recompute the HWE (exact-test p-value) INFO tag.
    pub hwe: bool,
}

impl Recompute {
    /// Recompute every supported tag (the conf default).
    pub const ALL: Self = Self { ac: true, an: true, ns: true, maf: true, hwe: true };
}

impl RecordView {
    /// Biallelic SNP test: matches bcftools `-v snps` intersected with
    /// plink2 `--snps-only just-acgt` (REF/ALT single base, ACGT,
    /// ALT != ".", REF != ALT; symbolic alleles and `*` fail the length and
    /// character checks).
    pub fn is_biallelic_snp(&self) -> bool {
        if self.ref_allele.len() != 1 || self.alt_alleles.len() != 1 {
            return false;
        }
        let r = self.ref_allele.as_bytes()[0].to_ascii_uppercase();
        let a = self.alt_alleles[0].as_bytes()[0].to_ascii_uppercase();
        matches!(r, b'A' | b'C' | b'G' | b'T')
            && matches!(a, b'A' | b'C' | b'G' | b'T')
            && r != a
    }

    /// Whether sample i is kept (the keep table).
    #[inline]
    fn kept(&self, keep: &[bool], i: usize) -> bool {
        keep.get(i).copied().unwrap_or(true)
    }

    /// Genotype mask: kept samples whose DP/GQ fall outside the rule
    /// thresholds become `./.`.
    ///
    /// Missing DP/GQ never masks (a bcftools numeric comparison is false on
    /// missing). Only [`RecordView::gt`] changes; dp/gq keep their original
    /// values. Recompute [`SiteCounts`] after masking.
    pub fn mask(&mut self, keep: &[bool], rule: &Option<GtMask>) {
        self.masked = vec![false; self.n_samples];
        let Some(rule) = rule else { return };
        for i in 0..self.n_samples {
            if !self.kept(keep, i) {
                continue;
            }
            if rule.fails(self.dp[i], self.gq[i]) {
                self.gt[2 * i] = -1;
                self.gt[2 * i + 1] = -1;
                self.masked[i] = true;
            }
        }
    }

    /// Counts over retained samples after masking.
    ///
    /// Delegates to `reed_core::estimators::site_counts::count_site`, the
    /// canonical counting primitive: dropped samples never enter any
    /// denominator (their output column is removed too), half-calls count
    /// their one present allele.
    pub fn counts(&self, retain: &[bool]) -> SiteCounts {
        count_site(&self.gt, Some(retain))
    }
}

/// Streams every record of `path` to `f` as a [`RecordBuf`], dispatched by
/// file extension: `.bcf`, `.vcf.gz` (or any `.gz`/`.bgz`), and plain
/// `.vcf`.
///
/// The reader consumes the header itself before the first record (BCF
/// records cannot be decoded before the header block is read; a VCF reader
/// that has not consumed its header text would parse header lines as
/// records). The BCF reader is built as
/// `bcf::io::Reader::from(bgzf::io::Reader::new(..))`: exactly one BGZF
/// layer over the file.
///
/// Reader failures (I/O, unsupported format, per-record decode) convert
/// into `E` through `From<ReadError>`; the callback's own error passes
/// through unchanged.
pub fn for_each_record<F, E>(path: &Path, mut f: F) -> Result<(), E>
where
    F: FnMut(RecordBuf) -> Result<(), E>,
    E: From<ReadError>,
{
    let name = path
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| ReadError::UnsupportedFormat(path.display().to_string()))?
        .to_string();
    let file = File::open(path).map_err(ReadError::Io)?;
    if name.ends_with(".bcf") {
        let mut reader =
            noodles::bcf::io::Reader::from(bgzf::io::Reader::new(BufReader::new(file)));
        let header = reader.read_header().map_err(ReadError::Io)?;
        let mut rec = noodles::bcf::Record::default();
        loop {
            let n = reader.read_record(&mut rec).map_err(ReadError::Io)?;
            if n == 0 {
                return Ok(());
            }
            let buf = to_record_buf(&rec, &header)?;
            f(buf)?;
        }
    } else if name.ends_with(".gz") || name.ends_with(".bgz") {
        let mut reader =
            noodles::vcf::io::Reader::new(bgzf::io::Reader::new(BufReader::new(file)));
        let header = reader.read_header().map_err(ReadError::Io)?;
        for result in reader.record_bufs(&header) {
            f(result.map_err(ReadError::Io)?)?;
        }
        Ok(())
    } else if name.ends_with(".vcf") {
        let mut reader = noodles::vcf::io::Reader::new(BufReader::new(file));
        let header = reader.read_header().map_err(ReadError::Io)?;
        for result in reader.record_bufs(&header) {
            f(result.map_err(ReadError::Io)?)?;
        }
        Ok(())
    } else {
        Err(ReadError::UnsupportedFormat(name).into())
    }
}

/// Reads only the header of a VCF/BCF file (record layer), dispatched by the
/// same extensions as [`for_each_record`].
///
/// The window layer's [`crate::header_info`] carries only a summary; the
/// record layer's consumers need the full header for output-side work
/// (sample-name pruning, INFO definitions, writing the filtered header).
pub fn record_header(path: &Path) -> Result<noodles::vcf::Header, ReadError> {
    let name = path
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| ReadError::UnsupportedFormat(path.display().to_string()))?
        .to_string();
    let file = File::open(path).map_err(ReadError::Io)?;
    if name.ends_with(".bcf") {
        let mut reader =
            noodles::bcf::io::Reader::from(bgzf::io::Reader::new(BufReader::new(file)));
        reader.read_header().map_err(ReadError::Io)
    } else if name.ends_with(".gz") || name.ends_with(".bgz") {
        let mut reader =
            noodles::vcf::io::Reader::new(bgzf::io::Reader::new(BufReader::new(file)));
        reader.read_header().map_err(ReadError::Io)
    } else if name.ends_with(".vcf") {
        let mut reader = noodles::vcf::io::Reader::new(BufReader::new(file));
        reader.read_header().map_err(ReadError::Io)
    } else {
        Err(ReadError::UnsupportedFormat(name))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Three samples: s0=0/1 s1=1/1 s2=./.
    fn site() -> RecordView {
        RecordView {
            contig: "Chr01".into(),
            pos: 100,
            qual: Some(50.0),
            is_pass: true,
            ref_allele: "A".into(),
            alt_alleles: vec!["G".into()],
            info_qd: None,
            info_fs: None,
            info_mq: None,
            info_sor: None,
            gt: vec![0, 1, 1, 1, -1, -1],
            dp: vec![Some(30), Some(30), Some(30)],
            gq: vec![Some(99), Some(99), Some(99)],
            n_samples: 3,
            masked: vec![false; 3],
        }
    }

    #[test]
    fn counts_full() {
        let s = site();
        let keep = vec![true; 3];
        let c = s.counts(&keep);
        assert_eq!((c.ac, c.an, c.ns), (3, 4, 2));
        assert_eq!((c.hom_ref, c.het, c.hom_alt), (0, 1, 1));
        assert_eq!(c.n_missing, 1);
        // MAF = min(3,1)/4 = 0.25
        assert!((c.maf() - 0.25).abs() < 1e-12);
    }

    #[test]
    fn counts_with_sample_drop() {
        let s = site();
        let keep = vec![true, true, false]; // drop the missing sample s2
        let c = s.counts(&keep);
        assert_eq!((c.ac, c.an, c.ns), (3, 4, 2));
        assert_eq!(c.n_missing, 0); // denominator and missing counts cover kept samples only
    }

    #[test]
    fn counts_half_call_counts_one_allele() {
        let mut s = site();
        s.gt = vec![0, -1, 1, 1, -1, -1]; // s0 is a half call with REF present
        let c = s.counts(&[true; 3]);
        assert_eq!((c.ac, c.an, c.ns), (2, 3, 2));
        assert_eq!(c.n_missing, 1);
    }

    #[test]
    fn mask_by_dp() {
        let mut s = site();
        s.dp = vec![Some(3), Some(30), Some(30)]; // s0 DP<5 -> masked
        let keep = vec![true; 3];
        s.mask(
            &keep,
            &Some(GtMask { dp_min: Some(5), dp_max: Some(100), gq_min: Some(20) }),
        );
        assert!(s.masked[0]);
        assert_eq!((s.gt[0], s.gt[1]), (-1, -1));
        let c = s.counts(&keep);
        assert_eq!((c.ac, c.an, c.ns), (2, 2, 1));
        assert_eq!(c.n_missing, 2);
    }

    /// Missing DP/GQ never masks, and a dropped sample is never masked even
    /// if its values fail the rule.
    #[test]
    fn mask_skips_missing_values_and_dropped_samples() {
        let mut s = site(); // s0=0/1 s1=1/1 s2=./.
        s.gt = vec![0, 1, 1, 1, 0, 0]; // give s2 a real call
        s.dp = vec![None, Some(3), Some(3)];
        let keep = vec![true, true, false];
        s.mask(
            &keep,
            &Some(GtMask { dp_min: Some(5), dp_max: Some(100), gq_min: Some(20) }),
        );
        assert_eq!(s.masked, vec![false, true, false]);
        assert_eq!((s.gt[0], s.gt[1]), (0, 1)); // s0: DP missing -> untouched
        assert_eq!((s.gt[2], s.gt[3]), (-1, -1)); // s1: DP below min -> masked
        assert_eq!((s.gt[4], s.gt[5]), (0, 0)); // s2 dropped: genotype preserved
    }

    #[test]
    fn mask_with_no_rule_is_a_noop() {
        let mut s = site();
        s.mask(&[true; 3], &None);
        assert!(s.masked.iter().all(|&m| !m));
        assert_eq!(s.gt, vec![0, 1, 1, 1, -1, -1]);
    }

    #[test]
    fn hwe_degenerate_sentinels() {
        let mut c = SiteCounts { ac: 0, an: 10, ..Default::default() };
        assert!(c.hwe_is_degenerate());
        c.ac = 10;
        assert!(c.hwe_is_degenerate());
        c.ac = 5;
        assert!(!c.hwe_is_degenerate());
        assert!((c.maf() - 0.5).abs() < 1e-12);
    }

    // ---- for_each_record ----

    const HEADER_TEXT: &str = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1
";

    fn test_records() -> Vec<RecordBuf> {
        let keys: noodles::vcf::variant::record_buf::samples::Keys =
            ["GT".to_string(), "DP".to_string()].into_iter().collect();
        (1..=3)
            .map(|pos| {
                RecordBuf::builder()
                    .set_reference_sequence_name("Chr01")
                    .set_variant_start(
                        noodles::core::Position::try_from(pos).expect("valid position"),
                    )
                    .set_reference_bases("A")
                    .set_alternate_bases(noodles::vcf::variant::record_buf::AlternateBases::from(
                        vec!["G".to_string()],
                    ))
                    .set_samples(noodles::vcf::variant::record_buf::Samples::new(
                        keys.clone(),
                        vec![
                            vec![
                                Some(noodles::vcf::variant::record_buf::samples::sample::Value::String("0/1".into())),
                                Some(noodles::vcf::variant::record_buf::samples::sample::Value::Integer(15)),
                            ],
                            vec![
                                Some(noodles::vcf::variant::record_buf::samples::sample::Value::String("./.".into())),
                                Some(noodles::vcf::variant::record_buf::samples::sample::Value::Integer(0)),
                            ],
                        ],
                    ))
                    .build()
            })
            .collect()
    }

    fn write_input(path: &Path) {
        let mut sink = crate::sink::open(Some(path), 1).expect("open sink");
        let header: noodles::vcf::Header = HEADER_TEXT.parse().expect("valid header");
        sink.write_header(&header).expect("write header");
        for rec in &test_records() {
            sink.write_record(&header, rec).expect("write record");
        }
        sink.finish().expect("finish");
    }

    #[test]
    fn for_each_record_reads_all_formats() {
        let dir = std::env::temp_dir().join(format!(
            "reed-io-record-fe-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("create scratch dir");
        for name in ["in.vcf", "in.vcf.gz", "in.bcf"] {
            let path = dir.join(name);
            write_input(&path);
            let mut n = 0;
            let mut poses = Vec::new();
            for_each_record(&path, |buf| {
                n += 1;
                poses.push(buf.variant_start().map(usize::from));
                assert_eq!(buf.reference_sequence_name(), "Chr01");
                Ok::<(), ReadError>(())
            })
            .unwrap_or_else(|e| panic!("{name}: {e}"));
            assert_eq!(n, 3, "{name}: record count");
            assert_eq!(
                poses,
                vec![Some(1), Some(2), Some(3)],
                "{name}: POS values"
            );
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn record_header_reads_all_formats() {
        let dir = std::env::temp_dir().join(format!(
            "reed-io-record-header-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("create scratch dir");
        for name in ["in.vcf", "in.vcf.gz", "in.bcf"] {
            let path = dir.join(name);
            write_input(&path);
            let header = record_header(&path).unwrap_or_else(|e| panic!("{name}: {e}"));
            let names: Vec<&str> = header.sample_names().iter().map(|s| s.as_str()).collect();
            assert_eq!(names, vec!["s0", "s1"], "{name}: sample names");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn for_each_record_rejects_unknown_extension() {
        let dir = std::env::temp_dir().join(format!(
            "reed-io-record-fe-bad-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("create scratch dir");
        let path = dir.join("in.txt");
        std::fs::write(&path, "not a vcf").expect("write scratch file");
        let err = for_each_record(&path, |_| Ok::<(), ReadError>(())).unwrap_err();
        assert!(matches!(err, ReadError::UnsupportedFormat(_)), "{err}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn for_each_record_propagates_callback_error() {
        let dir = std::env::temp_dir().join(format!(
            "reed-io-record-fe-cb-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("create scratch dir");
        let path = dir.join("in.vcf");
        write_input(&path);
        let err = for_each_record(&path, |_| {
            Err(ReadError::User("callback stopped".into()))
        })
        .unwrap_err();
        assert!(matches!(err, ReadError::User(_)), "{err}");
        let _ = std::fs::remove_dir_all(&dir);
    }
}
