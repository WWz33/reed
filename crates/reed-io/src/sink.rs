//! Output sinks for VCF/BCF records.
//!
//! [`open`] dispatches on the output path: None -> plain text VCF on stdout
//! (Unix filter semantics; SIGPIPE is left to the caller), `*.bcf` -> BCF,
//! `*.gz` -> BGZF-compressed VCF, anything else -> plain text VCF file. The
//! compressed paths compress on a multithreaded BGZF pool sized by `threads`
//! (a zero value clamps to one worker); stdout and plain files ignore it.
//!
//! Writer facts these choices rest on (noodles 0.117):
//! - `bgzf::io::MultithreadedWriter::with_worker_count(NonZeroUsize, file)`
//!   compresses concurrently on a worker pool;
//! - `bcf::io::Writer::from(mtw)` and `vcf::io::Writer::new(mtw)` are plain
//!   `W: Write` wrappers; a BCF file is itself BGZF, so the BCF sink wraps
//!   `File` directly;
//! - records are written uniformly through
//!   `vcf::variant::io::Write::write_variant_record` (both writers implement
//!   the trait);
//! - finish is `writer.into_inner()` then `MultithreadedWriter::finish()`
//!   (flush + BGZF EOF block);
//! - stdout is `vcf::io::Writer::new(io::stdout().lock())` (plain text).

use std::fs::File;
use std::io::{self, BufWriter, Write as _};
use std::num::NonZeroUsize;
use std::path::Path;

use noodles::bcf;
use noodles::bgzf;
use noodles::vcf;
use noodles::vcf::variant::io::Write as _;
use noodles::vcf::variant::RecordBuf;

use crate::ReadError;

/// Output sink abstraction.
pub trait VariantSink {
    /// Writes the output header (must precede the first record).
    fn write_header(&mut self, header: &vcf::Header) -> io::Result<()>;
    /// Writes one record rendered against `header`.
    fn write_record(&mut self, header: &vcf::Header, record: &RecordBuf) -> io::Result<()>;
    /// Finish: flush / MT writer finish (BGZF EOF block).
    ///
    /// Idempotent: a repeated call is a no-op.
    fn finish(&mut self) -> io::Result<()>;
}

/// The concrete writer stack behind the [`VariantSink`] impl.
enum Sink {
    /// Plain text VCF on stdout (the default output).
    Stdout(vcf::io::Writer<io::StdoutLock<'static>>),
    /// Plain text VCF file.
    PlainFile(vcf::io::Writer<BufWriter<File>>),
    /// BGZF-compressed VCF.
    VcfGz(vcf::io::Writer<bgzf::io::MultithreadedWriter<File>>),
    /// BCF (the file itself is BGZF).
    Bcf(bcf::io::Writer<bgzf::io::MultithreadedWriter<File>>),
}

/// Sink handle returned by [`open`].
struct Writer {
    /// Taken by finish; `None` afterwards is what keeps finish idempotent.
    sink: Option<Sink>,
    /// Header actually used to render text-VCF records (see [`text_header`]).
    text_header: Option<vcf::Header>,
}

impl VariantSink for Writer {
    fn write_header(&mut self, header: &vcf::Header) -> io::Result<()> {
        let stored = match self.sink.as_mut() {
            Some(Sink::Stdout(w)) => {
                let h = text_header(header);
                w.write_header(&h)?;
                Some(h)
            }
            Some(Sink::PlainFile(w)) => {
                let h = text_header(header);
                w.write_header(&h)?;
                Some(h)
            }
            Some(Sink::VcfGz(w)) => {
                let h = text_header(header);
                w.write_header(&h)?;
                Some(h)
            }
            Some(Sink::Bcf(w)) => {
                w.write_header(header)?;
                None
            }
            None => return Err(closed()),
        };
        self.text_header = stored;
        Ok(())
    }

    fn write_record(&mut self, header: &vcf::Header, record: &RecordBuf) -> io::Result<()> {
        // Text paths render records against the stored (version-clamped)
        // header; BCF uses the caller's header, from which its string maps
        // were built.
        let header = match self.text_header.as_ref() {
            Some(h) => h,
            None => header,
        };
        match self.sink.as_mut() {
            Some(Sink::Stdout(w)) => w.write_variant_record(header, record),
            Some(Sink::PlainFile(w)) => w.write_variant_record(header, record),
            Some(Sink::VcfGz(w)) => w.write_variant_record(header, record),
            Some(Sink::Bcf(w)) => w.write_variant_record(header, record),
            None => Err(closed()),
        }
    }

    fn finish(&mut self) -> io::Result<()> {
        match self.sink.take() {
            None => Ok(()),
            Some(Sink::Stdout(w)) => w.into_inner().flush(),
            Some(Sink::PlainFile(w)) => w.into_inner().flush(),
            Some(Sink::VcfGz(w)) => {
                let mut mtw = w.into_inner();
                mtw.finish().map(|_| ())
            }
            Some(Sink::Bcf(w)) => {
                let mut mtw = w.into_inner();
                mtw.finish().map(|_| ())
            }
        }
    }
}

fn closed() -> io::Error {
    io::Error::new(io::ErrorKind::BrokenPipe, "sink already finished")
}

/// Open a sink by output path:
/// - None -> plain text VCF on stdout (the default; Unix pipe semantics);
/// - `*.bcf` -> BCF (MultithreadedWriter(threads));
/// - `*.vcf.gz` (or any `*.gz`) -> VCF.gz (MultithreadedWriter(threads));
/// - anything else -> plain text VCF file.
///
/// `threads` only sizes the MT compression pool; the stdout and plain paths
/// ignore it (a zero value clamps to one worker). An existing file is
/// truncated.
pub fn open(path: Option<&Path>, threads: usize) -> Result<Box<dyn VariantSink>, ReadError> {
    let sink = match path {
        None => Sink::Stdout(vcf::io::Writer::new(io::stdout().lock())),
        Some(p) => {
            let name = p.as_os_str().to_string_lossy();
            if name.ends_with(".bcf") {
                let file = File::create(p)?;
                Sink::Bcf(bcf::io::Writer::from(multithreaded(file, threads)))
            } else if name.ends_with(".gz") {
                let file = File::create(p)?;
                Sink::VcfGz(vcf::io::Writer::new(multithreaded(file, threads)))
            } else {
                let file = File::create(p)?;
                Sink::PlainFile(vcf::io::Writer::new(BufWriter::new(file)))
            }
        }
    };
    Ok(Box::new(Writer {
        sink: Some(sink),
        text_header: None,
    }))
}

/// MT BGZF compression pool over `file` (at least one worker).
fn multithreaded(file: File, threads: usize) -> bgzf::io::MultithreadedWriter<File> {
    bgzf::io::MultithreadedWriter::with_worker_count(
        NonZeroUsize::new(threads).unwrap_or(NonZeroUsize::MIN),
        file,
    )
}

/// Clamp the declared VCF version on text output to at most VCFv4.3.
///
/// noodles renders GT with a leading phasing separator ("/0/1") once the
/// header declares VCFv4.4 or newer; bcftools rejects that form ("Couldn't
/// read GT data: value not a number or '.'"). Text output therefore declares
/// at most VCFv4.3, whose GT rendering ("0/1") every reader accepts. BCF
/// output is binary with packed genotypes and is unaffected. The gold
/// pipeline data declares VCFv4.2, so this is a guard, not a hot path.
fn text_header(header: &vcf::Header) -> vcf::Header {
    const VCF_4_4: vcf::header::FileFormat = vcf::header::FileFormat::new(4, 4);
    let mut header = header.clone();
    if header.file_format() >= VCF_4_4 {
        *header.file_format_mut() = vcf::header::FileFormat::new(4, 3);
    }
    header
}

#[cfg(test)]
mod tests {
    use std::io::{BufRead as _, BufReader};
    use std::path::PathBuf;
    use std::process::Command;

    use noodles::vcf::variant::record::samples::keys::key;
    use noodles::vcf::variant::record::samples::series::value::genotype::Phasing;
    use noodles::vcf::variant::record_buf::info::field::Value as InfoValue;
    use noodles::vcf::variant::record_buf::samples::sample::value::genotype::Allele;
    use noodles::vcf::variant::record_buf::samples::sample::value::Genotype;
    use noodles::vcf::variant::record_buf::samples::sample::Value;
    use noodles::vcf::variant::record_buf::{AlternateBases, Filters, Info, Samples};
    use noodles::core::Position;

    use super::*;

    const HEADER_TEXT: &str = "\
##fileformat=VCFv4.4
##contig=<ID=chr1,length=1000000>
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2
";

    fn header() -> vcf::Header {
        HEADER_TEXT.parse().expect("valid header text")
    }

    fn unphased(a: usize, b: usize) -> Value {
        Value::Genotype(Genotype::from_iter([
            Allele::new(Some(a), Phasing::Unphased),
            Allele::new(Some(b), Phasing::Unphased),
        ]))
    }

    fn samples() -> Samples {
        Samples::new(
            [String::from(key::GENOTYPE)].into_iter().collect(),
            vec![
                vec![Some(unphased(0, 1))],
                vec![Some(unphased(0, 0))],
            ],
        )
    }

    fn test_records() -> Vec<RecordBuf> {
        (1..=3)
            .map(|pos| {
                RecordBuf::builder()
                    .set_reference_sequence_name("chr1")
                    .set_variant_start(Position::try_from(pos).unwrap())
                    .set_reference_bases("A")
                    .set_alternate_bases(AlternateBases::from_iter([String::from("G")]))
                    .set_quality_score(30.0)
                    .set_info(Info::from_iter([(
                        String::from("DP"),
                        Some(InfoValue::Integer(57)),
                    )]))
                    .set_samples(samples())
                    .build()
            })
            .collect()
    }

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "reed-io-sink-{}-{name}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("create scratch dir");
        dir
    }

    fn write_with_sink(path: Option<&Path>) -> Result<(), ReadError> {
        let mut snk = open(path, 2)?;
        let header = header();
        snk.write_header(&header)?;
        for rec in &test_records() {
            snk.write_record(&header, rec)?;
        }
        Ok(snk.finish()?)
    }

    fn read_back_plain_vcf(path: &Path) -> io::Result<Vec<RecordBuf>> {
        let file = File::open(path)?;
        let mut reader = vcf::io::Reader::new(BufReader::new(file));
        let header = reader.read_header()?;
        reader.record_bufs(&header).collect()
    }

    fn read_back_vcf_gz(path: &Path) -> io::Result<Vec<RecordBuf>> {
        let file = File::open(path)?;
        let mut reader = vcf::io::Reader::new(bgzf::io::Reader::new(BufReader::new(file)));
        let header = reader.read_header()?;
        reader.record_bufs(&header).collect()
    }

    fn read_back_bcf(path: &Path) -> io::Result<Vec<RecordBuf>> {
        let file = File::open(path)?;
        // Reader::from over an explicit bgzf::io::Reader: one BGZF layer.
        // (Handing a bgzf::io::Reader to Reader::new would stack a second.)
        let mut reader = bcf::io::Reader::from(bgzf::io::Reader::new(BufReader::new(file)));
        let header = reader.read_header()?;
        reader.record_bufs(&header).collect()
    }

    /// BCF stores an empty FILTER column as filter id 0, which reads back as
    /// PASS; both mean "no filter applied".
    fn filters_equiv(a: &Filters, b: &Filters) -> bool {
        a == b
            || (a.as_ref().is_empty() && b.is_pass())
            || (b.as_ref().is_empty() && a.is_pass())
    }

    fn assert_same_records(expected: &[RecordBuf], got: &[RecordBuf]) {
        assert_eq!(expected.len(), got.len(), "record count");
        for (e, g) in expected.iter().zip(got) {
            assert_eq!(e.reference_sequence_name(), g.reference_sequence_name());
            assert_eq!(e.variant_start(), g.variant_start(), "POS");
            assert_eq!(e.reference_bases(), g.reference_bases(), "REF");
            assert_eq!(e.alternate_bases(), g.alternate_bases(), "ALT");
            assert_eq!(e.quality_score(), g.quality_score(), "QUAL");
            assert!(
                filters_equiv(e.filters(), g.filters()),
                "FILTER mismatch: {:?} vs {:?}",
                e.filters(),
                g.filters()
            );
            assert_eq!(e.info(), g.info(), "INFO");
            assert_eq!(e.samples(), g.samples(), "FORMAT/samples");
        }
    }

    #[test]
    fn roundtrip_plain_vcf_file() {
        let dir = scratch("plain");
        let path = dir.join("out.vcf");
        write_with_sink(Some(&path)).expect("write plain vcf");
        let got = read_back_plain_vcf(&path).expect("read back plain vcf");
        assert_same_records(&test_records(), &got);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn roundtrip_vcf_gz() {
        let dir = scratch("vcf-gz");
        let path = dir.join("out.vcf.gz");
        write_with_sink(Some(&path)).expect("write vcf.gz");
        let got = read_back_vcf_gz(&path).expect("read back vcf.gz");
        assert_same_records(&test_records(), &got);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn roundtrip_bcf() {
        let dir = scratch("bcf");
        let path = dir.join("out.bcf");
        write_with_sink(Some(&path)).expect("write bcf");
        let got = read_back_bcf(&path).expect("read back bcf");
        assert_same_records(&test_records(), &got);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn stdout_sink_smoke() {
        // Content cannot be asserted under the test harness; the point is that
        // the default (None) path opens, writes and finishes without error.
        write_with_sink(None).expect("write to stdout");
    }

    #[test]
    fn threads_zero_clamps_to_one_worker() {
        let dir = scratch("zero-threads");
        let path = dir.join("out.vcf.gz");
        let mut snk = open(Some(&path), 0).expect("open with threads=0");
        let header = header();
        snk.write_header(&header).expect("write header");
        snk.write_record(&header, &test_records()[0])
            .expect("write record");
        snk.finish().expect("finish");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn finish_is_idempotent_across_all_sinks() {
        let dir = scratch("idempotent");
        for name in ["out.vcf", "out.vcf.gz", "out.bcf"] {
            let path = dir.join(name);
            let mut snk = open(Some(&path), 1).expect("open");
            let header = header();
            snk.write_header(&header).expect("write header");
            snk.write_record(&header, &test_records()[0])
                .expect("write record");
            snk.finish().expect("first finish");
            snk.finish().expect("second finish is a no-op");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn write_after_finish_reports_error() {
        let dir = scratch("after-finish");
        let path = dir.join("out.vcf");
        let mut snk = open(Some(&path), 1).expect("open");
        let header = header();
        snk.write_header(&header).expect("write header");
        snk.finish().expect("finish");
        assert!(snk.write_header(&header).is_err());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn text_output_clamps_vcf_4_4_to_4_3() {
        let dir = scratch("clamp");
        let path = dir.join("out.vcf");
        write_with_sink(Some(&path)).expect("write plain vcf");
        let first = std::fs::read_to_string(&path)
            .expect("read output")
            .lines()
            .next()
            .expect("non-empty output")
            .to_string();
        // HEADER_TEXT declares VCFv4.4; noodles would render GT as "/0/1",
        // which bcftools rejects. The sink must downgrade the declaration.
        assert_eq!(first, "##fileformat=VCFv4.3");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn bcftools_reads_both_on_disk_formats() {
        // Optional sanity check against the system toolchain; skipped when
        // bcftools is not installed.
        if Command::new("bcftools").arg("--version").output().is_err() {
            eprintln!("bcftools not found; skipping");
            return;
        }
        let dir = scratch("bcftools");
        for name in ["out.vcf.gz", "out.bcf"] {
            let path = dir.join(name);
            write_with_sink(Some(&path)).expect("write output");
            let out = Command::new("bcftools")
                .arg("view")
                .arg("-H")
                .arg(&path)
                .output()
                .expect("run bcftools view");
            assert!(
                out.status.success(),
                "bcftools view failed on {name}: {}",
                String::from_utf8_lossy(&out.stderr)
            );
            let lines = out.stdout.as_slice().lines().count();
            assert_eq!(lines, 3, "bcftools view -H line count for {name}");
        }
        let _ = std::fs::remove_dir_all(&dir);
    }
}
