//! hapmap → VCF (docs/rounds/v13/02-architecture.md §3.12, §7.7).
//!
//! Unlike the matrix targets, this is a **record→record** translation: a
//! hapmap row maps to one VCF record, so it never folds into an alignment.
//! Both formats are site-major; the worker renders each record into a
//! worker-local [`BlockWriter`] and the ordered sink streams the blocks out,
//! reusing the filter's `-O` machinery (D13-8: no new compression path).
//!
//! The header is fully known once the hapmap column row is scanned, so no
//! prescan or `##contig` emission is needed (TASSEL itself writes no
//! `##contig`).

use std::io::{BufRead, Read};
use std::path::{Path, PathBuf};
use std::sync::Arc;

use reed_core::ChunkId;
use reed_core::iupac;
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{ChunkJob, OrderedSink, run};
use reed_io::chunk::{CHUNK_TARGET_COMPRESSED, ChunkBytes, ChunkPayload, ChunkSource};
use reed_io::error::IoError;
use reed_io::write::block_writer::{BlockWriter, ChunkOutput, IndexKey};
use reed_io::write::ordered_output::OrderedOutput;

use crate::convert::{ConvertCounters, ConvertReport};

/// [`CHUNK_TARGET_COMPRESSED`] in `usize` (the same constant the plain text
/// source uses).
#[allow(clippy::cast_possible_truncation)]
const PLAIN_CHUNK_BYTES: usize = CHUNK_TARGET_COMPRESSED as usize;

/// The hapmap column row: sample names in column order plus the genotype
/// cell width the file uses (one-letter IUPAC or two-letter diploid).
#[derive(Clone, Debug)]
pub(crate) struct HapmapHeader {
    pub(crate) samples: Vec<String>,
    /// The first data row's sample region is two letters per taxon
    /// (TASSEL's `isOneLetter` inverse, BuilderFromHapMap.java:112-124).
    pub(crate) two_letter: bool,
    /// The first data row, read while sniffing the cell width. The source
    /// replays it so the sniff does not drop a record.
    pub(crate) first_row: Vec<u8>,
}

/// Reads the hapmap header (`##` comment lines, then the `rs#` column row)
/// and decides one- vs two-letter coding from the first data row.
///
/// # Errors
///
/// [`IoError`] on read failures or a stream that is not hapmap.
pub(crate) fn read_header(reader: &mut impl BufRead) -> Result<HapmapHeader, IoError> {
    let mut line = Vec::new();
    // Skip `##` metadata (TASSEL's `##SAMPLE=` annotations included).
    loop {
        line.clear();
        let n = reader.read_until(b'\n', &mut line)?;
        if n == 0 {
            return Err(IoError::Header(
                "hapmap input has no column header row (a line starting with `rs#`)".into(),
            ));
        }
        if !line.starts_with(b"##") {
            break;
        }
    }
    if !line.starts_with(b"rs#") {
        return Err(IoError::Header(
            "not a hapmap file: the first non-`##` line must be the `rs#` column row".into(),
        ));
    }
    let header_fields = split_ws(&line);
    if header_fields.len() <= 11 {
        return Err(IoError::Header(
            "the hapmap column row names no samples".into(),
        ));
    }
    let samples: Vec<String> = header_fields[11..]
        .iter()
        .map(|f| String::from_utf8_lossy(f).into_owned())
        .collect();
    // The first data row decides the cell width (TASSEL's heuristic, but on
    // the *token* width rather than its average: hapmap cells are one token
    // each, `A` or `AG`).
    line.clear();
    let n = reader.read_until(b'\n', &mut line)?;
    if n == 0 {
        return Err(IoError::Header("the hapmap file has no data rows".into()));
    }
    let first = split_ws(&line);
    if first.len() <= 11 {
        return Err(IoError::Header(
            "the hapmap data row names no samples".into(),
        ));
    }
    let cell = first[11];
    let two_letter = match cell.len() {
        1 => false,
        2 => true,
        other => {
            return Err(IoError::Header(format!(
                "hapmap genotype cell {cell:?} is {other} letters; expected 1 (IUPAC) or 2 (diploid)"
            )));
        }
    };
    Ok(HapmapHeader {
        samples,
        two_letter,
        // `line` still holds the sniffed row: hand it back so the source can
        // replay it (the sniff must not consume a record).
        first_row: line,
    })
}

/// Splits on runs of ASCII whitespace (the hapmap delimiter is a tab, but
/// TASSEL accepts spaces too).
fn split_ws(line: &[u8]) -> Vec<&[u8]> {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    let line = line.strip_suffix(b"\r").unwrap_or(line);
    line.split(|&b| b == b'\t' || b == b' ')
        .filter(|f| !f.is_empty())
        .collect()
}

/// A hapmap text source: plain 1 MiB chunks, `plain: true` (the engine
/// passes the bytes through untouched).
pub(crate) struct HapmapSource {
    reader: Box<dyn BufRead + Send>,
    /// The first data row (replayed before the reader's stream).
    pending: Option<Vec<u8>>,
    next_chunk_id: u64,
    buf: Vec<u8>,
}

impl HapmapSource {
    /// Wraps a reader whose header has already been consumed plus the first
    /// data row the header sniff read through.
    #[must_use]
    pub(crate) fn new(reader: impl Read + Send + 'static, first_row: Vec<u8>) -> Self {
        Self {
            reader: Box::new(std::io::BufReader::with_capacity(1 << 20, reader)),
            pending: (!first_row.is_empty()).then_some(first_row),
            next_chunk_id: 0,
            buf: Vec::with_capacity(PLAIN_CHUNK_BYTES + (1 << 14)),
        }
    }
}

impl ChunkSource for HapmapSource {
    fn next_chunk(&mut self) -> Result<Option<ChunkBytes>, IoError> {
        self.buf.clear();
        if let Some(pending) = self.pending.take() {
            self.buf.extend_from_slice(&pending);
        }
        while self.buf.len() < PLAIN_CHUNK_BYTES {
            let available = self.reader.fill_buf()?;
            if available.is_empty() {
                break;
            }
            let take = available.len().min(PLAIN_CHUNK_BYTES - self.buf.len());
            self.buf.extend_from_slice(&available[..take]);
            self.reader.consume(take);
        }
        if self.buf.is_empty() {
            return Ok(None);
        }
        let id = ChunkId(self.next_chunk_id);
        self.next_chunk_id += 1;
        Ok(Some(ChunkBytes {
            id,
            compressed: ChunkPayload::from(Arc::<[u8]>::from(self.buf.as_slice())),
            first_in_block_offset: 0,
            plain: true,
            compressed_start: 0,
            blocks: Vec::new(),
        }))
    }
}

/// The translation job: one hapmap row → one VCF record into the worker's
/// block writer.
pub(crate) struct HapmapVcfJob {
    pub(crate) n_samples: usize,
    pub(crate) two_letter: bool,
    pub(crate) level: u8,
    pub(crate) compressed: bool,
}

/// Worker-local scratch: the block writer, a record buffer and the per-row
/// cell byte ranges.
pub(crate) struct HapmapVcfScratch {
    writer: BlockWriter,
    record: Vec<u8>,
    warned: bool,
}

#[derive(Default)]
pub(crate) struct HapmapVcfPartial {
    pub(crate) output: ChunkOutput,
    pub(crate) n_records: u64,
    pub(crate) n_malformed: u64,
}

impl ChunkJob for HapmapVcfJob {
    type Partial = HapmapVcfPartial;
    type Scratch = Option<HapmapVcfScratch>;

    fn run_chunk(
        &self,
        _chunk: &ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| HapmapVcfScratch {
            writer: BlockWriter::new(self.level, self.compressed)
                .expect("invariant: the level comes from the CLI clamp"),
            record: Vec::with_capacity(1 << 10),
            warned: false,
        });
        let mut partial = HapmapVcfPartial::default();
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            self.translate(line, scratch, &mut partial)?;
        }
        partial.output = scratch.writer.finish_chunk().map_err(EngineError::from)?;
        Ok(partial)
    }

    fn run_seam(
        &self,
        _id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let mut partial = HapmapVcfPartial::default();
        let scratch = scratch.get_or_insert_with(|| HapmapVcfScratch {
            writer: BlockWriter::new(self.level, self.compressed)
                .expect("invariant: the level comes from the CLI clamp"),
            record: Vec::with_capacity(1 << 10),
            warned: false,
        });
        self.translate(line, scratch, &mut partial)?;
        partial.output = scratch.writer.finish_chunk().map_err(EngineError::from)?;
        Ok(partial)
    }
}

impl HapmapVcfJob {
    /// Translates one hapmap row into a VCF record (doc_v13 §3.12).
    fn translate(
        &self,
        line: &[u8],
        scratch: &mut HapmapVcfScratch,
        partial: &mut HapmapVcfPartial,
    ) -> Result<(), EngineError> {
        let fields = split_ws(line);
        let expected = 11 + self.n_samples;
        if fields.len() != expected {
            return self.malformed(scratch, partial, "column count differs from the header");
        }
        let pos1: u64 = match std::str::from_utf8(fields[3])
            .ok()
            .and_then(|s| s.parse().ok())
        {
            Some(pos) => pos,
            None => return self.malformed(scratch, partial, "POS is not an integer"),
        };
        // Alleles: `A/G`, `A/G/T`, `A` or `NA` (TASSEL's no-data row).
        let alleles_field = fields[1];
        let alleles: Vec<&[u8]> = if alleles_field == b"NA" {
            Vec::new()
        } else {
            split_ws(alleles_field) // `A/G` has no whitespace: one token
                .iter()
                .flat_map(|t| t.split(|&b| b == b'/'))
                .filter(|a| !a.is_empty())
                .collect()
        };
        let out = &mut scratch.record;
        out.clear();
        // CHROM.
        out.extend_from_slice(fields[2]);
        out.push(b'\t');
        // POS.
        out.extend_from_slice(pos1.to_string().as_bytes());
        out.push(b'\t');
        // ID (verbatim; `rs#` carries a real site name or `chrom:pos`).
        out.extend_from_slice(fields[0]);
        out.push(b'\t');
        // REF/ALT (declared order, B13-11).
        if alleles.is_empty() {
            out.extend_from_slice(b".\t.\t.\tPASS\t.\tGT");
        } else {
            out.extend_from_slice(alleles[0]);
            out.push(b'\t');
            if alleles.len() == 1 {
                out.push(b'.');
            } else {
                for (i, alt) in alleles[1..].iter().enumerate() {
                    if i > 0 {
                        out.push(b',');
                    }
                    out.extend_from_slice(alt);
                }
            }
            out.extend_from_slice(b"\t.\tPASS\t.\tGT");
        }
        // Genotypes, one per kept-in-header sample.
        for cell in &fields[11..] {
            out.push(b'\t');
            if alleles.is_empty() {
                out.extend_from_slice(b"./.");
                continue;
            }
            match self.cell_gt(cell, &alleles) {
                Ok(gt) => out.extend_from_slice(&gt),
                Err(()) => {
                    return self.malformed(
                        scratch,
                        partial,
                        "a genotype cell names an undeclared allele",
                    );
                }
            }
        }
        out.push(b'\n');
        let key = IndexKey {
            // The index is not written for `--vcf` (no `-W`), so the contig
            // id is inert; the position keeps the record anchored anyway.
            contig: reed_core::ContigId(0),
            start1: pos1.max(1),
            ref_len: u32::try_from(alleles.first().map_or(1, |a| a.len())).unwrap_or(1),
        };
        scratch
            .writer
            .push_record(&scratch.record, key)
            .map_err(EngineError::from)?;
        partial.n_records += 1;
        Ok(())
    }

    #[allow(clippy::unnecessary_wraps)] // callers `return self.malformed(..)` to skip the row
    fn malformed(
        &self,
        scratch: &mut HapmapVcfScratch,
        partial: &mut HapmapVcfPartial,
        why: &str,
    ) -> Result<(), EngineError> {
        partial.n_malformed += 1;
        if !scratch.warned {
            scratch.warned = true;
            tracing::warn!("reed convert --vcf: skipping malformed hapmap rows ({why})");
        }
        Ok(())
    }

    /// The GT string for one hapmap cell (doc_v13 §3.12).
    fn cell_gt(&self, cell: &[u8], alleles: &[&[u8]]) -> Result<Vec<u8>, ()> {
        let index_of = |letter: u8| -> Option<usize> {
            alleles
                .iter()
                .position(|a| a.len() == 1 && a[0].eq_ignore_ascii_case(&letter))
        };
        let mut set: Vec<usize> = Vec::with_capacity(2);
        if self.two_letter {
            if cell.len() != 2 {
                return Err(());
            }
            // TASSEL reads (second, first) (BuilderFromHapMap.java:258-266);
            // the output is unphased either way, so keep the written order.
            for &c in cell {
                if iupac::is_hapmap_missing(c) {
                    continue;
                }
                if let Some(i) = index_of(c) {
                    set.push(i);
                } else if let Some(exp) = iupac::ambig_expand(c) {
                    // A two-letter cell should not carry an ambiguity code,
                    // but tolerate one only if it resolves to a single
                    // allele.
                    if exp.len() == 1
                        && let Some(i) = index_of(exp[0])
                    {
                        set.push(i);
                        continue;
                    }
                    return Err(());
                } else {
                    return Err(());
                }
            }
        } else {
            if cell.len() != 1 {
                return Err(());
            }
            let c = cell[0];
            if iupac::is_hapmap_missing(c) {
                return Ok(b"./.".to_vec());
            }
            let Some(exp) = iupac::ambig_expand(c) else {
                return Err(());
            };
            if exp.len() > 2 {
                // V/H/D/B: a three-base set is not a diploid genotype.
                return Err(());
            }
            for &letter in exp {
                match index_of(letter) {
                    Some(i) => set.push(i),
                    None => return Err(()),
                }
            }
        }
        if set.is_empty() {
            return Ok(b"./.".to_vec());
        }
        set.sort_unstable();
        set.dedup();
        if set.len() == 1 {
            let i = set[0];
            return Ok(format!("{i}/{i}").into_bytes());
        }
        Ok(format!("{}/{}", set[0], set[1]).into_bytes())
    }
}

/// Ordered sink: writes each chunk's blocks in order through
/// [`OrderedOutput`] (the filter's `-O` writer), accumulating the record and
/// malformed counters.
struct HapmapVcfSink {
    out: OrderedOutput,
    n_records: u64,
    n_malformed: u64,
    /// `(records, malformed)` published for the driver (the engine returns
    /// only the sink's `u64`).
    counters: Arc<std::sync::Mutex<Option<(u64, u64)>>>,
}

impl OrderedSink<HapmapVcfPartial> for HapmapVcfSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut HapmapVcfPartial) -> Result<(), EngineError> {
        self.n_records += partial.n_records;
        self.n_malformed += partial.n_malformed;
        self.out
            .accept_chunk(&partial.output)
            .map_err(EngineError::from)
    }

    fn finish(self) -> Result<u64, EngineError> {
        self.out.finish().map_err(EngineError::from)?;
        *self
            .counters
            .lock()
            .expect("invariant: hapmap slot not poisoned") =
            Some((self.n_records, self.n_malformed));
        Ok(self.n_records)
    }

    fn deflate_time(&self) -> std::time::Duration {
        self.out.deflate_time()
    }
}

/// What the CLI hands the driver: where to write and how.
pub(crate) struct VcfOutput {
    /// Output file; `None` = stdout.
    pub(crate) path: Option<PathBuf>,
    /// `-Oz` (BGZF) vs `-Ov` (plain text).
    pub(crate) compressed: bool,
    /// Deflate level (1–9).
    pub(crate) level: u8,
    /// The `##reed_command` value; empty = omit the line.
    pub(crate) command_line: String,
    /// `--no-version`: drop the `##reed_*` provenance lines.
    pub(crate) no_version: bool,
}

/// Runs the hapmap→VCF conversion (doc_v13 §7.7).
///
/// # Errors
///
/// Request (usage-class) failures — a non-hapmap input, a bad output path —
/// and transport failures.
pub(crate) fn run_hapmap_to_vcf(
    input: Option<&Path>,
    output: &VcfOutput,
    budget: &ThreadBudget,
) -> Result<ConvertReport, crate::Error> {
    // Open the input; gzip (including BGZF) is inflated serially like
    // `PlainChunkSource::open_gzip`.
    let mut reader: Box<dyn BufRead + Send> = match input {
        None => {
            // A gzipped hapmap must also work when piped (doc_v15 D15-29):
            // peek the magic, then hand the snooped reader to the decoder —
            // the peeked bytes are still in its buffer.
            let mut probe = std::io::BufReader::with_capacity(1 << 20, std::io::stdin());
            let gzipped = {
                let head = probe
                    .fill_buf()
                    .map_err(|e| crate::Error::usage(format!("cannot read stdin: {e}")))?;
                head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b
            };
            if gzipped {
                let gz = flate2::read::MultiGzDecoder::new(probe);
                Box::new(std::io::BufReader::with_capacity(1 << 20, gz))
            } else {
                Box::new(probe)
            }
        }
        Some(p) => {
            let file = std::fs::File::open(p)
                .map_err(|e| crate::Error::usage(format!("cannot open {}: {e}", p.display())))?;
            let mut probe = std::io::BufReader::with_capacity(1 << 20, file);
            let gzipped = {
                let head = probe.fill_buf().map_err(|e| {
                    crate::Error::usage(format!("cannot read {}: {e}", p.display()))
                })?;
                head.len() >= 2 && head[0] == 0x1f && head[1] == 0x8b
            };
            if gzipped {
                let file = std::fs::File::open(p).map_err(|e| {
                    crate::Error::usage(format!("cannot open {}: {e}", p.display()))
                })?;
                let gz = flate2::read::MultiGzDecoder::new(file);
                Box::new(std::io::BufReader::with_capacity(1 << 20, gz))
            } else {
                Box::new(probe)
            }
        }
    };
    let header = read_header(&mut reader)?;
    if header.samples.is_empty() {
        return Err(crate::Error::usage("the hapmap header names no samples"));
    }
    let n_samples = header.samples.len();

    // The VCF header: everything but `##contig` is known now (§3.12).
    let mut header_bytes = String::from("##fileformat=VCFv4.3\n");
    if !output.no_version {
        use std::fmt::Write as _;
        let _ = writeln!(header_bytes, "##reed_version={}", env!("CARGO_PKG_VERSION"));
        if !output.command_line.is_empty() {
            let _ = writeln!(header_bytes, "##reed_command={}", output.command_line);
        }
    }
    header_bytes.push_str(
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT",
    );
    for name in &header.samples {
        header_bytes.push('\t');
        header_bytes.push_str(name);
    }
    header_bytes.push('\n');

    let writer: Box<dyn std::io::Write + Send> = match &output.path {
        None => Box::new(std::io::BufWriter::with_capacity(
            1 << 20,
            std::io::stdout(),
        )),
        Some(p) => Box::new(std::io::BufWriter::with_capacity(
            1 << 20,
            std::fs::File::create(p)
                .map_err(|e| crate::Error::usage(format!("cannot create {}: {e}", p.display())))?,
        )),
    };
    let mut out = OrderedOutput::new(
        writer,
        output.compressed,
        output.level,
        None,
        None,
        Arc::from(Vec::<Vec<u8>>::new()),
    );
    out.write_header(header_bytes.as_bytes())
        .map_err(EngineError::from)?;

    let job = HapmapVcfJob {
        n_samples,
        two_letter: header.two_letter,
        level: output.level,
        compressed: output.compressed,
    };
    let source = HapmapSource::new(reader, header.first_row);
    let counters = Arc::new(std::sync::Mutex::new(None));
    let sink = HapmapVcfSink {
        out,
        n_records: 0,
        n_malformed: 0,
        counters: Arc::clone(&counters),
    };
    let abort = reed_engine::abort::abort_for(budget);
    let mut source = source;
    let report = run(&mut source, &job, sink, budget, &abort)?;
    let (n_records, n_malformed) = counters
        .lock()
        .expect("invariant: hapmap slot not poisoned")
        .take()
        .unwrap_or((report.sink_report, 0));
    let mut files = Vec::new();
    if let Some(p) = &output.path {
        files.push(p.clone());
    }
    Ok(ConvertReport {
        files,
        counters: Some(ConvertCounters {
            total_records: n_records + n_malformed,
            malformed: n_malformed,
            ..ConvertCounters::default()
        }),
        sfs: None,
        n_records,
        elapsed: report.elapsed,
        metrics: report.metrics,
    })
}
