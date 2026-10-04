//! VCF/BCF window reading: L1 (noodles) -> L2 ([`WindowBuffer`]).
//!
//! Rules encoded here (see `reed/DESIGN.md`):
//!
//! - v0.1 is diploid and biallelic: haploid, polyploid and multiallelic records
//!   are explicit errors, never silently altered.
//! - `contig_len` is required: a contig without a length in the header (Beagle
//!   strips `##contig`) is an error, not a silent fallback.
//! - Windows never cross contigs; the grid is `[start, start + window_size)`
//!   with `start` advancing by `step`; a trailing window is clamped to the
//!   contig length. Empty grid windows between data are emitted with zero sites.
//! - Genotypes are decoded once per site (column-major `series.iter()` — never
//!   per-sample random access, which is quadratic in noodles).

use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

use noodles::bcf as bcf;
use noodles::bgzf as bgzf;
use noodles::vcf as vcf;
use vcf::variant::record::samples::series::value::genotype::Phasing;
use vcf::variant::record::samples::series::Value;

use reed_core::{Bits, ContigId, SiteCompleteness, WindowBuffer, WindowError};

/// Failures while reading windows.
#[derive(Debug)]
pub enum ReadError {
    Io(io::Error),
    Core(WindowError),
    /// File name does not end in `.vcf`, `.vcf.gz`, `.vcf.bgz` or `.bcf`.
    UnsupportedFormat(String),
    /// A contig used by a record has no length in the VCF header.
    MissingContigLength { contig: String },
    /// A record has no POS.
    MissingPosition { contig: String },
    /// A genotype is not diploid (v0.1).
    NotDiploid { contig: String, pos: u64, observed: usize },
    /// An allele index above 1 (v0.1: biallelic).
    Multiallelic { contig: String, pos: u64 },
    /// The record has no GT series.
    MissingGt { contig: String, pos: u64 },
    /// Records are not sorted by position within a contig.
    OutOfOrder { contig: String, prev: u32, pos: u32 },
    /// A record's POS (1-based) exceeds the contig length declared in the
    /// header — the record would fall outside every window.
    PositionBeyondContigLength { contig: String, pos: u64, len: u64 },
    /// A contig's records are not contiguous (the contig re-appears after
    /// another contig); the window grid for it would be corrupted.
    ContigRevisit { contig: String },
    /// Invalid options or a user callback error.
    User(String),
}

impl fmt::Display for ReadError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            ReadError::Io(e) => write!(f, "I/O error: {e}"),
            ReadError::Core(e) => write!(f, "window invariant violated: {e}"),
            ReadError::UnsupportedFormat(name) => write!(
                f,
                "unsupported input {name:?} (v0.1: .vcf, .vcf.gz, .vcf.bgz, .bcf)"
            ),
            ReadError::MissingContigLength { contig } => write!(
                f,
                "contig {contig:?} has no length in the VCF header; window semantics require \
                 contig lengths (note: Beagle output strips ##contig). Add ##contig lines or \
                 supply lengths before reading"
            ),
            ReadError::MissingPosition { contig } => {
                write!(f, "record on contig {contig:?} has no POS")
            }
            ReadError::NotDiploid { contig, pos, observed } => write!(
                f,
                "non-diploid genotype ({observed} alleles) at {contig}:{pos}; v0.1 supports \
                 diploid calls only"
            ),
            ReadError::Multiallelic { contig, pos } => write!(
                f,
                "multiallelic genotype at {contig}:{pos}; v0.1 supports biallelic sites only"
            ),
            ReadError::MissingGt { contig, pos } => {
                write!(f, "record at {contig}:{pos} has no GT format series")
            }
            ReadError::OutOfOrder { contig, prev, pos } => write!(
                f,
                "records on {contig:?} are not position-sorted ({prev} followed by {pos})"
            ),
            ReadError::PositionBeyondContigLength { contig, pos, len } => write!(
                f,
                "record at {contig}:{pos} is beyond the declared contig length {len} \
                 (header/data contradiction)"
            ),
            ReadError::ContigRevisit { contig } => write!(
                f,
                "records for {contig:?} appear again after another contig; VCF requires \
                 each contig's records to be contiguous"
            ),
            ReadError::User(msg) => write!(f, "{msg}"),
        }
    }
}

impl std::error::Error for ReadError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ReadError::Io(e) => Some(e),
            ReadError::Core(e) => Some(e),
            _ => None,
        }
    }
}

impl From<io::Error> for ReadError {
    fn from(e: io::Error) -> Self {
        ReadError::Io(e)
    }
}

impl From<WindowError> for ReadError {
    fn from(e: WindowError) -> Self {
        ReadError::Core(e)
    }
}

/// How to cut the stream into windows and what the input is declared to contain.
#[derive(Clone, Copy, Debug)]
pub struct ReadOptions {
    pub window_size: u32,
    pub step: u32,
    pub manifest: SiteCompleteness,
}

impl Default for ReadOptions {
    fn default() -> Self {
        Self {
            window_size: 100_000,
            step: 100_000,
            manifest: SiteCompleteness::VariantOnly(reed_core::Ascertainment::Unknown),
        }
    }
}

/// Callback receiving each materialised window.
pub type WindowSink<'a> = &'a mut dyn FnMut(WindowBuffer) -> Result<(), ReadError>;

/// Stream summary.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct StreamSummary {
    pub windows: u64,
    pub sites: u64,
}

/// Streams every window to `f` in contig/position order.
pub fn for_each_window<F>(path: &Path, opts: &ReadOptions, f: F) -> Result<StreamSummary, ReadError>
where
    F: FnMut(WindowBuffer) -> Result<(), ReadError>,
{
    if opts.window_size == 0 || opts.step == 0 || opts.step > opts.window_size {
        return Err(ReadError::User(format!(
            "invalid window options: window_size={}, step={} (require 0 < step <= window_size)",
            opts.window_size, opts.step
        )));
    }
    let name = path
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| ReadError::UnsupportedFormat(path.display().to_string()))?
        .to_string();
    let file = File::open(path)?;
    if name.ends_with(".bcf") {
        stream_bcf(
            bcf::io::Reader::new(BufReader::new(file)),
            opts,
            f,
        )
    } else if name.ends_with(".vcf.gz") || name.ends_with(".vcf.bgz") {
        stream_vcf(
            vcf::io::Reader::new(bgzf::io::Reader::new(BufReader::new(file))),
            opts,
            f,
        )
    } else if name.ends_with(".vcf") {
        stream_vcf(vcf::io::Reader::new(BufReader::new(file)), opts, f)
    } else {
        Err(ReadError::UnsupportedFormat(name))
    }
}

/// Convenience: materialise all windows (small inputs and tests).
pub fn read_windows(path: &Path, opts: &ReadOptions) -> Result<Vec<WindowBuffer>, ReadError> {
    let mut out = Vec::new();
    for_each_window(path, opts, |w| {
        out.push(w);
        Ok(())
    })?;
    Ok(out)
}

/// Header facts useful for declaring [`ReadOptions::manifest`] and building a
/// [`reed_core::PopulationView`].
#[derive(Clone, Debug)]
pub struct HeaderInfo {
    pub n_samples: usize,
    /// Contigs in header order with their declared lengths (`None` = missing;
    /// the reader will refuse records on such contigs).
    pub contigs: Vec<(String, Option<u32>)>,
}

/// Reads only the header of a VCF/BCF file.
pub fn header_info(path: &Path) -> Result<HeaderInfo, ReadError> {
    let name = path
        .file_name()
        .and_then(|s| s.to_str())
        .ok_or_else(|| ReadError::UnsupportedFormat(path.display().to_string()))?
        .to_string();
    let file = File::open(path)?;
    let header = if name.ends_with(".bcf") {
        bcf::io::Reader::new(BufReader::new(file))
            .read_header()
            .map_err(ReadError::Io)?
    } else if name.ends_with(".vcf.gz") || name.ends_with(".vcf.bgz") {
        vcf::io::Reader::new(bgzf::io::Reader::new(BufReader::new(file)))
            .read_header()
            .map_err(ReadError::Io)?
    } else if name.ends_with(".vcf") {
        vcf::io::Reader::new(BufReader::new(file))
            .read_header()
            .map_err(ReadError::Io)?
    } else {
        return Err(ReadError::UnsupportedFormat(name));
    };
    let contigs = header
        .contigs()
        .iter()
        .map(|(name, contig)| (name.clone(), contig.length().map(|l| l as u32)))
        .collect();
    Ok(HeaderInfo {
        n_samples: header.sample_names().len(),
        contigs,
    })
}

// ---------------------------------------------------------------- assembler

struct PendingSite {
    pos: u32,
    alleles: Vec<i8>,
    missing: Vec<bool>,
    phased: Vec<bool>,
    is_variant: bool,
}

struct ContigState {
    name: String,
    id: ContigId,
    len: u32,
    next_start: u64,
    last_pos: u32,
    pending: Vec<PendingSite>,
}

struct Assembler {
    window_size: u64,
    step: u64,
    manifest: SiteCompleteness,
    n_samples: u32,
    contig_lens: HashMap<String, u32>,
    contig_ids: HashMap<String, ContigId>,
    finished_contigs: std::collections::HashSet<String>,
    next_contig_id: u32,
    state: Option<ContigState>,
    sites: u64,
}

impl Assembler {
    fn new(opts: &ReadOptions, header: &vcf::Header) -> Self {
        let mut contig_lens = HashMap::new();
        for (name, contig) in header.contigs() {
            if let Some(len) = contig.length() {
                contig_lens.insert(name.clone(), len as u32);
            }
        }
        Self {
            window_size: opts.window_size as u64,
            step: opts.step as u64,
            manifest: opts.manifest,
            n_samples: header.sample_names().len() as u32,
            contig_lens,
            contig_ids: HashMap::new(),
            finished_contigs: std::collections::HashSet::new(),
            next_contig_id: 0,
            state: None,
            sites: 0,
        }
    }

    fn push_site(
        &mut self,
        contig: &str,
        pos: u32,
        alleles: Vec<i8>,
        missing: Vec<bool>,
        phased: Vec<bool>,
        on_window: WindowSink<'_>,
    ) -> Result<(), ReadError> {
        if let Some(state) = &self.state {
            if state.name != contig {
                self.finish_contig(on_window)?;
            }
        }
        if self.state.is_none() {
            if self.finished_contigs.contains(contig) {
                // VCF requires a contig's records to be contiguous; a revisit
                // would corrupt the window grid for that contig.
                return Err(ReadError::ContigRevisit {
                    contig: contig.to_string(),
                });
            }
            let len = match self.contig_lens.get(contig) {
                Some(&len) => len,
                None => {
                    return Err(ReadError::MissingContigLength {
                        contig: contig.to_string(),
                    })
                }
            };
            let id = if let Some(&id) = self.contig_ids.get(contig) {
                id
            } else {
                let id = ContigId(self.next_contig_id);
                self.next_contig_id += 1;
                self.contig_ids.insert(contig.to_string(), id);
                id
            };
            let first_start = (pos as u64 / self.step) * self.step;
            self.state = Some(ContigState {
                name: contig.to_string(),
                id,
                len,
                next_start: first_start,
                last_pos: pos,
                pending: Vec::new(),
            });
        }

        let state = self.state.as_mut().expect("state initialised above");
        if pos < state.last_pos {
            return Err(ReadError::OutOfOrder {
                contig: contig.to_string(),
                prev: state.last_pos,
                pos,
            });
        }
        if pos as u64 >= state.len as u64 {
            return Err(ReadError::PositionBeyondContigLength {
                contig: contig.to_string(),
                pos: pos as u64 + 1,
                len: state.len as u64,
            });
        }
        let is_variant = alleles.iter().any(|&a| a > 0);
        while pos as u64 >= state.next_start + self.window_size {
            if let Some(window) = build_next_window(
                state,
                self.window_size,
                self.step,
                self.manifest,
                self.n_samples,
            )? {
                on_window(window)?;
            }
        }
        state.pending.push(PendingSite {
            pos,
            alleles,
            missing,
            phased,
            is_variant,
        });
        state.last_pos = pos;
        self.sites += 1;
        Ok(())
    }

    fn finish_contig(&mut self, on_window: WindowSink<'_>) -> Result<(), ReadError> {
        if let Some(state) = self.state.as_mut() {
            while state.next_start <= state.last_pos as u64 {
                if let Some(window) = build_next_window(
                    state,
                    self.window_size,
                    self.step,
                    self.manifest,
                    self.n_samples,
                )? {
                    on_window(window)?;
                }
            }
            self.finished_contigs.insert(std::mem::take(&mut state.name));
        }
        self.state = None;
        Ok(())
    }
}

fn build_next_window(
    state: &mut ContigState,
    window_size: u64,
    step: u64,
    manifest: SiteCompleteness,
    n_samples: u32,
) -> Result<Option<WindowBuffer>, ReadError> {
    let start = state.next_start;
    if start >= state.len as u64 {
        state.next_start += step;
        return Ok(None);
    }
    let end = (start + window_size).min(state.len as u64);

    let sites: Vec<&PendingSite> = state
        .pending
        .iter()
        .filter(|s| (s.pos as u64) >= start && (s.pos as u64) < end)
        .collect();

    let n = sites.len();
    let mut pos = Vec::with_capacity(n);
    let mut is_variant = Vec::with_capacity(n);
    let mut genotypes = Vec::with_capacity(n * n_samples as usize * 2);
    let mut missing = Vec::with_capacity(n * n_samples as usize);
    let mut phased = Vec::with_capacity(n * n_samples as usize);
    for site in &sites {
        pos.push(site.pos);
        is_variant.push(site.is_variant);
        genotypes.extend_from_slice(&site.alleles);
        missing.extend_from_slice(&site.missing);
        phased.extend_from_slice(&site.phased);
    }

    let window = WindowBuffer::new(
        manifest,
        state.id,
        state.len,
        start,
        end,
        n_samples,
        pos,
        Bits::from_bools(&is_variant),
        genotypes,
        Bits::from_bools(&missing),
        Bits::from_bools(&phased),
    )?;

    state.next_start += step;
    state.pending.retain(|s| (s.pos as u64) >= state.next_start);
    Ok(Some(window))
}

// ------------------------------------------------------------------ backends

fn record_key<R: vcf::variant::Record + ?Sized>(
    record: &R,
    header: &vcf::Header,
) -> Result<(String, u32), ReadError> {
    let contig = record
        .reference_sequence_name(header)
        .map_err(ReadError::Io)?
        .to_string();
    let start = record
        .variant_start()
        .transpose()
        .map_err(ReadError::Io)?
        .ok_or_else(|| ReadError::MissingPosition {
            contig: contig.clone(),
        })?;
    let pos0 = usize::from(start).saturating_sub(1) as u32;
    Ok((contig, pos0))
}

/// Decoded GT column: allele indices (2 per sample), missing bits, phasing bits.
type DecodedGt = (Vec<i8>, Vec<bool>, Vec<bool>);

fn decode_gt(
    series: &dyn vcf::variant::record::samples::Series,
    header: &vcf::Header,
    n_samples: usize,
    contig: &str,
    pos: u32,
) -> Result<DecodedGt, ReadError> {
    let at = || ReadError::MissingGt {
        contig: contig.to_string(),
        pos: pos as u64 + 1,
    };
    let mut alleles = Vec::with_capacity(n_samples * 2);
    let mut missing = Vec::with_capacity(n_samples);
    let mut phased = Vec::with_capacity(n_samples);

    for item in series.iter(header) {
        let value = item.map_err(ReadError::Io)?;
        // `None` means "no genotype for this sample" (a `.` sample column, or a
        // missing GT subfield) — valid per VCF spec; decode as ./. .
        let Some(Value::Genotype(genotype)) = value else {
            alleles.push(-1);
            alleles.push(-1);
            missing.push(true);
            phased.push(false);
            continue;
        };
        let mut n_alleles = 0usize;
        let mut phased_any = false;
        let mut is_missing = false;
        for pair in genotype.iter() {
            let (allele, phasing) = pair.map_err(ReadError::Io)?;
            n_alleles += 1;
            if phasing == Phasing::Phased {
                phased_any = true;
            }
            match allele {
                Some(index) => {
                    if index > 1 {
                        return Err(ReadError::Multiallelic {
                            contig: contig.to_string(),
                            pos: pos as u64 + 1,
                        });
                    }
                    alleles.push(index as i8);
                }
                None => {
                    alleles.push(-1);
                    is_missing = true;
                }
            }
        }
        if n_alleles == 0 {
            // BCF end-of-vector-only GT: spec-legal "no genotype".
            alleles.push(-1);
            alleles.push(-1);
            missing.push(true);
            phased.push(false);
            continue;
        }
        if n_alleles != 2 {
            return Err(ReadError::NotDiploid {
                contig: contig.to_string(),
                pos: pos as u64 + 1,
                observed: n_alleles,
            });
        }
        missing.push(is_missing);
        phased.push(phased_any);
    }

    if phased.len() != n_samples {
        return Err(at());
    }
    Ok((alleles, missing, phased))
}

fn stream_bcf<R, F>(
    mut reader: bcf::io::Reader<R>,
    opts: &ReadOptions,
    mut f: F,
) -> Result<StreamSummary, ReadError>
where
    R: BufRead,
    F: FnMut(WindowBuffer) -> Result<(), ReadError>,
{
    let header = reader.read_header().map_err(ReadError::Io)?;
    let mut assembler = Assembler::new(opts, &header);
    let n_samples = assembler.n_samples as usize;
    let mut windows = 0u64;
    let mut emit = |window: WindowBuffer| -> Result<(), ReadError> {
        windows += 1;
        f(window)
    };

    for result in reader.records() {
        let record = result.map_err(ReadError::Io)?;
        let (contig, pos) = record_key(&record, &header)?;
        let samples = record.samples().map_err(ReadError::Io)?;
        let mut gt_series = None;
        for series in samples.series() {
            let series = series.map_err(ReadError::Io)?;
            if series.name(&header).map_err(ReadError::Io)? == "GT" {
                gt_series = Some(series);
                break;
            }
        }
        let series = gt_series.ok_or_else(|| ReadError::MissingGt {
            contig: contig.clone(),
            pos: pos as u64 + 1,
        })?;
        let (alleles, missing, phased) = decode_gt(&series, &header, n_samples, &contig, pos)?;
        assembler.push_site(&contig, pos, alleles, missing, phased, &mut emit)?;
    }
    assembler.finish_contig(&mut emit)?;

    Ok(StreamSummary {
        windows,
        sites: assembler.sites,
    })
}

fn stream_vcf<R, F>(
    mut reader: vcf::io::Reader<R>,
    opts: &ReadOptions,
    mut f: F,
) -> Result<StreamSummary, ReadError>
where
    R: BufRead,
    F: FnMut(WindowBuffer) -> Result<(), ReadError>,
{
    let header = reader.read_header().map_err(ReadError::Io)?;
    let mut assembler = Assembler::new(opts, &header);
    let n_samples = assembler.n_samples as usize;
    let mut windows = 0u64;
    let mut emit = |window: WindowBuffer| -> Result<(), ReadError> {
        windows += 1;
        f(window)
    };

    for result in reader.records() {
        let record = result.map_err(ReadError::Io)?;
        let (contig, pos) = record_key(&record, &header)?;
        let samples = record.samples();
        let series = samples.select("GT").ok_or_else(|| ReadError::MissingGt {
            contig: contig.clone(),
            pos: pos as u64 + 1,
        })?;
        let (alleles, missing, phased) = decode_gt(&series, &header, n_samples, &contig, pos)?;
        assembler.push_site(&contig, pos, alleles, missing, phased, &mut emit)?;
    }
    assembler.finish_contig(&mut emit)?;

    Ok(StreamSummary {
        windows,
        sites: assembler.sites,
    })
}
