//! `run_filter`: the L4 filter driver (§11.4). Only this layer sees both
//! the engine and the strategies, so the strategy choice lives here. The
//! compiled [`FilterPlan`] executes as `SinglePass` (one chunk pipeline:
//! head → mask → tail → edit → worker-local deflate) or `TwoPassReparse`
//! (pass0 = GT-only per-sample missing reduction; pass1 = SinglePass over
//! the surviving columns). The ordered sink appends worker-compressed
//! BGZF blocks in chunk order and feeds the tabix indexer while writing
//! (§6.6) — no second read for `--write-index`.

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use reed_core::estimators::site_counts::SiteCounts;
use reed_core::sample::SampleMask;
use reed_core::Bitset;
use reed_core::ChunkId;
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{run, ChunkJob, OrderedSink};
use reed_filter::eval::edit_plan::fill_info_set;
use reed_filter::eval::head::{eval_head, HeadScratch};
use reed_filter::eval::mask::eval_mask;
use reed_filter::eval::tail::eval_tail;
use reed_filter::funnel::{FilterReport, FunnelPartial};
use reed_filter::plan::{ExecStrategy, FilterPlan, HeadRule};
use reed_filter::spec::compat::from_gold_conf;
use reed_filter::spec::{CohortRule as CohortRuleV2, FilterSpecV2};
use reed_filter::FilterError;
use reed_io::bgzf::block::EOF_BLOCK;
use reed_io::bgzf::deflate::Deflater;
use reed_io::chunk::{BgzfChunkSource, ChunkSource};
use noodles;
use reed_io::noodles_vcf as vcf;
use reed_io::vcf::fixed::parse_pos1;
use reed_io::Locus;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::line_edit::{edit_vcf_line, EditPlan, InfoValueOut};
use reed_io::vcf::sample_scan::SampleScanner;

/// BGZF maximum uncompressed payload per block (htslib `BGZF_BLOCK_SIZE`).
const BGZF_BLOCK_SIZE: usize = 65_280;

/// Run options.
#[derive(Clone, Debug, Default)]
pub struct FilterOptions {
    /// Worker threads (clamped to ≥ 1).
    pub threads: usize,
    /// Output path; `None` = stdout (plain VCF).
    pub output: Option<PathBuf>,
    /// Write a tabix index next to a `.vcf.gz` output.
    pub write_index: bool,
    /// BGZF compression level (1–9; htslib default 6).
    pub compress_level: u8,
    /// Test-only strategy override (§11.4: the chosen strategy lands in
    /// the report either way).
    pub strategy: Option<ExecStrategy>,
    /// The `##reed_command` line (captured argv).
    pub command_line: String,
    /// Suppress the human stderr table.
    pub quiet: bool,
    /// Write the JSON report here.
    pub report_json: Option<PathBuf>,
}

fn io_err(e: reed_io::error::IoError) -> FilterError {
    FilterError::Io(std::io::Error::other(e.to_string()))
}

fn ferr(err: EngineError) -> FilterError {
    FilterError::Io(std::io::Error::other(err.to_string()))
}

/// Runs the v2 filter.
///
/// # Errors
///
/// Spec compilation, transport, job and sink failures.
pub fn run_filter(
    input: &Path,
    spec: &FilterSpecV2,
    opts: &FilterOptions,
) -> Result<FilterReport, FilterError> {
    let t0 = Instant::now();
    let header = reed_io::record::record_header(input).map_err(FilterError::from)?;
    let n_samples_in = header.sample_names().len();
    let mut plan = reed_filter::plan::compile(spec, &header, &opts.command_line)?;
    let strategy = opts.strategy.unwrap_or(plan.strategy);
    plan.strategy = strategy;

    // ---- cohort resolution (§11.4) ----
    let (keep, keep_mask): (Arc<[bool]>, Option<SampleMask>) = match &spec.cohort {
        None => (Arc::from(vec![true; n_samples_in]), None),
        Some(CohortRuleV2::MaxMissing(threshold)) => {
            let partial = run_missing_pass(input, n_samples_in, opts)?;
            let keep_vec: Vec<bool> = if partial.n_records == 0 {
                vec![true; n_samples_in]
            } else {
                partial
                    .counts
                    .iter()
                    .map(|&m| m as f64 / partial.n_records as f64 <= *threshold)
                    .collect()
            };
            build_keep(keep_vec, &header, &mut plan)?
        }
        Some(CohortRuleV2::KeepList(path)) => {
            build_keep(keep_list_from_file(Path::new(path), &header, true)?, &header, &mut plan)?
        }
        Some(CohortRuleV2::RemoveList(path)) => {
            build_keep(keep_list_from_file(Path::new(path), &header, false)?, &header, &mut plan)?
        }
    };
    let n_samples_out = keep.iter().filter(|&&k| k).count();
    let subset = n_samples_out < n_samples_in;

    // ---- output header: prune samples, add recompute defs, provenance ----
    let out_header = if subset {
        reed_io::record::filter_header(&header, &keep)
    } else {
        header.clone()
    };
    let mut out_header = ensure_info_defs(out_header, &plan.recompute);
    let ascertainment = serde_json::to_string(&plan.ascertainment)
        .map_err(|e| FilterError::Spec(format!("ascertainment serialization: {e}")))?;
    for (key, value) in plan.provenance.header_lines(&ascertainment) {
        out_header
            .insert(key.as_str().parse().expect("invariant: valid key"), value.into())
            .map_err(|e| FilterError::Spec(format!("header insert {key}: {e}")))?;
    }
    let mut header_bytes = Vec::new();
    vcf::io::Writer::new(&mut header_bytes)
        .write_header(&out_header)
        .map_err(FilterError::from)?;

    // ---- output target ----
    let gz_out = opts.output.as_ref().is_some_and(|p| is_gz(p));
    let mut writer: Box<dyn Write + Send> = match &opts.output {
        None => Box::new(std::io::BufWriter::with_capacity(1 << 20, std::io::stdout())),
        Some(p) => Box::new(std::io::BufWriter::with_capacity(
            1 << 20,
            std::fs::File::create(p).map_err(FilterError::from)?,
        )),
    };
    let compressed = gz_out;

    // ---- the single pass ----
    let head_ids: Vec<String> = plan.head.iter().map(HeadRule::id).collect();
    let tail_ids: Vec<String> = plan.tail.iter().map(|t| t.id().to_string()).collect();
    let mask_ids: Vec<String> = plan
        .mask
        .iter()
        .map(|m| {
            let field = String::from_utf8_lossy(&m.field);
            format!("genotype.mask.{field}")
        })
        .collect();
    let contig_names: Arc<[Vec<u8>]> = Arc::from(
        header
            .contigs()
            .iter()
            .map(|(name, _)| name.as_bytes().to_vec())
            .collect::<Vec<_>>()
            .into_boxed_slice(),
    );
    let mut indexer = if gz_out && opts.write_index {
        use noodles::csi::binning_index::index::header::ReferenceSequenceNames;
        let names: ReferenceSequenceNames = contig_names
            .iter()
            .map(|v| bstr::BString::from(v.as_slice()))
            .collect();
        let header = noodles::csi::binning_index::index::Header::builder()
            .set_format(noodles::csi::binning_index::index::header::Format::Vcf)
            .set_reference_sequence_name_index(0)
            .set_start_position_index(1)
            .set_end_position_index(None)
            .set_line_comment_prefix(b'#')
            .set_reference_sequence_names(names)
            .build();
        let mut indexer = noodles::tabix::index::Indexer::default();
        indexer.set_header(header);
        Some(indexer)
    } else {
        None
    };
    let job = FilterJob {
        plan: Arc::new(plan),
        keep,
        keep_mask,
        n_samples_in,
        n_samples_out,
        compressed,
        contig_names: Arc::clone(&contig_names),
        level: opts.compress_level,
    };
    let slot: Arc<Mutex<Option<FunnelPartial>>> = Arc::new(Mutex::new(None));
    let index_path = opts
        .output
        .as_ref()
        .filter(|_| gz_out && opts.write_index)
        .map(|p| {
            let mut s = p.as_os_str().to_os_string();
            s.push(".tbi");
            PathBuf::from(s)
        });
    let mut sink = OutputSink {
        writer,
        indexer,
        index_path,
        coffset: 0,
        funnel: FunnelPartial::default(),
        slot: Arc::clone(&slot),
        compressed,
        contig_names: Arc::clone(&contig_names),
    };
    // Header first (as BGZF block(s) on the compressed path).
    sink.write_header_bytes(&header_bytes, compressed, opts.compress_level)?;

    let budget = ThreadBudget { threads: opts.threads.max(1), mem: reed_engine::budget::MemBudget::from_bytes(u64::MAX) };
    let (mut source, _facts) =
        BgzfChunkSource::open(std::io::BufReader::with_capacity(
            1 << 20,
            std::fs::File::open(input).map_err(FilterError::from)?,
        ))
        .map_err(|e: reed_io::error::IoError| FilterError::Io(std::io::Error::other(e.to_string())))?;
    let (_abort, abort_signal) = reed_engine::abort::abort_channel();
    run(&mut source as &mut dyn ChunkSource, &job, sink, &budget, &abort_signal).map_err(ferr)?;

    let funnel = slot
        .lock()
        .expect("invariant: funnel slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the funnel");

    // ---- report (§11.5) ----
    let mut funnel_stages: Vec<reed_filter::funnel::FunnelStage> = Vec::new();
    for (id, n) in head_ids.iter().zip(&funnel.head_hits) {
        funnel_stages.push(reed_filter::funnel::FunnelStage { rule: id.clone(), n_passed: *n });
    }
    for (id, n) in tail_ids.iter().zip(&funnel.tail_hits) {
        funnel_stages.push(reed_filter::funnel::FunnelStage { rule: id.clone(), n_passed: *n });
    }
    let masked_by_rule: Vec<reed_filter::funnel::MaskedBy> = mask_ids
        .iter()
        .zip(&funnel.mask_hits)
        .map(|(id, n)| reed_filter::funnel::MaskedBy { rule: id.clone(), n_masked: *n })
        .collect();
    let report = FilterReport {
        n_input_records: funnel.n_input,
        n_output_records: funnel.n_output,
        n_samples_in: n_samples_in as u32,
        n_samples_out: n_samples_out as u32,
        funnel: funnel_stages,
        n_masked_gt: funnel.masked_gt,
        masked_by_rule,
        policy_counters: funnel.policy,
        strategy: strategy.id().to_string(),
        elapsed_secs: t0.elapsed().as_secs_f64(),
    };
    if !opts.quiet {
        eprint_report(&report);
    }
    if let Some(path) = &opts.report_json {
        let json = serde_json::to_string_pretty(&report)
            .map_err(|e| FilterError::Spec(format!("report serialization: {e}")))?;
        std::fs::write(path, json).map_err(FilterError::from)?;
    }
    Ok(report)
}

fn is_gz(p: &Path) -> bool {
    p.extension().and_then(|e| e.to_str()) == Some("gz")
}

/// Builds the keep mask plus the editor's `SampleMask`, and completes the
/// ascertainment cohort dimension (§11.7).
fn build_keep(
    keep: Vec<bool>,
    header: &vcf::Header,
    plan: &mut FilterPlan,
) -> Result<(Arc<[bool]>, Option<SampleMask>), FilterError> {
    let n = header.sample_names().len();
    let ids: Vec<u32> = keep
        .iter()
        .enumerate()
        .filter(|(_, k)| **k)
        .map(|(i, _)| u32::try_from(i).expect("invariant: sample count < u32::MAX"))
        .collect();
    let mask = if ids.len() == n {
        None
    } else {
        Some(SampleMask::from_ids(
            n,
            ids.iter().copied().map(reed_core::SampleId),
        ).map_err(|e| FilterError::Spec(e.to_string()))?)
    };
    if ids.len() < n {
        use sha2::Digest as _;
        let mut hasher = sha2::Sha256::new();
        for (name, &k) in header.sample_names().iter().zip(&keep) {
            if k {
                hasher.update(name.as_bytes());
                hasher.update(b"\n");
            }
        }
        let keep_hash = format!("{:x}", hasher.finalize());
        plan.ascertainment = reed_filter::plan::ascertainment::with_cohort_subset(
            plan.ascertainment.clone(),
            &keep_hash,
            u32::try_from(ids.len()).expect("invariant: sample count < u32::MAX"),
        );
    }
    Ok((Arc::from(keep), mask))
}

/// Reads a keep/remove list (one sample per line, `#` comments).
fn keep_list_from_file(path: &Path, header: &vcf::Header, keep_mode: bool) -> Result<Vec<bool>, FilterError> {
    let text = std::fs::read_to_string(path).map_err(FilterError::from)?;
    let mut keep = vec![!keep_mode; header.sample_names().len()];
    let index: HashMap<&str, usize> = header
        .sample_names()
        .iter()
        .enumerate()
        .map(|(i, name)| (name.as_str(), i))
        .collect();
    for line in text.lines() {
        let name = line.trim();
        if name.is_empty() || name.starts_with('#') {
            continue;
        }
        let i = *index
            .get(name)
            .ok_or_else(|| FilterError::Spec(format!("sample {name:?} is not in the header")))?;
        keep[i] = keep_mode;
    }
    Ok(keep)
}

/// Adds header definitions for recomputed tags the input lacks (§11.6:
/// canonical Numbers — AC/HWE are `A`).
fn ensure_info_defs(mut header: vcf::Header, recompute: &reed_filter::plan::RecomputeFlags) -> vcf::Header {
    use noodles::vcf::header::record::value::map::info::{Number, Type};
    type InfoMap = noodles::vcf::header::record::value::Map<noodles::vcf::header::record::value::map::Info>;
    let mut add = |header: &mut vcf::Header, id: &str, number: Number, ty: Type, desc: &str| {
        if !header.infos().iter().any(|(k, _)| k.as_str() == id) {
            header
                .infos_mut()
                .insert(id.to_string(), InfoMap::new(number, ty, desc));
        }
    };
    if recompute.ac {
        add(&mut header, "AC", Number::A, Type::Integer, "Allele count in genotypes");
    }
    if recompute.an {
        add(&mut header, "AN", Number::Count(1), Type::Integer, "Total number of alleles in called genotypes");
    }
    if recompute.ns {
        add(&mut header, "NS", Number::Count(1), Type::Integer, "Number of samples with data");
    }
    if recompute.maf {
        add(&mut header, "MAF", Number::Count(1), Type::Float, "Frequency of the second most common allele");
    }
    if recompute.hwe {
        add(&mut header, "HWE", Number::A, Type::Float, "HWE test (PMID:15789306); 1=good, 0=bad");
    }
    header
}

// ---------------------------------------------------------------- pass0

/// Pass0 partial: per-sample fully-missing GT counts.
#[derive(Default)]
pub struct MissingPartial {
    n_records: u64,
    counts: Vec<u64>,
}

/// The GT-only missing-count job (`TwoPassReparse` pass0).
pub struct MissingJob {
    n_samples: usize,
}

impl ChunkJob for MissingJob {
    type Partial = MissingPartial;
    type Scratch = (Option<SampleScanner>, Vec<u64>);

    fn run_chunk(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let (scanner, miss) = scratch;
        let scanner = scanner.get_or_insert_with(|| SampleScanner::new(&[]));
        miss.clear();
        miss.resize(self.n_samples, 0);
        let mut n_records = 0_u64;
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            let Some(lr) = LineRef::parse(line) else {
                return Err(EngineError::Job("record line has fewer than 9 columns".into()));
            };
            let pos1 = parse_pos1(lr.pos_text(), &Locus { contig_name: String::new(), pos1: 0 }).map_err(|e| EngineError::Job(e.to_string()))?;
            let mut any = false;
            scanner
                .scan_fields(&lr, |i, gt, _| {
                    if let Some(gt) = gt {
                        if is_fully_missing(gt) {
                            miss[i] += 1;
                            any = true;
                        }
                    } else {
                        miss[i] += 1;
                        any = true;
                    }
                })
                .map_err(|e| EngineError::Job(e.to_string()))?;
            let _ = any;
            n_records += 1;
            let _ = pos1;
        }
        Ok(MissingPartial { n_records, counts: miss.clone() })
    }

    fn run_seam(&self, line: &[u8], scratch: &mut Self::Scratch) -> Result<Self::Partial, EngineError> {
        let (scanner, miss) = scratch;
        let scanner = scanner.get_or_insert_with(|| SampleScanner::new(&[]));
        miss.clear();
        miss.resize(self.n_samples, 0);
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::Job("seam record has fewer than 9 columns".into()));
        };
        scanner
            .scan_fields(&lr, |i, gt, _| {
                let missing = match gt {
                    Some(gt) => is_fully_missing(gt),
                    None => true,
                };
                if missing {
                    miss[i] += 1;
                }
            })
            .map_err(|e| EngineError::Job(e.to_string()))?;
        Ok(MissingPartial { n_records: 1, counts: miss.clone() })
    }
}

/// `./.` or `.` (both alleles missing).
fn is_fully_missing(gt: &[u8]) -> bool {
    gt.iter().all(|&b| b == b'.' || b == b'/' || b == b'|')
}

fn run_missing_pass(
    input: &Path,
    n_samples: usize,
    opts: &FilterOptions,
) -> Result<MissingPartial, FilterError> {
    let budget = ThreadBudget { threads: opts.threads.max(1), mem: reed_engine::budget::MemBudget::from_bytes(u64::MAX) };
    let (mut source, _facts) =
        BgzfChunkSource::open(std::io::BufReader::with_capacity(
            1 << 20,
            std::fs::File::open(input).map_err(FilterError::from)?,
        ))
        .map_err(|e: reed_io::error::IoError| FilterError::Io(std::io::Error::other(e.to_string())))?;
    let job = MissingJob { n_samples };
    let total: Arc<Mutex<Option<MissingPartial>>> = Arc::new(Mutex::new(None));
    let sink = SumSink { total: Arc::clone(&total) };
    let (_abort, abort_signal) = reed_engine::abort::abort_channel();
    run(&mut source as &mut dyn ChunkSource, &job, sink, &budget, &abort_signal).map_err(ferr);
    Ok(total
        .lock()
        .expect("invariant: missing slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the partial"))
}

/// Sums pass0 partials (integer sums are order-independent).
struct SumSink {
    total: Arc<Mutex<Option<MissingPartial>>>,
}

impl OrderedSink<MissingPartial> for SumSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut MissingPartial) -> Result<(), EngineError> {
        let mut guard = self
            .total
            .lock()
            .expect("invariant: missing slot not poisoned");
        match guard.as_mut() {
            Some(total) => {
                total.n_records += partial.n_records;
                for (t, m) in total.counts.iter_mut().zip(&partial.counts) {
                    *t += m;
                }
            }
            None => {
                let mut counts = vec![0_u64; partial.counts.len()];
                for (t, m) in counts.iter_mut().zip(&partial.counts) {
                    *t += m;
                }
                *guard = Some(MissingPartial { n_records: partial.n_records, counts });
            }
        }
        Ok(())
    }

    fn finish(self) -> Result<u64, EngineError> {
        let total = self.total.lock().expect("invariant: missing slot not poisoned");
        Ok(total.as_ref().map_or(0, |t| t.n_records))
    }
}

// ---------------------------------------------------------------- pass1

/// One compressed BGZF block plus its record index entries.
pub struct OutBlock {
    /// The deflated BGZF block (header + payload + footer).
    pub compressed: Vec<u8>,
    /// Uncompressed payload length (the block's virtual-offset span).
    pub uncompressed_len: u32,
    /// `(tid, pos1, ref_len, uoffset)` per output record, in order.
    pub records: Vec<(u16, u32, u32, u32)>,
}

/// The chunk partial: output bytes plus the funnel contribution.
#[derive(Default)]
pub struct FilterPartial {
    /// Compressed path: the chunk's BGZF blocks.
    pub blocks: Vec<OutBlock>,
    /// Plain path: the chunk's edited bytes.
    pub raw: Vec<u8>,
    /// The funnel contribution.
    pub funnel: FunnelPartial,
}

/// The pass1 chunk job: head → mask → tail → edit → worker-local deflate.
pub struct FilterJob {
    plan: Arc<FilterPlan>,
    keep: Arc<[bool]>,
    keep_mask: Option<SampleMask>,
    n_samples_in: usize,
    n_samples_out: usize,
    compressed: bool,
    contig_names: Arc<[Vec<u8>]>,
    level: u8,
}

/// Worker scratch (reused across chunks).
pub struct WorkerScratch {
    scanner: SampleScanner,
    head: HeadScratch,
    masked: Bitset,
    mask_values: Vec<Option<f64>>,
    counts: SiteCounts,
    info_set: Vec<(&'static [u8], InfoValueOut)>,
    edit_buf: Vec<u8>,
    block_buf: Vec<u8>,
    records: Vec<(u16, u32, u32, u32)>,
    blocks: Vec<OutBlock>,
    raw: Vec<u8>,
    deflater: Option<Deflater>,
    last_contig: (Vec<u8>, u16),
    funnel: FunnelPartial,
}

impl ChunkJob for FilterJob {
    type Partial = FilterPartial;
    type Scratch = Option<WorkerScratch>;

    fn run_chunk(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        scratch.funnel.reset(
            self.plan.head.len(),
            self.plan.tail.len(),
            self.plan.mask.len(),
        );
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            self.push_line(line, scratch)?;
        }
        self.flush_block(scratch)?;
        let partial = FilterPartial {
            blocks: std::mem::take(&mut scratch.blocks),
            raw: std::mem::take(&mut scratch.raw),
            funnel: std::mem::take(&mut scratch.funnel),
        };
        Ok(partial)
    }

    fn run_seam(&self, line: &[u8], scratch: &mut Self::Scratch) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        scratch.funnel.reset(
            self.plan.head.len(),
            self.plan.tail.len(),
            self.plan.mask.len(),
        );
        self.push_line(line, scratch)?;
        self.flush_block(scratch)?;
        Ok(FilterPartial {
            blocks: std::mem::take(&mut scratch.blocks),
            raw: std::mem::take(&mut scratch.raw),
            funnel: std::mem::take(&mut scratch.funnel),
        })
    }
}

impl FilterJob {
    fn build_scratch(&self) -> WorkerScratch {
        WorkerScratch {
            scanner: SampleScanner::new(
                &self
                    .plan
                    .int_fields
                    .iter()
                    .map(Vec::as_slice)
                    .collect::<Vec<_>>(),
            ),
            head: HeadScratch::new(),
            masked: Bitset::new(self.n_samples_in),
            mask_values: vec![None; self.plan.mask.len()],
            counts: SiteCounts::default(),
            info_set: Vec::new(),
            edit_buf: Vec::new(),
            block_buf: Vec::with_capacity(BGZF_BLOCK_SIZE + (1 << 14)),
            records: Vec::new(),
            blocks: Vec::new(),
            raw: Vec::new(),
            deflater: self.compressed.then(|| Deflater::new(i32::from(self.level)).expect("invariant: level")),
            last_contig: (Vec::new(), 0),
            funnel: FunnelPartial {
                head_hits: vec![0; self.plan.head.len()],
                tail_hits: vec![0; self.plan.tail.len()],
                mask_hits: vec![0; self.plan.mask.len()],
                ..FunnelPartial::default()
            },
        }
    }

    fn push_line(&self, line: &[u8], scratch: &mut WorkerScratch) -> Result<(), EngineError> {
        let plan = &self.plan;
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::Job("record line has fewer than 9 columns".into()));
        };
        let pos1 =
            parse_pos1(lr.pos_text(), &Locus { contig_name: String::new(), pos1: 0 }).map_err(|e| EngineError::Job(e.to_string()))?;

        // ---- head stage (§11.2: order is semantics) ----
        if !eval_head(&lr, pos1, &plan.head, &mut scratch.funnel.head_hits, &mut scratch.head) {
            scratch.funnel.n_input += 1;
            return Ok(());
        }

        // ---- scan + mask stage ----
        scratch.masked.clear_padding();
        for bit in 0..self.n_samples_in {
            scratch.masked.set(bit, false);
        }
        scratch.counts = SiteCounts::default();
        let keep = &self.keep;
        let mut masked_any = false;
        let mask_rules = &plan.mask;
        let scan_result = scratch.scanner.scan_fields(&lr, |i, gt, ints| {
            // One observation per mask RULE: rule i reads its own scanned
            // column, so several rules may share one field (DP<5 and
            // DP>100 both read DP).
            for (rule, v) in mask_rules.iter().zip(scratch.mask_values.iter_mut()) {
                *v = ints
                    .get(rule.int_column)
                    .and_then(|raw| raw.and_then(parse_f64_field));
            }
            let kept = keep.get(i).copied().unwrap_or(true);
            if !kept {
                return;
            }
            let hit = eval_mask(&scratch.mask_values, mask_rules);
            if let Some(rule_index) = hit {
                // Attribution: the first firing rule owns this sample
                // (called→missing transitions are counted below, so the
                // per-rule counters here include already-missing samples).
                scratch.funnel.mask_hits[rule_index] += 1;
                scratch.masked.set(i, true);
                masked_any = true;
            }
            // Called-ness from the GT text (half-calls are called).
            let (a, b) = match gt {
                Some(g) => parse_gt(g),
                None => (None, None),
            };
            let masked = hit.is_some();
            match (a, b) {
                (None, None) => {
                    // already missing: not a called→missing transition
                    scratch.counts.n_missing += 1;
                }
                (x, y) => {
                    if masked {
                        scratch.funnel.masked_gt += 1;
                        scratch.counts.n_missing += 1;
                    } else {
                        // count_site semantics, inlined
                        scratch.counts.ns += 1;
                        match (x, y) {
                            (Some(x), Some(y)) => {
                                scratch.counts.an += 2;
                                let n_alt = u64::from(x > 0) + u64::from(y > 0);
                                scratch.counts.ac += n_alt;
                                match n_alt {
                                    0 => scratch.counts.hom_ref += 1,
                                    1 => scratch.counts.het += 1,
                                    _ => scratch.counts.hom_alt += 1,
                                }
                            }
                            (Some(x), None) | (None, Some(x)) => {
                                scratch.counts.an += 1;
                                scratch.counts.ac += u64::from(x > 0);
                            }
                            (None, None) => unreachable!("handled above"),
                        }
                    }
                }
            }
        });
        if let Err(e) = scan_result {
            return Err(EngineError::Job(e.to_string()));
        }
        let _ = masked_any;

        // ---- tail stage ----
        if !eval_tail(&scratch.counts, self.n_samples_out, &plan.tail, &mut scratch.funnel.tail_hits) {
            scratch.funnel.n_input += 1;
            return Ok(());
        }

        // ---- edit + emit ----
        fill_info_set(&scratch.counts, &plan.recompute, &mut scratch.info_set);
        let stale: Vec<&[u8]> = plan.stale_tags.iter().map(Vec::as_slice).collect();
        let edit_plan = EditPlan {
            keep_samples: self.keep_mask.as_ref(),
            mask_gt: Some(&scratch.masked),
            info_remove: &stale,
            info_set: &scratch.info_set,
        };
        let gt_index = gt_subfield_index(lr.format_text());
        scratch.edit_buf.clear();
        edit_vcf_line(&lr, &edit_plan, gt_index, &mut scratch.edit_buf)
            .map_err(|e| EngineError::Job(e.to_string()))?;

        scratch.funnel.n_input += 1;
        scratch.funnel.n_output += 1;
        if self.compressed {
            // htslib inflates into a 64 KiB buffer: a member's uncompressed
            // payload must stay at or below BGZF_BLOCK_SIZE (§6.6), so cut
            // BEFORE appending a record that would overflow it.
            if scratch.block_buf.len() + scratch.edit_buf.len() > BGZF_BLOCK_SIZE {
                self.cut_block(scratch)?;
            }
            let uoffset = scratch.block_buf.len() as u32;
            let tid = self.tid_of(lr.chrom(), scratch);
            let ref_len = u32::try_from(lr.reference().len()).unwrap_or(1).max(1);
            scratch
                .records
                .push((tid, u32::try_from(pos1).unwrap_or(u32::MAX), ref_len, uoffset));
            scratch.block_buf.extend_from_slice(&scratch.edit_buf);
        } else {
            scratch.raw.extend_from_slice(&scratch.edit_buf);
        }
        Ok(())
    }

    fn tid_of(&self, chrom: &[u8], scratch: &mut WorkerScratch) -> u16 {
        if scratch.last_contig.0 == chrom {
            return scratch.last_contig.1;
        }
        let tid = self
            .contig_names
            .iter()
            .position(|name| name == chrom)
            .map(|i| i as u16)
            .unwrap_or(0);
        scratch.last_contig = (chrom.to_vec(), tid);
        tid
    }

    fn cut_block(&self, scratch: &mut WorkerScratch) -> Result<(), EngineError> {
        if scratch.block_buf.is_empty() {
            return Ok(());
        }
        let Some(deflater) = &mut scratch.deflater else {
            return Ok(());
        };
        let mut compressed = Vec::with_capacity(scratch.block_buf.len() / 2 + 64);
        deflater
            .deflate_block(&scratch.block_buf, &mut compressed)
            .map_err(|e| EngineError::Job(format!("bgzf deflate: {e}")))?;
        let uncompressed_len = u32::try_from(scratch.block_buf.len()).expect("invariant: < u32::MAX");
        scratch.blocks.push(OutBlock {
            compressed,
            uncompressed_len,
            records: std::mem::take(&mut scratch.records),
        });
        scratch.block_buf.clear();
        Ok(())
    }

    fn flush_block(&self, scratch: &mut WorkerScratch) -> Result<(), EngineError> {
        self.cut_block(scratch)
    }
}

fn gt_subfield_index(format_text: &[u8]) -> Option<u8> {
    format_text
        .split(|&b| b == b':')
        .position(|k| k == b"GT")
        .and_then(|i| u8::try_from(i).ok())
}

/// Parses a GT subfield: `(allele0, allele1)` with `None` = missing.
/// `None` alleles beyond two → policy error handled by the caller? The
/// polyploid case returns `(None, None)` + flag; for M5 the biallelic gate
/// makes it unreachable in practice.
fn parse_gt(gt: &[u8]) -> (Option<u16>, Option<u16>) {
    let mut out = [None, None];
    let mut seen = 0;
    for token in gt.split(|&b| b == b'/' || b == b'|') {
        if seen >= 2 {
            break;
        }
        if token != b"." && !token.is_empty() {
            let mut v = 0_u16;
            let ok = token.iter().all(|&b| b.is_ascii_digit());
            if ok {
                for &b in token {
                    v = v.saturating_mul(10).saturating_add(u16::from(b - b'0'));
                }
                out[seen] = Some(v);
            }
        }
        seen += 1;
    }
    if seen == 1 {
        // Haploid: normalize to a half-call.
        let only = out[0];
        out = [only, None];
    }
    (out[0], out[1])
}

fn parse_f64_field(raw: &[u8]) -> Option<f64> {
    if raw == b"." || raw.is_empty() {
        return None;
    }
    std::str::from_utf8(raw).ok().and_then(|s| s.parse().ok())
}

// ---------------------------------------------------------------- sink

/// The ordered output sink: appends worker blocks in chunk order, tracks
/// virtual offsets, and feeds the tabix indexer (§6.6).
struct OutputSink {
    writer: Box<dyn Write + Send>,
    indexer: Option<noodles::tabix::index::Indexer>,
    index_path: Option<PathBuf>,
    coffset: u64,
    funnel: FunnelPartial,
    slot: Arc<Mutex<Option<FunnelPartial>>>,
    compressed: bool,
    contig_names: Arc<[Vec<u8>]>,
}

impl OutputSink {
    /// Writes the (already rendered) header bytes as BGZF block(s).
    fn write_header_bytes(
        &mut self,
        header: &[u8],
        compressed: bool,
        level: u8,
    ) -> Result<(), FilterError> {
        if compressed {
            let mut deflater = Deflater::new(i32::from(level)).map_err(io_err)?;
            for chunk in header.chunks(BGZF_BLOCK_SIZE) {
                let mut block = Vec::with_capacity(chunk.len() / 2 + 64);
                deflater.deflate_block(chunk, &mut block).map_err(io_err)?;
                self.coffset += block.len() as u64;
                self.writer.write_all(&block).map_err(FilterError::from)?;
            }
        } else {
            self.writer.write_all(header).map_err(FilterError::from)?;
        }
        Ok(())
    }

    fn accept_block(&mut self, block: &OutBlock) -> Result<(), FilterError> {
        let coffset_before = self.coffset;
        self.writer.write_all(&block.compressed).map_err(FilterError::from)?;
        self.coffset += block.compressed.len() as u64;
        if let Some(indexer) = &mut self.indexer {
            // htslib tabix semantics, one add_record per output record:
            // span = [pos, pos + REF.len() - 1], chunk = [record voffset,
            // next record voffset] (block end for the last one).
            use noodles::bgzf::VirtualPosition;
            use noodles::core::Position as CorePosition;
            use noodles::csi::binning_index::index::reference_sequence::bin::Chunk;
            let block_end = (coffset_before << 16) | u64::from(block.uncompressed_len);
            for (i, &(tid, pos1, ref_len, uoffset)) in block.records.iter().enumerate() {
                let voff_start = (coffset_before << 16) | u64::from(uoffset);
                let voff_end = match block.records.get(i + 1) {
                    Some(&(_, _, _, next_uoffset)) => {
                        (coffset_before << 16) | u64::from(next_uoffset)
                    }
                    None => block_end,
                };
                let name = self
                    .contig_names
                    .get(usize::from(tid))
                    .map(|v| String::from_utf8_lossy(v).into_owned())
                    .unwrap_or_default();
                let pos = usize::try_from(pos1)
                    .ok()
                    .and_then(|p| CorePosition::try_from(p).ok())
                    .ok_or_else(|| FilterError::Spec(format!("index position {pos1}")))?;
                let end_pos1 = pos1.saturating_add(ref_len.saturating_sub(1));
                let end = usize::try_from(end_pos1)
                    .ok()
                    .and_then(|p| CorePosition::try_from(p).ok())
                    .ok_or_else(|| FilterError::Spec(format!("index end {end_pos1}")))?;
                indexer
                    .add_record(name.as_str(), pos, end, Chunk::new(
                        VirtualPosition::from(voff_start),
                        VirtualPosition::from(voff_end),
                    ))
                    .map_err(|e| FilterError::Spec(format!("index add_record: {e}")))?;
            }
        }
        Ok(())
    }
}

impl OrderedSink<FilterPartial> for OutputSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut FilterPartial) -> Result<(), EngineError> {
        self.funnel.merge(&partial.funnel);
        if self.compressed {
            for block in &partial.blocks {
                self.accept_block(block).map_err(|e| EngineError::Job(e.to_string()))?;
            }
        } else {
            self.writer.write_all(&partial.raw).map_err(|e| EngineError::Io(e.to_string()))?;
        }
        Ok(())
    }

    fn finish(mut self) -> Result<u64, EngineError> {
        if self.compressed {
            self.writer.write_all(&EOF_BLOCK).map_err(|e| EngineError::Io(e.to_string()))?;
        }
        self.writer.flush().map_err(|e| EngineError::Io(e.to_string()))?;
        if let (Some(indexer), Some(path)) = (self.indexer.take(), self.index_path.take()) {
            let index = indexer.build();
            noodles::tabix::fs::write(path, &index).map_err(|e| EngineError::Io(e.to_string()))?;
        }
        *self.slot.lock().expect("invariant: funnel slot not poisoned") =
            Some(std::mem::take(&mut self.funnel));
        Ok(self.funnel.n_output)
    }
}

fn eprint_report(report: &FilterReport) {
    eprintln!("reed filter: {} -> {} records ({} of {} samples), strategy {}, {:.1}s",
        report.n_input_records, report.n_output_records,
        report.n_samples_out, report.n_samples_in, report.strategy, report.elapsed_secs);
}

/// Reads a v2 spec from a TOML file.
///
/// # Errors
///
/// [`FilterError`] when the file cannot be read or parsed.
pub fn read_spec_toml(path: &Path) -> Result<FilterSpecV2, FilterError> {
    let text = std::fs::read_to_string(path).map_err(FilterError::from)?;
    FilterSpecV2::from_toml_str(&text)
}
