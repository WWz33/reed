//! `run_filter`: the L4 filter driver (§11.4). Only this layer sees both
//! the engine and the strategies, so the strategy choice lives here. The
//! compiled [`FilterPlan`] executes as `SinglePass` (one chunk pipeline:
//! head → mask → tail → edit → worker-local deflate) or `TwoPassReparse`
//! (pass0 = GT-only per-sample missing reduction; pass1 = SinglePass over
//! the surviving columns). The ordered sink appends worker-compressed
//! BGZF blocks in chunk order and feeds the tabix indexer while writing
//! (§6.6) — no second read for `--write-index`.
//!
//! The per-record stage sequence itself lives in
//! [`record_pipeline`](record_pipeline::RecordPipeline); the strategy files
//! only decide what to do with a record.

use std::collections::HashMap;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::Instant;

use reed_core::sample::SampleMask;
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{run, OrderedSink};
use reed_filter::funnel::{FilterReport, FunnelPartial};
use reed_filter::plan::{ExecStrategy, FilterPlan, HeadRule};
use reed_filter::spec::{FilterSpecV2, SampleList};
use reed_filter::FilterError;
use reed_core::ChunkId;
use reed_io::noodles_vcf as vcf;
use reed_io::write::index::Indexer;
use reed_io::write::ordered_output::OrderedOutput;

mod emit;
mod mid_stage;
mod record_pipeline;
mod reduce;
mod report;
mod single_pass;
mod sketch;

use emit::EmitJob;
use sketch::estimate_records;
use mid_stage::{mid_stage, ChunkDecision};
use reduce::run_missing_pass;
use single_pass::{FilterJob, FilterPartial};
use sketch::run_sketch_pass;

// Re-exported so `reed::filter::WriteIndex` keeps working (the type lives in
// reed-io's write layer, shared with `reed index`).
pub use reed_io::write::index::WriteIndex;

/// Run options.
#[derive(Clone, Debug, Default)]
pub struct FilterOptions {
    /// Worker threads (clamped to ≥ 1).
    pub threads: usize,
    /// Output path; `None` = stdout (plain VCF).
    pub output: Option<PathBuf>,
    /// Write a `.csi` (default) or `.tbi` index next to the output.
    pub write_index: WriteIndex,
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
    /// Sample lists naming samples absent from the header only warn
    /// (bcftools `--force-samples`).
    pub force_samples: bool,
    /// Index-accelerated regions (`-r/-R`): the reader seeks to the range the
    /// index maps them to; the spec's `targets` predicate drops records
    /// outside them.
    pub regions: Vec<reed_io::chunk::region::Region>,
    /// Directory for the stdin spill of multi-pass runs (default: the
    /// platform temp directory).
    pub tmp_dir: Option<PathBuf>,
    /// Explicit output compression (`-O v|z`); `None` infers it from the
    /// output path's extension.
    pub output_compressed: Option<bool>,
}

/// Pass1's source: the still-unconsumed first open (single-pass runs), else
/// a fresh open of the replay path (a file, or the stdin spill).
fn open_pass1(
    first: &mut Option<crate::input::OpenedInput>,
    replay_path: Option<&Path>,
    opts: &FilterOptions,
) -> Result<crate::input::OpenedInput, FilterError> {
    match first.take() {
        Some(opened) => Ok(opened),
        None => crate::input::open_regions(replay_path, &opts.regions).map_err(io_err),
    }
}

fn io_err(e: reed_io::error::IoError) -> FilterError {
    FilterError::Io(std::io::Error::other(e.to_string()))
}

/// Maps a transport error inside a chunk job.
fn job_err(e: reed_io::error::IoError) -> EngineError {
    EngineError::Job(e.to_string())
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
    input: Option<&Path>,
    spec: &FilterSpecV2,
    opts: &FilterOptions,
) -> Result<FilterReport, FilterError> {
    let t0 = Instant::now();
    if input.is_none() && !opts.regions.is_empty() {
        return Err(FilterError::Spec(
            "regions need a seekable file, not stdin".into(),
        ));
    }
    // ---- input: stdin streams straight into the single-pass pipeline;
    // multi-pass runs (`--mind`) spill it once and delete the spill on every
    // exit path (03-performance.md O5). The first open doubles as pass0's
    // source when nothing consumes it before pass1. ----
    let _spill: Option<crate::input::StdinSpill>;
    let mut replay_path: Option<PathBuf> = input.map(Path::to_path_buf);
    let mut first: Option<crate::input::OpenedInput>;
    match input {
        Some(p) => {
            first = Some(crate::input::open_regions(Some(p), &opts.regions).map_err(io_err)?);
            _spill = None;
        }
        None if spec.mind.is_some() => {
            let guard =
                crate::input::StdinSpill::new(opts.tmp_dir.as_deref()).map_err(io_err)?;
            let path = guard.path().to_path_buf();
            first = Some(crate::input::open_regions(Some(&path), &opts.regions).map_err(io_err)?);
            replay_path = Some(path);
            _spill = Some(guard);
        }
        None => {
            first = Some(crate::input::open(None).map_err(io_err)?);
            _spill = None;
        }
    }
    let facts = first
        .as_ref()
        .expect("invariant: the input was just opened")
        .facts
        .clone();
    let is_bcf = first.as_ref().is_some_and(|opened| opened.bcf.is_some());
    let header = facts.header.clone();
    let n_samples_in = header.sample_names().len();
    let mut plan = reed_filter::plan::compile(spec, &header, &opts.command_line)?;
    let strategy = match opts.strategy {
        Some(forced) => forced,
        // The sketch path's record-count estimate is line-based; BCF input
        // takes the two-pass path (doc_v3 02-bcf.md B1).
        None if is_bcf => ExecStrategy::TwoPassReparse,
        None => match &spec.mind {
            Some(_) if replay_path.is_some() => {
                // §11.4: the sketch path when the upper-bounded sketch fits
                // the memory budget (records extrapolated from the first
                // chunk's records/byte density).
                let estimated = estimate_records(
                    replay_path.as_deref().expect("invariant: guarded above"),
                )?;
                reed_filter::plan::strategy::choose_for_max_missing(
                    estimated,
                    n_samples_in,
                    plan.mask.len(),
                    reed_engine::budget::MemBudget::default_half_memory().bytes,
                )
            }
            _ => plan.strategy,
        },
    };
    plan.strategy = strategy;

    // ---- pass0 (the sketch path also captures the per-chunk sketches the
    // mid stage turns into emit decisions) ----
    let sketch_pass = if strategy == ExecStrategy::SketchThenEmit {
        let opened = first.take().expect("invariant: pass0 consumes the first open");
        Some(run_sketch_pass(opened, n_samples_in, &plan, opts)?)
    } else {
        None
    };

    // ---- static subset (`-s/-S/--keep/--remove`), then `--mind` (04 §C0b:
    // plink2's order; `--mind` counts missing calls after masking, over the
    // records the head rules keep) ----
    let mut keep_vec: Vec<bool> = match &spec.sample_select {
        None => vec![true; n_samples_in],
        Some(select) => {
            let mut keep = match &select.keep {
                None => vec![true; n_samples_in],
                Some(SampleList::Names(names)) => {
                    keep_list_from_names(names, &header, opts.force_samples)?
                }
                Some(SampleList::File(path)) => {
                    keep_list_from_file(Path::new(path), &header, true, opts.force_samples)?
                }
            };
            if let Some(remove) = &select.remove {
                let drop = match remove {
                    SampleList::Names(names) => {
                        keep_list_from_names(names, &header, opts.force_samples)?
                    }
                    SampleList::File(path) => {
                        keep_list_from_file(Path::new(path), &header, true, opts.force_samples)?
                    }
                };
                for (keep, drop) in keep.iter_mut().zip(&drop) {
                    *keep &= !*drop;
                }
            }
            keep
        }
    };
    if let Some(threshold) = spec.mind {
        let missing_partial;
        let (n_records, counts): (u64, &[u64]) = match &sketch_pass {
            Some(sketch) => (sketch.n_records, &sketch.cohort_miss),
            None => {
                let opened = first.take().expect("invariant: pass0 consumes the first open");
                missing_partial = run_missing_pass(opened, &plan, opts)?;
                (missing_partial.n_records, &missing_partial.counts)
            }
        };
        if n_records > 0 {
            for (keep, &miss) in keep_vec.iter_mut().zip(counts) {
                if miss as f64 / n_records as f64 > threshold {
                    *keep = false;
                }
            }
        }
    }
    let (keep, keep_mask) = build_keep(keep_vec, &header, &mut plan)?;
    let n_samples_out = keep.iter().filter(|&&k| k).count();
    let subset = n_samples_out < n_samples_in;
    // bcftools `view` updates INFO/AC,INFO/AN when samples are dropped; only
    // masking leaves them stale, like `bcftools filter -S .`
    // (02-cli.md §3.5). `-I/--no-update` suppresses the update.
    if subset && !plan.no_update {
        plan.recompute.ac = true;
        plan.recompute.an = true;
        for tag in ["AC", "AN"] {
            if !plan.stale_tags.iter().any(|t| t == tag.as_bytes()) {
                plan.stale_tags.push(tag.as_bytes().to_vec());
            }
        }
    }

    // ---- mid stage (§11.4): keep ∧ planes → per-chunk emit decisions ----
    let sketch_out: Option<(FunnelPartial, std::collections::HashMap<u64, Arc<ChunkDecision>>)> =
        sketch_pass.map(|sketch| {
            let mut funnel = sketch.funnel;
            let map = mid_stage(&sketch.sketches, &keep, &plan, n_samples_out, &mut funnel);
            (funnel, map)
        });

    // ---- output header: prune samples, add recompute defs, provenance ----
    let out_header = if subset {
        facts.pruned_header(&keep)
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
    // Extension inference (01-review D12): .vcf.gz/.vcf.bgz (or any .gz/
    // .bgz) compress; BCF output is refused until the BCF encoder exists —
    // a silent text-VCF-in-.bcf is worse than a loud error.
    let gz_out = match opts.output_compressed {
        Some(compressed) => compressed,
        None => match &opts.output {
            None => false,
            Some(p) => match p.extension().and_then(|e| e.to_str()) {
                Some("gz") | Some("bgz") => true,
                Some("bcf") => {
                    return Err(FilterError::Spec(
                        "BCF output is not supported yet; use a .vcf.gz output".into(),
                    ))
                }
                _ => false,
            },
        },
    };
    let writer: Box<dyn Write + Send> = match &opts.output {
        None => Box::new(std::io::BufWriter::with_capacity(1 << 20, std::io::stdout())),
        Some(p) => Box::new(std::io::BufWriter::with_capacity(
            1 << 20,
            std::fs::File::create(p).map_err(FilterError::from)?,
        )),
    };
    let compressed = gz_out;

    // ---- pass1 ----
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
    let indexer = if gz_out && opts.write_index.is_some() {
        Some(Indexer::new(opts.write_index, &contig_names).map_err(io_err)?)
    } else {
        None
    };
    let index_path = opts
        .output
        .as_ref()
        .filter(|_| gz_out && opts.write_index.is_some())
        .map(|p| opts.write_index.path_for(p));
    let budget = ThreadBudget { threads: opts.threads.max(1), mem: reed_engine::budget::MemBudget::default_half_memory() };

    // The report funnel: the FilterJob paths collect it in the sink slot;
    // the sketch path assembled it during the mid stage.
    let funnel = if let Some((sketch_funnel, decisions)) = sketch_out {
        let job = EmitJob {
            plan: Arc::new(plan),
            decisions: Arc::new(decisions),
            keep_mask,
            n_samples_in,
            compressed,
            contig_names: Arc::clone(&contig_names),
            level: opts.compress_level,
        };
        let slot: Arc<Mutex<Option<FunnelPartial>>> = Arc::new(Mutex::new(None));
        let mut sink = OutputSink {
            out: OrderedOutput::new(
                writer,
                compressed,
                opts.compress_level,
                indexer,
                index_path,
                Arc::clone(&contig_names),
            ),
            funnel: FunnelPartial::default(),
            slot: Arc::clone(&slot),
        };
        // Header first (as BGZF block(s) on the compressed path).
        sink.out.write_header(&header_bytes).map_err(io_err)?;
        let mut opened = open_pass1(&mut first, replay_path.as_deref(), opts)?;
        let abort = reed_engine::abort::abort();
        let job = crate::bcf_input::InputJob::new(job, &opened, &abort);
        let report = run(opened.source.as_mut(), &job, sink, &budget, &abort).map_err(ferr)?;
        job.finish(report.n_chunks, opts.regions.is_empty()).map_err(ferr)?;
        sketch_funnel
    } else {
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
        let mut sink = OutputSink {
            out: OrderedOutput::new(
                writer,
                compressed,
                opts.compress_level,
                indexer,
                index_path,
                Arc::clone(&contig_names),
            ),
            funnel: FunnelPartial::default(),
            slot: Arc::clone(&slot),
        };
        sink.out.write_header(&header_bytes).map_err(io_err)?;
        let mut opened = open_pass1(&mut first, replay_path.as_deref(), opts)?;
        let abort = reed_engine::abort::abort();
        let job = crate::bcf_input::InputJob::new(job, &opened, &abort);
        let report = run(opened.source.as_mut(), &job, sink, &budget, &abort).map_err(ferr)?;
        job.finish(report.n_chunks, opts.regions.is_empty()).map_err(ferr)?;
        slot
            .lock()
            .expect("invariant: funnel slot not poisoned")
            .take()
            .expect("invariant: sink finish stored the funnel")
    };

    // ---- report (§11.5) ----
    let funnel_stages = report::funnel_stages(&head_ids, &tail_ids, &funnel);
    let masked_by_rule = report::masked_by_rule(&mask_ids, &funnel);
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
        report::eprint_report(&report);
    }
    if let Some(path) = &opts.report_json {
        let json = serde_json::to_string_pretty(&report)
            .map_err(|e| FilterError::Spec(format!("report serialization: {e}")))?;
        std::fs::write(path, json).map_err(FilterError::from)?;
    }
    Ok(report)
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

/// Reads a keep/remove list (one sample per line, `#` comments). Unknown
/// samples are an error unless `force` (bcftools `--force-samples`).
fn keep_list_from_file(
    path: &Path,
    header: &vcf::Header,
    keep_mode: bool,
    force: bool,
) -> Result<Vec<bool>, FilterError> {
    let text = std::fs::read_to_string(path).map_err(FilterError::from)?;
    let names: Vec<String> = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(str::to_string)
        .collect();
    keep_list_from_names(&names, header, force).map(|mut keep| {
        if !keep_mode {
            for k in &mut keep {
                *k = !*k;
            }
        }
        keep
    })
}

/// Builds the keep mask from an inline sample list.
fn keep_list_from_names(
    names: &[String],
    header: &vcf::Header,
    force: bool,
) -> Result<Vec<bool>, FilterError> {
    let index: HashMap<&str, usize> = header
        .sample_names()
        .iter()
        .enumerate()
        .map(|(i, name)| (name.as_str(), i))
        .collect();
    let mut keep = vec![false; header.sample_names().len()];
    for name in names {
        match index.get(name.as_str()) {
            Some(&i) => keep[i] = true,
            None if force => eprintln!("reed filter: warning: sample {name:?} is not in the header"),
            None => {
                return Err(FilterError::Spec(format!("sample {name:?} is not in the header")))
            }
        }
    }
    Ok(keep)
}

/// Adds header definitions for recomputed tags the input lacks (§11.6:
/// canonical Numbers — AC/HWE are `A`).
fn ensure_info_defs(mut header: vcf::Header, recompute: &reed_filter::plan::RecomputeFlags) -> vcf::Header {
    use noodles::vcf::header::record::value::map::info::{Number, Type};
    type InfoMap = noodles::vcf::header::record::value::Map<noodles::vcf::header::record::value::map::Info>;
    let add = |header: &mut vcf::Header, id: &str, number: Number, ty: Type, desc: &str| {
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

/// The ordered output sink: [`OrderedOutput`] plus the funnel counters the
/// filter reports (§6.6, §11.5).
struct OutputSink {
    out: OrderedOutput,
    funnel: FunnelPartial,
    slot: Arc<Mutex<Option<FunnelPartial>>>,
}


impl OrderedSink<FilterPartial> for OutputSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut FilterPartial) -> Result<(), EngineError> {
        self.funnel.merge(&partial.funnel);
        self.out.accept_chunk(&partial.output).map_err(job_err)
    }

    fn finish(self) -> Result<u64, EngineError> {
        let OutputSink { out, mut funnel, slot } = self;
        out.finish().map_err(job_err)?;
        let n_output = funnel.n_output;
        *slot.lock().expect("invariant: funnel slot not poisoned") = Some(std::mem::take(&mut funnel));
        Ok(n_output)
    }
}
