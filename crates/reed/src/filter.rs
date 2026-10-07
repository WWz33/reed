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

use reed_core::count::SiteCounts;
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
use reed_filter::spec::{CohortRule as CohortRuleV2, FilterSpecV2};
use reed_filter::FilterError;
use reed_io::bgzf::block::EOF_BLOCK;
use reed_io::bgzf::deflate::Deflater;
use reed_io::chunk::ChunkSource;
use noodles;
use reed_io::noodles_vcf as vcf;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::line_edit::{edit_vcf_line, EditPlan, InfoValueOut};
use reed_core::genotype::{GenotypePolicies, PolicyCounters, STATE_ALT, STATE_MISSING};
use reed_io::vcf::genotype_scan::{INT_MISSING, GenotypeScanner};
use reed_io::vcf::site::{ContigTable, SiteView};

use crate::index::Indexer;
// Re-exported so `reed::filter::WriteIndex` keeps working (the type now
// lives with the `reed index` command it is shared with).
pub use crate::index::WriteIndex;

/// BGZF maximum uncompressed payload per block (htslib `BGZF_BLOCK_SIZE`).
const BGZF_BLOCK_SIZE: usize = 65_280;

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
    input: Option<&Path>,
    spec: &FilterSpecV2,
    opts: &FilterOptions,
) -> Result<FilterReport, FilterError> {
    let t0 = Instant::now();
    // One open for the header (the pass loops below open again per pass).
    let (_source, facts) = crate::input::open_regions(input, &opts.regions).map_err(io_err)?;
    let header = facts.header.clone();
    let n_samples_in = header.sample_names().len();
    let mut plan = reed_filter::plan::compile(spec, &header, &opts.command_line)?;
    let strategy = match opts.strategy {
        Some(forced) => forced,
        None => match &spec.cohort {
            Some(CohortRuleV2::MaxMissing(_)) if input.is_some() => {
                // §11.4: the sketch path when the upper-bounded sketch fits
                // the memory budget (records extrapolated from the first
                // chunk's records/byte density).
                let estimated = estimate_records(input.expect("invariant: guarded above"))?;
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

    // ---- pass0 (cohort resolution; the sketch path also captures the
    // per-chunk sketches the mid stage turns into emit decisions) ----
    let sketch_pass = if strategy == ExecStrategy::SketchThenEmit {
        Some(run_sketch_pass(input, n_samples_in, &plan, opts)?)
    } else {
        None
    };
    let (keep, keep_mask): (Arc<[bool]>, Option<SampleMask>) = match &spec.cohort {
        None => (Arc::from(vec![true; n_samples_in]), None),
        Some(CohortRuleV2::MaxMissing(threshold)) => {
            let missing_partial;
            let (n_records, counts): (u64, &[u64]) = match &sketch_pass {
                Some(sketch) => (sketch.n_records, &sketch.cohort_miss),
                None => {
                    missing_partial = run_missing_pass(input, n_samples_in, plan.policies, opts)?;
                    (missing_partial.n_records, &missing_partial.counts)
                }
            };
            let keep_vec: Vec<bool> = if n_records == 0 {
                vec![true; n_samples_in]
            } else {
                counts
                    .iter()
                    .map(|&m| m as f64 / n_records as f64 <= *threshold)
                    .collect()
            };
            build_keep(keep_vec, &header, &mut plan)?
        }
        Some(CohortRuleV2::KeepList(path)) => build_keep(
            keep_list_from_file(Path::new(path), &header, true, opts.force_samples)?,
            &header,
            &mut plan,
        )?,
        Some(CohortRuleV2::KeepNames(names)) => {
            build_keep(keep_list_from_names(names, &header, opts.force_samples)?, &header, &mut plan)?
        }
        Some(CohortRuleV2::RemoveList(path)) => build_keep(
            keep_list_from_file(Path::new(path), &header, false, opts.force_samples)?,
            &header,
            &mut plan,
        )?,
    };
    let n_samples_out = keep.iter().filter(|&&k| k).count();
    let subset = n_samples_out < n_samples_in;

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
    let gz_out = match &opts.output {
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
    };
    let mut writer: Box<dyn Write + Send> = match &opts.output {
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
    let mut indexer = if gz_out && opts.write_index.is_some() {
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
        match opts.write_index {
            WriteIndex::Tbi => {
                let mut indexer = noodles::tabix::index::Indexer::default();
                indexer.set_header(header);
                Some(Indexer::Tbi(Box::new(indexer)))
            }
            WriteIndex::Csi | WriteIndex::None => {
                use noodles::csi::binning_index::index::reference_sequence::index::BinnedIndex;
                // depth 6 (shift 32) so CSI covers every VCF-legal position;
                // the noodles default (depth 5) is TBI-sized (2^29).
                let indexer = noodles::csi::binning_index::indexer::Indexer::<BinnedIndex>::new(14, 6)
                    .expect("invariant: min_shift/depth are valid")
                    .set_header(header);
                Some(Indexer::Csi(Box::new(indexer), contig_names.len()))
            }
        }
    } else {
        None
    };
    let index_path = opts
        .output
        .as_ref()
        .filter(|_| gz_out && opts.write_index.is_some())
        .map(|p| {
            let suffix = match opts.write_index {
                WriteIndex::Tbi => ".tbi",
                WriteIndex::Csi | WriteIndex::None => ".csi",
            };
            let mut s = p.as_os_str().to_os_string();
            s.push(suffix);
            PathBuf::from(s)
        });
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
            writer,
            indexer,
            index_path,
            coffset: 0,
            funnel: FunnelPartial::default(),
            slot: Arc::clone(&slot),
            compressed,
            contig_names: Arc::clone(&contig_names),
            pending_record: None,
            last_block_end: 0,
        };
        // Header first (as BGZF block(s) on the compressed path).
        sink.write_header_bytes(&header_bytes, compressed, opts.compress_level)?;
        let (mut source, _facts) = crate::input::open_regions(input, &opts.regions).map_err(io_err)?;
        let abort = reed_engine::abort::abort();
        run(source.as_mut(), &job, sink, &budget, &abort).map_err(ferr)?;
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
            writer,
            indexer,
            index_path,
            coffset: 0,
            funnel: FunnelPartial::default(),
            slot: Arc::clone(&slot),
            compressed,
            contig_names: Arc::clone(&contig_names),
            pending_record: None,
            last_block_end: 0,
        };
        sink.write_header_bytes(&header_bytes, compressed, opts.compress_level)?;
        let (mut source, _facts) = crate::input::open_regions(input, &opts.regions).map_err(io_err)?;
        let abort = reed_engine::abort::abort();
        run(source.as_mut(), &job, sink, &budget, &abort).map_err(ferr)?;
        slot
            .lock()
            .expect("invariant: funnel slot not poisoned")
            .take()
            .expect("invariant: sink finish stored the funnel")
    };

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
    policies: GenotypePolicies,
}

impl ChunkJob for MissingJob {
    type Partial = MissingPartial;
    type Scratch = (Option<GenotypeScanner>, Vec<u64>, PolicyCounters);

    fn run_chunk(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let (scanner, miss, counters) = scratch;
        let n_samples = self.n_samples;
        let scanner = scanner.get_or_insert_with(|| {
            GenotypeScanner::new(n_samples, &[], &[], self.policies)
        });
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
            let site = SiteView::parse(&lr).map_err(|e| EngineError::Job(e.to_string()))?;
            *counters = PolicyCounters::default();
            let out = scanner
                .scan(&lr, &site, counters)
                .map_err(|e| EngineError::Job(e.to_string()))?;
            for (i, &code) in out.codes.iter().enumerate() {
                if code.is_fully_missing() {
                    miss[i] += 1;
                }
            }
            n_records += 1;
        }
        Ok(MissingPartial { n_records, counts: miss.clone() })
    }

    fn run_seam(&self, _id: ChunkId, line: &[u8], scratch: &mut Self::Scratch) -> Result<Self::Partial, EngineError> {
        let (scanner, miss, counters) = scratch;
        let n_samples = self.n_samples;
        let scanner = scanner.get_or_insert_with(|| {
            GenotypeScanner::new(n_samples, &[], &[], self.policies)
        });
        miss.clear();
        miss.resize(self.n_samples, 0);
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::Job("seam record has fewer than 9 columns".into()));
        };
        let site = SiteView::parse(&lr).map_err(|e| EngineError::Job(e.to_string()))?;
        *counters = PolicyCounters::default();
        let out = scanner
            .scan(&lr, &site, counters)
            .map_err(|e| EngineError::Job(e.to_string()))?;
        for (i, &code) in out.codes.iter().enumerate() {
            if code.is_fully_missing() {
                miss[i] += 1;
            }
        }
        Ok(MissingPartial { n_records: 1, counts: miss.clone() })
    }
}

fn run_missing_pass(
    input: Option<&Path>,
    n_samples: usize,
    policies: GenotypePolicies,
    opts: &FilterOptions,
) -> Result<MissingPartial, FilterError> {
    let budget = ThreadBudget { threads: opts.threads.max(1), mem: reed_engine::budget::MemBudget::default_half_memory() };
    let (mut source, _facts) = crate::input::open_regions(input, &opts.regions).map_err(io_err)?;
    let job = MissingJob { n_samples, policies };
    let total: Arc<Mutex<Option<MissingPartial>>> = Arc::new(Mutex::new(None));
    let sink = SumSink { total: Arc::clone(&total) };
    let abort = reed_engine::abort::abort();
    run(source.as_mut(), &job, sink, &budget, &abort).map_err(ferr);
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
    /// `(tid, pos1, ref_len, uoffset)` per output record, in order. `tid`
    /// is the header contig's index (u32 — no truncation); a record whose
    /// bytes span into the next block starts in *this* one, at `uoffset`.
    pub records: Vec<(u32, u64, u32, u32)>,
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
    scanner: GenotypeScanner,
    /// Contig table for the record's tid (last-hit cached).
    contigs: ContigTable,
    /// Per-record decode policy triggers, folded into the funnel.
    counters: PolicyCounters,
    head: HeadScratch,
    masked: Bitset,
    mask_values: Vec<Option<f64>>,
    counts: SiteCounts,
    info_set: Vec<(&'static [u8], InfoValueOut)>,
    edit_buf: Vec<u8>,
    block_buf: Vec<u8>,
    records: Vec<(u32, u64, u32, u32)>,
    blocks: Vec<OutBlock>,
    raw: Vec<u8>,
    deflater: Option<Deflater>,
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

    fn run_seam(&self, _id: ChunkId, line: &[u8], scratch: &mut Self::Scratch) -> Result<Self::Partial, EngineError> {
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
        let int_keys: Vec<&[u8]> = self.plan.int_fields.iter().map(Vec::as_slice).collect();
        let float_keys: Vec<&[u8]> = self.plan.float_fields.iter().map(Vec::as_slice).collect();
        WorkerScratch {
            scanner: GenotypeScanner::new(
                self.n_samples_in,
                &int_keys,
                &float_keys,
                self.plan.policies,
            ),
            contigs: ContigTable::new(self.contig_names.iter().cloned()),
            counters: PolicyCounters::default(),
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
        // The contig table is the index's tid space: a record on an
        // undeclared contig has no addressable tid and must fail loudly.
        let site = SiteView::parse_declared(&lr, &mut scratch.contigs)
            .map_err(|e| EngineError::Job(e.to_string()))?;
        let pos1 = site.pos1;

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
        scratch.counters = PolicyCounters::default();
        let scan_result = scratch.scanner.scan(&lr, &site, &mut scratch.counters);
        let out = match scan_result {
            Ok(out) => out,
            Err(e) => return Err(EngineError::Job(e.to_string())),
        };
        scratch.funnel.policy.add_decode(&scratch.counters);
        let keep = &self.keep;
        let mask_rules = &plan.mask;
        for (i, &code) in out.codes.iter().enumerate() {
            // One observation per mask RULE: rule i reads its own scanned
            // column, so several rules may share one field (DP<5 and
            // DP>100 both read DP).
            for (rule, v) in mask_rules.iter().zip(scratch.mask_values.iter_mut()) {
                *v = match rule.float_column {
                    Some(j) => out
                        .float_columns
                        .get(j)
                        .and_then(|column| column.get(i))
                        .copied()
                        .filter(|x| !x.is_nan())
                        .map(f64::from),
                    None => out
                        .int_columns
                        .get(rule.int_column)
                        .and_then(|column| column.get(i))
                        .copied()
                        .filter(|&x| x != INT_MISSING)
                        .map(f64::from),
                };
            }
            let kept = keep.get(i).copied().unwrap_or(true);
            if !kept {
                continue;
            }
            let hit = eval_mask(&scratch.mask_values, mask_rules);
            if let Some(rule_index) = hit {
                // Attribution: the first firing rule owns this sample
                // (called→missing transitions are counted below, so the
                // per-rule counters here include already-missing samples).
                scratch.funnel.mask_hits[rule_index] += 1;
                scratch.masked.set(i, true);
            }
            let masked = hit.is_some();
            if !code.is_called_any() {
                // already missing: not a called→missing transition
                scratch.counts.n_missing += 1;
            } else if masked {
                scratch.funnel.masked_gt += 1;
                scratch.counts.n_missing += 1;
            } else {
                // count_site semantics, from the decoded slots.
                let (s0, s1) = (code.slot_state(0), code.slot_state(1));
                let n_alt = u64::from(s0 == STATE_ALT) + u64::from(s1 == STATE_ALT);
                let n_miss = u64::from(s0 == STATE_MISSING) + u64::from(s1 == STATE_MISSING);
                scratch.counts.ns += 1;
                scratch.counts.an += 2 - n_miss;
                scratch.counts.ac += n_alt;
                match (n_miss, n_alt) {
                    (0, 0) => scratch.counts.hom_ref += 1,
                    (0, 1) => scratch.counts.het += 1,
                    (0, _) => scratch.counts.hom_alt += 1,
                    // Half-calls feed AN/AC/NS but no genotype class.
                    _ => {}
                }
            }
        }

        // ---- tail stage ----
        let tail_values = reed_filter::expr::SiteValues::from_line(
            &lr,
            pos1,
            Some(reed_filter::expr::GtCounts {
                ac: scratch.counts.ac,
                an: scratch.counts.an,
                n_missing: scratch.counts.n_missing,
                n_samples: self.n_samples_out as u64,
            }),
        );
        if !eval_tail(&tail_values, &scratch.counts, self.n_samples_out, &plan.tail, &mut scratch.funnel.tail_hits) {
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
            // Blocks are cut on the byte boundary, not the record boundary:
            // a record wider than the 65,280-byte member limit must span
            // blocks instead of corrupting one (htslib behavior). The
            // record's index entry anchors at its start offset.
            let tid = site
                .contig
                .expect("invariant: parse_declared resolved the contig");
            let ref_len = u32::try_from(site.reference.len()).unwrap_or(1).max(1);
            let uoffset = u32::try_from(scratch.block_buf.len())
                .expect("invariant: pending block below u32::MAX");
            scratch.records.push((tid, pos1, ref_len, uoffset));
            scratch.block_buf.extend_from_slice(&scratch.edit_buf);
            while scratch.block_buf.len() >= BGZF_BLOCK_SIZE {
                self.cut_block(scratch)?;
            }
        } else {
            scratch.raw.extend_from_slice(&scratch.edit_buf);
        }
        Ok(())
    }

    /// Deflates one block of at most [`BGZF_BLOCK_SIZE`] bytes. Records that
    /// start before the cut stay with this block; later starts are re-based
    /// into the continuation block.
    fn cut_block(&self, scratch: &mut WorkerScratch) -> Result<(), EngineError> {
        if scratch.block_buf.is_empty() {
            return Ok(());
        }
        let Some(deflater) = &mut scratch.deflater else {
            return Ok(());
        };
        let cut = scratch.block_buf.len().min(BGZF_BLOCK_SIZE);
        let cut32 = u32::try_from(cut).expect("invariant: cut ≤ 65,280");
        let mut compressed = Vec::with_capacity(cut / 2 + 64);
        deflater
            .deflate_block(&scratch.block_buf[..cut], &mut compressed)
            .map_err(|e| EngineError::Job(format!("bgzf deflate: {e}")))?;
        let boundary = scratch
            .records
            .partition_point(|&(_, _, _, uoffset)| uoffset < cut32);
        let mut records = std::mem::take(&mut scratch.records);
        let mut continuation = records.split_off(boundary);
        for record in continuation.iter_mut() {
            record.3 -= cut32;
        }
        scratch.records = continuation;
        scratch.blocks.push(OutBlock {
            compressed,
            uncompressed_len: cut32,
            records,
        });
        scratch.block_buf.drain(..cut);
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

// ---------------------------------------------------------------- sink

/// The ordered output sink: appends worker blocks in chunk order, tracks
/// virtual offsets, and feeds the indexer while writing (§6.6).
struct OutputSink {
    writer: Box<dyn Write + Send>,
    indexer: Option<Indexer>,
    index_path: Option<PathBuf>,
    coffset: u64,
    funnel: FunnelPartial,
    slot: Arc<Mutex<Option<FunnelPartial>>>,
    compressed: bool,
    contig_names: Arc<[Vec<u8>]>,
    /// The record whose index chunk is still open: its end is the next
    /// record's start (which may sit in a later block when records span
    /// block boundaries).
    pending_record: Option<(u32, u64, u32, u64)>,
    /// Virtual offset just past the last written block.
    last_block_end: u64,
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

    fn accept_block(&mut self, block: &OutBlock) -> Result<(), EngineError> {
        let coffset_before = self.coffset;
        self.writer
            .write_all(&block.compressed)
            .map_err(|e| EngineError::Io(e.to_string()))?;
        self.coffset += block.compressed.len() as u64;
        self.last_block_end = (coffset_before << 16) | u64::from(block.uncompressed_len);
        if self.indexer.is_some() {
            // Index chunks follow the record stream, not the block stream:
            // a record's chunk ends where the next record starts, even when
            // that is a later block.
            for &(tid, pos1, ref_len, uoffset) in &block.records {
                let voff_start = (coffset_before << 16) | u64::from(uoffset);
                if let Some(pending) = self.pending_record.take() {
                    self.close_index_record(pending, voff_start)?;
                }
                self.pending_record = Some((tid, pos1, ref_len, voff_start));
            }
        }
        Ok(())
    }

    /// Closes one record's index chunk `[start, end)` against the indexer.
    fn close_index_record(
        &mut self,
        pending: (u32, u64, u32, u64),
        voff_end: u64,
    ) -> Result<(), EngineError> {
        let Some(indexer) = &mut self.indexer else {
            return Ok(());
        };
        let (tid, pos1, ref_len, voff_start) = pending;
        let name = self
            .contig_names
            .get(tid as usize)
            .map(|v| String::from_utf8_lossy(v).into_owned())
            .unwrap_or_default();
        indexer.add(&name, tid, pos1, ref_len, voff_start, voff_end)
    }
}

impl OrderedSink<FilterPartial> for OutputSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut FilterPartial) -> Result<(), EngineError> {
        self.funnel.merge(&partial.funnel);
        if self.compressed {
            for block in &partial.blocks {
                self.accept_block(block)?;
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
        if let Some(pending) = self.pending_record.take() {
            // The last record's chunk runs to the end of the last block.
            self.close_index_record(pending, self.last_block_end)?;
        }
        if let (Some(indexer), Some(path)) = (self.indexer.take(), self.index_path.take()) {
            match indexer {
                Indexer::Csi(indexer, n_ref) => noodles::csi::fs::write(path, &indexer.build(n_ref))
                    .map_err(|e| EngineError::Io(e.to_string()))?,
                Indexer::Tbi(indexer) => {
                    noodles::tabix::fs::write(path, &indexer.build())
                        .map_err(|e| EngineError::Io(e.to_string()))?
                }
            }
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

// ------------------------------------------------------------ SketchThenEmit

/// One survivor's pass1 emit data (mid-stage output, §11.4).
#[derive(Clone, Debug)]
pub struct SiteEmit {
    /// Kept ∧ mask-fired sample bitmap (the EditPlan's `mask_gt`).
    pub mask_words: Vec<u64>,
    /// Recomputed INFO values in fill-tags order.
    pub infos: Vec<(&'static [u8], InfoValueOut)>,
}

/// A chunk's pass1 decisions, produced by the mid stage.
#[derive(Clone, Debug, Default)]
pub struct ChunkDecision {
    /// Records in the chunk (body + trailing seam).
    pub n_records: u64,
    /// Bit per record: passed the head rules.
    pub head_bits: Vec<u64>,
    /// Bit per site (head-passing ordinal): survived the tail rules.
    pub survivor_bits: Vec<u64>,
    /// Per-site emit data (`None` = tail-failed).
    pub survivors: Vec<Option<SiteEmit>>,
    /// The chunk's trailing seam record's emit (the seam site is appended
    /// to this chunk's sketch by the sink).
    pub seam: Option<SiteEmit>,
}

impl ChunkDecision {
    /// Whether record `r` passed the head rules.
    fn head_passed(&self, r: usize) -> bool {
        self.head_bits.get(r / 64).is_some_and(|word| word >> (r % 64) & 1 == 1)
    }

    /// Whether site `s` survived the tail rules.
    fn survivor(&self, s: usize) -> bool {
        self.survivor_bits.get(s / 64).is_some_and(|word| word >> (s % 64) & 1 == 1)
    }
}

/// Pass0 output of the sketch run.
pub struct SketchPassOutput {
    /// Per-chunk sketches in chunk order.
    pub sketches: Vec<ChunkSketch>,
    /// Genome-wide per-sample raw missing counts (all records).
    pub cohort_miss: Vec<u64>,
    /// Records observed (body + seams).
    pub n_records: u64,
    /// The pass0 funnel (head-rule hits; tail/mask counters still zero —
    /// the mid stage fills them).
    pub funnel: FunnelPartial,
}

/// One chunk's sketch plus where its trailing seam site landed.
pub struct ChunkSketch {
    pub id: u64,
    pub fragment: reed_filter::sketch::SketchFragment,
    /// Site index of the seam record within the fragment (the seam is the
    /// fragment's last record when it passed the head rules).
    pub seam_site: Option<usize>,
}

/// Worker scratch of the sketch pass.
struct SketchScratch {
    scanner: GenotypeScanner,
    head: HeadScratch,
    mask_values: Vec<Option<f64>>,
    p0a: Vec<u64>,
    p0m: Vec<u64>,
    p1a: Vec<u64>,
    p1m: Vec<u64>,
    rule_bits: Vec<Vec<u64>>,
    miss: Vec<u64>,
    counters: PolicyCounters,
    funnel: FunnelPartial,
}

/// The pass0 job: head rules → per-sample GT planes + per-rule mask
/// bitmaps (head survivors) + per-sample raw missing counts (all records,
/// the cohort's own semantics).
struct SketchJob {
    plan: Arc<FilterPlan>,
    n_samples: usize,
    w: usize,
}

/// The sketch pass's chunk partial. `pre` carries the seam record's
/// fragment — it belongs to the *previous* chunk's stream.
#[derive(Default)]
pub struct SketchPartial {
    pre: Option<reed_filter::sketch::SketchFragment>,
    fragment: reed_filter::sketch::SketchFragment,
    cohort_miss: Vec<u64>,
    funnel: FunnelPartial,
}

impl ChunkJob for SketchJob {
    type Partial = SketchPartial;
    type Scratch = Option<SketchScratch>;

    fn run_chunk(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        let mut fragment = reed_filter::sketch::SketchFragment::new(self.w, self.plan.mask.len());
        scratch.funnel.reset(
            self.plan.head.len(),
            self.plan.tail.len(),
            self.plan.mask.len(),
        );
        scratch.miss.clear();
        scratch.miss.resize(self.n_samples, 0);
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            self.sketch_line(line, scratch, &mut fragment)?;
        }
        Ok(SketchPartial {
            pre: None,
            fragment,
            cohort_miss: scratch.miss.clone(),
            funnel: std::mem::take(&mut scratch.funnel),
        })
    }

    fn run_seam(
        &self,
        _id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        let mut pre = reed_filter::sketch::SketchFragment::new(self.w, self.plan.mask.len());
        scratch.funnel.reset(
            self.plan.head.len(),
            self.plan.tail.len(),
            self.plan.mask.len(),
        );
        scratch.miss.clear();
        scratch.miss.resize(self.n_samples, 0);
        self.sketch_line(line, scratch, &mut pre)?;
        Ok(SketchPartial {
            pre: Some(pre),
            fragment: reed_filter::sketch::SketchFragment::new(self.w, self.plan.mask.len()),
            cohort_miss: scratch.miss.clone(),
            funnel: std::mem::take(&mut scratch.funnel),
        })
    }
}

impl SketchJob {
    fn build_scratch(&self) -> SketchScratch {
        SketchScratch {
            scanner: GenotypeScanner::new(
                self.n_samples,
                &self.plan.int_fields.iter().map(Vec::as_slice).collect::<Vec<_>>(),
                &self.plan.float_fields.iter().map(Vec::as_slice).collect::<Vec<_>>(),
                self.plan.policies,
            ),
            head: HeadScratch::new(),
            mask_values: vec![None; self.plan.mask.len()],
            p0a: vec![0; self.w],
            p0m: vec![0; self.w],
            p1a: vec![0; self.w],
            p1m: vec![0; self.w],
            rule_bits: vec![vec![0; self.w]; self.plan.mask.len()],
            miss: vec![0; self.n_samples],
            counters: PolicyCounters::default(),
            funnel: FunnelPartial {
                head_hits: vec![0; self.plan.head.len()],
                tail_hits: vec![0; self.plan.tail.len()],
                mask_hits: vec![0; self.plan.mask.len()],
                ..FunnelPartial::default()
            },
        }
    }

    /// Sketches one record: head rules first (the funnel counts every rule
    /// passed); head survivors contribute a site (planes + per-rule mask
    /// bitmaps over all samples — the mid stage intersects with `keep`).
    /// Every record feeds the cohort's raw missing counts.
    fn sketch_line(
        &self,
        line: &[u8],
        scratch: &mut SketchScratch,
        fragment: &mut reed_filter::sketch::SketchFragment,
    ) -> Result<(), EngineError> {
        let plan = &self.plan;
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::Job("record line has fewer than 9 columns".into()));
        };
        let site = SiteView::parse(&lr).map_err(|e| EngineError::Job(e.to_string()))?;
        let pos1 = site.pos1;
        let passed = eval_head(&lr, pos1, &plan.head, &mut scratch.funnel.head_hits, &mut scratch.head);
        if passed {
            for v in [&mut scratch.p0a, &mut scratch.p0m, &mut scratch.p1a, &mut scratch.p1m] {
                v.iter_mut().for_each(|word| *word = 0);
            }
            for v in &mut scratch.rule_bits {
                v.iter_mut().for_each(|word| *word = 0);
            }
        }
        // The scan runs for EVERY record: the cohort's raw missing counts
        // cover all records, while the plane/mask bitmaps are only filled
        // for head survivors (§11.4 pass0).
        let site_open = passed;
        let mask_rules = &plan.mask;
        scratch.counters = PolicyCounters::default();
        let out = scratch.scanner.scan(&lr, &site, &mut scratch.counters)?;
        scratch.funnel.policy.add_decode(&scratch.counters);
        for (i, &code) in out.codes.iter().enumerate() {
            if code.is_fully_missing()
                && let Some(slot) = scratch.miss.get_mut(i)
            {
                *slot += 1;
            }
            if !site_open {
                continue;
            }
            let (word, bit) = (i / 64, 1_u64 << (i % 64));
            if code.slot_alt(0) {
                scratch.p0a[word] |= bit;
            }
            if code.slot_missing(0) {
                scratch.p0m[word] |= bit;
            }
            if code.slot_alt(1) {
                scratch.p1a[word] |= bit;
            }
            if code.slot_missing(1) {
                scratch.p1m[word] |= bit;
            }
            for (rule, v) in mask_rules.iter().zip(scratch.mask_values.iter_mut()) {
                *v = match rule.float_column {
                    Some(j) => out
                        .float_columns
                        .get(j)
                        .and_then(|column| column.get(i))
                        .copied()
                        .filter(|x| !x.is_nan())
                        .map(f64::from),
                    None => out
                        .int_columns
                        .get(rule.int_column)
                        .and_then(|column| column.get(i))
                        .copied()
                        .filter(|&x| x != INT_MISSING)
                        .map(f64::from),
                };
            }
            if let Some(rule_index) = eval_mask(&scratch.mask_values, mask_rules) {
                scratch.rule_bits[rule_index][word] |= bit;
            }
        }
        if passed {
            fragment.push_site(
                &[&scratch.p0a, &scratch.p0m, &scratch.p1a, &scratch.p1m],
                &scratch.rule_bits.iter().map(Vec::as_slice).collect::<Vec<_>>(),
            );
        }
        fragment.push_record(passed);
        scratch.funnel.n_input += 1;
        Ok(())
    }
}

/// Collects chunk sketches in order; the seam fragment closes the
/// previous chunk (§11.4: the seam record belongs to the previous
/// chunk's stream).
struct SketchSink {
    sketches: Vec<ChunkSketch>,
    cohort: Vec<u64>,
    funnel: FunnelPartial,
    open: Option<(u64, reed_filter::sketch::SketchFragment)>,
    out: Arc<Mutex<Option<SketchPassOutput>>>,
}

impl OrderedSink<SketchPartial> for SketchSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut SketchPartial) -> Result<(), EngineError> {
        self.funnel.merge(&partial.funnel);
        if self.cohort.len() < partial.cohort_miss.len() {
            self.cohort.resize(partial.cohort_miss.len(), 0);
        }
        for (dst, src) in self.cohort.iter_mut().zip(&partial.cohort_miss) {
            *dst += src;
        }
        if let Some(pre) = partial.pre.take() {
            // The seam record closes the previously opened chunk.
            if let Some((cid, mut frag)) = self.open.take() {
                let seam_site = (pre.n_sites > 0).then_some(frag.n_sites);
                frag.append(&pre);
                self.sketches.push(ChunkSketch { id: cid, fragment: frag, seam_site });
            }
        }
        // Open the new chunk with its body fragment (a stale open without a
        // seam would mean a missing seam partial — finalize defensively).
        if let Some((cid, frag)) = self.open.take() {
            self.sketches.push(ChunkSketch { id: cid, fragment: frag, seam_site: None });
        }
        self.open = Some((_id.0, std::mem::take(&mut partial.fragment)));
        Ok(())
    }

    fn finish(mut self) -> Result<u64, EngineError> {
        if let Some((cid, frag)) = self.open.take() {
            self.sketches.push(ChunkSketch { id: cid, fragment: frag, seam_site: None });
        }
        let n_records = self.funnel.n_input;
        *self.out.lock().expect("invariant: sketch slot not poisoned") = Some(SketchPassOutput {
            sketches: self.sketches,
            cohort_miss: self.cohort,
            n_records,
            funnel: self.funnel,
        });
        Ok(n_records)
    }
}

/// The sketch pass (engine run #1).
fn run_sketch_pass(
    input: Option<&Path>,
    n_samples: usize,
    plan: &FilterPlan,
    opts: &FilterOptions,
) -> Result<SketchPassOutput, FilterError> {
    let budget = ThreadBudget { threads: opts.threads.max(1), mem: reed_engine::budget::MemBudget::default_half_memory() };
    let (mut source, _facts) = crate::input::open_regions(input, &opts.regions).map_err(io_err)?;
    let job = SketchJob {
        plan: Arc::new(plan.clone()),
        n_samples,
        w: n_samples.div_ceil(64),
    };
    let out: Arc<Mutex<Option<SketchPassOutput>>> = Arc::new(Mutex::new(None));
    let sink = SketchSink {
        sketches: Vec::new(),
        cohort: Vec::new(),
        funnel: sized_funnel(plan),
        open: None,
        out: Arc::clone(&out),
    };
    let abort = reed_engine::abort::abort();
    run(source.as_mut(), &job, sink, &budget, &abort).map_err(ferr)?;
    Ok(out
        .lock()
        .expect("invariant: sketch slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the sketch output"))
}

/// A sized funnel (the sketch sink's counters must exist before merges).
fn sized_funnel(plan: &FilterPlan) -> FunnelPartial {
    FunnelPartial {
        head_hits: vec![0; plan.head.len()],
        tail_hits: vec![0; plan.tail.len()],
        mask_hits: vec![0; plan.mask.len()],
        ..FunnelPartial::default()
    }
}

/// The mid stage (§11.4): per site, counts from `keep ∧ planes`, tail
/// rules, per-rule mask attribution, and the survivor's emit data.
fn mid_stage(
    sketches: &[ChunkSketch],
    keep: &[bool],
    plan: &FilterPlan,
    n_samples_out: usize,
    funnel: &mut FunnelPartial,
) -> std::collections::HashMap<u64, Arc<ChunkDecision>> {
    // Clean-padded keep words (the count math anchors on this).
    let w = keep.len().div_ceil(64);
    let mut keep_words = vec![0_u64; w];
    for (i, &k) in keep.iter().enumerate() {
        if k {
            keep_words[i / 64] |= 1 << (i % 64);
        }
    }
    let mut map = std::collections::HashMap::new();
    for sketch in sketches {
        let frag = &sketch.fragment;
        let mut survivor_bits = vec![0_u64; frag.n_sites.div_ceil(64)];
        let mut survivors: Vec<Option<SiteEmit>> = Vec::with_capacity(frag.n_sites);
        for s in 0..frag.n_sites {
            let c = reed_filter::sketch::count_site_sketch(frag, s, &keep_words);
            funnel.masked_gt += c.masked_gt;
            for (r, hits) in
                reed_filter::sketch::rule_hits_of(frag, s, &keep_words).into_iter().enumerate()
            {
                if let Some(slot) = funnel.mask_hits.get_mut(r) {
                    *slot += hits;
                }
            }
            let counts = SiteCounts {
                ac: c.ac,
                an: c.an,
                ns: c.ns,
                n_missing: c.n_missing,
                hom_ref: c.hom_ref,
                het: c.het,
                hom_alt: c.hom_alt,
            };
            let tail_values = reed_filter::expr::SiteValues::from_counts(reed_filter::expr::GtCounts {
                ac: counts.ac,
                an: counts.an,
                n_missing: counts.n_missing,
                n_samples: n_samples_out as u64,
            });
            if eval_tail(&tail_values, &counts, n_samples_out, &plan.tail, &mut funnel.tail_hits) {
                survivor_bits[s / 64] |= 1 << (s % 64);
                let mut infos = Vec::new();
                fill_info_set(&counts, &plan.recompute, &mut infos);
                let mask_words: Vec<u64> = frag
                    .mask_of(s)
                    .iter()
                    .zip(&keep_words)
                    .map(|(m, k)| m & k)
                    .collect();
                survivors.push(Some(SiteEmit { mask_words, infos }));
                funnel.n_output += 1;
            } else {
                survivors.push(None);
            }
        }
        let seam = sketch
            .seam_site
            .and_then(|s| survivors.get(s).cloned())
            .flatten();
        map.insert(
            sketch.id,
            Arc::new(ChunkDecision {
                n_records: frag.n_records,
                head_bits: frag.head_bits.clone(),
                survivor_bits,
                survivors,
                seam,
            }),
        );
    }
    map
}

/// The pass1 emit job: re-inflate the chunk, skip head-failed and
/// tail-failed records with a newline scan, and edit+compress survivors
/// from the mid-stage decisions (§11.4).
struct EmitJob {
    plan: Arc<FilterPlan>,
    decisions: Arc<std::collections::HashMap<u64, Arc<ChunkDecision>>>,
    keep_mask: Option<SampleMask>,
    n_samples_in: usize,
    compressed: bool,
    contig_names: Arc<[Vec<u8>]>,
    level: u8,
}

struct EmitScratch {
    edit_buf: Vec<u8>,
    block_buf: Vec<u8>,
    records: Vec<(u32, u64, u32, u32)>,
    blocks: Vec<OutBlock>,
    raw: Vec<u8>,
    deflater: Option<Deflater>,
    contigs: ContigTable,
    mask_bitset: Bitset,
}

impl ChunkJob for EmitJob {
    type Partial = FilterPartial;
    type Scratch = Option<EmitScratch>;

    fn run_chunk(
        &self,
        chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let dec = self
            .decisions
            .get(&chunk.id.0)
            .ok_or_else(|| EngineError::Job(format!("no emit decision for chunk {}", chunk.id.0)))?;
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        let mut record = 0_usize;
        let mut site = 0_usize;
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            let passed = dec.head_passed(record);
            record += 1;
            if !passed {
                continue; // head-failed: newline scan only
            }
            let s = site;
            site += 1;
            if !dec.survivor(s) {
                continue; // tail-failed: newline scan only
            }
            let Some(emit) = dec.survivors.get(s).and_then(|e| e.as_ref()) else {
                continue;
            };
            self.emit_line(line, emit, scratch)?;
        }
        self.cut_block(scratch)?;
        Ok(FilterPartial {
            blocks: std::mem::take(&mut scratch.blocks),
            raw: std::mem::take(&mut scratch.raw),
            funnel: FunnelPartial::default(),
        })
    }

    fn run_seam(
        &self,
        id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        // The seam record belongs to the previous chunk's stream; its
        // decision is that chunk's `seam` field.
        let emitted = id
            .0
            .checked_sub(1)
            .and_then(|prev| self.decisions.get(&prev))
            .and_then(|dec| dec.seam.clone());
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        if let Some(emit) = emitted {
            self.emit_line(line, &emit, scratch)?;
        }
        self.cut_block(scratch)?;
        Ok(FilterPartial {
            blocks: std::mem::take(&mut scratch.blocks),
            raw: std::mem::take(&mut scratch.raw),
            funnel: FunnelPartial::default(),
        })
    }
}

impl EmitJob {
    fn build_scratch(&self) -> EmitScratch {
        EmitScratch {
            edit_buf: Vec::new(),
            block_buf: Vec::with_capacity(BGZF_BLOCK_SIZE + (1 << 14)),
            records: Vec::new(),
            blocks: Vec::new(),
            raw: Vec::new(),
            deflater: self.compressed.then(|| Deflater::new(i32::from(self.level)).expect("invariant: level")),
            contigs: ContigTable::new(self.contig_names.iter().cloned()),
            mask_bitset: Bitset::new(self.n_samples_in),
        }
    }

    /// Edits one survivor line from its mid-stage emit and appends the
    /// bytes to the output block (or the raw buffer).
    fn emit_line(
        &self,
        line: &[u8],
        emit: &SiteEmit,
        scratch: &mut EmitScratch,
    ) -> Result<(), EngineError> {
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::Job("survivor line has fewer than 9 columns".into()));
        };
        let site = SiteView::parse_declared(&lr, &mut scratch.contigs)
            .map_err(|e| EngineError::Job(e.to_string()))?;
        let pos1 = site.pos1;
        let stale: Vec<&[u8]> = self.plan.stale_tags.iter().map(Vec::as_slice).collect();
        // Refresh the reusable mask bitset from the emit's words.
        for (word, src) in scratch.mask_bitset.words_mut().iter_mut().zip(&emit.mask_words) {
            *word = *src;
        }
        scratch.mask_bitset.clear_padding();
        let edit_plan = EditPlan {
            keep_samples: self.keep_mask.as_ref(),
            mask_gt: Some(&scratch.mask_bitset),
            info_remove: &stale,
            info_set: &emit.infos,
        };
        scratch.edit_buf.clear();
        edit_vcf_line(&lr, &edit_plan, gt_subfield_index(lr.format_text()), &mut scratch.edit_buf)
            .map_err(|e| EngineError::Job(e.to_string()))?;
        if self.compressed {
            // Byte-boundary cutting, records may span blocks (as in
            // FilterJob::push_line).
            let tid = site
                .contig
                .expect("invariant: parse_declared resolved the contig");
            let ref_len = u32::try_from(site.reference.len()).unwrap_or(1).max(1);
            let uoffset = u32::try_from(scratch.block_buf.len())
                .expect("invariant: pending block below u32::MAX");
            scratch.records.push((tid, pos1, ref_len, uoffset));
            scratch.block_buf.extend_from_slice(&scratch.edit_buf);
            while scratch.block_buf.len() >= BGZF_BLOCK_SIZE {
                self.cut_block(scratch)?;
            }
        } else {
            scratch.raw.extend_from_slice(&scratch.edit_buf);
        }
        Ok(())
    }


    fn cut_block(&self, scratch: &mut EmitScratch) -> Result<(), EngineError> {
        if scratch.block_buf.is_empty() {
            return Ok(());
        }
        let Some(deflater) = &mut scratch.deflater else {
            return Ok(());
        };
        let cut = scratch.block_buf.len().min(BGZF_BLOCK_SIZE);
        let cut32 = u32::try_from(cut).expect("invariant: cut ≤ 65,280");
        let mut compressed = Vec::with_capacity(cut / 2 + 64);
        deflater
            .deflate_block(&scratch.block_buf[..cut], &mut compressed)
            .map_err(|e| EngineError::Job(format!("bgzf deflate: {e}")))?;
        let boundary = scratch
            .records
            .partition_point(|&(_, _, _, uoffset)| uoffset < cut32);
        let mut records = std::mem::take(&mut scratch.records);
        let mut continuation = records.split_off(boundary);
        for record in continuation.iter_mut() {
            record.3 -= cut32;
        }
        scratch.records = continuation;
        scratch.blocks.push(OutBlock {
            compressed,
            uncompressed_len: cut32,
            records,
        });
        scratch.block_buf.drain(..cut);
        Ok(())
    }
}

/// Estimates the record count from the first chunk's records/byte density
/// extrapolated over the file size (§11.4, spec-sanctioned without an
/// index).
fn estimate_records(input: &Path) -> Result<u64, FilterError> {
    let file_len = std::fs::metadata(input).map_err(FilterError::from)?.len();
    let (mut source, _facts) = crate::input::open(Some(input)).map_err(io_err)?;
    let Some(chunk) = source.next_chunk().map_err(io_err)? else {
        return Ok(0);
    };
    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
    let mut scratch = Vec::new();
    let view = reed_io::chunk::inflate_chunk(&chunk, &mut inflater, &mut scratch).map_err(io_err)?;
    let body = view.body_bytes(view.text);
    let n_lines = u64::try_from(body.iter().filter(|&&b| b == b'\n').count()).unwrap_or(0);
    let compressed = u64::try_from(chunk.compressed.len()).unwrap_or(1).max(1);
    if n_lines == 0 {
        return Ok(0);
    }
    Ok(file_len.saturating_mul(n_lines) / compressed)
}
