//! Execution engine: pass0 (genome-wide sample missing-rate reduction) then
//! pass1 (filter + surgical edits + emit).
//!
//! Evaluation order (bit-for-bit aligned with the gold pipeline's bcftools
//! branch, docs/filter-architecture.md §4): RegionRule -> biallelic SNP ->
//! FILTER==PASS -> QUAL -> INFO (raw values) -> [GT mask] -> site missing
//! (post-mask) -> MAF (recomputed post-mask) -> HWE (recomputed post-mask).
//!
//! Execution strategy (v0.2, docs/filter-architecture.md §7): a batched
//! pipeline. A reader stage owns the noodles reader via
//! [`reed_io::record::for_each_record`] and hands fixed-size record batches
//! to a bounded crossbeam dispatch channel; `opts.threads` compute workers
//! pull batches from it and run the per-record pipeline (decode -> head
//! rules -> mask -> counts -> tail rules -> apply for survivors); the
//! calling thread reorders finished batches by batch id and emits them to
//! the sink in strict input order. The reorder plus order-independent
//! integer summation of the funnel partials make the output byte-identical
//! to a sequential run regardless of the thread count. Memory stays
//! bounded: only a few batches are ever queued or in flight. The two passes
//! run sequentially (the sample reduction must finish before site
//! filtering); decompression is still serial in reed-io's record transport
//! and compression concurrency comes from the caller's multithreaded BGZF
//! sink. No index appears at this layer: RegionRule is a streaming
//! predicate. All record I/O goes through reed-io
//! ([`reed_io::record::for_each_record`], [`reed_io::record::record_header`]);
//! this crate owns policy only.
//!
//! Error handling is first-error-wins: every stage shares an error slot; a
//! failing stage records its error and every other stage observes the flag
//! at its next blocking point and unwinds (bounded sends poll the flag),
//! so the pipeline always terminates and reports the first failure.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::time::{Duration, Instant};

use crossbeam_channel::{SendTimeoutError, Sender};

use reed_core::estimators::site_counts::SiteCounts;
use reed_io::record::{self, GtMask, RecordView, Recompute};
use reed_io::sink::VariantSink;
use reed_io::noodles_vcf as vcf;
use reed_io::noodles_vcf::variant::RecordBuf;

use crate::plan;
use crate::spec::legacy::{CohortRule, FilterSpec};
use crate::FilterError;

/// Records per batch handed to the compute stage
/// (docs/filter-architecture.md §7).
const BATCH_SIZE: usize = 4096;

/// Batches allowed to sit in a dispatch queue at once; together with the
/// bounded outcome channel this bounds pipeline memory to a few batches.
const MAX_BATCHES_IN_FLIGHT: usize = 4;

/// Finished batches allowed to wait for the emit stage.
const MAX_OUTCOMES_QUEUED: usize = 8;

/// Poll interval at which a stage blocked on a full channel re-checks the
/// shared error flag (the cancel path of every blocking send).
const CANCEL_POLL: Duration = Duration::from_millis(20);

/// Execution options.
#[derive(Clone, Debug)]
pub struct ExecOptions {
    /// Compute-stage worker count (the batched pipeline's filter workers).
    /// The output compression pool is sized separately by the caller when
    /// opening the sink; input decompression is currently serial (the next
    /// performance lever).
    pub threads: usize,
}

impl Default for ExecOptions {
    fn default() -> Self {
        Self {
            threads: std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(4),
        }
    }
}

/// Report of one execution (filter funnel plus counts).
#[derive(Clone, Debug)]
pub struct RunReport {
    /// Input path.
    pub input: PathBuf,
    /// Total records read by pass1 (including records excluded by the
    /// region rule).
    pub n_input_records: u64,
    /// Funnel: every activated rule in evaluation order with the cumulative
    /// count of records that passed it.
    pub funnel: Vec<(String, u64)>,
    /// Records written to the sink (the funnel survivors).
    pub n_output_records: u64,
    /// Sample columns before the cohort rule.
    pub n_samples_in: usize,
    /// Sample columns after the cohort rule.
    pub n_samples_out: usize,
    /// Genotypes set to ./. by the mask rule: called-to-missing transitions
    /// only (samples already missing before the mask are not counted).
    pub n_masked_gt: u64,
    /// Total wall-clock time of the run.
    pub elapsed: std::time::Duration,
}

/// First-error-wins slot shared by every pipeline stage. A failing stage
/// records its error here; all stages poll [`ErrorSlot::is_failed`] at
/// every blocking point and unwind, which is what keeps the pipeline
/// deadlock-free while still reporting the first failure.
struct ErrorSlot {
    failed: AtomicBool,
    first: std::sync::Mutex<Option<FilterError>>,
}

impl ErrorSlot {
    fn new() -> Self {
        Self {
            failed: AtomicBool::new(false),
            first: std::sync::Mutex::new(None),
        }
    }

    /// Records `err` unless an earlier error is already stored.
    fn record(&self, err: FilterError) {
        self.failed.store(true, Ordering::SeqCst);
        let mut first = self.first.lock().unwrap_or_else(|e| e.into_inner());
        if first.is_none() {
            *first = Some(err);
        }
    }

    fn is_failed(&self) -> bool {
        self.failed.load(Ordering::SeqCst)
    }

    /// Takes the recorded error (the pipeline failed) or None (success).
    fn take(&self) -> Option<FilterError> {
        let mut first = self.first.lock().unwrap_or_else(|e| e.into_inner());
        first.take()
    }
}

/// Reader-stage error: a real failure, or the sentinel that stops the read
/// once another stage has failed ([`record::for_each_record`] propagates
/// callback errors, which is how the reader breaks out of the stream).
enum ReadStageError {
    Failed(FilterError),
    Cancelled,
}

impl From<reed_io::ReadError> for ReadStageError {
    fn from(e: reed_io::ReadError) -> Self {
        ReadStageError::Failed(FilterError::from(e))
    }
}

/// Error-aware bounded send: retries `msg` until the channel accepts it,
/// the pipeline is cancelled (another stage failed), or the receiving end
/// is gone. The payload is restored to `msg` on failure so the caller keeps
/// ownership. Returns false when the send was abandoned.
fn send_bounded<T>(tx: &Sender<T>, err: &ErrorSlot, msg: &mut Option<T>) -> bool {
    loop {
        if err.is_failed() {
            return false;
        }
        let Some(payload) = msg.take() else {
            return true; // nothing to send; treat as delivered
        };
        match tx.send_timeout(payload, CANCEL_POLL) {
            Ok(()) => return true,
            Err(SendTimeoutError::Timeout(again)) => *msg = Some(again),
            Err(SendTimeoutError::Disconnected(again)) => {
                *msg = Some(again);
                return false;
            }
        }
    }
}

/// Reader stage: streams the input into fixed-size batches and dispatches
/// `(batch_id, batch)` pairs, with batch ids assigned in input order. Sends
/// are error-aware (see [`send_bounded`]) so a failure in any other stage
/// unblocks this one. The final partial batch keeps its sequence position.
fn read_batches(
    input: &Path,
    err: &ErrorSlot,
    dispatch: &Sender<(u64, Vec<RecordBuf>)>,
    label: &str,
) -> Result<(), ReadStageError> {
    let t_reader = Instant::now();
    let mut cb_ns: u128 = 0;
    let mut batch: Vec<RecordBuf> = Vec::with_capacity(BATCH_SIZE);
    let mut batch_id: u64 = 0;

    let result = record::for_each_record(input, |buf| -> Result<(), ReadStageError> {
        let t_cb = Instant::now();
        if err.is_failed() {
            return Err(ReadStageError::Cancelled);
        }
        batch.push(buf);
        if batch.len() == BATCH_SIZE {
            let mut msg = Some((
                batch_id,
                std::mem::replace(&mut batch, Vec::with_capacity(BATCH_SIZE)),
            ));
            if !send_bounded(dispatch, err, &mut msg) {
                return Err(ReadStageError::Cancelled);
            }
            batch_id += 1;
        }
        cb_ns += t_cb.elapsed().as_nanos();
        Ok(())
    });

    tracing::debug!(
        label,
        reader_wall_s = t_reader.elapsed().as_secs_f64(),
        callback_s = cb_ns as f64 / 1e9,
        "reader stage"
    );

    match result {
        // The stream ended (or was cancelled after another stage failed);
        // hand over the trailing partial batch, if any and if still wanted.
        Ok(()) | Err(ReadStageError::Cancelled) => {
            if !batch.is_empty() && !err.is_failed() {
                let mut msg = Some((batch_id, batch));
                send_bounded(dispatch, err, &mut msg);
            }
            Ok(())
        }
        Err(ReadStageError::Failed(e)) => Err(ReadStageError::Failed(e)),
    }
}

/// Pass0 as a batched pipeline: per-batch partial per-sample missing counts,
/// computed in parallel and summed afterwards. Integer sums are
/// order-independent, so the reduction matches a sequential scan exactly.
fn reduce_missing(
    input: &Path,
    n_samples: usize,
    threads: usize,
    err: &Arc<ErrorSlot>,
) -> Result<(u64, Vec<u64>), FilterError> {
    let (dispatch_tx, dispatch_rx) =
        crossbeam_channel::bounded::<(u64, Vec<RecordBuf>)>(MAX_BATCHES_IN_FLIGHT);
    let (partial_tx, partial_rx) = crossbeam_channel::unbounded::<(u64, Vec<u64>)>();

    let mut n_records: u64 = 0;
    let mut total: Vec<u64> = vec![0; n_samples];

    std::thread::scope(|scope| {
        // Reader stage.
        scope.spawn({
            let dispatch_tx = dispatch_tx.clone();
            let err = Arc::clone(err);
            move || {
                if let Err(ReadStageError::Failed(e)) = read_batches(input, &err, &dispatch_tx, "pass0") {
                    err.record(e);
                }
            }
        });
        // Compute stage: count fully-missing genotypes per batch.
        for w in 0..threads {
            scope.spawn({
                let dispatch_rx = dispatch_rx.clone();
                let partial_tx = partial_tx.clone();
                let err = Arc::clone(err);
                move || {
                    let mut recv = Duration::ZERO;
                    let mut work = Duration::ZERO;
                    loop {
                        let t_recv = Instant::now();
                        let Ok((_batch_id, batch)) = dispatch_rx.recv() else { break };
                        recv += t_recv.elapsed();
                        let t_work = Instant::now();
                        if err.is_failed() {
                            break;
                        }
                        let mut miss = vec![0u64; n_samples];
                        let mut failed = false;
                        for buf in &batch {
                            if let Err(e) = record::count_missing(buf, &mut miss) {
                                err.record(e.into());
                                failed = true;
                                break;
                            }
                        }
                        if failed {
                            break;
                        }
                        let len = batch.len() as u64;
                        if partial_tx.send((len, miss)).is_err() {
                            break;
                        }
                        work += t_work.elapsed();
                    }
                    tracing::debug!(worker = w, ?recv, ?work, "pass0 worker");
                }
            });
        }
        // Reduce stage (this thread): sum the partials as they arrive.
        // Dropping the original channel ends here is what turns the reader's
        // exit into the workers' end-of-stream signal.
        drop(dispatch_tx);
        drop(dispatch_rx);
        drop(partial_tx);
        let t_reduce = Instant::now();
        let mut recv = Duration::ZERO;
        for (len, miss) in partial_rx.iter() {
            let t_r = Instant::now();
            n_records += len;
            for (t, m) in total.iter_mut().zip(miss) {
                *t += m;
            }
            recv += t_r.elapsed();
        }
        tracing::debug!(reduce_wall_s = t_reduce.elapsed().as_secs_f64(), ?recv, "pass0 reduce");
    });

    // The scope joined every stage; a recorded error means the pass failed.
    match err.take() {
        Some(e) => Err(e),
        None => Ok((n_records, total)),
    }
}

/// One activated head rule: its funnel name and its predicate.
type HeadRule<'a> = (&'static str, Box<dyn Fn(&RecordView) -> bool + Send + Sync + 'a>);

/// Shared read-only context of the pass1 compute stage (borrowed; `Sync`
/// via the `Send + Sync` predicate bounds, so every worker can hold a
/// `&ScanCtx`).
struct ScanCtx<'a> {
    rules: &'a [HeadRule<'a>],
    genotypes: &'a Option<GtMask>,
    recompute: Recompute,
    keep: &'a [bool],
    n_samples_in: usize,
    n_samples_out: usize,
    site_missing: Option<f64>,
    min_maf: Option<f64>,
    min_hwe: Option<f64>,
}

/// Result of computing one batch in pass1: the surviving records plus this
/// batch's partials of every funnel counter. The emit stage sums the
/// partials in batch order, so the totals match a sequential scan exactly.
struct BatchOutcome {
    survivors: Vec<RecordBuf>,
    /// Cumulative pass counts of the head rules, aligned with the rule list.
    head_hits: Vec<u64>,
    site_missing_pass: u64,
    maf_pass: u64,
    hwe_pass: u64,
    /// Called-to-missing GT transitions the mask performed in this batch.
    masked_gt: u64,
    /// Records read for this batch (survivors are a subset of them).
    n_input: u64,
}

/// Computes one batch: the exact per-record pipeline of the sequential scan
/// (decode -> head rules -> mask -> counts -> tail rules -> apply for
/// survivors), collected into a [`BatchOutcome`].
fn process_batch(batch: Vec<RecordBuf>, ctx: &ScanCtx<'_>) -> Result<BatchOutcome, FilterError> {
    let mut outcome = BatchOutcome {
        survivors: Vec::new(),
        head_hits: vec![0; ctx.rules.len()],
        site_missing_pass: 0,
        maf_pass: 0,
        hwe_pass: 0,
        masked_gt: 0,
        n_input: batch.len() as u64,
    };

    for mut buf in batch {
        let mut site = record::decode(&buf, ctx.n_samples_in)?;
        // Head rules: increment every rule the record passes; the first
        // failure ends evaluation, so later head counters never see this
        // record.
        let mut passed_heads = true;
        for (hit, (_, ok)) in outcome.head_hits.iter_mut().zip(ctx.rules) {
            if !ok(&site) {
                passed_heads = false;
                break;
            }
            *hit += 1;
        }
        if !passed_heads {
            continue;
        }
        // The mask flags kept samples whose DP/GQ fail the rule, including
        // samples that were already missing; the report counts only true
        // called-to-missing transitions, so capture the called state
        // pre-mask.
        let called: Vec<bool> = (0..site.n_samples)
            .map(|i| site.gt[2 * i] >= 0 || site.gt[2 * i + 1] >= 0)
            .collect();
        site.mask(ctx.keep, ctx.genotypes);
        let counts: SiteCounts = site.counts(ctx.keep);
        // Tail rules run in order over the post-mask counts. Each counter
        // increments only for records that reach the rule and pass it; the
        // first tail failure stops evaluation, so later tail counters never
        // see this record.
        let mut ok = true;
        if let Some(m) = ctx.site_missing {
            if ctx.n_samples_out > 0 && counts.n_missing as f64 / ctx.n_samples_out as f64 > m {
                ok = false;
            } else {
                outcome.site_missing_pass += 1;
            }
        }
        if ok && let Some(m) = ctx.min_maf {
            if counts.maf() < m {
                ok = false;
            } else {
                outcome.maf_pass += 1;
            }
        }
        let mut hwe_p = None;
        if ok && let Some(threshold) = ctx.min_hwe {
            let p = counts.hwe_p();
            hwe_p = Some(p);
            if p < threshold {
                ok = false;
            } else {
                outcome.hwe_pass += 1;
            }
        }
        if !ok {
            continue;
        }
        record::apply(
            &mut buf,
            &site,
            &counts,
            ctx.keep,
            ctx.recompute,
            hwe_p,
            record::STALE_INFO,
        )?;
        outcome.masked_gt += site
            .masked
            .iter()
            .zip(called)
            .filter(|(masked, was_called)| **masked && *was_called)
            .count() as u64;
        outcome.survivors.push(buf);
    }
    Ok(outcome)
}

/// Run the filter. The caller opens the output sink (reed-io's
/// [`VariantSink`]); the library assumes nothing about where the records go.
///
/// `opts.threads` sizes the pass1 compute pool; the output compression pool
/// is applied by the caller when opening the sink.
pub fn execute(
    input: &Path,
    spec: &FilterSpec,
    emitter: &mut dyn VariantSink,
    opts: &ExecOptions,
) -> Result<RunReport, FilterError> {
    let t0 = Instant::now();
    let _plan = plan::legacy::compile(spec);

    let header = record::record_header(input)?;
    let n_samples_in = header.sample_names().len();
    let compute_threads = opts.threads.max(1);
    let err = Arc::new(ErrorSlot::new());

    // Stage timers feed the tracing events below.
    let t_header = t0.elapsed();
    let t_pass0 = Instant::now();

    // ---- pass0 (reduce): genome-wide sample missing-rate reduction (raw
    // counts, before masking), computed as a batched parallel reduction ----
    let (keep, _n_pass0_records) = match &spec.cohort {
        None => (Vec::new(), 0),
        Some(CohortRule::MaxMissing(threshold)) => {
            let (n_records, miss) = reduce_missing(input, n_samples_in, compute_threads, &err)?;
            let keep: Vec<bool> = if n_records == 0 {
                vec![true; n_samples_in]
            } else {
                miss.iter()
                    .map(|&m| (m as f64) / (n_records as f64) <= *threshold)
                    .collect()
            };
            (keep, n_records)
        }
    };
    tracing::info!(?t_header, pass0 = ?t_pass0.elapsed(), "filter phases");
    let t_pass1 = Instant::now();
    let n_samples_out = keep.iter().filter(|&&k| k).count();
    let all_kept = keep.is_empty() || keep.iter().all(|&k| k);

    // ---- pass1 (scan): filter + surgical edits + emit ----
    let out_header = if all_kept {
        header.clone()
    } else {
        record::filter_header(&header, &keep)
    };
    let out_header = ensure_info_defs(out_header, spec);
    emitter.write_header(&out_header)?;

    // Activated head rules in evaluation order; disabled rules stay out of
    // the funnel. The predicates are shared by every compute worker, hence
    // the Send + Sync bounds (all captures are Sync spec data or copies).
    let mut rules: Vec<HeadRule<'_>> = Vec::new();
    if !spec.regions.is_empty() {
        rules.push((
            "regions",
            Box::new(|s: &RecordView| spec.regions.contains(&s.contig, s.pos)),
        ));
    }
    if spec.sites.biallelic_snps {
        rules.push(("biallelic_snp", Box::new(|s: &RecordView| s.is_biallelic_snp())));
    }
    if spec.sites.require_pass {
        rules.push(("pass", Box::new(|s: &RecordView| s.is_pass)));
    }
    if let Some(q) = spec.sites.min_qual {
        rules.push(("qual", Box::new(move |s: &RecordView| {
            s.qual.is_some_and(|v| v >= q) // a missing QUAL fails (matches bcftools)
        })));
    }
    if spec.sites.info.any() {
        let info = spec.sites.info;
        rules.push(("info", Box::new(move |s: &RecordView| {
            // a sub-rule whose tag is missing passes (bcftools
            // missing-comparison semantics)
            info.min_qd.is_none_or(|m| s.info_qd.is_none_or(|v| v >= m))
                && info.max_fs.is_none_or(|m| s.info_fs.is_none_or(|v| v <= m))
                && info.min_mq.is_none_or(|m| s.info_mq.is_none_or(|v| v >= m))
                && info.max_sor.is_none_or(|m| s.info_sor.is_none_or(|v| v <= m))
        })));
    }

    let ctx = ScanCtx {
        rules: &rules,
        genotypes: &spec.genotypes,
        recompute: spec.recompute,
        keep: &keep,
        n_samples_in,
        n_samples_out,
        site_missing: spec.sites.max_missing,
        min_maf: spec.sites.min_maf,
        min_hwe: spec.sites.min_hwe_p,
    };

    let mut rule_hits: Vec<u64> = vec![0; rules.len()];
    let mut n_site_missing_pass: u64 = 0;
    let mut n_maf_pass: u64 = 0;
    let mut n_hwe_pass: u64 = 0;
    let mut n_input: u64 = 0;
    let mut n_output: u64 = 0;
    let mut n_masked_gt: u64 = 0;

    // ---- pass1 pipeline: read (its own thread) -> compute (worker pool) ->
    // emit (this thread, strictly in batch order) ----
    let (dispatch_tx, dispatch_rx) =
        crossbeam_channel::bounded::<(u64, Vec<RecordBuf>)>(MAX_BATCHES_IN_FLIGHT);
    let (outcome_tx, outcome_rx) =
        crossbeam_channel::bounded::<(u64, BatchOutcome)>(MAX_OUTCOMES_QUEUED);
    let ctx_ref = &ctx;

    std::thread::scope(|scope| {
        // Reader stage.
        scope.spawn({
            let dispatch_tx = dispatch_tx.clone();
            let err = Arc::clone(&err);
            move || {
                if let Err(ReadStageError::Failed(e)) = read_batches(input, &err, &dispatch_tx, "pass1") {
                    err.record(e);
                }
            }
        });
        // Compute stage: workers pull batches from the dispatch channel (the
        // work queue) and run the per-record pipeline; batch-level dispatch
        // balances the load by record count.
        for w in 0..compute_threads {
            let dispatch_rx = dispatch_rx.clone();
            let outcome_tx = outcome_tx.clone();
            let err = Arc::clone(&err);
            scope.spawn(move || {
                let mut recv = Duration::ZERO;
                let mut work = Duration::ZERO;
                loop {
                    let t_recv = Instant::now();
                    let Ok((batch_id, batch)) = dispatch_rx.recv() else { break };
                    recv += t_recv.elapsed();
                    let t_work = Instant::now();
                    if err.is_failed() {
                        break;
                    }
                    match process_batch(batch, ctx_ref) {
                        Ok(outcome) => {
                            let mut msg = Some((batch_id, outcome));
                            if !send_bounded(&outcome_tx, &err, &mut msg) {
                                break;
                            }
                        }
                        Err(e) => {
                            err.record(e);
                            break;
                        }
                    }
                    work += t_work.elapsed();
                }
                tracing::debug!(worker = w, ?recv, ?work, "pass1 worker");
            });
        }
        // Emit stage (this thread): reorder by batch id and hand the
        // survivors to the sink strictly in input order, summing the funnel
        // partials in batch order. Dropping the original channel ends here
        // is what turns the reader's exit into the workers' end-of-stream
        // signal (and the workers' exit into this loop's end-of-stream).
        drop(dispatch_tx);
        drop(dispatch_rx);
        drop(outcome_tx);
        let mut pending: BTreeMap<u64, BatchOutcome> = BTreeMap::new();
        let mut next_batch_id: u64 = 0;
        let t_emit = Instant::now();
        let mut emit_write = Duration::ZERO;
        let mut emit_drain = Duration::ZERO;
        while let Ok((batch_id, outcome)) = outcome_rx.recv() {
            let t_recv = Instant::now();
            pending.insert(batch_id, outcome);
            while let Some(outcome) = pending.remove(&next_batch_id) {
                n_input += outcome.n_input;
                for (total, hits) in rule_hits.iter_mut().zip(&outcome.head_hits) {
                    *total += hits;
                }
                n_site_missing_pass += outcome.site_missing_pass;
                n_maf_pass += outcome.maf_pass;
                n_hwe_pass += outcome.hwe_pass;
                n_masked_gt += outcome.masked_gt;
                for buf in &outcome.survivors {
                    // An output failure is this stage's error; recording it
                    // unblocks the reader and the workers.
                    let t_write = Instant::now();
                    let write_result = emitter.write_record(&out_header, buf);
                    emit_write += t_write.elapsed();
                    if let Err(e) = write_result {
                        err.record(FilterError::Io(e));
                        break;
                    }
                    n_output += 1;
                }
                next_batch_id += 1;
                if err.is_failed() {
                    break;
                }
            }
            emit_drain += t_recv.elapsed();
            if err.is_failed() {
                break;
            }
        }
        let emit_wall = t_emit.elapsed();
        tracing::debug!(
            ?emit_wall,
            ?emit_drain,
            write = ?emit_write,
            "pass1 emit"
        );
    });

    // The scope joined every stage; the first recorded error wins.
    if let Some(e) = err.take() {
        return Err(e);
    }
    tracing::info!(pass1 = ?t_pass1.elapsed(), "filter phases");

    // Funnel: activated rules in evaluation order, each with the cumulative
    // count of records that passed it. Head counters live in rule_hits; the
    // three tail rules accumulated their own counters in the pipeline.
    let mut funnel: Vec<(String, u64)> = Vec::new();
    for (i, (name, _)) in rules.iter().enumerate() {
        funnel.push(((*name).to_string(), rule_hits[i]));
    }
    if spec.sites.max_missing.is_some() {
        funnel.push(("site_missing".to_string(), n_site_missing_pass));
    }
    if spec.sites.min_maf.is_some() {
        funnel.push(("maf".to_string(), n_maf_pass));
    }
    if spec.sites.min_hwe_p.is_some() {
        funnel.push(("hwe".to_string(), n_hwe_pass));
    }

    emitter.finish()?;

    Ok(RunReport {
        input: input.to_path_buf(),
        n_input_records: n_input,
        funnel,
        n_output_records: n_output,
        n_samples_in,
        n_samples_out,
        n_masked_gt,
        elapsed: t0.elapsed(),
    })
}

/// Add header definitions for the recomputed tags the input lacks
/// (NS/MAF/HWE and friends); existing definitions stay untouched.
fn ensure_info_defs(mut header: vcf::Header, spec: &FilterSpec) -> vcf::Header {
    use reed_io::noodles_vcf::header::record::value::map::info::{Number, Type};
    type InfoMap = reed_io::noodles_vcf::header::record::value::Map<
        reed_io::noodles_vcf::header::record::value::map::Info,
    >;

    let add = |header: &mut vcf::Header, id: &str, number: Number, ty: Type, desc: &str| {
        let exists = header.infos().iter().any(|(k, _)| k.as_str() == id);
        if !exists {
            let key = id.to_string();
            header.infos_mut().insert(key, InfoMap::new(number, ty, desc));
        }
    };
    let r = spec.recompute;
    if r.ns {
        add(&mut header, "NS", Number::Count(1), Type::Integer, "Number of samples with data");
    }
    if r.maf {
        add(&mut header, "MAF", Number::Count(1), Type::Float,
            "Frequency of the second most common allele");
    }
    if r.hwe {
        add(&mut header, "HWE", Number::Count(1), Type::Float,
            "HWE test (PMID:15789306); 1=good, 0=bad");
    }
    if r.an {
        add(&mut header, "AN", Number::Count(1), Type::Integer,
            "Total number of alleles in called genotypes");
    }
    if r.ac {
        add(&mut header, "AC", Number::Count(1), Type::Integer,
            "Allele count in genotypes");
    }
    header
}

#[cfg(test)]
mod tests {
    //! End-to-end determinism gate for the batched pipeline: the same input
    //! filtered at several thread counts must produce byte-identical output
    //! and identical reports. The input is generated as plain VCF text (the
    //! record transport dispatches `.vcf` without any compression), with
    //! enough records to span multiple batches and rule hits at every stage
    //! of the funnel (mask, site missing, MAF, INFO, QUAL, PASS, biallelic).

    use super::*;
    use crate::spec::legacy::CohortRule;
    use reed_io::sink;

    const HEADER_TEXT: &str = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=100000000>
##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">
##INFO=<ID=FS,Number=1,Type=Float,Description=\"Fisher's exact test to filter strand-bias\">
##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">
##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\ts3
";

    /// Generates a plain-text VCF exercising every funnel rule. s3 is
    /// missing in 98% of the records, so the cohort rule drops it and the
    /// output columns are pruned.
    fn test_input_text(n_records: usize) -> String {
        let mut text = String::from(HEADER_TEXT);
        for i in 0..n_records {
            let pos = i + 1;
            let qual = if i % 50 == 25 { "10.0" } else { "50.0" }; // every 50th fails QUAL
            let filter = if i % 60 == 30 { "LowQual" } else { "PASS" }; // every 60th fails PASS
            let (ref_b, alt) = if i % 70 == 35 {
                ("AT", "AGG") // indel: fails the biallelic rule
            } else if i % 90 == 45 {
                ("A", "G,T") // multiallelic: fails the biallelic rule
            } else {
                ("A", "G")
            };
            let qd = if i % 80 == 40 { "2.0" } else { "10.0" }; // every 80th fails INFO
            let info = format!("QD={qd};FS=10.0;MQ=60.0;SOR=1.0");

            let s3_missing = i % 50 != 0;
            let gts: [&str; 3] = match i % 4 {
                0 => ["0/1", "0/1", "0/0"],
                1 => ["1/1", "0/1", "0/0"],
                2 => ["0/0", "0/1", "0/1"],
                _ => ["0/0", "0/0", "0/0"], // MAF 0 unless s3 carries the alt
            };
            let gts = if i % 40 == 20 { ["0/0", "0/0", "0/0"] } else { gts };
            let mut dps = [30u32; 3];
            if i % 10 == 5 {
                dps[1] = 3; // s1 below dp_min: masked
            }
            if i % 17 == 8 {
                dps[0] = 150; // s0 above dp_max: masked
            }
            let mut gqs = [99u32; 3];
            if i % 13 == 6 {
                gqs[2] = 5; // s2 below gq_min: masked
            }

            let sample = |gt: &str, dp: u32, gq: u32| format!("{gt}:{dp}:{gq}");
            let s3 = if s3_missing {
                "./.:0:0".to_string()
            } else {
                sample("1/1", 30, 99)
            };
            text.push_str(&format!(
                "Chr01\t{pos}\t.\t{ref_b}\t{alt}\t{qual}\t{filter}\t{info}\tGT:DP:GQ\t{}\t{}\t{}\t{s3}\n",
                sample(gts[0], dps[0], gqs[0]),
                sample(gts[1], dps[1], gqs[1]),
                sample(gts[2], dps[2], gqs[2]),
            ));
        }
        text
    }

    fn test_spec() -> FilterSpec {
        // s3 is missing in 98% of records; a 0.9 threshold drops it. With
        // only three kept samples, the default 0.10 site-missing cap would
        // reject every masked site (one missing sample of three is already
        // 0.33); 0.5 lets masked records survive so the masked transition
        // counter is exercised.
        FilterSpec {
            cohort: Some(CohortRule::MaxMissing(0.9)),
            sites: crate::spec::legacy::SiteRules {
                max_missing: Some(0.5),
                ..FilterSpec::default().sites
            },
            ..FilterSpec::default()
        }
    }

    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!(
            "reed-filter-exec-{name}-{}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("create scratch dir");
        dir
    }

    /// The batched pipeline must be byte-identical regardless of the thread
    /// count (including the threads=0 clamp), with identical funnels.
    #[test]
    fn pipeline_is_deterministic_across_thread_counts() {
        const N_RECORDS: usize = 10_000; // ~3 batches of 4096
        let dir = scratch("determinism");
        let input = dir.join("in.vcf");
        std::fs::write(&input, test_input_text(N_RECORDS)).expect("write test input");

        let mut reports = Vec::new();
        let mut out_paths = Vec::new();
        for (name, threads) in [("t0", 0usize), ("t1", 1), ("t3", 3), ("t7", 7)] {
            let out = dir.join(format!("out-{name}.vcf"));
            let mut out_sink = sink::open(Some(&out), threads).expect("open output sink");
            let report = execute(&input, &test_spec(), out_sink.as_mut(), &ExecOptions { threads })
                .unwrap_or_else(|e| panic!("execute {name}: {e}"));
            out_sink.finish().expect("finish output sink");
            reports.push(report);
            out_paths.push(out);
        }

        // Byte-identical output files.
        let reference = std::fs::read(&out_paths[0]).expect("read output");
        assert!(!reference.is_empty(), "output must not be empty");
        for out in &out_paths[1..] {
            assert_eq!(
                reference,
                std::fs::read(out).expect("read output"),
                "output differs for {}",
                out.display()
            );
        }

        // Identical reports (funnel, counts, masked transitions).
        let first = &reports[0];
        for r in &reports[1..] {
            assert_eq!(first.funnel, r.funnel);
            assert_eq!(first.n_input_records, r.n_input_records);
            assert_eq!(first.n_output_records, r.n_output_records);
            assert_eq!(first.n_masked_gt, r.n_masked_gt);
            assert_eq!(first.n_samples_in, r.n_samples_in);
            assert_eq!(first.n_samples_out, r.n_samples_out);
        }

        // The funnel actually fired at several stages and filtered records.
        assert_eq!(first.n_samples_out, 3, "cohort rule must drop s3");
        assert_eq!(first.n_input_records, N_RECORDS as u64);
        assert!(
            first.n_output_records > 0 && first.n_output_records < first.n_input_records,
            "the filter must keep and drop records: {} of {}",
            first.n_output_records,
            first.n_input_records
        );
        let hits: Vec<u64> = first.funnel.iter().map(|(_, n)| *n).collect();
        assert!(hits.windows(2).all(|w| w[0] >= w[1]), "funnel must be non-increasing");
        assert!(first.n_masked_gt > 0, "the mask rule must have fired");

        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The activated-rule list of the funnel is assembled deterministically:
    /// every run reports the same rule names in evaluation order.
    #[test]
    fn funnel_rule_list_is_stable() {
        let dir = scratch("funnel");
        let input = dir.join("in.vcf");
        std::fs::write(&input, test_input_text(100)).expect("write test input");

        let out = dir.join("out.vcf");
        let mut out_sink = sink::open(Some(&out), 2).expect("open output sink");
        let report = execute(&input, &test_spec(), out_sink.as_mut(), &ExecOptions { threads: 2 })
            .expect("execute");
        out_sink.finish().expect("finish output sink");

        let names: Vec<&str> = report.funnel.iter().map(|(name, _)| name.as_str()).collect();
        assert_eq!(
            names,
            vec!["biallelic_snp", "pass", "qual", "info", "site_missing", "maf"],
            "funnel order"
        );

        let _ = std::fs::remove_dir_all(&dir);
    }

    /// TEMP performance probe (remove after diagnosis): isolates why the
    /// batched reader stage parses several times slower than the v0.1
    /// sequential scan. V1 = the v0.1 pass0 shape (parse + count inline on
    /// the calling thread); V2 = the same work on a spawned thread; V3 =
    /// the batched reduce pipeline. Run with:
    ///   REED_PROBE_INPUT=<vcf.gz> cargo test --release -p reed-filter \
    ///     -- --ignored --nocapture probe_reader_variants
    #[test]
    #[ignore = "performance probe"]
    fn probe_reader_variants() {
        let Ok(input) = std::env::var("REED_PROBE_INPUT") else {
            eprintln!("[probe] set REED_PROBE_INPUT to a .vcf.gz; skipping");
            return;
        };
        let input = PathBuf::from(input);
        let n_samples = record::record_header(&input).expect("header").sample_names().len();

        let count = |input: &Path| -> (u64, Vec<u64>) {
            let mut n = 0u64;
            let mut miss = vec![0u64; n_samples];
            record::for_each_record(input, |buf| {
                n += 1;
                record::count_missing(&buf, &mut miss)
            })
            .expect("read");
            (n, miss)
        };

        let t = Instant::now();
        let (n1, m1) = count(&input);
        eprintln!("[probe] V1 sequential main-thread: {:.2?} n={n1}", t.elapsed());

        let t = Instant::now();
        std::thread::scope(|s| {
            s.spawn(|| count(&input));
        });
        eprintln!("[probe] V2 spawned thread: {:.2?}", t.elapsed());

        let err = Arc::new(ErrorSlot::new());
        let t = Instant::now();
        let (n3, _) = reduce_missing(&input, n_samples, 16, &err).expect("reduce");
        eprintln!("[probe] V3 batched pipeline: {:.2?} n={n3}", t.elapsed());

        let t = Instant::now();
        let (n4, m4) = count(&input);
        eprintln!("[probe] V1 again (warm): {:.2?} n={n4}", t.elapsed());
        assert_eq!(n1, n3);
        assert_eq!(m1, m4);
    }
}
