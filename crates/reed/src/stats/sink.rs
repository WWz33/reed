//! The ordered sink of the statistics run: merges chunk partials in input order,
//! emits window rows through the watermark, and closes contigs and the genome.

use crate::stats::request::{
    FStatsPartial, InvariantMode, LdPlan, Levels, SfsPartial, WindowLayoutOptions, WindowRow,
};
use crate::stats::run::{ContigAcc, GridMode, SegmentPartial, StatsPartial, decay_contig};
use reed_core::ChunkId;
use reed_core::coord::{Interval0, Pos0};
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::OrderedSink;
use reed_io::write::atomic::AtomicWriter;
use reed_stats::ld::LdDecayAcc;
use reed_stats::request::StatAcc;
use reed_stats::window::{BpBinning, ContigBins};
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

/// The LD analysis results the sink accumulated across contigs: the cohort
/// curve, one curve per named group, and the cohort pair rows.
pub(super) struct LdResult {
    pub(super) cohort: LdDecayAcc,
    pub(super) groups: Vec<LdDecayAcc>,
    pub(super) pairs: Vec<(Vec<u8>, u64, u64, f64)>,
}

/// Everything the sink hands back through the shared slot at `finish`.
pub(super) struct SinkOutput {
    pub(super) windows: Vec<WindowRow>,
    /// Genome-wide projected spectra.
    pub(super) sfs: Option<SfsPartial>,
    /// Genome-wide f-statistics.
    pub(super) fstats: Option<FStatsPartial>,
    /// LD results (the parallel matrix analysis).
    pub(super) ld: Option<LdResult>,
    /// The streaming `--ld-pairs` file, written and flushed but **not yet
    /// committed**. The driver commits it only after the input-completeness
    /// check passes; dropping it removes the temp and leaves the previous
    /// `.ld` untouched (doc_v27 V27-02).
    pub(super) pending_commit: Option<reed_io::write::atomic::AtomicFile>,
    /// Uncommitted streaming window files and their metadata. The driver
    /// publishes them only after input-integrity validation (V27-18).
    pub(super) window_session: Option<crate::stats::writers::session::WindowOutputSession>,
    /// The streamed pair table's final path, once it has been committed.
    /// Reported by the driver instead of being guessed by the CLI (V27-20).
    pub(super) pairs_path: Option<std::path::PathBuf>,
    pub(super) n_invariant_sites: u64,
    pub(super) n_no_alt_records: u64,
    pub(super) n_windows: u64,
    pub(super) n_records: u64,
    pub(super) n_multiallelic_skipped: u64,
}

/// The sink's per-contig accumulation state.
pub(super) struct CurrentContig {
    pub(super) name: Vec<u8>,
    pub(super) acc: ContigAcc,
    /// Contig-level totals (§8.7): every touched bin folded in.
    pub(super) totals: Vec<StatAcc>,
    pub(super) first_pos0: u64,
    pub(super) last_pos0: u64,
    /// Grid origin: windows start at 0 (1-based position 1), the
    /// the reference-tool grid (04-algorithms.md §4).
    pub(super) anchor: u64,
    /// Start of the next window to emit (the watermark frontier).
    pub(super) next_start: u64,
    /// End of the last emitted window (0-based).
    pub(super) span_end: u64,
}

/// Ordered consumer: merges same-contig segments, emits complete windows
/// through the watermark (§8.3), and closes contigs with a contig-level
/// row (§8.7). Revisited contigs are refused — unsorted input would
/// otherwise silently double-count windows.
pub(super) struct StatsSink {
    pub(super) grid: GridMode,
    pub(super) seed: Vec<StatAcc>,
    /// Request kinds, used to materialize scalar CLI output rows.
    pub(super) kinds: Arc<[reed_stats::request::StatKind]>,
    pub(super) levels: Levels,
    /// Whether completed window rows should drop large θ/Tajima vectors.
    /// Contig/genome rows always keep full accumulators for aggregation.
    pub(super) compact_window_accumulators: bool,
    pub(super) skip_empty: bool,
    /// Window-layout mode (structured rows + invariant policy).
    pub(super) layout: Option<WindowLayoutOptions>,
    /// Structured window rows; the flat table renders from these
    /// (`writers::windows::render_tsv`).
    pub(super) windows: Vec<WindowRow>,
    /// Reused full accumulator for one in-progress Bp window. It never
    /// enters `windows`; this avoids a fresh distribution allocation for
    /// every emitted window when the CLI compacts output (V27-18).
    pub(super) window_scratch: Vec<StatAcc>,
    /// Site-wide invariant sites seen.
    pub(super) n_invariant_sites: u64,
    /// Records with ALT = `.` (the all-sites check).
    pub(super) n_no_alt_records: u64,
    pub(super) contig_len: HashMap<Vec<u8>, u64>,
    pub(super) current: Option<CurrentContig>,
    pub(super) closed: HashSet<Vec<u8>>,
    pub(super) genome_totals: Vec<StatAcc>,
    /// Merged genome-wide projected spectra.
    pub(super) sfs: Option<SfsPartial>,
    /// Merged genome-wide f-statistics.
    pub(super) fstats: Option<FStatsPartial>,
    /// The LD plane layout (group-major); present only with an LD request.
    pub(super) ld_layout: Option<Arc<reed_stats::SampleLayout>>,
    /// The matrix being accumulated for the pending contig.
    pub(super) ld_matrix: Option<reed_stats::ContigMatrix>,
    /// One shared, abortable LD execution pool for the whole run.
    pub(super) ld_executor: Option<LdExecutor>,
    /// Decay results already joined mid-run by the in-flight cap — `finish`
    /// starts its merge from these (doc_v8 D8-2).
    pub(super) ld_merged: Option<LdResult>,
    /// Cap on unmerged contig matrices (not on OS threads).
    pub(super) max_ld_jobs: usize,
    /// The LD plan (the analysis needs the window and the bin edges).
    pub(super) ld: Option<LdPlan>,
    /// The streaming `--ld-pairs` table ([`LdPlan::pairs_path`]): the file is
    /// opened at construction and its header written even for an empty run,
    /// so a successful zero-pair run replaces the previous table instead of
    /// leaving stale rows behind (doc_v27 V27-05). Kept as a `(writer,
    /// atomic owner)` pair so publication still waits for every integrity
    /// check (doc_v27 V27-02).
    pub(super) pairs_file: Option<AtomicWriter>,
    pub(super) genome_span: u64,
    pub(super) genome_denom_bp: u64,
    pub(super) n_windows: u64,
    pub(super) n_records: u64,
    pub(super) n_multiallelic_skipped: u64,
    pub(super) out: Arc<std::sync::Mutex<Option<SinkOutput>>>,
    /// CLI streaming output; `None` retains the library's full rows.
    pub(super) output_session: Option<crate::stats::writers::session::WindowOutputSession>,
}

impl OrderedSink<StatsPartial> for StatsSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut StatsPartial) -> Result<(), EngineError> {
        self.n_records += partial.n_records;
        self.n_multiallelic_skipped += partial.n_multiallelic_skipped;
        self.n_invariant_sites += partial.n_invariant_sites;
        self.n_no_alt_records += partial.n_no_alt_records;
        if let Some(src) = &partial.sfs {
            match &mut self.sfs {
                Some(dst) => dst.merge(src),
                None => self.sfs = Some(src.clone()),
            }
        }
        if let Some(src) = &partial.fstats {
            match &mut self.fstats {
                Some(dst) => dst.merge(src),
                None => self.fstats = Some(src.clone()),
            }
        }
        if let Some(layout) = self.ld_layout.clone() {
            // The worker's runs append to the pending contig's matrix; a
            // contig change finalizes it (the analysis runs per contig, so
            // at most one contig's matrix is resident).
            if let Some(src) = partial.ld.take() {
                for run in src.runs {
                    let same = self
                        .ld_matrix
                        .as_ref()
                        .is_some_and(|m| m.contig() == run.contig.as_slice());
                    if !same {
                        self.finalize_ld()?;
                        self.ld_matrix = Some(reed_stats::ContigMatrix::new(
                            run.contig,
                            Arc::clone(&layout),
                        ));
                    }
                    if let Some(matrix) = self.ld_matrix.as_mut() {
                        matrix.push_run(run.flat, run.pos0);
                    }
                }
            }
        }
        for mut segment in std::mem::take(&mut partial.segments) {
            let mut current = match self.current.take() {
                Some(mut cur) if cur.name == segment.contig => {
                    // Cross-chunk seam (04-algorithms.md §6): a later chunk
                    // may not rewind below the last input site. `next_start`
                    // is a *window frontier*, not an input position: with
                    // `step > size`, sites legally fall between two windows,
                    // so comparing against it rejected sorted input at a
                    // chunk seam (doc_v27 V27-06). A same-POS duplicate at a
                    // seam is legal and lands in the still-open window
                    // (doc_v15 D15-19).
                    if segment.first_pos0 < cur.last_pos0 {
                        return Err(EngineError::from(
                            reed_core::CoreError::PositionNotAscending {
                                prev_pos0: cur.last_pos0,
                                pos0: segment.first_pos0,
                            },
                        ));
                    }
                    if matches!(cur.acc, ContigAcc::Bins(_)) {
                        segment.acc.fold(&mut cur.totals);
                    }
                    if !matches!(segment.acc, ContigAcc::Rows(_)) {
                        cur.acc.merge(&mut segment.acc);
                    }
                    cur.first_pos0 = cur.first_pos0.min(segment.first_pos0);
                    if !matches!(segment.acc, ContigAcc::Rows(_)) {
                        cur.last_pos0 = cur.last_pos0.max(segment.last_pos0);
                    }
                    cur
                }
                Some(cur) => {
                    self.close_contig(cur)?;
                    Self::open_and_absorb(&mut segment, &self.grid, &self.seed)
                }
                None => Self::open_and_absorb(&mut segment, &self.grid, &self.seed),
            };
            if let ContigAcc::Rows(rows) = &mut segment.acc {
                // Absorb one site and immediately advance the watermark.
                // Accumulating all rows into bins first recreated the old
                // per-site memory blow-up inside a single chunk (V27-18).
                let mut rows = rows.drain(..).peekable();
                while let Some(row) = rows.next() {
                    current
                        .acc
                        .push_site(row.pos0, &row.counts, row.all_variant);
                    current.first_pos0 = current.first_pos0.min(row.pos0);
                    current.last_pos0 = current.last_pos0.max(row.pos0);
                    // Do not close at the last position of a window until
                    // the next position differs: same-POS duplicates may be
                    // split across chunks and must share the window.
                    if rows.peek().is_none_or(|next| next.pos0 != row.pos0) {
                        self.emit_complete_windows(&mut current)?;
                    }
                }
            } else {
                self.emit_complete_windows(&mut current)?;
            }
            self.current = Some(current);
        }
        Ok(())
    }

    fn finish(mut self) -> Result<u64, EngineError> {
        self.finalize_ld()?;
        let mut ld = self.ld_merged.take();
        if let Some(executor) = &mut self.ld_executor {
            let pairs_path = self.ld.as_ref().and_then(|p| p.pairs_path.clone());
            let mut merge = |mut set: LdResult| -> Result<(), EngineError> {
                if let Some(path) = &pairs_path {
                    write_pairs(&mut self.pairs_file, path, &set.pairs)?;
                    set.pairs.clear();
                }
                match &mut ld {
                    Some(dst) => {
                        dst.cohort.merge(&set.cohort);
                        for (dst, src) in dst.groups.iter_mut().zip(&set.groups) {
                            dst.merge(src);
                        }
                        dst.pairs.extend(set.pairs);
                    }
                    None => ld = Some(set),
                }
                Ok(())
            };
            executor.finish(&mut merge)?;
        }
        if let Some(current) = self.current.take() {
            self.close_contig(current)?;
        }
        if self.levels.genome {
            let bp = self.genome_span;
            let denom_bp = self.genome_denom_bp;
            let totals = self.genome_totals.clone();
            self.write_row("genome", 0, bp, bp, denom_bp, &totals)?;
        }
        // The all-sites check (04-algorithms.md §3): a file with no
        // invariant sites must be declared with --variant-only or
        // --bypass-invariant-check.
        if matches!(
            self.layout,
            Some(WindowLayoutOptions {
                invariant: InvariantMode::Require
            })
        ) && self.n_no_alt_records == 0
        {
            return Err(EngineError::usage(format!(
                "input has no invariant sites ({} records, none with ALT = \".\"): declare \
                 --variant-only for an all-variant file, or --bypass-invariant-check to compute \
                 over present sites only (both overestimate π/dxy denominators)",
                self.n_records
            )));
        }
        let mut pending_commit = None;
        let pairs_path = self.ld.as_ref().and_then(|plan| plan.pairs_path.clone());
        // The streaming pair table is *not* committed here: the driver still
        // has to validate the input's completeness (`job.finish`, which
        // detects a record truncated at the end of a BCF). Committing before
        // that check published a truncated `.ld` over the previous run's
        // good file while the run then failed (doc_v27 V27-02). The
        // `AtomicFile` travels out in the slot; the driver commits it only
        // after every check passed, and dropping it unpublishes the temp.
        // doc_v28 V28-23: `finish` flushes, closes and hands back the
        // unpublished temporary; the driver commits it after every check.
        if let Some(writer) = self.pairs_file.take() {
            pending_commit = Some(writer.finish().map_err(EngineError::Io)?);
        }
        *self.out.lock().expect("invariant: stats slot not poisoned") = Some(SinkOutput {
            ld,
            windows: std::mem::take(&mut self.windows),
            sfs: self.sfs.take(),
            fstats: self.fstats.take(),
            pending_commit,
            window_session: self.output_session.take(),
            pairs_path,
            n_invariant_sites: self.n_invariant_sites,
            n_no_alt_records: self.n_no_alt_records,
            n_windows: self.n_windows,
            n_records: self.n_records,
            n_multiallelic_skipped: self.n_multiallelic_skipped,
        });
        Ok(self.n_windows)
    }
}

/// Opens the run's pair table and writes its header immediately. Publication
/// remains deferred: `AtomicFile` is committed only after input validation,
/// and dropping it removes the temporary on every failure path.
pub(super) fn open_pairs(path: &std::path::Path) -> Result<Option<AtomicWriter>, EngineError> {
    let mut writer = AtomicWriter::with_capacity(path, 1 << 16).map_err(EngineError::Io)?;
    std::io::Write::write_all(&mut writer, crate::stats::writers::ld::PAIRS_HEADER)
        .map_err(EngineError::Io)?;
    Ok(Some(writer))
}

/// Appends one contig's pair rows to the streaming table. The rows and the
/// header are the same bytes `writers::ld` writes for the in-memory report,
/// so a streamed run and an accumulated run produce identical files
/// (doc_v23 V23-14).
fn write_pairs(
    file: &mut Option<AtomicWriter>,
    path: &std::path::Path,
    rows: &[(Vec<u8>, u64, u64, f64)],
) -> Result<(), EngineError> {
    if rows.is_empty() {
        return Ok(());
    }
    if file.is_none() {
        *file = open_pairs(path)?;
    }
    let sink = file.as_mut().expect("invariant: opened above");
    let mut buf = String::with_capacity(rows.len() * 40);
    for (contig, a, b, r2) in rows {
        crate::stats::writers::ld::render_pair_row(
            &mut buf,
            &String::from_utf8_lossy(contig),
            *a,
            *b,
            *r2,
        );
    }
    std::io::Write::write_all(sink, buf.as_bytes()).map_err(EngineError::Io)
}

/// One LD executor per run. Contig analyses share this pool, so their
/// parallelism is additive rather than `(contig tasks) × (new pool per
/// task)` (doc_v27 V27-16). Results are buffered by contig sequence before
/// the sink merges them: a pool may finish contigs out of order, but the
/// `.ld` bytes remain in contig order.
pub(super) struct LdExecutor {
    pool: Option<rayon::ThreadPool>,
    results: std::sync::mpsc::Receiver<(usize, Result<LdResult, String>)>,
    sender: std::sync::mpsc::Sender<(usize, Result<LdResult, String>)>,
    abort: reed_engine::abort::AbortHandle,
    signal: reed_engine::abort::AbortSignal,
    next_sequence: usize,
    next_merge: usize,
    completed: std::collections::BTreeMap<usize, Result<LdResult, String>>,
    /// Set only after every result was merged successfully. Dropping an
    /// unfinished executor aborts; dropping a finished one must not flip the
    /// run's shared abort flag (V27-17).
    finished: bool,
}

impl LdExecutor {
    pub(super) fn new(
        worker_threads: usize,
        abort: reed_engine::abort::AbortHandle,
        signal: reed_engine::abort::AbortSignal,
    ) -> Result<Self, EngineError> {
        let (sender, results) = std::sync::mpsc::channel();
        // Sequential runs execute contigs inline: no extra thread at all
        // (doc_v27 V27-37). Parallel runs build exactly one auxiliary pool.
        let pool = if worker_threads == 0 {
            None
        } else {
            Some(
                rayon::ThreadPoolBuilder::new()
                    .num_threads(worker_threads)
                    .build()
                    .map_err(|e| {
                        EngineError::job(format!("cannot create the LD worker pool: {e}"))
                    })?,
            )
        };
        Ok(Self {
            pool,
            results,
            sender,
            abort,
            signal,
            next_sequence: 0,
            next_merge: 0,
            completed: std::collections::BTreeMap::new(),
            finished: false,
        })
    }

    /// Number of tasks that have not yet been merged by the sink.
    fn in_flight(&self) -> usize {
        self.next_sequence - self.next_merge
    }

    /// Queues one contig, or computes it inline in the sequential path.
    fn submit(
        &mut self,
        matrix: reed_stats::ContigMatrix,
        max_bp: u64,
        edges: Vec<u64>,
        collect: bool,
    ) {
        let sequence = self.next_sequence;
        self.next_sequence += 1;
        let Some(pool) = &self.pool else {
            let result = decay_contig(
                &matrix,
                max_bp,
                &edges,
                collect,
                &ThreadBudget {
                    threads: 1,
                    mem: reed_engine::budget::MemBudget::from_bytes(u64::MAX),
                },
            );
            self.completed.insert(
                sequence,
                Ok(LdResult {
                    cohort: result.cohort,
                    groups: result.groups,
                    pairs: result
                        .pairs
                        .into_iter()
                        .map(|(a, b, r2)| (matrix.contig().to_vec(), a, b, r2))
                        .collect(),
                }),
            );
            return;
        };
        let sender = self.sender.clone();
        let signal = self.signal.clone();
        pool.spawn(move || {
            let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                decay_contig_with_abort(&matrix, max_bp, &edges, collect, true, &signal)
            }));
            let result = match result {
                Ok(Ok(set)) => Ok(LdResult {
                    cohort: set.cohort,
                    groups: set.groups,
                    pairs: set
                        .pairs
                        .into_iter()
                        .map(|(a, b, r2)| (matrix.contig().to_vec(), a, b, r2))
                        .collect(),
                }),
                Ok(Err(err)) => Err(err),
                Err(_) => Err("the LD decay worker panicked".to_string()),
            };
            // If the sink is gone, nobody needs this result; dropping the
            // sender is sufficient and does not strand the worker.
            let _ = sender.send((sequence, result));
        });
    }

    /// Merges completed results in contig order until one is still running.
    /// This is the old `join_ld` behavior, without assuming task completion
    /// order is contig order.
    fn drain_ordered(
        &mut self,
        mut merge: impl FnMut(LdResult) -> Result<(), EngineError>,
    ) -> Result<(), EngineError> {
        while self.completed.contains_key(&self.next_merge) {
            let (_, result) = self.completed.pop_first().expect("invariant: key checked");
            self.next_merge += 1;
            let result = result.map_err(EngineError::job)?;
            merge(result)?;
        }
        Ok(())
    }

    /// Waits for the next still-running contig and merges it. The wait is
    /// abortable, so an error in the input pipeline cannot leave the sink
    /// blocked forever on a result channel.
    fn merge_one(
        &mut self,
        mut merge: impl FnMut(LdResult) -> Result<(), EngineError>,
    ) -> Result<(), EngineError> {
        while self.in_flight() > 0 {
            if let Some(result) = self.completed.remove(&self.next_merge) {
                self.next_merge += 1;
                merge(result.map_err(EngineError::job)?)?;
                return Ok(());
            }
            if self.signal.is_aborted() {
                return Err(EngineError::Aborted);
            }
            match self
                .results
                .recv_timeout(std::time::Duration::from_millis(10))
            {
                Ok((sequence, result)) => {
                    self.completed.insert(sequence, result);
                }
                Err(std::sync::mpsc::RecvTimeoutError::Timeout) => {}
                Err(std::sync::mpsc::RecvTimeoutError::Disconnected) => {
                    return Err(EngineError::job(
                        "the LD result channel closed unexpectedly",
                    ));
                }
            }
        }
        Ok(())
    }

    /// Drains all tasks in contig order. Used at successful finish.
    fn finish(
        &mut self,
        mut merge: impl FnMut(LdResult) -> Result<(), EngineError>,
    ) -> Result<(), EngineError> {
        while self.in_flight() > 0 {
            self.merge_one(&mut merge)?;
        }
        self.drain_ordered(merge)?;
        self.finished = true;
        Ok(())
    }
}

impl Drop for LdExecutor {
    fn drop(&mut self) {
        // An error returned by the library drops the sink before the caller
        // continues. Cancel every block loop first, then let pool drop join
        // its workers: no LD matrix scan survives the failed run (V27-17).
        if !self.finished {
            self.abort.trigger();
        }
        // Wake a receiver that is racing the abort before the pool drops.
        while self.results.try_recv().is_ok() {}
        self.completed.clear();
        self.pool.take();
    }
}

/// The contig worker with an abort check between 4096-site blocks. The
/// nested block parallelism stays inside the run's one Rayon pool.
fn decay_contig_with_abort(
    matrix: &reed_stats::ContigMatrix,
    max_bp: u64,
    edges: &[u64],
    collect: bool,
    parallel: bool,
    signal: &reed_engine::abort::AbortSignal,
) -> Result<reed_stats::LdDecaySet, String> {
    use rayon::prelude::*;
    let n = matrix.n_sites();
    let mut merged = reed_stats::LdDecaySet::new(edges, matrix.layout().n_groups());
    if n < 2 {
        return Ok(merged);
    }
    let ranges: Vec<std::ops::Range<usize>> = (0..n)
        .step_by(crate::stats::run::LD_BLOCK)
        .map(|start| start..(start + crate::stats::run::LD_BLOCK).min(n))
        .collect();
    let compute = |parallel: bool| -> Result<Vec<reed_stats::LdDecaySet>, String> {
        let one = |range: &std::ops::Range<usize>| -> Result<reed_stats::LdDecaySet, String> {
            if signal.is_aborted() {
                return Err("LD decay aborted".to_string());
            }
            Ok(reed_stats::decay_range(
                matrix,
                range.clone(),
                max_bp,
                edges,
                collect,
            ))
        };
        if parallel {
            ranges.par_iter().map(one).collect::<Result<Vec<_>, _>>()
        } else {
            ranges.iter().map(one).collect::<Result<Vec<_>, _>>()
        }
    };
    let parts = compute(parallel)?;
    for part in &parts {
        merged.merge(part);
    }
    Ok(merged)
}

impl StatsSink {
    /// Hands the pending contig's matrix to a decay worker: the analysis
    /// runs off the reducer thread (it overlaps the decode), and `finish`
    /// joins the jobs in contig order — the merge order the deterministic
    /// output requires (docs/rounds/v6/02-models.md §2.2). In-flight jobs
    /// are capped at `budget.auxiliary()`: M contigs cannot keep M matrices
    /// resident (doc_v8 D8-2, the same cap as the prune path's D7-19) —
    /// the oldest job is joined when the cap is reached.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] when a capped-out join hits a panicked worker
    /// (fail here, not at `finish`, so the run stops at the cause).
    pub(super) fn finalize_ld(&mut self) -> Result<(), EngineError> {
        let Some(matrix) = self.ld_matrix.take() else {
            return Ok(());
        };
        let Some(plan) = &self.ld else { return Ok(()) };
        let edges = plan.decay_edges();
        let max_bp = plan.decay_max_bp.unwrap_or(100_000);
        let collect = plan.pairs;
        // Bound resident matrices, not the number of OS threads. The oldest
        // contig is merged before another one is queued, so a run with many
        // contigs still has one pool and at most `max_ld_jobs` matrices.
        if self
            .ld_executor
            .as_ref()
            .is_some_and(|executor| executor.in_flight() >= self.max_ld_jobs)
        {
            let mut merged = self.ld_merged.take();
            self.join_ld(&mut merged)?;
            self.ld_merged = merged;
        }
        if let Some(executor) = &mut self.ld_executor {
            executor.submit(matrix, max_bp, edges, collect);
        }
        Ok(())
    }

    /// Merges completed LD results in contig order. `out` starts from and
    /// returns to the caller's accumulator, so mid-run capped merges and
    /// `finish` share exactly the old join semantics.
    pub(super) fn join_ld(&mut self, out: &mut Option<LdResult>) -> Result<(), EngineError> {
        let pairs_path = self.ld.as_ref().and_then(|p| p.pairs_path.clone());
        let mut merge = |mut set: LdResult| -> Result<(), EngineError> {
            if let Some(path) = &pairs_path {
                write_pairs(&mut self.pairs_file, path, &set.pairs)?;
                set.pairs.clear();
            }
            match out {
                Some(dst) => {
                    dst.cohort.merge(&set.cohort);
                    for (dst, src) in dst.groups.iter_mut().zip(&set.groups) {
                        dst.merge(src);
                    }
                    dst.pairs.extend(set.pairs);
                }
                None => *out = Some(set),
            }
            Ok(())
        };
        if let Some(executor) = &mut self.ld_executor {
            executor.merge_one(&mut merge)?;
        }
        Ok(())
    }

    /// Opens a contig's state and absorbs its first segment (folding the
    /// segment's bins into the contig totals and the bin map).
    pub(super) fn open_and_absorb(
        segment: &mut SegmentPartial,
        grid: &GridMode,
        seed: &[StatAcc],
    ) -> CurrentContig {
        // Grid from position 1: `[1, w]`, `[1+s, w+s]`, … (04 §4).
        let anchor = 0;
        let mut cur = CurrentContig {
            name: segment.contig.clone(),
            acc: match grid {
                GridMode::Bp { width, .. } => {
                    ContigAcc::Bins(ContigBins::new(BpBinning::new(*width), seed))
                }
                GridMode::Regions(_regions) => {
                    // Workers share the contig's accumulator by Arc. Adopt
                    // it here rather than allocating a second copy: later
                    // partials with the same Arc must not be merged again.
                    let ContigAcc::Regions(acc) = &segment.acc else {
                        unreachable!("invariant: Regions grid implies region accumulators");
                    };
                    ContigAcc::Regions(Arc::clone(acc))
                }
                GridMode::None => ContigAcc::None,
            },
            totals: seed.to_vec(),
            first_pos0: segment.first_pos0,
            last_pos0: if matches!(segment.acc, ContigAcc::Rows(_)) {
                // Raw rows advance the watermark one position at a time;
                // using the whole segment's last position here emits the
                // window before its earlier rows are absorbed (V27-18).
                segment.first_pos0
            } else {
                segment.last_pos0
            },
            anchor,
            next_start: anchor,
            span_end: 0,
        };
        if matches!(cur.acc, ContigAcc::Bins(_)) {
            segment.acc.fold(&mut cur.totals);
        }
        if !matches!(segment.acc, ContigAcc::Rows(_)) {
            cur.acc.merge(&mut segment.acc);
        }
        cur
    }

    /// Emits every window that the watermark proves complete: with sorted
    /// input, no future site can land in `[s, s+size)` once a position at or
    /// beyond `s+size` has been seen (§8.3). The check demands a *strictly
    /// later* position rather than merely reaching the window's last one
    /// (doc_v15 D15-19): the VCF specification permits multiple records at
    /// the same POS ("Positions are sorted numerically, in increasing
    /// order ... It is permitted to have multiple records with the same
    /// POS", VCFv4.3 Fixed fields), and a same-POS duplicate that straddles
    /// a BGZF chunk seam must still land in the window it belongs to.
    /// Consumed bins are freed. The contig length does not gate this loop —
    /// only the trailing windows in `close_contig` need it. Regions grids
    /// emit nothing here: their windows are the BED intervals, known up
    /// front.
    #[allow(clippy::unnecessary_wraps)] // callers share the Result-returning shape of the window writers
    pub(super) fn emit_complete_windows(
        &mut self,
        cur: &mut CurrentContig,
    ) -> Result<(), EngineError> {
        let GridMode::Bp { size, step, width } = self.grid else {
            return Ok(());
        };
        let watermark = cur.last_pos0;
        if !self.levels.window && !self.levels.contig {
            // Genome-only output has no window or contig row to emit. Only
            // the observed span is needed; walking every possible window is
            // wasted work, especially with --window 1 (doc_v27 V27-30).
            cur.span_end = watermark.saturating_add(1);
            cur.next_start = cur.span_end;
            return Ok(());
        }
        let mut s = cur.next_start;
        let ContigAcc::Bins(bins) = &mut cur.acc else {
            return Ok(());
        };
        while s.saturating_add(size) <= watermark {
            if self.skip_empty {
                let Some(next_bin) = bins.next_bin_at_or_after(s / width) else {
                    break;
                };
                // The next touched bin determines the earliest window that
                // can be nonempty; windows before it need no iteration.
                let candidate = next_bin.saturating_mul(width);
                let earliest = candidate.saturating_sub(size - 1);
                if earliest > s {
                    s = earliest.div_ceil(step) * step;
                }
            }
            let end = s + size;
            let range = (s / width)..end.div_ceil(width);
            let empty = self.skip_empty && !bins.any_bin_in(range.clone());
            // The window span advances whether or not the row is written, so
            // the contig row's bp stays correct when `--levels` drops the
            // window rows (doc_v7 D7-13).
            cur.span_end = end;
            if self.levels.window && !empty {
                let mut accs = std::mem::take(&mut self.window_scratch);
                let owned;
                let accs_ref = if self.compact_window_accumulators {
                    bins.finish_bins_into(range, &mut accs);
                    &accs
                } else {
                    owned = bins.finish_bins(range);
                    &owned
                };
                self.write_window_row(
                    &String::from_utf8_lossy(&cur.name),
                    s + 1,
                    end,
                    end - s,
                    end - s,
                    accs_ref,
                )?;
                self.window_scratch = accs;
                self.n_windows += 1;
            }
            s = s.saturating_add(step);
        }
        cur.next_start = s;
        bins.drop_bins_below(s / width);
        Ok(())
    }

    /// Closes a contig: emits the trailing windows (start at or before the
    /// last site, clipped to the contig length — the contig length is
    /// required here and is never fabricated, 01-review D7), then the
    /// contig-level row (§8.7), then folds into the genome total.
    pub(super) fn close_contig(&mut self, mut cur: CurrentContig) -> Result<(), EngineError> {
        if !self.closed.insert(cur.name.clone()) {
            return Err(EngineError::job(format!(
                "contig {} revisited after a later contig: input is not contig-sorted",
                String::from_utf8_lossy(&cur.name)
            )));
        }
        let contig = String::from_utf8_lossy(&cur.name).into_owned();
        match self.grid.clone() {
            GridMode::Bp { size, step, width } => {
                let len = self.contig_len.get(&cur.name).copied().ok_or_else(|| {
                    EngineError::job(format!(
                        "contig {} has no ##contig length; window clipping needs it (pass --contig-lengths)",
                        String::from_utf8_lossy(&cur.name)
                    ))
                })?;
                let ContigAcc::Bins(bins) = &cur.acc else {
                    unreachable!("invariant: Bp grid implies bin accumulators");
                };
                let mut s = cur.next_start;
                while s < len && s <= cur.last_pos0 {
                    if self.skip_empty {
                        let Some(next_bin) = bins.next_bin_at_or_after(s / width) else {
                            break;
                        };
                        let candidate = next_bin.saturating_mul(width);
                        let earliest = candidate.saturating_sub(size - 1);
                        if earliest > s {
                            s = earliest.div_ceil(step) * step;
                        }
                    }
                    let end = s.saturating_add(size).min(len);
                    let range = (s / width)..end.div_ceil(width);
                    let empty = self.skip_empty && !bins.any_bin_in(range.clone());
                    cur.span_end = end;
                    if self.levels.window && !empty {
                        let mut accs = std::mem::take(&mut self.window_scratch);
                        let owned;
                        let accs_ref = if self.compact_window_accumulators {
                            bins.finish_bins_into(range, &mut accs);
                            &accs
                        } else {
                            owned = bins.finish_bins(range);
                            &owned
                        };
                        self.write_window_row(&contig, s + 1, end, end - s, end - s, accs_ref)?;
                        self.window_scratch = accs;
                        self.n_windows += 1;
                    }
                    s = s.saturating_add(step);
                }
            }
            GridMode::None => {
                // No window rows exist for an SFS/f-statistics/LD-only
                // request. There is consequently no trailing window to clip
                // and no contig length to require (doc_v27 V27-07); the
                // span is the observed input extent.
                cur.span_end = cur.last_pos0.saturating_add(1);
            }
            GridMode::Regions(_) => {
                // The windows are the BED intervals in input order; every
                // region is emitted (skip-empty drops untouched ones).
                let ContigAcc::Regions(regions) = &cur.acc else {
                    unreachable!("invariant: Regions grid implies region accumulators");
                };
                let regions = regions
                    .lock()
                    .expect("invariant: region accumulator not poisoned");
                let intervals: Vec<reed_core::coord::Interval0> =
                    regions.regions().windows().to_vec();
                if intervals.is_empty() {
                    // A contig absent from BED has zero BED windows. Do not
                    // emit the old `start1 = 1, end = 0` row: its end was
                    // before its start (doc_v27 V27-23). The observed span
                    // keeps any contig/genome row well-formed; |S| stays 0.
                    cur.span_end = cur.last_pos0.saturating_add(1);
                }
                for (idx, interval) in intervals.iter().enumerate() {
                    let start0 = interval.start.0;
                    let end0 = interval.end.0;
                    // The span covers every region even when the row is
                    // dropped by `--levels` (doc_v7 D7-13).
                    cur.span_end = cur.span_end.max(end0);
                    if !self.levels.window || (self.skip_empty && !regions.touched(idx)) {
                        continue;
                    }
                    let bp = end0.saturating_sub(start0);
                    let accs = regions.accs(idx);
                    self.write_window_row(&contig, start0 + 1, end0, bp, bp, accs)?;
                    self.n_windows += 1;
                }
            }
        }
        if matches!(cur.acc, ContigAcc::Regions(_)) {
            // Region partials share their accumulator across chunks; fold
            // its full state exactly once at contig close, never once per
            // chunk (doc_v27 V27-19).
            cur.acc.fold(&mut cur.totals);
        }
        if self.levels.contig {
            let bp = cur.span_end.saturating_sub(cur.anchor);
            let denom_bp = self.denominator_for_contig(&cur);
            self.write_row(
                &contig,
                cur.anchor + 1,
                cur.span_end,
                bp,
                denom_bp,
                &cur.totals,
            )?;
        }
        for (dst, src) in self.genome_totals.iter_mut().zip(&cur.totals) {
            dst.merge(src);
        }
        self.genome_span += cur.span_end.saturating_sub(cur.anchor);
        self.genome_denom_bp += self.denominator_for_contig(&cur);
        Ok(())
    }

    /// The variant-only denominator of a contig row. BED windows use the
    /// union of their intervals clipped to the declared contig length, not
    /// the whole contig: an unsampled interval is not an absent monomorphic
    /// position (doc_v27 V27-11). Without a declared length the observed
    /// span is the only honest fallback; no length is fabricated.
    fn denominator_for_contig(&self, cur: &CurrentContig) -> u64 {
        let span = cur.span_end.saturating_sub(cur.anchor);
        match &self.grid {
            GridMode::Regions(regions) => {
                let end = self
                    .contig_len
                    .get(&cur.name)
                    .copied()
                    .unwrap_or(cur.span_end);
                regions.get(&cur.name).map_or(0, |set| {
                    set.accessible_bp(Interval0::new(Pos0(0), Pos0(end)))
                })
            }
            GridMode::Bp { .. } | GridMode::None => {
                self.contig_len.get(&cur.name).copied().unwrap_or(span)
            }
        }
    }

    /// Emits a completed window and releases its large distribution vectors
    /// when the CLI requested compact output. The contig/genome totals are
    /// separate full accumulators, so this cannot affect their merges.
    pub(super) fn write_window_row(
        &mut self,
        contig: &str,
        start1: u64,
        end0: u64,
        bp: u64,
        denom_bp: u64,
        accs: &[StatAcc],
    ) -> Result<(), EngineError> {
        if self.compact_window_accumulators {
            let variant_only = self
                .layout
                .as_ref()
                .is_some_and(|layout| layout.invariant == InvariantMode::VariantOnly);
            let outputs = accs
                .iter()
                .zip(self.kinds.iter())
                .map(|(acc, kind)| {
                    reed_stats::request::StatOutput::from_acc(acc, kind, bp, denom_bp, variant_only)
                })
                .collect();
            // Do not call `write_row` here: cloning the full accumulator
            // and then clearing it would still churn the heap once per
            // window. Only scalar outputs enter the CLI's output session.
            let row = WindowRow {
                contig: contig.to_string(),
                start1,
                end: end0,
                bp,
                denom_bp,
                accs: Vec::new(),
                outputs,
            };
            if let Some(session) = &mut self.output_session {
                session.write_row(&row)?;
            } else {
                self.windows.push(row);
            }
            Ok(())
        } else {
            self.write_row(contig, start1, end0, bp, denom_bp, accs)
        }
    }

    /// Records one row at any level. `bp` is the public span; `denom_bp`
    /// is independent and is used only for variant-only π/dxy. Ordinary
    /// windows pass the same value for both; contig/genome rows pass their
    /// declared contig length as the denominator (doc_v27 V27-11).
    pub(super) fn write_row(
        &mut self,
        contig: &str,
        start1: u64,
        end0: u64,
        bp: u64,
        denom_bp: u64,
        accs: &[StatAcc],
    ) -> Result<(), EngineError> {
        if let Some(session) = &mut self.output_session {
            let variant_only = self
                .layout
                .as_ref()
                .is_some_and(|layout| layout.invariant == InvariantMode::VariantOnly);
            let outputs = accs
                .iter()
                .zip(self.kinds.iter())
                .map(|(acc, kind)| {
                    reed_stats::request::StatOutput::from_acc(acc, kind, bp, denom_bp, variant_only)
                })
                .collect();
            let row = WindowRow {
                contig: contig.to_string(),
                start1,
                end: end0,
                bp,
                denom_bp,
                accs: Vec::new(),
                outputs,
            };
            session.write_row(&row)
        } else {
            self.windows.push(WindowRow {
                contig: contig.to_string(),
                start1,
                end: end0,
                bp,
                denom_bp,
                accs: accs.to_vec(),
                outputs: Vec::new(),
            });
            Ok(())
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use reed_core::count::SiteGroupCounts;
    use reed_stats::request::{StatKind, build_accs};

    fn site_row(n_ref: u32, n_alt: u32) -> Vec<SiteGroupCounts> {
        vec![SiteGroupCounts {
            n_ref_alleles: n_ref,
            n_alt_alleles: n_alt,
            ..SiteGroupCounts::default()
        }]
    }

    /// One worker segment: one contig, one site, bins built exactly like the
    /// worker builds them (`BpBinning::new(width)`).
    fn segment(pos0: u64, n_ref: u32, n_alt: u32, seed: &[StatAcc]) -> SegmentPartial {
        let mut bins = ContigBins::new(BpBinning::new(100), seed);
        bins.push_site(pos0, &site_row(n_ref, n_alt), n_alt > 0);
        SegmentPartial {
            contig: b"chr1".to_vec(),
            acc: ContigAcc::Bins(bins),
            first_pos0: pos0,
            last_pos0: pos0,
        }
    }

    /// Builds a sink without going through the public driver, so a chunk
    /// seam can be placed at an exact position rather than guessed from the
    /// compressor's block layout (the original V27-06 repro depends on that
    /// layout and therefore cannot be a stable unit fixture).
    fn new_sink(
        seed: &[StatAcc],
        grid: GridMode,
        contig_lengths: HashMap<Vec<u8>, u64>,
    ) -> (StatsSink, Arc<std::sync::Mutex<Option<SinkOutput>>>) {
        let out = Arc::new(std::sync::Mutex::new(None));
        let sink = StatsSink {
            sfs: None,
            fstats: None,
            ld_layout: None,
            ld_matrix: None,
            ld_executor: None,
            ld_merged: None,
            max_ld_jobs: 1,
            ld: None,
            pairs_file: None,
            grid,
            kinds: Vec::new().into(),
            seed: seed.to_vec(),
            levels: Levels {
                window: true,
                contig: true,
                genome: true,
            },
            compact_window_accumulators: false,
            skip_empty: false,
            layout: None,
            windows: Vec::new(),
            window_scratch: Vec::new(),
            output_session: None,
            n_invariant_sites: 0,
            n_no_alt_records: 0,
            contig_len: contig_lengths,
            current: None,
            closed: HashSet::new(),
            genome_totals: seed.to_vec(),
            genome_span: 0,
            genome_denom_bp: 0,
            n_windows: 0,
            n_records: 0,
            n_multiallelic_skipped: 0,
            out: Arc::clone(&out),
        };
        (sink, out)
    }

    fn accept_site(sink: &mut StatsSink, id: u64, pos0: u64, seed: &[StatAcc]) {
        let mut partial = StatsPartial {
            segments: vec![segment(pos0, 1, 1, seed)],
            n_records: 1,
            ..StatsPartial::default()
        };
        sink.accept(ChunkId(id), &mut partial)
            .expect("sorted input must be accepted");
    }

    /// doc_v27 V27-06: with `step > size`, a sorted position between two
    /// windows is not an ordering error, including at a chunk seam. A
    /// subsequent position in the next window must still be accumulated.
    #[test]
    fn a_site_between_sparse_windows_is_legal_across_a_chunk_seam() {
        let seed = build_accs(
            &[StatKind::Pi {
                group: 0,
                per_bp: false,
            }],
            &[2],
        )
        .expect("seed");
        let (mut sink, out) = new_sink(
            &seed,
            GridMode::Bp {
                size: 100,
                step: 1_000_000,
                width: 100,
            },
            HashMap::from([(b"chr1".to_vec(), 1_100_000_u64)]),
        );
        // [1,100] closes, then the legal inter-window sites straddle a seam.
        for (id, pos0) in [50_u64, 150, 200, 1_000_050].into_iter().enumerate() {
            accept_site(&mut sink, id as u64, pos0, &seed);
        }
        sink.finish().expect("finish");
        let output = out.lock().expect("slot").take().expect("sink output");
        let windows: Vec<(u64, u64)> = output
            .windows
            .iter()
            .filter(|w| w.start1 != 1 || w.end != 1_100_000)
            .filter(|w| w.contig == "chr1")
            .map(|w| (w.start1, w.end))
            .collect();
        assert!(windows.contains(&(1, 100)), "first window: {windows:?}");
        assert!(
            windows.contains(&(1_000_001, 1_000_100)),
            "the later window must not be lost: {windows:?}"
        );
    }

    /// Removing `next_start` from the seam check must not drop a site that
    /// belongs to an already-emitted window. The two records straddle the
    /// seam at exactly the window-closing position.
    #[test]
    fn a_same_pos_duplicate_still_lands_in_the_closed_window_when_step_exceeds_size() {
        let seed = build_accs(
            &[StatKind::Pi {
                group: 0,
                per_bp: false,
            }],
            &[2],
        )
        .expect("seed");
        let (mut sink, out) = new_sink(
            &seed,
            GridMode::Bp {
                size: 100,
                step: 1_000_000,
                width: 100,
            },
            HashMap::from([(b"chr1".to_vec(), 1_100_000_u64)]),
        );
        for (id, pos0) in [50_u64, 99, 99].into_iter().enumerate() {
            accept_site(&mut sink, id as u64, pos0, &seed);
        }
        sink.finish().expect("finish");
        let output = out.lock().expect("slot").take().expect("sink output");
        let first = output
            .windows
            .iter()
            .find(|w| w.start1 == 1 && w.end == 100)
            .expect("first window");
        let StatAcc::Pi(pi) = &first.accs[0] else {
            panic!("pi acc");
        };
        assert_eq!(pi.site_tally().0, 3, "all three records belong to [1,100]");
    }

    /// The seam still rejects a genuine backwards move; the fix only
    /// separates the two cursors, it does not disable input validation.
    #[test]
    fn a_backwards_move_is_still_rejected_when_step_exceeds_size() {
        let seed = build_accs(
            &[StatKind::Pi {
                group: 0,
                per_bp: false,
            }],
            &[2],
        )
        .expect("seed");
        let (mut sink, _out) = new_sink(
            &seed,
            GridMode::Bp {
                size: 100,
                step: 1_000_000,
                width: 100,
            },
            HashMap::from([(b"chr1".to_vec(), 1_100_000_u64)]),
        );
        for (id, pos0) in [50_u64, 200].into_iter().enumerate() {
            accept_site(&mut sink, id as u64, pos0, &seed);
        }
        let mut bad = StatsPartial {
            segments: vec![segment(150, 1, 1, &seed)],
            n_records: 1,
            ..StatsPartial::default()
        };
        let err = sink.accept(ChunkId(2), &mut bad).expect_err("must reject");
        assert!(err.to_string().contains("strictly ascending"), "{err}");
    }

    /// doc_v15 D15-19: the ordered sink receives one segment per chunk, so a
    /// *legal* same-POS duplicate (VCFv4.3: "It is permitted to have
    /// multiple records with the same POS") can straddle a seam at a
    /// window's last position. Window [1,100] closes only once a strictly
    /// later position is seen, so both copies must land inside it; the old
    /// rule emitted the window at POS 100 and silently dropped the second
    /// copy from every window row.
    #[test]
    fn same_pos_seam_duplicate_lands_in_its_window() {
        let seed = build_accs(
            &[StatKind::Pi {
                group: 0,
                per_bp: false,
            }],
            &[2],
        )
        .expect("seed");
        let out = Arc::new(std::sync::Mutex::new(None));
        let mut sink = StatsSink {
            sfs: None,
            fstats: None,
            ld_layout: None,
            ld_matrix: None,
            ld_executor: None,
            ld_merged: None,
            max_ld_jobs: 1,
            ld: None,
            pairs_file: None,
            grid: GridMode::Bp {
                size: 100,
                step: 100,
                width: 100,
            },
            kinds: Vec::new().into(),
            seed: seed.clone(),
            levels: Levels {
                window: true,
                contig: true,
                genome: false,
            },
            compact_window_accumulators: false,
            skip_empty: false,
            layout: None,
            windows: Vec::new(),
            window_scratch: Vec::new(),
            output_session: None,
            n_invariant_sites: 0,
            n_no_alt_records: 0,
            contig_len: HashMap::from([(b"chr1".to_vec(), 1000_u64)]),
            current: None,
            closed: HashSet::new(),
            genome_totals: seed.clone(),
            genome_span: 0,
            genome_denom_bp: 0,
            n_windows: 0,
            n_records: 0,
            n_multiallelic_skipped: 0,
            out: Arc::clone(&out),
        };
        // Chunk 0 ends with the first copy at POS 100 (pos0 99).
        let mut chunk0 = StatsPartial {
            segments: vec![segment(99, 1, 1, &seed)],
            n_records: 1,
            ..StatsPartial::default()
        };
        sink.accept(ChunkId(0), &mut chunk0)
            .expect("accept chunk 0");
        // Chunk 1 opens with the duplicate: same POS, still sorted.
        let mut chunk1 = StatsPartial {
            segments: vec![segment(99, 0, 2, &seed)],
            n_records: 1,
            ..StatsPartial::default()
        };
        sink.accept(ChunkId(1), &mut chunk1)
            .expect("accept chunk 1 (same POS across the seam)");
        // Chunk 2 pushes past the window so [1,100] closes through the
        // watermark, then the run finishes.
        let mut chunk2 = StatsPartial {
            segments: vec![segment(249, 1, 1, &seed)],
            n_records: 1,
            ..StatsPartial::default()
        };
        sink.accept(ChunkId(2), &mut chunk2)
            .expect("accept chunk 2");
        sink.finish().expect("finish");
        let output = out.lock().expect("slot").take().expect("sink output");

        let window = output
            .windows
            .iter()
            .find(|w| w.contig == "chr1" && w.start1 == 1 && w.end == 100)
            .expect("window [1,100]");
        let StatAcc::Pi(pi) = &window.accs[0] else {
            panic!("pi acc");
        };
        assert_eq!(
            pi.site_tally().0,
            2,
            "both same-POS copies must be inside window [1,100]"
        );
        assert_eq!(
            pi.counts().diffs,
            1,
            "the het copy contributes one REF x ALT pair, the hom-alt one none"
        );

        // Window rows sum to the contig row.
        let contig = output
            .windows
            .iter()
            .find(|w| w.contig == "chr1" && w.start1 == 1 && w.end == 300)
            .expect("contig row spans the emitted windows");
        let StatAcc::Pi(contig_pi) = &contig.accs[0] else {
            panic!("pi acc");
        };
        let window_sites: u64 = output
            .windows
            .iter()
            .filter(|w| w.contig == "chr1" && !(w.start1 == 1 && w.end == 300))
            .map(|w| {
                let StatAcc::Pi(pi) = &w.accs[0] else {
                    panic!("pi acc")
                };
                pi.site_tally().0
            })
            .sum();
        assert_eq!(
            window_sites,
            contig_pi.site_tally().0,
            "window sums == contig row"
        );
        assert_eq!(
            contig_pi.site_tally().0,
            3,
            "two seam copies plus the later site"
        );
    }
}
