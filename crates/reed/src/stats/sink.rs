//! The ordered sink of the statistics run: merges chunk partials in input order,
//! emits window rows through the watermark, and closes contigs and the genome.

use crate::stats::request::{
    FStatsPartial, InvariantMode, LdPlan, Levels, SfsPartial, WindowLayoutOptions, WindowRow,
};
use crate::stats::run::{ContigAcc, GridMode, SegmentPartial, StatsPartial, decay_contig};
use reed_core::ChunkId;
use reed_core::grid::RegionSet;
use reed_engine::error::EngineError;
use reed_engine::pipeline::OrderedSink;
use reed_stats::ld::LdDecayAcc;
use reed_stats::request::StatAcc;
use reed_stats::window::{BpBinning, ContigBins, RegionAccs};
use std::collections::HashMap;
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
    pub(super) levels: Levels,
    pub(super) skip_empty: bool,
    /// Window-layout mode (structured rows + invariant policy).
    pub(super) layout: Option<WindowLayoutOptions>,
    /// Structured window rows; the flat table renders from these
    /// (`writers::windows::render_tsv`).
    pub(super) windows: Vec<WindowRow>,
    /// Site-wide invariant sites seen.
    pub(super) n_invariant_sites: u64,
    /// Records with ALT = `.` (the all-sites check).
    pub(super) n_no_alt_records: u64,
    pub(super) contig_len: HashMap<Vec<u8>, u64>,
    pub(super) current: Option<CurrentContig>,
    pub(super) closed: Vec<Vec<u8>>,
    pub(super) genome_totals: Vec<StatAcc>,
    /// Merged genome-wide projected spectra.
    pub(super) sfs: Option<SfsPartial>,
    /// Merged genome-wide f-statistics.
    pub(super) fstats: Option<FStatsPartial>,
    /// The LD plane layout (group-major); present only with an LD request.
    pub(super) ld_layout: Option<Arc<reed_stats::SampleLayout>>,
    /// The matrix being accumulated for the pending contig.
    pub(super) ld_matrix: Option<reed_stats::ContigMatrix>,
    /// Per-contig decay jobs in contig order (joined at `finish`, so the
    /// analysis overlaps the decode instead of stalling the reducer).
    pub(super) ld_jobs: Vec<std::thread::JoinHandle<LdResult>>,
    /// The LD plan (the analysis needs the window and the bin edges).
    pub(super) ld: Option<LdPlan>,
    pub(super) genome_span: u64,
    pub(super) n_windows: u64,
    pub(super) n_records: u64,
    pub(super) n_multiallelic_skipped: u64,
    pub(super) out: Arc<std::sync::Mutex<Option<SinkOutput>>>,
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
                        self.finalize_ld();
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
                    // may not rewind into positions this sink already
                    // emitted windows for.
                    if segment.first_pos0 < cur.last_pos0 {
                        return Err(EngineError::from(
                            reed_core::CoreError::PositionNotAscending {
                                prev_pos0: cur.last_pos0,
                                pos0: segment.first_pos0,
                            },
                        ));
                    }
                    segment.acc.fold(&mut cur.totals);
                    cur.acc.merge(&mut segment.acc);
                    cur.first_pos0 = cur.first_pos0.min(segment.first_pos0);
                    cur.last_pos0 = cur.last_pos0.max(segment.last_pos0);
                    cur
                }
                Some(cur) => {
                    self.close_contig(cur)?;
                    Self::open_and_absorb(&mut segment, &self.grid, &self.seed)
                }
                None => Self::open_and_absorb(&mut segment, &self.grid, &self.seed),
            };
            self.emit_complete_windows(&mut current)?;
            self.current = Some(current);
        }
        Ok(())
    }

    fn finish(mut self) -> Result<u64, EngineError> {
        self.finalize_ld();
        let mut ld = None;
        self.join_ld(&mut ld)?;
        if let Some(current) = self.current.take() {
            self.close_contig(current)?;
        }
        if self.levels.genome {
            let bp = self.genome_span;
            let totals = self.genome_totals.clone();
            self.write_row("genome", 0, bp, bp, &totals);
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
            return Err(EngineError::job(format!(
                "input has no invariant sites ({} records, none with ALT = \".\"): declare \
                 --variant-only for an all-variant file, or --bypass-invariant-check to compute \
                 over present sites only (both overestimate π/dxy denominators)",
                self.n_records
            )));
        }
        *self.out.lock().expect("invariant: stats slot not poisoned") = Some(SinkOutput {
            ld,
            windows: std::mem::take(&mut self.windows),
            sfs: self.sfs.take(),
            fstats: self.fstats.take(),
            n_invariant_sites: self.n_invariant_sites,
            n_no_alt_records: self.n_no_alt_records,
            n_windows: self.n_windows,
            n_records: self.n_records,
            n_multiallelic_skipped: self.n_multiallelic_skipped,
        });
        Ok(self.n_windows)
    }
}

impl StatsSink {
    /// Hands the pending contig's matrix to a decay worker: the analysis
    /// runs off the reducer thread (it overlaps the decode), and `finish`
    /// joins the jobs in contig order — the merge order the deterministic
    /// output requires (docs/rounds/v6/02-models.md §2.2).
    pub(super) fn finalize_ld(&mut self) {
        let Some(matrix) = self.ld_matrix.take() else {
            return;
        };
        let Some(plan) = &self.ld else { return };
        let edges = plan.decay_edges();
        let max_bp = plan.decay_max_bp.unwrap_or(100_000);
        let collect = plan.pairs;
        self.ld_jobs.push(std::thread::spawn(move || {
            let set = decay_contig(&matrix, max_bp, &edges, collect);
            LdResult {
                cohort: set.cohort,
                groups: set.groups,
                pairs: set
                    .pairs
                    .into_iter()
                    .map(|(a, b, r2)| (matrix.contig().to_vec(), a, b, r2))
                    .collect(),
            }
        }));
    }

    /// Joins the per-contig decay jobs in contig order.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] when a decay worker panicked.
    pub(super) fn join_ld(&mut self, out: &mut Option<LdResult>) -> Result<(), EngineError> {
        for handle in self.ld_jobs.drain(..) {
            let set = handle
                .join()
                .map_err(|_| EngineError::job("the LD decay worker panicked"))?;
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
                GridMode::Regions(regions) => ContigAcc::Regions(RegionAccs::new(
                    regions
                        .get(&segment.contig)
                        .unwrap_or_else(|| Arc::new(RegionSet::default())),
                    seed,
                )),
            },
            totals: seed.to_vec(),
            first_pos0: segment.first_pos0,
            last_pos0: segment.last_pos0,
            anchor,
            next_start: anchor,
            span_end: 0,
        };
        segment.acc.fold(&mut cur.totals);
        cur.acc.merge(&mut segment.acc);
        cur
    }

    /// Emits every window that the watermark proves complete: with sorted
    /// input, no future site can land in `[s, s+size)` once the last seen
    /// position reaches `s+size−1` (§8.3). Consumed bins are freed. The
    /// contig length does not gate this loop — only the trailing windows
    /// in `close_contig` need it. Regions grids emit nothing here: their
    /// windows are the BED intervals, known up front.
    #[allow(clippy::unnecessary_wraps)] // callers share the Result-returning shape of the window writers
    pub(super) fn emit_complete_windows(
        &mut self,
        cur: &mut CurrentContig,
    ) -> Result<(), EngineError> {
        let GridMode::Bp { size, step, width } = self.grid else {
            return Ok(());
        };
        let watermark = cur.last_pos0.saturating_add(1);
        let mut s = cur.next_start;
        let ContigAcc::Bins(bins) = &mut cur.acc else {
            return Ok(());
        };
        while s.saturating_add(size) <= watermark {
            let end = s + size;
            let range = (s / width)..end.div_ceil(width);
            let skip = self.skip_empty && !bins.any_bin_in(range.clone());
            if !skip {
                let accs = bins.finish_bins(range);
                self.write_row(
                    &String::from_utf8_lossy(&cur.name),
                    s + 1,
                    end,
                    end - s,
                    &accs,
                );
                cur.span_end = end;
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
        if self.closed.iter().any(|closed| closed == &cur.name) {
            return Err(EngineError::job(format!(
                "contig {} revisited after a later contig: input is not contig-sorted",
                String::from_utf8_lossy(&cur.name)
            )));
        }
        self.closed.push(cur.name.clone());
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
                    let end = s.saturating_add(size).min(len);
                    let range = (s / width)..end.div_ceil(width);
                    let skip = self.skip_empty && !bins.any_bin_in(range.clone());
                    if !skip {
                        let accs = bins.finish_bins(range);
                        self.write_row(&contig, s + 1, end, end - s, &accs);
                        cur.span_end = end;
                        self.n_windows += 1;
                    }
                    s = s.saturating_add(step);
                }
            }
            GridMode::Regions(_) => {
                // The windows are the BED intervals in input order; every
                // region is emitted (skip-empty drops untouched ones).
                let ContigAcc::Regions(regions) = &cur.acc else {
                    unreachable!("invariant: Regions grid implies region accumulators");
                };
                let intervals: Vec<reed_core::coord::Interval0> =
                    regions.regions().windows().to_vec();
                for (idx, interval) in intervals.iter().enumerate() {
                    if self.skip_empty && !regions.touched(idx) {
                        continue;
                    }
                    let start0 = interval.start.0;
                    let end0 = interval.end.0;
                    let bp = end0.saturating_sub(start0);
                    let accs = regions.accs(idx);
                    self.write_row(&contig, start0 + 1, end0, bp, accs);
                    cur.span_end = cur.span_end.max(end0);
                    self.n_windows += 1;
                }
            }
        }
        if self.levels.contig {
            let bp = cur.span_end.saturating_sub(cur.anchor);
            self.write_row(&contig, cur.anchor + 1, cur.span_end, bp, &cur.totals);
        }
        for (dst, src) in self.genome_totals.iter_mut().zip(&cur.totals) {
            dst.merge(src);
        }
        self.genome_span += cur.span_end.saturating_sub(cur.anchor);
        Ok(())
    }

    /// Records one row at any level; `bp` is the span's base-pair count
    /// (the AccessibleBp/per-bp denominator at that level). The text forms
    /// render from these rows (`writers::windows`).
    pub(super) fn write_row(
        &mut self,
        contig: &str,
        start1: u64,
        end0: u64,
        bp: u64,
        accs: &[StatAcc],
    ) {
        self.windows.push(WindowRow {
            contig: contig.to_string(),
            start1,
            end0,
            bp,
            accs: accs.to_vec(),
        });
    }
}
