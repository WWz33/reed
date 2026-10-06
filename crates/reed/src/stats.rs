//! `run_stats`: the L4 statistics driver (`docs/architecture-spec.md` §8).
//! Chunk workers count genotypes straight from the VCF text into per-bin
//! additive accumulators ([`reed_stats::window::ContigBins`]); the ordered
//! sink merges chunk partials in input order, emits each window as soon as
//! the watermark proves no future site can reach it (§8.3), and folds
//! contig/genome-level totals alongside (§8.7). The TSV output is
//! byte-identical across thread counts (§7.6: chunk boundaries come from
//! the input only).
//!
//! Policy notes: multiallelic records are skipped (counted separately) —
//! the fast GT path is biallelic-only, keyed off the `ALT` column.
//! Unassigned samples still feed the site-wide variant flag but no group's
//! counts. Output levels: window rows use 1-based `[start, end]` spans;
//! the contig level spans the contig's emitted windows; the genome level
//! is a single `genome` row whose `end` is the summed span in bp (not a
//! coordinate).

use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;
use std::time::Duration;

use reed_core::count::SiteGroupCounts;
use reed_core::grid::WindowGrid;
use reed_core::ChunkId;
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{run, ChunkJob, OrderedSink};
use reed_io::chunk::{BgzfChunkSource, ChunkSource};
use reed_io::vcf::fixed::parse_pos1;
use reed_io::vcf::gt::tally_gt;
use reed_io::vcf::line::LineRef;
use reed_io::Locus;
use reed_stats::fst::FstKind;
use reed_stats::pbs::PbsKind;
use reed_stats::request::{build_accs, StatAcc, StatKey, StatKind};
use reed_stats::window::{BpBinning, ContigBins};

/// The full per-group/per-pair statistic set over `n_groups` groups: π and
/// dxy in both conventions, Fst in all three estimators, PBS in both
/// flavours (three groups required), and `SFS`/`θ_W`/`Tajima's D` per group —
/// the M4 acceptance set, in output-column order.
#[must_use]
pub fn default_kinds(n_groups: usize) -> Vec<StatKind> {
    let mut kinds = Vec::new();
    for group in 0..n_groups {
        kinds.push(StatKind::Pi { group, per_bp: false });
        kinds.push(StatKind::Pi { group, per_bp: true });
    }
    for a in 0..n_groups {
        for b in (a + 1)..n_groups {
            kinds.push(StatKind::Dxy { a, b, per_bp: false });
            kinds.push(StatKind::Dxy { a, b, per_bp: true });
            for kind in [FstKind::Hudson, FstKind::WeirCockerham84, FstKind::Nei86] {
                kinds.push(StatKind::Fst { a, b, kind });
            }
        }
    }
    if n_groups >= 3 {
        kinds.push(StatKind::Pbs { a: 0, b: 1, c: 2, kind: PbsKind::Yi2010 });
        kinds.push(StatKind::Pbs { a: 0, b: 1, c: 2, kind: PbsKind::PbsN1 });
    }
    for group in 0..n_groups {
        kinds.push(StatKind::Sfs { group, folded: false });
        kinds.push(StatKind::Sfs { group, folded: true });
        kinds.push(StatKind::ThetaW { group, per_bp: false });
        kinds.push(StatKind::ThetaW { group, per_bp: true });
        kinds.push(StatKind::TajimaD { group, min_sites: 0 });
    }
    kinds
}

/// Which output levels to emit (§8.7). Window rows always carry the same
/// columns; contig rows aggregate each contig's bins, the genome row
/// aggregates everything — all three are plain accumulator merges, so PBS
/// stays re-aggregable (its accumulator keeps per-pair num/den).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Levels {
    /// One row per window (the default).
    pub window: bool,
    /// One row per contig, emitted when the contig closes.
    pub contig: bool,
    /// One `genome` row, emitted at the end.
    pub genome: bool,
}

impl Levels {
    /// Parses a comma-separated `window,contig,genome` list.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] on an unknown level token.
    pub fn parse(text: &str) -> Result<Self, EngineError> {
        let mut levels = Self::default();
        for token in text.split(',').map(str::trim).filter(|t| !t.is_empty()) {
            match token {
                "window" => levels.window = true,
                "contig" => levels.contig = true,
                "genome" => levels.genome = true,
                other => {
                    return Err(EngineError::Job(format!(
                        "unknown level {other:?} (expected window, contig or genome)"
                    )))
                }
            }
        }
        Ok(levels)
    }
}

/// One `run_stats` request: statistics, group assignment, window grid.
#[derive(Clone, Debug)]
pub struct StatsRequest {
    /// Statistics to compute, in output-column order.
    pub kinds: Vec<StatKind>,
    /// Sample index → group index; `None` samples feed only the site-wide
    /// variant flag. Length must equal the header's sample count.
    pub gt_groups: Vec<Option<u32>>,
    /// Window width in bp.
    pub window_size: u64,
    /// Window stride in bp (`step ≤ size` for non-overlapping windows).
    pub step: u64,
    /// Output levels (default: window rows only).
    pub levels: Levels,
    /// Skip windows no site fell into (§8.5 `--skip-empty`).
    pub skip_empty: bool,
    /// Extra contig lengths (`--contig-lengths`, .fai-style), merged over
    /// the header's `##contig` declarations. Trailing-window clipping needs
    /// a length for every contig with records and is never fabricated
    /// (01-review D7).
    pub contig_lengths: HashMap<Vec<u8>, u64>,
}

/// The run's output and counters.
#[derive(Clone, Debug)]
pub struct StatsReport {
    /// TSV: a header line, then one row per selected level.
    pub tsv: String,
    /// Window-level rows emitted.
    pub n_windows: u64,
    /// Records absorbed into windows.
    pub n_records: u64,
    /// Records skipped for a multi-allelic ALT column.
    pub n_multiallelic_skipped: u64,
    /// Chunks processed.
    pub n_chunks: u64,
    /// Wall time.
    pub elapsed: Duration,
}

/// Parses a `sample<TAB>population` groups file — pixy's populations-file
/// convention: two columns, no header row, `#` comments allowed
/// (reed_core::GroupTable::from_tsv). A file row naming a sample absent
/// from the VCF header is an error (pixy behavior) unless `force_samples`,
/// which skips it with a warning. Group ids follow first-appearance order.
///
/// # Errors
///
/// [`EngineError::Job`] for an unreadable file, a malformed row, or an
/// unknown sample (without `force_samples`).
pub fn assign_groups(
    groups_path: &Path,
    sample_names: &[String],
    force_samples: bool,
) -> Result<(Vec<Option<u32>>, Vec<String>), EngineError> {
    let samples = reed_core::sample::SampleTable::from_names(sample_names.iter().map(String::as_str))
        .map_err(|err| EngineError::Job(err.to_string()))?;
    let table = reed_core::group::GroupTable::from_tsv(
        groups_path,
        &samples,
        reed_core::group::GroupFileOptions { ignore_missing_samples: force_samples },
    )
    .map_err(|err| EngineError::Job(err.to_string()))?;
    let mut gt_groups: Vec<Option<u32>> = vec![None; sample_names.len()];
    let mut names: Vec<String> = Vec::with_capacity(table.n_groups());
    for (group_idx, group) in table.iter().enumerate() {
        names.push(group.name.clone());
        let group_idx = u32::try_from(group_idx)
            .map_err(|_| EngineError::Job("more than u32::MAX groups".into()))?;
        for (word_idx, &word) in group.mask.words().iter().enumerate() {
            let mut word = word;
            while word != 0 {
                let bit = word.trailing_zeros() as usize;
                gt_groups[word_idx * 64 + bit] = Some(group_idx);
                word &= word - 1;
            }
        }
    }
    Ok((gt_groups, names))
}

/// One contig's run of sites: bin accumulators plus the span it covers.
struct SegmentPartial {
    contig: Vec<u8>,
    bins: ContigBins,
    first_pos0: u64,
    last_pos0: u64,
}

/// Worker → sink partial: per-contig bin accumulators in record order.
#[derive(Default)]
pub struct StatsPartial {
    segments: Vec<SegmentPartial>,
    n_records: u64,
    n_multiallelic_skipped: u64,
}

/// The chunk job: line → `LineRef` → GT tallies → group count row → the
/// site's bin accumulators.
pub struct StatsJob {
    width: u64,
    n_groups: usize,
    gt_groups: Arc<[Option<u32>]>,
    seed: Arc<[StatAcc]>,
    /// Declared contig lengths (header + `--contig-lengths`); positions
    /// beyond the declared length are rejected (01-review D8).
    contig_len: Arc<HashMap<Vec<u8>, u64>>,
}

impl ChunkJob for StatsJob {
    type Partial = StatsPartial;
    type Scratch = Vec<SiteGroupCounts>;

    fn run_chunk(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let mut partial = StatsPartial::default();
        let mut current: Option<SegmentPartial> = None;
        self.reset_row(scratch);
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            self.push_line(line, &mut current, &mut partial, scratch)?;
        }
        if let Some(segment) = current.take() {
            partial.segments.push(segment);
        }
        Ok(partial)
    }

    fn run_seam(
        &self,
        _id: reed_core::ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let mut partial = StatsPartial::default();
        let mut current: Option<SegmentPartial> = None;
        self.reset_row(scratch);
        self.push_line(line, &mut current, &mut partial, scratch)?;
        if let Some(segment) = current.take() {
            partial.segments.push(segment);
        }
        Ok(partial)
    }
}

impl StatsJob {
    fn reset_row(&self, row: &mut Vec<SiteGroupCounts>) {
        row.clear();
        row.resize(self.n_groups, SiteGroupCounts::default());
    }

    fn push_line(
        &self,
        line: &[u8],
        current: &mut Option<SegmentPartial>,
        partial: &mut StatsPartial,
        row: &mut [SiteGroupCounts],
    ) -> Result<(), EngineError> {
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::Job("record line has fewer than 9 columns".into()));
        };
        let chrom = lr.chrom();
        match current {
            Some(segment) if segment.contig == chrom => {}
            _ => {
                if let Some(segment) = current.take() {
                    partial.segments.push(segment);
                }
                *current = Some(SegmentPartial {
                    contig: chrom.to_vec(),
                    bins: ContigBins::new(BpBinning::new(self.width), &self.seed),
                    first_pos0: u64::MAX,
                    last_pos0: 0,
                });
            }
        }
        let locus = |pos1: u64| Locus {
            contig_name: String::from_utf8_lossy(chrom).into_owned(),
            pos1,
        };
        let pos1 =
            parse_pos1(lr.pos_text(), &locus(0)).map_err(|err| EngineError::Job(err.to_string()))?;
        let Some(pos0) = pos1.checked_sub(1) else {
            let here = locus(0);
            return Err(EngineError::Job(format!("{here}: POS is zero")));
        };
        // Input contract (04-algorithms.md §6): positions ascending within a
        // contig, never beyond the declared length.
        if let Some(&len) = self.contig_len.get(chrom) {
            if pos1 > len {
                return Err(EngineError::Job(
                    reed_core::CoreError::PositionBeyondContigLength {
                        contig_name: String::from_utf8_lossy(chrom).into_owned(),
                        pos1,
                        len_bp: len,
                    }
                    .to_string(),
                ));
            }
        }
        if memchr::memchr(b',', lr.alts()).is_some() {
            partial.n_multiallelic_skipped += 1;
            return Ok(());
        }
        for counts in row.iter_mut() {
            *counts = SiteGroupCounts::default();
        }
        let mut any_alt = false;
        for (sample_idx, sample) in lr.samples_text().split(|&b| b == b'\t').enumerate() {
            let Some(t) = tally_gt(sample) else {
                let here = locus(pos1);
                return Err(EngineError::Job(format!(
                    "{here}: GT has more than two alleles"
                )));
            };
            any_alt |= t.n_alt > 0;
            if let Some(group) = self.gt_groups.get(sample_idx).copied().flatten() {
                let counts = &mut row[group as usize];
                counts.n_ref_alleles += t.n_ref;
                counts.n_alt_alleles += t.n_alt;
                counts.n_hom_ref += t.n_hom_ref;
                counts.n_het += t.n_het;
                counts.n_hom_alt += t.n_hom_alt;
                counts.n_half_called += t.n_half;
            }
        }
        let segment = current
            .as_mut()
            .expect("invariant: segment opened for this contig above");
        if pos0 < segment.last_pos0 {
            return Err(EngineError::Job(
                reed_core::CoreError::PositionNotAscending {
                    prev_pos0: segment.last_pos0,
                    pos0,
                }
                .to_string(),
            ));
        }
        segment.bins.push_site(pos0, row, any_alt);
        segment.first_pos0 = segment.first_pos0.min(pos0);
        segment.last_pos0 = pos0;
        partial.n_records += 1;
        Ok(())
    }
}

/// Everything the sink hands back through the shared slot at `finish`.
struct SinkOutput {
    rows: String,
    n_windows: u64,
    n_records: u64,
    n_multiallelic_skipped: u64,
}

/// The sink's per-contig accumulation state.
struct CurrentContig {
    name: Vec<u8>,
    bins: ContigBins,
    /// Contig-level totals (§8.7): every touched bin folded in.
    totals: Vec<StatAcc>,
    first_pos0: u64,
    last_pos0: u64,
    /// First window start (a multiple of the stride, legacy anchoring §5.9).
    anchor: u64,
    /// Start of the next window to emit (the watermark frontier).
    next_start: u64,
    /// End of the last emitted window (0-based).
    span_end: u64,
}

/// Ordered consumer: merges same-contig segments, emits complete windows
/// through the watermark (§8.3), and closes contigs with a contig-level
/// row (§8.7). Revisited contigs are refused — unsorted input would
/// otherwise silently double-count windows.
struct StatsSink {
    size: u64,
    step: u64,
    width: u64,
    kinds: Arc<[StatKind]>,
    seed: Vec<StatAcc>,
    levels: Levels,
    skip_empty: bool,
    contig_len: HashMap<Vec<u8>, u64>,
    current: Option<CurrentContig>,
    closed: Vec<Vec<u8>>,
    genome_totals: Vec<StatAcc>,
    genome_span: u64,
    rows: String,
    n_windows: u64,
    n_records: u64,
    n_multiallelic_skipped: u64,
    out: Arc<std::sync::Mutex<Option<SinkOutput>>>,
}

impl OrderedSink<StatsPartial> for StatsSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut StatsPartial) -> Result<(), EngineError> {
        self.n_records += partial.n_records;
        self.n_multiallelic_skipped += partial.n_multiallelic_skipped;
        for mut segment in std::mem::take(&mut partial.segments) {
            let mut current = match self.current.take() {
                Some(mut cur) if cur.name == segment.contig => {
                    // Cross-chunk seam (04-algorithms.md §6): a later chunk
                    // may not rewind into positions this sink already
                    // emitted windows for.
                    if segment.first_pos0 < cur.last_pos0 {
                        return Err(EngineError::Job(
                            reed_core::CoreError::PositionNotAscending {
                                prev_pos0: cur.last_pos0,
                                pos0: segment.first_pos0,
                            }
                            .to_string(),
                        ));
                    }
                    segment.bins.fold(&mut cur.totals);
                    cur.bins.merge_bins(&mut segment.bins);
                    cur.first_pos0 = cur.first_pos0.min(segment.first_pos0);
                    cur.last_pos0 = cur.last_pos0.max(segment.last_pos0);
                    cur
                }
                Some(cur) => {
                    self.close_contig(cur)?;
                    Self::open_and_absorb(&mut segment, self.step, self.width, &self.seed)
                }
                None => Self::open_and_absorb(&mut segment, self.step, self.width, &self.seed),
            };
            self.emit_complete_windows(&mut current)?;
            self.current = Some(current);
        }
        Ok(())
    }

    fn finish(mut self) -> Result<u64, EngineError> {
        if let Some(current) = self.current.take() {
            self.close_contig(current)?;
        }
        if self.levels.genome {
            let bp = self.genome_span;
            let totals = self.genome_totals.clone();
            self.write_row("genome", 0, bp, bp, &totals);
        }
        *self
            .out
            .lock()
            .expect("invariant: stats slot not poisoned") = Some(SinkOutput {
            rows: self.rows,
            n_windows: self.n_windows,
            n_records: self.n_records,
            n_multiallelic_skipped: self.n_multiallelic_skipped,
        });
        Ok(self.n_windows)
    }
}

impl StatsSink {
    /// Opens a contig's state and absorbs its first segment (folding the
    /// segment's bins into the contig totals and the bin map).
    fn open_and_absorb(
        segment: &mut SegmentPartial,
        step: u64,
        width: u64,
        seed: &[StatAcc],
    ) -> CurrentContig {
        let anchor = (segment.first_pos0 / step) * step;
        let mut cur = CurrentContig {
            name: segment.contig.clone(),
            bins: ContigBins::new(BpBinning::new(width), seed),
            totals: seed.to_vec(),
            first_pos0: segment.first_pos0,
            last_pos0: segment.last_pos0,
            anchor,
            next_start: anchor,
            span_end: 0,
        };
        segment.bins.fold(&mut cur.totals);
        cur.bins.merge_bins(&mut segment.bins);
        cur
    }

    /// Emits every window that the watermark proves complete: with sorted
    /// input, no future site can land in `[s, s+size)` once the last seen
    /// position reaches `s+size−1` (§8.3). Consumed bins are freed. The
    /// contig length does not gate this loop — only the trailing windows
    /// in `close_contig` need it.
    fn emit_complete_windows(&mut self, cur: &mut CurrentContig) -> Result<(), EngineError> {
        let watermark = cur.last_pos0.saturating_add(1);
        let mut s = cur.next_start;
        while s.saturating_add(self.size) <= watermark {
            let end = s + self.size;
            let bins = (s / self.width)..end.div_ceil(self.width);
            if !(self.skip_empty && !cur.bins.any_bin_in(bins.clone())) {
                let accs = cur.bins.finish_bins(bins);
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
            s = s.saturating_add(self.step);
        }
        cur.next_start = s;
        cur.bins.drop_bins_below(s / self.width);
        Ok(())
    }

    /// Closes a contig: emits the trailing windows (start at or before the
    /// last site, clipped to the contig length — the contig length is
    /// required here and is never fabricated, 01-review D7), then the
    /// contig-level row (§8.7), then folds into the genome total.
    fn close_contig(&mut self, mut cur: CurrentContig) -> Result<(), EngineError> {
        if self.closed.iter().any(|closed| closed == &cur.name) {
            return Err(EngineError::Job(format!(
                "contig {} revisited after a later contig: input is not contig-sorted",
                String::from_utf8_lossy(&cur.name)
            )));
        }
        self.closed.push(cur.name.clone());
        let len = self.contig_len.get(&cur.name).copied().ok_or_else(|| {
            EngineError::Job(format!(
                "contig {} has no ##contig length; window clipping needs it (pass --contig-lengths)",
                String::from_utf8_lossy(&cur.name)
            ))
        })?;
        let contig = String::from_utf8_lossy(&cur.name).into_owned();
        let mut s = cur.next_start;
        while s < len && s <= cur.last_pos0 {
            let end = s.saturating_add(self.size).min(len);
            let bins = (s / self.width)..end.div_ceil(self.width);
            if !(self.skip_empty && !cur.bins.any_bin_in(bins.clone())) {
                let accs = cur.bins.finish_bins(bins);
                self.write_row(&contig, s + 1, end, end - s, &accs);
                cur.span_end = end;
                self.n_windows += 1;
            }
            s = s.saturating_add(self.step);
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

    /// Writes one TSV row at any level; `bp` is the span's base-pair count
    /// (the AccessibleBp/per-bp denominator at that level).
    fn write_row(&mut self, contig: &str, start1: u64, end0: u64, bp: u64, accs: &[StatAcc]) {
        self.rows.push_str(contig);
        self.rows.push('\t');
        self.rows.push_str(&start1.to_string());
        self.rows.push('\t');
        self.rows.push_str(&end0.to_string());
        for (kind, acc) in self.kinds.iter().zip(accs) {
            self.rows.push('\t');
            match acc.sfs_bins() {
                Some(bins) => {
                    for (i, bin) in bins.iter().enumerate() {
                        if i > 0 {
                            self.rows.push(';');
                        }
                        self.rows.push_str(&bin.to_string());
                    }
                }
                None => match acc.value(kind, bp) {
                    Some(value) => self.rows.push_str(&format!("{value}")),
                    None => self.rows.push('.'),
                },
            }
        }
        self.rows.push('\n');
    }
}

/// Runs windowed statistics over a BGZF VCF (`docs/architecture-spec.md`
/// milestone M4: the full 3-population statistic set).
///
/// # Errors
///
/// Transport, job, and request failures — including a sample-count mismatch
/// against the header, an invalid grid, or a group index without samples.
pub fn run_stats(
    path: Option<&Path>,
    request: StatsRequest,
    budget: &ThreadBudget,
) -> Result<StatsReport, EngineError> {
    let grid = WindowGrid::Bp { size: request.window_size, step: request.step };
    grid.validate().map_err(|err| EngineError::Job(err.to_string()))?;
    let binning = BpBinning::from_grid(&grid)
        .ok_or_else(|| EngineError::Job("stats driver requires a Bp grid".into()))?;

    let n_groups = request
        .gt_groups
        .iter()
        .flatten()
        .copied()
        .max()
        .map_or(0, |g| g as usize + 1);
    let mut group_sizes = vec![0_u32; n_groups];
    for group in request.gt_groups.iter().flatten() {
        group_sizes[*group as usize] += 1;
    }
    let accs = build_accs(&request.kinds, &group_sizes)
        .map_err(|err| EngineError::Job(format!("invalid stats request: {err}")))?;
    let keys: Vec<StatKey> = request.kinds.iter().map(StatKind::key).collect();

    let (mut source, header) = crate::input::open(path).map_err(EngineError::from)?;
    if header.n_samples() != request.gt_groups.len() {
        return Err(EngineError::Job(format!(
            "request assigns {} samples but the header declares {}",
            request.gt_groups.len(),
            header.n_samples()
        )));
    }
    let contig_len: Arc<HashMap<Vec<u8>, u64>> = Arc::new({
        let mut map: HashMap<Vec<u8>, u64> = request.contig_lengths;
        for (name, len) in header
            .contig_names
            .iter()
            .zip(header.contig_lengths.iter().copied())
        {
            if let Some(len) = len {
                map.insert(name.as_bytes().to_vec(), len);
            }
        }
        map
    });

    let job = StatsJob {
        width: binning.width(),
        n_groups,
        gt_groups: request.gt_groups.into(),
        seed: accs.clone().into(),
        contig_len: Arc::clone(&contig_len),
    };
    let out = Arc::new(std::sync::Mutex::new(None));
    let sink = StatsSink {
        size: request.window_size,
        step: request.step,
        width: binning.width(),
        kinds: Arc::from(request.kinds.into_boxed_slice()),
        seed: accs.clone(),
        levels: request.levels,
        skip_empty: request.skip_empty,
        contig_len: (*contig_len).clone(),
        current: None,
        closed: Vec::new(),
        genome_totals: accs,
        genome_span: 0,
        rows: String::new(),
        n_windows: 0,
        n_records: 0,
        n_multiallelic_skipped: 0,
        out: Arc::clone(&out),
    };
    let (_abort_handle, abort_signal) = reed_engine::abort::abort_channel();
    let report = run(source.as_mut(), &job, sink, budget, &abort_signal)?;
    let output = out
        .lock()
        .expect("invariant: stats slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the output");

    let mut tsv = String::from("contig\tstart\tend");
    for key in &keys {
        tsv.push('\t');
        tsv.push_str(key.as_str());
    }
    tsv.push('\n');
    tsv.push_str(&output.rows);
    Ok(StatsReport {
        tsv,
        n_windows: output.n_windows,
        n_records: output.n_records,
        n_multiallelic_skipped: output.n_multiallelic_skipped,
        n_chunks: report.n_chunks,
        elapsed: report.elapsed,
    })
}
