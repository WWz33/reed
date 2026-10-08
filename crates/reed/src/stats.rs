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
use reed_core::count::code_has_alt;
use reed_core::genotype::{GenotypePolicies, MultiallelicPolicy, PloidyPolicy, PolicyCounters};
use reed_io::vcf::genotype_scan::GenotypeScanner;
use reed_io::vcf::site::SiteView;
use reed_io::vcf::line::LineRef;
pub use reed_stats::fst::FstKind;
pub use reed_stats::pbs::PbsKind;
use reed_stats::request::{build_accs, StatAcc, StatKey};
pub use reed_stats::request::StatKind;
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
            for kind in [FstKind::Hudson, FstKind::WeirCockerham84] {
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

/// The pixy-mode statistic set: π per group and/or dxy per pair, in output
/// column order (02-cli.md §5.2).
#[must_use]
pub fn pixy_kinds(
    n_groups: usize,
    pi: bool,
    dxy: bool,
    theta: bool,
    tajima: bool,
    fst: Option<FstKind>,
    pbs: Option<PbsKind>,
) -> Vec<StatKind> {
    let mut kinds = Vec::new();
    if pi {
        for group in 0..n_groups {
            kinds.push(StatKind::Pi { group, per_bp: false });
        }
    }
    if dxy {
        for a in 0..n_groups {
            for b in (a + 1)..n_groups {
                kinds.push(StatKind::Dxy { a, b, per_bp: false });
            }
        }
    }
    if theta {
        for group in 0..n_groups {
            kinds.push(StatKind::ThetaW { group, per_bp: false });
        }
    }
    if tajima {
        for group in 0..n_groups {
            kinds.push(StatKind::TajimaD { group, min_sites: 0 });
        }
    }
    if let Some(kind) = fst {
        for a in 0..n_groups {
            for b in (a + 1)..n_groups {
                kinds.push(StatKind::Fst { a, b, kind });
            }
        }
    }
    if let Some(kind) = pbs {
        // One triple: group 0 is the branch of interest (Yi et al. 2010).
        if n_groups >= 3 {
            kinds.push(StatKind::Pbs { a: 0, b: 1, c: 2, kind });
        }
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

/// How the input's invariant sites are treated (04-algorithms.md §3,
/// 02-cli.md §5.3). pixy's rule: an all-sites input is the default, and
/// anything else must be declared.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InvariantMode {
    /// The input must contain invariant sites (pixy's check); a file with
    /// only variant sites fails at the end of the run.
    Require,
    /// `--variant-only`: the file has only variant sites; absent window
    /// positions count as monomorphic with all samples called, so π/dxy
    /// denominators are completed by `(L − |S|) · pairs`.
    VariantOnly,
    /// `--bypass-invariant-check`: compute over the sites the file contains
    /// (overestimates π and dxy; the caller warns).
    Bypass,
}

/// pixy-layout window output (02-cli.md §5.4): π/dxy files with the
/// integer audit columns.
#[derive(Clone, Copy, Debug)]
pub struct PixyOptions {
    /// Invariant-site policy.
    pub invariant: InvariantMode,
}

/// One window's row, structured for the writers.
#[derive(Clone, Debug)]
pub struct WindowRow {
    /// Contig name.
    pub contig: String,
    /// 1-based inclusive start.
    pub start1: u64,
    /// 1-based inclusive end.
    pub end0: u64,
    /// Window span in bp (the per-bp denominator).
    pub bp: u64,
    /// Accumulators in `kinds` order.
    pub accs: Vec<StatAcc>,
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
    /// `Some` = pixy-layout output mode (structured rows, invariant-site
    /// policy); `None` = the wide debug TSV.
    pub pixy: Option<PixyOptions>,
}

/// The run's output and counters.
#[derive(Clone, Debug)]
pub struct StatsReport {
    /// TSV: a header line, then one row per selected level.
    pub tsv: String,
    /// Structured window rows (pixy mode only).
    pub windows: Vec<WindowRow>,
    /// Site-wide invariant sites seen in the input.
    pub n_invariant_sites: u64,
    /// Records with ALT = `.` (pixy's all-sites check).
    pub n_no_alt_records: u64,
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
    /// Site-wide invariant sites (no ALT called by any sample).
    n_invariant_sites: u64,
    /// Records whose ALT column is `.` (pixy's all-sites marker).
    n_no_alt_records: u64,
}

/// The chunk job: line → `LineRef` → GT tallies → group count row → the
/// site's bin accumulators.
pub struct StatsJob {
    width: u64,
    n_groups: usize,
    gt_groups: Arc<[Option<u32>]>,
    n_samples: usize,
    seed: Arc<[StatAcc]>,
    /// Declared contig lengths (header + `--contig-lengths`); positions
    /// beyond the declared length are rejected (01-review D8).
    contig_len: Arc<HashMap<Vec<u8>, u64>>,
}

/// Worker scratch: the group count row, the shared record decoder and its
/// per-record policy counters.
#[derive(Default)]
pub struct StatsScratch {
    row: Vec<SiteGroupCounts>,
    scanner: Option<GenotypeScanner>,
    counters: PolicyCounters,
}

/// The statistics decode policies: haploid calls are half calls and ALT
/// indices above 1 collapse onto the single ALT allele (pixy's biallelic
/// view; multi-allelic records are skipped before the scan anyway).
fn stats_policies() -> GenotypePolicies {
    GenotypePolicies {
        ploidy: PloidyPolicy::HaploidAsHalfCall,
        multiallelic: MultiallelicPolicy::CollapseAlt,
        ..GenotypePolicies::default()
    }
}

impl ChunkJob for StatsJob {
    type Partial = StatsPartial;
    type Scratch = StatsScratch;

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
    fn reset_row(&self, scratch: &mut StatsScratch) {
        scratch.row.clear();
        scratch.row.resize(self.n_groups, SiteGroupCounts::default());
    }

    fn push_line(
        &self,
        line: &[u8],
        current: &mut Option<SegmentPartial>,
        partial: &mut StatsPartial,
        scratch: &mut StatsScratch,
    ) -> Result<(), EngineError> {
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::Job("record line has fewer than 9 columns".into()));
        };
        let site = SiteView::parse(&lr).map_err(|err| EngineError::Job(err.to_string()))?;
        let chrom = site.contig_name;
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
        let pos1 = site.pos1;
        let Some(pos0) = pos1.checked_sub(1) else {
            return Err(EngineError::Job(format!("{}: POS is zero", site.locus())));
        };
        // Input contract (04-algorithms.md §6): positions ascending within a
        // contig, never beyond the declared length.
        if let Some(&len) = self.contig_len.get(chrom)
            && pos1 > len
        {
            return Err(EngineError::Job(
                reed_core::CoreError::PositionBeyondContigLength {
                    contig_name: String::from_utf8_lossy(chrom).into_owned(),
                    pos1,
                    len_bp: len,
                }
                .to_string(),
            ));
        }
        if site.is_invariant() {
            partial.n_no_alt_records += 1;
        }
        if site.n_alt() > 1 {
            partial.n_multiallelic_skipped += 1;
            return Ok(());
        }
        let StatsScratch { row, scanner, counters } = scratch;
        for counts in row.iter_mut() {
            *counts = SiteGroupCounts::default();
        }
        *counters = PolicyCounters::default();
        let n_samples = self.n_samples;
        let scanner = scanner
            .get_or_insert_with(|| GenotypeScanner::new(n_samples, &[], &[], stats_policies()));
        let out = scanner
            .scan(&lr, &site, counters)
            .map_err(|err| EngineError::Job(err.to_string()))?;
        let mut any_alt = false;
        for (sample_idx, &code) in out.codes.iter().enumerate() {
            any_alt |= code_has_alt(code);
            if let Some(group) = self.gt_groups.get(sample_idx).copied().flatten() {
                row[group as usize].add_code(code);
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
        if !any_alt {
            partial.n_invariant_sites += 1;
        }
        partial.n_records += 1;
        Ok(())
    }
}

/// Everything the sink hands back through the shared slot at `finish`.
struct SinkOutput {
    rows: String,
    windows: Vec<WindowRow>,
    n_invariant_sites: u64,
    n_no_alt_records: u64,
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
    /// Grid origin: windows start at 0 (1-based position 1), the
    /// pixy/vcftools grid (04-algorithms.md §4).
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
    /// pixy-layout mode (structured rows + invariant policy).
    pixy: Option<PixyOptions>,
    /// Structured rows (pixy mode only).
    windows: Vec<WindowRow>,
    /// Site-wide invariant sites seen.
    n_invariant_sites: u64,
    /// Records with ALT = `.` (pixy's all-sites check).
    n_no_alt_records: u64,
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
        self.n_invariant_sites += partial.n_invariant_sites;
        self.n_no_alt_records += partial.n_no_alt_records;
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
        // pixy's all-sites check (04-algorithms.md §3): a file with no
        // invariant sites must be declared with --variant-only or
        // --bypass-invariant-check.
        if matches!(self.pixy, Some(PixyOptions { invariant: InvariantMode::Require }))
            && self.n_no_alt_records == 0
        {
            return Err(EngineError::Job(format!(
                "input has no invariant sites ({} records, none with ALT = \".\"): declare \
                 --variant-only for an all-variant file, or --bypass-invariant-check to compute \
                 over present sites only (both overestimate π/dxy denominators)",
                self.n_records
            )));
        }
        *self
            .out
            .lock()
            .expect("invariant: stats slot not poisoned") = Some(SinkOutput {
            rows: self.rows,
            windows: std::mem::take(&mut self.windows),
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
    /// Opens a contig's state and absorbs its first segment (folding the
    /// segment's bins into the contig totals and the bin map).
    fn open_and_absorb(
        segment: &mut SegmentPartial,
        step: u64,
        width: u64,
        seed: &[StatAcc],
    ) -> CurrentContig {
        // Grid from position 1: `[1, w]`, `[1+s, w+s]`, … (04 §4).
        let anchor = 0;
        let _ = step;
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
        let variant_only = matches!(
            self.pixy,
            Some(PixyOptions { invariant: InvariantMode::VariantOnly })
        );
        if self.pixy.is_some() {
            self.windows.push(WindowRow {
                contig: contig.to_string(),
                start1,
                end0,
                bp,
                accs: accs.to_vec(),
            });
        }
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
                None => {
                    let value = if variant_only {
                        acc.value_variant_only(kind, bp)
                    } else {
                        acc.value(kind, bp)
                    };
                    match value {
                        Some(value) => self.rows.push_str(&format!("{value}")),
                        None => self.rows.push('.'),
                    }
                }
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
    let opened = crate::input::open(path).map_err(EngineError::from)?;
    run_stats_opened(opened, request, budget)
}

/// Runs windowed statistics over an already-opened input (the CLI opens
/// once so stdin streams; doc_v3 03-performance.md O5).
///
/// # Errors
///
/// As [`run_stats`].
pub fn run_stats_opened(
    opened: crate::input::OpenedInput,
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

    let abort = reed_engine::abort::abort();
    let header = opened.facts.clone();
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
        gt_groups: request.gt_groups.clone().into(),
        n_samples: request.gt_groups.len(),
        seed: accs.clone().into(),
        contig_len: Arc::clone(&contig_len),
    };
    let job = crate::bcf_input::InputJob::new(job, &opened, &abort);
    let crate::input::OpenedInput { mut source, .. } = opened;
    let out = Arc::new(std::sync::Mutex::new(None));
    let sink = StatsSink {
        size: request.window_size,
        step: request.step,
        width: binning.width(),
        kinds: Arc::from(request.kinds.into_boxed_slice()),
        seed: accs.clone(),
        levels: request.levels,
        skip_empty: request.skip_empty,
        pixy: request.pixy,
        windows: Vec::new(),
        n_invariant_sites: 0,
        n_no_alt_records: 0,
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
    let report = run(source.as_mut(), &job, sink, budget, &abort)?;
    job.finish(report.n_chunks, true)?;
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
        windows: output.windows,
        n_invariant_sites: output.n_invariant_sites,
        n_no_alt_records: output.n_no_alt_records,
        n_windows: output.n_windows,
        n_records: output.n_records,
        n_multiallelic_skipped: output.n_multiallelic_skipped,
        n_chunks: report.n_chunks,
        elapsed: report.elapsed,
    })
}
