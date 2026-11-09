//! The statistics run driver: chunk job (`StatsJob`), per-chunk partials, the
//! grid/contig accumulation state, and `run_popgen` / `run_popgen_opened`.

use crate::stats::request::{FStatsPartial, FStatsPlan, LdPartial, LdRun, SfsPartial, SfsPlan};
use reed_core::count::SiteGroupCounts;
use reed_core::count::code_has_alt;
use reed_core::genotype::{
    GenotypePolicies, GtCode, MultiallelicPolicy, PloidyPolicy, PolicyCounters,
};
use reed_core::grid::{ContigRegions, RegionSet};
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::ChunkJob;
use reed_io::bcf::BcfDict;
use reed_io::vcf::genotype_scan::GenotypeScanner;
use reed_io::vcf::site::SiteView;
use reed_stats::request::StatAcc;
use reed_stats::sfs::{ProjectionCache, WeightsArc};
use reed_stats::window::{BpBinning, ContigBins, RegionAccs};
use smallvec::SmallVec;
use std::collections::HashMap;
use std::sync::Arc;

/// A contig's shared BED accumulator. Partial segments carry this `Arc`;
/// the sink merges it once, rather than once per chunk (doc_v27 V27-19).
type SharedRegionAcc = Arc<std::sync::Mutex<RegionAccs>>;
type SharedRegionMap = Arc<std::sync::Mutex<HashMap<Vec<u8>, SharedRegionAcc>>>;

/// One contig's run of sites: bin accumulators plus the span it covers.
#[derive(Clone)]
pub(super) struct SegmentPartial {
    pub(super) contig: Vec<u8>,
    pub(super) acc: ContigAcc,
    pub(super) first_pos0: u64,
    pub(super) last_pos0: u64,
}

/// One raw site in the streaming Bp path. The row is small compared with a
/// full `StatAcc` and is consumed by the sink before the next chunk arrives.
#[derive(Clone, Debug)]
pub(super) struct SiteRow {
    pub(super) pos0: u64,
    pub(super) counts: SmallVec<[SiteGroupCounts; 4]>,
    pub(super) all_variant: bool,
}

/// The per-contig accumulation state, one variant per grid kind: `Bp` bins
/// (watermark window emission) or per-region accumulators (`--bed-file`).
#[derive(Clone, Debug)]
pub(super) enum ContigAcc {
    Bins(ContigBins),
    /// Raw per-site rows for the CLI's streaming window path. Keeping one
    /// `StatAcc` per site is what made `--window 1` allocate hundreds of MiB;
    /// rows are merged into the sink's reusable bins instead (V27-18).
    Rows(Vec<SiteRow>),
    Regions(SharedRegionAcc),
    /// No window bookkeeping at all: the SFS-only composition the convert
    /// driver builds (`StatsJob::for_sfs`, doc_v13 §7.5 plan A).
    None,
}

impl ContigAcc {
    pub(super) fn push_site(&mut self, pos0: u64, row: &[SiteGroupCounts], all_variant: bool) {
        match self {
            Self::Bins(bins) => bins.push_site(pos0, row, all_variant),
            Self::Rows(rows) => rows.push(SiteRow {
                pos0,
                counts: SmallVec::from_slice(row),
                all_variant,
            }),
            Self::Regions(regions) => regions
                .lock()
                .expect("invariant: region accumulator not poisoned")
                .push_site(pos0, row, all_variant),
            Self::None => {}
        }
    }

    pub(super) fn fold(&self, totals: &mut [StatAcc]) {
        match self {
            Self::Bins(bins) => bins.fold(totals),
            Self::Rows(rows) => {
                for row in rows {
                    for acc in totals.iter_mut() {
                        acc.push(&row.counts, row.all_variant);
                    }
                }
            }
            Self::Regions(regions) => regions
                .lock()
                .expect("invariant: region accumulator not poisoned")
                .fold(totals),
            Self::None => {}
        }
    }

    /// Merges another same-kind partial's state in.
    pub(super) fn merge(&mut self, other: &mut Self) {
        match (self, other) {
            (Self::Bins(dst), Self::Bins(src)) => dst.merge_bins(src),
            (Self::Bins(dst), Self::Rows(src)) => {
                for row in src {
                    dst.push_site(row.pos0, &row.counts, row.all_variant);
                }
            }
            (Self::Rows(dst), Self::Rows(src)) => dst.append(src),
            (Self::Regions(dst), Self::Regions(src)) => {
                if !Arc::ptr_eq(dst, src) {
                    dst.lock()
                        .expect("invariant: region accumulator not poisoned")
                        .merge_regions(
                            &src.lock()
                                .expect("invariant: region accumulator not poisoned"),
                        );
                }
            }
            (Self::None, Self::None) => {}
            _ => unreachable!("invariant: one grid kind per run"),
        }
    }
}

/// The run's grid, resolved for accumulation: fixed bp windows or explicit
/// per-contig regions (`--bed-file`).
#[derive(Clone, Debug)]
pub(super) enum GridMode {
    Bp {
        size: u64,
        step: u64,
        width: u64,
    },
    Regions(Arc<ContigRegions>),
    /// No grid at all (the SFS-only composition; its output is genome-wide
    /// spectra, nothing renders a window row).
    None,
}

/// Worker → sink partial: per-contig bin accumulators in record order.
#[derive(Default)]
pub struct StatsPartial {
    pub(super) segments: Vec<SegmentPartial>,
    pub(super) n_records: u64,
    pub(super) n_multiallelic_skipped: u64,
    /// Site-wide invariant sites (no ALT called by any sample).
    pub(super) n_invariant_sites: u64,
    /// Records whose ALT column is `.` (the all-sites marker).
    pub(super) n_no_alt_records: u64,
    /// Genome-wide projected spectra (docs/rounds/v5/05-algorithms.md §3.1).
    pub(super) sfs: Option<SfsPartial>,
    /// Genome-wide Patterson f-statistics.
    pub(super) fstats: Option<FStatsPartial>,
    /// LD partials.
    pub(super) ld: Option<LdPartial>,
}

/// The chunk job: line → `LineRef` → GT tallies → group count row → the
/// site's bin accumulators.
pub struct StatsJob {
    pub(super) grid: GridMode,
    pub(super) n_groups: usize,
    pub(super) gt_groups: Arc<[Option<u32>]>,
    pub(super) n_samples: usize,
    pub(super) seed: Arc<[StatAcc]>,
    /// The record-level site selection (`-r/-R` ∧ `-t/-T` ∧
    /// `--sites-file`; doc_v27 V27-14/V27-15). Shared by every worker, so
    /// it sits behind an `Arc`; an empty selector keeps everything.
    pub(super) select: Arc<reed_filter::spec::region::SiteSelector>,
    /// Declared contig lengths (header + `--contig-lengths`); positions
    /// beyond the declared length are rejected (01-review D8).
    pub(super) contig_len: Arc<HashMap<Vec<u8>, u64>>,
    /// Use raw rows instead of full per-site accumulators for Bp windows.
    pub(super) stream_windows: bool,
    /// One shared accumulator per contig for BED windows. Every worker and
    /// seam uses the same `Arc`; the sink merges each contig's state once,
    /// while each partial keeps its own local position span (V27-19).
    pub(super) region_accs: SharedRegionMap,
    /// The BCF-native fast path (doc_v3 02-bcf.md §2.4/§2.5): present only
    /// when the input is BCF.
    pub(super) bcf: Option<BcfNative>,
    /// The genome-wide projected-SFS plan; `None` = not requested.
    pub(super) sfs: Option<SfsPlan>,
    /// The genome-wide f-statistics plan; `None` = not requested.
    pub(super) fstats: Option<FStatsPlan>,
    /// The LD plane layout (group-major); present only with an LD request.
    pub(super) ld_layout: Option<Arc<reed_stats::SampleLayout>>,
}

/// The BCF fast path's shared state: the dictionaries, the FORMAT index of
/// `GT` (resolved once) and the byte-pair GT table.
pub(super) struct BcfNative {
    pub(super) dict: Arc<BcfDict>,
    pub(super) gt_key: Option<u32>,
    pub(super) table: reed_io::bcf::GtTable,
}

/// Worker-local LD state: the run being filled for the current contig.
#[derive(Default)]
pub struct LdScratch {
    pub(super) contig: Vec<u8>,
    /// The run's site-major planes (4·w words per site).
    pub(super) flat: Vec<u64>,
    /// The run's site positions.
    pub(super) pos0: Vec<u64>,
}

impl LdScratch {
    /// Moves the open run into the chunk partial (a no-op on an empty run).
    pub(super) fn flush(&mut self, partial: &mut StatsPartial) {
        if self.pos0.is_empty() {
            return;
        }
        let run = LdRun {
            contig: std::mem::take(&mut self.contig),
            flat: std::mem::take(&mut self.flat),
            pos0: std::mem::take(&mut self.pos0),
        };
        partial
            .ld
            .get_or_insert_with(LdPartial::default)
            .runs
            .push(run);
    }
}

/// The positional-block size of the parallel decay scan: coarse enough that
/// the per-block lookback span stays a small multiple of the block, fine
/// enough to fill every thread.
pub(super) const LD_BLOCK: usize = 4096;

/// Positional-block parallel decay over one contig's matrix
/// (docs/rounds/v6/02-models.md §2.2). Each block pairs its sites with
/// every predecessor within `max_bp` — read-only, so the blocks are
/// independent; the results merge in block order, which is the serial
/// `(later site, earlier site)` order, keeping the output deterministic
/// across thread counts. The pool comes from the budget, never the rayon
/// global pool (doc_v8 D8-1).
pub(super) fn decay_contig(
    m: &reed_stats::ContigMatrix,
    max_bp: u64,
    edges: &[u64],
    collect_pairs: bool,
    budget: &ThreadBudget,
) -> reed_stats::LdDecaySet {
    use rayon::prelude::*;
    let n = m.n_sites();
    if n < 2 {
        return reed_stats::LdDecaySet::new(edges, m.layout().n_groups());
    }
    let ranges: Vec<std::ops::Range<usize>> = (0..n)
        .step_by(LD_BLOCK)
        .map(|s| s..(s + LD_BLOCK).min(n))
        .collect();
    let parts: Vec<reed_stats::LdDecaySet> =
        reed_engine::budget::run_auxiliary(budget, |parallel| {
            if parallel {
                ranges
                    .par_iter()
                    .map(|r| reed_stats::decay_range(m, r.clone(), max_bp, edges, collect_pairs))
                    .collect()
            } else {
                ranges
                    .iter()
                    .map(|r| reed_stats::decay_range(m, r.clone(), max_bp, edges, collect_pairs))
                    .collect()
            }
        })
        .expect("invariant: the auxiliary pool builds");
    let mut set = reed_stats::LdDecaySet::new(edges, m.layout().n_groups());
    for part in &parts {
        set.merge(part);
    }
    set
}

/// Worker scratch: the group count row, the shared record decoder and its
/// per-record policy counters.
#[derive(Default)]
pub struct StatsScratch {
    pub(super) row: Vec<SiteGroupCounts>,
    pub(super) scanner: Option<GenotypeScanner>,
    pub(super) counters: PolicyCounters,
    /// The native BCF path's decoded GT codes (one per sample, per record).
    pub(super) codes: Vec<GtCode>,
    /// Per-group projection-weight memos for the SFS path.
    pub(super) sfs_caches: Vec<ProjectionCache>,
    /// The per-site projection-weight row the SFS joint accumulation reads
    /// (reused across sites instead of allocating per site; doc_v15 D15-35).
    pub(super) sfs_weights: Vec<Option<WeightsArc>>,
    /// Worker-local LD state.
    pub(super) ld: LdScratch,
}

/// The statistics decode policies: haploid calls are half calls and ALT
/// indices above 1 collapse onto the single ALT allele. `StatsJob` rejects
/// multiallelic records at its site gate before this decoder runs (V27-39).
/// The `stats` QC reports also use these policies and retain their existing
/// collapsed-ALT convention; they do not pass through that popg-only gate.
pub(super) fn stats_policies() -> GenotypePolicies {
    GenotypePolicies {
        ploidy: PloidyPolicy::HaploidAsHalfCall,
        multiallelic: MultiallelicPolicy::CollapseAlt,
        ..GenotypePolicies::default()
    }
}

impl StatsJob {
    /// Builds the SFS-only job the convert driver composes (doc_v13 §7.5
    /// plan A): the statistics chunk path restricted to the genome-wide SFS
    /// fold. The same line/site/scan stages run — so `reed convert --sfs`
    /// and `reed popg --sfs` cannot drift — but no window accumulator is
    /// ever allocated (`GridMode::None`).
    pub(crate) fn for_sfs(
        n_samples: usize,
        gt_groups: Vec<Option<u32>>,
        plan: SfsPlan,
        contig_len: Arc<HashMap<Vec<u8>, u64>>,
    ) -> Self {
        Self {
            grid: GridMode::None,
            stream_windows: false,
            n_groups: plan.projections.len(),
            gt_groups: gt_groups.into(),
            n_samples,
            seed: Vec::<StatAcc>::new().into(),
            select: Arc::new(reed_filter::spec::region::SiteSelector::default()),
            contig_len,
            bcf: None,
            region_accs: Arc::new(std::sync::Mutex::new(HashMap::new())),
            sfs: Some(plan),
            fstats: None,
            ld_layout: None,
        }
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
            self.push_line(line, &mut current, &mut partial, scratch)?;
        }
        if let Some(segment) = current.take() {
            partial.segments.push(segment);
        }
        scratch.ld.flush(&mut partial);
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
        scratch.ld.flush(&mut partial);
        Ok(partial)
    }
}

impl StatsPartial {
    /// Absorbs another partial's SFS state in chunk order (the convert
    /// sink's merge; doc_v13 §7.5). Every other field is dropped: the
    /// convert driver only consumes the spectra.
    pub(crate) fn merge_sfs(&mut self, other: &Self) {
        if let Some(src) = &other.sfs {
            match &mut self.sfs {
                Some(dst) => dst.merge(src),
                None => self.sfs = Some(src.clone()),
            }
        }
    }

    /// The merged SFS partial (convert's drain).
    #[must_use]
    pub(crate) fn sfs_partial(&self) -> Option<&SfsPartial> {
        self.sfs.as_ref()
    }

    /// Records the partial absorbed (convert's report counter).
    #[must_use]
    pub(crate) fn records(&self) -> u64 {
        self.n_records
    }
}

impl StatsJob {
    fn reset_row(&self, scratch: &mut StatsScratch) {
        scratch.row.clear();
        scratch
            .row
            .resize(self.n_groups, SiteGroupCounts::default());
    }

    /// The site stage both record sources share (text and native BCF): the
    /// streamed-target predicate, the biallelic gate, contig segment
    /// bookkeeping, the POS contract and invariant accounting. Multiallelic
    /// records are rejected here, before any segment or accumulator is
    /// touched (doc_v27 V27-39). Returns the 0-based position when the
    /// record continues to the genotype stage.
    fn site_prefix(
        &self,
        contig: &[u8],
        pos1: u64,
        n_alt: u32,
        is_invariant: bool,
        current: &mut Option<SegmentPartial>,
        partial: &mut StatsPartial,
    ) -> Result<Option<u64>, EngineError> {
        // The site selection gates every record before it reaches an
        // accumulator: `-r/-R`, `-t/-T` and `--sites-file` are three
        // families that intersect (doc_v27 V27-14/V27-15).
        if !self.select.is_empty() && !self.select.allows(contig, pos1) {
            return Ok(None);
        }
        // Position validation is independent of the multiallelic filter.
        // A filtered record must not turn an otherwise descending input
        // into a successful run (doc_v27 V27-28).
        if let Some(segment) = current.as_ref()
            && segment.contig == contig
            && pos1
                .checked_sub(1)
                .is_some_and(|pos0| pos0 < segment.last_pos0)
        {
            return Err(EngineError::from(
                reed_core::CoreError::PositionNotAscending {
                    prev_pos0: segment.last_pos0,
                    pos0: pos1.saturating_sub(1),
                },
            ));
        }
        // The default `-M 2` gate lives before segment creation. A filtered
        // multiallelic record therefore cannot advance the trailing window,
        // cannot enter a denominator, and cannot bypass the ascending-POS
        // check of a later biallelic record (doc_v27 V27-39/V27-28).
        if n_alt > 1 {
            partial.n_multiallelic_skipped += 1;
            return Ok(None);
        }
        match current {
            Some(segment) if segment.contig == contig => {}
            _ => {
                if let Some(segment) = current.take() {
                    partial.segments.push(segment);
                }
                let acc = match &self.grid {
                    GridMode::Bp { width, .. } => {
                        if self.stream_windows {
                            ContigAcc::Rows(Vec::new())
                        } else {
                            ContigAcc::Bins(ContigBins::new(BpBinning::new(*width), &self.seed))
                        }
                    }
                    GridMode::Regions(regions) => {
                        let set = regions
                            .get(contig)
                            .unwrap_or_else(|| Arc::new(RegionSet::default()));
                        let mut shared = self
                            .region_accs
                            .lock()
                            .expect("invariant: region accumulator map not poisoned");
                        let acc = shared
                            .entry(contig.to_vec())
                            .or_insert_with(|| {
                                Arc::new(std::sync::Mutex::new(RegionAccs::new(set, &self.seed)))
                            })
                            .clone();
                        ContigAcc::Regions(acc)
                    }
                    GridMode::None => ContigAcc::None,
                };
                *current = Some(SegmentPartial {
                    contig: contig.to_vec(),
                    acc,
                    first_pos0: u64::MAX,
                    last_pos0: 0,
                });
            }
        }
        let Some(pos0) = pos1.checked_sub(1) else {
            return Err(EngineError::job(format!(
                "{}:{}: POS is zero",
                String::from_utf8_lossy(contig),
                pos1
            )));
        };
        // Input contract (04-algorithms.md §6): positions ascending within a
        // contig, never beyond the declared length.
        if let Some(&len) = self.contig_len.get(contig)
            && pos1 > len
        {
            return Err(EngineError::from(
                reed_core::CoreError::PositionBeyondContigLength {
                    contig_name: String::from_utf8_lossy(contig).into_owned(),
                    pos1,
                    len_bp: len,
                },
            ));
        }
        if is_invariant {
            partial.n_no_alt_records += 1;
        }
        Ok(Some(pos0))
    }

    /// The genotype stage both record sources share: the per-group counts,
    /// the ascending-position check and the bin accumulation.
    // The arguments are distinct caller-owned scratch buffers (row, caches,
    // segment, partial, LD scratch); bundling them into a struct would only
    // move the same borrows around.
    #[allow(clippy::too_many_arguments)]
    fn push_row(
        &self,
        pos0: u64,
        codes: &[GtCode],
        row: &mut [SiteGroupCounts],
        sfs_caches: &mut Vec<ProjectionCache>,
        sfs_weights: &mut Vec<Option<WeightsArc>>,
        current: &mut Option<SegmentPartial>,
        partial: &mut StatsPartial,
        ld_scratch: &mut LdScratch,
    ) -> Result<(), EngineError> {
        for counts in row.iter_mut() {
            *counts = SiteGroupCounts::default();
        }
        let mut any_alt = false;
        for (sample_idx, &code) in codes.iter().enumerate() {
            any_alt |= code_has_alt(code);
            if let Some(group) = self.gt_groups.get(sample_idx).copied().flatten() {
                row[group as usize].add_code(code);
            }
        }
        let segment = current
            .as_mut()
            .expect("invariant: segment opened for this contig above");
        if pos0 < segment.last_pos0 {
            return Err(EngineError::from(
                reed_core::CoreError::PositionNotAscending {
                    prev_pos0: segment.last_pos0,
                    pos0,
                },
            ));
        }
        segment.acc.push_site(pos0, row, any_alt);
        segment.first_pos0 = segment.first_pos0.min(pos0);
        segment.last_pos0 = pos0;
        if !any_alt {
            partial.n_invariant_sites += 1;
        }
        self.push_sfs(row, partial, sfs_caches, sfs_weights);
        self.push_fstats(&segment.contig, pos0, row, partial);
        self.push_ld(&segment.contig, pos0, codes, partial, ld_scratch);
        partial.n_records += 1;
        Ok(())
    }

    /// Feeds the LD path for one site (docs/rounds/v6/02-models.md §2.2):
    /// the site's planes join the worker's contig run in the group-major
    /// layout; the pairs are computed position-block parallel after the
    /// stream ends, so no carry or head/tail windows exist here.
    fn push_ld(
        &self,
        contig: &[u8],
        pos0: u64,
        codes: &[GtCode],
        partial: &mut StatsPartial,
        ld: &mut LdScratch,
    ) {
        let Some(layout) = &self.ld_layout else {
            return;
        };
        if ld.contig != contig {
            ld.flush(partial);
            ld.contig.clear();
            ld.contig.extend_from_slice(contig);
        }
        let w = layout.w();
        let base = ld.flat.len();
        ld.flat.resize(base + 4 * w, 0);
        for (sample, &code) in codes.iter().enumerate() {
            let Some(bit) = layout.bit_of(sample) else {
                continue;
            };
            let (word, mask) = (bit / 64, 1_u64 << (bit % 64));
            if code.slot_alt(0) {
                ld.flat[base + word] |= mask;
            }
            if code.slot_missing(0) {
                ld.flat[base + w + word] |= mask;
            }
            if code.slot_alt(1) {
                ld.flat[base + 2 * w + word] |= mask;
            }
            if code.slot_missing(1) {
                ld.flat[base + 3 * w + word] |= mask;
            }
        }
        ld.pos0.push(pos0);
    }

    /// Feeds the genome-wide f-statistics for one site (docs/rounds/v5/05
    /// §3.3). The segment's contig names the jackknife block key.
    fn push_fstats(
        &self,
        contig: &[u8],
        pos0: u64,
        row: &[SiteGroupCounts],
        partial: &mut StatsPartial,
    ) {
        let Some(plan) = &self.fstats else { return };
        let stats = partial
            .fstats
            .get_or_insert_with(|| FStatsPartial::new(plan));
        stats.push_grid(contig, pos0);
        for acc in &mut stats.per_test {
            acc.push(contig, pos0, row);
        }
    }

    /// Feeds the genome-wide projected spectra for one site (docs/rounds/v5/
    /// 05-algorithms.md §3.1): `m < n` sites are skipped and counted. The
    /// joint (2D) accumulators share the same weight vectors, combined as an
    /// outer product per pair (easySFS `_from_count_dict`).
    fn push_sfs(
        &self,
        row: &[SiteGroupCounts],
        partial: &mut StatsPartial,
        sfs_caches: &mut Vec<ProjectionCache>,
        sfs_weights: &mut Vec<Option<WeightsArc>>,
    ) {
        let Some(plan) = &self.sfs else { return };
        let sfs = partial.sfs.get_or_insert_with(|| {
            SfsPartial::new(plan).expect("invariant: projections validated at request time")
        });
        if sfs_caches.len() != sfs.per_group.len() {
            *sfs_caches = plan
                .projections
                .iter()
                .map(|&n| ProjectionCache::new(n))
                .collect();
        }
        sfs_weights.clear();
        for (g, acc) in sfs.per_group.iter_mut().enumerate() {
            let counts = &row[g];
            let m = counts.n_ref_alleles + counts.n_alt_alleles;
            if m >= acc.n() {
                let w = sfs_caches[g].weights(m, counts.n_alt_alleles);
                acc.push_weighted(Some(&w));
                sfs_weights.push(Some(w));
            } else {
                acc.push_weighted(None);
                sfs_weights.push(None);
            }
        }
        if let Some(joint) = &mut sfs.joint {
            for ((g1, g2), acc) in plan.joint_pairs.iter().zip(joint.iter_mut()) {
                acc.push_weighted(
                    sfs_weights[*g1].as_deref().map(|w| w.as_slice()),
                    sfs_weights[*g2].as_deref().map(|w| w.as_slice()),
                );
            }
        }
    }

    /// The one line-entry point for the text path (doc_v25 V25-03): blank
    /// lines are skipped (a seam can deliver `"\n"` or `"\r\n"` as the
    /// stitched line), malformed lines stay fatal — stats' strict policy is
    /// unchanged. `classify_line` owns the judgement.
    fn push_line(
        &self,
        line: &[u8],
        current: &mut Option<SegmentPartial>,
        partial: &mut StatsPartial,
        scratch: &mut StatsScratch,
    ) -> Result<(), EngineError> {
        let lr = match reed_io::vcf::line::classify_line(line, self.n_samples) {
            reed_io::vcf::line::RecordClass::Record(lr) => lr,
            reed_io::vcf::line::RecordClass::Blank => return Ok(()),
            reed_io::vcf::line::RecordClass::Malformed(why) => {
                return Err(EngineError::job(format!(
                    "record line is malformed ({})",
                    why.label()
                )));
            }
        };
        let site = SiteView::parse(&lr).map_err(EngineError::from)?;
        let Some(pos0) = self.site_prefix(
            site.contig_name,
            site.pos1,
            site.n_alt(),
            site.is_invariant(),
            current,
            partial,
        )?
        else {
            return Ok(());
        };
        let StatsScratch {
            row,
            scanner,
            counters,
            sfs_caches,
            sfs_weights,
            ld,
            ..
        } = scratch;
        *counters = PolicyCounters::default();
        let n_samples = self.n_samples;
        let scanner = scanner
            .get_or_insert_with(|| GenotypeScanner::new(n_samples, &[], &[], stats_policies()));
        let out = scanner
            .scan(&lr, &site, counters)
            .map_err(EngineError::from)?;
        self.push_row(
            pos0,
            out.codes,
            row,
            sfs_caches,
            sfs_weights,
            current,
            partial,
            ld,
        )
    }

    /// The record sink for a filtering wrapper (doc_v4 O-2): the wrapper
    /// already decoded and masked the genotypes, so this path parses only
    /// the fixed columns again and skips both the GT decode and the
    /// rewritten-text round-trip.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] for malformed records.
    pub fn run_records(
        &self,
        records: &[crate::filtered::ScannedRecord<'_>],
        scratch: &mut StatsScratch,
    ) -> Result<StatsPartial, EngineError> {
        let mut partial = StatsPartial::default();
        let mut current: Option<SegmentPartial> = None;
        self.reset_row(scratch);
        for rec in records {
            let lr = match reed_io::vcf::line::classify_line(rec.line, self.n_samples) {
                reed_io::vcf::line::RecordClass::Record(lr) => lr,
                reed_io::vcf::line::RecordClass::Blank => {
                    return Err(EngineError::job("blank record line"));
                }
                reed_io::vcf::line::RecordClass::Malformed(why) => {
                    return Err(EngineError::job(format!(
                        "record line is malformed ({})",
                        why.label()
                    )));
                }
            };
            let site = SiteView::parse(&lr).map_err(EngineError::from)?;
            let Some(pos0) = self.site_prefix(
                site.contig_name,
                site.pos1,
                site.n_alt(),
                site.is_invariant(),
                &mut current,
                &mut partial,
            )?
            else {
                continue;
            };
            self.push_row(
                pos0,
                rec.codes,
                &mut scratch.row,
                &mut scratch.sfs_caches,
                &mut scratch.sfs_weights,
                &mut current,
                &mut partial,
                &mut scratch.ld,
            )?;
        }
        if let Some(segment) = current.take() {
            partial.segments.push(segment);
        }
        scratch.ld.flush(&mut partial);
        Ok(partial)
    }

    /// The native BCF record sink (doc_v3 02-bcf.md §2.4/§2.5): each framed
    /// record is decoded straight into the shared site/row stages — the
    /// shared block past its header and every unplanned FORMAT field are
    /// skipped by descriptor length, so no text is ever produced.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] for malformed records and policy rejections.
    pub fn run_bcf_records(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        records: &[&[u8]],
        scratch: &mut StatsScratch,
    ) -> Result<StatsPartial, EngineError> {
        let Some(bcf) = &self.bcf else {
            return Err(EngineError::job(
                "internal: BCF records without a BCF input",
            ));
        };
        let mut partial = StatsPartial::default();
        let mut current: Option<SegmentPartial> = None;
        self.reset_row(scratch);
        for bytes in records {
            let rec = reed_io::bcf::record_ref(bytes).map_err(EngineError::from)?;
            self.push_bcf_record(&rec, bcf, &mut current, &mut partial, scratch)?;
        }
        if let Some(segment) = current.take() {
            partial.segments.push(segment);
        }
        scratch.ld.flush(&mut partial);
        Ok(partial)
    }

    fn push_bcf_record(
        &self,
        rec: &reed_io::bcf::BcfRecordRef<'_>,
        bcf: &BcfNative,
        current: &mut Option<SegmentPartial>,
        partial: &mut StatsPartial,
        scratch: &mut StatsScratch,
    ) -> Result<(), EngineError> {
        let header = reed_io::bcf::site_header(rec).map_err(EngineError::from)?;
        let contig = bcf.dict.contig(header.chrom).ok_or_else(|| {
            EngineError::job(format!(
                "BCF contig index {} is not in the header dictionary",
                header.chrom
            ))
        })?;
        let Some(pos1) = u64::try_from(header.pos0)
            .ok()
            .and_then(|p| p.checked_add(1))
        else {
            return Err(EngineError::job("BCF POS is negative"));
        };
        let n_alt = u32::try_from(header.n_allele.saturating_sub(1)).unwrap_or(u32::MAX);
        let Some(pos0) = self.site_prefix(contig, pos1, n_alt, n_alt == 0, current, partial)?
        else {
            return Ok(());
        };
        scratch.counters = PolicyCounters::default();
        let query = reed_io::bcf::GtQuery {
            gt_key: bcf.gt_key,
            n_samples: self.n_samples,
            contig_name: contig,
            pos1,
        };
        reed_io::bcf::decode_record_gt(
            rec,
            &header,
            &query,
            &bcf.table,
            &mut scratch.codes,
            &mut scratch.counters,
        )
        .map_err(EngineError::from)?;
        self.push_row(
            pos0,
            &scratch.codes,
            &mut scratch.row,
            &mut scratch.sfs_caches,
            &mut scratch.sfs_weights,
            current,
            partial,
            &mut scratch.ld,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use reed_io::chunk::{ChunkBytes, ChunkPayload};

    fn sfs_job() -> StatsJob {
        let plan = crate::stats::request::SfsPlan {
            projections: vec![2].into(),
            folded: false,
            joint_pairs: Vec::new().into(),
        };
        StatsJob::for_sfs(2, vec![Some(0), Some(0)], plan, Arc::new(HashMap::new()))
    }

    fn plain_chunk() -> ChunkBytes {
        ChunkBytes {
            id: reed_core::ChunkId(0),
            compressed: ChunkPayload::from(Vec::new()),
            first_in_block_offset: 0,
            plain: true,
            compressed_start: 0,
            blocks: Vec::new(),
            truncate_to: None,
        }
    }

    /// doc_v27 V27-04: the real 4096-site block partition must retain the
    /// cross-block pairs of a single-site tail. These are pair counts, not
    /// merely a smoke test: 4097 sites used to lose 4096 pairs, and the
    /// whole 4096k+1 family was affected.
    #[test]
    fn ld_block_boundaries_keep_all_pairs_including_a_single_site_tail() {
        let budget = ThreadBudget {
            threads: 2,
            mem: reed_engine::budget::MemBudget::from_bytes(u64::MAX),
        };
        for n in [4095_usize, 4096, 4097, 4098, 8193] {
            let layout = Arc::new(reed_stats::SampleLayout::from_assignments(
                &[Some(0); 4],
                1,
                None,
            ));
            let mut matrix = reed_stats::ContigMatrix::new(b"chr1".to_vec(), layout);
            for pos0 in 0..n as u64 {
                // Four samples at dosages 0/1/2/1: every pair has r² = 1.
                matrix.push_run(vec![0b1010, 0, 0b1100, 0], vec![pos0]);
            }
            let set = decay_contig(&matrix, u64::MAX, &[0, u64::MAX], true, &budget);
            let want = (n * (n - 1) / 2) as u64;
            assert_eq!(set.pairs.len() as u64, want, "n={n}");
            assert_eq!(set.cohort.counts().iter().sum::<u64>(), want, "n={n}");
        }
    }

    /// The stats line policy (doc_v25 V25-03): a blank seam record (`"\n"`)
    /// or a CRLF blank line in the body is skipped; a malformed line stays
    /// fatal, with the old diagnostic still in the message.
    #[test]
    fn blank_lines_are_skipped_and_malformed_stays_fatal() {
        let job = sfs_job();
        let mut scratch = StatsScratch::default();
        for blank in [b"\n".as_slice(), b"\r\n"] {
            let partial = job
                .run_seam(reed_core::ChunkId(1), blank, &mut scratch)
                .expect("a blank seam line is not a record");
            assert_eq!(partial.n_records, 0);
        }
        // A CRLF blank (`b"\r"` after the split) between two records.
        let body = b"chr1\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n\r\n\
                     chr1\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\n";
        let partial = job
            .run_chunk(&plain_chunk(), body, &mut scratch)
            .expect("a blank body line is not a record");
        assert_eq!(partial.n_records, 2);
        let Err(err) = job.run_seam(reed_core::ChunkId(1), b"chr1\t1\t.\n", &mut scratch) else {
            panic!("a malformed line must stay fatal");
        };
        assert!(err.to_string().contains("fewer than 8 columns"), "{err}");
    }
}
