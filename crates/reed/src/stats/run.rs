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
use reed_io::vcf::line::LineRef;
use reed_io::vcf::site::SiteView;
use reed_stats::request::StatAcc;
use reed_stats::sfs::{ProjectionCache, WeightsArc};
use reed_stats::window::{BpBinning, ContigBins, RegionAccs};
use std::collections::HashMap;
use std::sync::Arc;

/// One contig's run of sites: bin accumulators plus the span it covers.
pub(super) struct SegmentPartial {
    pub(super) contig: Vec<u8>,
    pub(super) acc: ContigAcc,
    pub(super) first_pos0: u64,
    pub(super) last_pos0: u64,
}

/// The per-contig accumulation state, one variant per grid kind: `Bp` bins
/// (watermark window emission) or per-region accumulators (`--bed-file`).
#[derive(Clone, Debug)]
pub(super) enum ContigAcc {
    Bins(ContigBins),
    Regions(RegionAccs),
    /// No window bookkeeping at all: the SFS-only composition the convert
    /// driver builds (`StatsJob::for_sfs`, doc_v13 §7.5 plan A).
    None,
}

impl ContigAcc {
    pub(super) fn push_site(&mut self, pos0: u64, row: &[SiteGroupCounts], all_variant: bool) {
        match self {
            Self::Bins(bins) => bins.push_site(pos0, row, all_variant),
            Self::Regions(regions) => regions.push_site(pos0, row, all_variant),
            Self::None => {}
        }
    }

    /// Registers a present-but-uncountable site (a skipped multiallelic
    /// record) so the variant-only pi/dxy denominator counts it as present
    /// (doc_v7 D7-11).
    pub(super) fn push_skipped(&mut self, pos0: u64) {
        match self {
            Self::Bins(bins) => bins.push_skipped(pos0),
            Self::Regions(regions) => regions.push_skipped(pos0),
            Self::None => {}
        }
    }

    pub(super) fn fold(&self, totals: &mut [StatAcc]) {
        match self {
            Self::Bins(bins) => bins.fold(totals),
            Self::Regions(regions) => regions.fold(totals),
            Self::None => {}
        }
    }

    /// Merges another same-kind partial's state in.
    pub(super) fn merge(&mut self, other: &mut Self) {
        match (self, other) {
            (Self::Bins(dst), Self::Bins(src)) => dst.merge_bins(src),
            (Self::Regions(dst), Self::Regions(src)) => dst.merge_regions(src),
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
    /// Streamed target regions (`-t/-T`): kept when inside any include and
    /// no exclude; empty lists keep everything.
    pub(super) targets: Arc<[reed_filter::spec::region::RegionTriple]>,
    pub(super) targets_exclude: Arc<[reed_filter::spec::region::RegionTriple]>,
    /// Declared contig lengths (header + `--contig-lengths`); positions
    /// beyond the declared length are rejected (01-review D8).
    pub(super) contig_len: Arc<HashMap<Vec<u8>, u64>>,
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
const LD_BLOCK: usize = 4096;

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
    /// Worker-local LD state.
    pub(super) ld: LdScratch,
}

/// The statistics decode policies: haploid calls are half calls and ALT
/// indices above 1 collapse onto the single ALT allele (the biallelic
/// view; multi-allelic records are skipped before the scan anyway).
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
            n_groups: plan.projections.len(),
            gt_groups: gt_groups.into(),
            n_samples,
            seed: Vec::<StatAcc>::new().into(),
            targets: Vec::new().into(),
            targets_exclude: Vec::new().into(),
            contig_len,
            bcf: None,
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
            if line.is_empty() {
                continue;
            }
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
    /// streamed-target predicate, the contig segment bookkeeping, the POS
    /// contract and the invariant/multiallelic accounting. Returns the
    /// 0-based position when the record continues to the genotype stage.
    fn site_prefix(
        &self,
        contig: &[u8],
        pos1: u64,
        n_alt: u32,
        is_invariant: bool,
        current: &mut Option<SegmentPartial>,
        partial: &mut StatsPartial,
    ) -> Result<Option<u64>, EngineError> {
        // Streamed targets (`-t/-T`): the region predicate gates every record
        // before it reaches an accumulator.
        if (!self.targets.is_empty() || !self.targets_exclude.is_empty())
            && !reed_filter::spec::region::region_allows(
                &self.targets,
                &self.targets_exclude,
                contig,
                pos1,
            )
        {
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
                        ContigAcc::Bins(ContigBins::new(BpBinning::new(*width), &self.seed))
                    }
                    GridMode::Regions(regions) => {
                        let set = regions
                            .get(contig)
                            .unwrap_or_else(|| Arc::new(RegionSet::default()));
                        ContigAcc::Regions(RegionAccs::new(set, &self.seed))
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
        if n_alt > 1 {
            partial.n_multiallelic_skipped += 1;
            // The file contains this site even though pi/dxy skip it, so the
            // variant-only denominator must not treat it as absent
            // (doc_v7 D7-11).
            if let Some(segment) = current.as_mut() {
                segment.acc.push_skipped(pos0);
            }
            return Ok(None);
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
        self.push_sfs(row, partial, sfs_caches);
        self.push_fstats(&segment.contig.clone(), pos0, row, partial);
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
        let mut weights: Vec<Option<WeightsArc>> = Vec::with_capacity(row.len());
        for (g, acc) in sfs.per_group.iter_mut().enumerate() {
            let counts = &row[g];
            let m = counts.n_ref_alleles + counts.n_alt_alleles;
            if m >= acc.n() {
                let w = sfs_caches[g].weights(m, counts.n_alt_alleles);
                acc.push_weighted(Some(&w));
                weights.push(Some(w));
            } else {
                acc.push_weighted(None);
                weights.push(None);
            }
        }
        if let Some(joint) = &mut sfs.joint {
            for ((g1, g2), acc) in plan.joint_pairs.iter().zip(joint.iter_mut()) {
                acc.push_weighted(
                    weights[*g1].as_deref().map(|w| w.as_slice()),
                    weights[*g2].as_deref().map(|w| w.as_slice()),
                );
            }
        }
    }

    fn push_line(
        &self,
        line: &[u8],
        current: &mut Option<SegmentPartial>,
        partial: &mut StatsPartial,
        scratch: &mut StatsScratch,
    ) -> Result<(), EngineError> {
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::job("record line has fewer than 9 columns"));
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
        self.push_row(pos0, out.codes, row, sfs_caches, current, partial, ld)
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
            let Some(lr) = LineRef::parse(rec.line) else {
                return Err(EngineError::job("record line has fewer than 9 columns"));
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
            current,
            partial,
            &mut scratch.ld,
        )
    }
}
