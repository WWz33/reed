//! The unified pass0 (`docs/rounds/v6/03-architecture.md` §7): ONE chunk
//! job fills the site-major bit-plane matrices every global analysis reads.
//! Before Y1e three jobs normalized the same stages three ways
//! (`MissingJob`, `SketchJob`, `LdPruneJob`); this job replaces them.
//!
//! For every record the shared [`RecordPipeline`] runs head → decode →
//! mask; every head survivor contributes one site of `(4 + R) · w` words —
//! the raw `[slot0_alt, slot0_miss, slot1_alt, slot1_miss]` planes plus one
//! firing bitmap per mask rule. From those planes the analyses derive
//! everything else (doc_v6 02 §2.2):
//!
//! * `--mind` (doc_v4 02 回退版 §3.2): a sample's called alleles over the
//!   head survivors are `2·sites − colpop(slot0_miss ∨ mask) −
//!   colpop(slot1_miss ∨ mask)` — the per-column popcount of the post-mask
//!   miss planes, accumulated here so no second counting loop can drift.
//! * `--indep-pairwise`: the per-contig matrices feed the pruner, which
//!   sees only the sites the tail rules kept (the site rules run before
//!   the pruning, per 02-cli §6 step 6b).
//! * The sketch path's mid stage: the planes and rule bitmaps replace the
//!   chunk sketch — the decisions are recomputed over the final keep set
//!   from the same bytes.

use std::sync::{Arc, Mutex};

use reed_core::genotype::GtCode;
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{ChunkJob, OrderedSink, run};
use reed_filter::funnel::FunnelPartial;
use reed_filter::plan::FilterPlan;

use reed_core::ChunkId;

use crate::filter::record_pipeline::{NativeScratch, PipelineOptions, RecordPipeline};
use crate::filter::{FilterOptions, ferr};
use reed_filter::FilterError;

/// What this pass0 run collects (decided by the spec, `pass0::plan_mode`).
#[derive(Clone, Copy, Debug)]
pub(crate) struct MatrixMode {
    /// Retain the per-contig site matrices (the pruner or the mid stage
    /// reads them after the pass).
    pub retain: bool,
    /// Record one firing bitmap per mask rule (the mid stage's per-rule
    /// attribution needs them; the planes alone carry the effective mask).
    pub rule_bits: bool,
    /// Evaluate the tail rules inline (the keep set is already final: no
    /// `--mind`); the matrix then holds exactly the sites the pruner sees.
    pub tail_inline: bool,
    /// Record the per-site effective mask bitmap (the pruner excludes
    /// masked samples, doc_v7 D7-4). Set whenever the pruner runs.
    pub mask_plane: bool,
}

/// One pass0 partial: a chunk body's or a seam record's stream unit. The
/// seam record belongs to the previous chunk's stream (§6.2), so the sink
/// merges it into the chunk it closes.
#[derive(Default)]
pub(crate) struct MatrixPartial {
    /// This partial is a seam record (belongs to the previous chunk).
    pub seam: bool,
    /// Bit per processed record: the head rules passed.
    pub head_bits: Vec<u64>,
    /// Records processed (the head bitmap's population).
    pub n_records: usize,
    /// Head-surviving records (the `--mind` denominator, `2L`).
    pub n_head: u64,
    /// Per-sample post-mask called alleles over the head survivors.
    pub called: Vec<u64>,
    /// The head survivors' site-major planes, `(4 + R) · w` words each, in
    /// record order (empty when the mode retains nothing).
    pub flat: Vec<u64>,
    /// The sites' 0-based positions, aligned with `flat`.
    pub pos0: Vec<u64>,
    /// Whether each site's record carried a GT field, aligned with `flat`
    /// (the mid stage's `-g/--genotype` short-circuit needs it: a GT-less
    /// record and an all-`./.` record have identical planes).
    pub has_gt: Vec<bool>,
    /// Contig runs over the sites: `(first site offset, contig name)`.
    pub contigs: Vec<(usize, Vec<u8>)>,
    /// This unit's funnel (head hits, decode policies).
    pub funnel: FunnelPartial,
}

/// Worker scratch: the shared record pipeline plus the per-site buffers.
pub(crate) struct MatrixScratch {
    pipeline: RecordPipeline,
    funnel: FunnelPartial,
    head_bits: Vec<u64>,
    n_records: usize,
    /// Per-sample post-mask slot-missing tallies over the unit's sites
    /// (the `--mind` column popcounts; the per-unit fold is
    /// `called = 2·sites − m0 − m1`).
    m0: Vec<u64>,
    m1: Vec<u64>,
    n_head: u64,
    /// The current site's plane buffer, `(4 + R) · w` words, reused.
    site: Vec<u64>,
    flat: Vec<u64>,
    pos0: Vec<u64>,
    /// Whether each retained site's record carried a GT field, aligned
    /// with `pos0` (the mid stage's `-g` short-circuit).
    has_gt: Vec<bool>,
    contigs: Vec<(usize, Vec<u8>)>,
    /// Native BCF input buffers (doc_v4 O-4).
    native: Option<NativeScratch>,
}

/// The unified pass0 job.
pub(crate) struct MatrixJob {
    pub plan: Arc<FilterPlan>,
    pub n_samples: usize,
    /// Words per sample bitmap (`ceil(n_samples / 64)`, at least 1).
    pub w: usize,
    pub mode: MatrixMode,
    /// Accumulate the per-sample called alleles (`--mind` requested).
    pub want_counts: bool,
    /// The kept-sample count the inline tail rules take as their
    /// denominator (`mode.tail_inline` only; `None` defers the tail).
    pub tail_n_samples: Option<u64>,
    /// The static sample subset (the inline tail counts and the mask stage
    /// cover kept samples exactly like pass1's single pass does; the
    /// dynamic `--mind` subset does not exist yet in pass0).
    pub keep: Arc<[bool]>,
    /// BCF input: when set, the job decodes the records natively instead of
    /// rendering them (doc_v4 O-4).
    pub bcf_input: Option<Arc<reed_io::bcf::BcfDict>>,
    /// The GT decode table for the native path, built once per run and
    /// shared by every worker (A6-10). `Some` exactly when `bcf_input` is.
    pub gt_table: Option<Arc<reed_io::bcf::GtTable>>,
    /// The plan's Integer FORMAT keys.
    pub int_keys: Arc<[Vec<u8>]>,
    /// The plan's Float FORMAT keys.
    pub float_keys: Arc<[Vec<u8>]>,
}

impl MatrixJob {
    fn build_scratch(&self) -> MatrixScratch {
        let int_keys: Vec<&[u8]> = self.plan.int_fields.iter().map(Vec::as_slice).collect();
        let float_keys: Vec<&[u8]> = self.plan.float_fields.iter().map(Vec::as_slice).collect();
        let mut pipeline = RecordPipeline::new(
            self.n_samples,
            &int_keys,
            &float_keys,
            self.plan.policies,
            &[],
        );
        pipeline.resize_for(&self.plan);
        MatrixScratch {
            pipeline,
            funnel: FunnelPartial {
                head_hits: vec![0; self.plan.head.len()],
                tail_hits: vec![0; self.plan.tail.len()],
                mask_hits: vec![0; self.plan.mask.len()],
                ..FunnelPartial::default()
            },
            head_bits: Vec::new(),
            n_records: 0,
            m0: vec![0; self.n_samples],
            m1: vec![0; self.n_samples],
            n_head: 0,
            site: vec![0; self.stride()],
            flat: Vec::new(),
            pos0: Vec::new(),
            has_gt: Vec::new(),
            contigs: Vec::new(),
            native: match (&self.bcf_input, &self.gt_table) {
                (Some(_), Some(table)) => Some(NativeScratch::new(
                    std::sync::Arc::clone(table),
                    self.n_samples,
                    &self.int_keys,
                    &self.float_keys,
                )),
                _ => None,
            },
        }
    }

    /// Words one site occupies in the flat buffer.
    fn stride(&self) -> usize {
        (4 + usize::from(self.mode.rule_bits) * self.plan.mask.len()
            + usize::from(self.mode.mask_plane))
            * self.w
    }

    /// The pipeline options this mode asks of the shared stages: the mask
    /// runs for every sample unless the keep set is final and only the
    /// local counts need it (the prune pass0); the tail runs inline only
    /// when no dynamic rule can shift it.
    fn opts(&self) -> PipelineOptions {
        PipelineOptions {
            mask_all_samples: self.mode.rule_bits || self.want_counts,
            require_contig: false,
            tail_n_samples: if self.mode.tail_inline {
                self.tail_n_samples
            } else {
                None
            },
            counts: self.mode.tail_inline,
        }
    }

    /// Resets the unit's accumulators (a chunk body or a seam record).
    fn reset_unit(&self, scratch: &mut MatrixScratch) {
        scratch.funnel.reset(
            self.plan.head.len(),
            self.plan.tail.len(),
            self.plan.mask.len(),
        );
        scratch.head_bits.clear();
        scratch.n_records = 0;
        for slot in scratch.m0.iter_mut().chain(scratch.m1.iter_mut()) {
            *slot = 0;
        }
        scratch.n_head = 0;
        scratch.flat.clear();
        scratch.pos0.clear();
        scratch.has_gt.clear();
        scratch.contigs.clear();
    }

    /// Records one processed record's head outcome after its body ran.
    fn record_head_bit(&self, scratch: &mut MatrixScratch, head_passed: bool) {
        let r = scratch.n_records;
        let word = r / 64;
        if word >= scratch.head_bits.len() {
            scratch.head_bits.resize(word + 1, 0);
        }
        if head_passed {
            scratch.head_bits[word] |= 1 << (r % 64);
            scratch.n_head += 1;
        }
        scratch.n_records += 1;
    }

    /// The survivor body both record paths share: the site's planes go
    /// into the reusable buffer, the `--mind` popcounts fold from them,
    /// and the retained modes append the site to the unit's run.
    #[allow(clippy::too_many_arguments)]
    fn survivor_body(
        &self,
        site: &mut [u64],
        m0: &mut [u64],
        m1: &mut [u64],
        flat: &mut Vec<u64>,
        pos0: &mut Vec<u64>,
        contigs: &mut Vec<(usize, Vec<u8>)>,
        has_gt: &mut Vec<bool>,
        codes: &[GtCode],
        mask_rule: &[u16],
        masked: &reed_core::Bitset,
        contig: &[u8],
        pos1: u64,
        record_has_gt: bool,
    ) {
        let w = self.w;
        site.iter_mut().for_each(|word| *word = 0);
        for (i, &code) in codes.iter().enumerate() {
            let (word, bit) = (i / 64, 1_u64 << (i % 64));
            if code.slot_alt(0) {
                site[word] |= bit;
            }
            if code.slot_missing(0) {
                site[w + word] |= bit;
            }
            if code.slot_alt(1) {
                site[2 * w + word] |= bit;
            }
            if code.slot_missing(1) {
                site[3 * w + word] |= bit;
            }
            if self.mode.rule_bits
                && let Some(rule) = mask_rule.get(i).copied().filter(|&r| r != u16::MAX)
                && let Some(slot) = site.get_mut(4 * w + usize::from(rule) * w + word)
            {
                *slot |= bit;
            }
        }
        // The per-site effective mask bitmap the pruner excludes over
        // (doc_v7 D7-4), stored past the four genotype planes and the rule
        // bitmaps.
        if self.mode.mask_plane {
            let base = (4 + usize::from(self.mode.rule_bits) * self.plan.mask.len()) * w;
            for wi in 0..w {
                let m = masked.words().get(wi).copied().unwrap_or(0);
                if let Some(slot) = site.get_mut(base + wi) {
                    *slot = m;
                }
            }
        }
        // --mind (doc_v4 02 回退版 §3.2): the post-mask miss planes' moving
        // per-column popcounts. A masked sample contributes no called
        // alleles, a half call one, a full call two — exactly
        // `2·sites − m0 − m1` per sample once the unit folds.
        if self.want_counts {
            for wi in 0..w {
                let mask_word = masked.words().get(wi).copied().unwrap_or(0);
                let mut post0 = site[w + wi] | mask_word;
                while post0 != 0 {
                    let b = post0.trailing_zeros() as usize;
                    if let Some(slot) = m0.get_mut(wi * 64 + b) {
                        *slot += 1;
                    }
                    post0 &= post0 - 1;
                }
                let mut post1 = site[3 * w + wi] | mask_word;
                while post1 != 0 {
                    let b = post1.trailing_zeros() as usize;
                    if let Some(slot) = m1.get_mut(wi * 64 + b) {
                        *slot += 1;
                    }
                    post1 &= post1 - 1;
                }
            }
        }
        if self.mode.retain {
            if contigs.last().is_none_or(|(_, name)| name != contig) {
                contigs.push((pos0.len(), contig.to_vec()));
            }
            flat.extend_from_slice(site);
            pos0.push(pos1.saturating_sub(1));
            has_gt.push(record_has_gt);
        }
    }

    /// Scans one text record through the shared pipeline.
    fn scan_line(&self, line: &[u8], scratch: &mut MatrixScratch) -> Result<(), EngineError> {
        let plan = &self.plan;
        let opts = self.opts();
        let MatrixScratch {
            pipeline,
            funnel,
            site,
            m0,
            m1,
            flat,
            pos0,
            contigs,
            has_gt,
            ..
        } = scratch;
        let processed = pipeline.process(line, plan, &self.keep, opts, funnel, |view| {
            self.survivor_body(
                site,
                m0,
                m1,
                flat,
                pos0,
                contigs,
                has_gt,
                view.scan.codes,
                view.mask_rule,
                view.masked,
                view.site.contig_name,
                view.pos1,
                view.scan.has_gt,
            );
            Ok(())
        })?;
        self.record_head_bit(scratch, processed.head_passed);
        Ok(())
    }

    /// Scans one native BCF record (doc_v4 O-4): the wrapper hands the
    /// framed record over decoded.
    fn scan_bcf_record(
        &self,
        dict: &reed_io::bcf::BcfDict,
        bytes: &[u8],
        scratch: &mut MatrixScratch,
        native: &mut NativeScratch,
    ) -> Result<(), EngineError> {
        native.decode(
            dict,
            bytes,
            &self.int_keys,
            &self.float_keys,
            self.n_samples,
        )?;
        let plan = &self.plan;
        let opts = self.opts();
        let MatrixScratch {
            pipeline,
            funnel,
            site,
            m0,
            m1,
            flat,
            pos0,
            contigs,
            has_gt,
            ..
        } = scratch;
        let processed = pipeline.process_scanned(
            &native.site_line,
            plan,
            &self.keep,
            opts,
            funnel,
            crate::filter::record_pipeline::NativeScan {
                codes: &native.codes,
                int_columns: &native.int_columns,
                float_columns: &native.float_columns,
                counters: &native.counters,
                has_gt: native.has_gt,
            },
            |view| {
                self.survivor_body(
                    site,
                    m0,
                    m1,
                    flat,
                    pos0,
                    contigs,
                    has_gt,
                    view.scan.codes,
                    view.mask_rule,
                    view.masked,
                    view.site.contig_name,
                    view.pos1,
                    view.scan.has_gt,
                );
                Ok(())
            },
        )?;
        self.record_head_bit(scratch, processed.head_passed);
        Ok(())
    }

    /// Folds the unit's accumulators into a partial (`called = 2·sites −
    /// m0 − m1` per sample).
    fn take_partial(&self, scratch: &mut MatrixScratch, seam: bool) -> MatrixPartial {
        let mut called = vec![0_u64; self.n_samples];
        if self.want_counts {
            let sites = scratch.n_head.saturating_mul(2);
            for (i, slot) in called.iter_mut().enumerate() {
                *slot = sites
                    .saturating_sub(*scratch.m0.get(i).unwrap_or(&0))
                    .saturating_sub(*scratch.m1.get(i).unwrap_or(&0));
            }
        }
        MatrixPartial {
            seam,
            head_bits: std::mem::take(&mut scratch.head_bits),
            n_records: scratch.n_records,
            n_head: scratch.n_head,
            called,
            flat: std::mem::take(&mut scratch.flat),
            pos0: std::mem::take(&mut scratch.pos0),
            has_gt: std::mem::take(&mut scratch.has_gt),
            contigs: std::mem::take(&mut scratch.contigs),
            funnel: std::mem::take(&mut scratch.funnel),
        }
    }
    /// The native BCF record entry point (doc_v4 O-4): one chunk's framed
    /// records in, one chunk partial out.
    ///
    /// # Errors
    ///
    /// As [`Self::scan_bcf_record`].
    pub(crate) fn run_bcf_records(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        records: &[&[u8]],
        scratch: &mut Option<MatrixScratch>,
    ) -> Result<MatrixPartial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        self.reset_unit(scratch);
        let dict = self
            .bcf_input
            .as_ref()
            .expect("invariant: the native path needs the input dictionary");
        for bytes in records {
            let mut native = scratch
                .native
                .take()
                .expect("invariant: the native scratch exists with the dictionary");
            let result = self.scan_bcf_record(dict, bytes, scratch, &mut native);
            scratch.native = Some(native);
            result?;
        }
        Ok(self.take_partial(scratch, false))
    }
}

impl ChunkJob for MatrixJob {
    type Partial = MatrixPartial;
    type Scratch = Option<MatrixScratch>;

    fn run_chunk(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        self.reset_unit(scratch);
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            self.scan_line(line, scratch)?;
        }
        Ok(self.take_partial(scratch, false))
    }

    fn run_seam(
        &self,
        _id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        self.reset_unit(scratch);
        self.scan_line(line, scratch)?;
        Ok(self.take_partial(scratch, true))
    }
}

/// One contig's analysed sites, in position order, as flat site-major
/// blocks — the filter-side twin of `reed_stats::ContigMatrix` with the
/// sketch's `(4 + R) · w` stride (the rule bitmaps ride along when the mid
/// stage needs them).
#[derive(Clone, Debug, Default)]
pub(crate) struct SiteMatrix {
    /// Contig name bytes.
    pub contig: Vec<u8>,
    /// Words per plane row of a site (`4 + R` rows of `w` words).
    pub w: usize,
    /// Number of rule bitmaps (`R` in the stride; 0 without the mid stage).
    pub n_rules: usize,
    /// One effective mask bitmap per site rides past the rule bitmaps
    /// (the pruner's sample exclusion, doc_v7 D7-4).
    pub has_mask: bool,
    /// Flat site-major blocks, one per contributing stream unit.
    blocks: Vec<Vec<u64>>,
    /// First site index of each block; `base.len() == blocks.len() + 1`.
    base: Vec<usize>,
    /// Every site's 0-based position (ascending).
    pos0: Vec<u64>,
    /// Whether each site's record carried a GT field, aligned with `pos0`
    /// (the mid stage's `-g/--genotype` short-circuit; empty when the mode
    /// retains nothing).
    has_gt: Vec<bool>,
}

impl SiteMatrix {
    /// An empty matrix for one contig.
    pub fn new(contig: Vec<u8>, w: usize, n_rules: usize, has_mask: bool) -> Self {
        Self {
            contig,
            w,
            n_rules,
            has_mask,
            blocks: Vec::new(),
            base: vec![0],
            pos0: Vec::new(),
            has_gt: Vec::new(),
        }
    }

    /// Appends one stream unit's sites: `flat` is
    /// `positions.len() · (4 + R) · w` words.
    pub fn push_block(&mut self, flat: Vec<u64>, pos0: Vec<u64>, has_gt: Vec<bool>) {
        if pos0.is_empty() {
            return;
        }
        debug_assert_eq!(flat.len(), pos0.len() * self.stride());
        debug_assert_eq!(has_gt.len(), pos0.len());
        let n = self.pos0.len() + pos0.len();
        self.blocks.push(flat);
        self.pos0.extend(pos0);
        self.has_gt.extend(has_gt);
        self.base.push(n);
    }

    /// Words one site occupies.
    pub const fn stride(&self) -> usize {
        (4 + self.n_rules + self.has_mask as usize) * self.w
    }

    /// Sites in the matrix.
    pub fn n_sites(&self) -> usize {
        self.pos0.len()
    }

    /// Every site's 0-based position (ascending).
    pub fn pos0(&self) -> &[u64] {
        &self.pos0
    }

    /// Whether one site's record carried a GT field (`true` when the
    /// vector is absent — a matrix built without retention is never asked).
    pub fn has_gt(&self, local: usize) -> bool {
        self.has_gt.get(local).copied().unwrap_or(true)
    }

    /// Estimated resident size in bytes (the `--memory` judgment, Y1f).
    pub fn size_bytes(&self) -> u64 {
        let words: u64 = self.blocks.iter().map(|b| b.len() as u64).sum();
        words
            .saturating_mul(8)
            .saturating_add(self.pos0.len() as u64 * 8)
            .saturating_add(self.base.len() as u64 * 8)
    }

    fn block_of(&self, i: usize) -> usize {
        self.base.partition_point(|&b| b <= i) - 1
    }

    /// The four raw genotype plane slices of site `i`.
    ///
    /// # Panics
    ///
    /// Panics if `i >= n_sites()`.
    pub fn planes(&self, i: usize) -> [&[u64]; 4] {
        let b = self.block_of(i);
        let row = self.stride();
        let k = (i - self.base[b]) * row;
        let w = self.w;
        let blk = &self.blocks[b];
        [
            &blk[k..k + w],
            &blk[k + w..k + 2 * w],
            &blk[k + 2 * w..k + 3 * w],
            &blk[k + 3 * w..k + 4 * w],
        ]
    }

    /// The firing bitmap of mask rule `r` at site `i` (empty without
    /// `mode.rule_bits`).
    pub fn rule_mask(&self, i: usize, r: usize) -> &[u64] {
        let b = self.block_of(i);
        let row = self.stride();
        let k = (i - self.base[b]) * row + (4 + r) * self.w;
        &self.blocks[b][k..k + self.w]
    }

    /// The effective mask bitmap at site `i` (all-zero without
    /// `mode.mask_plane`): the samples the mask turned missing (doc_v7 D7-4).
    pub fn mask_plane(&self, i: usize) -> &[u64] {
        let b = self.block_of(i);
        let row = self.stride();
        let k = (i - self.base[b]) * row + (4 + self.n_rules) * self.w;
        &self.blocks[b][k..k + self.w]
    }
}

/// A chunk's stream bookkeeping for the mid stage: the head bitmap over
/// the chunk's records (body plus its seam) and the global site ordinals
/// the chunk's survivors occupy.
#[derive(Clone, Debug)]
pub(crate) struct ChunkSites {
    /// Chunk id.
    pub id: u64,
    /// Bit per record (body records first, then the seam record): head
    /// passed.
    pub head_bits: Vec<u64>,
    /// Records covered (the bitmap's population).
    pub n_records: usize,
    /// First global site ordinal of this chunk's sites.
    pub start: usize,
    /// Sites covered (body sites plus the seam site when it survived).
    pub n_sites: usize,
    /// The seam site's ordinal within the chunk (0-based among the
    /// chunk's sites); `None` when the seam record failed the head rules.
    pub seam_site: Option<usize>,
}

/// The open chunk the sink is still extending (its seam record may arrive
/// with the next chunk's partial).
struct OpenChunk {
    sites: ChunkSites,
}

/// Everything the pass0 hands back.
#[derive(Default)]
pub(crate) struct MatrixOutput {
    /// Per-sample post-mask called alleles over the head survivors.
    pub called: Vec<u64>,
    /// Head-surviving records (the `--mind` denominator).
    pub n_head: u64,
    /// The pass0 funnel (head hits, decode policies; the mid stage fills
    /// the tail/mask counters on the sketch path).
    pub funnel: FunnelPartial,
    /// Per-contig matrices in first-appearance order (empty when nothing
    /// was retained or the pruner consumed them).
    pub matrices: Vec<SiteMatrix>,
    /// First global site ordinal of each matrix, aligned with `matrices`.
    pub contig_base: Vec<usize>,
    /// Per-chunk stream bookkeeping, in chunk order (the mid stage's
    /// input; empty unless rule bitmaps were recorded).
    pub chunks: Vec<ChunkSites>,
    /// Total sites recorded.
    pub total_sites: usize,
    /// The per-contig prune results (the mode that prunes during the
    /// stream; `None` when the driver prunes afterwards).
    pub prune: Option<super::prune::PruneOutcome>,
    /// Pass0's engine metrics.
    pub metrics: reed_engine::metrics::StageMetrics,
    /// Chunks delivered (the timing report's chunk count).
    pub n_chunks: u64,
}

/// The sink: appends stream units in chunk order, folds the funnel and the
/// `--mind` counters, extends the per-contig matrices, and (when the mode
/// prunes during the stream) hands each closed contig to its own prune
/// thread so only one matrix is resident at a time.
struct MatrixSink {
    w: usize,
    n_rules: usize,
    /// Record the per-site effective mask bitmap (pruner exclusion).
    mask_plane: bool,
    want_counts: bool,
    n_samples: usize,
    /// Record the per-chunk stream bookkeeping (the mid stage's input;
    /// only the sketch shape reads it back).
    track_chunks: bool,
    /// The prune configuration (streaming prune mode only).
    prune_config: Option<super::prune::PruneConfig>,
    /// The memory budget the per-contig matrices are judged against.
    budget_bytes: u64,
    /// Upper bound on concurrently running per-contig prune threads (doc_v7
    /// D7-19: without a cap, M contigs keep M matrices resident).
    max_prune: usize,
    /// The final keep bitmap the streaming pruner excludes over (doc_v7 D7-4).
    keep_words: Arc<[u64]>,
    called: Vec<u64>,
    n_head: u64,
    funnel: FunnelPartial,
    matrices: Vec<SiteMatrix>,
    contig_base: Vec<usize>,
    total_sites: usize,
    pending: Option<SiteMatrix>,
    open: Option<OpenChunk>,
    chunks: Vec<ChunkSites>,
    prune_jobs: Vec<std::thread::JoinHandle<super::prune::PruneOutcome>>,
    /// Prune results already joined to respect `max_prune` (doc_v7 D7-19).
    prune_done: Vec<super::prune::PruneOutcome>,
    out: Arc<Mutex<Option<MatrixOutput>>>,
}

impl OrderedSink<MatrixPartial> for MatrixSink {
    fn accept(&mut self, id: ChunkId, partial: &mut MatrixPartial) -> Result<(), EngineError> {
        self.funnel.merge(&partial.funnel);
        if self.want_counts {
            if self.called.len() < partial.called.len() {
                self.called.resize(partial.called.len(), 0);
            }
            for (dst, src) in self.called.iter_mut().zip(&partial.called) {
                *dst += src;
            }
        }
        self.n_head += partial.n_head;
        if partial.seam {
            self.close_seam(partial)?;
        } else {
            self.open_chunk(id.0, partial)?;
        }
        Ok(())
    }

    fn finish(mut self) -> Result<u64, EngineError> {
        if self.track_chunks
            && let Some(open) = self.open.take()
        {
            self.chunks.push(open.sites);
        }
        self.close_pending()?;
        let mut prune = None;
        if self.prune_config.is_some() {
            let mut seen = 0_u64;
            let mut pruned = 0_u64;
            let mut per_contig = Vec::new();
            let started = std::time::Instant::now();
            for outcome in self.prune_done.drain(..) {
                seen += outcome.n_seen;
                pruned += outcome.n_pruned;
                per_contig.extend(outcome.keep.per_contig);
            }
            for handle in self.prune_jobs.drain(..) {
                let outcome = handle
                    .join()
                    .map_err(|_| EngineError::job("the per-contig prune worker panicked"))?;
                seen += outcome.n_seen;
                pruned += outcome.n_pruned;
                per_contig.extend(outcome.keep.per_contig);
            }
            prune = Some(super::prune::PruneOutcome {
                keep: super::prune::SiteKeepSet { per_contig },
                n_seen: seen,
                n_pruned: pruned,
                prune_secs: started.elapsed().as_secs_f64(),
            });
        }
        let out = MatrixOutput {
            called: std::mem::take(&mut self.called),
            n_head: self.n_head,
            funnel: std::mem::take(&mut self.funnel),
            matrices: std::mem::take(&mut self.matrices),
            contig_base: std::mem::take(&mut self.contig_base),
            chunks: std::mem::take(&mut self.chunks),
            total_sites: self.total_sites,
            prune,
            metrics: reed_engine::metrics::StageMetrics::default(),
            n_chunks: 0,
        };
        let report = out.total_sites as u64;
        *self
            .out
            .lock()
            .expect("invariant: matrix slot not poisoned") = Some(out);
        Ok(report)
    }
}

impl MatrixSink {
    /// Closes the previously opened chunk with this seam record: the seam
    /// belongs to the previous chunk's stream (§6.2).
    fn close_seam(&mut self, partial: &mut MatrixPartial) -> Result<(), EngineError> {
        if !self.track_chunks {
            // The seam site still joins the matrices; no bookkeeping.
            if partial.n_head > 0 {
                self.append_sites(partial)?;
            }
            return Ok(());
        }
        let bit = partial.head_bits.first().is_some_and(|word| word & 1 == 1);
        if let Some(mut open) = self.open.take() {
            let r = open.sites.n_records;
            let word = r / 64;
            if word >= open.sites.head_bits.len() {
                open.sites.head_bits.resize(word + 1, 0);
            }
            if bit {
                open.sites.head_bits[word] |= 1 << (r % 64);
            }
            open.sites.n_records += 1;
            if partial.n_head > 0 {
                self.append_sites(partial)?;
                open.sites.seam_site = Some(open.sites.n_sites);
                open.sites.n_sites += 1;
            }
            self.chunks.push(open.sites);
        }
        Ok(())
    }

    /// Opens a new chunk unit and appends its body sites.
    fn open_chunk(&mut self, id: u64, partial: &mut MatrixPartial) -> Result<(), EngineError> {
        if !self.track_chunks {
            return self.append_sites(partial);
        }
        // A stale open without a seam would mean a missing seam partial —
        // finalize defensively (the sketch sink did the same).
        if let Some(open) = self.open.take() {
            self.chunks.push(open.sites);
        }
        let start = self.total_sites;
        self.append_sites(partial)?;
        self.open = Some(OpenChunk {
            sites: ChunkSites {
                id,
                head_bits: std::mem::take(&mut partial.head_bits),
                n_records: partial.n_records,
                start,
                n_sites: partial.pos0.len(),
                seam_site: None,
            },
        });
        Ok(())
    }

    /// Appends a partial's sites to the per-contig matrices (the pending
    /// contig's matrix is closed when a new contig starts).
    fn append_sites(&mut self, partial: &mut MatrixPartial) -> Result<(), EngineError> {
        if partial.pos0.is_empty() {
            return Ok(());
        }
        for (run_index, (offset, contig)) in
            partial.contigs.iter().map(|(o, c)| (*o, c)).enumerate()
        {
            let end = partial
                .contigs
                .get(run_index + 1)
                .map_or(partial.pos0.len(), |(next, _)| *next);
            let same = self.pending.as_ref().is_some_and(|m| m.contig == *contig);
            if !same {
                self.close_pending()?;
                self.pending = Some(SiteMatrix::new(
                    contig.clone(),
                    self.w,
                    self.n_rules,
                    self.mask_plane,
                ));
                self.contig_base.push(self.total_sites);
            }
            let stride = self.pending.as_ref().map_or(1, SiteMatrix::stride);
            let flat = partial.flat[offset * stride..end * stride].to_vec();
            let pos0 = partial.pos0[offset..end].to_vec();
            let has_gt = partial.has_gt[offset..end].to_vec();
            self.total_sites += pos0.len();
            if let Some(matrix) = self.pending.as_mut() {
                matrix.push_block(flat, pos0, has_gt);
            }
        }
        Ok(())
    }

    /// Closes the pending contig's matrix: either parks it for the
    /// driver's post-pass analyses or hands it to its own prune thread
    /// (one contig resident at a time, the Y1f segmentation).
    fn close_pending(&mut self) -> Result<(), EngineError> {
        let Some(matrix) = self.pending.take() else {
            return Ok(());
        };
        if matrix.n_sites() == 0 {
            self.contig_base.pop();
            return Ok(());
        }
        let Some(config) = self.prune_config else {
            self.matrices.push(matrix);
            return Ok(());
        };
        // Streaming per-contig pruning (one matrix resident at a time)
        // only fits the mode whose tail ran inline; the sketch shape keeps
        // its matrices for the post-`--mind` mid stage.
        if self.track_chunks {
            self.matrices.push(matrix);
            return Ok(());
        }
        if matrix.size_bytes() > self.budget_bytes {
            return Err(EngineError::MemBudgetExceeded {
                peak_bytes: matrix.size_bytes(),
                budget_bytes: self.budget_bytes,
            });
        }
        // Bound the number of concurrently running prune threads so M contigs
        // cannot keep M matrices resident (doc_v7 D7-19): join the oldest
        // when the cap is reached, keeping its outcome for `finish`.
        if self.prune_jobs.len() >= self.max_prune {
            let handle = self.prune_jobs.remove(0);
            let outcome = handle
                .join()
                .map_err(|_| EngineError::job("the per-contig prune worker panicked"))?;
            self.prune_done.push(outcome);
        }
        let n_samples = self.n_samples;
        let keep_words = Arc::clone(&self.keep_words);
        self.prune_jobs.push(std::thread::spawn(move || {
            super::prune::prune_one(matrix, &config, None, 0, n_samples, &keep_words)
        }));
        Ok(())
    }
}

/// Runs the unified pass0 over an opened input.
///
/// # Errors
///
/// [`FilterError`] for transport, decode and job failures.
#[allow(clippy::too_many_arguments)]
pub(crate) fn run_matrix_pass(
    opened: crate::input::OpenedInput,
    plan: &FilterPlan,
    mode: MatrixMode,
    want_counts: bool,
    tail_n_samples: Option<u64>,
    keep: Arc<[bool]>,
    prune_config: Option<super::prune::PruneConfig>,
    budget_bytes: u64,
    opts: &FilterOptions,
    native_bcf: bool,
) -> Result<MatrixOutput, FilterError> {
    let budget = ThreadBudget {
        threads: opts.threads.max(1),
        mem: opts.mem_bytes.map_or_else(
            reed_engine::budget::MemBudget::default_half_memory,
            reed_engine::budget::MemBudget::from_bytes,
        ),
    };
    let n_samples = opened.facts.header.sample_names().len();
    let abort = reed_engine::abort::abort_for(&budget);
    let inner = MatrixJob {
        plan: Arc::new(plan.clone()),
        n_samples,
        w: n_samples.div_ceil(64).max(1),
        mode,
        want_counts,
        tail_n_samples,
        keep: Arc::clone(&keep),
        bcf_input: native_bcf.then(|| {
            opened
                .bcf
                .clone()
                .expect("invariant: a BCF input carries its dictionary")
        }),
        gt_table: native_bcf.then(|| Arc::new(reed_io::bcf::GtTable::new(plan.policies))),
        int_keys: Arc::from(plan.int_fields.clone().into_boxed_slice()),
        float_keys: Arc::from(plan.float_fields.clone().into_boxed_slice()),
    };
    let job = if native_bcf {
        crate::bcf_input::InputJob::with_bcf_records(
            inner,
            &opened,
            &abort,
            MatrixJob::run_bcf_records,
        )
    } else {
        crate::bcf_input::InputJob::new(inner, &opened, &abort)
    };
    let crate::input::OpenedInput { mut source, .. } = opened;
    let out: Arc<Mutex<Option<MatrixOutput>>> = Arc::new(Mutex::new(None));
    let w = n_samples.div_ceil(64).max(1);
    let mut keep_words = vec![0_u64; w];
    for (i, &k) in keep.iter().enumerate() {
        if k {
            keep_words[i / 64] |= 1 << (i % 64);
        }
    }
    let sink = MatrixSink {
        w,
        n_rules: usize::from(mode.rule_bits) * plan.mask.len(),
        mask_plane: mode.mask_plane,
        want_counts,
        n_samples,
        track_chunks: mode.rule_bits,
        prune_config,
        budget_bytes,
        max_prune: opts.threads.max(1),
        keep_words: keep_words.into(),
        called: Vec::new(),
        n_head: 0,
        funnel: FunnelPartial {
            head_hits: vec![0; plan.head.len()],
            tail_hits: vec![0; plan.tail.len()],
            mask_hits: vec![0; plan.mask.len()],
            ..FunnelPartial::default()
        },
        matrices: Vec::new(),
        contig_base: Vec::new(),
        total_sites: 0,
        pending: None,
        open: None,
        chunks: Vec::new(),
        prune_jobs: Vec::new(),
        prune_done: Vec::new(),
        out: Arc::clone(&out),
    };
    let report = run(source.as_mut(), &job, sink, &budget, &abort).map_err(ferr)?;
    job.finish(report.n_chunks, opts.regions.is_empty())
        .map_err(ferr)?;
    let mut output = out
        .lock()
        .expect("invariant: matrix slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the output");
    output.metrics = report.metrics;
    output.n_chunks = report.n_chunks;
    Ok(output)
}
