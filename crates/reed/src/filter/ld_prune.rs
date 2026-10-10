//! `--indep-pairwise` pass0 (docs/rounds/v5/05-algorithms.md §3.4.1): the
//! site-level pruning rule. Workers decode each head-and-tail-surviving
//! record's masked genotype planes and its major-allele frequency; the
//! reducer feeds each contig's sites, in order, to one [`LdPruner`] and
//! returns the global ordinals to keep. Pass1 then emits only those sites.
//!
//! The pruning is order-dependent, so it cannot run per chunk. Today one
//! [`LdPruner`] runs inside the ordered reducer (`LdPruneSink`), i.e. all
//! contigs are pruned serially on the reduce thread — the 66 s pass0 cost
//! measured in docs/rounds/v6/04-performance.md §6. Per-contig parallel
//! pruning is the Y1 matrix work (docs/rounds/v6/06-roadmap.md); the pruner
//! itself is already contig-agnostic. Only the current window's sites are
//! resident (P5-5: memory does not grow with contig length).

use std::sync::{Arc, Mutex};

use reed_core::genotype::PlaneBlock;
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{ChunkJob, OrderedSink, run};
use reed_filter::FilterError;
use reed_filter::funnel::FunnelPartial;
use reed_filter::plan::FilterPlan;
use reed_stats::ld::{LdPruner, LdWindow};

use reed_core::ChunkId;

use super::record_pipeline::{PipelineOptions, RecordPipeline};
use super::{FilterOptions, ferr};

/// One surviving site handed to the reducer.
pub(super) struct PruneSite {
    /// Contig name bytes.
    pub(super) contig: Vec<u8>,
    /// 0-based position.
    pub(super) pos0: u64,
    /// The masked genotype planes (four `w`-word slices, flattened).
    pub(super) planes: Vec<u64>,
    /// plink2's per-site aggregates (`--freq` major frequency and the
    /// homozygote/heterozygote counts its `FillVaggs` produces).
    pub(super) stats: reed_stats::ld::SitePruneStats,
}

/// Pass0 partial: the chunk's surviving sites, in order.
#[derive(Default)]
pub(super) struct LdPrunePartial {
    /// The chunk's sites, in record order.
    pub(super) sites: Vec<PruneSite>,
    /// Sites in this chunk.
    pub(super) n_survivors: u64,
}

/// The pass0 job: head → mask → tail survival, then planes + major freq.
pub(super) struct LdPruneJob {
    pub(super) plan: Arc<FilterPlan>,
    pub(super) n_samples: usize,
    pub(super) keep: Arc<[bool]>,
}

/// Worker scratch.
pub(super) struct LdPruneScratch {
    pipeline: RecordPipeline,
    funnel: FunnelPartial,
    planes: PlaneBlock,
}

/// The same survival contract pass1 uses, with the tail rules evaluated
/// inline (the pruning happens after every site rule, per 02 §6 step 6b).
const OPTS: PipelineOptions = PipelineOptions {
    scan_head_failures: false,
    mask_all_samples: false,
    require_contig: false,
    // The tail rules decide survival; `--mind`'s keep set is already fixed.
    tail_n_samples: None,
    counts: true,
};

impl LdPruneJob {
    fn build_scratch(&self) -> LdPruneScratch {
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
        LdPruneScratch {
            pipeline,
            funnel: FunnelPartial {
                head_hits: vec![0; self.plan.head.len()],
                tail_hits: vec![0; self.plan.tail.len()],
                mask_hits: vec![0; self.plan.mask.len()],
                ..FunnelPartial::default()
            },
            planes: PlaneBlock::new(1, self.n_samples, false),
        }
    }

    /// Scans one record; appends the surviving site to `out`.
    fn scan_line(
        &self,
        line: &[u8],
        scratch: &mut LdPruneScratch,
        out: &mut LdPrunePartial,
    ) -> Result<(), EngineError> {
        let plan = &self.plan;
        let keep = &self.keep;
        let LdPruneScratch {
            pipeline,
            funnel,
            planes,
            ..
        } = scratch;
        let mut site = None;
        let processed = pipeline.process(line, plan, keep, OPTS, funnel, |view| {
            planes.push_codes(0, view.scan.codes);
            let sp = planes.site_planes(0);
            let mut flat = Vec::with_capacity(sp.slot0_alt.len() * 4);
            flat.extend_from_slice(sp.slot0_alt);
            flat.extend_from_slice(sp.slot0_miss);
            flat.extend_from_slice(sp.slot1_alt);
            flat.extend_from_slice(sp.slot1_miss);
            let n = sp.slot0_alt.len();
            let stats = {
                let planes4 = [
                    &flat[0..n],
                    &flat[n..2 * n],
                    &flat[2 * n..3 * n],
                    &flat[3 * n..4 * n],
                ];
                site_prune_stats(planes4, self.n_samples)
            };
            site = Some(PruneSite {
                contig: view.site.contig_name.to_vec(),
                pos0: view.pos1.saturating_sub(1),
                planes: flat,
                stats,
            });
            Ok(())
        })?;
        // `emitted` is exactly "head passed and the inline tail passed".
        if processed.emitted
            && let Some(site) = site
        {
            out.n_survivors += 1;
            out.sites.push(site);
        }
        Ok(())
    }
}

impl ChunkJob for LdPruneJob {
    type Partial = LdPrunePartial;
    type Scratch = Option<LdPruneScratch>;

    fn run_chunk(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        let mut out = LdPrunePartial::default();
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            self.scan_line(line, scratch, &mut out)?;
        }
        Ok(out)
    }

    fn run_seam(
        &self,
        _id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        let mut out = LdPrunePartial::default();
        self.scan_line(line, scratch, &mut out)?;
        Ok(out)
    }
}

/// One site's plink2 aggregates over all samples, from its four planes:
/// the major-allele frequency and the hom-ref / hom-alt / het counts
/// (`FillVaggs`).
pub(super) fn site_prune_stats(
    planes: [&[u64]; 4],
    n_samples: usize,
) -> reed_stats::ld::SitePruneStats {
    let mut alt = 0_u64;
    let mut called = 0_u64;
    let mut n_hom_ref = 0_u64;
    let mut n_hom_alt = 0_u64;
    let mut n_het = 0_u64;
    for (w, _) in planes[0].iter().enumerate() {
        let real = if (w + 1) * 64 <= n_samples {
            !0_u64
        } else {
            (1_u64 << (n_samples % 64)) - 1
        };
        let miss = (planes[1][w] | planes[3][w]) & real;
        let called_w = real & !miss;
        let alt0 = planes[0][w] & called_w;
        let alt1 = planes[2][w] & called_w;
        alt += alt0.count_ones() as u64 + alt1.count_ones() as u64;
        called += 2 * called_w.count_ones() as u64;
        // Both slots called and agreeing / disagreeing.
        let both_called = called_w & !(planes[1][w] ^ planes[3][w]) & called_w;
        let _ = both_called;
        let hom = called_w & !(planes[1][w] | planes[3][w]) & !(alt0 ^ alt1);
        let hom_alt = hom & alt0;
        let hom_ref = hom & !alt0;
        n_hom_alt += hom_alt.count_ones() as u64;
        n_hom_ref += hom_ref.count_ones() as u64;
        n_het += (called_w & !(planes[1][w] | planes[3][w]) & (alt0 ^ alt1)).count_ones() as u64;
    }
    let maj_freq = if called == 0 {
        1.0
    } else {
        let alt_freq = alt as f64 / called as f64;
        alt_freq.max(1.0 - alt_freq)
    };
    reed_stats::ld::SitePruneStats {
        maj_freq,
        n_hom_ref,
        n_hom_alt,
        n_het,
    }
}

/// The prune configuration the sink needs.
#[derive(Clone, Copy, Debug)]
pub(super) struct PruneConfig {
    pub(super) threshold: f64,
    pub(super) window: LdWindow,
    pub(super) step: usize,
}

/// The kept sites, per contig, sorted by position (pass1 binary-searches
/// this: worker chunks are unordered, so an ordinal bitmap cannot be used).
#[derive(Clone, Debug, Default)]
pub(super) struct SiteKeepSet {
    /// `(contig, sorted 0-based positions)`.
    per_contig: Vec<(Vec<u8>, Vec<u64>)>,
}

impl SiteKeepSet {
    /// Whether the site at `(contig, pos0)` was kept.
    pub(super) fn contains(&self, contig: &[u8], pos0: u64) -> bool {
        self.per_contig
            .iter()
            .find(|(c, _)| c == contig)
            .is_some_and(|(_, positions)| positions.binary_search(&pos0).is_ok())
    }
}

/// The pass0 result: the kept sites and the counters.
pub(super) struct LdPruneResult {
    /// The kept sites.
    pub(super) keep: SiteKeepSet,
    /// Sites that survived the filter rules.
    pub(super) n_survivors: u64,
    /// Sites the pruner removed (plus monomorphic ones).
    pub(super) n_pruned: u64,
    /// Pass0's engine metrics.
    pub(super) metrics: reed_engine::metrics::StageMetrics,
}

/// Runs the pruning pass0 over an opened input.
///
/// # Errors
///
/// [`FilterError`] for transport, decode and job failures.
pub(super) fn run_ld_prune_pass(
    opened: crate::input::OpenedInput,
    plan: &FilterPlan,
    keep: &Arc<[bool]>,
    config: PruneConfig,
    opts: &FilterOptions,
    native_bcf: bool,
) -> Result<LdPruneResult, FilterError> {
    let budget = ThreadBudget {
        threads: opts.threads.max(1),
        mem: opts.mem_bytes.map_or_else(
            reed_engine::budget::MemBudget::default_half_memory,
            reed_engine::budget::MemBudget::from_bytes,
        ),
    };
    let n_samples = opened.facts.header.sample_names().len();
    let abort = reed_engine::abort::abort();
    let inner = LdPruneJob {
        plan: Arc::new(plan.clone()),
        n_samples,
        keep: Arc::clone(keep),
    };
    // BCF input renders through the text path (the native prune path is not
    // wired yet); `native_bcf` is accepted for signature parity.
    let _ = native_bcf;
    let job = crate::bcf_input::InputJob::new(inner, &opened, &abort);
    let crate::input::OpenedInput { mut source, .. } = opened;
    let total: Arc<Mutex<Option<LdPruneResult>>> = Arc::new(Mutex::new(None));
    let sink = LdPruneSink {
        base: 0,
        pruner: None,
        current_contig: Vec::new(),
        keys: Vec::new(),
        kept: Vec::new(),
        n_survivors: 0,
        n_pruned: 0,
        config,
        w: n_samples.div_ceil(64),
        n_samples,
        out: Arc::clone(&total),
    };
    let report = run(source.as_mut(), &job, sink, &budget, &abort).map_err(ferr)?;
    job.finish(report.n_chunks, opts.regions.is_empty())
        .map_err(ferr)?;
    let mut result = total
        .lock()
        .expect("invariant: prune slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the result");
    result.metrics = report.metrics;
    Ok(result)
}

/// The reducer: feeds each contig's sites, in order, to one pruner.
struct LdPruneSink {
    /// Ordinal base of the sites seen so far.
    base: u64,
    pruner: Option<LdPruner>,
    /// The contig the current pruner is working on.
    current_contig: Vec<u8>,
    /// Every surviving site's key, in arrival order (parallel to the
    /// pruner's ordinals).
    keys: Vec<(Vec<u8>, u64)>,
    /// The kept keys, filled at each contig's end.
    kept: Vec<(Vec<u8>, u64)>,
    n_survivors: u64,
    n_pruned: u64,
    config: PruneConfig,
    w: usize,
    n_samples: usize,
    out: Arc<Mutex<Option<LdPruneResult>>>,
}

impl OrderedSink<LdPrunePartial> for LdPruneSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut LdPrunePartial) -> Result<(), EngineError> {
        for site in partial.sites.drain(..) {
            // A new contig: close the previous one's window.
            if site.contig != self.current_contig {
                if let Some(mut pruner) = self.pruner.take() {
                    pruner.finish_contig();
                    self.absorb(&pruner);
                }
                self.current_contig.clone_from(&site.contig);
                self.pruner = Some(LdPruner::new(
                    self.config.threshold,
                    self.config.window,
                    self.config.step,
                    self.w,
                    self.n_samples,
                ));
            }
            let pruner = self.pruner.as_mut().expect("invariant: just set");
            let planes: [&[u64]; 4] = [
                &site.planes[0..self.w],
                &site.planes[self.w..2 * self.w],
                &site.planes[2 * self.w..3 * self.w],
                &site.planes[3 * self.w..4 * self.w],
            ];
            let ordinal = usize::try_from(self.n_survivors).unwrap_or(usize::MAX);
            pruner.push(ordinal, site.pos0, planes, site.stats);
            self.keys.push((site.contig, site.pos0));
            self.n_survivors += 1;
        }
        self.base += partial.n_survivors;
        Ok(())
    }

    fn finish(mut self) -> Result<u64, EngineError> {
        if let Some(mut pruner) = self.pruner.take() {
            pruner.finish_contig();
            self.absorb(&pruner);
        }
        self.n_pruned = self.n_survivors - self.kept.len() as u64;
        // Group by contig; the keys arrive in contig/position order, so each
        // contig's positions are already ascending.
        let mut per_contig: Vec<(Vec<u8>, Vec<u64>)> = Vec::new();
        for (contig, pos0) in &self.kept {
            match per_contig.last_mut() {
                Some((c, positions)) if c == contig => positions.push(*pos0),
                _ => per_contig.push((contig.clone(), vec![*pos0])),
            }
        }
        *self.out.lock().expect("invariant: prune slot not poisoned") = Some(LdPruneResult {
            keep: SiteKeepSet { per_contig },
            n_survivors: self.n_survivors,
            n_pruned: self.n_pruned,
            metrics: reed_engine::metrics::StageMetrics::default(),
        });
        Ok(self.n_survivors)
    }
}

impl LdPruneSink {
    /// Records one closed contig's keep decisions.
    fn absorb(&mut self, pruner: &LdPruner) {
        for ordinal in pruner.kept_sites() {
            if let Some(key) = self.keys.get(ordinal) {
                self.kept.push(key.clone());
            }
        }
    }
}
