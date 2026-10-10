//! `--indep-pairwise` pass0 (docs/rounds/v5/05-algorithms.md §3.4.1): the
//! site-level pruning rule. Workers decode each head-and-tail-surviving
//! record's masked genotype planes and its major-allele frequency; the
//! reducer feeds each contig's sites, in order, to one [`LdPruner`] and
//! returns the global ordinals to keep. Pass1 then emits only those sites.
//!
//! The pruning is order-dependent **within a contig** and independent
//! **between contigs** (plink2's subcontig parallelism), so the reducer only
//! collects each chunk's surviving sites in arrival order and the per-contig
//! pruners run in parallel at the end (`rayon`). Before Y1 one pruner ran on
//! the reduce thread for the whole file: 66 s pass0, 64 s of it in the
//! reducer (docs/rounds/v6/04-performance.md §6). Each pruner is
//! contig-agnostic and keeps only its window's sites resident (P5-5).

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
    /// Wall time of the parallel per-contig prune (Y1 measurement).
    pub(super) prune_secs: f64,
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
    let collected: Arc<Mutex<Option<Vec<PruneContig>>>> = Arc::new(Mutex::new(None));
    let sink = LdPruneSink {
        contigs: Vec::new(),
        n_survivors: 0,
        out: Arc::clone(&collected),
    };
    let report = run(source.as_mut(), &job, sink, &budget, &abort).map_err(ferr)?;
    job.finish(report.n_chunks, opts.regions.is_empty())
        .map_err(ferr)?;
    let contigs = collected
        .lock()
        .expect("invariant: prune slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the sites");
    let prune_start = std::time::Instant::now();
    let mut result = prune_contigs(contigs, config, n_samples.div_ceil(64).max(1), n_samples);
    result.prune_secs = prune_start.elapsed().as_secs_f64();
    result.metrics = report.metrics;
    Ok(result)
}

/// The reducer: collects every surviving site in arrival order, grouped by
/// contig. No pruning happens here — the per-contig pruners run in parallel
/// once the stream is exhausted (`prune_contigs`), which is the Y1 fix for
/// the 64 s serial reducer (docs/rounds/v6/04-performance.md §6).
struct LdPruneSink {
    /// Surviving sites per contig, in arrival order.
    contigs: Vec<PruneContig>,
    n_survivors: u64,
    out: Arc<Mutex<Option<Vec<PruneContig>>>>,
}

/// One contig's surviving sites, in arrival order.
struct PruneContig {
    name: Vec<u8>,
    sites: Vec<PruneSite>,
}

impl OrderedSink<LdPrunePartial> for LdPruneSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut LdPrunePartial) -> Result<(), EngineError> {
        for site in partial.sites.drain(..) {
            self.n_survivors += 1;
            match self.contigs.last_mut() {
                Some(cur) if cur.name == site.contig => cur.sites.push(site),
                _ => {
                    let name = site.contig.clone();
                    self.contigs.push(PruneContig {
                        name,
                        sites: vec![site],
                    });
                }
            }
        }
        Ok(())
    }

    fn finish(self) -> Result<u64, EngineError> {
        let n = self.n_survivors;
        *self.out.lock().expect("invariant: prune slot not poisoned") = Some(self.contigs);
        Ok(n)
    }
}

/// Prunes every contig in parallel (plink2's subcontig parallelism) and
/// returns the kept sites plus the counters.
///
/// Each contig gets its own [`LdPruner`] fed in arrival order, so the
/// within-contig decision sequence is exactly the serial one; contigs share
/// nothing, so the parallel result is deterministic.
fn prune_contigs(
    contigs: Vec<PruneContig>,
    config: PruneConfig,
    w: usize,
    n_samples: usize,
) -> LdPruneResult {
    use rayon::prelude::*;
    let n_survivors: u64 = contigs.iter().map(|c| c.sites.len() as u64).sum();
    let per_contig: Vec<(Vec<u8>, Vec<u64>)> = contigs
        .into_par_iter()
        .map(|c| {
            let mut pruner =
                LdPruner::new(config.threshold, config.window, config.step, w, n_samples);
            for (i, site) in c.sites.iter().enumerate() {
                let planes: [&[u64]; 4] = [
                    &site.planes[0..w],
                    &site.planes[w..2 * w],
                    &site.planes[2 * w..3 * w],
                    &site.planes[3 * w..4 * w],
                ];
                pruner.push(i, site.pos0, planes, site.stats);
            }
            pruner.finish_contig();
            let kept: Vec<u64> = pruner
                .kept_sites()
                .into_iter()
                .filter_map(|ordinal| c.sites.get(ordinal).map(|s| s.pos0))
                .collect();
            (c.name, kept)
        })
        .collect();
    let n_kept: u64 = per_contig.iter().map(|(_, p)| p.len() as u64).sum();
    LdPruneResult {
        keep: SiteKeepSet { per_contig },
        n_survivors,
        n_pruned: n_survivors - n_kept,
        metrics: reed_engine::metrics::StageMetrics::default(),
        prune_secs: 0.0,
    }
}
