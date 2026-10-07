//! `--mind` pass0 (04 §C0, `doc_v2/03-architecture.md` §8): the per-sample
//! missing reduction that resolves the dynamic sample rule.
//!
//! plink2 semantics: a sample's missing rate is counted **after** the
//! genotype mask and only over the records that pass the head rules (the
//! genotype-independent site rules); the denominator is that record count.
//! This runs the same pipeline the single pass runs (`RecordPipeline`), so
//! the three strategies cannot drift.

use std::sync::{Arc, Mutex};

use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{run, ChunkJob, OrderedSink};
use reed_filter::funnel::FunnelPartial;
use reed_filter::plan::FilterPlan;
use reed_filter::FilterError;

use reed_core::ChunkId;

use super::record_pipeline::{PipelineOptions, RecordPipeline};
use super::{ferr, FilterOptions};

/// Pass0 partial: per-sample fully-missing GT counts after masking, over
/// head-surviving records.
#[derive(Default)]
pub(super) struct MissingPartial {
    /// Records that passed the head rules (the `--mind` denominator).
    pub(super) n_records: u64,
    /// Per-sample fully-missing genotype counts (post-mask).
    pub(super) counts: Vec<u64>,
}

/// Worker scratch: the shared record pipeline plus the per-sample counters.
pub(super) struct MissingScratch {
    pipeline: RecordPipeline,
    miss: Vec<u64>,
    funnel: FunnelPartial,
}

/// The GT-only missing-count job (`--mind` pass0).
pub(super) struct MissingJob {
    pub(super) plan: Arc<FilterPlan>,
    pub(super) n_samples: usize,
}

const OPTS: PipelineOptions = PipelineOptions {
    // The denominator is the head-surviving record count: records that fail
    // the head rules never reach the cohort.
    scan_head_failures: false,
    // The mask runs for every sample (the keep set is not decided yet).
    mask_all_samples: true,
    require_contig: false,
    // The tail rules are not evaluated in pass0 (`--mind` runs before them).
    tail_n_samples: None,
    counts: false,
};

impl MissingJob {
    fn build_scratch(&self) -> MissingScratch {
        // The mask rules read the plan's FORMAT columns, so pass0 scans the
        // same fields the single pass does (only GT is otherwise needed).
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
        MissingScratch {
            pipeline,
            miss: vec![0; self.n_samples],
            funnel: FunnelPartial {
                head_hits: vec![0; self.plan.head.len()],
                tail_hits: vec![0; self.plan.tail.len()],
                mask_hits: vec![0; self.plan.mask.len()],
                ..FunnelPartial::default()
            },
        }
    }

    /// Runs one record and returns whether it survived the head rules; the
    /// body counts the post-mask fully-missing calls of such records.
    fn scan_line(
        &self,
        line: &[u8],
        scratch: &mut MissingScratch,
    ) -> Result<bool, EngineError> {
        let plan = &self.plan;
        let MissingScratch { pipeline, miss, funnel } = scratch;
        let miss = &mut *miss;
        let processed = pipeline.process(line, plan, &[], OPTS, funnel, |view| {
            for (i, &code) in view.scan.codes.iter().enumerate() {
                if (code.is_fully_missing() || view.masked.get(i))
                    && let Some(slot) = miss.get_mut(i)
                {
                    *slot += 1;
                }
            }
            Ok(())
        })?;
        Ok(processed.head_passed)
    }
}

impl ChunkJob for MissingJob {
    type Partial = MissingPartial;
    type Scratch = Option<MissingScratch>;

    fn run_chunk(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        scratch.miss.clear();
        scratch.miss.resize(self.n_samples, 0);
        let mut n_records = 0_u64;
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            if self.scan_line(line, scratch)? {
                n_records += 1;
            }
        }
        Ok(MissingPartial { n_records, counts: scratch.miss.clone() })
    }

    fn run_seam(
        &self,
        _id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        scratch.miss.clear();
        scratch.miss.resize(self.n_samples, 0);
        let n_records = u64::from(self.scan_line(line, scratch)?);
        Ok(MissingPartial { n_records, counts: scratch.miss.clone() })
    }
}

pub(super) fn run_missing_pass(
    opened: crate::input::OpenedInput,
    plan: &FilterPlan,
    opts: &FilterOptions,
) -> Result<MissingPartial, FilterError> {
    let budget = ThreadBudget { threads: opts.threads.max(1), mem: reed_engine::budget::MemBudget::default_half_memory() };
    let crate::input::OpenedInput { mut source, facts } = opened;
    let job = MissingJob {
        plan: Arc::new(plan.clone()),
        n_samples: facts.header.sample_names().len(),
    };
    let total: Arc<Mutex<Option<MissingPartial>>> = Arc::new(Mutex::new(None));
    let sink = SumSink { total: Arc::clone(&total) };
    let abort = reed_engine::abort::abort();
    run(source.as_mut(), &job, sink, &budget, &abort).map_err(ferr)?;
    Ok(total
        .lock()
        .expect("invariant: missing slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the partial"))
}

/// Sums pass0 partials (integer sums are order-independent).
pub(super) struct SumSink {
    total: Arc<Mutex<Option<MissingPartial>>>,
}

impl OrderedSink<MissingPartial> for SumSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut MissingPartial) -> Result<(), EngineError> {
        let mut guard = self
            .total
            .lock()
            .expect("invariant: missing slot not poisoned");
        match guard.as_mut() {
            Some(total) => {
                total.n_records += partial.n_records;
                if total.counts.len() < partial.counts.len() {
                    total.counts.resize(partial.counts.len(), 0);
                }
                for (t, m) in total.counts.iter_mut().zip(&partial.counts) {
                    *t += m;
                }
            }
            None => {
                *guard = Some(MissingPartial {
                    n_records: partial.n_records,
                    counts: partial.counts.clone(),
                });
            }
        }
        Ok(())
    }

    fn finish(self) -> Result<u64, EngineError> {
        let total = self.total.lock().expect("invariant: missing slot not poisoned");
        Ok(total.as_ref().map_or(0, |t| t.n_records))
    }
}
