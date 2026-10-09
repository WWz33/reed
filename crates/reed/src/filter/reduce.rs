//! `--max-missing-sample` pass0 (doc_v4 02-cli-semantics.md §3.2): the
//! per-sample reduction that resolves the dynamic sample rule.
//!
//! vcftools semantics: a sample's call rate is counted **after** the
//! genotype mask and only over the records that pass the head rules (the
//! genotype-independent site rules); the numerator is its called alleles
//! (a half call counts one, a masked call none) and the denominator is
//! `2 × that record count`. This runs the same pipeline the single pass
//! runs (`RecordPipeline`), so the three strategies cannot drift.

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

/// Pass0 partial: per-sample called-allele counts after masking, over
/// head-surviving records.
#[derive(Default)]
pub(super) struct MissingPartial {
    /// Records that passed the head rules (the `2L` denominator).
    pub(super) n_records: u64,
    /// Per-sample called-allele counts (post-mask).
    pub(super) called: Vec<u64>,
    /// Pass0's engine metrics (doc_v4 03 A4-5: the driver sums the passes).
    pub(super) metrics: reed_engine::metrics::StageMetrics,
}

/// Worker scratch: the shared record pipeline plus the per-sample counters.
pub(super) struct MissingScratch {
    pipeline: RecordPipeline,
    called: Vec<u64>,
    funnel: FunnelPartial,
}

/// The GT-only call-rate job (`--max-missing-sample` pass0).
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
    // The tail rules are not evaluated in pass0 (the sample rule runs
    // before them).
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
            called: vec![0; self.n_samples],
            funnel: FunnelPartial {
                head_hits: vec![0; self.plan.head.len()],
                tail_hits: vec![0; self.plan.tail.len()],
                mask_hits: vec![0; self.plan.mask.len()],
                ..FunnelPartial::default()
            },
        }
    }

    /// Runs one record and returns whether it survived the head rules; the
    /// body adds the post-mask called alleles of such records.
    fn scan_line(
        &self,
        line: &[u8],
        scratch: &mut MissingScratch,
    ) -> Result<bool, EngineError> {
        let plan = &self.plan;
        let MissingScratch { pipeline, called, funnel } = scratch;
        let called = &mut *called;
        let processed = pipeline.process(line, plan, &[], OPTS, funnel, |view| {
            for (i, &code) in view.scan.codes.iter().enumerate() {
                let n = called_alleles(code, view.masked.get(i));
                if let Some(slot) = called.get_mut(i) {
                    *slot += n;
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
        scratch.called.clear();
        scratch.called.resize(self.n_samples, 0);
        let mut n_records = 0_u64;
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            if self.scan_line(line, scratch)? {
                n_records += 1;
            }
        }
        Ok(MissingPartial {
            n_records,
            called: scratch.called.clone(),
            metrics: reed_engine::metrics::StageMetrics::default(),
        })
    }

    fn run_seam(
        &self,
        _id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        scratch.called.clear();
        scratch.called.resize(self.n_samples, 0);
        let n_records = u64::from(self.scan_line(line, scratch)?);
        Ok(MissingPartial {
            n_records,
            called: scratch.called.clone(),
            metrics: reed_engine::metrics::StageMetrics::default(),
        })
    }
}

pub(super) fn run_missing_pass(
    opened: crate::input::OpenedInput,
    plan: &FilterPlan,
    opts: &FilterOptions,
) -> Result<MissingPartial, FilterError> {
    let budget = ThreadBudget {
        threads: opts.threads.max(1),
        mem: opts.mem_bytes.map_or_else(
            reed_engine::budget::MemBudget::default_half_memory,
            reed_engine::budget::MemBudget::from_bytes,
        ),
    };
    let n_samples = opened.facts.header.sample_names().len();
    let abort = reed_engine::abort::abort();
    let job = crate::bcf_input::InputJob::new(
        MissingJob { plan: Arc::new(plan.clone()), n_samples },
        &opened,
        &abort,
    );
    let crate::input::OpenedInput { mut source, .. } = opened;
    let total: Arc<Mutex<Option<MissingPartial>>> = Arc::new(Mutex::new(None));
    let sink = SumSink { total: Arc::clone(&total) };
    let report = run(source.as_mut(), &job, sink, &budget, &abort).map_err(ferr)?;
    job.finish(report.n_chunks, opts.regions.is_empty()).map_err(ferr)?;
    let mut partial = total
        .lock()
        .expect("invariant: missing slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the partial");
    partial.metrics = report.metrics;
    Ok(partial)
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
                if total.called.len() < partial.called.len() {
                    total.called.resize(partial.called.len(), 0);
                }
                for (t, m) in total.called.iter_mut().zip(&partial.called) {
                    *t += m;
                }
            }
            None => {
                *guard = Some(MissingPartial {
                    n_records: partial.n_records,
                    called: partial.called.clone(),
                    metrics: reed_engine::metrics::StageMetrics::default(),
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

/// Post-mask called alleles of one sample's genotype (doc_v4 02 §3.1: a
/// half call contributes one allele, a masked call none).
pub(super) fn called_alleles(code: reed_core::genotype::GtCode, masked: bool) -> u64 {
    if masked || !code.is_called_any() {
        return 0;
    }
    2 - u64::from(code.slot_missing(0)) - u64::from(code.slot_missing(1))
}
