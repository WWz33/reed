//! A filtering wrapper job (04-cleanup §2.1): applies the compiled filter
//! plan to the text stream before an inner consumer sees it, so `reed stats`
//! can accept the filter's site/threshold flags in the *same pass* instead
//! of running the filter API into an intermediate file.
//!
//! The wrapper runs the filter's own [`RecordPipeline`] (head → mask →
//! tail) and rewrites the surviving records with the filter's
//! [`edit_vcf_line`]: masked genotypes become `./.`, exactly what the inner
//! consumer would see after a physical `reed filter` + `reed stats`. INFO
//! updates stay filter-only, so the edit plan never touches INFO.

use std::sync::Arc;

use reed_core::ChunkId;
use reed_engine::error::EngineError;
use reed_engine::pipeline::ChunkJob;
use reed_filter::funnel::FunnelPartial;
use reed_filter::plan::FilterPlan;
use reed_io::chunk::ChunkBytes;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::line_edit::{edit_vcf_line, EditPlan};

use crate::filter::record_pipeline::{PipelineOptions, RecordPipeline};

/// Wraps an inner job with the filter's per-record stages; `plan = None`
/// passes every record through untouched (no buffering, no pipeline).
pub struct FilteredJob<J> {
    plan: Option<Arc<FilterPlan>>,
    /// Samples the filter keeps (all when no static subset was requested).
    keep: Arc<[bool]>,
    /// The retained sample count, the tail rules' missing-rate denominator.
    n_out: u64,
    inner: J,
}

/// Worker scratch: the record pipeline, the rewritten-record buffer and the
/// inner job's scratch. The pipeline is built lazily (it needs the plan).
pub struct FilteredScratch<S> {
    pipeline: Option<RecordPipeline>,
    /// Survivors, rewritten and newline-terminated (the inner job's body).
    buffer: Vec<u8>,
    funnel: FunnelPartial,
    inner: S,
}

impl<S: Default> Default for FilteredScratch<S> {
    fn default() -> Self {
        Self {
            pipeline: None,
            buffer: Vec::new(),
            funnel: FunnelPartial::default(),
            inner: S::default(),
        }
    }
}

impl<J: ChunkJob> FilteredJob<J> {
    /// Wraps `inner` with `plan` (the same compiler `reed filter` uses) over
    /// the sample keep mask.
    #[must_use]
    pub fn new(plan: Option<Arc<FilterPlan>>, keep: Arc<[bool]>, inner: J) -> Self {
        let n_out = u64::try_from(keep.iter().filter(|k| **k).count()).unwrap_or(u64::MAX);
        Self { plan, keep, n_out, inner }
    }

    fn build_scratch(&self) -> FilteredScratch<J::Scratch> {
        let Some(plan) = &self.plan else {
            return FilteredScratch::default();
        };
        let int_keys: Vec<&[u8]> = plan.int_fields.iter().map(Vec::as_slice).collect();
        let float_keys: Vec<&[u8]> = plan.float_fields.iter().map(Vec::as_slice).collect();
        let mut pipeline = RecordPipeline::new(
            self.keep.len(),
            &int_keys,
            &float_keys,
            plan.policies,
            &[],
        );
        pipeline.resize_for(plan);
        FilteredScratch {
            pipeline: Some(pipeline),
            buffer: Vec::with_capacity(1 << 16),
            funnel: FunnelPartial::default(),
            inner: J::Scratch::default(),
        }
    }

    /// Applies the filter plan to one record; a survivor is appended to
    /// `buffer` with its masked genotypes rewritten.
    fn filter_line(
        &self,
        line: &[u8],
        scratch: &mut FilteredScratch<J::Scratch>,
    ) -> Result<(), EngineError> {
        let plan = self.plan.as_ref().expect("invariant: guarded by the caller");
        let FilteredScratch { pipeline, buffer, funnel, .. } = scratch;
        let pipeline = pipeline.as_mut().expect("invariant: start_chunk built it");
        let opts = PipelineOptions {
            scan_head_failures: false,
            mask_all_samples: true,
            require_contig: false,
            // The tail's missing-rate denominator is the retained sample set.
            tail_n_samples: Some(self.n_out),
            counts: true,
        };
        pipeline.process(line, plan, &self.keep, opts, funnel, |view| {
            let edit = EditPlan {
                // Every column stays: the inner consumer selects samples
                // through its own group assignment.
                keep_samples: None,
                mask_gt: Some(view.masked),
                // INFO updates are filter-only.
                info_remove: &[],
                info_set: &[],
            };
            let gt_index = view
                .line
                .format_text()
                .split(|&b| b == b':')
                .position(|k| k == b"GT")
                .and_then(|i| u8::try_from(i).ok());
            edit_vcf_line(view.line, &edit, gt_index, buffer).map_err(EngineError::from)
        })?;
        Ok(())
    }

    fn start_chunk(&self, scratch: &mut FilteredScratch<J::Scratch>) {
        let Some(plan) = &self.plan else { return };
        if scratch.pipeline.is_none() {
            scratch.pipeline = Some(self.build_scratch().pipeline.expect("invariant: built"));
            scratch.buffer = Vec::with_capacity(1 << 16);
        }
        scratch.funnel.reset(plan.head.len(), plan.tail.len(), plan.mask.len());
        scratch.buffer.clear();
    }
}

impl<J: ChunkJob> ChunkJob for FilteredJob<J> {
    type Partial = J::Partial;
    type Scratch = FilteredScratch<J::Scratch>;

    fn run_chunk(
        &self,
        chunk: &ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        if self.plan.is_none() {
            return self.inner.run_chunk(chunk, body, &mut scratch.inner);
        }
        self.start_chunk(scratch);
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            self.filter_line(line, scratch)?;
        }
        let FilteredScratch { buffer, inner, .. } = scratch;
        self.inner.run_chunk(chunk, buffer, inner)
    }

    fn run_seam(
        &self,
        id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        if self.plan.is_none() {
            return self.inner.run_seam(id, line, &mut scratch.inner);
        }
        self.start_chunk(scratch);
        if LineRef::parse(line).is_none() {
            return Ok(J::Partial::default());
        }
        self.filter_line(line, scratch)?;
        let FilteredScratch { buffer, inner, .. } = scratch;
        if buffer.is_empty() {
            // The filter dropped the seam record.
            return Ok(J::Partial::default());
        }
        self.inner.run_seam(id, buffer, inner)
    }
}
