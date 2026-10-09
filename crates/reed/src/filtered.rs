//! A filtering wrapper job (04-cleanup §2.1): applies the compiled filter
//! plan to the text stream before an inner consumer sees it, so `reed stats`
//! can accept the filter's site/threshold flags in the *same pass* instead
//! of running the filter API into an intermediate file.
//!
//! The wrapper runs the filter's own [`RecordPipeline`] (head → mask →
//! tail) and hands the survivors to the inner job. Two hand-over paths:
//!
//! - text ([`FilteredJob::new`]): the survivors are rewritten with the
//!   filter's [`edit_vcf_line`] (masked genotypes become `./.`) and the
//!   inner job re-parses them — exactly what a physical `reed filter` +
//!   `reed stats` would produce;
//! - records ([`FilteredJob::with_records`], doc_v4 O-2): the wrapper hands
//!   over the original line plus the already-masked genotype codes, and the
//!   inner job skips the text rewrite and the re-scan. Only jobs that
//!   declare the entry point use it.
//!
//! INFO updates stay filter-only, so neither path touches INFO.

use std::sync::Arc;

use reed_core::genotype::GtCode;
use reed_core::ChunkId;
use reed_engine::error::EngineError;
use reed_engine::pipeline::ChunkJob;
use reed_filter::funnel::FunnelPartial;
use reed_filter::plan::FilterPlan;
use reed_io::chunk::ChunkBytes;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::line_edit::{edit_vcf_line, EditPlan};

use crate::filter::record_pipeline::{PipelineOptions, RecordPipeline};

/// One survivor handed to a record-consuming inner job (doc_v4 O-2): the
/// original line plus the decoded genotypes with the filter's mask applied.
pub struct ScannedRecord<'a> {
    /// The record's original bytes (INFO/FORMAT untouched).
    pub line: &'a [u8],
    /// One masked [`GtCode`] per sample column.
    pub codes: &'a [GtCode],
}

/// The record-level entry point of a wrapped job (the
/// [`FilteredJob::with_records`] counterpart of [`ChunkJob::run_chunk`]).
/// Both the chunk body and the seam arrive through it, so there is no chunk
/// parameter — the consumer must not key state on chunk identity.
pub type ScannedRecordFn<J> =
    fn(&J, &[ScannedRecord<'_>], &mut <J as ChunkJob>::Scratch) -> Result<<J as ChunkJob>::Partial, EngineError>;

/// Wraps an inner job with the filter's per-record stages; `plan = None`
/// passes every record through untouched (no buffering, no pipeline).
pub struct FilteredJob<J: ChunkJob> {
    plan: Option<Arc<FilterPlan>>,
    /// Samples the filter keeps (all when no static subset was requested).
    keep: Arc<[bool]>,
    /// The retained sample count, the tail rules' missing-rate denominator.
    n_out: u64,
    /// The record path's entry point; `None` uses rewritten text instead.
    records: Option<ScannedRecordFn<J>>,
    inner: J,
}

/// Worker scratch: the record pipeline, the rewritten-record buffer, the
/// record path's line/code buffers and the inner job's scratch. The pipeline
/// is built lazily (it needs the plan).
pub struct FilteredScratch<S> {
    pipeline: Option<RecordPipeline>,
    /// Survivors, newline-terminated (the text path's inner-job body).
    buffer: Vec<u8>,
    /// The record path: masked codes and `(line offset, line len, code
    /// offset)` spans into the two buffers.
    codes: Vec<GtCode>,
    spans: Vec<(u32, u32, u32)>,
    funnel: FunnelPartial,
    inner: S,
}

impl<S: Default> Default for FilteredScratch<S> {
    fn default() -> Self {
        Self {
            pipeline: None,
            buffer: Vec::new(),
            codes: Vec::new(),
            spans: Vec::new(),
            funnel: FunnelPartial::default(),
            inner: S::default(),
        }
    }
}

impl<J: ChunkJob> FilteredJob<J> {
    /// Wraps `inner` with `plan` (the same compiler `reed filter` uses) over
    /// the sample keep mask; survivors reach the inner job as rewritten text.
    #[must_use]
    pub fn new(plan: Option<Arc<FilterPlan>>, keep: Arc<[bool]>, inner: J) -> Self {
        Self::build(plan, keep, None, inner)
    }

    /// Wraps `inner` with the record path: survivors reach `records` as
    /// [`ScannedRecord`]s, so the inner job neither re-parses the rewritten
    /// record nor decodes GT again (doc_v4 O-2).
    #[must_use]
    pub fn with_records(
        plan: Option<Arc<FilterPlan>>,
        keep: Arc<[bool]>,
        records: ScannedRecordFn<J>,
        inner: J,
    ) -> Self {
        Self::build(plan, keep, Some(records), inner)
    }

    fn build(
        plan: Option<Arc<FilterPlan>>,
        keep: Arc<[bool]>,
        records: Option<ScannedRecordFn<J>>,
        inner: J,
    ) -> Self {
        let n_out = u64::try_from(keep.iter().filter(|k| **k).count()).unwrap_or(u64::MAX);
        Self { plan, keep, n_out, records, inner }
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
            codes: Vec::new(),
            spans: Vec::new(),
            funnel: FunnelPartial::default(),
            inner: J::Scratch::default(),
        }
    }

    /// Stages one record for the text path: survivors are appended to
    /// `buffer` with their masked genotypes rewritten.
    fn filter_line(
        &self,
        line: &[u8],
        scratch: &mut FilteredScratch<J::Scratch>,
    ) -> Result<(), EngineError> {
        let plan = self.plan.as_ref().expect("invariant: guarded by the caller");
        let FilteredScratch { pipeline, buffer, funnel, .. } = scratch;
        let pipeline = pipeline.as_mut().expect("invariant: start_chunk built it");
        let opts = self.pipeline_options();
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

    /// Stages one record for the record path: the original line plus the
    /// masked codes land in the buffers, with a span for the view below.
    fn scan_line(
        &self,
        line: &[u8],
        scratch: &mut FilteredScratch<J::Scratch>,
    ) -> Result<(), EngineError> {
        let plan = self.plan.as_ref().expect("invariant: guarded by the caller");
        let FilteredScratch { pipeline, buffer, codes, spans, funnel, .. } = scratch;
        let pipeline = pipeline.as_mut().expect("invariant: start_chunk built it");
        let opts = self.pipeline_options();
        let line_start = buffer.len();
        let code_start = codes.len();
        let processed = pipeline.process(line, plan, &self.keep, opts, funnel, |view| {
            for (i, &code) in view.scan.codes.iter().enumerate() {
                codes.push(if view.masked.get(i) { code.masked() } else { code });
            }
            buffer.extend_from_slice(line);
            Ok(())
        })?;
        if processed.emitted {
            spans.push((
                u32::try_from(line_start).unwrap_or(u32::MAX),
                u32::try_from(buffer.len() - line_start).unwrap_or(u32::MAX),
                u32::try_from(code_start).unwrap_or(u32::MAX),
            ));
        } else {
            // The tail dropped the record before the body ran; keep the
            // buffers exactly as they were.
            buffer.truncate(line_start);
            codes.truncate(code_start);
        }
        Ok(())
    }

    fn pipeline_options(&self) -> PipelineOptions {
        PipelineOptions {
            scan_head_failures: false,
            mask_all_samples: true,
            require_contig: false,
            // The tail's missing-rate denominator is the retained sample set.
            tail_n_samples: Some(self.n_out),
            counts: true,
        }
    }

    /// Builds the record views over the scratch buffers and runs the inner
    /// entry point.
    fn run_records(
        &self,
        scratch: &mut FilteredScratch<J::Scratch>,
    ) -> Result<J::Partial, EngineError> {
        let f = self.records.expect("invariant: the record path is guarded");
        let FilteredScratch { buffer, codes, spans, inner, .. } = scratch;
        let n_samples = self.keep.len();
        let mut views: Vec<ScannedRecord<'_>> = Vec::with_capacity(spans.len());
        for &(line_off, line_len, code_off) in spans.iter() {
            let (line_off, line_len, code_off) =
                (line_off as usize, line_len as usize, code_off as usize);
            views.push(ScannedRecord {
                line: &buffer[line_off..line_off + line_len],
                codes: &codes[code_off..code_off + n_samples],
            });
        }
        f(&self.inner, &views, inner)
    }

    fn start_chunk(&self, scratch: &mut FilteredScratch<J::Scratch>) {
        let Some(plan) = &self.plan else { return };
        if scratch.pipeline.is_none() {
            scratch.pipeline = Some(self.build_scratch().pipeline.expect("invariant: built"));
            scratch.buffer = Vec::with_capacity(1 << 16);
            scratch.codes = Vec::new();
            scratch.spans = Vec::new();
        }
        scratch.funnel.reset(plan.head.len(), plan.tail.len(), plan.mask.len());
        scratch.buffer.clear();
        scratch.codes.clear();
        scratch.spans.clear();
    }

    /// Whether the record path is active (the inner job must expect
    /// [`ScannedRecord`]s instead of text).
    #[must_use]
    pub const fn uses_records(&self) -> bool {
        self.records.is_some()
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
            if self.records.is_some() {
                self.scan_line(line, scratch)?;
            } else {
                self.filter_line(line, scratch)?;
            }
        }
        if self.records.is_some() {
            return self.run_records(scratch);
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
        if self.records.is_some() {
            self.scan_line(line, scratch)?;
            if scratch.spans.is_empty() {
                // The filter dropped the seam record.
                return Ok(J::Partial::default());
            }
            let _ = id;
            return self.run_records(scratch);
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
