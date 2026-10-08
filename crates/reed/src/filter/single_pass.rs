//! The single-pass strategy (`doc_v2/03-architecture.md` §8): one chunk
//! pipeline that runs the shared record pipeline per line and edits+writes
//! the survivors through a worker-local [`BlockWriter`].

use std::sync::Arc;

use reed_core::sample::SampleMask;
use reed_core::ChunkId;
use reed_engine::error::EngineError;
use reed_engine::pipeline::ChunkJob;
use reed_filter::eval::edit_plan::fill_info_set;
use reed_filter::funnel::FunnelPartial;
use reed_filter::plan::FilterPlan;
use reed_io::vcf::line_edit::{edit_vcf_line, EditPlan, InfoValueOut};
use reed_io::bcf::{BcfEncoder, NativeEncoder};
use reed_io::write::block_writer::{BlockWriter, ChunkOutput, IndexKey};

use super::record_pipeline::{gt_subfield_index, PipelineOptions, RecordPipeline};
use super::job_err;

/// The BCF encoders a worker may use: the native B4 encoder by default,
/// noodles' B2 encoder behind `REED_BCF_ENCODER=noodles` (the two must
/// produce equivalent records; the comparison is part of the B4 gate).
pub(super) struct BcfEncoders {
    header: Arc<noodles::vcf::Header>,
    native: NativeEncoder,
    noodles: Option<BcfEncoder>,
}

impl BcfEncoders {
    fn new(header: Arc<noodles::vcf::Header>) -> Result<Self, reed_io::IoError> {
        let native = NativeEncoder::new(&header)?;
        let noodles = if std::env::var("REED_BCF_ENCODER").is_ok_and(|v| v == "noodles") {
            let mut encoder = BcfEncoder::new()?;
            encoder.write_header(&header)?;
            Some(encoder)
        } else {
            None
        };
        Ok(Self { header, native, noodles })
    }

    fn encode(&mut self, line: &[u8], out: &mut Vec<u8>) -> Result<(), reed_io::IoError> {
        match &mut self.noodles {
            Some(encoder) => encoder.encode_record(&self.header, line, out),
            None => self.native.encode(line, out),
        }
    }
}

/// The chunk partial: output bytes plus the funnel contribution.
#[derive(Default)]
pub(super) struct FilterPartial {
    /// The chunk's output blocks (compressed) or raw bytes (plain).
    pub output: ChunkOutput,
    /// The funnel contribution.
    pub funnel: FunnelPartial,
}

/// The pass1 chunk job: head → mask → tail → edit → worker-local deflate.
pub(super) struct FilterJob {
    pub(super) plan: Arc<FilterPlan>,
    pub(super) keep: Arc<[bool]>,
    pub(super) keep_mask: Option<SampleMask>,
    pub(super) n_samples_in: usize,
    pub(super) n_samples_out: usize,
    pub(super) compressed: bool,
    pub(super) contig_names: Arc<[Vec<u8>]>,
    pub(super) level: u8,
    /// BCF output (`-Ob/-Ou`): the header the dictionary keys come from.
    /// Each worker builds its own encoder from it (and the environment
    /// selects the native encoder over noodles' B2 fallback for A/B tests).
    pub(super) bcf: Option<Arc<noodles::vcf::Header>>,
}

/// Worker scratch (reused across chunks).
pub(super) struct WorkerScratch {
    /// The shared record pipeline (head → decode → mask → count → tail).
    pipeline: RecordPipeline,
    info_set: Vec<(&'static [u8], InfoValueOut)>,
    edit_buf: Vec<u8>,
    /// Encoded BCF bytes before they enter the block writer.
    bcf_buf: Vec<u8>,
    /// This worker's BCF encoder (native by default, noodles on
    /// `REED_BCF_ENCODER=noodles` for the B2/B4 byte comparison).
    bcf: Option<BcfEncoders>,
    writer: BlockWriter,
    funnel: FunnelPartial,
}

impl ChunkJob for FilterJob {
    type Partial = FilterPartial;
    type Scratch = Option<WorkerScratch>;

    fn run_chunk(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        scratch.funnel.reset(
            self.plan.head.len(),
            self.plan.tail.len(),
            self.plan.mask.len(),
        );
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            self.push_line(line, scratch)?;
        }
        let output = scratch.writer.finish_chunk().map_err(job_err)?;
        let partial = FilterPartial {
            output,
            funnel: std::mem::take(&mut scratch.funnel),
        };
        Ok(partial)
    }

    fn run_seam(&self, _id: ChunkId, line: &[u8], scratch: &mut Self::Scratch) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        scratch.funnel.reset(
            self.plan.head.len(),
            self.plan.tail.len(),
            self.plan.mask.len(),
        );
        self.push_line(line, scratch)?;
        let output = scratch.writer.finish_chunk().map_err(job_err)?;
        Ok(FilterPartial {
            output,
            funnel: std::mem::take(&mut scratch.funnel),
        })
    }
}

impl FilterJob {
    fn build_scratch(&self) -> WorkerScratch {
        let int_keys: Vec<&[u8]> = self.plan.int_fields.iter().map(Vec::as_slice).collect();
        let float_keys: Vec<&[u8]> = self.plan.float_fields.iter().map(Vec::as_slice).collect();
        let mut pipeline = RecordPipeline::new(
            self.n_samples_in,
            &int_keys,
            &float_keys,
            self.plan.policies,
            &self.contig_names,
        );
        pipeline.resize_for(&self.plan);
        WorkerScratch {
            pipeline,
            info_set: Vec::new(),
            edit_buf: Vec::new(),
            bcf_buf: Vec::new(),
            bcf: self.bcf.as_ref().map(|header| {
                BcfEncoders::new(Arc::clone(header))
                    .expect("invariant: the output header encodes")
            }),
            writer: if self.compressed {
                BlockWriter::new(self.level, true)
            } else {
                BlockWriter::new(self.level, false)
            }
            .expect("invariant: the level comes from the CLI clamp"),
            funnel: FunnelPartial {
                head_hits: vec![0; self.plan.head.len()],
                tail_hits: vec![0; self.plan.tail.len()],
                mask_hits: vec![0; self.plan.mask.len()],
                ..FunnelPartial::default()
            },
        }
    }

    /// Runs the shared record pipeline and writes the survivors: head →
    /// decode → mask → count → tail → edit → worker-local deflate.
    fn push_line(&self, line: &[u8], scratch: &mut WorkerScratch) -> Result<(), EngineError> {
        let plan = &self.plan;
        let opts = PipelineOptions {
            // Head-failed records are never scanned on this path.
            scan_head_failures: false,
            mask_all_samples: false,
            // The contig table is the index's tid space: a record on an
            // undeclared contig has no addressable tid and must fail loudly.
            require_contig: true,
            tail_n_samples: Some(self.n_samples_out as u64),
            counts: true,
        };
        let WorkerScratch { pipeline, info_set, edit_buf, bcf_buf, bcf, writer, funnel } = scratch;
        pipeline.process(line, plan, &self.keep, opts, funnel, |view| {
            let site = view.site;
            let counts = view.counts.expect("invariant: counts requested");

            // ---- edit + emit ----
            fill_info_set(counts, &plan.recompute, info_set);
            let stale: Vec<&[u8]> = plan.stale_tags.iter().map(Vec::as_slice).collect();
            let edit_plan = EditPlan {
                keep_samples: self.keep_mask.as_ref(),
                mask_gt: Some(view.masked),
                info_remove: &stale,
                info_set,
            };
            let gt_index = gt_subfield_index(view.line.format_text());
            edit_buf.clear();
            edit_vcf_line(view.line, &edit_plan, gt_index, edit_buf).map_err(job_err)?;

            view.funnel.n_output += 1;
            // Block boundaries are byte-based, not record-based (htslib): the
            // record's index entry anchors at the block and offset where its
            // bytes start, even when they span into the next block.
            let key = IndexKey {
                contig: site
                    .contig
                    .expect("invariant: parse_declared resolved the contig"),
                start1: view.pos1,
                ref_len: u32::try_from(site.reference.len()).unwrap_or(1).max(1),
            };
            match bcf {
                None => writer.push_record(edit_buf, key).map_err(job_err),
                // The edited text is encoded to BCF; the record's bytes then
                // flow through the block writer like any other payload.
                Some(encoders) => {
                    bcf_buf.clear();
                    encoders.encode(edit_buf, bcf_buf).map_err(job_err)?;
                    writer.push_record(bcf_buf, key).map_err(job_err)
                }
            }
        })?;
        Ok(())
    }
}
