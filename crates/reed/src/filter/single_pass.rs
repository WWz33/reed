//! The single-pass strategy (`docs/rounds/v2/03-architecture.md` §8): one chunk
//! pipeline that runs the shared record pipeline per line and edits+writes
//! the survivors through a worker-local [`BlockWriter`].

use std::sync::Arc;

use reed_core::ChunkId;
use reed_core::math::hwe::HweCache;
use reed_core::sample::SampleMask;
use reed_engine::error::EngineError;
use reed_engine::pipeline::ChunkJob;
use reed_filter::eval::edit_plan::fill_info_set;
use reed_filter::funnel::FunnelPartial;
use reed_filter::plan::FilterPlan;
use reed_io::bcf::{BcfDict, record_ref, render_record};
use reed_io::bcf::{BcfEncoder, NativeEncoder};
use reed_io::vcf::line::LineRef;
use reed_io::vcf::line_edit::{EditPlan, InfoValueOut, edit_vcf_line};
use reed_io::write::block_writer::{BlockWriter, ChunkOutput, IndexKey};

use super::job_err;
use super::record_pipeline::{
    NativeScan, NativeScratch, PipelineOptions, RecordPipeline, RecordView, gt_subfield_index,
};

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
        Ok(Self {
            header,
            native,
            noodles,
        })
    }

    /// Encodes the edited text through noodles (the oracle A/B path).
    fn encode_noodles(&mut self, line: &[u8], out: &mut Vec<u8>) -> Result<(), reed_io::IoError> {
        self.noodles
            .as_mut()
            .expect("invariant: the noodles branch checks the encoder exists")
            .encode_record(&self.header, line, out)
    }
}

/// The record the emitter encodes/edits from (doc_v6 03 §8).
enum EmitSource<'a, 'b> {
    /// Text input: the record's own line.
    Text(&'a LineRef<'b>),
    /// BCF input: the framed record bytes and the input dictionary. Only the
    /// native `-Ob` branch takes it; text output (and the noodles oracle)
    /// renders the record first.
    Bcf { bytes: &'a [u8], dict: &'a BcfDict },
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
    /// BCF *input* (doc_v4 O-4): when set, the job takes the native record
    /// path — the site line plus decoded GT/FORMAT columns — instead of
    /// rendered text.
    pub(super) bcf_input: Option<Arc<BcfDict>>,
    /// The GT decode table for the native path, built once per run and
    /// shared by every worker (A6-10). `Some` exactly when `bcf_input` is.
    pub(super) gt_table: Option<Arc<reed_io::bcf::GtTable>>,
    /// The plan's Integer/Float FORMAT keys (the mask's columns).
    pub(super) int_keys: Arc<[Vec<u8>]>,
    /// See [`FilterJob::int_keys`].
    pub(super) float_keys: Arc<[Vec<u8>]>,
    /// `--indep-pairwise`: the kept sites, per contig, from pass0
    /// (docs/rounds/v5/05-algorithms.md §3.4.1). Pass1 drops a surviving
    /// record the pruner removed.
    pub(super) site_keep: Option<Arc<super::pass0::SiteKeepSet>>,
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
    /// Native BCF input buffers (doc_v4 O-4).
    native: Option<NativeScratch>,
    /// The HWE memo behind `--fill-tags HWE` (A6-16).
    hwe: HweCache,
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
        // One build per chunk, not per surviving record (A6-15).
        let stale = self.stale_refs();
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            self.push_line(line, scratch, &stale)?;
        }
        let output = scratch.writer.finish_chunk().map_err(job_err)?;
        let partial = FilterPartial {
            output,
            funnel: std::mem::take(&mut scratch.funnel),
        };
        Ok(partial)
    }

    fn run_seam(
        &self,
        _id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        scratch.funnel.reset(
            self.plan.head.len(),
            self.plan.tail.len(),
            self.plan.mask.len(),
        );
        let stale = self.stale_refs();
        self.push_line(line, scratch, &stale)?;
        let output = scratch.writer.finish_chunk().map_err(job_err)?;
        Ok(FilterPartial {
            output,
            funnel: std::mem::take(&mut scratch.funnel),
        })
    }
}

/// The scratch pieces the survivor emitter touches. They are borrowed out
/// of the worker scratch as one bundle (the record pipeline holds the rest),
/// so the emit step is one call in both the text and the native path (A6-20).
struct EmitScratchMut<'a> {
    info_set: &'a mut Vec<(&'static [u8], InfoValueOut)>,
    edit_buf: &'a mut Vec<u8>,
    bcf_buf: &'a mut Vec<u8>,
    bcf: &'a mut Option<BcfEncoders>,
    writer: &'a mut BlockWriter,
    /// The HWE memo behind `--fill-tags HWE` (A6-16).
    hwe: &'a mut HweCache,
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
                BcfEncoders::new(Arc::clone(header)).expect("invariant: the output header encodes")
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
            native: match (&self.bcf_input, &self.gt_table) {
                (Some(_), Some(table)) => Some(NativeScratch::new(
                    std::sync::Arc::clone(table),
                    self.n_samples_in,
                    &self.int_keys,
                    &self.float_keys,
                )),
                _ => None,
            },
            hwe: HweCache::new(),
        }
    }

    /// The native BCF record path (doc_v4 O-4): the framed records arrive
    /// decoded from the BCF side, so the pipeline decides on the site line
    /// plus typed GT/FORMAT values and the full text is rendered only for
    /// survivors.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] for malformed records, policy rejections and
    /// write failures.
    pub(super) fn push_bcf_records(
        &self,
        records: &[&[u8]],
        scratch: &mut WorkerScratch,
        stale: &[&[u8]],
    ) -> Result<(), EngineError> {
        let dict = self
            .bcf_input
            .as_ref()
            .expect("invariant: the native path needs the input dictionary");
        for bytes in records {
            // Take the native buffers out so the pipeline's scratch and the
            // decode buffers are borrowed disjointly.
            let mut native = scratch
                .native
                .take()
                .expect("invariant: the native scratch exists with the dictionary");
            let result = self.push_bcf_record(dict, bytes, scratch, &mut native, stale);
            scratch.native = Some(native);
            result?;
        }
        Ok(())
    }

    /// The [`crate::bcf_input::BcfRecordFn`] entry point (doc_v4 O-4): one
    /// chunk's framed records in, one chunk partial out.
    ///
    /// # Errors
    ///
    /// As [`Self::push_bcf_records`].
    pub(super) fn run_bcf_records(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        records: &[&[u8]],
        scratch: &mut Option<WorkerScratch>,
    ) -> Result<FilterPartial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        scratch.funnel.reset(
            self.plan.head.len(),
            self.plan.tail.len(),
            self.plan.mask.len(),
        );
        let stale = self.stale_refs();
        self.push_bcf_records(records, scratch, &stale)?;
        let output = scratch.writer.finish_chunk().map_err(job_err)?;
        Ok(FilterPartial {
            output,
            funnel: std::mem::take(&mut scratch.funnel),
        })
    }

    /// The plan's stale INFO tags as borrowed slices; built once per chunk
    /// by the callers (A6-15 — it used to be one `Vec` per surviving
    /// record).
    fn stale_refs(&self) -> Vec<&[u8]> {
        self.plan.stale_tags.iter().map(Vec::as_slice).collect()
    }

    /// The one emit step shared by the text and native paths (A6-20): the
    /// pass0 pruning rule, the recomputed INFO fill, the edit and the write.
    /// `source` is the record itself — its text line, or the framed BCF
    /// bytes for BCF input. The native `-Ob` branch walks the source with
    /// the edit applied (A6-9: no edited text is materialised and re-parsed);
    /// text output and the noodles oracle still edit the line first.
    fn emit_survivor(
        &self,
        emit: &mut EmitScratchMut<'_>,
        stale: &[&[u8]],
        view: &mut RecordView<'_, '_>,
        source: EmitSource<'_, '_>,
    ) -> Result<(), EngineError> {
        // The LD-pruning rule (§3.4.1): pass0 removed this site.
        if let Some(keep) = self.site_keep.as_deref()
            && !keep.contains(view.site.contig_name, view.pos1.saturating_sub(1))
        {
            view.funnel.n_ld_pruned += 1;
            return Ok(());
        }
        let counts = view.counts.expect("invariant: counts requested");
        fill_info_set(counts, &self.plan.recompute, emit.hwe, emit.info_set);
        // A GT-less record keeps the input's AC/AN: the automatic update has
        // no value to put back, so only user-requested removals apply
        // (doc_v12 D12-2).
        let auto: Vec<&[u8]> = self.plan.auto_stale_tags.iter().map(Vec::as_slice).collect();
        let stale_record = reed_filter::eval::edit_plan::stale_for_record(
            stale,
            &auto,
            counts.has_gt,
        );
        let masked = view.masked;
        let site = view.site;
        let pos1 = view.pos1;
        let edit_plan = EditPlan {
            keep_samples: self.keep_mask.as_ref(),
            mask_gt: Some(masked),
            info_remove: &stale_record,
            info_set: emit.info_set,
        };
        view.funnel.n_output += 1;
        // Block boundaries are byte-based, not record-based (htslib): the
        // record's index entry anchors at the block and offset where its
        // bytes start, even when they span into the next block.
        let key = IndexKey {
            contig: site
                .contig
                .expect("invariant: parse_declared resolved the contig"),
            start1: pos1,
            ref_len: u32::try_from(site.reference.len()).unwrap_or(1).max(1),
        };
        // Native BCF output: encode from the source record with the edit
        // applied in the same walk (no edited-text round-trip).
        if let Some(encoders) = emit.bcf.as_mut().filter(|e| e.noodles.is_none()) {
            emit.bcf_buf.clear();
            match source {
                EmitSource::Text(line) => encoders
                    .native
                    .encode_text(line, site, self.n_samples_in, &edit_plan, emit.bcf_buf)
                    .map_err(job_err)?,
                EmitSource::Bcf { bytes, dict } => {
                    let rec = record_ref(bytes).map_err(EngineError::from)?;
                    encoders
                        .native
                        .encode_bcf(&rec, dict, &edit_plan, emit.bcf_buf)
                        .map_err(job_err)?;
                }
            }
            return emit.writer.push_record(emit.bcf_buf, key).map_err(job_err);
        }
        let edit_line = match source {
            EmitSource::Text(line) => line,
            EmitSource::Bcf { .. } => {
                unreachable!("invariant: the BCF caller renders for text output")
            }
        };
        let gt_index = gt_subfield_index(edit_line.format_text());
        emit.edit_buf.clear();
        edit_vcf_line(edit_line, &edit_plan, gt_index, emit.edit_buf).map_err(job_err)?;
        match emit.bcf {
            None => emit.writer.push_record(emit.edit_buf, key).map_err(job_err),
            Some(encoders) => {
                emit.bcf_buf.clear();
                encoders
                    .encode_noodles(emit.edit_buf, emit.bcf_buf)
                    .map_err(job_err)?;
                emit.writer.push_record(emit.bcf_buf, key).map_err(job_err)
            }
        }
    }

    fn push_bcf_record(
        &self,
        dict: &BcfDict,
        bytes: &[u8],
        scratch: &mut WorkerScratch,
        native: &mut NativeScratch,
        stale: &[&[u8]],
    ) -> Result<(), EngineError> {
        native.decode(
            dict,
            bytes,
            &self.int_keys,
            &self.float_keys,
            self.n_samples_in,
        )?;
        let plan = &self.plan;
        let opts = PipelineOptions {
            mask_all_samples: false,
            require_contig: true,
            tail_n_samples: Some(self.n_samples_out as u64),
            counts: true,
        };
        let NativeScratch {
            site_line,
            full,
            codes,
            int_columns,
            float_columns,
            counters,
            has_gt,
            ..
        } = native;
        // The survivor render buffer is moved out for the closure and put
        // back afterwards so its capacity survives across records.
        let full_slot: &mut Vec<u8> = full;
        let mut full = std::mem::take(full_slot);
        let WorkerScratch {
            pipeline,
            info_set,
            edit_buf,
            bcf_buf,
            bcf,
            writer,
            funnel,
            ..
        } = scratch;
        let WorkerScratch { hwe, .. } = scratch;
        let mut emit = EmitScratchMut {
            info_set,
            edit_buf,
            bcf_buf,
            bcf,
            writer,
            hwe,
        };
        pipeline.process_scanned(
            site_line,
            plan,
            &self.keep,
            opts,
            funnel,
            NativeScan {
                codes: &codes[..],
                int_columns: &int_columns[..],
                float_columns: &float_columns[..],
                counters,
                has_gt: *has_gt,
            },
            |view| {
                // The native `-Ob` branch encodes the framed record directly
                // (doc_v6 03 §8); text output and the noodles oracle render
                // the full record and edit the text.
                if emit.bcf.as_ref().is_some_and(|e| e.noodles.is_none()) {
                    return self.emit_survivor(
                        &mut emit,
                        stale,
                        view,
                        EmitSource::Bcf { bytes, dict },
                    );
                }
                let rec = record_ref(bytes).map_err(EngineError::from)?;
                full.clear();
                render_record(&rec, dict, &mut full).map_err(EngineError::from)?;
                let full_lr = LineRef::parse(&full)
                    .ok_or_else(|| EngineError::job("rendered record has fewer than 9 columns"))?;
                self.emit_survivor(&mut emit, stale, view, EmitSource::Text(&full_lr))
            },
        )?;
        *full_slot = full;
        Ok(())
    }

    /// Runs the shared record pipeline and writes the survivors: head →
    /// decode → mask → count → tail → edit → worker-local deflate.
    fn push_line(
        &self,
        line: &[u8],
        scratch: &mut WorkerScratch,
        stale: &[&[u8]],
    ) -> Result<(), EngineError> {
        let plan = &self.plan;
        let opts = PipelineOptions {
            // Head-failed records are never scanned on this path.
            mask_all_samples: false,
            // The contig table is the index's tid space: a record on an
            // undeclared contig has no addressable tid and must fail loudly.
            require_contig: true,
            tail_n_samples: Some(self.n_samples_out as u64),
            counts: true,
        };
        let WorkerScratch {
            pipeline,
            info_set,
            edit_buf,
            bcf_buf,
            bcf,
            writer,
            funnel,
            ..
        } = scratch;
        let WorkerScratch { hwe, .. } = scratch;
        let mut emit = EmitScratchMut {
            info_set,
            edit_buf,
            bcf_buf,
            bcf,
            writer,
            hwe,
        };
        pipeline.process(line, plan, &self.keep, opts, funnel, |view| {
            let edit_line = view.line;
            self.emit_survivor(&mut emit, stale, view, EmitSource::Text(edit_line))
        })?;
        Ok(())
    }
}
