//! The sketch path's pass1 emit (`docs/rounds/v2/03-architecture.md` §8, §11.4):
//! re-inflate the chunk, skip head-failed and tail-failed records with a
//! newline scan, and edit+compress survivors from the mid-stage decisions.

use std::sync::Arc;

use reed_core::sample::SampleMask;
use reed_core::{Bitset, ChunkId};
use reed_engine::error::EngineError;
use reed_engine::pipeline::ChunkJob;
use reed_filter::funnel::FunnelPartial;
use reed_filter::plan::FilterPlan;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::line_edit::{EditPlan, edit_vcf_line};
use reed_io::vcf::site::{ContigTable, SiteView};
use reed_io::write::block_writer::{BlockWriter, IndexKey};

use super::job_err;
use super::mid_stage::{ChunkDecision, SiteEmit};
use super::record_pipeline::gt_subfield_index;
use super::single_pass::FilterPartial;

pub(super) struct EmitJob {
    pub(super) plan: Arc<FilterPlan>,
    pub(super) decisions: Arc<std::collections::HashMap<u64, Arc<ChunkDecision>>>,
    pub(super) keep_mask: Option<SampleMask>,
    pub(super) n_samples_in: usize,
    pub(super) compressed: bool,
    pub(super) contig_names: Arc<[Vec<u8>]>,
    pub(super) level: u8,
}

pub(super) struct EmitScratch {
    edit_buf: Vec<u8>,
    writer: BlockWriter,
    contigs: ContigTable,
    mask_bitset: Bitset,
}

impl ChunkJob for EmitJob {
    type Partial = FilterPartial;
    type Scratch = Option<EmitScratch>;

    fn run_chunk(
        &self,
        chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let dec = self.decisions.get(&chunk.id.0).ok_or_else(|| {
            EngineError::job(format!("no emit decision for chunk {}", chunk.id.0))
        })?;
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        let mut record = 0_usize;
        let mut site = 0_usize;
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            let passed = dec.head_passed(record);
            record += 1;
            if !passed {
                continue; // head-failed: newline scan only
            }
            let s = site;
            site += 1;
            if !dec.survivor(s) {
                continue; // tail-failed: newline scan only
            }
            let Some(emit) = dec.survivors.get(s).and_then(|e| e.as_ref()) else {
                continue;
            };
            self.emit_line(line, emit, scratch)?;
        }
        let output = scratch.writer.finish_chunk().map_err(job_err)?;
        Ok(FilterPartial {
            output,
            funnel: FunnelPartial::default(),
        })
    }

    fn run_seam(
        &self,
        id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        // The seam record belongs to the previous chunk's stream; its
        // decision is that chunk's `seam` field.
        let emitted =
            id.0.checked_sub(1)
                .and_then(|prev| self.decisions.get(&prev))
                .and_then(|dec| dec.seam.clone());
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        if let Some(emit) = emitted {
            self.emit_line(line, &emit, scratch)?;
        }
        let output = scratch.writer.finish_chunk().map_err(job_err)?;
        Ok(FilterPartial {
            output,
            funnel: FunnelPartial::default(),
        })
    }
}

impl EmitJob {
    fn build_scratch(&self) -> EmitScratch {
        EmitScratch {
            edit_buf: Vec::new(),
            writer: BlockWriter::new(self.level, self.compressed)
                .expect("invariant: the level comes from the CLI clamp"),
            contigs: ContigTable::new(self.contig_names.iter().cloned()),
            mask_bitset: Bitset::new(self.n_samples_in),
        }
    }

    /// Edits one survivor line from its mid-stage emit and appends the
    /// bytes to the output block (or the raw buffer).
    fn emit_line(
        &self,
        line: &[u8],
        emit: &SiteEmit,
        scratch: &mut EmitScratch,
    ) -> Result<(), EngineError> {
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::job("survivor line has fewer than 9 columns"));
        };
        let site =
            SiteView::parse_declared(&lr, &mut scratch.contigs).map_err(EngineError::from)?;
        let pos1 = site.pos1;
        let stale: Vec<&[u8]> = self.plan.stale_tags.iter().map(Vec::as_slice).collect();
        // Refresh the reusable mask bitset from the emit's words.
        for (word, src) in scratch
            .mask_bitset
            .words_mut()
            .iter_mut()
            .zip(&emit.mask_words)
        {
            *word = *src;
        }
        scratch.mask_bitset.clear_padding();
        let edit_plan = EditPlan {
            keep_samples: self.keep_mask.as_ref(),
            mask_gt: Some(&scratch.mask_bitset),
            info_remove: &stale,
            info_set: &emit.infos,
        };
        scratch.edit_buf.clear();
        edit_vcf_line(
            &lr,
            &edit_plan,
            gt_subfield_index(lr.format_text()),
            &mut scratch.edit_buf,
        )
        .map_err(EngineError::from)?;
        // Byte-boundary cutting, records may span blocks (as in
        // FilterJob::push_line).
        let key = IndexKey {
            contig: site
                .contig
                .expect("invariant: parse_declared resolved the contig"),
            start1: pos1,
            ref_len: u32::try_from(site.reference.len()).unwrap_or(1).max(1),
        };
        scratch
            .writer
            .push_record(&scratch.edit_buf, key)
            .map_err(job_err)?;
        Ok(())
    }
}
