//! BCF input through the shared text pipeline (doc_v3 02-bcf.md B1/B3a).
//!
//! A BCF chunk has no newline framing, so the engine hands the job the whole
//! inflated chunk ([`RecordFraming::Raw`]) and this wrapper does the rest:
//! it takes the carried bytes of the record that started in the previous
//! chunk from the [`FrameChain`], walks the length prefixes, renders every
//! complete record to VCF text, and hands that text to the command's own
//! job unchanged. Every consumer therefore keeps its single implementation:
//! the only format-specific code is the renderer in `reed_io::bcf`.

// Casts in this file are bounded by construction: Carried record prefixes are bounded by the 1 GiB record ceiling.
#![allow(clippy::cast_possible_truncation)]

use std::sync::Arc;

use reed_engine::abort::Abort;
use reed_engine::error::EngineError;
use reed_engine::frame_chain::{Carried, FrameChain};
use reed_engine::pipeline::ChunkJob;
use reed_io::chunk::ChunkBytes;
use reed_io::bcf::{declared_len, record_len, record_ref, render_record, BcfDict, VirtualMap};
use reed_io::chunk::framing::RecordFraming;

use crate::input::OpenedInput;

/// The largest record the chain will carry; a runaway length prefix must not
/// allocate.
const MAX_RECORD_BYTES: usize = 1 << 30;

/// The BCF side of the wrapper: the frame chain plus the tag dictionaries.
struct BcfState {
    chain: Arc<FrameChain>,
    dict: Arc<BcfDict>,
}

/// The native BCF record hook (doc_v3 02-bcf.md §2.3–§2.5): a job entry
/// point that takes the chunk's framed record bytes instead of rendered
/// text, so a consumer that only needs a few fields never pays for the text
/// round-trip.
pub type BcfRecordFn<J> = fn(
    &J,
    &ChunkBytes,
    &[&[u8]],
    &mut <J as ChunkJob>::Scratch,
) -> Result<<J as ChunkJob>::Partial, EngineError>;

/// Wraps a text job so it also runs on BCF input.
pub struct InputJob<J: ChunkJob> {
    inner: J,
    bcf: Option<BcfState>,
    /// When set and the input is BCF, framed records go here instead of
    /// being rendered to text.
    native: Option<BcfRecordFn<J>>,
}

/// Scratch: the wrapped job's scratch plus the render buffer.
pub struct InputScratch<S> {
    inner: S,
    rendered: Vec<u8>,
}

impl<S: Default> Default for InputScratch<S> {
    fn default() -> Self {
        Self { inner: S::default(), rendered: Vec::new() }
    }
}

impl<J: ChunkJob> InputJob<J> {
    /// The chain's cumulative wait (engine metric, doc_v4 03 A4-5).
    fn chain_wait(&self) -> std::time::Duration {
        self.bcf.as_ref().map_or(std::time::Duration::ZERO, |bcf| bcf.chain.wait_time())
    }

    /// Wraps `inner` for the input `opened` opened; BCF inputs get a chain
    /// cell 0 (chunk 0 never waits) bound to the run's abort signal.
    #[must_use]
    pub fn new(inner: J, opened: &OpenedInput, abort: &Abort) -> Self {
        let bcf = opened.bcf.as_ref().map(|dict| {
            let chain = Arc::new(FrameChain::new(abort.signal.clone()));
            chain.publish(0, Carried::default());
            BcfState { chain, dict: Arc::clone(dict) }
        });
        Self { inner, bcf, native: None }
    }

    /// Wraps `inner` with a native BCF record entry point: on BCF input the
    /// job receives the chunk's framed record bytes ([`BcfRecordFn`])
    /// instead of rendered text; text input is unaffected.
    #[must_use]
    pub fn with_bcf_records(
        inner: J,
        opened: &OpenedInput,
        abort: &Abort,
        records: BcfRecordFn<J>,
    ) -> Self {
        let mut job = Self::new(inner, opened, abort);
        job.native = Some(records);
        job
    }

    /// Ends the run: releases chain waiters and rejects a file that stopped
    /// inside a record.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] when the last cell still carries bytes.
    pub fn finish(&self, n_chunks: u64, expect_complete: bool) -> Result<(), EngineError> {
        let Some(bcf) = &self.bcf else { return Ok(()) };
        bcf.chain.finish();
        if expect_complete && bcf.chain.peek(n_chunks).is_some_and(|tail| !tail.is_empty()) {
            return Err(EngineError::Job(
                "BCF input ends inside a record (truncated file)".into(),
            ));
        }
        Ok(())
    }

    /// Frames one chunk (02-bcf.md §2.2): walks the length prefixes, updates
    /// `carried` in place, and publishes the successor's carried bytes
    /// *before* any decoding, so a successor never waits for a predecessor's
    /// render or decode.
    ///
    /// Returns `Some((body_start, body_end))` when the previous chunk's
    /// record completes here (`carried.head` then holds its bytes) plus the
    /// complete in-chunk records at `text[body_start..body_end]`; `None`
    /// means this chunk holds no complete record and the carried bytes grew.
    fn frame_chunk(
        &self,
        chunk: &ChunkBytes,
        text: &[u8],
        carried: &mut Carried,
        chain: &FrameChain,
    ) -> Result<Option<(usize, usize)>, EngineError> {
        let chunk_id = chunk.id.get();
        let mut guard = chain.publish_guard(chunk_id + 1);
        let skip = if chunk_id == 0 {
            usize::try_from(chunk.first_in_block_offset).unwrap_or(0)
        } else {
            0
        };
        let map = VirtualMap { chunk_start: chunk.compressed_start, blocks: &chunk.blocks };
        if skip >= text.len() {
            // Defensive: the header always ends inside chunk 0.
            guard.publish(Carried::default());
            return Ok(None);
        }
        let crossing_len = carried.prefix_len as usize;
        let body_start;
        let body_end;
        if carried.is_empty() {
            body_start = skip;
            body_end = walk_records(text, skip);
        } else {
            let Some(total) = declared_len_across(&carried.head, text) else {
                // Even the length prefix is still incomplete.
                carried.head.extend_from_slice(text);
                carried.prefix_len += text.len() as u64;
                guard.publish(carried.clone());
                return Ok(None);
            };
            if total > MAX_RECORD_BYTES {
                return Err(EngineError::job(format!(
                    "BCF record of {total} bytes exceeds the {MAX_RECORD_BYTES}-byte ceiling"
                )));
            }
            let need = total.saturating_sub(crossing_len);
            if need > text.len() {
                // The record spans this whole chunk.
                carried.head.extend_from_slice(text);
                carried.prefix_len += text.len() as u64;
                guard.publish(carried.clone());
                return Ok(None);
            }
            body_start = need;
            body_end = walk_records(text, need);
            carried.head.extend_from_slice(&text[..body_start]);
        }
        // Publish the crossing record of the *next* chunk before decoding.
        guard.publish(Carried {
            head: text[body_end..].to_vec(),
            prefix_len: (text.len() - body_end) as u64,
            voff: map.voff(body_end),
        });
        Ok(Some((body_start, body_end)))
    }
}

/// Walks the length prefixes from `from` and returns the offset after the
/// last complete record (the record starting there crosses the boundary).
fn walk_records(text: &[u8], from: usize) -> usize {
    let mut p = from;
    while let Some(len) = record_len(&text[p..]) {
        p += len;
    }
    p
}

impl<J: ChunkJob + Sync> ChunkJob for InputJob<J> {
    type Partial = J::Partial;
    type Scratch = InputScratch<J::Scratch>;

    fn wait_time(&self) -> std::time::Duration {
        self.chain_wait() + self.inner.wait_time()
    }

    fn framing(&self) -> RecordFraming {
        if self.bcf.is_some() {
            RecordFraming::Raw
        } else {
            RecordFraming::Newline
        }
    }

    fn run_chunk(
        &self,
        chunk: &ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let Some(bcf) = &self.bcf else {
            return self.inner.run_chunk(chunk, body, &mut scratch.inner);
        };
        let InputScratch { inner, rendered } = scratch;
        let mut carried = match bcf.chain.wait(chunk.id.get()) {
            Ok(carried) => carried,
            Err(err) => {
                bcf.chain.poison();
                return Err(err);
            }
        };
        let framed = match self.frame_chunk(chunk, body, &mut carried, &bcf.chain) {
            Ok(framed) => framed,
            Err(err) => {
                bcf.chain.poison();
                return Err(err);
            }
        };
        let Some((body_start, body_end)) = framed else {
            // No complete record in this chunk; the inner job gets an empty
            // batch (exactly what the text path would hand it).
            return match self.native {
                Some(records) => (records)(&self.inner, chunk, &[], inner),
                None => self.inner.run_chunk(chunk, &[], inner),
            };
        };
        if let Some(records_fn) = self.native {
            let mut records: Vec<&[u8]> = Vec::new();
            if !carried.head.is_empty() {
                records.push(&carried.head);
            }
            let mut q = body_start;
            while q < body_end {
                let len = record_len(&body[q..]).expect("invariant: walked records are complete");
                records.push(&body[q..q + len]);
                q += len;
            }
            return (records_fn)(&self.inner, chunk, &records, inner);
        }
        rendered.clear();
        // The record that started earlier completes at `body_start`.
        if !carried.head.is_empty() {
            let rec = record_ref(&carried.head).map_err(EngineError::from)?;
            render_record(&rec, &bcf.dict, rendered).map_err(EngineError::from)?;
        }
        let mut q = body_start;
        while q < body_end {
            let len = record_len(&body[q..]).expect("invariant: walked records are complete");
            let rec = record_ref(&body[q..q + len]).map_err(EngineError::from)?;
            render_record(&rec, &bcf.dict, rendered).map_err(EngineError::from)?;
            q += len;
        }
        self.inner.run_chunk(chunk, rendered, inner)
    }

    fn run_seam(
        &self,
        id: reed_core::ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        match &self.bcf {
            // Text framing seams through the engine as usual.
            None => self.inner.run_seam(id, line, &mut scratch.inner),
            // Raw framing never seams: this wrapper stitches records itself.
            Some(_) => Err(EngineError::Job(
                "internal: the BCF wrapper never receives seam lines".into(),
            )),
        }
    }
}

/// The record's declared total length, taking the missing prefix bytes from
/// `text` (a chunk may end inside the 8-byte prefix).
fn declared_len_across(carried: &[u8], text: &[u8]) -> Option<usize> {
    if carried.len() >= 8 {
        return declared_len(carried);
    }
    let need = 8 - carried.len();
    if text.len() < need {
        return None;
    }
    let mut head = [0_u8; 8];
    head[..carried.len()].copy_from_slice(carried);
    head[carried.len()..].copy_from_slice(&text[..need]);
    declared_len(&head)
}
