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
use reed_io::bcf::{
    BcfDict, FrameStep, VirtualMap, declared_len, frame_chunk, record_len, record_ref,
    render_record,
};
use reed_io::chunk::ChunkBytes;
use reed_io::chunk::framing::RecordFraming;

use crate::input::OpenedInput;

/// The largest record the chain will carry; a runaway length prefix must not
/// allocate. Shared with the engine's seam stitcher (A6-2).
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
    /// The run's effective record ceiling: the smaller of the engine's hard
    /// limit and `--memory` (doc_v26 V26-07). Carried bytes count against it,
    /// so a small budget cannot accumulate a prefix the pipeline itself could
    /// never admit.
    record_ceiling: usize,
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
        Self {
            inner: S::default(),
            rendered: Vec::new(),
        }
    }
}

impl<J: ChunkJob> InputJob<J> {
    /// The chain's cumulative wait (engine metric, doc_v4 03 A4-5).
    fn chain_wait(&self) -> std::time::Duration {
        self.bcf
            .as_ref()
            .map_or(std::time::Duration::ZERO, |bcf| bcf.chain.wait_time())
    }

    /// Wraps `inner` for the input `opened` opened; BCF inputs get a chain
    /// cell 0 (chunk 0 never waits) bound to the run's abort signal.
    ///
    /// `budget` sets the carried-record ceiling together with the engine's
    /// hard limit (`reed_engine::pipeline::record_ceiling`).
    #[must_use]
    pub fn new(
        inner: J,
        opened: &OpenedInput,
        abort: &Abort,
        budget: &reed_engine::budget::ThreadBudget,
    ) -> Self {
        let record_ceiling = reed_engine::pipeline::record_ceiling(budget);
        let bcf = opened.bcf.as_ref().map(|dict| {
            let chain = Arc::new(FrameChain::new(abort.signal.clone()));
            chain.publish(0, Carried::default());
            BcfState {
                chain,
                dict: Arc::clone(dict),
            }
        });
        Self {
            inner,
            record_ceiling,
            bcf,
            native: None,
        }
    }

    /// Wraps `inner` with a native BCF record entry point: on BCF input the
    /// job receives the chunk's framed record bytes ([`BcfRecordFn`])
    /// instead of rendered text; text input is unaffected.
    #[must_use]
    pub fn with_bcf_records(
        inner: J,
        opened: &OpenedInput,
        abort: &Abort,
        budget: &reed_engine::budget::ThreadBudget,
        records: BcfRecordFn<J>,
    ) -> Self {
        let mut job = Self::new(inner, opened, abort, budget);
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
        if expect_complete
            && bcf
                .chain
                .peek(n_chunks)
                .is_some_and(|tail| !tail.is_empty())
        {
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
        let map = VirtualMap {
            chunk_start: chunk.compressed_start,
            blocks: &chunk.blocks,
        };
        if skip >= text.len() {
            // Defensive: the header always ends inside chunk 0.
            guard.publish(Carried::default());
            return Ok(None);
        }
        let crossing_len = carried.prefix_len as usize;
        // doc_v28 V28-25: the seam arithmetic itself is shared with the index
        // job; this job keeps the carried bytes, the ceiling and the decode.
        let (body_start, body_end) = match frame_chunk(&carried.head, crossing_len, text, skip) {
            FrameStep::PrefixIncomplete => {
                // Even the length prefix is still incomplete.
                self.check_carried(crossing_len, text.len(), None)?;
                carried.head.extend_from_slice(text);
                carried.prefix_len += text.len() as u64;
                guard.publish(carried.clone());
                return Ok(None);
            }
            FrameStep::Spans { total } => {
                if total > self.record_ceiling {
                    return Err(EngineError::job(format!(
                        "BCF record of {total} bytes exceeds the {}-byte ceiling (the smaller of \
                         the hard limit and the --memory budget)",
                        self.record_ceiling
                    )));
                }
                // The record spans this whole chunk.
                self.check_carried(crossing_len, text.len(), Some(total))?;
                carried.head.extend_from_slice(text);
                carried.prefix_len += text.len() as u64;
                guard.publish(carried.clone());
                return Ok(None);
            }
            FrameStep::Complete {
                body_start,
                body_end,
            } => (body_start, body_end),
        };
        if crossing_len == 0 {
            // doc_v26 review item 2: the first-fragment branch publishes the
            // crossing record below without the ceiling check the has-carry
            // branch performs. `text[body_end..]` can be the whole chunk (a
            // record that starts in this chunk and only ends after it), so a
            // small --memory could admit a carried prefix far above its
            // budget. The check must cover this branch too — and, when the
            // length prefix is already complete in this chunk, the declared
            // total is the tighter of the two bounds.
            let crossing = text.len() - body_end;
            self.check_carried(0, crossing, declared_len(&text[body_end..]))?;
        } else {
            carried.head.extend_from_slice(&text[..body_start]);
        }
        // Publish the crossing record of the *next* chunk before decoding.
        // The ceiling was checked above on whichever branch produced
        // `body_end` (doc_v26 review item 2: the no-carry branch had no
        // check and could publish an over-budget prefix).
        guard.publish(Carried {
            head: text[body_end..].to_vec(),
            prefix_len: (text.len() - body_end) as u64,
            voff: map.voff(body_end),
        });
        Ok(Some((body_start, body_end)))
    }

    /// Rejects a carried prefix that would pass the run's ceiling, before the
    /// caller appends `added` bytes to it (doc_v26 V26-07).
    ///
    /// `declared` is the record's total length when the prefix already holds
    /// a complete length field; otherwise the accumulated length is the only
    /// bound available.
    fn check_carried(
        &self,
        have: usize,
        added: usize,
        declared: Option<usize>,
    ) -> Result<(), EngineError> {
        let next = have.saturating_add(added);
        let limit = self.record_ceiling;
        let over = match declared {
            Some(total) => total > limit,
            None => next > limit,
        };
        if over {
            return Err(EngineError::job(format!(
                "BCF carried record grew to {next} bytes, past the {limit}-byte ceiling (the \
                 smaller of the hard limit and the --memory budget)"
            )));
        }
        Ok(())
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use reed_engine::budget::{MemBudget, ThreadBudget};
    use reed_engine::pipeline::ChunkJob;

    /// A job that records nothing: these tests only exercise framing.
    struct NoopJob;

    impl ChunkJob for NoopJob {
        type Partial = ();
        type Scratch = ();
        fn framing(&self) -> RecordFraming {
            RecordFraming::Newline
        }
        fn run_chunk(
            &self,
            _chunk: &ChunkBytes,
            _body: &[u8],
            _scratch: &mut Self::Scratch,
        ) -> Result<Self::Partial, EngineError> {
            Ok(())
        }
        fn run_seam(
            &self,
            _id: reed_core::ChunkId,
            _seam: &[u8],
            _scratch: &mut Self::Scratch,
        ) -> Result<Self::Partial, EngineError> {
            Ok(())
        }
    }

    /// A BCF record's bytes: an 8-byte length prefix declaring `total`, with
    /// only `present` bytes actually materialised. `present < total` models a
    /// record that starts in this chunk and continues in the next.
    fn record_prefixed(total: usize, present: usize) -> Vec<u8> {
        let mut out = Vec::new();
        // `declared_len` = 8 + l_shared + l_indiv; l_indiv counts the sample
        // block, which is empty for these framing-only fixtures.
        let l_shared = u32::try_from(total - 8).expect("invariant: fits");
        out.extend_from_slice(&l_shared.to_le_bytes());
        out.extend_from_slice(&0_u32.to_le_bytes());
        out.resize(present.min(total), 0x00);
        out
    }

    /// Builds an `InputJob` whose BCF state exists but whose dictionary is
    /// never consulted (framing only), with an explicit ceiling.
    fn job_with_ceiling(ceiling: usize) -> InputJob<NoopJob> {
        let abort = reed_engine::abort::abort();
        let chain = Arc::new(FrameChain::new(abort.signal.clone()));
        chain.publish(0, Carried::default());
        let dict = Arc::new(
            BcfDict::from_header(
                b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
            )
            .expect("invariant: minimal BCF header"),
        );
        InputJob {
            inner: NoopJob,
            record_ceiling: ceiling,
            bcf: Some(BcfState { chain, dict }),
            native: None,
        }
    }

    fn chunk(id: u64, compressed_start: u64, blocks: Vec<(u32, u32)>) -> ChunkBytes {
        ChunkBytes {
            id: reed_core::ChunkId(id),
            compressed: reed_io::chunk::ChunkPayload::from(Vec::new()),
            first_in_block_offset: 0,
            plain: true,
            compressed_start,
            blocks,
            truncate_to: None,
        }
    }

    /// Review item 2: when a record starts inside this chunk and only ends in
    /// a later one, the first fragment is published by the no-carry branch —
    /// which used to skip the ceiling check entirely, so a small `--memory`
    /// could publish a carried prefix of arbitrary size.
    #[test]
    fn the_first_crossing_fragment_is_checked_against_the_ceiling() {
        // A record whose prefix declares 4096 bytes but only 64 of them are
        // in this chunk: the carried fragment is 64 bytes, while the declared
        // total (4096) is what the ceiling must judge.
        let text = record_prefixed(4096, 64);
        let job = job_with_ceiling(1024);
        let chain = Arc::new(FrameChain::new(reed_engine::abort::abort().signal));
        chain.publish(0, Carried::default());
        let mut carried = Carried::default();
        let c = chunk(0, 0, vec![]);
        let err = job
            .frame_chunk(&c, &text, &mut carried, &chain)
            .expect_err("a declared 4096-byte record past a 1024-byte ceiling must fail");
        let msg = err.to_string();
        assert!(msg.contains("ceiling"), "{msg}");

        // The control: the same shape under a ceiling that admits it passes
        // and publishes the crossing fragment.
        let roomy = job_with_ceiling(1 << 20);
        let chain = Arc::new(FrameChain::new(reed_engine::abort::abort().signal));
        chain.publish(0, Carried::default());
        let mut carried = Carried::default();
        let out = roomy
            .frame_chunk(&c, &text, &mut carried, &chain)
            .expect("invariant: within the ceiling");
        // The record does not end in this chunk, so the framed range is empty
        // (`walk_records` returns the start when the first record does not
        // fit); its whole text becomes the carried prefix.
        assert_eq!(out, Some((0, 0)), "no complete record in this chunk");
        // The successor's cell carries the rest of the record; that is where
        // the published fragment is observable (the caller's `carried` is the
        // *predecessor's*, which is not mutated in this branch).
        let next = chain
            .peek(1)
            .expect("invariant: the successor cell was published");
        assert_eq!(next.prefix_len, 64, "the whole chunk is carried");
    }

    /// The real budget path: a tiny `--memory` must bound the first carried
    /// fragment too, not just the later accumulates.
    #[test]
    fn a_small_memory_budget_bounds_the_first_carried_fragment() {
        let budget = ThreadBudget {
            threads: 1,
            mem: MemBudget::from_bytes(1),
        };
        let ceiling = reed_engine::pipeline::record_ceiling(&budget);
        assert_eq!(ceiling, 1, "the ceiling follows the budget");

        // The declared total is 136 bytes, over the 1-byte ceiling, with 64
        // of them present: framing must refuse rather than carry it.
        let text = record_prefixed(136, 64);
        let job = job_with_ceiling(ceiling);
        let chain = Arc::new(FrameChain::new(reed_engine::abort::abort().signal));
        chain.publish(0, Carried::default());
        let mut carried = Carried::default();
        let c = chunk(0, 0, vec![]);
        assert!(
            job.frame_chunk(&c, &text, &mut carried, &chain).is_err(),
            "a record past the --memory-derived ceiling must be refused"
        );
    }
}
