//! The sketch path's pass0 (`docs/rounds/v2/03-architecture.md` §8, §11.4):
//! head rules → per-sample GT planes + per-rule mask bitmaps (head
//! survivors) + per-sample post-mask called-allele counts over head
//! survivors (`--mind`'s vcftools semantics, doc_v4 02 回退版 §3.2).

use std::path::Path;
use std::sync::{Arc, Mutex};

use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{run, ChunkJob, OrderedSink};
use reed_filter::funnel::FunnelPartial;
use reed_filter::plan::FilterPlan;
use reed_filter::FilterError;

use reed_core::ChunkId;

use super::record_pipeline::{PipelineOptions, RecordPipeline};
use super::{ferr, io_err, FilterOptions};

pub(super) struct SketchPassOutput {
    /// Per-chunk sketches in chunk order.
    pub sketches: Vec<ChunkSketch>,
    /// Genome-wide per-sample post-mask called-allele counts (head survivors).
    pub cohort_called: Vec<u64>,
    /// Head-surviving records (the `2L` denominator, body + seams).
    pub n_records: u64,
    /// The pass0 funnel (head-rule hits; tail/mask counters still zero —
    /// the mid stage fills them).
    pub funnel: FunnelPartial,
    /// Pass0's engine metrics (doc_v4 03 A4-5: the driver sums the passes).
    pub metrics: reed_engine::metrics::StageMetrics,
}

/// One chunk's sketch plus where its trailing seam site landed.
pub(super) struct ChunkSketch {
    pub id: u64,
    pub fragment: reed_filter::sketch::SketchFragment,
    /// Site index of the seam record within the fragment (the seam is the
    /// fragment's last record when it passed the head rules).
    pub seam_site: Option<usize>,
}

/// Worker scratch of the sketch pass: the shared record pipeline plus the
/// per-site planes and per-rule mask bitmaps the sketch carries.
pub(super) struct SketchScratch {
    pipeline: RecordPipeline,
    p0a: Vec<u64>,
    p0m: Vec<u64>,
    p1a: Vec<u64>,
    p1m: Vec<u64>,
    rule_bits: Vec<Vec<u64>>,
    called: Vec<u64>,
    funnel: FunnelPartial,
}

/// The pass0 job: head rules → per-sample GT planes + per-rule mask
/// bitmaps (head survivors) + per-sample post-mask called-allele counts
/// (head survivors, the `--mind` numerator).
pub(super) struct SketchJob {
    pub(super) plan: Arc<FilterPlan>,
    pub(super) n_samples: usize,
    pub(super) w: usize,
    pub(super) contig_names: Arc<[Vec<u8>]>,
}

/// The sketch pass's chunk partial. `pre` carries the seam record's
/// fragment — it belongs to the *previous* chunk's stream.
#[derive(Default)]
pub(super) struct SketchPartial {
    pre: Option<reed_filter::sketch::SketchFragment>,
    fragment: reed_filter::sketch::SketchFragment,
    cohort_called: Vec<u64>,
    /// Head-surviving records in this chunk.
    n_head_passed: u64,
    funnel: FunnelPartial,
}
impl ChunkJob for SketchJob {
    type Partial = SketchPartial;
    type Scratch = Option<SketchScratch>;

    fn run_chunk(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        let mut fragment = reed_filter::sketch::SketchFragment::new(self.w, self.plan.mask.len());
        scratch.funnel.reset(
            self.plan.head.len(),
            self.plan.tail.len(),
            self.plan.mask.len(),
        );
        scratch.called.clear();
        scratch.called.resize(self.n_samples, 0);
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            self.sketch_line(line, scratch, &mut fragment)?;
        }
        let n_head_passed = head_survivors(self.plan.head.len(), &scratch.funnel);
        Ok(SketchPartial {
            pre: None,
            fragment,
            cohort_called: scratch.called.clone(),
            n_head_passed,
            funnel: std::mem::take(&mut scratch.funnel),
        })
    }

    fn run_seam(
        &self,
        _id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| self.build_scratch());
        let mut pre = reed_filter::sketch::SketchFragment::new(self.w, self.plan.mask.len());
        scratch.funnel.reset(
            self.plan.head.len(),
            self.plan.tail.len(),
            self.plan.mask.len(),
        );
        scratch.called.clear();
        scratch.called.resize(self.n_samples, 0);
        self.sketch_line(line, scratch, &mut pre)?;
        let n_head_passed = head_survivors(self.plan.head.len(), &scratch.funnel);
        Ok(SketchPartial {
            pre: Some(pre),
            fragment: reed_filter::sketch::SketchFragment::new(self.w, self.plan.mask.len()),
            cohort_called: scratch.called.clone(),
            n_head_passed,
            funnel: std::mem::take(&mut scratch.funnel),
        })
    }
}

impl SketchJob {
    fn build_scratch(&self) -> SketchScratch {
        let int_keys: Vec<&[u8]> = self.plan.int_fields.iter().map(Vec::as_slice).collect();
        let float_keys: Vec<&[u8]> = self.plan.float_fields.iter().map(Vec::as_slice).collect();
        let mut pipeline = RecordPipeline::new(
            self.n_samples,
            &int_keys,
            &float_keys,
            self.plan.policies,
            &self.contig_names,
        );
        pipeline.resize_for(&self.plan);
        SketchScratch {
            pipeline,
            p0a: vec![0; self.w],
            p0m: vec![0; self.w],
            p1a: vec![0; self.w],
            p1m: vec![0; self.w],
            rule_bits: vec![vec![0; self.w]; self.plan.mask.len()],
            called: vec![0; self.n_samples],
            funnel: FunnelPartial {
                head_hits: vec![0; self.plan.head.len()],
                tail_hits: vec![0; self.plan.tail.len()],
                mask_hits: vec![0; self.plan.mask.len()],
                ..FunnelPartial::default()
            },
        }
    }

    /// Sketches one record through the shared pipeline: the head rules run
    /// first (the funnel counts every rule passed); head survivors
    /// contribute a site (planes + per-rule mask bitmaps over all samples —
    /// the mid stage intersects with `keep`); every head survivor feeds the
    /// cohort's called-allele counts (§11.4 pass0).
    fn sketch_line(
        &self,
        line: &[u8],
        scratch: &mut SketchScratch,
        fragment: &mut reed_filter::sketch::SketchFragment,
    ) -> Result<(), EngineError> {
        let plan = &self.plan;
        let SketchScratch { pipeline, p0a, p0m, p1a, p1m, rule_bits, called, funnel } = scratch;
        // The scan runs for EVERY record (the denominator counts every
        // head-surviving record); the plane/mask bitmaps and the cohort
        // counts are only filled for head survivors. The tail is the mid
        // stage's job.
        let opts = PipelineOptions {
            scan_head_failures: true,
            mask_all_samples: true,
            require_contig: false,
            tail_n_samples: None,
            counts: false,
        };
        let keep_all: Vec<bool> = Vec::new();
        let processed = pipeline.process(line, plan, &keep_all, opts, funnel, |view| {
            let head_passed = view.head_passed;
            if !head_passed {
                return Ok(());
            }
            for v in [&mut *p0a, &mut *p0m, &mut *p1a, &mut *p1m] {
                v.iter_mut().for_each(|word| *word = 0);
            }
            for v in rule_bits.iter_mut() {
                v.iter_mut().for_each(|word| *word = 0);
            }
            for (i, &code) in view.scan.codes.iter().enumerate() {
                // `--mind` (doc_v4 02 回退版 §3.2): post-mask called
                // alleles over the head survivors.
                let n = super::reduce::called_alleles(code, view.masked.get(i));
                if let Some(slot) = called.get_mut(i) {
                    *slot += n;
                }
                let (word, bit) = (i / 64, 1_u64 << (i % 64));
                if code.slot_alt(0) {
                    p0a[word] |= bit;
                }
                if code.slot_missing(0) {
                    p0m[word] |= bit;
                }
                if code.slot_alt(1) {
                    p1a[word] |= bit;
                }
                if code.slot_missing(1) {
                    p1m[word] |= bit;
                }
                if let Some(rule) = view.mask_rule.get(i).copied().filter(|&r| r != u16::MAX) {
                    rule_bits[usize::from(rule)][word] |= bit;
                }
            }
            Ok(())
        })?;
        if processed.head_passed {
            fragment.push_site(
                &[p0a, p0m, p1a, p1m],
                &rule_bits.iter().map(Vec::as_slice).collect::<Vec<_>>(),
            );
        }
        fragment.push_record(processed.head_passed);
        Ok(())
    }
}

/// Records that survived every head rule (the funnel's cumulative
/// semantics: the last rule's count, or the input count with no head rules).
fn head_survivors(n_head: usize, funnel: &FunnelPartial) -> u64 {
    if n_head == 0 {
        funnel.n_input
    } else {
        funnel.head_hits.get(n_head - 1).copied().unwrap_or(0)
    }
}

/// Collects chunk sketches in order; the seam fragment closes the
/// previous chunk (§11.4: the seam record belongs to the previous
pub(super) struct SketchSink {
    sketches: Vec<ChunkSketch>,
    cohort: Vec<u64>,
    n_records: u64,
    funnel: FunnelPartial,
    open: Option<(u64, reed_filter::sketch::SketchFragment)>,
    out: Arc<Mutex<Option<SketchPassOutput>>>,
}

impl OrderedSink<SketchPartial> for SketchSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut SketchPartial) -> Result<(), EngineError> {
        self.funnel.merge(&partial.funnel);
        self.n_records += partial.n_head_passed;
        if self.cohort.len() < partial.cohort_called.len() {
            self.cohort.resize(partial.cohort_called.len(), 0);
        }
        for (dst, src) in self.cohort.iter_mut().zip(&partial.cohort_called) {
            *dst += src;
        }
        if let Some(pre) = partial.pre.take() {
            // The seam record closes the previously opened chunk.
            if let Some((cid, mut frag)) = self.open.take() {
                let seam_site = (pre.n_sites > 0).then_some(frag.n_sites);
                frag.append(&pre);
                self.sketches.push(ChunkSketch { id: cid, fragment: frag, seam_site });
            }
        }
        // Open the new chunk with its body fragment (a stale open without a
        // seam would mean a missing seam partial — finalize defensively).
        if let Some((cid, frag)) = self.open.take() {
            self.sketches.push(ChunkSketch { id: cid, fragment: frag, seam_site: None });
        }
        self.open = Some((_id.0, std::mem::take(&mut partial.fragment)));
        Ok(())
    }

    fn finish(mut self) -> Result<u64, EngineError> {
        if let Some((cid, frag)) = self.open.take() {
            self.sketches.push(ChunkSketch { id: cid, fragment: frag, seam_site: None });
        }
        let n_records = self.n_records;
        *self.out.lock().expect("invariant: sketch slot not poisoned") = Some(SketchPassOutput {
            sketches: self.sketches,
            cohort_called: self.cohort,
            n_records,
            funnel: self.funnel,
            // Filled by `run_sketch_pass` from the engine's RunReport.
            metrics: reed_engine::metrics::StageMetrics::default(),
        });
        Ok(n_records)
    }
}

/// The sketch pass (engine run #1).
pub(super) fn run_sketch_pass(
    opened: crate::input::OpenedInput,
    n_samples: usize,
    plan: &FilterPlan,
    opts: &FilterOptions,
) -> Result<SketchPassOutput, FilterError> {
    let budget = ThreadBudget {
        threads: opts.threads.max(1),
        mem: opts.mem_bytes.map_or_else(
            reed_engine::budget::MemBudget::default_half_memory,
            reed_engine::budget::MemBudget::from_bytes,
        ),
    };
    let crate::input::OpenedInput { mut source, facts, .. } = opened;
    let contig_names: Arc<[Vec<u8>]> = Arc::from(
        facts
            .header
            .contigs()
            .iter()
            .map(|(name, _)| name.as_bytes().to_vec())
            .collect::<Vec<_>>()
            .into_boxed_slice(),
    );
    let job = SketchJob {
        plan: Arc::new(plan.clone()),
        n_samples,
        w: n_samples.div_ceil(64),
        contig_names,
    };
    let out: Arc<Mutex<Option<SketchPassOutput>>> = Arc::new(Mutex::new(None));
    let sink = SketchSink {
        sketches: Vec::new(),
        cohort: Vec::new(),
        n_records: 0,
        funnel: sized_funnel(plan),
        open: None,
        out: Arc::clone(&out),
    };
    let abort = reed_engine::abort::abort();
    run(source.as_mut(), &job, sink, &budget, &abort).map_err(ferr)?;
    Ok(out
        .lock()
        .expect("invariant: sketch slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the sketch output"))
}

/// A sized funnel (the sketch sink's counters must exist before merges).
pub(super) fn sized_funnel(plan: &FilterPlan) -> FunnelPartial {
    FunnelPartial {
        head_hits: vec![0; plan.head.len()],
        tail_hits: vec![0; plan.tail.len()],
        mask_hits: vec![0; plan.mask.len()],
        ..FunnelPartial::default()
    }
}

pub(super) fn estimate_records(input: &Path) -> Result<u64, FilterError> {
    let file_len = std::fs::metadata(input).map_err(FilterError::from)?.len();
    let mut opened = crate::input::open(Some(input)).map_err(io_err)?;
    let Some(chunk) = opened.source.next_chunk().map_err(io_err)? else {
        return Ok(0);
    };
    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
    let mut scratch = Vec::new();
    let view = reed_io::chunk::inflate_chunk(&chunk, &mut inflater, &mut scratch).map_err(io_err)?;
    let body = view.body_bytes(view.text);
    let n_lines = u64::try_from(body.iter().filter(|&&b| b == b'\n').count()).unwrap_or(0);
    let compressed = u64::try_from(chunk.compressed.len()).unwrap_or(1).max(1);
    if n_lines == 0 {
        return Ok(0);
    }
    Ok(file_len.saturating_mul(n_lines) / compressed)
}
