//! `run_qc`: the first engine consumer — the QC `summary` section over a
//! BGZF VCF (`docs/architecture-spec.md` §10.1, M3 first consumer).

use std::path::Path;

use reed_core::site::classify_variant;
use reed_engine::budget::ThreadBudget;
use reed_engine::pipeline::{run, ChunkJob, OrderedSink};
use reed_io::chunk::BgzfChunkSource;
use reed_qc::summary::SummaryAcc;

/// Result of one `run_qc` call.
#[derive(Clone, Debug)]
pub struct SummaryReport {
    /// Samples in the input header.
    pub n_samples: usize,
    /// The `summary` section counts.
    pub summary: SummaryAcc,
    /// Chunks processed.
    pub n_chunks: u64,
    /// Wall time.
    pub elapsed: std::time::Duration,
}

/// The QC summary job: fixed columns only — the sample region is never
/// touched (§6.3: head predicates cost no sample-region parsing).
pub struct QcJob;

impl ChunkJob for QcJob {
    type Partial = SummaryAcc;
    type Scratch = ();

    fn run_chunk(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        _scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, reed_engine::error::EngineError> {
        let mut acc = Self::Partial::default();
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            push_line(&mut acc, line)?;
        }
        Ok(acc)
    }

    fn run_seam(
        &self,
        _id: reed_core::ChunkId,
        line: &[u8],
        _scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, reed_engine::error::EngineError> {
        let mut acc = Self::Partial::default();
        push_line(&mut acc, line)?;
        Ok(acc)
    }
}

fn push_line(
    acc: &mut SummaryAcc,
    line: &[u8],
) -> Result<(), reed_engine::error::EngineError> {
    let Some(line_ref) = reed_io::vcf::line::LineRef::parse(line) else {
        return Err(reed_engine::error::EngineError::Job(
            "record line has fewer than 9 columns".into(),
        ));
    };
    let reference = line_ref.reference();
    let alts = line_ref.alts();
    let (primary_alt, n_alt) = match alts {
        b"." => (alts, 0_u8),
        _ => {
            let first_end = memchr::memchr(b',', alts).unwrap_or(alts.len());
            let n = alts.iter().filter(|&&b| b == b',').count() + 1;
            (&alts[..first_end], u8::try_from(n).unwrap_or(u8::MAX))
        }
    };
    acc.push(classify_variant(reference, primary_alt), n_alt);
    Ok(())
}

/// Sink collecting the merged summary; hands the accumulator back through a
/// shared slot at `finish` (the engine consumes the sink by value).
struct SummarySink {
    acc: SummaryAcc,
    out: std::sync::Arc<std::sync::Mutex<Option<SummaryAcc>>>,
}

impl OrderedSink<SummaryAcc> for SummarySink {
    fn accept(
        &mut self,
        _id: reed_core::ChunkId,
        partial: &mut SummaryAcc,
    ) -> Result<(), reed_engine::error::EngineError> {
        self.acc.merge(partial);
        Ok(())
    }

    fn finish(mut self) -> Result<u64, reed_engine::error::EngineError> {
        *self
            .out
            .lock()
            .expect("invariant: summary slot not poisoned") = Some(self.acc);
        Ok(self.acc.n_records)
    }
}

/// Runs the QC summary section over a BGZF-compressed VCF.
///
/// # Errors
///
/// Anything the transport or the job rejects, wrapped from
/// [`reed_engine::error::EngineError`].
pub fn run_qc(
    input: &Path,
    budget: &ThreadBudget,
) -> Result<SummaryReport, Box<dyn std::error::Error + Send + Sync>> {
    let file = std::fs::File::open(input)?;
    let (mut source, header) = BgzfChunkSource::open(std::io::BufReader::with_capacity(1 << 20, file))?;
    let n_samples = header.n_samples();
    let abort = reed_engine::abort::abort();
    let slot = std::sync::Arc::new(std::sync::Mutex::new(None));
    let report = run(
        &mut source as &mut dyn reed_io::chunk::ChunkSource,
        &QcJob,
        SummarySink {
            acc: SummaryAcc::default(),
            out: std::sync::Arc::clone(&slot),
        },
        budget,
        &abort,
    )?;
    let summary = slot
        .lock()
        .expect("invariant: summary slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the accumulator");
    Ok(SummaryReport {
        n_samples,
        summary,
        n_chunks: report.n_chunks,
        elapsed: report.elapsed,
    })
}
