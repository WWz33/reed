//! `TwoPassReparse` pass0 (`doc_v2/03-architecture.md` §8): the GT-only
//! per-sample fully-missing reduction that resolves the cohort rule.

use std::path::Path;
use std::sync::{Arc, Mutex};

use reed_core::genotype::{GenotypePolicies, PolicyCounters};
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{run, ChunkJob, OrderedSink};
use reed_io::vcf::genotype_scan::GenotypeScanner;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::site::SiteView;

use reed_core::ChunkId;
use reed_filter::FilterError;

use super::{ferr, io_err, FilterOptions};

/// Pass0 partial: per-sample fully-missing GT counts.
#[derive(Default)]
pub(super) struct MissingPartial {
    /// Records scanned (all of them, no head filtering).
    pub(super) n_records: u64,
    /// Per-sample fully-missing genotype counts.
    pub(super) counts: Vec<u64>,
}

/// The GT-only missing-count job (`TwoPassReparse` pass0).
pub(super) struct MissingJob {
    pub(super) n_samples: usize,
    pub(super) policies: GenotypePolicies,
}

impl ChunkJob for MissingJob {
    type Partial = MissingPartial;
    type Scratch = (Option<GenotypeScanner>, Vec<u64>, PolicyCounters);

    fn run_chunk(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let (scanner, miss, counters) = scratch;
        let n_samples = self.n_samples;
        let scanner = scanner.get_or_insert_with(|| {
            GenotypeScanner::new(n_samples, &[], &[], self.policies)
        });
        miss.clear();
        miss.resize(self.n_samples, 0);
        let mut n_records = 0_u64;
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            let Some(lr) = LineRef::parse(line) else {
                return Err(EngineError::Job("record line has fewer than 9 columns".into()));
            };
            let site = SiteView::parse(&lr).map_err(|e| EngineError::Job(e.to_string()))?;
            *counters = PolicyCounters::default();
            let out = scanner
                .scan(&lr, &site, counters)
                .map_err(|e| EngineError::Job(e.to_string()))?;
            for (i, &code) in out.codes.iter().enumerate() {
                if code.is_fully_missing() {
                    miss[i] += 1;
                }
            }
            n_records += 1;
        }
        Ok(MissingPartial { n_records, counts: miss.clone() })
    }

    fn run_seam(&self, _id: ChunkId, line: &[u8], scratch: &mut Self::Scratch) -> Result<Self::Partial, EngineError> {
        let (scanner, miss, counters) = scratch;
        let n_samples = self.n_samples;
        let scanner = scanner.get_or_insert_with(|| {
            GenotypeScanner::new(n_samples, &[], &[], self.policies)
        });
        miss.clear();
        miss.resize(self.n_samples, 0);
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::Job("seam record has fewer than 9 columns".into()));
        };
        let site = SiteView::parse(&lr).map_err(|e| EngineError::Job(e.to_string()))?;
        *counters = PolicyCounters::default();
        let out = scanner
            .scan(&lr, &site, counters)
            .map_err(|e| EngineError::Job(e.to_string()))?;
        for (i, &code) in out.codes.iter().enumerate() {
            if code.is_fully_missing() {
                miss[i] += 1;
            }
        }
        Ok(MissingPartial { n_records: 1, counts: miss.clone() })
    }
}

pub(super) fn run_missing_pass(
    input: Option<&Path>,
    n_samples: usize,
    policies: GenotypePolicies,
    opts: &FilterOptions,
) -> Result<MissingPartial, FilterError> {
    let budget = ThreadBudget { threads: opts.threads.max(1), mem: reed_engine::budget::MemBudget::default_half_memory() };
    let (mut source, _facts) = crate::input::open_regions(input, &opts.regions).map_err(io_err)?;
    let job = MissingJob { n_samples, policies };
    let total: Arc<Mutex<Option<MissingPartial>>> = Arc::new(Mutex::new(None));
    let sink = SumSink { total: Arc::clone(&total) };
    let abort = reed_engine::abort::abort();
    run(source.as_mut(), &job, sink, &budget, &abort).map_err(ferr)?;
    Ok(total
        .lock()
        .expect("invariant: missing slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the partial"))
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
                for (t, m) in total.counts.iter_mut().zip(&partial.counts) {
                    *t += m;
                }
            }
            None => {
                let mut counts = vec![0_u64; partial.counts.len()];
                for (t, m) in counts.iter_mut().zip(&partial.counts) {
                    *t += m;
                }
                *guard = Some(MissingPartial { n_records: partial.n_records, counts });
            }
        }
        Ok(())
    }

    fn finish(self) -> Result<u64, EngineError> {
        let total = self.total.lock().expect("invariant: missing slot not poisoned");
        Ok(total.as_ref().map_or(0, |t| t.n_records))
    }
}
