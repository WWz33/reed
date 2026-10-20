//! Reading a VCF into a dense genotype matrix (doc_v17).
//!
//! Both `reed pca` and `reed distance` need the same thing: every biallelic
//! site's dosage per sample, in input order. Sharing one implementation
//! avoids the "two parsers drift apart" failure mode.


use reed_io::vcf::line::LineRef;
use reed_engine::pipeline::{ChunkJob, OrderedSink};
use reed_core::ChunkId;
use reed_engine::error::EngineError;
use reed_engine::abort::Abort;
use reed_engine::budget::ThreadBudget;
/// The sites and dosages collected from the input.
///
/// Dosages are stored site-major (`n_sites x n_samples`), which is the layout
/// the facade expects after a transpose at the end.
#[derive(Default)]
pub(crate) struct GenotypeBlock {
    /// `(chrom, pos1)` per site, in input order.
    pub(crate) sites: Vec<(Vec<u8>, u64)>,
    /// Row-major dosages: `sites.len() * n_samples`, `NaN` for missing.
    pub(crate) dosages: Vec<f64>,
}

impl GenotypeBlock {
    pub(crate) fn merge(&mut self, other: Self) {
        self.sites.extend(other.sites);
        self.dosages.extend(other.dosages);
    }
}

#[derive(Default)]
pub(crate) struct PcaScratch {
    /// Per-sample GT token of the current record.
    gt: Vec<Option<u32>>,
    /// Allele count of the current record's ALT field.
    n_alt: usize,
}

pub(crate) struct PcaJob {
    pub(crate) n_samples: usize,
    pub(crate) keep: Vec<bool>,
    /// Only biallelic records contribute; multiallelic sites are skipped.
    pub(crate) allow_multiallelic: bool,
}

impl PcaJob {
    /// Parses the GT column of one sample into `0/1/2` or `NaN`.
    ///
    /// Returns `None` when the genotype is not a diploid biallelic call
    /// (missing, haploid, or an allele index beyond ALT).
    pub(crate) fn dosage(raw: &[u8], n_alt: usize) -> Option<f64> {
        let mut a = 0u32;
        let mut b = 0u32;
        let mut seen = 0usize;
        for token in raw.split(|&c| c == b'/' || c == b'|') {
            if seen >= 2 {
                return None;
            }
            if token == b"." || token.is_empty() {
                return None;
            }
            let mut value = 0u32;
            for &c in token {
                if !c.is_ascii_digit() {
                    return None;
                }
                value = value.checked_mul(10)?.checked_add(u32::from(c - b'0'))?;
            }
            if seen == 0 {
                a = value;
            } else {
                b = value;
            }
            seen += 1;
        }
        if seen != 2 || a as usize > n_alt || b as usize > n_alt {
            return None;
        }
        Some(f64::from(a + b))
    }
}

impl ChunkJob for PcaJob {
    type Partial = GenotypeBlock;
    type Scratch = PcaScratch;

    fn run_chunk(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let mut out = GenotypeBlock::default();
        for line in body.split(|&c| c == b'\n') {
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            let Some(lr) = LineRef::parse(line) else {
                continue;
            };
            if lr.n_sample_cols() != self.n_samples {
                continue;
            }
            let Ok(pos1) = std::str::from_utf8(lr.pos_text())
                .ok()
                .and_then(|s| s.parse::<u64>().ok())
                .ok_or(())
            else {
                continue;
            };
            // Allele count gate: only biallelic sites have a 0/1/2 dosage.
            let alts = lr.alts();
            let n_alt = if alts == b"." {
                0
            } else {
                alts.split(|&c| c == b',').count()
            };
            if n_alt == 0 || (n_alt > 1 && !self.allow_multiallelic) {
                continue;
            }
            scratch.n_alt = n_alt;

            let samples = lr.samples_text();
            scratch.gt.clear();
            let mut start = 0usize;
            let mut values = Vec::with_capacity(self.n_samples);
            for i in 0..self.n_samples {
                let begin = start.min(samples.len());
                let end = samples[begin..]
                    .iter()
                    .position(|&c| c == b'\t')
                    .map_or(samples.len(), |rel| begin + rel);
                start = end + 1;
                if !self.keep[i] {
                    continue;
                }
                let column = &samples[begin..end];
                // The GT is the first colon-separated field.
                let colon = column.iter().position(|&c| c == b':').unwrap_or(column.len());
                let raw = &column[..colon];
                values.push(PcaJob::dosage(raw, n_alt).unwrap_or(f64::NAN));
            }
            if values.is_empty() {
                continue;
            }
            // Skip sites with no variation at all: they contribute nothing and
            // would make the standardization degenerate.
            let mut any_alt = false;
            for v in &values {
                if !v.is_nan() && *v > 0.0 {
                    any_alt = true;
                    break;
                }
            }
            if !any_alt {
                continue;
            }
            out.sites.push((lr.chrom().to_vec(), pos1));
            out.dosages.extend_from_slice(&values);
        }
        Ok(out)
    }

    fn run_seam(
        &self,
        _id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        // A record spanning a chunk boundary: process it like any other line.
        let mut out = GenotypeBlock::default();
        if line.is_empty() || line[0] == b'#' {
            return Ok(out);
        }
        let Some(lr) = LineRef::parse(line) else {
            return Ok(out);
        };
        if lr.n_sample_cols() != self.n_samples {
            return Ok(out);
        }
        let Ok(pos1) = std::str::from_utf8(lr.pos_text())
            .ok()
            .and_then(|s| s.parse::<u64>().ok())
            .ok_or(())
        else {
            return Ok(out);
        };
        let alts = lr.alts();
        let n_alt = if alts == b"." {
            0
        } else {
            alts.split(|&c| c == b',').count()
        };
        if n_alt == 0 || (n_alt > 1 && !self.allow_multiallelic) {
            return Ok(out);
        }
        scratch.n_alt = n_alt;
        let samples = lr.samples_text();
        let mut start = 0usize;
        let mut values = Vec::with_capacity(self.n_samples);
        for i in 0..self.n_samples {
            let begin = start.min(samples.len());
            let end = samples[begin..]
                .iter()
                .position(|&c| c == b'\t')
                .map_or(samples.len(), |rel| begin + rel);
            start = end + 1;
            if !self.keep[i] {
                continue;
            }
            let column = &samples[begin..end];
            let colon = column.iter().position(|&c| c == b':').unwrap_or(column.len());
            values.push(PcaJob::dosage(&column[..colon], n_alt).unwrap_or(f64::NAN));
        }
        let mut any_alt = false;
        for v in &values {
            if !v.is_nan() && *v > 0.0 {
                any_alt = true;
                break;
            }
        }
        if values.is_empty() || !any_alt {
            return Ok(out);
        }
        out.sites.push((lr.chrom().to_vec(), pos1));
        out.dosages.extend_from_slice(&values);
        Ok(out)
    }
}

/// Collects the blocks in chunk order into one matrix.
pub(crate) struct CollectSink {
    pub(crate) acc: GenotypeBlock,
}

impl OrderedSink<GenotypeBlock> for CollectSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut GenotypeBlock) -> Result<(), EngineError> {
        self.acc.merge(std::mem::take(partial));
        Ok(())
    }

    fn finish(self) -> Result<u64, EngineError> {
        Ok(self.acc.sites.len() as u64)
    }
}

pub(crate) fn collect_via_engine(
    opened: &mut crate::OpenedInput,
    job: &PcaJob,
    sink: CollectSink,
    budget: &ThreadBudget,
    abort: &Abort,
) -> Result<GenotypeBlock, crate::Error> {
    // The engine's sink consumes the accumulator, so stash the result in a
    // shared cell while `finish` runs.
    use std::sync::{Arc, Mutex};
    let slot: Arc<Mutex<Option<GenotypeBlock>>> = Arc::new(Mutex::new(None));

    struct StashSink {
        inner: CollectSink,
        slot: Arc<Mutex<Option<GenotypeBlock>>>,
    }
    impl OrderedSink<GenotypeBlock> for StashSink {
        fn accept(
            &mut self,
            id: ChunkId,
            partial: &mut GenotypeBlock,
        ) -> Result<(), EngineError> {
            self.inner.accept(id, partial)
        }
        fn finish(self) -> Result<u64, EngineError> {
            let n = self.inner.acc.sites.len() as u64;
            *self.slot.lock().expect("invariant: sink slot") = Some(self.inner.acc);
            Ok(n)
        }
    }

    let stash = StashSink { inner: sink, slot: Arc::clone(&slot) };
    reed_engine::pipeline::run(opened.source_mut(), job, stash, budget, abort)
        .map_err(crate::Error::from)?;
    let block = slot
        .lock()
        .expect("invariant: sink slot")
        .take()
        .unwrap_or_default();
    Ok(block)
}
