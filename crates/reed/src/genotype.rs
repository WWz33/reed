//! Reading a VCF into a dense genotype matrix (doc_v17).
//!
//! Both `reed pca` and `reed distance` need the same thing: every biallelic
//! site's dosage per sample, in input order. Sharing one implementation
//! avoids the "two parsers drift apart" failure mode.


use reed_io::vcf::line::LineRef;
use reed_io::{IoError, Locus};
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
    /// Sites collected, in input order. Only the count is consumed
    /// downstream (doc_v20 R6: the per-site `(chrom, pos)` vectors were
    /// dead allocations).
    pub(crate) n_sites: usize,
    /// Row-major dosages: `n_sites * n_samples`, `NaN` for missing.
    pub(crate) dosages: Vec<f64>,
}

impl GenotypeBlock {
    pub(crate) fn merge(&mut self, other: Self) {
        self.n_sites += other.n_sites;
        self.dosages.extend(other.dosages);
    }
}

/// A malformed record line, mirroring the scanner's `IoError::MalformedRecord`
/// (`genotype_scan.rs`). The locus is built from whatever columns are
/// present, with `pos1 = 0` when POS does not parse (the `parse_pos1`
/// convention) — the whole line is not re-parsed just to name it.
fn malformed_line(line: &[u8], reason: &'static str) -> EngineError {
    let chrom_end = line.iter().position(|&b| b == b'\t').unwrap_or(line.len());
    let chrom = &line[..chrom_end];
    let pos_start = (chrom_end + 1).min(line.len());
    let pos_end = line[pos_start..]
        .iter()
        .position(|&b| b == b'\t')
        .map_or(line.len(), |rel| pos_start + rel);
    let pos1 = std::str::from_utf8(&line[pos_start..pos_end])
        .ok()
        .and_then(|t| t.parse().ok())
        .unwrap_or(0);
    EngineError::from(IoError::MalformedRecord {
        locus: Locus {
            contig_name: String::from_utf8_lossy(chrom).into_owned(),
            pos1,
        },
        reason,
    })
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
            // A malformed record is an error, not a skip: silently dropping
            // it yields a wrong matrix on a truncated or mismatched file
            // (doc_v20 R6).
            let Some(lr) = LineRef::parse(line) else {
                return Err(malformed_line(line, "record line has fewer than 9 columns"));
            };
            if lr.n_sample_cols() != self.n_samples {
                return Err(malformed_line(
                    line,
                    "sample column count disagrees with the header",
                ));
            }
            // doc_v20 review of R6: the position is not stored any more, but
            // a non-numeric POS is still a malformed record and must not pass
            // silently here while `filter`/`stats` reject the same line.
            if !lr.pos_text().iter().all(u8::is_ascii_digit) || lr.pos_text().is_empty() {
                return Err(malformed_line(line, "POS is not a positive integer"));
            }
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
            out.n_sites += 1;
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
        // Malformed seam records are errors too, like the body path
        // (doc_v20 R6).
        let Some(lr) = LineRef::parse(line) else {
            return Err(malformed_line(line, "record line has fewer than 9 columns"));
        };
        if lr.n_sample_cols() != self.n_samples {
            return Err(malformed_line(
                line,
                "sample column count disagrees with the header",
            ));
        }
        // doc_v20 review of R6: same POS contract as `run_chunk`.
        if !lr.pos_text().iter().all(u8::is_ascii_digit) || lr.pos_text().is_empty() {
            return Err(malformed_line(line, "POS is not a positive integer"));
        }
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
        out.n_sites += 1;
        out.dosages.extend_from_slice(&values);
        Ok(out)
    }
}

/// Collects the blocks in chunk order into one matrix.
pub(crate) struct CollectSink {
    pub(crate) acc: GenotypeBlock,
    /// Upper bound on the collected dosage bytes (doc_v20 review of S76).
    /// The block is accumulated site by site, so checking only after the
    /// scan (as the first version of the fix did) still lets the full
    /// allocation happen before the refusal. `None` means unbounded.
    pub(crate) max_bytes: Option<u64>,
    /// Samples per site in the collected layout.
    pub(crate) n_samples: usize,
    /// Set when the cap was hit; the caller turns it into a usage error
    /// after the scan aborts (an `EngineError` from the sink would be
    /// classified as a runtime failure, not a command-line one).
    pub(crate) exceeded: std::sync::Arc<std::sync::atomic::AtomicBool>,
}

impl CollectSink {
    /// Bytes the collected dosages occupy once `n_sites` sites are in.
    fn bytes_for(&self, n_sites: usize) -> u64 {
        (n_sites as u64)
            .saturating_mul(self.n_samples as u64)
            .saturating_mul(std::mem::size_of::<f64>() as u64)
    }

    fn over_budget(&self, n_sites: usize) -> bool {
        self.max_bytes.is_some_and(|max| self.bytes_for(n_sites) > max)
    }
}

impl OrderedSink<GenotypeBlock> for CollectSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut GenotypeBlock) -> Result<(), EngineError> {
        self.acc.merge(std::mem::take(partial));
        if self.over_budget(self.acc.n_sites) {
            // Stop the scan now; `collect_via_engine` reports the usage
            // error after the abort unwinds.
            self.exceeded
                .store(true, std::sync::atomic::Ordering::Release);
            return Err(EngineError::job("dosage matrix exceeded --memory"));
        }
        Ok(())
    }

    fn finish(self) -> Result<u64, EngineError> {
        Ok(self.acc.n_sites as u64)
    }
}

pub(crate) fn collect_via_engine(
    opened: &mut crate::OpenedInput,
    job: PcaJob,
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
            let n = self.inner.acc.n_sites as u64;
            *self.slot.lock().expect("invariant: sink slot") = Some(self.inner.acc);
            Ok(n)
        }
    }

    let stash = StashSink { inner: sink, slot: Arc::clone(&slot) };
    // doc_v21 M1: every other entry point wraps its job in `InputJob`, which
    // renders BCF records to text; a bare `PcaJob` has `RecordFraming::Newline`
    // and splits raw BCF bytes on '\n' (documented-input command failure). Wrap
    // unconditionally — on text input `InputJob` is a transparent passthrough.
    let wrapped = crate::bcf_input::InputJob::new(job, opened, abort);
    let report = reed_engine::pipeline::run(opened.source_mut(), &wrapped, stash, budget, abort)
        .map_err(crate::Error::from)?;
    wrapped
        .finish(report.n_chunks, true)
        .map_err(crate::Error::from)?;
    let block = slot
        .lock()
        .expect("invariant: sink slot")
        .take()
        .unwrap_or_default();
    Ok(block)
}

#[cfg(test)]
mod tests {
    use super::*;
    use reed_io::chunk::{ChunkBytes, ChunkPayload};

    fn job(n_samples: usize) -> PcaJob {
        PcaJob {
            n_samples,
            keep: vec![true; n_samples],
            allow_multiallelic: false,
        }
    }

    /// `run_chunk` ignores the chunk bytes themselves; a placeholder is all
    /// the in-process call needs.
    fn chunk_bytes() -> ChunkBytes {
        ChunkBytes {
            id: ChunkId(0),
            compressed: ChunkPayload::Owned(vec![].into()),
            first_in_block_offset: 0,
            plain: true,
            compressed_start: 0,
            blocks: Vec::new(),
            truncate_to: None,
        }
    }

    /// doc_v20 R6: a record whose sample-column count disagrees with the
    /// header must fail the run (mirroring the scanner's
    /// `IoError::MalformedRecord`), not be skipped — the skip silently
    /// shrinks the PCA/distance matrix on a truncated or mismatched file.
    #[test]
    fn r6_run_chunk_errors_on_sample_column_mismatch() {
        let job = job(2);
        let mut scratch = PcaScratch::default();
        // One sample column where the header declared two.
        let body = b"chr1\t101\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n";
        let Err(err) = job.run_chunk(&chunk_bytes(), body, &mut scratch) else {
            panic!("a sample-column mismatch must fail the chunk, not be skipped (doc_v20 R6)")
        };
        let msg = err.to_string();
        assert!(msg.contains("chr1:101"), "error must name the locus: {msg}");
        assert!(
            msg.contains("sample column count"),
            "error must name the mismatch: {msg}"
        );
    }

    /// doc_v20 R6: a truncated record line (fewer than nine tabs) must fail
    /// the run, not be skipped.
    #[test]
    fn r6_run_chunk_errors_on_unparseable_line() {
        let job = job(2);
        let mut scratch = PcaScratch::default();
        let body = b"chr1\t102\t.\tA\n";
        let Err(err) = job.run_chunk(&chunk_bytes(), body, &mut scratch) else {
            panic!("an unparseable line must fail the chunk (doc_v20 R6)")
        };
        let msg = err.to_string();
        assert!(msg.contains("chr1:102"), "error must name the locus: {msg}");
        assert!(msg.contains("9 columns"), "error must name the reason: {msg}");
    }

    /// doc_v20 R6: the seam path (a record spanning a chunk boundary)
    /// rejects the same malformed lines as the body path.
    #[test]
    fn r6_run_chunk_errors_on_unparseable_pos() {
        // doc_v20 review of R6: dropping the POS parse also dropped the
        // check. `filter`/`stats` reject this line; pca/distance must too.
        let job = job(1);
        let body = b"chr1\tX\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\n";
        let Err(err) = job.run_chunk(&chunk_bytes(), body, &mut PcaScratch::default()) else {
            panic!("a non-numeric POS must fail the chunk");
        };
        assert!(err.to_string().contains("POS"), "{err}");
    }

    #[test]
    fn r6_run_seam_errors_on_sample_column_mismatch() {
        let job = job(2);
        let mut scratch = PcaScratch::default();
        let line = b"chr1\t103\t.\tA\tG\t50\tPASS\t.\tGT\t0/1";
        let Err(err) = job.run_seam(ChunkId(1), line, &mut scratch) else {
            panic!("a sample-column mismatch must fail the seam (doc_v20 R6)")
        };
        assert!(
            err.to_string().contains("chr1:103"),
            "error must name the locus: {err}"
        );
    }
}
