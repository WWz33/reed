//! The convert chunk job (`ConvertJob`), its partials, the ordered sink and
//! the `run_convert` / `run_convert_opened` driver entry points
//! (docs/rounds/v13/02-architecture.md §2, §7).
//!
//! One line parse per record feeds every enabled target: the phylogenetic
//! column fold (vcf2phylip semantics), the binary-NEXUS stream, the `--012`
//! count matrix, the TASSEL hapmap rows and the reference-driven MSA sites.
//! The projected-SFS target composes the statistics chunk path (`StatsJob`,
//! v13 §7.5 plan A), so the SFS fold keeps exactly one implementation.

// Casts in this file are bounded by construction: sample/allele counts and
// offsets inside one chunk body are far below u32::MAX (the engine's record
// ceiling is 1 GiB).
#![allow(clippy::cast_possible_truncation)]

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};
use std::time::Duration;

use reed_core::iupac;
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::metrics::StageMetrics;
use reed_engine::pipeline::{ChunkJob, OrderedSink, run};
use reed_io::vcf::line::{LineRef, RecordClass, classify_line};

use crate::convert::writers;
use crate::convert::{ConvertRequest, ConvertTarget, MsaRegion, PhyloFormats, PhyloOptions};
use crate::stats::driver::{build_sfs_plan, resolve_groups, sfs_report};
use crate::stats::request::{SfsPartial, SfsReport};
use crate::stats::run::{StatsJob, StatsPartial, StatsScratch};

pub(crate) use super::spill::{RowsBlob, SPILL_THRESHOLD, SPILL_TOTAL_BUDGET, SpilledMsa};

/// One conversion run's record counters (the vcf2phylip stdout report,
/// v13 §3.1; reed adds `malformed`, B13-2).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct ConvertCounters {
    /// Every processed record line (vcf2phylip's `snp_num`).
    pub total_records: u64,
    /// Records below the `--min-samples-locus` gate (`snp_shallow`).
    pub excess_missing: u64,
    /// Records failing the SNP gate (`mnp_num`).
    pub mnps: u64,
    /// Nucleotide columns accepted (`snp_accepted`).
    pub accepted: u64,
    /// Binary-NEXUS columns accepted (`snp_biallelic`).
    pub biallelic_for_bin: u64,
    /// Skipped malformed records: column-count mismatches, bad genotype
    /// tokens, unknown allele letters or indices (B13-2/B13-9).
    pub malformed: u64,
}

impl ConvertCounters {
    fn absorb(&mut self, other: &Self) {
        self.total_records += other.total_records;
        self.excess_missing += other.excess_missing;
        self.mnps += other.mnps;
        self.accepted += other.accepted;
        self.biallelic_for_bin += other.biallelic_for_bin;
        self.malformed += other.malformed;
    }
}

/// The run's output description: the files written and the counters.
#[derive(Debug)]
pub struct ConvertReport {
    /// Every file written, in write order.
    pub files: Vec<PathBuf>,
    /// The record counters; `None` when no row target ran (an SFS-only run
    /// carries statistics instead).
    pub counters: Option<ConvertCounters>,
    /// The projected spectra when the SFS target ran.
    pub sfs: Option<SfsReport>,
    /// Records processed (the statistics count when SFS ran, else the
    /// row-target total).
    pub n_records: u64,
    /// Total wall time.
    pub elapsed: Duration,
    /// Pipeline stage metrics.
    pub metrics: StageMetrics,
}

// ---------------------------------------------------------------------------
// Plan structures (shared by the fold and the render)
// ---------------------------------------------------------------------------

/// The phylogenetic stream (nucleotide formats and/or binary NEXUS).
#[derive(Clone, Debug)]
pub(crate) struct PhyloPlan {
    pub(crate) formats: PhyloFormats,
    /// `--nexus-binary`.
    pub(crate) bin: bool,
    pub(crate) options: PhyloOptions,
    /// The effective sample gate: `--min-samples-locus` (default 4) clamped
    /// to the kept sample count (vcf2phylip's startup clamp, v13 §3.7).
    pub(crate) gate: u32,
}

/// Per-contig MSA regions, 0-based half-open.
pub(crate) type MsaRegions = Arc<HashMap<Vec<u8>, Vec<(u64, u64)>>>;

/// The MSA plan: the in-region predicate (0-based half-open per contig).
#[derive(Clone, Debug)]
pub(crate) struct MsaPlan {
    pub(crate) regions: MsaRegions,
}

impl MsaPlan {
    fn allows(&self, contig: &[u8], pos0: u64) -> bool {
        self.regions
            .get(contig)
            .is_some_and(|runs| runs.iter().any(|&(s, e)| pos0 >= s && pos0 < e))
    }
}

/// The row-target job: one parse, per-target folds (crate-internal; the
/// facade exports only the request and report types).
pub(crate) struct ConvertJob {
    pub(crate) n_samples: usize,
    pub(crate) keep: Arc<[bool]>,
    pub(crate) n_kept: usize,
    pub(crate) phylo: Option<PhyloPlan>,
    /// The `--012` sample gate (explicit `--min-samples-locus` only).
    pub(crate) count_gate: Option<u32>,
    /// The hapmap sample gate (explicit only).
    pub(crate) hmp_gate: Option<u32>,
    pub(crate) count: bool,
    pub(crate) hmp: bool,
    pub(crate) msa: Option<MsaPlan>,
    /// The composed statistics job for the SFS target (v13 §7.5).
    pub(crate) sfs: Option<StatsJob>,
    pub(crate) seed: u64,
}

/// The chunk identity key for the deterministic `--resolve-iupac` draw
/// (v13 §7.4.3): chunk boundaries are input-determined, so
/// `(chunk id, line ordinal)` is stable across thread counts; the seam
/// record of a chunk boundary keys off the chunk it belongs to (`id - 1`)
/// with the sentinel ordinal.
#[derive(Clone, Copy)]
pub(crate) struct DrawKey {
    pub(crate) chunk: u64,
    pub(crate) line: u32,
}

/// One sample's parsed genotype tokens (the single row-level view every
/// target consumes, D13-7).
#[derive(Clone, Copy, Default)]
pub(crate) struct GtSlots {
    /// Allele indices; `None` is the `.` missing allele. Up to four tokens
    /// are stored (more is malformed for reed).
    pub(crate) slots: [Option<u32>; 4],
    /// Token count (ploidy), capped at four.
    pub(crate) n: u8,
    /// More than four tokens.
    pub(crate) overflow: bool,
    /// A token that is neither `.` nor an unsigned integer.
    pub(crate) bad: bool,
    /// An allele token the *strict statistics scanner* rejects but convert
    /// itself accepts: more than five digits (`0000000`) or an index that
    /// does not fit the scanner's `u16`. Row targets parse it like the old
    /// code did (bcftools accepts zero-padded alleles), but the composed
    /// `--sfs` target must not hand the record to the strict scanner.
    pub(crate) strict_bad: bool,
    /// The raw GT text does not start with `.` (vcf2phylip's presence rule).
    pub(crate) present: bool,
    /// The raw GT field range within the samples region (for GEN_BIN).
    pub(crate) raw_start: u32,
    pub(crate) raw_len: u32,
}

/// Worker-local buffers.
#[derive(Default)]
pub struct ConvertScratch {
    gts: Vec<GtSlots>,
    /// Processed allele texts (flat) and their `(start, len)` ranges; range
    /// 0 is REF, range `i + 1` is `ALT[i]`.
    alleles: Vec<u8>,
    allele_ranges: Vec<(u32, u32)>,
    /// The reusable column / row buffer.
    column: Vec<u8>,
    /// The sanitized chunk body for the composed SFS job (mixed runs drop
    /// the structurally malformed lines the row targets skip, B13-2).
    sanitized: Vec<u8>,
    /// MSA per-record value buffers.
    msa_vals: Vec<u8>,
    msa_starts: Vec<u32>,
    msa_lens: Vec<u32>,
    msa_padded: Vec<u8>,
    /// Once-flags for the skip warnings.
    warned: WarnFlags,
    /// The composed statistics scratch.
    stats: StatsScratch,
}

#[derive(Default)]
struct WarnFlags {
    malformed: bool,
    multiallelic_012: bool,
    padded_msa: bool,
}

/// The per-chunk fold of every target.
#[derive(Default)]
pub(crate) struct ConvertPartial {
    pub(crate) sfs: Option<StatsPartial>,
    pub(crate) phylo: Option<PhyloPartial>,
    pub(crate) count: Option<CountPartial>,
    pub(crate) hmp: Option<HmpPartial>,
    pub(crate) msa: Option<MsaPartial>,
    pub(crate) counters: ConvertCounters,
}

/// The nucleotide + binary-NEXUS columns of one chunk.
#[derive(Default)]
pub(crate) struct PhyloPartial {
    pub(crate) n_sites: u64,
    pub(crate) n_bin_sites: u64,
    /// Per kept sample, this chunk's nucleotide letters.
    pub(crate) cols: Vec<Vec<u8>>,
    /// Per kept sample, the binary-NEXUS states.
    pub(crate) bin_cols: Vec<Vec<u8>>,
    /// Rendered `used_sites` rows (`CHROM\tPOS\tN\n` each).
    pub(crate) used: Vec<u8>,
}

/// The `--012` chunk: per kept sample values and the `.pos` lines.
#[derive(Default)]
pub(crate) struct CountPartial {
    pub(crate) n_sites: u64,
    pub(crate) cells: Vec<Vec<i8>>,
    pub(crate) pos: Vec<u8>,
}

/// The hapmap chunk: rendered site rows.
#[derive(Default)]
pub(crate) struct HmpPartial {
    pub(crate) n_sites: u64,
    pub(crate) rows: Vec<u8>,
}

/// One contig run of in-region MSA sites (positions ascending).
pub(crate) struct MsaRun {
    pub(crate) contig: Vec<u8>,
    pub(crate) sites: Vec<MsaSite>,
}

/// One MSA site: every kept sample's value, all `width` bytes long.
pub(crate) struct MsaSite {
    pub(crate) pos0: u64,
    pub(crate) vals: Vec<u8>,
    pub(crate) width: u32,
}

#[derive(Default)]
pub(crate) struct MsaPartial {
    pub(crate) runs: Vec<MsaRun>,
}

impl ConvertPartial {
    fn absorb(&mut self, other: &mut Self) {
        self.counters.absorb(&other.counters);
        if let Some(src) = &mut other.sfs {
            match &mut self.sfs {
                Some(dst) => dst.merge_sfs(src),
                None => self.sfs = Some(std::mem::take(src)),
            }
        }
    }
}

/// Ordered sink: merges every chunk partial in input order into one state.
pub(crate) struct ConvertSink {
    merged: ConvertPartial,
    /// Per-target blobs taken out of each partial before the merge and
    /// spilled past the threshold instead of held whole (doc_v15 D15-7).
    phylo_cols: Option<Vec<RowsBlob>>,
    phylo_bin_cols: Option<Vec<RowsBlob>>,
    phylo_used: RowsBlob,
    phylo_sites: u64,
    phylo_bin_sites: u64,
    count_cells: Option<Vec<RowsBlob>>,
    count_pos: RowsBlob,
    msa_runs: Option<Vec<(Vec<u8>, RowsBlob)>>,
    hmp_rows: RowsBlob,
    /// Per-chunk conversion scratch (cells are `i8`, the blobs hold bytes).
    cell_scratch: Vec<u8>,
    msa_scratch: Vec<u8>,
    out: Arc<Mutex<Option<MergedOutput>>>,
}

/// The sink's finished state: the merged partial plus the (possibly
/// spilled) per-target blobs.
pub(crate) struct MergedOutput {
    pub(crate) partial: ConvertPartial,
    pub(crate) phylo: Option<SpilledPhylo>,
    pub(crate) count: Option<SpilledCount>,
    pub(crate) msa: Option<SpilledMsa>,
    pub(crate) hmp_rows: RowsBlob,
}

/// The phylogenetic streams with their site tallies.
pub(crate) struct SpilledPhylo {
    pub(crate) n_sites: u64,
    pub(crate) n_bin_sites: u64,
    pub(crate) cols: Vec<RowsBlob>,
    pub(crate) bin_cols: Vec<RowsBlob>,
    pub(crate) used: RowsBlob,
}

/// The `--012` streams.
pub(crate) struct SpilledCount {
    pub(crate) cells: Vec<RowsBlob>,
    pub(crate) pos: RowsBlob,
}

impl OrderedSink<ConvertPartial> for ConvertSink {
    fn accept(
        &mut self,
        _id: reed_core::ChunkId,
        partial: &mut ConvertPartial,
    ) -> Result<(), EngineError> {
        if let Some(mut phylo) = partial.phylo.take() {
            let cols = self.phylo_cols.get_or_insert_with(|| {
                (0..phylo.cols.len())
                    .map(|_| RowsBlob::default())
                    .collect()
            });
            debug_assert_eq!(cols.len(), phylo.cols.len(), "invariant: kept-sample count");
            for (dst, src) in cols.iter_mut().zip(phylo.cols.iter_mut()) {
                dst.append(src, SPILL_THRESHOLD)?;
            }
            let bin_cols = self.phylo_bin_cols.get_or_insert_with(|| {
                (0..phylo.bin_cols.len())
                    .map(|_| RowsBlob::default())
                    .collect()
            });
            debug_assert_eq!(
                bin_cols.len(),
                phylo.bin_cols.len(),
                "invariant: kept-sample count"
            );
            for (dst, src) in bin_cols.iter_mut().zip(phylo.bin_cols.iter_mut()) {
                dst.append(src, SPILL_THRESHOLD)?;
            }
            self.phylo_used.append(&mut phylo.used, SPILL_THRESHOLD)?;
            self.phylo_sites += phylo.n_sites;
            self.phylo_bin_sites += phylo.n_bin_sites;
        }
        if let Some(mut count) = partial.count.take() {
            let cells = self.count_cells.get_or_insert_with(|| {
                (0..count.cells.len())
                    .map(|_| RowsBlob::default())
                    .collect()
            });
            debug_assert_eq!(cells.len(), count.cells.len(), "invariant: kept-sample count");
            for (dst, src) in cells.iter_mut().zip(count.cells.iter_mut()) {
                self.cell_scratch.clear();
                self.cell_scratch.extend(src.iter().map(|&cell| cell as u8));
                dst.append(&mut self.cell_scratch, SPILL_THRESHOLD)?;
            }
            self.count_pos.append(&mut count.pos, SPILL_THRESHOLD)?;
        }
        if let Some(msa) = partial.msa.take() {
            let runs = self.msa_runs.get_or_insert_with(Vec::new);
            for run in msa.runs {
                // The `runs` merge rule: the last run continues while the
                // contig does not change.
                let new_entry = !matches!(runs.last(), Some((name, _)) if *name == run.contig);
                if new_entry {
                    runs.push((run.contig, RowsBlob::default()));
                }
                let (_, blob) = runs.last_mut().expect("invariant: just ensured");
                self.msa_scratch.clear();
                for site in &run.sites {
                    self.msa_scratch
                        .extend_from_slice(&site.pos0.to_le_bytes());
                    self.msa_scratch
                        .extend_from_slice(&site.width.to_le_bytes());
                    self.msa_scratch.extend_from_slice(&site.vals);
                }
                blob.append(&mut self.msa_scratch, SPILL_THRESHOLD)?;
            }
        }
        if let Some(mut hmp) = partial.hmp.take() {
            self.hmp_rows.append(&mut hmp.rows, SPILL_THRESHOLD)?;
        }
        self.merged.absorb(partial);
        self.enforce_spill_budget()?;
        Ok(())
    }

    fn finish(self) -> Result<u64, EngineError> {
        *self
            .out
            .lock()
            .expect("invariant: convert slot not poisoned") = Some(MergedOutput {
            partial: self.merged,
            phylo: self.phylo_cols.map(|cols| SpilledPhylo {
                n_sites: self.phylo_sites,
                n_bin_sites: self.phylo_bin_sites,
                cols,
                bin_cols: self.phylo_bin_cols.unwrap_or_default(),
                used: self.phylo_used,
            }),
            count: self.count_cells.map(|cells| SpilledCount {
                cells,
                pos: self.count_pos,
            }),
            msa: self.msa_runs.map(|runs| SpilledMsa { runs }),
            hmp_rows: self.hmp_rows,
        });
        Ok(0)
    }
}

impl ConvertSink {
    /// Bounds the sink's buffered blob bytes: past [`SPILL_TOTAL_BUDGET`]
    /// the in-memory buffers are spilled whole in a fixed order (doc_v15
    /// D15-7; the bytes are unaffected, so the output is too).
    fn enforce_spill_budget(&mut self) -> Result<(), EngineError> {
        let mut total = self.buffered_bytes();
        if total <= SPILL_TOTAL_BUDGET {
            return Ok(());
        }
        let mut all: Vec<&mut RowsBlob> = Vec::new();
        if let Some(cols) = &mut self.phylo_cols {
            all.extend(cols.iter_mut());
        }
        if let Some(cols) = &mut self.phylo_bin_cols {
            all.extend(cols.iter_mut());
        }
        if let Some(cells) = &mut self.count_cells {
            all.extend(cells.iter_mut());
        }
        if let Some(runs) = &mut self.msa_runs {
            all.extend(runs.iter_mut().map(|(_, blob)| blob));
        }
        all.push(&mut self.phylo_used);
        all.push(&mut self.count_pos);
        all.push(&mut self.hmp_rows);
        for blob in all {
            if total <= SPILL_TOTAL_BUDGET {
                break;
            }
            let len = blob.mem_len();
            if len > 0 {
                blob.spill_now()?;
                total -= len;
            }
        }
        Ok(())
    }

    /// The sink's current in-memory blob bytes.
    fn buffered_bytes(&self) -> usize {
        let cols = |blobs: &Option<Vec<RowsBlob>>| -> usize {
            blobs.as_ref().map_or(0, |v| v.iter().map(RowsBlob::mem_len).sum())
        };
        cols(&self.phylo_cols)
            + cols(&self.phylo_bin_cols)
            + cols(&self.count_cells)
            + self.msa_runs.as_ref().map_or(0, |runs| {
                runs.iter().map(|(_, blob)| blob.mem_len()).sum()
            })
            + self.phylo_used.mem_len()
            + self.count_pos.mem_len()
            + self.hmp_rows.mem_len()
    }
}

impl ChunkJob for ConvertJob {
    type Partial = ConvertPartial;
    type Scratch = ConvertScratch;

    fn run_chunk(
        &self,
        chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let mut partial = ConvertPartial::default();
        // The composed SFS body is rebuilt only when a record must be
        // dropped: structural blanks need no dropping (the statistics job
        // skips them itself), but malformed lines (doc_v25 V25-03) and
        // records whose GT the strict scanner would reject — out-of-domain
        // alleles, non-digit tokens, polyploids (doc_v25 V25-01) — must not
        // reach it. Convert's policy for those is "count malformed, exit 0"
        // (B13-9), so handing them to the strict scanner used to turn the
        // same file into a fatal `--sfs` run. A clean chunk still hands the
        // original slice over (zero copy).
        let want_sfs = self.sfs.is_some();
        let mut dropped: Vec<usize> = Vec::new();
        let mut line_no: u32 = 0;
        for (index, line) in body.split(|&b| b == b'\n').enumerate() {
            let key = DrawKey {
                chunk: chunk.id.get(),
                line: line_no,
            };
            // Blank lines are not records: they do not consume a draw slot
            // (the key must stay aligned with the emitted rows).
            let outcome = self.push_line(line, key, &mut partial, scratch);
            if outcome != LineOutcome::Blank {
                line_no += 1;
            }
            if want_sfs {
                let keep = match outcome {
                    LineOutcome::Record => sfs_accepts_gt(scratch),
                    // `Blank` passes through: the statistics job skips it.
                    LineOutcome::Blank => true,
                    LineOutcome::Malformed => false,
                };
                if !keep {
                    dropped.push(index);
                }
            }
        }
        if let Some(stats) = &self.sfs {
            let ConvertScratch {
                sanitized,
                stats: stats_scratch,
                ..
            } = scratch;
            let sfs_body = if dropped.is_empty() {
                body
            } else {
                sanitized.clear();
                let mut next_drop = dropped.iter().copied().peekable();
                for (index, line) in body.split(|&b| b == b'\n').enumerate() {
                    if next_drop.peek() == Some(&index) {
                        next_drop.next();
                        continue;
                    }
                    if !sanitized.is_empty() {
                        sanitized.push(b'\n');
                    }
                    sanitized.extend_from_slice(line);
                }
                sanitized.as_slice()
            };
            partial.sfs = Some(stats.run_chunk(chunk, sfs_body, stats_scratch)?);
        }
        Ok(partial)
    }

    fn run_seam(
        &self,
        id: reed_core::ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let mut partial = ConvertPartial::default();
        // The seam record belongs to the previous chunk's stream.
        let key = DrawKey {
            chunk: id.get().saturating_sub(1),
            line: u32::MAX,
        };
        let outcome = self.push_line(line, key, &mut partial, scratch);
        // Feed the SFS target only what convert accepts (the same policy as
        // the body path: structural blanks/malformed and GT-domain rejects
        // never reach the strict statistics scanner).
        if let Some(stats) = &self.sfs
            && outcome == LineOutcome::Record
            && sfs_accepts_gt(scratch)
        {
            let ConvertScratch {
                stats: stats_scratch,
                ..
            } = scratch;
            partial.sfs = Some(stats.run_seam(id, line, stats_scratch)?);
        }
        Ok(partial)
    }
}

/// Whether the genotype just parsed by [`ConvertJob::parse_gts`] passes the
/// strict statistics scanner's grammar (doc_v25 V25-01): every sample's GT is
/// diploid-or-haploid, digit-only within the scanner's five-digit token cap,
/// and inside the record's allele domain.
fn sfs_accepts_gt(scratch: &ConvertScratch) -> bool {
    !scratch
        .gts
        .iter()
        .any(|gt| gt.bad || gt.strict_bad || gt.overflow || gt.n > 2)
}

impl ConvertJob {
    fn has_row_targets(&self) -> bool {
        self.phylo.is_some() || self.count || self.hmp || self.msa.is_some()
    }

    /// Folds one line into every row target. Returns whether the line was a
    /// record slot (non-blank): malformed lines still consume one, exactly as
    /// vcf2phylip counts them; blank lines do not exist for it.
    fn push_line(
        &self,
        line: &[u8],
        key: DrawKey,
        partial: &mut ConvertPartial,
        scratch: &mut ConvertScratch,
    ) -> LineOutcome {
        let lr = match classify_line(line, self.n_samples) {
            RecordClass::Record(lr) => lr,
            // `line_is_sane` used to call a blank seam line sane and hand it
            // to the statistics job; a blank line is now simply not a record
            // (doc_v25 V25-03).
            RecordClass::Blank => return LineOutcome::Blank,
            RecordClass::Malformed(why) => {
                // vcf2phylip counts every non-empty, non-`#` line as
                // processed, malformed ones included (`snp_num` increments
                // before its checks).
                partial.counters.total_records += 1;
                self.bump_malformed(partial, scratch, why.label());
                return LineOutcome::Malformed;
            }
        };
        partial.counters.total_records += 1;
        let pos1 = reed_io::vcf::fixed::parse_pos1(lr.pos_text(), lr.chrom())
            .expect("invariant: classify_line validated POS");
        let pos0 = pos1
            .checked_sub(1)
            .expect("invariant: classify_line rejects POS 0");

        // The single per-sample GT parse (D13-7).
        self.parse_gts(&lr, scratch);

        // Presence count over the kept samples (vcf2phylip counts samples
        // whose GT text does not start with `.`; v13 §3.1 step 2).
        let present = scratch
            .gts
            .iter()
            .zip(self.keep.iter())
            .filter(|(gt, kept)| **kept && gt.present)
            .count() as u32;

        // --- The phylogenetic stream (nucleotide + bin, vcf2phylip order) --
        // Its gates skip *this stream only*: a phylo early exit must not
        // starve the count/hapmap/MSA folds below (doc_v15 D15-4).
        'phylo: {
            let Some(plan) = &self.phylo else {
                break 'phylo;
            };
            if present < plan.gate {
                partial.counters.excess_missing += 1;
                break 'phylo;
            }
            self.build_alleles(&lr, true, scratch);
            if !is_snp(&lr, &scratch.allele_ranges) {
                partial.counters.mnps += 1;
                break 'phylo;
            }
            if plan.formats.any() {
                match self.phylo_column(scratch, key, plan.options.resolve_iupac) {
                    Ok(()) => {
                        let write_used = plan.options.write_used_sites;
                        let phy = Self::phylo_partial(partial, plan, self.n_kept);
                        for (i, byte) in scratch.column.iter().copied().enumerate() {
                            phy.cols[i].push(byte);
                        }
                        phy.n_sites += 1;
                        if write_used {
                            let used = &mut phy.used;
                            used.extend_from_slice(lr.chrom());
                            used.push(b'\t');
                            used.extend_from_slice(lr.pos_text());
                            used.push(b'\t');
                            used.extend_from_slice(present.to_string().as_bytes());
                            used.push(b'\n');
                        }
                        partial.counters.accepted += 1;
                    }
                    Err(Malformed) => {
                        // vcf2phylip's malformed skip also drops the record
                        // from the bin stream (the nucleotide check runs
                        // first and `continue`s).
                        self.bump_malformed(
                            partial,
                            scratch,
                            "a genotype with an unknown allele letter or an index beyond the record's alleles",
                        );
                        break 'phylo;
                    }
                }
            }
            // --- Binary NEXUS (nested under the SNP gate, v13 §3.4) -------
            if plan.bin && lr.alts().len() == 1 {
                let samples = lr.samples_text();
                let mut cells = Vec::with_capacity(self.n_kept);
                for (gt, &kept) in scratch.gts.iter().zip(self.keep.iter()) {
                    if !kept {
                        continue;
                    }
                    let raw = &samples
                        [gt.raw_start as usize..(gt.raw_start as usize + gt.raw_len as usize)];
                    cells.push(iupac::gen_bin_lookup(raw).unwrap_or(b'?'));
                }
                {
                    let phy = Self::phylo_partial(partial, plan, self.n_kept);
                    for (i, byte) in cells.into_iter().enumerate() {
                        phy.bin_cols[i].push(byte);
                    }
                    phy.n_bin_sites += 1;
                }
                partial.counters.biallelic_for_bin += 1;
            }
        }

        self.push_count_matrix(&lr, present, partial, scratch);
        self.push_hapmap(&lr, present, partial, scratch);
        self.push_msa(&lr, pos0, partial, scratch);
        LineOutcome::Record
    }

    /// The phylogenetic partial, created with the shape both branches agree
    /// on (nucleotide columns and/or binary columns, per the plan).
    fn phylo_partial<'a>(
        partial: &'a mut ConvertPartial,
        plan: &PhyloPlan,
        n_kept: usize,
    ) -> &'a mut PhyloPartial {
        partial.phylo.get_or_insert_with(|| PhyloPartial {
            n_sites: 0,
            n_bin_sites: 0,
            cols: if plan.formats.any() {
                vec![Vec::new(); n_kept]
            } else {
                Vec::new()
            },
            bin_cols: if plan.bin {
                vec![Vec::new(); n_kept]
            } else {
                Vec::new()
            },
            used: Vec::new(),
        })
    }

    fn bump_malformed(
        &self,
        partial: &mut ConvertPartial,
        scratch: &mut ConvertScratch,
        why: &str,
    ) {
        partial.counters.malformed += 1;
        if !scratch.warned.malformed {
            scratch.warned.malformed = true;
            tracing::warn!(
                "reed convert: skipping a malformed record ({why}); further ones are counted in PREFIX.log"
            );
        }
    }

    /// Parses every sample's GT tokens into `scratch.gts` (the D13-7 view).
    /// An allele index above the record's ALT count marks the genotype bad
    /// (doc_v25 V25-01): every target then counts the record as malformed —
    /// the `--012` path used to code a haploid out-of-domain call as `-1`
    /// silently.
    fn parse_gts(&self, lr: &LineRef<'_>, scratch: &mut ConvertScratch) {
        let samples = lr.samples_text();
        // The record's ALT count (0 for `ALT=.`): legal indices are 0..=n_alt.
        let n_alt: u32 = if lr.alts() == b"." {
            0
        } else {
            u32::try_from(memchr::memchr_iter(b',', lr.alts()).count() + 1).unwrap_or(u32::MAX)
        };
        // One GT-position lookup per record (doc_v25 V25-02): the old
        // `memchr(b':')` took the first subfield, so `FORMAT=DP:GT` handed
        // the DP value to the genotype parser.
        let gt_index = reed_io::vcf::format_layout::gt_index(lr.format_text());
        scratch.gts.clear();
        scratch.gts.resize(self.n_samples, GtSlots::default());
        let mut start = 0_usize;
        for gt in &mut scratch.gts {
            let begin = start.min(samples.len());
            let end =
                memchr::memchr(b'\t', &samples[begin..]).map_or(samples.len(), |rel| begin + rel);
            let column = &samples[begin..end];
            start = end + 1;
            // A column shorter than the GT position has no GT text; treat it
            // as missing, never as the first subfield.
            let raw = reed_io::vcf::format_layout::gt_subfield(column, gt_index)
                .unwrap_or(&b"."[..]);
            gt.raw_start = begin as u32;
            gt.raw_len = raw.len() as u32;
            gt.present = !raw.starts_with(b".");
            let mut tokens = 0_usize;
            for token in raw.split(|&b| b == b'/' || b == b'|') {
                if tokens >= 4 {
                    gt.overflow = true;
                    break;
                }
                if token == b"." {
                    gt.slots[tokens] = None;
                } else if !token.is_empty() && token.iter().all(u8::is_ascii_digit) {
                    // The row targets' grammar: any digit run (bcftools
                    // accepts zero-padded alleles), but the strict scanner
                    // caps a token at five digits — record that separately
                    // instead of failing this record for every target.
                    if token.len() > 5 {
                        gt.strict_bad = true;
                    }
                    match std::str::from_utf8(token)
                        .ok()
                        .and_then(|t| t.parse::<u32>().ok())
                    {
                        Some(idx) if idx <= n_alt => gt.slots[tokens] = Some(idx),
                        _ => gt.bad = true,
                    }
                } else {
                    gt.bad = true;
                }
                tokens += 1;
            }
            gt.n = u8::try_from(tokens.min(4)).unwrap_or(4);
        }
    }

    /// Builds the processed allele texts: range 0 = REF, `i + 1` = `ALT[i]`.
    /// `uppercase` selects the phylogenetic alphabet (uppercase, `-` → `*`,
    /// `<NON_REF>` → REF) versus MSA's literal alleles (`*` → `-`).
    fn build_alleles(&self, lr: &LineRef<'_>, uppercase: bool, scratch: &mut ConvertScratch) {
        scratch.alleles.clear();
        scratch.allele_ranges.clear();
        let push = |scratch: &mut ConvertScratch, bytes: &[u8]| {
            let start = scratch.alleles.len() as u32;
            scratch.alleles.extend_from_slice(bytes);
            scratch.allele_ranges.push((start, bytes.len() as u32));
        };
        if uppercase {
            let mut reference: Vec<u8> = lr
                .reference()
                .iter()
                .map(|&b| {
                    if b == b'-' {
                        b'*'
                    } else {
                        b.to_ascii_uppercase()
                    }
                })
                .collect();
            push(scratch, &reference);
            let ref_processed = std::mem::take(&mut reference);
            for alt in lr.alts().split(|&b| b == b',') {
                let mut processed: Vec<u8> = Vec::with_capacity(alt.len());
                let mut i = 0;
                while i < alt.len() {
                    if alt[i] == b'-' {
                        processed.push(b'*');
                        i += 1;
                    } else if alt[i..].starts_with(b"<NON_REF>") {
                        processed.extend_from_slice(&ref_processed);
                        i += b"<NON_REF>".len();
                    } else {
                        processed.push(alt[i].to_ascii_uppercase());
                        i += 1;
                    }
                }
                push(scratch, &processed);
            }
        } else {
            push(scratch, lr.reference());
            for alt in lr.alts().split(|&b| b == b',') {
                if alt == b"*" {
                    push(scratch, b"-");
                } else {
                    push(scratch, alt);
                }
            }
        }
    }

    /// vcf2phylip's SNP gate on the *processed* allele table: REF is one
    /// character and every ALT allele is one character.
    fn phylo_column(
        &self,
        scratch: &mut ConvertScratch,
        key: DrawKey,
        resolve: bool,
    ) -> Result<(), Malformed> {
        let ConvertScratch {
            gts,
            allele_ranges,
            alleles,
            column,
            ..
        } = scratch;
        column.clear();
        for (sample, gt) in gts.iter().enumerate() {
            if !self.keep[sample] {
                continue;
            }
            if gt.bad || gt.overflow {
                return Err(Malformed);
            }
            let letter_at = |index: u32| -> Result<&[u8], Malformed> {
                let &(start, len) = allele_ranges.get(index as usize).ok_or(Malformed)?;
                Ok(&alleles[start as usize..(start + len) as usize])
            };
            if resolve {
                if gt.n == 0 {
                    return Err(Malformed);
                }
                let n = usize::from(gt.n);
                let token = resolve_pick(key, sample as u32, n, self.seed);
                let letter = match gt.slots[token] {
                    None => b'N',
                    Some(index) => {
                        let allele = letter_at(index)?;
                        if allele.len() != 1 {
                            return Err(Malformed);
                        }
                        if allele[0] == b'*' { b'-' } else { allele[0] }
                    }
                };
                column.push(letter);
                continue;
            }
            let mut set = [0_u8; 8];
            let n = usize::from(gt.n).min(8);
            for (t, slot) in set.iter_mut().enumerate().take(n) {
                *slot = match gt.slots[t] {
                    None => b'N',
                    Some(index) => {
                        let allele = letter_at(index)?;
                        if allele.len() != 1 {
                            return Err(Malformed);
                        }
                        allele[0]
                    }
                };
            }
            set[..n].sort_unstable();
            let mut deduped = 0;
            for j in 0..n {
                if j == 0 || set[j] != set[j - 1] {
                    set[deduped] = set[j];
                    deduped += 1;
                }
            }
            match iupac::ambig_lookup(&set[..deduped]) {
                Some(letter) => column.push(letter),
                // The reference implementation raises here (its AMBIG lookup
                // is outside the malformed guard); reed skips + counts
                // (B13-9).
                None => return Err(Malformed),
            }
        }
        Ok(())
    }

    fn push_count_matrix(
        &self,
        lr: &LineRef<'_>,
        present: u32,
        partial: &mut ConvertPartial,
        scratch: &mut ConvertScratch,
    ) {
        if !self.count {
            return;
        }
        // vcftools' biallelic gate: `>2` alleles is skipped; `ALT='.'`
        // counts one allele and is kept.
        if raw_alt_count(lr) > 1 {
            if !scratch.warned.multiallelic_012 {
                scratch.warned.multiallelic_012 = true;
                tracing::warn!(
                    "reed convert --012: skipping multiallelic loci (vcftools behavior)"
                );
            }
            return;
        }
        if self.count_gate.is_some_and(|gate| present < gate) {
            return;
        }
        let mut cells = Vec::with_capacity(self.n_kept);
        for (gt, &kept) in scratch.gts.iter().zip(self.keep.iter()) {
            if !kept {
                continue;
            }
            let Some(cell) = cell_012(gt) else {
                self.bump_malformed(partial, scratch, "a genotype index beyond the ALT alleles");
                return;
            };
            cells.push(cell);
        }
        let out = partial.count.get_or_insert_with(|| CountPartial {
            n_sites: 0,
            cells: vec![Vec::new(); self.n_kept],
            pos: Vec::new(),
        });
        for (i, cell) in cells.into_iter().enumerate() {
            out.cells[i].push(cell);
        }
        out.n_sites += 1;
        out.pos.extend_from_slice(lr.chrom());
        out.pos.push(b'\t');
        out.pos.extend_from_slice(lr.pos_text());
        out.pos.push(b'\n');
    }

    fn push_hapmap(
        &self,
        lr: &LineRef<'_>,
        present: u32,
        partial: &mut ConvertPartial,
        scratch: &mut ConvertScratch,
    ) {
        if !self.hmp {
            return;
        }
        if self.hmp_gate.is_some_and(|gate| present < gate) {
            return;
        }
        // The SNP gate (shared with the phylogenetic stream, v13 §3.9).
        self.build_alleles(lr, true, scratch);
        if !is_snp(lr, &scratch.allele_ranges) {
            return;
        }
        match self.hapmap_row(lr, scratch) {
            Ok(row) => {
                let out = partial.hmp.get_or_insert_with(HmpPartial::default);
                out.rows.extend_from_slice(&row);
                out.n_sites += 1;
            }
            Err(Malformed) => {
                self.bump_malformed(partial, scratch, "a hapmap row without an IUPAC letter");
            }
        }
    }

    /// Renders one hapmap site row. Cells use TASSEL's diploid alphabet
    /// (v13 §3.9): the IUPAC letter when TASSEL reads it (A/C/G/T/R/Y/S/W/
    /// K/M/N and `-`), else `N` (its reader has no lowercase letters and no
    /// three-base codes V/H/D/B).
    fn hapmap_row(
        &self,
        lr: &LineRef<'_>,
        scratch: &mut ConvertScratch,
    ) -> Result<Vec<u8>, Malformed> {
        self.phylo_column(scratch, DrawKey { chunk: 0, line: 0 }, false)?;
        for letter in &mut scratch.column {
            *letter = match *letter {
                b'A' | b'C' | b'G' | b'T' | b'N' | b'R' | b'Y' | b'S' | b'W' | b'K' | b'M'
                | b'-' => *letter,
                _ => b'N',
            };
        }
        let mut row = Vec::with_capacity(64 + self.n_kept);
        // rs#: the VCF ID, or `chrom:pos` when absent.
        if lr.id() == b"." {
            row.extend_from_slice(lr.chrom());
            row.push(b':');
            row.extend_from_slice(lr.pos_text());
        } else {
            row.extend_from_slice(lr.id());
        }
        row.push(b'\t');
        // alleles: the declared REF/ALT set (B13-5).
        let mut first = true;
        let mut seen: Vec<&[u8]> = Vec::new();
        for (i, &(start, len)) in scratch.allele_ranges.iter().enumerate() {
            let mut text = &scratch.alleles[start as usize..(start + len) as usize];
            // Skip an ALT `.` (no alternative allele).
            if i > 0 && text == b"." {
                continue;
            }
            // The GATK deletion symbol is TASSEL's gap.
            if text == b"*" {
                text = b"-";
            }
            if seen.contains(&text) {
                continue;
            }
            seen.push(text);
            if !first {
                row.push(b'/');
            }
            row.extend_from_slice(text);
            first = false;
        }
        row.push(b'\t');
        row.extend_from_slice(lr.chrom());
        row.push(b'\t');
        row.extend_from_slice(lr.pos_text());
        row.extend_from_slice(b"\t+\tNA\tNA\tNA\tNA\tNA\tNA\t");
        for (i, &letter) in scratch.column.iter().enumerate() {
            row.push(letter);
            if i + 1 < self.n_kept {
                row.push(b'\t');
            }
        }
        row.push(b'\n');
        Ok(row)
    }

    fn push_msa(
        &self,
        lr: &LineRef<'_>,
        pos0: u64,
        partial: &mut ConvertPartial,
        scratch: &mut ConvertScratch,
    ) {
        let Some(plan) = &self.msa else { return };
        if !plan.allows(lr.chrom(), pos0) {
            return;
        }
        // Literal allele texts (vcf2msa's genotype_resolve domain).
        self.build_alleles(lr, false, scratch);
        match self.msa_values(scratch) {
            Ok((vals, width)) => {
                let out = partial.msa.get_or_insert_with(MsaPartial::default);
                let site = MsaSite { pos0, vals, width };
                match out.runs.last_mut() {
                    Some(run) if run.contig == lr.chrom() => run.sites.push(site),
                    _ => out.runs.push(MsaRun {
                        contig: lr.chrom().to_vec(),
                        sites: vec![site],
                    }),
                }
            }
            Err(Malformed) => {
                self.bump_malformed(partial, scratch, "an MSA genotype with an unknown allele");
            }
        }
    }

    /// Resolves every kept sample's MSA value for one record (vcf2msa's
    /// `genotype_resolve` minus the external aligner: values are right-padded
    /// with `-` to the site's maximum length, B13-6).
    fn msa_values(&self, scratch: &mut ConvertScratch) -> Result<(Vec<u8>, u32), Malformed> {
        let ConvertScratch {
            gts,
            allele_ranges,
            alleles,
            msa_vals,
            msa_starts,
            msa_lens,
            msa_padded,
            warned,
            ..
        } = scratch;
        msa_vals.clear();
        msa_starts.clear();
        msa_lens.clear();
        for (sample, gt) in gts.iter().enumerate() {
            if !self.keep[sample] {
                continue;
            }
            if gt.bad || gt.overflow {
                return Err(Malformed);
            }
            msa_starts.push(msa_vals.len() as u32);
            resolve_msa_value(gt, allele_ranges, alleles, msa_vals)?;
            let start = *msa_starts.last().expect("invariant: just pushed");
            msa_lens.push(msa_vals.len() as u32 - start);
        }
        let width = msa_lens.iter().copied().max().unwrap_or(1).max(1);
        if width == 1 {
            return Ok((std::mem::take(msa_vals), 1));
        }
        if !warned.padded_msa {
            warned.padded_msa = true;
            tracing::warn!(
                "reed convert --msa: indel sites are gap-padded deterministically \
                 (no external aligner; B13-6)"
            );
        }
        msa_padded.clear();
        for (i, &len) in msa_lens.iter().enumerate() {
            let start = msa_starts[i] as usize;
            msa_padded.extend_from_slice(&msa_vals[start..start + len as usize]);
            for _ in len..width {
                msa_padded.push(b'-');
            }
        }
        Ok((std::mem::take(msa_padded), width))
    }
}

/// What one input line was for the row targets (`push_line`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum LineOutcome {
    /// Only line terminators: not a record, consumes no draw slot.
    Blank,
    /// A non-blank line that fails the structural checks: counted as
    /// malformed (vcf2phylip's `snp_num` slot) and skipped.
    Malformed,
    /// A structurally valid record; GT-domain rejection is decided per
    /// target from the parsed `gts`.
    Record,
}

/// A malformed record for one target (skipped with a counter).
struct Malformed;

/// The raw ALT allele count (`.` counts zero; else commas + 1).
fn raw_alt_count(lr: &LineRef<'_>) -> u32 {
    let alts = lr.alts();
    if alts == b"." {
        return 0;
    }
    memchr::memchr_iter(b',', alts).count() as u32 + 1
}

/// vcf2phylip's SNP gate on the processed allele table: REF (raw) is one
/// character and every ALT allele is one character.
fn is_snp(lr: &LineRef<'_>, ranges: &[(u32, u32)]) -> bool {
    if lr.reference().len() != 1 {
        return false;
    }
    ranges[1..].iter().all(|&(_, len)| len == 1)
}

/// The `--012` cell for one sample (v13 §3.8, B13-8): `0`/`1`/`2` for a full
/// diploid call, `-1` for missing, half calls and haploids.
fn cell_012(gt: &GtSlots) -> Option<i8> {
    if gt.bad || gt.overflow || gt.n == 0 {
        return None;
    }
    if gt.n != 2 {
        // Haploid (and any other ploidy): reed emits -1; vcftools emits
        // 1/0/1 (B13-8).
        return Some(-1);
    }
    match (gt.slots[0], gt.slots[1]) {
        (Some(a), Some(b)) => {
            if a > 1 || b > 1 {
                return None;
            }
            Some(match (a, b) {
                (0, 0) => 0,
                (1, 1) => 2,
                _ => 1,
            })
        }
        // Missing or half call.
        _ => Some(-1),
    }
}

/// Resolves one sample's MSA value (vcf2msa `genotype_resolve`): single
/// characters join into an IUPAC letter; indel alleles use the literal text
/// (the shortest on conflict, ties lexicographic); a fully missing genotype
/// yields the REF allele (vcf2msa's step 4).
fn resolve_msa_value(
    gt: &GtSlots,
    ranges: &[(u32, u32)],
    flat: &[u8],
    out: &mut Vec<u8>,
) -> Result<(), Malformed> {
    let allele = |index: u32| -> Option<&[u8]> {
        let &(start, len) = ranges.get(index as usize)?;
        Some(&flat[start as usize..(start + len) as usize])
    };
    let mut called = 0;
    let mut var: Vec<u8> = Vec::new();
    let mut indel: Vec<Vec<u8>> = Vec::new();
    for t in 0..usize::from(gt.n) {
        let Some(index) = gt.slots[t] else { continue };
        let text = allele(index).ok_or(Malformed)?;
        called += 1;
        if text.len() == 1 {
            if !var.contains(&text[0]) {
                var.push(text[0]);
            }
        } else {
            let value = text.to_vec();
            if !indel.contains(&value) {
                indel.push(value);
            }
        }
    }
    if called == 0 {
        // No called allele: the REF string (vcf2msa step 4).
        out.extend_from_slice(allele(0).ok_or(Malformed)?);
        return Ok(());
    }
    if !var.is_empty() {
        if var.len() == 1 {
            out.push(var[0]);
            return Ok(());
        }
        var.sort_unstable();
        var.dedup();
        let letter = iupac::ambig_lookup(&var).ok_or(Malformed)?;
        out.push(letter);
        return Ok(());
    }
    // Indel-only: the shortest literal, ties lexicographic.
    indel.sort();
    indel.dedup();
    let shortest = indel.iter().min_by_key(|v| v.len()).ok_or(Malformed)?;
    out.extend_from_slice(shortest);
    Ok(())
}

/// The deterministic `--resolve-iupac` draw (splitmix64 over
/// `(chunk, line, sample, seed)`, v13 §7.4.3).
fn resolve_pick(key: DrawKey, sample: u32, n: usize, seed: u64) -> usize {
    let mut x = seed
        ^ key.chunk.wrapping_mul(0x9E37_79B9_7F4A_7C15)
        ^ u64::from(key.line).wrapping_mul(0xC2B2_AE3D_27D4_EB4F)
        ^ u64::from(sample).wrapping_mul(0x1656_67B1_9E37_79F9);
    x ^= x >> 30;
    x = x.wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^= x >> 31;
    (x % n as u64) as usize
}

// ---------------------------------------------------------------------------
// Driver
// ---------------------------------------------------------------------------

/// One MSA output unit (a region or a whole reference contig).
pub(crate) struct MsaOutput {
    pub(crate) contig: Vec<u8>,
    pub(crate) start1: u64,
    pub(crate) end1: u64,
    /// The file name carries the coordinates only for explicit regions.
    pub(crate) explicit: bool,
}

/// Reads a FASTA into per-contig sequences (the MSA reference).
///
/// Streams the file line by line: `fs::read` held the whole reference (a
/// human genome is ~3 GB) and then copied every contig into an `Arc<[u8]>`,
/// so the peak was twice the file plus all the sequence data (doc_v23
/// V23-11). Only the current contig's sequence is buffered twice (the
/// `Vec` and the `Arc` it is moved into), and the `Arc` is the only copy
/// that outlives a record.
fn read_reference(path: &PathBuf) -> Result<HashMap<Vec<u8>, Arc<[u8]>>, crate::Error> {
    use std::io::BufRead as _;
    let file = std::fs::File::open(path).map_err(|e| {
        crate::Error::usage(format!("cannot read --reference {}: {e}", path.display()))
    })?;
    let reader = std::io::BufReader::with_capacity(1 << 20, file);
    let mut map: HashMap<Vec<u8>, Arc<[u8]>> = HashMap::new();
    let mut name: Option<Vec<u8>> = None;
    let mut seq: Vec<u8> = Vec::new();
    for line in reader.split(b'\n') {
        let line = line.map_err(|e| {
            crate::Error::usage(format!("cannot read --reference {}: {e}", path.display()))
        })?;
        let line = line.strip_suffix(b"\r").unwrap_or(&line);
        if let Some(header) = line.strip_prefix(b">") {
            if let Some(n) = name.take() {
                map.insert(n, Arc::from(std::mem::take(&mut seq)));
            }
            let token = header
                .split(|&b| b == b' ' || b == b'\t')
                .next()
                .unwrap_or(header);
            if token.is_empty() {
                return Err(crate::Error::usage("--reference: empty FASTA record name"));
            }
            name = Some(token.to_vec());
        } else if !line.is_empty() {
            if name.is_none() {
                return Err(crate::Error::usage(
                    "--reference: FASTA sequence before the first `>` header",
                ));
            }
            seq.extend_from_slice(line);
        }
    }
    if let Some(n) = name.take() {
        map.insert(n, Arc::from(seq));
    }
    if map.is_empty() {
        return Err(crate::Error::usage("--reference: no FASTA records"));
    }
    Ok(map)
}

/// Builds the MSA plan and output list from the reference and the request
/// regions (whole contigs when none were given).
fn build_msa(
    reference: &HashMap<Vec<u8>, Arc<[u8]>>,
    regions: &[MsaRegion],
    flank: u64,
) -> Result<(MsaPlan, Vec<MsaOutput>), crate::Error> {
    let mut predicates: HashMap<Vec<u8>, Vec<(u64, u64)>> = HashMap::new();
    let mut outputs: Vec<MsaOutput> = Vec::new();
    if regions.is_empty() {
        // Deterministic file order: the reference sorted by name.
        let mut names: Vec<&Vec<u8>> = reference.keys().collect();
        names.sort();
        for name in names {
            let len = reference[name].len() as u64;
            predicates.insert(name.clone(), vec![(0, len)]);
            outputs.push(MsaOutput {
                contig: name.clone(),
                start1: 1,
                end1: len,
                explicit: false,
            });
        }
    } else {
        for region in regions {
            let contig = region.contig.as_bytes().to_vec();
            let Some(sequence) = reference.get(&contig) else {
                return Err(crate::Error::usage(format!(
                    "--region names contig {} which the reference FASTA does not contain",
                    region.contig
                )));
            };
            let len = sequence.len() as u64;
            let start1 = region.start1.unwrap_or(1).max(1);
            let end1 = region.end1.unwrap_or(len).min(len);
            if start1 > len || end1 < start1 {
                return Err(crate::Error::usage(format!(
                    "--region {}:{start1} exceeds the reference ({len} bp)",
                    region.contig
                )));
            }
            let start0 = (start1 - 1).saturating_sub(flank);
            let end0 = (end1 + flank).min(len);
            predicates
                .entry(contig.clone())
                .or_default()
                .push((start0, end0));
            outputs.push(MsaOutput {
                contig,
                start1: start0 + 1,
                end1: end0,
                explicit: true,
            });
        }
    }
    Ok((
        MsaPlan {
            regions: Arc::new(predicates),
        },
        outputs,
    ))
}

/// Runs the conversion over a path (or stdin when `None`).
///
/// # Errors
///
/// Request validation (Usage) and transport/render failures (Runtime).
pub fn run_convert(
    path: Option<&std::path::Path>,
    request: &ConvertRequest,
    budget: &ThreadBudget,
) -> Result<ConvertReport, crate::Error> {
    // The hapmap→VCF target reads a hapmap input, which the VCF opener would
    // reject; it is exclusive with every VCF-input target (doc_v13 §1.1
    // rule 6, §7.7).
    if request
        .targets
        .iter()
        .any(|t| matches!(t, ConvertTarget::Vcf { .. }))
    {
        let [
            ConvertTarget::Vcf {
                compressed,
                level,
                output,
                command_line,
                no_version,
            },
        ] = request.targets.as_slice()
        else {
            return Err(crate::Error::usage(
                "cannot mix --vcf (a hapmap input) with VCF-input targets",
            ));
        };
        let out = crate::convert::hapmap::VcfOutput {
            path: output.clone(),
            compressed: *compressed,
            level: *level,
            command_line: command_line.clone(),
            no_version: *no_version,
        };
        return crate::convert::hapmap::run_hapmap_to_vcf(path, &out, budget);
    }
    let opened = crate::input::open(path).map_err(EngineError::from)?;
    run_convert_opened(opened, request, budget)
}

/// Runs the conversion over an already-opened input (the CLI opens once so
/// stdin streams).
///
/// # Errors
///
/// As [`run_convert`].
pub fn run_convert_opened(
    opened: crate::input::OpenedInput,
    request: &ConvertRequest,
    budget: &ThreadBudget,
) -> Result<ConvertReport, crate::Error> {
    let sample_names = opened.facts.sample_names.clone();
    let n_samples = sample_names.len();
    if n_samples == 0 {
        return Err(crate::Error::usage(
            "the header declares no samples; nothing to convert",
        ));
    }
    if request.targets.is_empty() {
        return Err(crate::Error::usage(
            "nothing to convert: give --sfs, a phylo format (--phylip/--fasta/--nexus), \
             --012, --hapmap or --msa",
        ));
    }

    // The kept-sample view (deterministic header order).
    let keep = reed_core::sample::KeepMask::try_new(
        n_samples,
        request
            .keep_samples
            .clone()
            .unwrap_or_else(|| vec![true; n_samples]),
    )
    .map_err(|e| crate::Error::usage(e.to_string()))?;
    if !keep.iter().any(|&k| k) {
        return Err(crate::Error::usage("the sample subset selects no samples"));
    }
    let n_kept = keep.n_selected();
    let keep: Arc<[bool]> = Arc::from(keep.as_slice());

    // Resolve the targets.
    let mut phylo: Option<PhyloPlan> = None;
    let mut count = false;
    let mut hmp = false;
    let mut msa: Option<MsaPlan> = None;
    let mut msa_outputs: Vec<MsaOutput> = Vec::new();
    let mut reference: Option<HashMap<Vec<u8>, Arc<[u8]>>> = None;
    let mut sfs_request: Option<(crate::stats::request::SfsRequest, bool)> = None;
    let mut nexus_bin = false;
    let mut msa_collapse: Option<String> = None;
    let explicit_gate = request.min_samples_locus;
    for target in &request.targets {
        match target {
            ConvertTarget::Sfs {
                request: sfs,
                int_dtype,
            } => sfs_request = Some((sfs.clone(), *int_dtype)),
            ConvertTarget::Phylo { formats, options } => {
                // A formats-empty target is allowed when it carries an
                // outgroup: the bin-only run reads it from this plan
                // (doc_v15 D15-6).
                if !formats.any() && options.outgroup.is_none() {
                    return Err(crate::Error::usage(
                        "the phylo target selects no format (--phylip/--fasta/--nexus)",
                    ));
                }
                if phylo.is_some() {
                    return Err(crate::Error::usage("repeated phylo target"));
                }
                let gate = explicit_gate.unwrap_or(4).min(n_kept as u32);
                phylo = Some(PhyloPlan {
                    formats: *formats,
                    bin: false,
                    options: options.clone(),
                    gate,
                });
            }
            ConvertTarget::NexusBin => nexus_bin = true,
            ConvertTarget::CountMatrix => count = true,
            ConvertTarget::Hapmap => hmp = true,
            ConvertTarget::Vcf { .. } => {
                // Guarded in `run_convert`; `run_convert_opened` cannot open a
                // hapmap (its `OpenedInput` is VCF-shaped).
                return Err(crate::Error::usage(
                    "the --vcf target reads a hapmap file; it cannot run from an already-opened VCF input",
                ));
            }
            ConvertTarget::Msa {
                reference: path,
                regions,
                flank,
                collapse_suffix,
            } => {
                let map = read_reference(path)?;
                let (plan, outputs) = build_msa(&map, regions, *flank)?;
                msa = Some(plan);
                msa_outputs = outputs;
                reference = Some(map);
                msa_collapse.clone_from(collapse_suffix);
            }
        }
    }
    // A bin-only run still carries the phylo record counters and gate.
    if nexus_bin {
        match &mut phylo {
            Some(plan) => plan.bin = true,
            None => {
                phylo = Some(PhyloPlan {
                    formats: PhyloFormats::default(),
                    bin: true,
                    options: PhyloOptions {
                        outgroup: None,
                        resolve_iupac: false,
                        write_used_sites: false,
                    },
                    gate: explicit_gate.unwrap_or(4).min(n_kept as u32),
                });
            }
        }
    }
    // The outgroup must name a kept sample (B13-4: usage error).
    if let Some(plan) = &phylo
        && let Some(outgroup) = &plan.options.outgroup
    {
        {
            let found = keep
                .iter()
                .enumerate()
                .any(|(i, &k)| k && sample_names[i] == *outgroup);
            if !found {
                return Err(crate::Error::usage(format!(
                    "--outgroup {outgroup:?} is not among the samples (after -s/-S/--keep/--remove)"
                )));
            }
        }
    }

    // The output sample names (MSA may collapse theirs at a suffix; two
    // samples collapsing onto one name is refused instead of silently
    // merging, B13-7).
    let kept_names: Vec<String> = keep
        .iter()
        .enumerate()
        .filter(|&(_, &k)| k)
        .map(|(i, _)| sample_names[i].clone())
        .collect();
    let msa_names: Vec<String> = match &msa_collapse {
        None => kept_names.clone(),
        Some(suffix) => kept_names
            .iter()
            .map(|n| n.split(suffix.as_str()).next().unwrap_or(n).to_string())
            .collect(),
    };
    if msa_collapse.is_some() {
        let mut seen: Vec<&str> = Vec::new();
        for name in &msa_names {
            if seen.contains(&name.as_str()) {
                return Err(crate::Error::usage(format!(
                    "--collapse-suffix maps two samples onto {name:?};                      reed refuses to merge them (vcf2msa would)"
                )));
            }
            seen.push(name);
        }
    }

    // The SFS plan (the same validation as stats, D13-6).
    let (gt_groups, group_names) = resolve_groups(
        &request.populations,
        request.force_samples,
        &sample_names,
        request.keep_samples.as_deref(),
    )?;
    let n_groups = group_names.len();
    let mut group_sizes = vec![0_u32; n_groups];
    for group in gt_groups.iter().flatten() {
        group_sizes[*group as usize] += 1;
    }
    let sfs_plan = sfs_request
        .as_ref()
        .map(|(sfs, _)| build_sfs_plan(sfs, &group_names, &group_sizes))
        .transpose()?;

    let abort = reed_engine::abort::abort_for(budget);
    let contig_len: Arc<HashMap<Vec<u8>, u64>> = Arc::new(
        opened
            .facts
            .contig_names
            .iter()
            .zip(opened.facts.contig_lengths.iter().copied())
            .filter_map(|(name, len)| len.map(|len| (name.as_bytes().to_vec(), len)))
            .collect(),
    );
    let job = ConvertJob {
        n_samples,
        keep: Arc::clone(&keep),
        n_kept,
        phylo,
        count,
        count_gate: explicit_gate.map(|g| g.min(n_kept as u32)),
        hmp,
        hmp_gate: explicit_gate.map(|g| g.min(n_kept as u32)),
        msa,
        sfs: sfs_plan.as_ref().map(|plan| {
            StatsJob::for_sfs(
                n_samples,
                gt_groups.clone(),
                plan.clone(),
                Arc::clone(&contig_len),
            )
        }),
        seed: request.seed,
    };
    // The render decisions are captured before the job moves into the
    // input wrapper.
    let has_row_targets = job.has_row_targets();
    let render_phylo = job.phylo.clone();
    let render_count = job.count;
    let render_hmp = job.hmp;
    let out = Arc::new(Mutex::new(None));
    let sink = ConvertSink {
        merged: ConvertPartial::default(),
        phylo_cols: None,
        phylo_bin_cols: None,
        phylo_used: RowsBlob::default(),
        phylo_sites: 0,
        phylo_bin_sites: 0,
        count_cells: None,
        count_pos: RowsBlob::default(),
        msa_runs: None,
        hmp_rows: RowsBlob::default(),
        cell_scratch: Vec::new(),
        msa_scratch: Vec::new(),
        out: Arc::clone(&out),
    };
    let report = {
        // The text path (the InputJob wrapper renders BCF records to text):
        // the native record path is a statistics-only optimization.
        let job = crate::bcf_input::InputJob::new(job, &opened, &abort, budget);
        let mut source = opened.source;
        let report = run(source.as_mut(), &job, sink, budget, &abort)?;
        job.finish(report.n_chunks, true)?;
        report
    };
    let MergedOutput {
        partial: merged,
        phylo: spilled_phylo,
        count: spilled_count,
        msa: spilled_msa,
        hmp_rows,
    } = out
        .lock()
        .expect("invariant: convert slot not poisoned")
        .take()
        .expect("invariant: the sink stored the merged state");

    // ---- Render and write (drain) ---------------------------------------
    let mut files: Vec<PathBuf> = Vec::new();
    let prefix = request.prefix.as_path();
    let mut counters = None;
    if has_row_targets {
        counters = Some(merged.counters);
    }
    if let Some(plan) = &render_phylo {
        files.extend(writers::write_phylo(
            prefix,
            plan,
            &kept_names,
            spilled_phylo.as_ref(),
        )?);
    }
    if render_count {
        files.extend(writers::write_count_matrix(
            prefix,
            &kept_names,
            &keep,
            spilled_count.as_ref(),
        )?);
    }
    if render_hmp {
        files.extend(writers::write_hapmap(prefix, &kept_names, Some(&hmp_rows))?);
    }
    if let Some(reference) = &reference {
        files.extend(writers::write_msa(
            prefix,
            reference,
            &msa_outputs,
            &msa_names,
            spilled_msa.as_ref(),
        )?);
    }
    let mut sfs_report_out = None;
    let mut n_records = merged.counters.total_records;
    if let (Some((_, int_dtype)), Some(plan), Some(partial)) =
        (&sfs_request, &sfs_plan, &merged.sfs)
    {
        // `SfsPartial` is seeded lazily by the first site that reaches the
        // projector, so an input where every record is skipped (all
        // multiallelic, ALT `.` with a GT call, out-of-domain GT, polyploid)
        // leaves the accumulator unseeded. That is a legitimate empty
        // spectrum, not a broken invariant: write the zero spectrum like
        // `stats` does for a 0-site request instead of panicking (the
        // `--sfs` panic is what made v25's stricter GT rejection abort runs
        // the old build completed).
        let empty;
        let partial_sfs = match partial.sfs_partial() {
            Some(sfs) => sfs,
            None => {
                empty = SfsPartial::new(plan)
                    .ok_or_else(|| crate::Error::usage("invalid SFS projections"))?;
                &empty
            }
        };
        let report = sfs_report(partial_sfs, plan);
        files.extend(crate::stats::writers::sfs::write_projected_sfs(
            prefix,
            &group_names,
            &report,
            *int_dtype,
        )?);
        sfs_report_out = Some(report);
        n_records = partial.records();
    }
    Ok(ConvertReport {
        files,
        counters,
        sfs: sfs_report_out,
        n_records,
        elapsed: report.elapsed,
        metrics: report.metrics,
    })
}

#[cfg(test)]
mod spill_tests {
    use super::*;
    use std::io::Write as _;

    fn spill_path_of(blob: &RowsBlob) -> PathBuf {
        match blob {
            RowsBlob::Spilled { spill, .. } => spill.path().to_path_buf(),
            RowsBlob::Memory(_) => panic!("expected a spilled blob"),
        }
    }

    #[test]
    fn rows_blob_spills_and_cleans_up() {
        let mut blob = RowsBlob::default();
        let mut first = b"aaaa".to_vec();
        blob.append(&mut first, 8).expect("append");
        let mut second = b"bbbb".to_vec();
        blob.append(&mut second, 8).expect("append"); // 8 bytes: still memory
        assert!(matches!(blob, RowsBlob::Memory(_)));
        let mut third = b"cccc".to_vec();
        blob.append(&mut third, 8).expect("append"); // 12 > 8: spills
        let path = spill_path_of(&blob);
        assert!(path.exists(), "the spill file must exist");
        let mut fourth = b"dddd".to_vec();
        blob.append(&mut fourth, 8).expect("append"); // appended to the file
        let mut out = Vec::new();
        blob.write_to(&mut out).expect("write_to");
        assert_eq!(out, b"aaaabbbbccccdddd");
        drop(blob);
        assert!(!path.exists(), "the spill file must be removed on drop");
    }

    #[test]
    fn write_hapmap_streams_a_spilled_blob() {
        let dir = std::env::temp_dir().join(format!("reed-convert-spill-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("invariant: temp dir");
        let row: &[u8] = b"chr1:1\tA/G\tchr1\t1\t+\tNA\tNA\tNA\tNA\tNA\tNA\tA\n";
        // A real TempSpill, so the drop path is the production one.
        let (mut file, spill) =
            crate::spill::TempSpill::create(Some(&dir), "convert", "spill")
                .expect("invariant: spill file");
        file.write_all(row).expect("invariant: spill write");
        let spill_path = spill.path().to_path_buf();
        let rows = RowsBlob::Spilled {
            spill,
            len: row.len() as u64,
        };
        let names = vec!["s0".to_string()];
        let files = writers::write_hapmap(&dir.join("out"), &names, Some(&rows)).expect("write");
        let text = std::fs::read_to_string(&files[0]).expect("invariant: read");
        let expected = "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\ts0\nchr1:1\tA/G\tchr1\t1\t+\tNA\tNA\tNA\tNA\tNA\tNA\tA\n";
        assert_eq!(text, expected);
        drop(rows);
        assert!(
            !spill_path.exists(),
            "the spill file must be removed on drop"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn rows_blob_visit_sees_the_whole_stream() {
        let mut blob = RowsBlob::default();
        let mut head = vec![1_u8, 2, 255];
        blob.append(&mut head, 2).expect("append"); // 3 > 2: spills at once
        let mut tail = vec![0_u8; 200_000];
        blob.append(&mut tail, 2).expect("append");
        let mut seen = Vec::new();
        blob.visit(|chunk| {
            seen.extend_from_slice(chunk);
            Ok(())
        })
        .expect("visit");
        assert_eq!(seen.len(), 200_003, "visit must see every byte");
        assert_eq!(&seen[..3], &[1, 2, 255]);
    }

    #[test]
    fn write_count_matrix_streams_spilled_cells() {
        let dir = std::env::temp_dir().join(format!("reed-convert-cells-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("invariant: temp dir");
        let raw: Vec<u8> = vec![1, 2, 0xFF, 0]; // 0/1, 1/1, ./., 0/0
        let (mut cell_file, cells_spill) =
            crate::spill::TempSpill::create(Some(&dir), "convert", "spill")
                .expect("invariant: spill file");
        cell_file.write_all(&raw).expect("invariant: spill write");
        let count = SpilledCount {
            cells: vec![RowsBlob::Spilled {
                spill: cells_spill,
                len: raw.len() as u64,
            }],
            pos: RowsBlob::Memory(b"chr1\t1\nchr1\t2\nchr1\t3\nchr1\t4\n".to_vec()),
        };
        let names = vec!["s0".to_string()];
        writers::write_count_matrix(&dir.join("out"), &names, &[true], Some(&count))
            .expect("write");
        assert_eq!(
            std::fs::read_to_string(dir.join("out.012")).expect("invariant: read"),
            "0\t1\t2\t-1\t0\n"
        );
        assert_eq!(
            std::fs::read_to_string(dir.join("out.012.pos")).expect("invariant: read"),
            "chr1\t1\nchr1\t2\nchr1\t3\nchr1\t4\n"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }
}
