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

use crate::context::RunContext;
use reed_engine::error::EngineError;
use reed_engine::metrics::StageMetrics;
use reed_engine::pipeline::{ChunkJob, OrderedSink, run};
use reed_io::vcf::line::{LineRef, RecordClass, classify_line};

use crate::convert::writers;
use crate::convert::{ConvertRequest, ConvertTarget, PhyloFormats, PhyloOptions};
use crate::stats::driver::{build_sfs_plan, resolve_groups, sfs_report};
use crate::stats::request::{SfsPartial, SfsReport};
use crate::stats::run::{StatsJob, StatsPartial, StatsScratch};

mod convert_job;
mod msa;

pub(crate) use msa::MsaOutput;

pub(crate) use super::spill::{RowsBlob, SPILL_THRESHOLD, SPILL_TOTAL_BUDGET, SpilledMsa};
use msa::{build_msa, read_reference};

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
    /// Records with a polyploid kept genotype that entered an SFS-only or
    /// `--012` target. Polyploidy is outside the supported domain (diploid
    /// plus haploid-as-half-call): `--sfs` skips the record and `--012`
    /// writes its documented `-1`, both with one warning and this count
    /// (v29 D29-06).
    pub polyploid: u64,
}

impl ConvertCounters {
    fn absorb(&mut self, other: &Self) {
        self.total_records += other.total_records;
        self.excess_missing += other.excess_missing;
        self.mnps += other.mnps;
        self.accepted += other.accepted;
        self.biallelic_for_bin += other.biallelic_for_bin;
        self.malformed += other.malformed;
        self.polyploid += other.polyploid;
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
    /// Records skipped because a kept genotype is polyploid, regardless of
    /// which targets ran (also in [`ConvertCounters::polyploid`] when a row
    /// target ran). v29 D29-06.
    pub polyploid: u64,
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

/// The reference FASTA: contig name → sequence bytes.
pub(crate) type ReferenceMap = HashMap<Vec<u8>, Arc<[u8]>>;

/// The MSA plan: the in-region predicate (0-based half-open per contig) and
/// the reference sequence the fold validates each REF span against.
#[derive(Clone, Debug)]
pub(crate) struct MsaPlan {
    pub(crate) regions: MsaRegions,
    pub(crate) reference: Arc<ReferenceMap>,
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
#[derive(Clone, Default)]
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
    /// At least one allele is called (half-calls count; doc_v30 R01). Used by
    /// the count-family gates (`--012`, `--hapmap`, `--msa`).
    pub(crate) present: bool,
    /// vcf2phylip's presence rule for the phylo `--min-samples-locus` gate:
    /// the raw GT text does not start with `.`, so `./1` is missing and
    /// `0/.` is present (doc_v30 §4.6(a)). The two rules differ only on
    /// half calls, and the phylo targets are documented as vcf2phylip
    /// replacements, so they keep its rule.
    pub(crate) phylo_present: bool,
    /// The raw GT subfield range within the samples region (for GEN_BIN).
    /// `None` when the sample column has no GT subfield (`FORMAT=DP:GT` with
    /// a `30`-only column); callers must treat that as missing, never slice
    /// the sample column itself.
    pub(crate) raw: Option<std::ops::Range<u32>>,
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
    /// Whether the current record's malformed skip was already counted. A
    /// bad record is one input line, not one failure per enabled target
    /// (v29 V29-09 / D29-07).
    record_malformed: bool,
    /// The composed statistics scratch.
    stats: StatsScratch,
}

#[derive(Default)]
struct WarnFlags {
    malformed: bool,
    multiallelic_012: bool,
    padded_msa: bool,
    polyploid: bool,
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

/// One MSA site: every kept sample's value, all `width` bytes long, plus the
/// reference span the record consumes (`ref_len`; `width >= ref_len`, the
/// extra columns are insertion gaps).
pub(crate) struct MsaSite {
    pub(crate) pos0: u64,
    pub(crate) ref_len: u32,
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
    /// The run's blob budget (doc_v30 R09): a quarter of `--memory`, with
    /// the historical 8 MiB as the floor rather than the cap.
    spill_budget: usize,
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
                (0..phylo.cols.len()).map(|_| RowsBlob::default()).collect()
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
            debug_assert_eq!(
                cells.len(),
                count.cells.len(),
                "invariant: kept-sample count"
            );
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
                    self.msa_scratch.extend_from_slice(&site.pos0.to_le_bytes());
                    self.msa_scratch
                        .extend_from_slice(&site.width.to_le_bytes());
                    self.msa_scratch
                        .extend_from_slice(&site.ref_len.to_le_bytes());
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
    /// The sink's blob budget for `budget`: a quarter of the run's memory
    /// budget, never below [`SPILL_TOTAL_BUDGET`] (doc_v30 R09 — the old
    /// fixed 8 MiB ignored `--memory`).
    fn spill_budget(budget: &ThreadBudget) -> usize {
        let share = usize::try_from(budget.mem.bytes / 4).unwrap_or(usize::MAX);
        share.max(SPILL_TOTAL_BUDGET)
    }

    /// Bounds the sink's buffered blob bytes: past the run's spill budget
    /// the in-memory buffers are spilled whole in a fixed order (doc_v15
    /// D15-7; the bytes are unaffected, so the output is too).
    fn enforce_spill_budget(&mut self) -> Result<(), EngineError> {
        let mut total = self.buffered_bytes();
        if total <= self.spill_budget {
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
            if total <= self.spill_budget {
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
            blobs
                .as_ref()
                .map_or(0, |v| v.iter().map(RowsBlob::mem_len).sum())
        };
        cols(&self.phylo_cols)
            + cols(&self.phylo_bin_cols)
            + cols(&self.count_cells)
            + self
                .msa_runs
                .as_ref()
                .map_or(0, |runs| runs.iter().map(|(_, blob)| blob.mem_len()).sum())
            + self.phylo_used.mem_len()
            + self.count_pos.mem_len()
            + self.hmp_rows.mem_len()
    }
}

/// Runs the conversion over a path (or stdin when `None`).
///
/// # Errors
///
/// Request validation (Usage) and transport/render failures (Runtime).
pub fn run_convert(
    path: Option<&std::path::Path>,
    request: &ConvertRequest,
    ctx: &RunContext<'_>,
) -> Result<ConvertReport, crate::Error> {
    request.validate()?;
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
        // `-o` is written verbatim by this route: refuse it naming the
        // hapmap input before the writer truncates it (v29 V29-03).
        if let Some(output) = output {
            let mut inputs: Vec<(&std::path::Path, &str)> = Vec::new();
            if let Some(input) = path {
                inputs.push((input, "input"));
            }
            super::plan::validate_outputs(std::slice::from_ref(output), &inputs)?;
        }
        return crate::convert::hapmap::run_hapmap_to_vcf(path, &out, ctx);
    }
    let opened = crate::input::open(path).map_err(EngineError::from)?;
    run_convert_opened(opened, request, ctx)
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
    ctx: &RunContext<'_>,
) -> Result<ConvertReport, crate::Error> {
    let budget = ctx.budget;
    request.validate()?;
    let sample_names = opened.facts().sample_names.clone();
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
    let mut reference: Option<Arc<ReferenceMap>> = None;
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

    // The complete output set is known before the first writer opens: refuse
    // two targets writing the same file, and any output naming the input, the
    // MSA reference or the populations file (v29 V29-03 / D29-04 A).
    let mut planned = writers::planned_outputs(
        request.prefix.as_path(),
        phylo.as_ref(),
        count,
        hmp,
        &msa_outputs,
    );
    if let (Some(_), Some(plan)) = (&sfs_request, &sfs_plan) {
        planned.extend(crate::stats::writers::sfs::planned_paths(
            request.prefix.as_path(),
            &group_names,
            plan,
        ));
    }
    let mut guarded: Vec<(&std::path::Path, &str)> = Vec::new();
    if let Some(path) = opened.path() {
        guarded.push((path, "input"));
    }
    for target in &request.targets {
        if let ConvertTarget::Msa { reference, .. } = target {
            guarded.push((reference.as_path(), "reference FASTA"));
        }
    }
    if let crate::stats::request::Populations::File(path) = &request.populations {
        guarded.push((path.as_path(), "populations file"));
    }
    super::plan::validate_outputs(&planned, &guarded)?;

    let abort = ctx.abort;
    let contig_len: Arc<HashMap<Vec<u8>, u64>> = Arc::new(
        opened
            .facts()
            .contig_names
            .iter()
            .zip(opened.facts().contig_lengths.iter().copied())
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
        spill_budget: ConvertSink::spill_budget(budget),
        out: Arc::clone(&out),
    };
    let report = {
        // The text path (the InputJob wrapper renders BCF records to text):
        // the native record path is a statistics-only optimization.
        let job = crate::bcf_input::InputJob::new(job, &opened, abort, budget);
        let (mut source, _facts, _meta) = opened.into_parts();
        let report = run(source.as_mut(), &job, sink, budget, abort)?;
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
    let polyploid = merged.counters.polyploid;
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
    if let (Some((_, int_dtype)), Some(plan)) = (&sfs_request, &sfs_plan) {
        // A requested SFS target always writes its outputs (v27-05), even
        // when no chunk ever reached the sink: a header-only file has no
        // records, so the spectrum is all zero. `SfsPartial` is seeded lazily
        // by the first projector site, and an input where every record is
        // skipped (all multiallelic, ALT `.` with a GT call, out-of-domain
        // GT, polyploid) leaves it unseeded. Both shapes are the legitimate
        // zero spectrum `stats` writes for a 0-site request (v29 V29-05).
        let empty;
        let partial_sfs = match merged.sfs.as_ref().and_then(StatsPartial::sfs_partial) {
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
        n_records = merged.sfs.as_ref().map_or(0, StatsPartial::records);
    }
    Ok(ConvertReport {
        files,
        counters,
        sfs: sfs_report_out,
        n_records,
        polyploid,
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
        let (mut file, spill) = crate::spill::TempSpill::create(Some(&dir), "convert", "spill")
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
