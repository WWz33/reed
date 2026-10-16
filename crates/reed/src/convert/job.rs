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
use reed_io::vcf::line::LineRef;

use crate::convert::writers;
use crate::convert::{ConvertRequest, ConvertTarget, MsaRegion, PhyloFormats, PhyloOptions};
use crate::stats::driver::{build_sfs_plan, resolve_groups, sfs_report};
use crate::stats::request::SfsReport;
use crate::stats::run::{StatsJob, StatsPartial, StatsScratch};

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
    msa_lens: Vec<u8>,
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
    pub(crate) width: u8,
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
        if let Some(src) = &mut other.phylo {
            let dst = self.phylo.get_or_insert_with(|| PhyloPartial {
                n_sites: 0,
                n_bin_sites: 0,
                cols: vec![Vec::new(); src.cols.len()],
                bin_cols: vec![Vec::new(); src.bin_cols.len()],
                used: Vec::new(),
            });
            dst.n_sites += src.n_sites;
            dst.n_bin_sites += src.n_bin_sites;
            for (d, s) in dst.cols.iter_mut().zip(&mut src.cols) {
                d.append(s);
            }
            for (d, s) in dst.bin_cols.iter_mut().zip(&mut src.bin_cols) {
                d.append(s);
            }
            dst.used.append(&mut src.used);
        }
        if let Some(src) = &mut other.count {
            let dst = self.count.get_or_insert_with(|| CountPartial {
                n_sites: 0,
                cells: vec![Vec::new(); src.cells.len()],
                pos: Vec::new(),
            });
            dst.n_sites += src.n_sites;
            for (d, s) in dst.cells.iter_mut().zip(&mut src.cells) {
                d.append(s);
            }
            dst.pos.append(&mut src.pos);
        }
        if let Some(src) = &mut other.hmp {
            let dst = self.hmp.get_or_insert_with(HmpPartial::default);
            dst.n_sites += src.n_sites;
            dst.rows.append(&mut src.rows);
        }
        if let Some(src) = &mut other.msa {
            let dst = self.msa.get_or_insert_with(MsaPartial::default);
            for run in src.runs.drain(..) {
                match dst.runs.last_mut() {
                    Some(last) if last.contig == run.contig => last.sites.extend(run.sites),
                    _ => dst.runs.push(run),
                }
            }
        }
    }
}

/// Ordered sink: merges every chunk partial in input order into one state.
pub(crate) struct ConvertSink {
    merged: ConvertPartial,
    out: Arc<Mutex<Option<ConvertPartial>>>,
}

impl OrderedSink<ConvertPartial> for ConvertSink {
    fn accept(
        &mut self,
        _id: reed_core::ChunkId,
        partial: &mut ConvertPartial,
    ) -> Result<(), EngineError> {
        self.merged.absorb(partial);
        Ok(())
    }

    fn finish(self) -> Result<u64, EngineError> {
        *self
            .out
            .lock()
            .expect("invariant: convert slot not poisoned") = Some(self.merged);
        Ok(0)
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
        if let Some(stats) = &self.sfs {
            let ConvertScratch {
                sanitized,
                stats: stats_scratch,
                ..
            } = scratch;
            // A pure SFS run hands the body through untouched (exact
            // `stats --sfs` behavior); a mixed run drops the structurally
            // malformed lines the row targets skip (B13-2).
            let sfs_body = if self.has_row_targets() {
                sanitize_body(body, self.n_samples, sanitized)
            } else {
                body
            };
            partial.sfs = Some(stats.run_chunk(chunk, sfs_body, stats_scratch)?);
        }
        let mut line_no: u32 = 0;
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            let key = DrawKey {
                chunk: chunk.id.get(),
                line: line_no,
            };
            self.push_line(line, key, &mut partial, scratch);
            line_no += 1;
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
        if let Some(stats) = &self.sfs {
            let ConvertScratch {
                sanitized,
                stats: stats_scratch,
                ..
            } = scratch;
            let sfs_line = if self.has_row_targets() && !line_is_sane(line, self.n_samples) {
                b"".as_slice()
            } else if self.has_row_targets() {
                let _ = &mut *sanitized;
                line
            } else {
                line
            };
            if !sfs_line.is_empty() {
                partial.sfs = Some(stats.run_seam(id, sfs_line, stats_scratch)?);
            }
        }
        // The seam record belongs to the previous chunk's stream.
        let key = DrawKey {
            chunk: id.get().saturating_sub(1),
            line: u32::MAX,
        };
        self.push_line(line, key, &mut partial, scratch);
        Ok(partial)
    }
}

impl ConvertJob {
    fn has_row_targets(&self) -> bool {
        self.phylo.is_some() || self.count || self.hmp || self.msa.is_some()
    }

    fn push_line(
        &self,
        line: &[u8],
        key: DrawKey,
        partial: &mut ConvertPartial,
        scratch: &mut ConvertScratch,
    ) {
        // vcf2phylip counts every non-empty, non-`#` line as processed,
        // malformed ones included (`snp_num` increments before its checks).
        partial.counters.total_records += 1;
        let Some(lr) = LineRef::parse(line) else {
            self.bump_malformed(partial, scratch, "fewer than 9 columns");
            return;
        };
        if lr.n_sample_cols() != self.n_samples {
            self.bump_malformed(
                partial,
                scratch,
                "a sample-column count differing from the header",
            );
            return;
        }
        let Ok(pos1) = reed_io::vcf::fixed::parse_pos1(lr.pos_text(), lr.chrom()) else {
            self.bump_malformed(partial, scratch, "an unparsable POS");
            return;
        };
        let pos0 = pos1 - 1;

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
        if let Some(plan) = &self.phylo {
            if present < plan.gate {
                partial.counters.excess_missing += 1;
                return;
            }
            self.build_alleles(&lr, true, scratch);
            if !is_snp(&lr, &scratch.allele_ranges) {
                partial.counters.mnps += 1;
                return;
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
                        return;
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
    fn parse_gts(&self, lr: &LineRef<'_>, scratch: &mut ConvertScratch) {
        let samples = lr.samples_text();
        scratch.gts.clear();
        scratch.gts.resize(self.n_samples, GtSlots::default());
        let mut start = 0_usize;
        for gt in &mut scratch.gts {
            let begin = start.min(samples.len());
            let end =
                memchr::memchr(b'\t', &samples[begin..]).map_or(samples.len(), |rel| begin + rel);
            let column = &samples[begin..end];
            start = end + 1;
            let colon = memchr::memchr(b':', column).unwrap_or(column.len());
            let raw = &column[..colon];
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
                    match std::str::from_utf8(token)
                        .ok()
                        .and_then(|t| t.parse::<u32>().ok())
                    {
                        Some(idx) => gt.slots[tokens] = Some(idx),
                        None => gt.bad = true,
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
        if raw_alt_count(lr) > 2 {
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
    fn msa_values(&self, scratch: &mut ConvertScratch) -> Result<(Vec<u8>, u8), Malformed> {
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
            msa_lens.push((msa_vals.len() as u32 - start) as u8);
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

/// A malformed record for one target (skipped with a counter).
struct Malformed;

/// Whether one complete record line is structurally sane for the composed
/// statistics job (9+ columns and a sample-column count matching the
/// header).
fn line_is_sane(line: &[u8], n_samples: usize) -> bool {
    if line.iter().all(|&b| b == b'\n' || b == b'\r') {
        return true;
    }
    LineRef::parse(line).is_some_and(|lr| {
        lr.n_sample_cols() == n_samples
            && reed_io::vcf::fixed::parse_pos1(lr.pos_text(), lr.chrom()).is_ok()
    })
}

/// Copies `body`, dropping lines [`line_is_sane`] rejects; returns the input
/// untouched when nothing must be dropped (the common case).
fn sanitize_body<'a>(body: &'a [u8], n_samples: usize, out: &'a mut Vec<u8>) -> &'a [u8] {
    let mut clean = true;
    for line in body.split_inclusive(|&b| b == b'\n') {
        if !line_is_sane(line, n_samples) {
            clean = false;
            break;
        }
    }
    if clean {
        return body;
    }
    out.clear();
    for line in body.split_inclusive(|&b| b == b'\n') {
        if line_is_sane(line, n_samples) {
            out.extend_from_slice(line);
        }
    }
    out.as_slice()
}

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
fn read_reference(path: &PathBuf) -> Result<HashMap<Vec<u8>, Arc<[u8]>>, crate::Error> {
    let text = std::fs::read(path).map_err(|e| {
        crate::Error::usage(format!("cannot read --reference {}: {e}", path.display()))
    })?;
    let mut map: HashMap<Vec<u8>, Arc<[u8]>> = HashMap::new();
    let mut name: Option<Vec<u8>> = None;
    let mut seq: Vec<u8> = Vec::new();
    for line in text.split(|&b| b == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
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
    let keep: Vec<bool> = request
        .keep_samples
        .clone()
        .unwrap_or_else(|| vec![true; n_samples]);
    if keep.len() != n_samples {
        return Err(crate::Error::usage(format!(
            "sample keep mask has {} entries but the header declares {n_samples}",
            keep.len()
        )));
    }
    if !keep.iter().any(|&k| k) {
        return Err(crate::Error::usage("the sample subset selects no samples"));
    }
    let n_kept = keep.iter().filter(|&&k| k).count();
    let keep: Arc<[bool]> = keep.into();

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
                if !formats.any() {
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
    let sfs_plan = sfs_request
        .as_ref()
        .map(|(sfs, _)| build_sfs_plan(sfs, n_groups))
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
        out: Arc::clone(&out),
    };
    let report = {
        // The text path (the InputJob wrapper renders BCF records to text):
        // the native record path is a statistics-only optimization.
        let job = crate::bcf_input::InputJob::new(job, &opened, &abort);
        let mut source = opened.source;
        let report = run(source.as_mut(), &job, sink, budget, &abort)?;
        job.finish(report.n_chunks, true)?;
        report
    };
    let merged = out
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
            merged.phylo.as_ref(),
        )?);
    }
    if render_count {
        files.extend(writers::write_count_matrix(
            prefix,
            &kept_names,
            &keep,
            merged.count.as_ref(),
        )?);
    }
    if render_hmp {
        files.extend(writers::write_hapmap(
            prefix,
            &kept_names,
            merged.hmp.as_ref(),
        )?);
    }
    if let Some(reference) = &reference {
        files.extend(writers::write_msa(
            prefix,
            reference,
            &msa_outputs,
            &msa_names,
            merged.msa.as_ref(),
        )?);
    }
    let mut sfs_report_out = None;
    let mut n_records = merged.counters.total_records;
    if let (Some((_, int_dtype)), Some(plan), Some(partial)) =
        (&sfs_request, &sfs_plan, &merged.sfs)
    {
        let partial_sfs = partial
            .sfs_partial()
            .expect("invariant: the SFS job produced a spectra partial");
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
