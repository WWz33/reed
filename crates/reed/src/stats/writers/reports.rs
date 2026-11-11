//! Per-site and per-sample reports (02-cli.md §5.4): `--freq`, `--hardy`,
//! `--het`, `--missing`, `--depth`, written in plink2's column layouts.
//!
//! Workers render each requested row into per-file byte buffers and the
//! ordered sink only appends them (03-performance.md O1) — HWE, float
//! formatting and the row text never touch the reducer thread.
//!
//! The job assembles itself from the request (O2): `--depth` alone decodes
//! DP, `--het` alone builds the genotype planes, and per-site rows are only
//! rendered for the reports that print them. Genotypes are decoded by the
//! shared scanner and counted by the shared kernel.

// Casts in this file are bounded by construction: The %g exponent is bounded by the f64 range.
#![allow(clippy::cast_possible_truncation)]

use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use reed_core::ChunkId;
use reed_core::count::SiteGroupCounts;
use reed_core::genotype::{GenotypePolicies, PlaneBlock, PolicyCounters};
use reed_core::math::hwe::HweCache;
use reed_core::sample::SampleMask;
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{ChunkJob, OrderedSink, run};
use reed_io::vcf::genotype_scan::{GenotypeScanner, INT_MISSING, ScanOutput};
use reed_io::vcf::line::LineRef;
use reed_io::vcf::site::SiteView;
use reed_io::write::atomic::AtomicWriter;
use reed_stats::IndividualFAcc;
use reed_stats::summary::report::QcReportAcc;

use crate::stats::driver::resolve_groups;
use crate::stats::request::Populations;

/// Which reports to write.
#[derive(Clone, Debug, Default)]
pub struct StatsReportRequest {
    /// `--freq`: `PREFIX.afreq`.
    pub freq: bool,
    /// `--hardy`: `PREFIX.hardy`.
    pub hardy: bool,
    /// `--het`: `PREFIX.het` (per-sample inbreeding F, plink's moment form).
    pub het: bool,
    /// `--missing`: `PREFIX.smiss` and `PREFIX.vmiss`.
    pub missing: bool,
    /// `--depth`: `PREFIX.idepth` (per sample) and `PREFIX.ldepth.mean`
    /// (per site) in vcftools' layouts.
    pub depth: bool,
    /// `--summary`: also accumulate the bcftools-stats QC summary in this
    /// same pass (doc_v4 O-1: the two families used to read twice).
    pub summary: bool,
    /// `--hist LIST` (doc_v14 02 §3): built-in metrics (`qual`, `dp`,
    /// `dp-site`, `gq`, `idl`, `af`) and `INFO/TAG`/`FMT/TAG` requests, in
    /// the order given. Empty = no histograms.
    pub hist: Vec<crate::stats::writers::histogram::HistMetric>,
    /// `--singletons`: one row per singleton/doubleton allele
    /// (`PREFIX.singletons`, the vcftools layout).
    pub singletons: bool,
    /// Population assignment (`--pop FILE`, doc_v14 02 §2.2): stratifies
    /// `--freq`/`--missing` (`PREFIX.afreq.pop`, `PREFIX.vmiss.pop`) and adds
    /// the POP column to `PREFIX.singletons`. `None` = the whole cohort (the
    /// 0.3.x behavior; `Populations::All` also means that here).
    pub populations: Option<Populations>,
    /// `--force-samples`: population-file rows naming samples absent from the
    /// header are skipped with a warning instead of failing.
    pub force_samples: bool,
    /// The record-level site selection (`-r/-R`, `-t/-T`, `--sites-file`),
    /// three families that **AND** (doc_v27 V27-14/V27-15). Report mode used
    /// to take only the streamed target pair, which is why `--sites-file`
    /// never reached it (doc_v7 D7-6) and why `-r` could not be combined
    /// with `-t` the way `reed filter` combines them.
    pub select: reed_filter::spec::region::SiteSelector,
    /// Output prefix (`-o PREFIX`).
    pub prefix: PathBuf,
}

impl StatsReportRequest {
    /// Whether any report is requested.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        !(self.freq
            || self.hardy
            || self.het
            || self.missing
            || self.depth
            || self.summary
            || self.singletons)
            && self.hist.is_empty()
    }
}

/// Validates the header types the requested reports read (V22-R2-6).
///
/// The filter already types-checks its depth mask at compile time (V22-R7);
/// the reports used to accept a Float `DP`/`AD`/`GQ` and silently treat
/// every value as missing. The policies now match the readers:
///
/// - depth consumers (`--depth`, `--summary`'s DP/PSC sections, `--hist dp`)
///   need an Integer `FORMAT/DP`, or an Integer `FORMAT/AD` when the header
///   declares no DP. `--hist dp` additionally requires one of them
///   (otherwise it would write an empty histogram);
/// - `--hist gq` needs an Integer `FORMAT/GQ`;
/// - `--hist dp-site` needs an Integer `INFO/DP`.
///
/// `--summary`/`--depth` without either depth tag stay accepted: the depth
/// sections are simply empty, as vcftools/bcftools report for such a
/// callset.
///
/// # Errors
///
/// A usage error naming the tag and its declared type.
pub fn validate_report_tags(
    request: &StatsReportRequest,
    tags: &reed_core::tag_table::TagTable,
) -> Result<(), crate::error::Error> {
    use crate::error::Error;
    use reed_core::tag_table::TagType;
    let text = |ty: TagType| format!("{ty:?}");
    let needs_depth = request.depth || request.summary;
    let hist_dp = request
        .hist
        .iter()
        .any(crate::stats::writers::histogram::HistMetric::needs_dp);
    if needs_depth || hist_dp {
        match tags.format("DP") {
            Some(def) if def.ty != TagType::Integer => {
                return Err(Error::usage(format!(
                    "FORMAT field DP is {}; depth reports need an Integer DP (or an Integer AD \
                     fallback)",
                    text(def.ty)
                )));
            }
            Some(_) => {}
            None => {
                if let Some(def) = tags.format("AD")
                    && def.ty != TagType::Integer
                {
                    return Err(Error::usage(format!(
                        "FORMAT field DP is not declared and FORMAT/AD is {}; the depth fallback \
                         needs Integer AD",
                        text(def.ty)
                    )));
                }
                if hist_dp && tags.format("AD").is_none() {
                    return Err(Error::usage(
                        "--hist \"dp\": the header declares no Integer FORMAT/DP or FORMAT/AD",
                    ));
                }
            }
        }
    }
    for metric in &request.hist {
        match metric {
            crate::stats::writers::histogram::HistMetric::Gq => match tags.format("GQ") {
                Some(def) if def.ty != TagType::Integer => {
                    return Err(Error::usage(format!(
                        "--hist \"gq\": FORMAT/GQ is {}; the GQ histogram needs an Integer GQ",
                        text(def.ty)
                    )));
                }
                Some(_) => {}
                None => {
                    return Err(Error::usage(
                        "--hist \"gq\": the header declares no FORMAT/GQ",
                    ));
                }
            },
            crate::stats::writers::histogram::HistMetric::DpSite => match tags.info("DP") {
                Some(def) if def.ty != TagType::Integer => {
                    return Err(Error::usage(format!(
                        "--hist \"dp-site\": INFO/DP is {}; the site-depth histogram needs an \
                         Integer DP",
                        text(def.ty)
                    )));
                }
                Some(_) => {}
                None => {
                    return Err(Error::usage(
                        "--hist \"dp-site\": the header declares no INFO/DP",
                    ));
                }
            },
            _ => {}
        }
    }
    Ok(())
}

/// Files written and records seen.
#[derive(Clone, Debug)]
pub struct StatsOutcome {
    /// Paths written, in the order of the flags.
    pub files: Vec<PathBuf>,
    /// Records scanned.
    pub n_records: u64,
    /// The merged QC summary when `summary` was requested (`render_bcftools_stats`
    /// turns it into `PREFIX.stats`).
    pub qc: Option<QcReportAcc>,
    /// The merged histograms when `--hist` was requested: one entry per
    /// metric in request order.
    pub hist: Vec<crate::stats::writers::histogram::HistOut>,
    /// The engine's stage metrics (doc_v4 03 A4-5).
    pub metrics: reed_engine::metrics::StageMetrics,
}

/// The reports that need per-site rows, and the fields each needs.
#[derive(Clone, Copy, Debug, Default)]
struct ReportPlan {
    freq: bool,
    hardy: bool,
    vmiss: bool,
    ldepth: bool,
    smiss: bool,
    het: bool,
    depth: bool,
    summary: bool,
    hist: bool,
    singletons: bool,
    /// `--pop`: per-population afreq/vmiss rows (doc_v14 02 §2.2).
    pop_freq: bool,
    pop_vmiss: bool,
}

impl ReportPlan {
    fn from_request(request: &StatsReportRequest, n_groups: usize) -> Self {
        let pop = n_groups > 0;
        Self {
            freq: request.freq,
            hardy: request.hardy,
            vmiss: request.missing,
            ldepth: request.depth,
            smiss: request.missing,
            het: request.het,
            depth: request.depth,
            summary: request.summary,
            hist: !request.hist.is_empty(),
            singletons: request.singletons,
            pop_freq: pop && request.freq,
            pop_vmiss: pop && request.missing,
        }
    }

    /// Whether any per-site row is printed.
    const fn sites(&self) -> bool {
        self.freq || self.hardy || self.vmiss || self.ldepth || self.pop_freq || self.pop_vmiss
    }
}

/// `--pop`: the per-group counts at one site (stratified afreq/vmiss rows).
#[derive(Default)]
struct PopSiteCounts {
    /// ALT allele count per group.
    alt: Vec<u32>,
    /// Called allele count per group (the OBS_CT denominator).
    an: Vec<u32>,
    /// Fully-missing genotype count per group.
    missing: Vec<u32>,
}

/// Worker → sink partial: rendered rows plus the per-sample accumulators.
#[derive(Default)]
struct ReportPartial {
    /// `--hist`: the generic-tag and GQ accumulators for this chunk.
    hist: Option<crate::stats::writers::histogram::HistAcc>,
    /// The QC summary pieces this chunk accumulated (`--summary`).
    qc: Option<QcReportAcc>,
    /// The `--hist` pieces this chunk accumulated.
    freq: Vec<u8>,
    /// `--singletons`: rendered `PREFIX.singletons` rows.
    singletons: Vec<u8>,
    hardy: Vec<u8>,
    vmiss: Vec<u8>,
    ldepth: Vec<u8>,
    /// `--pop --freq`: rendered `PREFIX.afreq.pop` rows.
    freq_pop: Vec<u8>,
    /// `--pop --missing`: rendered `PREFIX.vmiss.pop` rows.
    vmiss_pop: Vec<u8>,
    /// Per sample: fully missing genotype count (`--missing`).
    missing: Vec<u64>,
    /// Per sample: plink `--het` moment components.
    individual: Option<IndividualFAcc>,
    /// Per sample: `(Σ DP, sites with a called genotype and a DP)` (`--depth`).
    depth_sums: Vec<(f64, u64)>,
    n_records: u64,
}

/// The scan job: one record → rendered rows + per-sample updates (+ the QC
/// summary sections when `--summary` shares the pass).
struct ReportJob {
    n_samples: usize,
    policies: GenotypePolicies,
    plan: ReportPlan,
    /// The `--hist` request (empty = no histograms).
    hist: Vec<crate::stats::writers::histogram::HistMetric>,
    /// Sample names (the `--singletons` INDV column).
    sample_names: Vec<String>,
    /// `--pop`: each sample's group slot (`None` = not in the populations
    /// file) and the group names in first-appearance order (empty = no
    /// `--pop`).
    groups: Option<(Vec<Option<u32>>, Vec<String>)>,
    /// The record-level site selection (`-r/-R` ∧ `-t/-T` ∧
    /// `--sites-file`, doc_v27 V27-14/V27-15).
    select: reed_filter::spec::region::SiteSelector,
}

/// Reused per-worker buffers.
struct ReportScratch {
    scanner: GenotypeScanner,
    /// Genotype planes: only when `--het` needs the per-sample moments.
    planes: Option<PlaneBlock>,
    counters: PolicyCounters,
    mask: SampleMask,
    /// The QC summary accumulator (`--summary`, doc_v4 O-1).
    qc: Option<crate::stats::writers::summary::QcAccumulator>,
    /// The `--hist` accumulators (`None` unless requested).
    hist: Option<crate::stats::writers::histogram::HistAcc>,
    /// HWE memo behind the `--hardy` p column (A6-16).
    hwe: HweCache,
}

impl ReportScratch {
    // The flags mirror `ReportPlan`'s; grouping them into a struct would
    // just restate it.
    #[allow(clippy::fn_params_excessive_bools)]
    fn new(
        n_samples: usize,
        policies: GenotypePolicies,
        het: bool,
        depth: bool,
        summary: bool,
        hist: Option<&[crate::stats::writers::histogram::HistMetric]>,
        singletons: bool,
    ) -> Self {
        use crate::stats::writers::histogram::HistMetric;
        let qc = summary || hist.is_some_and(HistMetric::any_needs_qc);
        // `--singletons` reads the raw allele indices like the QC summary.
        let needs_alleles = qc || singletons;
        // The QC accumulator's DP sections and `--hist dp` both read the
        // depth column: FMT/DP, or the AD sum when the record has no DP key
        // (bcftools `calc_sample_depth`).
        let needs_dp = depth || summary || hist.is_some_and(|m| m.iter().any(HistMetric::needs_dp));
        Self {
            // DP is decoded for `--depth`, `--summary` (whose DP sections
            // read the genotype depth) and `--hist dp`; the QC accumulator
            // also needs the raw allele indices for its per-ALT counts
            // (doc_v5 03 §5).
            scanner: {
                let scanner = if needs_dp {
                    GenotypeScanner::new_with_specs(
                        n_samples,
                        &[reed_io::vcf::format_layout::IntSpec::depth()],
                        &[],
                        policies,
                    )
                } else {
                    GenotypeScanner::new(n_samples, &[], &[], policies)
                };
                // The QC accumulator (which `--hist`'s built-in metrics
                // also ride) and `--singletons` read the raw allele indices
                // (doc_v5 03 §5).
                if needs_alleles {
                    scanner.tracking_alleles()
                } else {
                    scanner
                }
            },
            planes: het.then(|| PlaneBlock::new(1, n_samples, false)),
            counters: PolicyCounters::default(),
            mask: SampleMask::all(n_samples),
            qc: qc.then(|| crate::stats::writers::summary::QcAccumulator::new(n_samples)),
            hist: hist.map(crate::stats::writers::histogram::HistAcc::new),
            hwe: HweCache::new(),
        }
    }
}

impl ChunkJob for ReportJob {
    type Partial = ReportPartial;
    type Scratch = Option<ReportScratch>;

    fn run_chunk(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| {
            ReportScratch::new(
                self.n_samples,
                self.policies,
                self.plan.het,
                self.plan.depth,
                self.plan.summary,
                self.plan.hist.then_some(&self.hist),
                self.plan.singletons,
            )
        });
        let mut partial = self.new_partial();
        // `push_line` classifies: blank lines are skipped (a seam can hand
        // the report job `"\n"`), malformed lines stay fatal (doc_v25 V25-03).
        for line in body.split(|&b| b == b'\n') {
            self.push_line(line, scratch, &mut partial)?;
        }
        self.drain_qc(scratch, &mut partial);
        Ok(partial)
    }

    fn run_seam(
        &self,
        _id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| {
            ReportScratch::new(
                self.n_samples,
                self.policies,
                self.plan.het,
                self.plan.depth,
                self.plan.summary,
                self.plan.hist.then_some(&self.hist),
                self.plan.singletons,
            )
        });
        let mut partial = self.new_partial();
        self.push_line(line, scratch, &mut partial)?;
        self.drain_qc(scratch, &mut partial);
        Ok(partial)
    }
}

impl ReportJob {
    /// Moves the worker's QC accumulator into the chunk partial, leaving a
    /// fresh one behind.
    fn drain_qc(&self, scratch: &mut ReportScratch, partial: &mut ReportPartial) {
        if let Some(qc) = &mut scratch.qc {
            partial.qc = Some(qc.take_acc());
        }
        if let Some(hist) = &mut scratch.hist {
            partial.hist = Some(std::mem::replace(
                hist,
                crate::stats::writers::histogram::HistAcc::new(&self.hist),
            ));
        }
    }

    /// Allocates the partial's per-sample accumulators for this request.
    fn new_partial(&self) -> ReportPartial {
        ReportPartial {
            missing: if self.plan.smiss {
                vec![0; self.n_samples]
            } else {
                Vec::new()
            },
            individual: self.plan.het.then(|| IndividualFAcc::new(self.n_samples)),
            depth_sums: if self.plan.depth {
                vec![(0.0, 0); self.n_samples]
            } else {
                Vec::new()
            },
            ..ReportPartial::default()
        }
    }

    fn push_line(
        &self,
        line: &[u8],
        scratch: &mut ReportScratch,
        partial: &mut ReportPartial,
    ) -> Result<(), EngineError> {
        let lr = match reed_io::vcf::line::classify_line(line, self.n_samples) {
            reed_io::vcf::line::RecordClass::Record(lr) => lr,
            reed_io::vcf::line::RecordClass::Blank => return Ok(()),
            reed_io::vcf::line::RecordClass::Malformed(why) => {
                return Err(EngineError::job(format!(
                    "record line is malformed ({})",
                    why.label()
                )));
            }
        };
        let site = SiteView::parse(&lr).map_err(EngineError::from)?;
        let pos1 = site.pos1;
        // The site selection gates every record before any accumulator (the
        // same predicate the windowed path applies); `-r/-R`, `-t/-T` and
        // `--sites-file` intersect (doc_v27 V27-14/V27-15).
        if !self.select.is_empty() && !self.select.allows(site.contig_name, pos1) {
            return Ok(());
        }
        // The scanner adds its policy triggers; reset the per-record slot.
        scratch.counters = PolicyCounters::default();
        let out = scratch
            .scanner
            .scan(&lr, &site, &mut scratch.counters)
            .map_err(EngineError::from)?;

        // Counts straight from the codes — the planes are built only when
        // `--het` needs the per-sample moments (O2).
        let mut counts = SiteGroupCounts::default();
        let mut n_missing = 0_u32;
        for &code in out.codes {
            counts.add_code(code);
            if code.is_fully_missing() {
                n_missing += 1;
            }
        }
        if self.plan.smiss {
            for (slot, &code) in partial.missing.iter_mut().zip(out.codes) {
                if code.is_fully_missing() {
                    *slot += 1;
                }
            }
        }

        // vcftools depth: a sample contributes when its genotype is called
        // and the DP value is present.
        let mut site_depth = (0.0_f64, 0.0_f64, 0_u64);
        if self.plan.depth {
            let dp_column = out.int_columns.first().map_or(&[][..], Vec::as_slice);
            let mut sum = 0.0_f64;
            let mut sq = 0.0_f64;
            let mut n = 0_u64;
            for (i, &code) in out.codes.iter().enumerate() {
                let dp = dp_column.get(i).copied().unwrap_or(INT_MISSING);
                if code.is_fully_missing() || dp == INT_MISSING {
                    continue;
                }
                let dp = f64::from(dp);
                partial.depth_sums[i].0 += dp;
                partial.depth_sums[i].1 += 1;
                sum += dp;
                sq += dp * dp;
                n += 1;
            }
            site_depth = (sum, sq, n);
        }

        // `--pop`: the per-group alt-allele, called-allele and fully-missing
        // counts at this site (the stratified afreq/vmiss rows). Samples
        // outside every group are invisible to the stratified rows, as in a
        // per-population run.
        let mut pop = PopSiteCounts::default();
        if self.plan.pop_freq || self.plan.pop_vmiss {
            let Some((slots, names)) = &self.groups else {
                return Ok(());
            };
            let n_groups = names.len();
            pop = PopSiteCounts {
                alt: vec![0; n_groups],
                an: vec![0; n_groups],
                missing: vec![0; n_groups],
            };
            for (i, &code) in out.codes.iter().enumerate() {
                let Some(slot) = slots.get(i).copied().flatten() else {
                    continue;
                };
                let g = slot as usize;
                if g >= n_groups {
                    continue;
                }
                // The same counting class `SiteCounts::from_codes` uses,
                // scoped to one group (doc_v25 V25-07a).
                use reed_core::genotype::GtClass;
                if code.class() == GtClass::Missing {
                    pop.missing[g] += 1;
                    continue;
                }
                pop.an[g] += code.n_ref() + code.n_alt();
                pop.alt[g] += code.n_alt();
            }
        }
        // plink `--het` moment components at biallelic sites.
        let n_group = u32::try_from(self.n_samples).expect("invariant: sample count < u32::MAX");
        if self.plan.het {
            let planes = scratch
                .planes
                .as_mut()
                .expect("invariant: --het builds the planes");
            planes.push_codes(0, out.codes);
            let n_called = counts.n_called_alleles(n_group);
            if n_called > 0 {
                let p_alt = f64::from(counts.n_alt_alleles) / f64::from(n_called);
                if let Some(individual) = &mut partial.individual {
                    individual.push_site(planes, 0, &scratch.mask, p_alt);
                }
            }
        }

        if self.plan.sites() {
            self.render_site(
                &lr,
                &site,
                pos1,
                &counts,
                n_missing,
                site_depth,
                &pop,
                &mut scratch.hwe,
                partial,
            );
        }
        // The QC summary consumes the scan above (doc_v5 03 §5): the reports
        // scanner already carries the raw allele indices and the DP column,
        // so there is no second GT parse.
        if let Some(qc) = &mut scratch.qc {
            qc.push_scan(&out, &lr, &site)?;
        }
        // The scan output borrows the scratch scanner, so the hist
        // accumulator moves out for the call and back afterwards.
        if let Some(mut acc) = scratch.hist.take() {
            self.push_hist(&lr, &mut acc);
            scratch.hist = Some(acc);
        }
        if self.plan.singletons {
            self.render_singletons(&lr, &out, &mut partial.singletons);
        }
        partial.n_records += 1;
        Ok(())
    }

    /// Renders the `--singletons` rows for one record (vcftools layout:
    /// `CHROM POS SINGLETON/DOUBLETON ALLELE INDV`). Semantics follow
    /// vcftools `output_singletons` (variant_file_output.cpp:2715): per
    /// allele index including REF, `S` for a count of 1 and `D` for a count
    /// of 2 with a homozygous carrier. With `--pop` a sixth column carries
    /// the carrier's population (`.` for a sample outside every group);
    /// counting still runs over the whole cohort (doc_v14 02 §2.2).
    fn render_singletons(&self, lr: &LineRef<'_>, out: &ScanOutput<'_>, buf: &mut Vec<u8>) {
        let Some(alleles) = out.alleles else {
            return;
        };
        let alts = lr.alts();
        let n_alt = if alts == b"." {
            0
        } else {
            memchr::memchr_iter(b',', alts).count() + 1
        };
        let n_alleles = n_alt + 1;
        let mut counts = vec![0_u32; n_alleles];
        for index in alleles.iter().flatten() {
            if let Some(slot) = counts.get_mut(usize::from(*index)) {
                *slot += 1;
            }
        }
        for (index, &count) in counts.iter().enumerate() {
            if count != 1 && count != 2 {
                continue;
            }
            let want_hom = count == 2;
            for (sample, token) in alleles.chunks(2).enumerate() {
                let index = u16::try_from(index).unwrap_or(u16::MAX);
                let carries = if want_hom {
                    token == [Some(index), Some(index)]
                } else {
                    token.iter().flatten().any(|a| *a == index)
                };
                if !carries {
                    continue;
                }
                let allele: &[u8] = if index == 0 {
                    lr.reference()
                } else {
                    alt_at(alts, usize::from(index) - 1)
                };
                let kind = if want_hom { b'D' } else { b'S' };
                let name = self.sample_names.get(sample).map_or("?", String::as_str);
                let _ = match &self.groups {
                    Some((slots, names)) => writeln!(
                        buf,
                        "{}\t{}\t{}\t{}\t{}\t{}",
                        String::from_utf8_lossy(lr.chrom()),
                        String::from_utf8_lossy(lr.pos_text()),
                        kind as char,
                        String::from_utf8_lossy(allele),
                        name,
                        slots
                            .get(sample)
                            .copied()
                            .flatten()
                            .and_then(|g| names.get(g as usize))
                            .map_or_else(|| ".".to_string(), String::clone)
                    ),
                    None => writeln!(
                        buf,
                        "{}\t{}\t{}\t{}\t{}",
                        String::from_utf8_lossy(lr.chrom()),
                        String::from_utf8_lossy(lr.pos_text()),
                        kind as char,
                        String::from_utf8_lossy(allele),
                        name
                    ),
                };
                break;
            }
        }
    }

    /// Feeds the `--hist` accumulators for one record (doc_v14 02 §3.3): the
    /// built-in metrics ride the QC accumulator (same pass); `gq` and the
    /// generic `FMT/`/`INFO/` tags decode their own field here — the FORMAT
    /// keys are read once per record, so a tag costs one split per sample.
    fn push_hist(&self, lr: &LineRef<'_>, acc: &mut crate::stats::writers::histogram::HistAcc) {
        // FORMAT key offsets, computed lazily (only when a FMT metric asks).
        let mut fmt_keys: Option<()> = None;
        // doc_v28 V28-24: the FORMAT key position comes from the shared
        // colon walk (itself the `key_index` primitive).
        let key_index = |_keys: &mut Option<()>, tag: &[u8]| -> Option<usize> {
            reed_io::vcf::format_layout::key_index(lr.format_text(), tag)
        };
        for (idx, metric) in self.hist.iter().enumerate() {
            let tag_bytes: &[u8];
            let (is_fmt, kind, is_gq) = match metric {
                crate::stats::writers::histogram::HistMetric::Gq => {
                    tag_bytes = b"GQ";
                    (true, crate::stats::writers::histogram::HistKind::Int, true)
                }
                crate::stats::writers::histogram::HistMetric::Tag { format, tag, kind } => {
                    tag_bytes = tag.as_bytes();
                    (*format, *kind, false)
                }
                _ => continue,
            };
            if !is_fmt {
                // INFO/<TAG>: read the record's INFO field.
                let value = reed_io::vcf::info::info_fragments(lr.info_text())
                    .find(|(t, _)| *t == tag_bytes)
                    .and_then(|(_, v)| v);
                match (kind, value) {
                    (crate::stats::writers::histogram::HistKind::Int, Some(text)) => {
                        for token in text.split(|&b| b == b',') {
                            match parse_int_token(token) {
                                Some(v) => {
                                    if let Some(h) = acc.int.get_mut(idx).and_then(Option::as_mut) {
                                        h.push(v);
                                    }
                                }
                                None => bump_missing(acc, idx),
                            }
                        }
                    }
                    (crate::stats::writers::histogram::HistKind::Float, Some(text)) => {
                        for token in text.split(|&b| b == b',') {
                            match parse_float_token(token) {
                                Some(v) => {
                                    if let Some(f) = acc.float.get_mut(idx).and_then(Option::as_mut)
                                    {
                                        f.push(Some(v));
                                    }
                                }
                                None => bump_missing(acc, idx),
                            }
                        }
                    }
                    (_, None) => bump_missing(acc, idx),
                }
                continue;
            }
            // FMT/<TAG> (and `gq`): subfield `key_index` of every sample. A
            // key the record's FORMAT lacks is one missing observation per
            // sample (the opportunity count), not one per record.
            let Some(key) = key_index(&mut fmt_keys, tag_bytes) else {
                bump_missing_n(acc, idx, u64::try_from(self.n_samples).unwrap_or(u64::MAX));
                continue;
            };
            // The scanner accepts trailing empty columns (doc_v22 V22-5);
            // they are not samples, so never count them as missing
            // observations (V22-R4).
            for span in lr.sample_ranges().take(self.n_samples) {
                let sample = &lr.line[span];
                let token = reed_io::vcf::format_layout::subfield(sample, key).unwrap_or(b".");
                // A Number=A/G/R subfield carries several comma-separated
                // values; each is one observation (the INFO path's rule).
                for value in token.split(|&b| b == b',') {
                    match kind {
                        crate::stats::writers::histogram::HistKind::Int => {
                            match parse_int_token(value) {
                                Some(v) if is_gq => acc.gq.push(v),
                                Some(v) => {
                                    if let Some(h) = acc.int.get_mut(idx).and_then(Option::as_mut) {
                                        h.push(v);
                                    }
                                }
                                None => bump_missing(acc, idx),
                            }
                        }
                        crate::stats::writers::histogram::HistKind::Float => {
                            match parse_float_token(value) {
                                Some(v) => {
                                    if let Some(f) = acc.float.get_mut(idx).and_then(Option::as_mut)
                                    {
                                        f.push(Some(v));
                                    }
                                }
                                None => bump_missing(acc, idx),
                            }
                        }
                    }
                }
            }
        }
    }

    /// Renders every requested per-site row for one record (O1: the worker
    /// formats; the reducer only appends the bytes).
    // The columns come from the record's own fields; a column-group struct
    // would just restate the parameters.
    #[allow(clippy::too_many_arguments)]
    fn render_site(
        &self,
        lr: &LineRef<'_>,
        site: &SiteView<'_>,
        pos1: u64,
        counts: &SiteGroupCounts,
        n_missing: u32,
        depth: (f64, f64, u64),
        pop: &PopSiteCounts,
        hwe: &mut HweCache,
        partial: &mut ReportPartial,
    ) {
        let chrom = site.contig_name;
        let reference = site.reference;
        let alt = site.alts;
        let an = counts.n_called_alleles(
            u32::try_from(self.n_samples).expect("invariant: sample count < u32::MAX"),
        );
        // plink2's layout needs an identifiable row: when the input ID is
        // missing it synthesizes `chrom_pos` (same rule here).
        let id_missing = lr.id() == b"." || lr.id().is_empty();
        let write_prefix = |out: &mut Vec<u8>| {
            out.extend_from_slice(chrom);
            out.push(b'\t');
            if id_missing {
                out.extend_from_slice(chrom);
                out.push(b'_');
                let _ = write!(out, "{pos1}");
            } else {
                out.extend_from_slice(lr.id());
            }
        };

        if self.plan.freq {
            let out = &mut partial.freq;
            write_prefix(out);
            out.push(b'\t');
            out.extend_from_slice(reference);
            out.push(b'\t');
            out.extend_from_slice(alt);
            out.push(b'\t');
            if an > 0 {
                g6_into(out, f64::from(counts.n_alt_alleles) / f64::from(an));
            } else {
                out.extend_from_slice(b"NA");
            }
            let _ = write!(out, "\t{an}");
            out.push(b'\n');
        }
        if self.plan.hardy {
            let out = &mut partial.hardy;
            write_prefix(out);
            out.push(b'\t');
            out.extend_from_slice(reference);
            out.push(b'\t');
            out.extend_from_slice(alt);
            // A1 is the input's REF allele, AX the ALT — matching plink2's
            // VCF-import output (verified against plink2 2.x on this data:
            // A1 follows REF regardless of which allele is minor).
            let called = counts.n_full_called();
            let het = counts.n_het;
            let _ = write!(out, "\t{}\t{}\t{}", counts.n_hom_ref, het, counts.n_hom_alt);
            out.push(b'\t');
            if called > 0 {
                g6_into(out, f64::from(het) / f64::from(called));
            } else {
                out.extend_from_slice(b"NA");
            }
            out.push(b'\t');
            let p_a1 = if an > 0 {
                f64::from(counts.n_ref_alleles) / f64::from(an)
            } else {
                f64::NAN
            };
            if p_a1.is_nan() {
                out.extend_from_slice(b"NA");
            } else {
                g6_into(out, 2.0 * p_a1 * (1.0 - p_a1));
            }
            out.push(b'\t');
            let mac = 2 * counts.n_hom_ref.min(counts.n_hom_alt) + het;
            if called == 0 {
                out.extend_from_slice(b"NA");
            } else {
                g6_into(out, hwe.p_exact(called, mac, het));
            }
            out.push(b'\n');
        }
        if self.plan.vmiss {
            let out = &mut partial.vmiss;
            write_prefix(out);
            let obs = u64::try_from(self.n_samples).unwrap_or(u64::MAX);
            let _ = write!(out, "\t{n_missing}\t{obs}\t");
            if obs > 0 {
                g6_into(out, f64::from(n_missing) / obs as f64);
            } else {
                out.extend_from_slice(b"NA");
            }
            out.push(b'\n');
        }
        if self.plan.ldepth {
            let out = &mut partial.ldepth;
            out.extend_from_slice(chrom);
            let _ = write!(out, "\t{pos1}\t");
            let (sum, sq, n) = depth;
            if n > 0 {
                let mean = sum / n as f64;
                // Sample variance (n − 1), vcftools' VAR_DEPTH.
                let var = if n > 1 {
                    (sq - sum * sum / n as f64) / (n - 1) as f64
                } else {
                    0.0
                };
                g6_into(out, mean);
                out.push(b'\t');
                g6_into(out, var.max(0.0));
            } else {
                out.extend_from_slice(b"NA\tNA");
            }
            out.push(b'\n');
        }
        // `--pop`: one stratified row per group, in group order, after the
        // cohort files' row (the group order is the file's first-appearance
        // order, so the output is deterministic).
        if self.plan.pop_freq {
            let Some((_, names)) = &self.groups else {
                return;
            };
            for (g, name) in names.iter().enumerate() {
                let out = &mut partial.freq_pop;
                write_prefix(out);
                let _ = write!(out, "\t{name}\t");
                out.extend_from_slice(reference);
                out.push(b'\t');
                out.extend_from_slice(alt);
                out.push(b'\t');
                let an = pop.an[g];
                if an > 0 {
                    g6_into(out, f64::from(pop.alt[g]) / f64::from(an));
                } else {
                    out.extend_from_slice(b"NA");
                }
                let _ = write!(out, "\t{an}");
                out.push(b'\n');
            }
        }
        if self.plan.pop_vmiss {
            let Some((slots, names)) = &self.groups else {
                return;
            };
            let mut group_n = vec![0_u32; names.len()];
            for slot in slots.iter().flatten() {
                if let Some(n) = group_n.get_mut(*slot as usize) {
                    *n += 1;
                }
            }
            for (g, name) in names.iter().enumerate() {
                let out = &mut partial.vmiss_pop;
                write_prefix(out);
                let _ = write!(out, "\t{name}\t");
                let obs = group_n[g];
                let _ = write!(out, "{}\t{obs}\t", pop.missing[g]);
                if obs > 0 {
                    g6_into(out, f64::from(pop.missing[g]) / f64::from(obs));
                } else {
                    out.extend_from_slice(b"NA");
                }
                out.push(b'\n');
            }
        }
    }
}

/// Ordered consumer: appends the worker-rendered rows, folds per-sample
/// vectors.
struct ReportSink {
    request: StatsReportRequest,
    plan: ReportPlan,
    n_samples: usize,
    sample_names: Vec<String>,
    afreq: Option<AtomicSink>,
    singletons: Option<AtomicSink>,
    hardy: Option<AtomicSink>,
    vmiss: Option<AtomicSink>,
    ldepth: Option<AtomicSink>,
    /// `--pop --freq`: `PREFIX.afreq.pop`.
    freq_pop: Option<AtomicSink>,
    /// `--pop --missing`: `PREFIX.vmiss.pop`.
    vmiss_pop: Option<AtomicSink>,
    missing: Vec<u64>,
    depth_sums: Vec<(f64, u64)>,
    individual: IndividualFAcc,
    n_records: u64,
    /// The `--summary` accumulator (`None` unless it shares the pass).
    qc: Option<QcReportAcc>,
    /// The merged `--hist` accumulators.
    hist: crate::stats::writers::histogram::HistAcc,
    /// The rendered histograms, filled at `finish`.
    hist_out: Vec<crate::stats::writers::histogram::HistOut>,
    /// Prepared report files, published by the driver only after the input
    /// completeness check passes (doc_v28 V28-06).
    pending: Arc<Mutex<Vec<(PathBuf, reed_io::write::atomic::AtomicFile)>>>,
    out: Arc<Mutex<Option<StatsOutcome>>>,
}

impl ReportSink {
    /// Opens a report file atomically (doc_v26 V26-08): the returned sink
    /// writes to a sibling temporary and publishes on
    /// [`AtomicWriter::commit`].
    fn open(prefix: &Path, suffix: &str) -> Result<AtomicSink, EngineError> {
        let mut path = prefix.as_os_str().to_os_string();
        path.push(suffix);
        AtomicWriter::create(&PathBuf::from(&path)).map_err(EngineError::Io)
    }
}

/// One report file: a buffered writer over an atomic temporary.
///
/// The report sink opens its files *before* the run and streams rows into
/// them, so a plain `File::create` truncated the user's previous report at
/// start-up and left it truncated if the run failed. reed-io's
/// [`AtomicWriter`] owns the create→flush→close sequence (doc_v28 V28-23);
/// `prepare_commit` returns the unpublished temporary, which the driver
/// publishes only after `job.finish` (doc_v28 V28-06).
pub(super) type AtomicSink = AtomicWriter;

/// Closes one report file into its unpublished temporary. The driver
/// publishes the collected owners after `job.finish` (doc_v28 V28-06).
fn stage_report_file(
    pending: &Mutex<Vec<(PathBuf, reed_io::write::atomic::AtomicFile)>>,
    files: &mut Vec<PathBuf>,
    sink: AtomicSink,
    path: PathBuf,
) -> Result<(), EngineError> {
    let atomic = sink.finish().map_err(EngineError::Io)?;
    pending
        .lock()
        .expect("invariant: pending-report slot not poisoned")
        .push((path.clone(), atomic));
    files.push(path);
    Ok(())
}

impl OrderedSink<ReportPartial> for ReportSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut ReportPartial) -> Result<(), EngineError> {
        self.n_records += partial.n_records;
        let append = |writer: &mut Option<AtomicSink>, bytes: &[u8]| {
            if let Some(w) = writer {
                w.write_all(bytes).map_err(EngineError::Io)?;
            }
            Ok::<(), EngineError>(())
        };
        append(&mut self.afreq, &partial.freq)?;
        append(&mut self.singletons, &partial.singletons)?;
        append(&mut self.hardy, &partial.hardy)?;
        append(&mut self.vmiss, &partial.vmiss)?;
        append(&mut self.ldepth, &partial.ldepth)?;
        append(&mut self.freq_pop, &partial.freq_pop)?;
        append(&mut self.vmiss_pop, &partial.vmiss_pop)?;
        if self.plan.depth {
            for (dst, src) in self.depth_sums.iter_mut().zip(&partial.depth_sums) {
                dst.0 += src.0;
                dst.1 += src.1;
            }
        }
        if self.plan.smiss {
            for (dst, src) in self.missing.iter_mut().zip(&partial.missing) {
                *dst += src;
            }
        }
        if let Some(other) = partial.individual.take() {
            self.individual.merge(&other);
        }
        if let (Some(dst), Some(src)) = (&mut self.qc, partial.qc.take()) {
            dst.merge(&src);
        }
        if let Some(src) = partial.hist.take() {
            self.hist.merge(&src);
        }
        Ok(())
    }

    fn finish(mut self) -> Result<u64, EngineError> {
        let mut files = Vec::new();
        let path = |suffix: &str| {
            let mut p = self.request.prefix.as_os_str().to_os_string();
            p.push(suffix);
            PathBuf::from(p)
        };
        // doc_v28 V28-06: close every file into its temporary, then let the
        // driver publish after `job.finish` proves the input is complete.
        if let Some(w) = self.afreq.take() {
            stage_report_file(&self.pending, &mut files, w, path(".afreq"))?;
        }
        if let Some(w) = self.singletons.take() {
            stage_report_file(&self.pending, &mut files, w, path(".singletons"))?;
        }
        if let Some(w) = self.hardy.take() {
            stage_report_file(&self.pending, &mut files, w, path(".hardy"))?;
        }
        if let Some(w) = self.vmiss.take() {
            stage_report_file(&self.pending, &mut files, w, path(".vmiss"))?;
        }
        if let Some(w) = self.ldepth.take() {
            stage_report_file(&self.pending, &mut files, w, path(".ldepth.mean"))?;
        }
        if let Some(w) = self.freq_pop.take() {
            stage_report_file(&self.pending, &mut files, w, path(".afreq.pop"))?;
        }
        if let Some(w) = self.vmiss_pop.take() {
            stage_report_file(&self.pending, &mut files, w, path(".vmiss.pop"))?;
        }
        if self.request.depth {
            let mut w = ReportSink::open(&self.request.prefix, ".idepth")?;
            writeln!(w, "INDV\tN_SITES\tMEAN_DEPTH").map_err(EngineError::Io)?;
            for (i, &(sum, n)) in self.depth_sums.iter().enumerate() {
                let name = self.sample_names.get(i).map_or("?", String::as_str);
                let mean = if n > 0 {
                    g6(sum / n as f64)
                } else {
                    "NA".to_string()
                };
                writeln!(w, "{name}\t{n}\t{mean}").map_err(EngineError::Io)?;
            }
            stage_report_file(&self.pending, &mut files, w, path(".idepth"))?;
        }
        if self.request.het {
            let f = self.individual.f_values();
            let mut w = ReportSink::open(&self.request.prefix, ".het")?;
            writeln!(w, "#IID\tO(HOM)\tE(HOM)\tOBS_CT\tF").map_err(EngineError::Io)?;
            for i in 0..self.n_samples {
                let obs = self.individual.n_called[i];
                let o_hom = self.individual.o_hom[i];
                let e_hom = obs as f64 - self.individual.hetero_expected[i].value();
                let f_cell = match f.get(i).copied().flatten() {
                    Some(v) => g6(v),
                    None => "NA".to_string(),
                };
                let name = self.sample_names.get(i).map_or("?", String::as_str);
                writeln!(w, "{name}\t{o_hom}\t{}\t{obs}\t{f_cell}", g6(e_hom))
                    .map_err(EngineError::Io)?;
            }
            stage_report_file(&self.pending, &mut files, w, path(".het"))?;
        }
        if self.request.missing {
            // plink2's .smiss OBS_CT is the number of variants (the sample's
            // opportunity count), not the sample count.
            let obs = self.n_records;
            let mut w = ReportSink::open(&self.request.prefix, ".smiss")?;
            writeln!(w, "#IID\tMISSING_CT\tOBS_CT\tF_MISS").map_err(EngineError::Io)?;
            for (i, &missing) in self.missing.iter().enumerate() {
                let f_miss = if obs > 0 {
                    g6(missing as f64 / obs as f64)
                } else {
                    "NA".into()
                };
                let name = self.sample_names.get(i).map_or("?", String::as_str);
                writeln!(w, "{name}\t{missing}\t{obs}\t{f_miss}").map_err(EngineError::Io)?;
            }
            stage_report_file(&self.pending, &mut files, w, path(".smiss"))?;
        }
        // `--hist`: one `PREFIX.<m>.hist` per metric (doc_v14 02 §3.2).
        if self.request.hist.is_empty() {
            // nothing to write
        } else {
            let merged = crate::stats::writers::histogram::merge_hist(
                &self.request.hist,
                self.qc.as_ref(),
                &self.hist,
            );
            for out in &merged {
                let suffix = format!(".{}.hist", out.metric.infix());
                let mut w = ReportSink::open(&self.request.prefix, &suffix)?;
                let text = crate::stats::writers::histogram::render_hist(out);
                w.write_all(text.as_bytes()).map_err(EngineError::Io)?;
                stage_report_file(&self.pending, &mut files, w, path(&suffix))?;
            }
            // The library caller gets the same histograms the files carry
            // (`StatsOutcome::hist`); it used to be left empty.
            self.hist_out = merged;
        }
        // The QC accumulator may also run for `--hist`'s built-in metrics
        // (its data is consumed above); `PREFIX.stats` and `outcome.qc` exist
        // only when `--summary` asked for them.
        if !self.request.summary {
            self.qc = None;
        }
        let n_records = self.n_records;
        *self
            .out
            .lock()
            .expect("invariant: reports slot not poisoned") = Some(StatsOutcome {
            files,
            n_records,
            qc: self.qc.take(),
            hist: std::mem::take(&mut self.hist_out),
            // The engine metrics come from the run, not the sink; the
            // driver fills this after `run` returns.
            metrics: reed_engine::metrics::StageMetrics::default(),
        });
        Ok(n_records)
    }
}

/// The `index`-th comma-separated ALT allele of `alts`.
fn alt_at(alts: &[u8], index: usize) -> &[u8] {
    alts.split(|&b| b == b',').nth(index).unwrap_or(alts)
}

/// One `--hist` integer token, using the same parser as the filter's
/// integer columns (V22-R2-5): optional `+`/`-`, `.`/unparsable/negative/
/// outside-`i32` is missing. The old `str::parse::<u32>` accepted `+50` but
/// kept values above `i32::MAX`, diverging from the scanner.
fn parse_int_token(token: &[u8]) -> Option<u32> {
    let value = reed_io::vcf::genotype_scan::parse_int(token);
    (value != INT_MISSING).then_some(value)
}

/// One `--hist` float token: `.`, non-finite text (`NaN`, `inf`) and
/// unparsable text are missing — a non-finite value has no bin.
fn parse_float_token(token: &[u8]) -> Option<f32> {
    std::str::from_utf8(token)
        .ok()
        .and_then(|t| t.trim().parse::<f32>().ok())
        .filter(|v| v.is_finite())
}

/// Counts one unreadable observation for the metric at `idx`.
fn bump_missing(acc: &mut crate::stats::writers::histogram::HistAcc, idx: usize) {
    bump_missing_n(acc, idx, 1);
}

/// Counts `n` unreadable observations for the metric at `idx`.
fn bump_missing_n(acc: &mut crate::stats::writers::histogram::HistAcc, idx: usize, n: u64) {
    if let Some(m) = acc.missing.get_mut(idx) {
        *m += n;
    }
}

/// Appends `%g` with 6 significant digits, matching plink2's print precision.
fn g6_into(out: &mut Vec<u8>, x: f64) {
    if !x.is_finite() {
        out.extend_from_slice(b"NA");
        return;
    }
    if x == 0.0 {
        out.push(b'0');
        return;
    }
    let exponent = x.abs().log10().floor() as i32;
    let exponential = !(-4..=5).contains(&exponent);
    if exponential {
        // Rust prints `2.18867e6` / `1.23456e-3`; printf %g prints
        // `2.18867e+06` / `1.23456e-03`. The exponential form is rare, so a
        // small temporary is fine here.
        let raw = format!("{x:.5e}");
        if let Some((mantissa, exp)) = raw.split_once('e') {
            let exp: i32 = exp.parse().unwrap_or(0);
            let mantissa = trim_trailing_zeros(mantissa);
            let sign = if exp < 0 { '-' } else { '+' };
            let _ = write!(out, "{}e{sign}{:02}", mantissa, exp.abs());
        } else {
            out.extend_from_slice(raw.as_bytes());
        }
    } else {
        let decimals = (5 - exponent).max(0) as usize;
        let start = out.len();
        let _ = write!(out, "{x:.decimals$}");
        trim_trailing_zeros_in_place(out, start);
    }
}

/// `%g` with 6 significant digits as a `String` (per-sample files).
fn g6(x: f64) -> String {
    let mut buf = Vec::with_capacity(24);
    g6_into(&mut buf, x);
    String::from_utf8(buf).expect("invariant: g6 renders ASCII")
}

/// Strips trailing zeros (and a trailing dot) from a fixed-point rendering
/// appended at `start` of `out` (integers keep their zeros: no dot present).
fn trim_trailing_zeros_in_place(out: &mut Vec<u8>, start: usize) {
    if !out[start..].contains(&b'.') {
        return;
    }
    while out.len() > start && out.last() == Some(&b'0') {
        out.pop();
    }
    if out.len() > start && out.last() == Some(&b'.') {
        out.pop();
    }
}

/// Strips trailing zeros (and a trailing dot) from a fixed-point rendering.
fn trim_trailing_zeros(s: &str) -> String {
    if let Some((mantissa, exp)) = s.split_once('e') {
        return format!("{}e{exp}", trim_trailing_zeros(mantissa));
    }
    if s.contains('.') {
        let trimmed = s.trim_end_matches('0').trim_end_matches('.');
        return trimmed.to_string();
    }
    s.to_string()
}

/// Runs the requested reports over a BGZF VCF.
///
/// # Errors
///
/// Transport, job and output failures.
pub fn run_stats(
    input: Option<&Path>,
    request: &StatsReportRequest,
    budget: &ThreadBudget,
) -> Result<StatsOutcome, crate::error::Error> {
    let opened = crate::input::open(input).map_err(EngineError::from)?;
    run_stats_opened(opened, request, budget)
}

/// Runs the requested reports over an already-opened input (the CLI opens
/// once so stdin streams; doc_v3 03-performance.md O5).
///
/// # Errors
///
/// As [`run_stats`].
pub fn run_stats_opened(
    opened: crate::input::OpenedInput,
    request: &StatsReportRequest,
    budget: &ThreadBudget,
) -> Result<StatsOutcome, crate::error::Error> {
    let n_samples = opened.facts.n_samples();
    // Header types are a request property, checked before the first record
    // (V22-R2-6): a Float depth/GQ tag used to decode as all-missing.
    validate_report_tags(request, &reed_io::header::tag_table(&opened.facts.header))?;
    // `--pop` (doc_v14 02 §2.2): resolve the assignment against the header
    // before the run; `None` keeps the whole-cohort behavior. An empty group
    // set (every row skipped) is a usage error — the stratified files would
    // be all headers.
    let groups = match &request.populations {
        Some(populations) => {
            let (slots, names) = resolve_groups(
                populations,
                request.force_samples,
                &opened.facts.sample_names,
                None,
            )?;
            if names.is_empty() {
                return Err(crate::error::Error::usage(
                    "the populations file assigns no sample from the header",
                ));
            }
            Some((slots, names))
        }
        None => None,
    };
    let n_groups = groups.as_ref().map_or(0, |(_, names)| names.len());
    let plan = ReportPlan::from_request(request, n_groups);
    let job = ReportJob {
        n_samples,
        plan,
        hist: request.hist.clone(),
        sample_names: opened.facts.sample_names.clone(),
        groups,
        // The one stats decode policy (doc_v25 V25-07): the report job used
        // to hand-copy it, so a policy change had to touch two files.
        policies: crate::stats::run::stats_policies(),
        select: request.select.clone(),
    };
    let mut afreq = request
        .freq
        .then(|| ReportSink::open(&request.prefix, ".afreq"))
        .transpose()?;
    if let Some(w) = &mut afreq {
        writeln!(w, "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT").map_err(EngineError::Io)?;
    }
    let mut singletons = request
        .singletons
        .then(|| ReportSink::open(&request.prefix, ".singletons"))
        .transpose()?;
    if let Some(w) = &mut singletons {
        // With `--pop` the vcftools five columns gain the POP column
        // (doc_v14 02 §2.2: `.` for a carrier outside every group).
        if n_groups > 0 {
            writeln!(w, "CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV\tPOP")
                .map_err(EngineError::Io)?;
        } else {
            writeln!(w, "CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV")
                .map_err(EngineError::Io)?;
        }
    }
    let mut hardy = request
        .hardy
        .then(|| ReportSink::open(&request.prefix, ".hardy"))
        .transpose()?;
    if let Some(w) = &mut hardy {
        writeln!(
            w,
            "#CHROM\tID\tA1\tAX\tHOM_A1_CT\tHET_A1_CT\tTWO_AX_CT\tO(HET_A1)\tE(HET_A1)\tP"
        )
        .map_err(EngineError::Io)?;
    }
    let mut vmiss = request
        .missing
        .then(|| ReportSink::open(&request.prefix, ".vmiss"))
        .transpose()?;
    if let Some(w) = &mut vmiss {
        writeln!(w, "#CHROM\tID\tMISSING_CT\tOBS_CT\tF_MISS").map_err(EngineError::Io)?;
    }
    let mut ldepth = request
        .depth
        .then(|| ReportSink::open(&request.prefix, ".ldepth.mean"))
        .transpose()?;
    if let Some(w) = &mut ldepth {
        writeln!(w, "CHROM\tPOS\tMEAN_DEPTH\tVAR_DEPTH").map_err(EngineError::Io)?;
    }
    let mut freq_pop = plan
        .pop_freq
        .then(|| ReportSink::open(&request.prefix, ".afreq.pop"))
        .transpose()?;
    if let Some(w) = &mut freq_pop {
        writeln!(w, "#CHROM\tID\tPOP\tREF\tALT\tALT_FREQS\tOBS_CT").map_err(EngineError::Io)?;
    }
    let mut vmiss_pop = plan
        .pop_vmiss
        .then(|| ReportSink::open(&request.prefix, ".vmiss.pop"))
        .transpose()?;
    if let Some(w) = &mut vmiss_pop {
        writeln!(w, "#CHROM\tID\tPOP\tMISSING_CT\tOBS_CT\tF_MISS").map_err(EngineError::Io)?;
    }
    let abort = reed_engine::abort::abort_for(budget);
    let job = crate::bcf_input::InputJob::new(job, &opened, &abort, budget);
    let crate::input::OpenedInput {
        mut source, facts, ..
    } = opened;
    let pending_reports: Arc<Mutex<Vec<(PathBuf, reed_io::write::atomic::AtomicFile)>>> =
        Arc::new(Mutex::new(Vec::new()));
    let sink = ReportSink {
        request: request.clone(),
        plan,
        n_samples,
        sample_names: facts.sample_names.clone(),
        afreq,
        singletons,
        hardy,
        vmiss,
        ldepth,
        freq_pop,
        vmiss_pop,
        missing: vec![0; n_samples],
        depth_sums: vec![(0.0, 0); n_samples],
        individual: IndividualFAcc::new(n_samples),
        n_records: 0,
        qc: (request.summary
            || crate::stats::writers::histogram::HistMetric::any_needs_qc(&request.hist))
        .then(|| {
            reed_stats::summary::report::QcReportAcc::new(
                u32::try_from(n_samples).unwrap_or(u32::MAX),
            )
        }),
        hist: crate::stats::writers::histogram::HistAcc::new(&request.hist),
        hist_out: Vec::new(),
        pending: Arc::clone(&pending_reports),
        out: Arc::new(Mutex::new(None)),
    };
    let out = Arc::clone(&sink.out);
    let report = run(source.as_mut(), &job, sink, budget, &abort)?;
    job.finish(report.n_chunks, true)?;
    // Publish only after the BCF tail check: a failed run leaves every
    // previous report file byte-identical (doc_v28 V28-06).
    for (_path, atomic) in pending_reports
        .lock()
        .expect("invariant: pending-report slot not poisoned")
        .drain(..)
    {
        atomic.commit().map_err(EngineError::Io)?;
    }
    let mut outcome = out
        .lock()
        .expect("invariant: reports slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the outcome");
    outcome.metrics = report.metrics;
    Ok(outcome)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The buffer renderer keeps `%g`-6-digits semantics (float formatting
    /// must not drift: plink2 parity depends on it).
    #[test]
    fn g6_into_matches_the_string_form() {
        for x in [
            0.0,
            1.0,
            0.5,
            0.333_333_333_3,
            1e-5,
            1.234_567_89,
            123_456.0,
            1_234_567.0,
            0.000_123_456,
            2.18867e6,
            1.23456e-3,
            f64::NAN,
            1.0 / 3.0,
        ] {
            let mut buf = Vec::new();
            g6_into(&mut buf, x);
            assert_eq!(String::from_utf8(buf).expect("ascii"), g6(x), "x={x}");
        }
        assert_eq!(g6(123_456.0), "123456");
        assert_eq!(g6(0.1), "0.1");
        assert_eq!(g6(2.18867e6), "2.18867e+06");
    }
}
