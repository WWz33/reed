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

use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use reed_core::count::SiteGroupCounts;
use reed_core::genotype::{GenotypePolicies, MultiallelicPolicy, PlaneBlock, PloidyPolicy, PolicyCounters};
use reed_core::math::hwe::hwe_exact_p;
use reed_core::sample::SampleMask;
use reed_core::ChunkId;
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{run, ChunkJob, OrderedSink};
use reed_io::vcf::genotype_scan::{INT_MISSING, GenotypeScanner};
use reed_io::vcf::line::LineRef;
use reed_io::vcf::site::SiteView;
use reed_stats::summary::report::QcReportAcc;
use reed_stats::IndividualFAcc;

/// Which reports to write.
#[derive(Clone, Debug, Default)]
pub struct ReportRequest {
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
    /// Output prefix (`-o PREFIX`).
    pub prefix: PathBuf,
}

impl ReportRequest {
    /// Whether any report is requested.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        !(self.freq || self.hardy || self.het || self.missing || self.depth)
    }
}

/// Files written and records seen.
#[derive(Clone, Debug)]
pub struct ReportsOutcome {
    /// Paths written, in the order of the flags.
    pub files: Vec<PathBuf>,
    /// Records scanned.
    pub n_records: u64,
    /// The merged QC summary when `summary` was requested (`render_bcftools_stats`
    /// turns it into `PREFIX.stats`).
    pub qc: Option<QcReportAcc>,
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
}

impl ReportPlan {
    fn from_request(request: &ReportRequest) -> Self {
        Self {
            freq: request.freq,
            hardy: request.hardy,
            vmiss: request.missing,
            ldepth: request.depth,
            smiss: request.missing,
            het: request.het,
            depth: request.depth,
            summary: request.summary,
        }
    }

    /// Whether any per-site row is printed.
    const fn sites(&self) -> bool {
        self.freq || self.hardy || self.vmiss || self.ldepth
    }
}

/// Worker → sink partial: rendered rows plus the per-sample accumulators.
#[derive(Default)]
struct ReportPartial {
    /// The QC summary pieces this chunk accumulated (`--summary`).
    qc: Option<QcReportAcc>,
    freq: Vec<u8>,
    hardy: Vec<u8>,
    vmiss: Vec<u8>,
    ldepth: Vec<u8>,
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
}

/// Reused per-worker buffers.
struct ReportScratch {
    scanner: GenotypeScanner,
    /// Genotype planes: only when `--het` needs the per-sample moments.
    planes: Option<PlaneBlock>,
    counters: PolicyCounters,
    mask: SampleMask,
    /// The QC summary accumulator (`--summary`, doc_v4 O-1).
    qc: Option<crate::stats::writers::qc_report::QcAccumulator>,
}

impl ReportScratch {
    fn new(
        n_samples: usize,
        policies: GenotypePolicies,
        het: bool,
        depth: bool,
        summary: bool,
    ) -> Self {
        Self {
            // DP is decoded only for `--depth` (O2).
            scanner: if depth {
                GenotypeScanner::new(n_samples, &[b"DP"], &[], policies)
            } else {
                GenotypeScanner::new(n_samples, &[], &[], policies)
            },
            planes: het.then(|| PlaneBlock::new(1, n_samples, false)),
            counters: PolicyCounters::default(),
            mask: SampleMask::all(n_samples),
            qc: summary.then(|| {
                crate::stats::writers::qc_report::QcAccumulator::new(n_samples)
            }),
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
            )
        });
        let mut partial = self.new_partial();
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
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
    }

    /// Allocates the partial's per-sample accumulators for this request.
    fn new_partial(&self) -> ReportPartial {
        ReportPartial {
            missing: if self.plan.smiss { vec![0; self.n_samples] } else { Vec::new() },
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
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::job("record line has fewer than 9 columns"));
        };
        let site = SiteView::parse(&lr).map_err(EngineError::from)?;
        let pos1 = site.pos1;
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
            let dp_column = out
                .int_columns
                .first()
                .map(Vec::as_slice)
                .unwrap_or(&[]);
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
            self.render_site(&lr, &site, pos1, &counts, n_missing, site_depth, partial);
        }
        // The QC summary sections read the same record in the same pass
        // (doc_v4 O-1); its scanner has its own policies (haploid as half
        // call, per-ALT tracking), so it decodes independently.
        if let Some(qc) = &mut scratch.qc {
            qc.push_line(line)?;
        }
        partial.n_records += 1;
        Ok(())
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
            let _ = write!(out, "\t{}", an); out.push(b'\n');
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
                g6_into(out, hwe_exact_p(called, mac, het));
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
    }
}

/// Ordered consumer: appends the worker-rendered rows, folds per-sample
/// vectors.
struct ReportSink {
    request: ReportRequest,
    plan: ReportPlan,
    n_samples: usize,
    sample_names: Vec<String>,
    afreq: Option<BufWriter<std::fs::File>>,
    hardy: Option<BufWriter<std::fs::File>>,
    vmiss: Option<BufWriter<std::fs::File>>,
    ldepth: Option<BufWriter<std::fs::File>>,
    missing: Vec<u64>,
    depth_sums: Vec<(f64, u64)>,
    individual: IndividualFAcc,
    n_records: u64,
    /// The `--summary` accumulator (`None` unless it shares the pass).
    qc: Option<QcReportAcc>,
    out: Arc<Mutex<Option<ReportsOutcome>>>,
}

impl ReportSink {
    fn open(prefix: &Path, suffix: &str) -> Result<BufWriter<std::fs::File>, EngineError> {
        let mut path = prefix.as_os_str().to_os_string();
        path.push(suffix);
        let file = std::fs::File::create(PathBuf::from(&path))
            .map_err(EngineError::Io)?;
        Ok(BufWriter::new(file))
    }
}

impl OrderedSink<ReportPartial> for ReportSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut ReportPartial) -> Result<(), EngineError> {
        self.n_records += partial.n_records;
        let append = |writer: &mut Option<BufWriter<std::fs::File>>, bytes: &[u8]| {
            if let Some(w) = writer {
                w.write_all(bytes).map_err(EngineError::Io)?;
            }
            Ok::<(), EngineError>(())
        };
        append(&mut self.afreq, &partial.freq)?;
        append(&mut self.hardy, &partial.hardy)?;
        append(&mut self.vmiss, &partial.vmiss)?;
        append(&mut self.ldepth, &partial.ldepth)?;
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
        Ok(())
    }

    fn finish(mut self) -> Result<u64, EngineError> {
        let mut files = Vec::new();
        let path = |suffix: &str| {
            let mut p = self.request.prefix.as_os_str().to_os_string();
            p.push(suffix);
            PathBuf::from(p)
        };
        if let Some(mut w) = self.afreq.take() {
            w.flush().map_err(EngineError::Io)?;
            files.push(path(".afreq"));
        }
        if let Some(mut w) = self.hardy.take() {
            w.flush().map_err(EngineError::Io)?;
            files.push(path(".hardy"));
        }
        if let Some(mut w) = self.vmiss.take() {
            w.flush().map_err(EngineError::Io)?;
            files.push(path(".vmiss"));
        }
        if let Some(mut w) = self.ldepth.take() {
            w.flush().map_err(EngineError::Io)?;
            files.push(path(".ldepth.mean"));
        }
        if self.request.depth {
            let mut w = ReportSink::open(&self.request.prefix, ".idepth")?;
            writeln!(w, "INDV\tN_SITES\tMEAN_DEPTH").map_err(EngineError::Io)?;
            for (i, &(sum, n)) in self.depth_sums.iter().enumerate() {
                let name = self.sample_names.get(i).map_or("?", String::as_str);
                let mean = if n > 0 { g6(sum / n as f64) } else { "NA".to_string() };
                writeln!(w, "{name}\t{n}\t{mean}").map_err(EngineError::Io)?;
            }
            w.flush().map_err(EngineError::Io)?;
            files.push(path(".idepth"));
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
            w.flush().map_err(EngineError::Io)?;
            files.push(path(".het"));
        }
        if self.request.missing {
            // plink2's .smiss OBS_CT is the number of variants (the sample's
            // opportunity count), not the sample count.
            let obs = self.n_records;
            let mut w = ReportSink::open(&self.request.prefix, ".smiss")?;
            writeln!(w, "#IID\tMISSING_CT\tOBS_CT\tF_MISS").map_err(EngineError::Io)?;
            for (i, &missing) in self.missing.iter().enumerate() {
                let f_miss = if obs > 0 { g6(missing as f64 / obs as f64) } else { "NA".into() };
                let name = self.sample_names.get(i).map_or("?", String::as_str);
                writeln!(w, "{name}\t{missing}\t{obs}\t{f_miss}")
                    .map_err(EngineError::Io)?;
            }
            w.flush().map_err(EngineError::Io)?;
            files.push(path(".smiss"));
        }
        let n_records = self.n_records;
        *self.out.lock().expect("invariant: reports slot not poisoned") = Some(ReportsOutcome {
            files,
            n_records,
            qc: self.qc.take(),
            // The engine metrics come from the run, not the sink; the
            // driver fills this after `run` returns.
            metrics: reed_engine::metrics::StageMetrics::default(),
        });
        Ok(n_records)
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
    if !(-4..6).contains(&exponent) {
        // Rust prints `2.18867e6` / `1.23456e-3`; printf %g prints
        // `2.18867e+06` / `1.23456e-03`. The exponential form is rare, so a
        // small temporary is fine here.
        let raw = format!("{x:.5e}");
        match raw.split_once('e') {
            Some((mantissa, exp)) => {
                let exp: i32 = exp.parse().unwrap_or(0);
                let mantissa = trim_trailing_zeros(mantissa);
                let _ = write!(
                    out,
                    "{}e{}{:02}",
                    mantissa,
                    if exp < 0 { '-' } else { '+' },
                    exp.abs()
                );
            }
            None => out.extend_from_slice(raw.as_bytes()),
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
pub fn run_reports(
    input: Option<&Path>,
    request: &ReportRequest,
    budget: &ThreadBudget,
) -> Result<ReportsOutcome, EngineError> {
    let opened = crate::input::open(input).map_err(EngineError::from)?;
    run_reports_opened(opened, request, budget)
}

/// Runs the requested reports over an already-opened input (the CLI opens
/// once so stdin streams; doc_v3 03-performance.md O5).
///
/// # Errors
///
/// As [`run_reports`].
pub fn run_reports_opened(
    opened: crate::input::OpenedInput,
    request: &ReportRequest,
    budget: &ThreadBudget,
) -> Result<ReportsOutcome, EngineError> {
    let n_samples = opened.facts.n_samples();
    let plan = ReportPlan::from_request(request);
    let job = ReportJob {
        n_samples,
        plan,
        policies: GenotypePolicies {
            ploidy: PloidyPolicy::HaploidAsHalfCall,
            multiallelic: MultiallelicPolicy::CollapseAlt,
            ..GenotypePolicies::default()
        },
    };
    let mut afreq = request.freq.then(|| ReportSink::open(&request.prefix, ".afreq")).transpose()?;
    if let Some(w) = &mut afreq {
        writeln!(w, "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT")
            .map_err(EngineError::Io)?;
    }
    let mut hardy = request.hardy.then(|| ReportSink::open(&request.prefix, ".hardy")).transpose()?;
    if let Some(w) = &mut hardy {
        writeln!(w, "#CHROM\tID\tA1\tAX\tHOM_A1_CT\tHET_A1_CT\tTWO_AX_CT\tO(HET_A1)\tE(HET_A1)\tP")
            .map_err(EngineError::Io)?;
    }
    let mut vmiss = request.missing.then(|| ReportSink::open(&request.prefix, ".vmiss")).transpose()?;
    if let Some(w) = &mut vmiss {
        writeln!(w, "#CHROM\tID\tMISSING_CT\tOBS_CT\tF_MISS")
            .map_err(EngineError::Io)?;
    }
    let mut ldepth = request.depth.then(|| ReportSink::open(&request.prefix, ".ldepth.mean")).transpose()?;
    if let Some(w) = &mut ldepth {
        writeln!(w, "CHROM\tPOS\tMEAN_DEPTH\tVAR_DEPTH")
            .map_err(EngineError::Io)?;
    }
    let abort = reed_engine::abort::abort();
    let job = crate::bcf_input::InputJob::new(job, &opened, &abort);
    let crate::input::OpenedInput { mut source, facts, .. } = opened;
    let sink = ReportSink {
        request: request.clone(),
        plan,
        n_samples,
        sample_names: facts.sample_names.clone(),
        afreq,
        hardy,
        vmiss,
        ldepth,
        missing: vec![0; n_samples],
        depth_sums: vec![(0.0, 0); n_samples],
        individual: IndividualFAcc::new(n_samples),
        n_records: 0,
        qc: request.summary.then(|| {
            reed_stats::summary::report::QcReportAcc::new(
                u32::try_from(n_samples).unwrap_or(u32::MAX),
            )
        }),
        out: Arc::new(Mutex::new(None)),
    };
    let out = Arc::clone(&sink.out);
    let report = run(source.as_mut(), &job, sink, budget, &abort)?;
    job.finish(report.n_chunks, true)?;
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
            0.0, 1.0, 0.5, 0.3333333333, 1e-5, 1.23456789, 123456.0, 1234567.0, 0.000123456,
            2.18867e6, 1.23456e-3, f64::NAN, 1.0 / 3.0,
        ] {
            let mut buf = Vec::new();
            g6_into(&mut buf, x);
            assert_eq!(String::from_utf8(buf).expect("ascii"), g6(x), "x={x}");
        }
        assert_eq!(g6(123456.0), "123456");
        assert_eq!(g6(0.1), "0.1");
        assert_eq!(g6(2.18867e6), "2.18867e+06");
    }
}
