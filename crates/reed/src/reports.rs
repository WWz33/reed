//! Per-site and per-sample reports (02-cli.md §5.4): `--freq`, `--hardy`,
//! `--het`, `--missing`, written in plink2's column layouts.
//!
//! Rows stream out of the ordered sink as chunks arrive, so memory stays
//! bounded; the per-sample vectors fold at the end. Genotypes are decoded by
//! the shared `SampleScanner` and counted by the shared counting kernel
//! (`count_site_group`) — no report-specific genotype parser.

use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use reed_core::count::{count_site_group, SiteGroupCounts};
use reed_core::genotype::{GenotypePolicies, MultiallelicPolicy, PlaneBlock, PloidyPolicy, PolicyCounters};
use reed_core::math::hwe::hwe_exact_p;
use reed_core::sample::SampleMask;
use reed_core::ChunkId;
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{run, ChunkJob, OrderedSink};
use reed_io::vcf::line::LineRef;
use reed_io::vcf::genotype_scan::{INT_MISSING, GenotypeScanner};
use reed_io::vcf::site::SiteView;
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
}

/// One site's report row (chrom/pos/id/ref/alt plus the counted classes).
#[derive(Clone, Debug)]
struct SiteRow {
    chrom: Vec<u8>,
    pos1: u64,
    id: Vec<u8>,
    reference: Vec<u8>,
    alt: Vec<u8>,
    counts: SiteGroupCounts,
    /// Fully missing genotypes at the site.
    n_missing: u32,
    /// `(Σ DP, Σ DP², samples contributing)` for `--depth`.
    depth: (f64, f64, u64),
}

/// Worker → sink partial.
#[derive(Default)]
struct ReportPartial {
    sites: Vec<SiteRow>,
    /// Per sample: fully missing genotype count.
    missing: Vec<u64>,
    /// Per sample: plink `--het` moment components.
    individual: Option<IndividualFAcc>,
    /// Per sample: `(Σ DP, sites with a called genotype and a DP)`.
    depth_sums: Vec<(f64, u64)>,
    n_records: u64,
}

/// The scan job: one record → one site row + per-sample updates.
struct ReportJob {
    n_samples: usize,
    policies: GenotypePolicies,
}

/// Reused per-worker buffers.
struct ReportScratch {
    scanner: GenotypeScanner,
    planes: PlaneBlock,
    counters: PolicyCounters,
    mask: SampleMask,
}

/// The DP FORMAT column index in the scanner's int output.
const DP_COLUMN: usize = 0;

impl ReportScratch {
    fn new(n_samples: usize, policies: GenotypePolicies) -> Self {
        Self {
            scanner: GenotypeScanner::new(n_samples, &[b"DP"], &[], policies),
            planes: PlaneBlock::new(1, n_samples, false),
            counters: PolicyCounters::default(),
            mask: SampleMask::all(n_samples),
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
        let scratch = scratch.get_or_insert_with(|| ReportScratch::new(self.n_samples, self.policies));
        let mut partial = ReportPartial {
            missing: vec![0; self.n_samples],
            individual: Some(IndividualFAcc::new(self.n_samples)),
            depth_sums: vec![(0.0, 0); self.n_samples],
            ..ReportPartial::default()
        };
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            self.push_line(line, scratch, &mut partial)?;
        }
        Ok(partial)
    }

    fn run_seam(
        &self,
        _id: ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let scratch = scratch.get_or_insert_with(|| ReportScratch::new(self.n_samples, self.policies));
        let mut partial = ReportPartial {
            missing: vec![0; self.n_samples],
            individual: Some(IndividualFAcc::new(self.n_samples)),
            depth_sums: vec![(0.0, 0); self.n_samples],
            ..ReportPartial::default()
        };
        self.push_line(line, scratch, &mut partial)?;
        Ok(partial)
    }
}

impl ReportJob {
    fn push_line(
        &self,
        line: &[u8],
        scratch: &mut ReportScratch,
        partial: &mut ReportPartial,
    ) -> Result<(), EngineError> {
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::Job("record line has fewer than 9 columns".into()));
        };
        let site = SiteView::parse(&lr).map_err(|e| EngineError::Job(e.to_string()))?;
        let pos1 = site.pos1;
        // The scanner adds its policy triggers; reset the per-record slot.
        scratch.counters = PolicyCounters::default();
        let out = scratch
            .scanner
            .scan(&lr, &site, &mut scratch.counters)
            .map_err(|e| EngineError::Job(e.to_string()))?;
        scratch.planes.push_codes(0, out.codes);

        let mut counts = SiteGroupCounts::default();
        count_site_group(
            &scratch.planes.site_planes(0),
            scratch.mask.words(),
            u32::try_from(self.n_samples).expect("invariant: sample count < u32::MAX"),
            &mut counts,
        );

        // Per-sample fully-missing genotypes, straight from the codes.
        let mut n_missing = 0_u32;
        for (i, &code) in out.codes.iter().enumerate() {
            if code.is_fully_missing() {
                partial.missing[i] += 1;
                n_missing += 1;
            }
        }

        // vcftools depth: a sample contributes when its genotype is called
        // and the DP value is present.
        let dp_column = out.int_columns.get(DP_COLUMN).map(Vec::as_slice).unwrap_or(&[]);
        let mut site_dp_sum = 0.0_f64;
        let mut site_dp_sq = 0.0_f64;
        let mut site_dp_n = 0_u64;
        for (i, &code) in out.codes.iter().enumerate() {
            let dp = dp_column.get(i).copied().unwrap_or(INT_MISSING);
            if code.is_fully_missing() || dp == INT_MISSING {
                continue;
            }
            let dp = f64::from(dp);
            partial.depth_sums[i].0 += dp;
            partial.depth_sums[i].1 += 1;
            site_dp_sum += dp;
            site_dp_sq += dp * dp;
            site_dp_n += 1;
        }
        let site_depth = (site_dp_sum, site_dp_sq, site_dp_n);

        // plink `--het` moment components at biallelic sites.
        let n_called_alleles = counts.n_called_alleles(
            u32::try_from(self.n_samples).expect("invariant: sample count < u32::MAX"),
        );
        if n_called_alleles > 0 {
            let p_alt = f64::from(counts.n_alt_alleles) / f64::from(n_called_alleles);
            if let Some(individual) = &mut partial.individual {
                individual.push_site(&scratch.planes, 0, &scratch.mask, p_alt);
            }
        }

        partial.sites.push(SiteRow {
            chrom: lr.chrom().to_vec(),
            pos1,
            id: lr.id().to_vec(),
            reference: lr.reference().to_vec(),
            alt: lr.alts().to_vec(),
            counts,
            n_missing,
            depth: site_depth,
        });
        partial.n_records += 1;
        Ok(())
    }
}

/// Ordered consumer: streams per-site rows, folds per-sample vectors.
struct ReportSink {
    request: ReportRequest,
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
    out: Arc<Mutex<Option<ReportsOutcome>>>,
}

impl ReportSink {
    fn open(prefix: &Path, suffix: &str) -> Result<BufWriter<std::fs::File>, EngineError> {
        let mut path = prefix.as_os_str().to_os_string();
        path.push(suffix);
        let file = std::fs::File::create(PathBuf::from(&path))
            .map_err(|e| EngineError::Io(e.to_string()))?;
        Ok(BufWriter::new(file))
    }

    fn write_site_rows(&mut self, sites: &[SiteRow]) -> Result<(), EngineError> {
        for row in sites {
            let n_group = u32::try_from(self.n_samples).expect("invariant: sample count < u32::MAX");
            let an = row.counts.n_called_alleles(n_group);
            let chrom = String::from_utf8_lossy(&row.chrom);
            // plink2's layout needs an identifiable row: when the input ID is
            // missing it synthesizes `chrom_pos` (same rule here).
            let id = if row.id == b"." || row.id.is_empty() {
                format!("{chrom}_{}", row.pos1)
            } else {
                String::from_utf8_lossy(&row.id).into_owned()
            };
            if let Some(w) = &mut self.afreq {
                let freq = if an > 0 { g6(f64::from(row.counts.n_alt_alleles) / f64::from(an)) } else { "NA".into() };
                writeln!(
                    w,
                    "{chrom}\t{id}\t{}\t{}\t{freq}\t{an}",
                    String::from_utf8_lossy(&row.reference),
                    String::from_utf8_lossy(&row.alt)
                )
                .map_err(|e| EngineError::Io(e.to_string()))?;
            }
            if let Some(w) = &mut self.hardy {
                let called = row.counts.n_full_called();
                // A1 is the input's REF allele, AX the ALT — matching plink2's
                // VCF-import output (verified against plink2 2.x on this data:
                // A1 follows REF regardless of which allele is minor).
                let a1 = String::from_utf8_lossy(&row.reference).into_owned();
                let ax = String::from_utf8_lossy(&row.alt).into_owned();
                let (hom_a1, two_ax) = (row.counts.n_hom_ref, row.counts.n_hom_alt);
                let het = row.counts.n_het;
                let o_het = if called > 0 { g6(f64::from(het) / f64::from(called)) } else { "NA".into() };
                let p_a1 = if an > 0 {
                    f64::from(row.counts.n_ref_alleles) / f64::from(an)
                } else {
                    f64::NAN
                };
                let e_het = if p_a1.is_nan() { "NA".into() } else { g6(2.0 * p_a1 * (1.0 - p_a1)) };
                let mac = 2 * row.counts.n_hom_ref.min(row.counts.n_hom_alt) + het;
                let p = if called == 0 { "NA".into() } else { g6(hwe_exact_p(called, mac, het)) };
                writeln!(w, "{chrom}\t{id}\t{a1}\t{ax}\t{hom_a1}\t{het}\t{two_ax}\t{o_het}\t{e_het}\t{p}")
                    .map_err(|e| EngineError::Io(e.to_string()))?;
            }
            if let Some(w) = &mut self.vmiss {
                let obs = u64::try_from(self.n_samples).unwrap_or(u64::MAX);
                let f_miss = if obs > 0 { g6(f64::from(row.n_missing) / obs as f64) } else { "NA".into() };
                writeln!(w, "{chrom}\t{id}\t{}\t{obs}\t{f_miss}", row.n_missing)
                    .map_err(|e| EngineError::Io(e.to_string()))?;
            }
            if let Some(w) = &mut self.ldepth {
                let (sum, sq, n) = row.depth;
                let (mean, var) = if n > 0 {
                    let mean = sum / n as f64;
                    // Sample variance (n − 1), vcftools' VAR_DEPTH.
                    let var = if n > 1 { (sq - sum * sum / n as f64) / (n - 1) as f64 } else { 0.0 };
                    (g6(mean), g6(var.max(0.0)))
                } else {
                    ("NA".to_string(), "NA".to_string())
                };
                writeln!(w, "{chrom}\t{}\t{mean}\t{var}", row.pos1)
                    .map_err(|e| EngineError::Io(e.to_string()))?;
            }
        }
        Ok(())
    }
}

impl OrderedSink<ReportPartial> for ReportSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut ReportPartial) -> Result<(), EngineError> {
        self.n_records += partial.n_records;
        if self.request.freq || self.request.hardy || self.request.missing || self.request.depth {
            self.write_site_rows(&partial.sites)?;
        }
        if self.request.depth {
            for (dst, src) in self.depth_sums.iter_mut().zip(&partial.depth_sums) {
                dst.0 += src.0;
                dst.1 += src.1;
            }
        }
        if self.request.missing {
            for (dst, src) in self.missing.iter_mut().zip(&partial.missing) {
                *dst += src;
            }
        }
        if let Some(other) = partial.individual.take() {
            self.individual.merge(&other);
        }
        Ok(())
    }

    fn finish(mut self) -> Result<u64, EngineError> {
        let mut files = Vec::new();
        let mut path = |suffix: &str| {
            let mut p = self.request.prefix.as_os_str().to_os_string();
            p.push(suffix);
            PathBuf::from(p)
        };
        if let Some(mut w) = self.afreq.take() {
            w.flush().map_err(|e| EngineError::Io(e.to_string()))?;
            files.push(path(".afreq"));
        }
        if let Some(mut w) = self.hardy.take() {
            w.flush().map_err(|e| EngineError::Io(e.to_string()))?;
            files.push(path(".hardy"));
        }
        if let Some(mut w) = self.vmiss.take() {
            w.flush().map_err(|e| EngineError::Io(e.to_string()))?;
            files.push(path(".vmiss"));
        }
        if let Some(mut w) = self.ldepth.take() {
            w.flush().map_err(|e| EngineError::Io(e.to_string()))?;
            files.push(path(".ldepth.mean"));
        }
        if self.request.depth {
            let mut w = ReportSink::open(&self.request.prefix, ".idepth")?;
            writeln!(w, "INDV\tN_SITES\tMEAN_DEPTH").map_err(|e| EngineError::Io(e.to_string()))?;
            for (i, &(sum, n)) in self.depth_sums.iter().enumerate() {
                let name = self.sample_names.get(i).map_or("?", String::as_str);
                let mean = if n > 0 { g6(sum / n as f64) } else { "NA".to_string() };
                writeln!(w, "{name}\t{n}\t{mean}").map_err(|e| EngineError::Io(e.to_string()))?;
            }
            w.flush().map_err(|e| EngineError::Io(e.to_string()))?;
            files.push(path(".idepth"));
        }
        if self.request.het {
            let f = self.individual.f_values();
            let mut w = ReportSink::open(&self.request.prefix, ".het")?;
            writeln!(w, "#IID\tO(HOM)\tE(HOM)\tOBS_CT\tF").map_err(|e| EngineError::Io(e.to_string()))?;
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
                    .map_err(|e| EngineError::Io(e.to_string()))?;
            }
            w.flush().map_err(|e| EngineError::Io(e.to_string()))?;
            files.push(path(".het"));
        }
        if self.request.missing {
            // plink2's .smiss OBS_CT is the number of variants (the sample's
            // opportunity count), not the sample count.
            let obs = self.n_records;
            let mut w = ReportSink::open(&self.request.prefix, ".smiss")?;
            writeln!(w, "#IID\tMISSING_CT\tOBS_CT\tF_MISS").map_err(|e| EngineError::Io(e.to_string()))?;
            for (i, &missing) in self.missing.iter().enumerate() {
                let f_miss = if obs > 0 { g6(missing as f64 / obs as f64) } else { "NA".into() };
                let name = self.sample_names.get(i).map_or("?", String::as_str);
                writeln!(w, "{name}\t{missing}\t{obs}\t{f_miss}")
                    .map_err(|e| EngineError::Io(e.to_string()))?;
            }
            w.flush().map_err(|e| EngineError::Io(e.to_string()))?;
            files.push(path(".smiss"));
        }
        let n_records = self.n_records;
        *self.out.lock().expect("invariant: reports slot not poisoned") = Some(ReportsOutcome { files, n_records });
        Ok(n_records)
    }
}

/// `%g` with 6 significant digits, matching plink2's print precision.
fn g6(x: f64) -> String {
    if !x.is_finite() {
        return "NA".to_string();
    }
    if x == 0.0 {
        return "0".to_string();
    }
    let exponent = x.abs().log10().floor() as i32;
    let s = if !(-4..6).contains(&exponent) {
        // Rust prints `2.18867e6` / `1.23456e-3`; printf %g prints
        // `2.18867e+06` / `1.23456e-03`.
        let raw = format!("{x:.5e}");
        match raw.split_once('e') {
            Some((mantissa, exp)) => {
                let exp: i32 = exp.parse().unwrap_or(0);
                format!(
                    "{}e{}{:02}",
                    trim_trailing_zeros(mantissa),
                    if exp < 0 { '-' } else { '+' },
                    exp.abs()
                )
            }
            None => raw,
        }
    } else {
        let decimals = (5 - exponent).max(0) as usize;
        format!("{x:.decimals$}")
    };
    trim_trailing_zeros(&s)
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
    let (mut source, facts) = crate::input::open(input).map_err(EngineError::from)?;
    let n_samples = facts.n_samples();
    let job = ReportJob {
        n_samples,
        policies: GenotypePolicies {
            ploidy: PloidyPolicy::HaploidAsHalfCall,
            multiallelic: MultiallelicPolicy::CollapseAlt,
            ..GenotypePolicies::default()
        },
    };
    let mut afreq = request.freq.then(|| ReportSink::open(&request.prefix, ".afreq")).transpose()?;
    if let Some(w) = &mut afreq {
        writeln!(w, "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT")
            .map_err(|e| EngineError::Io(e.to_string()))?;
    }
    let mut hardy = request.hardy.then(|| ReportSink::open(&request.prefix, ".hardy")).transpose()?;
    if let Some(w) = &mut hardy {
        writeln!(w, "#CHROM\tID\tA1\tAX\tHOM_A1_CT\tHET_A1_CT\tTWO_AX_CT\tO(HET_A1)\tE(HET_A1)\tP")
            .map_err(|e| EngineError::Io(e.to_string()))?;
    }
    let mut vmiss = request.missing.then(|| ReportSink::open(&request.prefix, ".vmiss")).transpose()?;
    if let Some(w) = &mut vmiss {
        writeln!(w, "#CHROM\tID\tMISSING_CT\tOBS_CT\tF_MISS")
            .map_err(|e| EngineError::Io(e.to_string()))?;
    }
    let mut ldepth = request.depth.then(|| ReportSink::open(&request.prefix, ".ldepth.mean")).transpose()?;
    if let Some(w) = &mut ldepth {
        writeln!(w, "CHROM\tPOS\tMEAN_DEPTH\tVAR_DEPTH")
            .map_err(|e| EngineError::Io(e.to_string()))?;
    }
    let sink = ReportSink {
        request: request.clone(),
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
        out: Arc::new(Mutex::new(None)),
    };
    let out = Arc::clone(&sink.out);
    let abort = reed_engine::abort::abort();
    run(source.as_mut(), &job, sink, budget, &abort)?;
    Ok(out
        .lock()
        .expect("invariant: reports slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the outcome"))
}
