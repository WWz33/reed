//! `run_stats`: the L4 statistics driver (`docs/spec/architecture.md` §8).
//! Chunk workers count genotypes straight from the VCF text into per-bin
//! additive accumulators ([`reed_stats::window::ContigBins`]); the ordered
//! sink merges chunk partials in input order, emits each window as soon as
//! the watermark proves no future site can reach it (§8.3), and folds
//! contig/genome-level totals alongside (§8.7). The TSV output is
//! byte-identical across thread counts (§7.6: chunk boundaries come from
//! the input only).
//!
//! Policy notes: multiallelic records are skipped (counted separately) —
//! the fast GT path is biallelic-only, keyed off the `ALT` column.
//! Unassigned samples still feed the site-wide variant flag but no group's
//! counts. Output levels: window rows use 1-based `[start, end]` spans;
//! the contig level spans the contig's emitted windows; the genome level
//! is a single `genome` row whose `end` is the summed span in bp (not a
//! coordinate).

use std::collections::HashMap;
use std::fmt::Write as _;
use std::path::Path;
use std::sync::Arc;
use std::time::Duration;

use reed_core::ChunkId;
use reed_core::count::SiteGroupCounts;
use reed_core::count::code_has_alt;

use reed_core::genotype::{
    GenotypePolicies, GtCode, MultiallelicPolicy, PloidyPolicy, PolicyCounters,
};
use reed_core::grid::WindowGrid;
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::{ChunkJob, OrderedSink, run};
use reed_io::bcf::BcfDict;
use reed_io::vcf::genotype_scan::GenotypeScanner;
use reed_io::vcf::line::LineRef;
use reed_io::vcf::site::SiteView;
pub use reed_stats::fst::FstKind;
use reed_stats::jackknife::BpBlockJackknife;
use reed_stats::ld::LdDecayAcc;
use reed_stats::patterson::{PattersonKind, f2_site, f3_site, f4_site};
pub use reed_stats::pbs::PbsKind;
pub use reed_stats::request::StatKind;
use reed_stats::request::{StatAcc, StatKey, build_accs};
use reed_stats::sfs::{ProjectedSfsAcc, ProjectionCache};
use reed_stats::window::{BpBinning, ContigBins};

/// The full per-group/per-pair statistic set over `n_groups` groups: π and
/// dxy in both conventions, Fst in all three estimators, PBS in both
/// flavours (three groups required), and `SFS`/`θ_W`/`Tajima's D` per group —
/// the M4 acceptance set, in output-column order.
#[must_use]
pub fn default_kinds(n_groups: usize) -> Vec<StatKind> {
    let mut kinds = Vec::new();
    for group in 0..n_groups {
        kinds.push(StatKind::Pi {
            group,
            per_bp: false,
        });
        kinds.push(StatKind::Pi {
            group,
            per_bp: true,
        });
    }
    for a in 0..n_groups {
        for b in (a + 1)..n_groups {
            kinds.push(StatKind::Dxy {
                a,
                b,
                per_bp: false,
            });
            kinds.push(StatKind::Dxy { a, b, per_bp: true });
            for kind in [FstKind::Hudson, FstKind::WeirCockerham84] {
                kinds.push(StatKind::Fst { a, b, kind });
            }
        }
    }
    if n_groups >= 3 {
        kinds.push(StatKind::Pbs {
            a: 0,
            b: 1,
            c: 2,
            kind: PbsKind::Yi2010,
        });
        kinds.push(StatKind::Pbs {
            a: 0,
            b: 1,
            c: 2,
            kind: PbsKind::PbsN1,
        });
    }
    for group in 0..n_groups {
        kinds.push(StatKind::Sfs {
            group,
            folded: false,
        });
        kinds.push(StatKind::Sfs {
            group,
            folded: true,
        });
        kinds.push(StatKind::ThetaW {
            group,
            per_bp: false,
        });
        kinds.push(StatKind::ThetaW {
            group,
            per_bp: true,
        });
        kinds.push(StatKind::TajimaD {
            group,
            min_sites: 0,
        });
    }
    kinds
}

/// The window-mode statistic set: π per group and/or dxy per pair, in
/// output column order (02-cli.md §5.2).
// The booleans are the mode's independent switches; they map 1:1 to the
// user-facing flags, and a struct would not carry more meaning.
#[allow(clippy::fn_params_excessive_bools)]
#[must_use]
pub fn window_kinds(
    n_groups: usize,
    pi: bool,
    dxy: bool,
    theta: bool,
    tajima: bool,
    fst: Option<FstKind>,
    pbs: Option<PbsKind>,
) -> Vec<StatKind> {
    let mut kinds = Vec::new();
    if pi {
        for group in 0..n_groups {
            kinds.push(StatKind::Pi {
                group,
                per_bp: false,
            });
        }
    }
    if dxy {
        for a in 0..n_groups {
            for b in (a + 1)..n_groups {
                kinds.push(StatKind::Dxy {
                    a,
                    b,
                    per_bp: false,
                });
            }
        }
    }
    if theta {
        for group in 0..n_groups {
            kinds.push(StatKind::ThetaW {
                group,
                per_bp: false,
            });
        }
    }
    if tajima {
        for group in 0..n_groups {
            kinds.push(StatKind::TajimaD {
                group,
                min_sites: 0,
            });
        }
    }
    if let Some(kind) = fst {
        for a in 0..n_groups {
            for b in (a + 1)..n_groups {
                kinds.push(StatKind::Fst { a, b, kind });
            }
        }
    }
    if let Some(kind) = pbs {
        // One triple: group 0 is the branch of interest (Yi et al. 2010).
        if n_groups >= 3 {
            kinds.push(StatKind::Pbs {
                a: 0,
                b: 1,
                c: 2,
                kind,
            });
        }
    }
    kinds
}

/// Which output levels to emit (§8.7). Window rows always carry the same
/// columns; contig rows aggregate each contig's bins, the genome row
/// aggregates everything — all three are plain accumulator merges, so PBS
/// stays re-aggregable (its accumulator keeps per-pair num/den).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Levels {
    /// One row per window (the default).
    pub window: bool,
    /// One row per contig, emitted when the contig closes.
    pub contig: bool,
    /// One `genome` row, emitted at the end.
    pub genome: bool,
}

impl Levels {
    /// Parses a comma-separated `window,contig,genome` list.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] on an unknown level token.
    pub fn parse(text: &str) -> Result<Self, EngineError> {
        let mut levels = Self::default();
        for token in text.split(',').map(str::trim).filter(|t| !t.is_empty()) {
            match token {
                "window" => levels.window = true,
                "contig" => levels.contig = true,
                "genome" => levels.genome = true,
                other => {
                    return Err(EngineError::job(format!(
                        "unknown level {other:?} (expected window, contig or genome)"
                    )));
                }
            }
        }
        Ok(levels)
    }
}

/// How the input's invariant sites are treated (04-algorithms.md §3,
/// 02-cli.md §5.3). The reference rule: an all-sites input is the default, and
/// anything else must be declared.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InvariantMode {
    /// The input must contain invariant sites (the reference check); a file with
    /// only variant sites fails at the end of the run.
    Require,
    /// `--variant-only`: the file has only variant sites; absent window
    /// positions count as monomorphic with all samples called, so π/dxy
    /// denominators are completed by `(L − |S|) · pairs`.
    VariantOnly,
    /// `--bypass-invariant-check`: compute over the sites the file contains
    /// (overestimates π and dxy; the caller warns).
    Bypass,
}

/// Window-layout output (02-cli.md §5.4): π/dxy files with the
/// integer audit columns.
#[derive(Clone, Copy, Debug)]
pub struct WindowLayoutOptions {
    /// Invariant-site policy.
    pub invariant: InvariantMode,
}

/// One window's row, structured for the writers.
#[derive(Clone, Debug)]
pub struct WindowRow {
    /// Contig name.
    pub contig: String,
    /// 1-based inclusive start.
    pub start1: u64,
    /// 1-based inclusive end.
    pub end0: u64,
    /// Window span in bp (the per-bp denominator).
    pub bp: u64,
    /// Accumulators in `kinds` order.
    pub accs: Vec<StatAcc>,
}

/// One `run_stats` request: statistics, group assignment, window grid.
#[derive(Clone, Debug, Default)]
pub struct StatsRequest {
    /// Statistics to compute, in output-column order.
    pub kinds: Vec<StatKind>,
    /// Sample index → group index; `None` samples feed only the site-wide
    /// variant flag. Length must equal the header's sample count.
    pub gt_groups: Vec<Option<u32>>,
    /// Window width in bp.
    pub window_size: u64,
    /// Window stride in bp (`step ≤ size` for non-overlapping windows).
    pub step: u64,
    /// Output levels (default: window rows only).
    pub levels: Levels,
    /// Skip windows no site fell into (§8.5 `--skip-empty`).
    pub skip_empty: bool,
    /// Extra contig lengths (`--contig-lengths`, .fai-style), merged over
    /// the header's `##contig` declarations. Trailing-window clipping needs
    /// a length for every contig with records and is never fabricated
    /// (01-review D7).
    pub contig_lengths: HashMap<Vec<u8>, u64>,
    /// `Some` = window-layout output mode (structured rows, invariant-site
    /// policy); `None` = the wide debug TSV.
    pub layout: Option<WindowLayoutOptions>,
    /// Indexed regions (`-r/-R`): the reader seeks to the mapped range and
    /// records outside the regions are dropped.
    pub regions: Vec<reed_io::chunk::region::Region>,
    /// Streamed target regions (`-t/-T`, parsed): kept when inside any
    /// include and no exclude region.
    pub targets: Vec<reed_filter::spec::region::RegionTriple>,
    /// Excluded target regions (`^` entries of `-t/-T`).
    pub targets_exclude: Vec<reed_filter::spec::region::RegionTriple>,
    /// Per-sample keep mask (`-s/-S/--keep/--remove`); `None` = all kept.
    pub keep_samples: Option<Vec<bool>>,
    /// The compiled filter plan (`-v/-m/-f/-i/-e/--min-dp/--max-missing/--maf/...`,
    /// 04-cleanup §2.1): applied to every record in the same pass by the
    /// filtered job. `None` = no site-level filtering.
    pub filter: Option<Arc<reed_filter::plan::FilterPlan>>,
    /// Projected SFS, genome-wide (docs/rounds/v5/05-algorithms.md §3.1);
    /// `None` = not requested. Not a window statistic.
    pub sfs: Option<SfsRequest>,
    /// Patterson f-statistics, genome-wide with a block jackknife
    /// (docs/rounds/v5/05-algorithms.md §3.3); `None` = not requested.
    pub fstats: Option<FStatsRequest>,
    /// LD: `--ld-decay` / `--indep-pairwise`
    /// (docs/rounds/v5/05-algorithms.md §3.4); `None` = not requested.
    pub ld: Option<LdRequest>,
}

/// An LD request (docs/rounds/v5/05-algorithms.md §3.4).
#[derive(Clone, Debug, PartialEq)]
pub struct LdRequest {
    /// `--ld-decay MAX_BP`: pairwise r² within this span, binned by
    /// `decay_bin_bp`. `None` = no decay table.
    pub decay_max_bp: Option<u64>,
    /// Decay bin width in bp (default 1000).
    pub decay_bin_bp: u64,
    /// Write the per-pair table (`--ld-pairs`).
    pub pairs: bool,
}

/// One decay bin: `(bin_start_bp, bin_end_bp, n_pairs, mean_r2)`.
pub type DecayBin = (u64, u64, u64, Option<f64>);

/// The LD output (`--ld-decay` / `--ld-pairs`).
///
/// The cohort curve pools every stored sample (one r² per pair); the
/// per-population curves restrict each pair to that population's samples.
/// `.ld.decay` writes the cohort under `pop = all` and then, when more than
/// one population is named, one block per population.
#[derive(Clone, Debug, Default)]
pub struct LdReport {
    /// The cohort's decay bins.
    pub decay: Vec<DecayBin>,
    /// One decay curve per named group, in group order (`--pop`): bins in
    /// the cohort's shape. Empty without a populations file.
    pub group_decay: Vec<Vec<DecayBin>>,
    /// Per-pair rows: `(contig, pos1_a, pos1_b, r2)`.
    pub pairs: Vec<(String, u64, u64, f64)>,
}

/// The validated LD plan shared by the job and the sink.
#[derive(Clone, Debug)]
struct LdPlan {
    decay_max_bp: Option<u64>,
    decay_bin_bp: u64,
    pairs: bool,
}

impl LdPlan {
    /// Decay bin edges: `[0, bin), [bin, 2·bin), …` up to `decay_max_bp`.
    fn decay_edges(&self) -> Vec<u64> {
        let bin = self.decay_bin_bp.max(1);
        let max = self.decay_max_bp.unwrap_or(bin);
        let mut edges = Vec::new();
        let mut e = 0;
        while e < max {
            edges.push(e);
            e += bin;
        }
        edges
    }
}

/// One pass0 run: a contig's consecutive surviving sites as flat site-major
/// planes (4·w words per site, `[slot0_alt, slot0_miss, slot1_alt,
/// slot1_miss]`) plus their positions.
#[derive(Clone, Debug)]
struct LdRun {
    contig: Vec<u8>,
    flat: Vec<u64>,
    pos0: Vec<u64>,
}

/// Per-chunk LD partial: the runs the worker filled. The pairs are computed
/// position-block parallel from the accumulated matrices after the stream
/// ends (docs/rounds/v6/02-models.md §2.2), so no head/tail carry exists.
#[derive(Clone, Debug, Default)]
struct LdPartial {
    runs: Vec<LdRun>,
}

/// A projected site-frequency-spectrum request (docs/rounds/v5/05-algorithms.md
/// §3.1): one genome-wide spectrum per population.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SfsRequest {
    /// Projection size in chromosomes, one entry per group (the CLI
    /// broadcasts a single `--proj` value to every group).
    pub projections: Vec<u32>,
    /// Fold the spectrum (default: unfolded — REF is the ancestral state).
    pub folded: bool,
}

/// The spectra an [`SfsRequest`] produced.
#[derive(Clone, Debug)]
pub struct SfsReport {
    /// Per group: ξ_k for k = 0..=n (folded: 0..=n/2).
    pub spectra: Vec<Vec<f64>>,
    /// Per group: the projection size in chromosomes.
    pub projections: Vec<u32>,
    /// Whether the spectra are folded.
    pub folded: bool,
    /// Per group: sites absorbed.
    pub n_sites: Vec<u64>,
    /// Per group: sites skipped for `m < n`.
    pub n_skipped_m_lt_n: Vec<u64>,
}

/// One Patterson f-statistic test (docs/rounds/v5/05-algorithms.md §3.3).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FStatTest {
    /// `f2(A, B)`.
    F2 {
        /// First group.
        a: usize,
        /// Second group.
        b: usize,
    },
    /// `f3(C; A, B)`.
    F3 {
        /// Test population (the branch of interest).
        c: usize,
        /// First source.
        a: usize,
        /// Second source.
        b: usize,
    },
    /// `f4(A, B; C, D)`.
    F4 {
        /// First group.
        a: usize,
        /// Second group.
        b: usize,
        /// Third group.
        c: usize,
        /// Fourth group.
        d: usize,
    },
    /// `D(A, B; C, D)` (the normalized f4).
    Dstat {
        /// First group.
        a: usize,
        /// Second group.
        b: usize,
        /// Third group.
        c: usize,
        /// Fourth group.
        d: usize,
    },
    /// `f4-ratio = f4(A, O; X, C) / f4(A, O; B, C)`.
    F4Ratio {
        /// Outgroup O.
        o: usize,
        /// Sister A.
        a: usize,
        /// Test X.
        x: usize,
        /// Reference B.
        b: usize,
        /// Outgroup-side C.
        c: usize,
    },
}

impl FStatTest {
    /// The stable statistic id (`f2`, `f3`, `f4`, `dstat`, `f4ratio`).
    #[must_use]
    pub const fn id(&self) -> &'static str {
        match self {
            Self::F2 { .. } => "f2",
            Self::F3 { .. } => "f3",
            Self::F4 { .. } => "f4",
            Self::Dstat { .. } => "dstat",
            Self::F4Ratio { .. } => "f4ratio",
        }
    }

    /// The populations in the order the output file prints them.
    #[must_use]
    pub fn pops(&self) -> Vec<usize> {
        match *self {
            Self::F2 { a, b } => vec![a, b],
            Self::F3 { c, a, b } => vec![c, a, b],
            Self::F4 { a, b, c, d } | Self::Dstat { a, b, c, d } => vec![a, b, c, d],
            Self::F4Ratio { o, a, x, b, c } => vec![a, o, x, b, c],
        }
    }
}

/// A Patterson f-statistics request: the tests plus the block-jackknife
/// block size in bp (default 5M).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FStatsRequest {
    /// Tests, in output order.
    pub tests: Vec<FStatTest>,
    /// Block size in bp for the delete-one-block jackknife.
    pub block_bp: u64,
}

/// One test's result.
#[derive(Clone, Debug)]
pub struct FStatRow {
    /// The test.
    pub test: FStatTest,
    /// The point estimate (ratio of sums for f3*/D/f4-ratio).
    pub estimate: f64,
    /// Block-jackknife standard error.
    pub se: f64,
    /// `estimate / se`.
    pub z: f64,
    /// Blocks contributing.
    pub n_blocks: usize,
    /// Sites contributing.
    pub n_sites: u64,
}

/// The f-statistics report.
#[derive(Clone, Debug)]
pub struct FStatsReport {
    /// One row per test, in request order.
    pub rows: Vec<FStatRow>,
}

/// Per-test accumulator state: the per-site sums and the bp jackknife.
#[derive(Clone, Debug)]
struct FStatAcc {
    test: FStatTest,
    /// Σnum, Σden (the jackknife's inputs).
    num: reed_core::math::sum::NeumaierSum,
    den: reed_core::math::sum::NeumaierSum,
    n_sites: u64,
    jk: BpBlockJackknife,
}

impl FStatAcc {
    fn new(test: FStatTest, block_bp: u64) -> Self {
        Self {
            test,
            num: reed_core::math::sum::NeumaierSum::new(),
            den: reed_core::math::sum::NeumaierSum::new(),
            n_sites: 0,
            jk: BpBlockJackknife::new(block_bp),
        }
    }

    /// The per-site `(num, den)` pair of this test at one site; `None` when
    /// the site is skipped. `den` is 1 for the additive statistics, so the
    /// ratio path and the sum path share one accumulator.
    fn site(&self, row: &[SiteGroupCounts]) -> Option<(f64, f64)> {
        match self.test {
            FStatTest::F2 { a, b } => Some((f2_site(row, a, b, PattersonKind::Corrected)?, 1.0)),
            // The raw Patterson f3 = ΣT (docs/rounds/v5/05 §3.3).
            FStatTest::F3 { c, a, b } => {
                Some((f3_site(row, c, a, b, PattersonKind::Corrected)?.0, 1.0))
            }
            FStatTest::F4 { a, b, c, d } => Some((f4_site(row, a, b, c, d)?.0, 1.0)),
            FStatTest::Dstat { a, b, c, d } => Some(f4_site(row, a, b, c, d)?),
            // f4-ratio = Σf4(A,O;X,C) / Σf4(A,O;B,C): the denominator f4's
            // numerator sum, exactly like Dstat's but with the other f4.
            FStatTest::F4Ratio { o, a, x, b, c } => {
                let (num_n, _) = f4_site(row, a, o, x, c)?;
                let (num_d, _) = f4_site(row, a, o, b, c)?;
                Some((num_n, num_d))
            }
        }
    }

    fn push(&mut self, contig: &[u8], pos0: u64, row: &[SiteGroupCounts]) {
        let Some((num, den)) = self.site(row) else {
            return;
        };
        self.num.add(num);
        self.den.add(den);
        self.n_sites += 1;
        self.jk.push(contig, pos0, num, den);
    }

    fn merge(&mut self, other: &Self) {
        self.num.merge(&other.num);
        self.den.merge(&other.den);
        self.n_sites += other.n_sites;
        self.jk.merge(&other.jk);
    }

    fn row(&self) -> Option<FStatRow> {
        // Additive statistics (f2/f3/f4) report Σnum, with the jackknife
        // deleting one block from that same sum; the ratios (D, f4-ratio)
        // report Σnum/Σden over the pair.
        let (estimate, _se) = match self.test {
            FStatTest::Dstat { .. } | FStatTest::F4Ratio { .. } => {
                let den = self.den.value();
                if den == 0.0 {
                    return None;
                }
                (self.num.value() / den, f64::NAN)
            }
            _ => (self.num.value(), f64::NAN),
        };
        let (est_jk, se_jk, n_blocks) = self.jk.jackknife().unwrap_or((estimate, f64::NAN, 0));
        let z = if se_jk.is_finite() && se_jk > 0.0 {
            est_jk / se_jk
        } else {
            f64::NAN
        };
        Some(FStatRow {
            test: self.test.clone(),
            estimate,
            se: se_jk,
            z,
            n_blocks,
            n_sites: self.n_sites,
        })
    }
}

/// The validated f-statistics plan shared by the job and the sink.
#[derive(Clone, Debug)]
struct FStatsPlan {
    tests: Arc<[FStatTest]>,
    block_bp: u64,
}

/// Per-chunk f-statistics partials.
#[derive(Clone, Debug)]
struct FStatsPartial {
    per_test: Vec<FStatAcc>,
}

impl FStatsPartial {
    fn new(plan: &FStatsPlan) -> Self {
        Self {
            per_test: plan
                .tests
                .iter()
                .map(|test| FStatAcc::new(test.clone(), plan.block_bp))
                .collect(),
        }
    }

    fn merge(&mut self, other: &Self) {
        for (dst, src) in self.per_test.iter_mut().zip(&other.per_test) {
            dst.merge(src);
        }
    }
}

/// The validated SFS plan shared by the job and the sink.
#[derive(Clone, Debug)]
struct SfsPlan {
    projections: Arc<[u32]>,
    folded: bool,
}

/// Per-chunk SFS partials: one projected spectrum per group.
#[derive(Clone, Debug)]
struct SfsPartial {
    per_group: Vec<ProjectedSfsAcc>,
}

impl SfsPartial {
    fn new(plan: &SfsPlan) -> Option<Self> {
        let mut per_group = Vec::with_capacity(plan.projections.len());
        for &n in plan.projections.iter() {
            per_group.push(ProjectedSfsAcc::new(n, plan.folded)?);
        }
        Some(Self { per_group })
    }

    fn merge(&mut self, other: &Self) {
        for (dst, src) in self.per_group.iter_mut().zip(&other.per_group) {
            dst.merge(src);
        }
    }
}

/// The run's output and counters.
#[derive(Clone, Debug)]
pub struct StatsReport {
    /// TSV: a header line, then one row per selected level.
    pub tsv: String,
    /// Structured window rows (window-layout mode only).
    pub windows: Vec<WindowRow>,
    /// Site-wide invariant sites seen in the input.
    pub n_invariant_sites: u64,
    /// Records with ALT = `.` (the all-sites check).
    pub n_no_alt_records: u64,
    /// Window-level rows emitted.
    pub n_windows: u64,
    /// Records absorbed into windows.
    pub n_records: u64,
    /// Records skipped for a multi-allelic ALT column.
    pub n_multiallelic_skipped: u64,
    /// The projected spectra when [`StatsRequest::sfs`] was set.
    pub sfs: Option<SfsReport>,
    /// The Patterson f-statistics when [`StatsRequest::fstats`] was set.
    pub fstats: Option<FStatsReport>,
    /// The LD tables when [`StatsRequest::ld`] was set.
    pub ld: Option<LdReport>,
    /// Chunks processed.
    pub n_chunks: u64,
    /// The engine's stage metrics (doc_v4 03 A4-5): reader/reducer/chain
    /// cumulative times, printed at `--verbosity 2`.
    pub metrics: reed_engine::metrics::StageMetrics,
    /// Wall time.
    pub elapsed: Duration,
}

/// Parses a `sample<TAB>population` groups file — the reference populations-file
/// convention: two columns, no header row, `#` comments allowed
/// (reed_core::GroupTable::from_tsv). A file row naming a sample absent
/// from the VCF header is an error (the reference behavior) unless `force_samples`,
/// which skips it with a warning. Group ids follow first-appearance order.
///
/// # Errors
///
/// [`EngineError::Job`] for an unreadable file, a malformed row, or an
/// unknown sample (without `force_samples`).
pub fn assign_groups(
    groups_path: &Path,
    sample_names: &[String],
    force_samples: bool,
) -> Result<(Vec<Option<u32>>, Vec<String>), EngineError> {
    let samples =
        reed_core::sample::SampleTable::from_names(sample_names.iter().map(String::as_str))
            .map_err(EngineError::from)?;
    let table = reed_core::group::GroupTable::from_tsv(
        groups_path,
        &samples,
        reed_core::group::GroupFileOptions {
            ignore_missing_samples: force_samples,
        },
    )
    .map_err(EngineError::from)?;
    let mut gt_groups: Vec<Option<u32>> = vec![None; sample_names.len()];
    let mut names: Vec<String> = Vec::with_capacity(table.n_groups());
    for (group_idx, group) in table.iter().enumerate() {
        names.push(group.name.clone());
        let group_idx =
            u32::try_from(group_idx).map_err(|_| EngineError::job("more than u32::MAX groups"))?;
        for (word_idx, &word) in group.mask.words().iter().enumerate() {
            let mut word = word;
            while word != 0 {
                let bit = word.trailing_zeros() as usize;
                gt_groups[word_idx * 64 + bit] = Some(group_idx);
                word &= word - 1;
            }
        }
    }
    Ok((gt_groups, names))
}

/// One contig's run of sites: bin accumulators plus the span it covers.
struct SegmentPartial {
    contig: Vec<u8>,
    bins: ContigBins,
    first_pos0: u64,
    last_pos0: u64,
}

/// Worker → sink partial: per-contig bin accumulators in record order.
#[derive(Default)]
pub struct StatsPartial {
    segments: Vec<SegmentPartial>,
    n_records: u64,
    n_multiallelic_skipped: u64,
    /// Site-wide invariant sites (no ALT called by any sample).
    n_invariant_sites: u64,
    /// Records whose ALT column is `.` (the all-sites marker).
    n_no_alt_records: u64,
    /// Genome-wide projected spectra (docs/rounds/v5/05-algorithms.md §3.1).
    sfs: Option<SfsPartial>,
    /// Genome-wide Patterson f-statistics.
    fstats: Option<FStatsPartial>,
    /// LD partials.
    ld: Option<LdPartial>,
}

/// The chunk job: line → `LineRef` → GT tallies → group count row → the
/// site's bin accumulators.
pub struct StatsJob {
    width: u64,
    n_groups: usize,
    gt_groups: Arc<[Option<u32>]>,
    n_samples: usize,
    seed: Arc<[StatAcc]>,
    /// Streamed target regions (`-t/-T`): kept when inside any include and
    /// no exclude; empty lists keep everything.
    targets: Arc<[reed_filter::spec::region::RegionTriple]>,
    targets_exclude: Arc<[reed_filter::spec::region::RegionTriple]>,
    /// Declared contig lengths (header + `--contig-lengths`); positions
    /// beyond the declared length are rejected (01-review D8).
    contig_len: Arc<HashMap<Vec<u8>, u64>>,
    /// The BCF-native fast path (doc_v3 02-bcf.md §2.4/§2.5): present only
    /// when the input is BCF.
    bcf: Option<BcfNative>,
    /// The genome-wide projected-SFS plan; `None` = not requested.
    sfs: Option<SfsPlan>,
    /// The genome-wide f-statistics plan; `None` = not requested.
    fstats: Option<FStatsPlan>,
    /// The LD plane layout (group-major); present only with an LD request.
    ld_layout: Option<Arc<reed_stats::SampleLayout>>,
}

/// The BCF fast path's shared state: the dictionaries, the FORMAT index of
/// `GT` (resolved once) and the byte-pair GT table.
struct BcfNative {
    dict: Arc<BcfDict>,
    gt_key: Option<u32>,
    table: reed_io::bcf::GtTable,
}

/// Worker-local LD state: the run being filled for the current contig.
#[derive(Default)]
pub struct LdScratch {
    contig: Vec<u8>,
    /// The run's site-major planes (4·w words per site).
    flat: Vec<u64>,
    /// The run's site positions.
    pos0: Vec<u64>,
}

impl LdScratch {
    /// Moves the open run into the chunk partial (a no-op on an empty run).
    fn flush(&mut self, partial: &mut StatsPartial) {
        if self.pos0.is_empty() {
            return;
        }
        let run = LdRun {
            contig: std::mem::take(&mut self.contig),
            flat: std::mem::take(&mut self.flat),
            pos0: std::mem::take(&mut self.pos0),
        };
        partial
            .ld
            .get_or_insert_with(LdPartial::default)
            .runs
            .push(run);
    }
}

/// The positional-block size of the parallel decay scan: coarse enough that
/// the per-block lookback span stays a small multiple of the block, fine
/// enough to fill every thread.
const LD_BLOCK: usize = 4096;

/// Positional-block parallel decay over one contig's matrix
/// (docs/rounds/v6/02-models.md §2.2). Each block pairs its sites with
/// every predecessor within `max_bp` — read-only, so the blocks are
/// independent; the results merge in block order, which is the serial
/// `(later site, earlier site)` order, keeping the output deterministic
/// across thread counts.
fn decay_contig(
    m: &reed_stats::ContigMatrix,
    max_bp: u64,
    edges: &[u64],
    collect_pairs: bool,
) -> reed_stats::LdDecaySet {
    use rayon::prelude::*;
    let n = m.n_sites();
    if n < 2 {
        return reed_stats::LdDecaySet::new(edges, m.layout().n_groups());
    }
    let ranges: Vec<std::ops::Range<usize>> = (0..n)
        .step_by(LD_BLOCK)
        .map(|s| s..(s + LD_BLOCK).min(n))
        .collect();
    let parts: Vec<reed_stats::LdDecaySet> = ranges
        .par_iter()
        .map(|r| reed_stats::decay_range(m, r.clone(), max_bp, edges, collect_pairs))
        .collect();
    let mut set = reed_stats::LdDecaySet::new(edges, m.layout().n_groups());
    for part in &parts {
        set.merge(part);
    }
    set
}

/// Worker scratch: the group count row, the shared record decoder and its
/// per-record policy counters.
#[derive(Default)]
pub struct StatsScratch {
    row: Vec<SiteGroupCounts>,
    scanner: Option<GenotypeScanner>,
    counters: PolicyCounters,
    /// The native BCF path's decoded GT codes (one per sample, per record).
    codes: Vec<GtCode>,
    /// Per-group projection-weight memos for the SFS path.
    sfs_caches: Vec<ProjectionCache>,
    /// Worker-local LD state.
    ld: LdScratch,
}

/// The statistics decode policies: haploid calls are half calls and ALT
/// indices above 1 collapse onto the single ALT allele (the biallelic
/// view; multi-allelic records are skipped before the scan anyway).
fn stats_policies() -> GenotypePolicies {
    GenotypePolicies {
        ploidy: PloidyPolicy::HaploidAsHalfCall,
        multiallelic: MultiallelicPolicy::CollapseAlt,
        ..GenotypePolicies::default()
    }
}

impl ChunkJob for StatsJob {
    type Partial = StatsPartial;
    type Scratch = StatsScratch;

    fn run_chunk(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        body: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let mut partial = StatsPartial::default();
        let mut current: Option<SegmentPartial> = None;
        self.reset_row(scratch);
        for line in body.split(|&b| b == b'\n') {
            if line.is_empty() {
                continue;
            }
            self.push_line(line, &mut current, &mut partial, scratch)?;
        }
        if let Some(segment) = current.take() {
            partial.segments.push(segment);
        }
        scratch.ld.flush(&mut partial);
        Ok(partial)
    }

    fn run_seam(
        &self,
        _id: reed_core::ChunkId,
        line: &[u8],
        scratch: &mut Self::Scratch,
    ) -> Result<Self::Partial, EngineError> {
        let mut partial = StatsPartial::default();
        let mut current: Option<SegmentPartial> = None;
        self.reset_row(scratch);
        self.push_line(line, &mut current, &mut partial, scratch)?;
        if let Some(segment) = current.take() {
            partial.segments.push(segment);
        }
        scratch.ld.flush(&mut partial);
        Ok(partial)
    }
}

impl StatsJob {
    fn reset_row(&self, scratch: &mut StatsScratch) {
        scratch.row.clear();
        scratch
            .row
            .resize(self.n_groups, SiteGroupCounts::default());
    }

    /// The site stage both record sources share (text and native BCF): the
    /// streamed-target predicate, the contig segment bookkeeping, the POS
    /// contract and the invariant/multiallelic accounting. Returns the
    /// 0-based position when the record continues to the genotype stage.
    fn site_prefix(
        &self,
        contig: &[u8],
        pos1: u64,
        n_alt: u32,
        is_invariant: bool,
        current: &mut Option<SegmentPartial>,
        partial: &mut StatsPartial,
    ) -> Result<Option<u64>, EngineError> {
        // Streamed targets (`-t/-T`): the region predicate gates every record
        // before it reaches an accumulator.
        if (!self.targets.is_empty() || !self.targets_exclude.is_empty())
            && !reed_filter::spec::region::region_allows(
                &self.targets,
                &self.targets_exclude,
                contig,
                pos1,
            )
        {
            return Ok(None);
        }
        match current {
            Some(segment) if segment.contig == contig => {}
            _ => {
                if let Some(segment) = current.take() {
                    partial.segments.push(segment);
                }
                *current = Some(SegmentPartial {
                    contig: contig.to_vec(),
                    bins: ContigBins::new(BpBinning::new(self.width), &self.seed),
                    first_pos0: u64::MAX,
                    last_pos0: 0,
                });
            }
        }
        let Some(pos0) = pos1.checked_sub(1) else {
            return Err(EngineError::job(format!(
                "{}:{}: POS is zero",
                String::from_utf8_lossy(contig),
                pos1
            )));
        };
        // Input contract (04-algorithms.md §6): positions ascending within a
        // contig, never beyond the declared length.
        if let Some(&len) = self.contig_len.get(contig)
            && pos1 > len
        {
            return Err(EngineError::from(
                reed_core::CoreError::PositionBeyondContigLength {
                    contig_name: String::from_utf8_lossy(contig).into_owned(),
                    pos1,
                    len_bp: len,
                },
            ));
        }
        if is_invariant {
            partial.n_no_alt_records += 1;
        }
        if n_alt > 1 {
            partial.n_multiallelic_skipped += 1;
            return Ok(None);
        }
        Ok(Some(pos0))
    }

    /// The genotype stage both record sources share: the per-group counts,
    /// the ascending-position check and the bin accumulation.
    // The arguments are distinct caller-owned scratch buffers (row, caches,
    // segment, partial, LD scratch); bundling them into a struct would only
    // move the same borrows around.
    #[allow(clippy::too_many_arguments)]
    fn push_row(
        &self,
        pos0: u64,
        codes: &[GtCode],
        row: &mut [SiteGroupCounts],
        sfs_caches: &mut Vec<ProjectionCache>,
        current: &mut Option<SegmentPartial>,
        partial: &mut StatsPartial,
        ld_scratch: &mut LdScratch,
    ) -> Result<(), EngineError> {
        for counts in row.iter_mut() {
            *counts = SiteGroupCounts::default();
        }
        let mut any_alt = false;
        for (sample_idx, &code) in codes.iter().enumerate() {
            any_alt |= code_has_alt(code);
            if let Some(group) = self.gt_groups.get(sample_idx).copied().flatten() {
                row[group as usize].add_code(code);
            }
        }
        let segment = current
            .as_mut()
            .expect("invariant: segment opened for this contig above");
        if pos0 < segment.last_pos0 {
            return Err(EngineError::from(
                reed_core::CoreError::PositionNotAscending {
                    prev_pos0: segment.last_pos0,
                    pos0,
                },
            ));
        }
        segment.bins.push_site(pos0, row, any_alt);
        segment.first_pos0 = segment.first_pos0.min(pos0);
        segment.last_pos0 = pos0;
        if !any_alt {
            partial.n_invariant_sites += 1;
        }
        self.push_sfs(row, partial, sfs_caches);
        self.push_fstats(&segment.contig.clone(), pos0, row, partial);
        self.push_ld(&segment.contig, pos0, codes, partial, ld_scratch);
        partial.n_records += 1;
        Ok(())
    }

    /// Feeds the LD path for one site (docs/rounds/v6/02-models.md §2.2):
    /// the site's planes join the worker's contig run in the group-major
    /// layout; the pairs are computed position-block parallel after the
    /// stream ends, so no carry or head/tail windows exist here.
    fn push_ld(
        &self,
        contig: &[u8],
        pos0: u64,
        codes: &[GtCode],
        partial: &mut StatsPartial,
        ld: &mut LdScratch,
    ) {
        let Some(layout) = &self.ld_layout else {
            return;
        };
        if ld.contig != contig {
            ld.flush(partial);
            ld.contig.clear();
            ld.contig.extend_from_slice(contig);
        }
        let w = layout.w();
        let base = ld.flat.len();
        ld.flat.resize(base + 4 * w, 0);
        for (sample, &code) in codes.iter().enumerate() {
            let Some(bit) = layout.bit_of(sample) else {
                continue;
            };
            let (word, mask) = (bit / 64, 1_u64 << (bit % 64));
            if code.slot_alt(0) {
                ld.flat[base + word] |= mask;
            }
            if code.slot_missing(0) {
                ld.flat[base + w + word] |= mask;
            }
            if code.slot_alt(1) {
                ld.flat[base + 2 * w + word] |= mask;
            }
            if code.slot_missing(1) {
                ld.flat[base + 3 * w + word] |= mask;
            }
        }
        ld.pos0.push(pos0);
    }

    /// Feeds the genome-wide f-statistics for one site (docs/rounds/v5/05
    /// §3.3). The segment's contig names the jackknife block key.
    fn push_fstats(
        &self,
        contig: &[u8],
        pos0: u64,
        row: &[SiteGroupCounts],
        partial: &mut StatsPartial,
    ) {
        let Some(plan) = &self.fstats else { return };
        let stats = partial
            .fstats
            .get_or_insert_with(|| FStatsPartial::new(plan));
        for acc in &mut stats.per_test {
            acc.push(contig, pos0, row);
        }
    }

    /// Feeds the genome-wide projected spectra for one site (docs/rounds/v5/
    /// 05-algorithms.md §3.1): `m < n` sites are skipped and counted.
    fn push_sfs(
        &self,
        row: &[SiteGroupCounts],
        partial: &mut StatsPartial,
        sfs_caches: &mut Vec<ProjectionCache>,
    ) {
        let Some(plan) = &self.sfs else { return };
        let sfs = partial.sfs.get_or_insert_with(|| {
            SfsPartial::new(plan).expect("invariant: projections validated at request time")
        });
        if sfs_caches.len() != sfs.per_group.len() {
            *sfs_caches = plan
                .projections
                .iter()
                .map(|&n| ProjectionCache::new(n))
                .collect();
        }
        for (g, acc) in sfs.per_group.iter_mut().enumerate() {
            let counts = &row[g];
            let m = counts.n_ref_alleles + counts.n_alt_alleles;
            if m >= acc.n() {
                let weights = sfs_caches[g].weights(m, counts.n_alt_alleles);
                acc.push_weighted(Some(&weights));
            } else {
                acc.push_weighted(None);
            }
        }
    }

    fn push_line(
        &self,
        line: &[u8],
        current: &mut Option<SegmentPartial>,
        partial: &mut StatsPartial,
        scratch: &mut StatsScratch,
    ) -> Result<(), EngineError> {
        let Some(lr) = LineRef::parse(line) else {
            return Err(EngineError::job("record line has fewer than 9 columns"));
        };
        let site = SiteView::parse(&lr).map_err(EngineError::from)?;
        let Some(pos0) = self.site_prefix(
            site.contig_name,
            site.pos1,
            site.n_alt(),
            site.is_invariant(),
            current,
            partial,
        )?
        else {
            return Ok(());
        };
        let StatsScratch {
            row,
            scanner,
            counters,
            sfs_caches,
            ld,
            ..
        } = scratch;
        *counters = PolicyCounters::default();
        let n_samples = self.n_samples;
        let scanner = scanner
            .get_or_insert_with(|| GenotypeScanner::new(n_samples, &[], &[], stats_policies()));
        let out = scanner
            .scan(&lr, &site, counters)
            .map_err(EngineError::from)?;
        self.push_row(pos0, out.codes, row, sfs_caches, current, partial, ld)
    }

    /// The record sink for a filtering wrapper (doc_v4 O-2): the wrapper
    /// already decoded and masked the genotypes, so this path parses only
    /// the fixed columns again and skips both the GT decode and the
    /// rewritten-text round-trip.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] for malformed records.
    pub fn run_records(
        &self,
        records: &[crate::filtered::ScannedRecord<'_>],
        scratch: &mut StatsScratch,
    ) -> Result<StatsPartial, EngineError> {
        let mut partial = StatsPartial::default();
        let mut current: Option<SegmentPartial> = None;
        self.reset_row(scratch);
        for rec in records {
            let Some(lr) = LineRef::parse(rec.line) else {
                return Err(EngineError::job("record line has fewer than 9 columns"));
            };
            let site = SiteView::parse(&lr).map_err(EngineError::from)?;
            let Some(pos0) = self.site_prefix(
                site.contig_name,
                site.pos1,
                site.n_alt(),
                site.is_invariant(),
                &mut current,
                &mut partial,
            )?
            else {
                continue;
            };
            self.push_row(
                pos0,
                rec.codes,
                &mut scratch.row,
                &mut scratch.sfs_caches,
                &mut current,
                &mut partial,
                &mut scratch.ld,
            )?;
        }
        if let Some(segment) = current.take() {
            partial.segments.push(segment);
        }
        scratch.ld.flush(&mut partial);
        Ok(partial)
    }

    /// The native BCF record sink (doc_v3 02-bcf.md §2.4/§2.5): each framed
    /// record is decoded straight into the shared site/row stages — the
    /// shared block past its header and every unplanned FORMAT field are
    /// skipped by descriptor length, so no text is ever produced.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] for malformed records and policy rejections.
    pub fn run_bcf_records(
        &self,
        _chunk: &reed_io::chunk::ChunkBytes,
        records: &[&[u8]],
        scratch: &mut StatsScratch,
    ) -> Result<StatsPartial, EngineError> {
        let Some(bcf) = &self.bcf else {
            return Err(EngineError::job(
                "internal: BCF records without a BCF input",
            ));
        };
        let mut partial = StatsPartial::default();
        let mut current: Option<SegmentPartial> = None;
        self.reset_row(scratch);
        for bytes in records {
            let rec = reed_io::bcf::record_ref(bytes).map_err(EngineError::from)?;
            self.push_bcf_record(&rec, bcf, &mut current, &mut partial, scratch)?;
        }
        if let Some(segment) = current.take() {
            partial.segments.push(segment);
        }
        scratch.ld.flush(&mut partial);
        Ok(partial)
    }

    fn push_bcf_record(
        &self,
        rec: &reed_io::bcf::BcfRecordRef<'_>,
        bcf: &BcfNative,
        current: &mut Option<SegmentPartial>,
        partial: &mut StatsPartial,
        scratch: &mut StatsScratch,
    ) -> Result<(), EngineError> {
        let header = reed_io::bcf::site_header(rec).map_err(EngineError::from)?;
        let contig = bcf.dict.contig(header.chrom).ok_or_else(|| {
            EngineError::job(format!(
                "BCF contig index {} is not in the header dictionary",
                header.chrom
            ))
        })?;
        let Some(pos1) = u64::try_from(header.pos0)
            .ok()
            .and_then(|p| p.checked_add(1))
        else {
            return Err(EngineError::job("BCF POS is negative"));
        };
        let n_alt = u32::try_from(header.n_allele.saturating_sub(1)).unwrap_or(u32::MAX);
        let Some(pos0) = self.site_prefix(contig, pos1, n_alt, n_alt == 0, current, partial)?
        else {
            return Ok(());
        };
        scratch.counters = PolicyCounters::default();
        let query = reed_io::bcf::GtQuery {
            gt_key: bcf.gt_key,
            n_samples: self.n_samples,
            contig_name: contig,
            pos1,
        };
        reed_io::bcf::decode_record_gt(
            rec,
            &header,
            &query,
            &bcf.table,
            &mut scratch.codes,
            &mut scratch.counters,
        )
        .map_err(EngineError::from)?;
        self.push_row(
            pos0,
            &scratch.codes,
            &mut scratch.row,
            &mut scratch.sfs_caches,
            current,
            partial,
            &mut scratch.ld,
        )
    }
}

/// The LD analysis results the sink accumulated across contigs: the cohort
/// curve, one curve per named group, and the cohort pair rows.
struct LdResult {
    cohort: LdDecayAcc,
    groups: Vec<LdDecayAcc>,
    pairs: Vec<(Vec<u8>, u64, u64, f64)>,
}

/// Everything the sink hands back through the shared slot at `finish`.
struct SinkOutput {
    rows: String,
    windows: Vec<WindowRow>,
    /// Genome-wide projected spectra.
    sfs: Option<SfsPartial>,
    /// Genome-wide f-statistics.
    fstats: Option<FStatsPartial>,
    /// LD results (the parallel matrix analysis).
    ld: Option<LdResult>,
    n_invariant_sites: u64,
    n_no_alt_records: u64,
    n_windows: u64,
    n_records: u64,
    n_multiallelic_skipped: u64,
}

/// The sink's per-contig accumulation state.
struct CurrentContig {
    name: Vec<u8>,
    bins: ContigBins,
    /// Contig-level totals (§8.7): every touched bin folded in.
    totals: Vec<StatAcc>,
    first_pos0: u64,
    last_pos0: u64,
    /// Grid origin: windows start at 0 (1-based position 1), the
    /// the reference-tool grid (04-algorithms.md §4).
    anchor: u64,
    /// Start of the next window to emit (the watermark frontier).
    next_start: u64,
    /// End of the last emitted window (0-based).
    span_end: u64,
}

/// Ordered consumer: merges same-contig segments, emits complete windows
/// through the watermark (§8.3), and closes contigs with a contig-level
/// row (§8.7). Revisited contigs are refused — unsorted input would
/// otherwise silently double-count windows.
struct StatsSink {
    size: u64,
    step: u64,
    width: u64,
    kinds: Arc<[StatKind]>,
    seed: Vec<StatAcc>,
    levels: Levels,
    skip_empty: bool,
    /// Window-layout mode (structured rows + invariant policy).
    layout: Option<WindowLayoutOptions>,
    /// Structured rows (window-layout mode only).
    windows: Vec<WindowRow>,
    /// Site-wide invariant sites seen.
    n_invariant_sites: u64,
    /// Records with ALT = `.` (the all-sites check).
    n_no_alt_records: u64,
    contig_len: HashMap<Vec<u8>, u64>,
    current: Option<CurrentContig>,
    closed: Vec<Vec<u8>>,
    genome_totals: Vec<StatAcc>,
    /// Merged genome-wide projected spectra.
    sfs: Option<SfsPartial>,
    /// Merged genome-wide f-statistics.
    fstats: Option<FStatsPartial>,
    /// The LD plane layout (group-major); present only with an LD request.
    ld_layout: Option<Arc<reed_stats::SampleLayout>>,
    /// The matrix being accumulated for the pending contig.
    ld_matrix: Option<reed_stats::ContigMatrix>,
    /// Per-contig decay jobs in contig order (joined at `finish`, so the
    /// analysis overlaps the decode instead of stalling the reducer).
    ld_jobs: Vec<std::thread::JoinHandle<LdResult>>,
    /// The LD plan (the analysis needs the window and the bin edges).
    ld: Option<LdPlan>,
    genome_span: u64,
    rows: String,
    n_windows: u64,
    n_records: u64,
    n_multiallelic_skipped: u64,
    out: Arc<std::sync::Mutex<Option<SinkOutput>>>,
}

impl OrderedSink<StatsPartial> for StatsSink {
    fn accept(&mut self, _id: ChunkId, partial: &mut StatsPartial) -> Result<(), EngineError> {
        self.n_records += partial.n_records;
        self.n_multiallelic_skipped += partial.n_multiallelic_skipped;
        self.n_invariant_sites += partial.n_invariant_sites;
        self.n_no_alt_records += partial.n_no_alt_records;
        if let Some(src) = &partial.sfs {
            match &mut self.sfs {
                Some(dst) => dst.merge(src),
                None => self.sfs = Some(src.clone()),
            }
        }
        if let Some(src) = &partial.fstats {
            match &mut self.fstats {
                Some(dst) => dst.merge(src),
                None => self.fstats = Some(src.clone()),
            }
        }
        if let Some(layout) = self.ld_layout.clone() {
            // The worker's runs append to the pending contig's matrix; a
            // contig change finalizes it (the analysis runs per contig, so
            // at most one contig's matrix is resident).
            if let Some(src) = partial.ld.take() {
                for run in src.runs {
                    let same = self
                        .ld_matrix
                        .as_ref()
                        .is_some_and(|m| m.contig() == run.contig.as_slice());
                    if !same {
                        self.finalize_ld();
                        self.ld_matrix = Some(reed_stats::ContigMatrix::new(
                            run.contig,
                            Arc::clone(&layout),
                        ));
                    }
                    if let Some(matrix) = self.ld_matrix.as_mut() {
                        matrix.push_run(run.flat, run.pos0);
                    }
                }
            }
        }
        for mut segment in std::mem::take(&mut partial.segments) {
            let mut current = match self.current.take() {
                Some(mut cur) if cur.name == segment.contig => {
                    // Cross-chunk seam (04-algorithms.md §6): a later chunk
                    // may not rewind into positions this sink already
                    // emitted windows for.
                    if segment.first_pos0 < cur.last_pos0 {
                        return Err(EngineError::from(
                            reed_core::CoreError::PositionNotAscending {
                                prev_pos0: cur.last_pos0,
                                pos0: segment.first_pos0,
                            },
                        ));
                    }
                    segment.bins.fold(&mut cur.totals);
                    cur.bins.merge_bins(&mut segment.bins);
                    cur.first_pos0 = cur.first_pos0.min(segment.first_pos0);
                    cur.last_pos0 = cur.last_pos0.max(segment.last_pos0);
                    cur
                }
                Some(cur) => {
                    self.close_contig(cur)?;
                    Self::open_and_absorb(&mut segment, self.step, self.width, &self.seed)
                }
                None => Self::open_and_absorb(&mut segment, self.step, self.width, &self.seed),
            };
            self.emit_complete_windows(&mut current)?;
            self.current = Some(current);
        }
        Ok(())
    }

    fn finish(mut self) -> Result<u64, EngineError> {
        self.finalize_ld();
        let mut ld = None;
        self.join_ld(&mut ld)?;
        if let Some(current) = self.current.take() {
            self.close_contig(current)?;
        }
        if self.levels.genome {
            let bp = self.genome_span;
            let totals = self.genome_totals.clone();
            self.write_row("genome", 0, bp, bp, &totals);
        }
        // The all-sites check (04-algorithms.md §3): a file with no
        // invariant sites must be declared with --variant-only or
        // --bypass-invariant-check.
        if matches!(
            self.layout,
            Some(WindowLayoutOptions {
                invariant: InvariantMode::Require
            })
        ) && self.n_no_alt_records == 0
        {
            return Err(EngineError::job(format!(
                "input has no invariant sites ({} records, none with ALT = \".\"): declare \
                 --variant-only for an all-variant file, or --bypass-invariant-check to compute \
                 over present sites only (both overestimate π/dxy denominators)",
                self.n_records
            )));
        }
        *self.out.lock().expect("invariant: stats slot not poisoned") = Some(SinkOutput {
            ld,
            rows: self.rows,
            windows: std::mem::take(&mut self.windows),
            sfs: self.sfs.take(),
            fstats: self.fstats.take(),
            n_invariant_sites: self.n_invariant_sites,
            n_no_alt_records: self.n_no_alt_records,
            n_windows: self.n_windows,
            n_records: self.n_records,
            n_multiallelic_skipped: self.n_multiallelic_skipped,
        });
        Ok(self.n_windows)
    }
}

impl StatsSink {
    /// Hands the pending contig's matrix to a decay worker: the analysis
    /// runs off the reducer thread (it overlaps the decode), and `finish`
    /// joins the jobs in contig order — the merge order the deterministic
    /// output requires (docs/rounds/v6/02-models.md §2.2).
    fn finalize_ld(&mut self) {
        let Some(matrix) = self.ld_matrix.take() else {
            return;
        };
        let Some(plan) = &self.ld else { return };
        let edges = plan.decay_edges();
        let max_bp = plan.decay_max_bp.unwrap_or(100_000);
        let collect = plan.pairs;
        self.ld_jobs.push(std::thread::spawn(move || {
            let set = decay_contig(&matrix, max_bp, &edges, collect);
            LdResult {
                cohort: set.cohort,
                groups: set.groups,
                pairs: set
                    .pairs
                    .into_iter()
                    .map(|(a, b, r2)| (matrix.contig().to_vec(), a, b, r2))
                    .collect(),
            }
        }));
    }

    /// Joins the per-contig decay jobs in contig order.
    ///
    /// # Errors
    ///
    /// [`EngineError::Job`] when a decay worker panicked.
    fn join_ld(&mut self, out: &mut Option<LdResult>) -> Result<(), EngineError> {
        for handle in self.ld_jobs.drain(..) {
            let set = handle
                .join()
                .map_err(|_| EngineError::job("the LD decay worker panicked"))?;
            match out {
                Some(dst) => {
                    dst.cohort.merge(&set.cohort);
                    for (dst, src) in dst.groups.iter_mut().zip(&set.groups) {
                        dst.merge(src);
                    }
                    dst.pairs.extend(set.pairs);
                }
                None => *out = Some(set),
            }
        }
        Ok(())
    }

    /// Opens a contig's state and absorbs its first segment (folding the
    /// segment's bins into the contig totals and the bin map).
    fn open_and_absorb(
        segment: &mut SegmentPartial,
        step: u64,
        width: u64,
        seed: &[StatAcc],
    ) -> CurrentContig {
        // Grid from position 1: `[1, w]`, `[1+s, w+s]`, … (04 §4).
        let anchor = 0;
        let _ = step;
        let mut cur = CurrentContig {
            name: segment.contig.clone(),
            bins: ContigBins::new(BpBinning::new(width), seed),
            totals: seed.to_vec(),
            first_pos0: segment.first_pos0,
            last_pos0: segment.last_pos0,
            anchor,
            next_start: anchor,
            span_end: 0,
        };
        segment.bins.fold(&mut cur.totals);
        cur.bins.merge_bins(&mut segment.bins);
        cur
    }

    /// Emits every window that the watermark proves complete: with sorted
    /// input, no future site can land in `[s, s+size)` once the last seen
    /// position reaches `s+size−1` (§8.3). Consumed bins are freed. The
    /// contig length does not gate this loop — only the trailing windows
    /// in `close_contig` need it.
    #[allow(clippy::unnecessary_wraps)] // callers share the Result-returning shape of the window writers
    fn emit_complete_windows(&mut self, cur: &mut CurrentContig) -> Result<(), EngineError> {
        let watermark = cur.last_pos0.saturating_add(1);
        let mut s = cur.next_start;
        while s.saturating_add(self.size) <= watermark {
            let end = s + self.size;
            let bins = (s / self.width)..end.div_ceil(self.width);
            let skip = self.skip_empty && !cur.bins.any_bin_in(bins.clone());
            if !skip {
                let accs = cur.bins.finish_bins(bins);
                self.write_row(
                    &String::from_utf8_lossy(&cur.name),
                    s + 1,
                    end,
                    end - s,
                    &accs,
                );
                cur.span_end = end;
                self.n_windows += 1;
            }
            s = s.saturating_add(self.step);
        }
        cur.next_start = s;
        cur.bins.drop_bins_below(s / self.width);
        Ok(())
    }

    /// Closes a contig: emits the trailing windows (start at or before the
    /// last site, clipped to the contig length — the contig length is
    /// required here and is never fabricated, 01-review D7), then the
    /// contig-level row (§8.7), then folds into the genome total.
    fn close_contig(&mut self, mut cur: CurrentContig) -> Result<(), EngineError> {
        if self.closed.iter().any(|closed| closed == &cur.name) {
            return Err(EngineError::job(format!(
                "contig {} revisited after a later contig: input is not contig-sorted",
                String::from_utf8_lossy(&cur.name)
            )));
        }
        self.closed.push(cur.name.clone());
        let len = self.contig_len.get(&cur.name).copied().ok_or_else(|| {
            EngineError::job(format!(
                "contig {} has no ##contig length; window clipping needs it (pass --contig-lengths)",
                String::from_utf8_lossy(&cur.name)
            ))
        })?;
        let contig = String::from_utf8_lossy(&cur.name).into_owned();
        let mut s = cur.next_start;
        while s < len && s <= cur.last_pos0 {
            let end = s.saturating_add(self.size).min(len);
            let bins = (s / self.width)..end.div_ceil(self.width);
            let skip = self.skip_empty && !cur.bins.any_bin_in(bins.clone());
            if !skip {
                let accs = cur.bins.finish_bins(bins);
                self.write_row(&contig, s + 1, end, end - s, &accs);
                cur.span_end = end;
                self.n_windows += 1;
            }
            s = s.saturating_add(self.step);
        }
        if self.levels.contig {
            let bp = cur.span_end.saturating_sub(cur.anchor);
            self.write_row(&contig, cur.anchor + 1, cur.span_end, bp, &cur.totals);
        }
        for (dst, src) in self.genome_totals.iter_mut().zip(&cur.totals) {
            dst.merge(src);
        }
        self.genome_span += cur.span_end.saturating_sub(cur.anchor);
        Ok(())
    }

    /// Writes one TSV row at any level; `bp` is the span's base-pair count
    /// (the AccessibleBp/per-bp denominator at that level).
    fn write_row(&mut self, contig: &str, start1: u64, end0: u64, bp: u64, accs: &[StatAcc]) {
        let variant_only = matches!(
            self.layout,
            Some(WindowLayoutOptions {
                invariant: InvariantMode::VariantOnly
            })
        );
        if self.layout.is_some() {
            self.windows.push(WindowRow {
                contig: contig.to_string(),
                start1,
                end0,
                bp,
                accs: accs.to_vec(),
            });
        }
        self.rows.push_str(contig);
        self.rows.push('\t');
        self.rows.push_str(&start1.to_string());
        self.rows.push('\t');
        self.rows.push_str(&end0.to_string());
        for (kind, acc) in self.kinds.iter().zip(accs) {
            self.rows.push('\t');
            match acc.sfs_bins() {
                Some(bins) => {
                    for (i, bin) in bins.iter().enumerate() {
                        if i > 0 {
                            self.rows.push(';');
                        }
                        self.rows.push_str(&bin.to_string());
                    }
                }
                None => {
                    let value = if variant_only {
                        acc.value_variant_only(kind, bp)
                    } else {
                        acc.value(kind, bp)
                    };
                    match value {
                        Some(value) => {
                            write!(self.rows, "{value}").expect("invariant: writing to String");
                        }
                        None => self.rows.push('.'),
                    }
                }
            }
        }
        self.rows.push('\n');
    }
}

/// Runs windowed statistics over a BGZF VCF (`docs/spec/architecture.md`
/// milestone M4: the full 3-population statistic set).
///
/// # Errors
///
/// Transport, job, and request failures — including a sample-count mismatch
/// against the header, an invalid grid, or a group index without samples.
pub fn run_stats(
    path: Option<&Path>,
    request: StatsRequest,
    budget: &ThreadBudget,
) -> Result<StatsReport, EngineError> {
    let opened = crate::input::open(path).map_err(EngineError::from)?;
    run_stats_opened(opened, request, budget)
}

/// Runs windowed statistics over an already-opened input (the CLI opens
/// once so stdin streams; doc_v3 03-performance.md O5).
///
/// # Errors
///
/// As [`run_stats`].
pub fn run_stats_opened(
    opened: crate::input::OpenedInput,
    request: StatsRequest,
    budget: &ThreadBudget,
) -> Result<StatsReport, EngineError> {
    let grid = WindowGrid::Bp {
        size: request.window_size,
        step: request.step,
    };
    grid.validate().map_err(EngineError::from)?;
    let binning = BpBinning::from_grid(&grid)
        .ok_or_else(|| EngineError::job("stats driver requires a Bp grid"))?;

    let n_groups = request
        .gt_groups
        .iter()
        .flatten()
        .copied()
        .max()
        .map_or(0, |g| g as usize + 1);
    // The SFS is genome-wide and per group: validate the projections once
    // (docs/rounds/v5/05-algorithms.md §3.1).
    let fstats_plan = match &request.fstats {
        None => None,
        Some(fstats) => {
            if fstats.tests.is_empty() {
                return Err(EngineError::job("the f-statistics request lists no tests"));
            }
            for test in &fstats.tests {
                for g in test.pops() {
                    if g >= n_groups {
                        return Err(EngineError::job(format!(
                            "f-statistic {} names group {g} but the sample assignment declares                              {n_groups} group(s)",
                            test.id()
                        )));
                    }
                }
            }
            Some(FStatsPlan {
                tests: fstats.tests.clone().into(),
                block_bp: fstats.block_bp.max(1),
            })
        }
    };
    let ld_plan = request.ld.as_ref().map(|ld| LdPlan {
        decay_max_bp: ld.decay_max_bp,
        decay_bin_bp: ld.decay_bin_bp.max(1),
        pairs: ld.pairs,
    });
    let sfs_plan = match &request.sfs {
        None => None,
        Some(sfs) => {
            if sfs.projections.len() != n_groups {
                return Err(EngineError::job(format!(
                    "the SFS requests {} group projection(s) but the sample assignment declares                      {n_groups} group(s)",
                    sfs.projections.len()
                )));
            }
            for &n in &sfs.projections {
                if n < 2 {
                    return Err(EngineError::job(format!(
                        "projection size must be at least 2 chromosomes, got {n}"
                    )));
                }
            }
            Some(SfsPlan {
                projections: sfs.projections.clone().into(),
                folded: sfs.folded,
            })
        }
    };
    // `-s/-S/--keep/--remove`: excluded samples lose their group before the
    // group sizes are derived, so the denominators shrink with them.
    let mut gt_groups = request.gt_groups.clone();
    if let Some(keep) = &request.keep_samples {
        if keep.len() != gt_groups.len() {
            return Err(EngineError::job(format!(
                "sample keep mask has {} entries but the header declares {}",
                keep.len(),
                gt_groups.len()
            )));
        }
        for (slot, kept) in gt_groups.iter_mut().zip(keep) {
            if !kept {
                *slot = None;
            }
        }
    }
    let mut group_sizes = vec![0_u32; n_groups];
    for group in gt_groups.iter().flatten() {
        group_sizes[*group as usize] += 1;
    }
    // The LD plane layout: samples packed group-major (word-aligned), so a
    // group's pair sums read one contiguous slice and the cohort sums are
    // the exact sum of the group sums (docs/rounds/v6/02-models.md §2.2).
    let ld_layout = ld_plan.is_some().then(|| {
        Arc::new(reed_stats::SampleLayout::from_assignments(
            &gt_groups, n_groups,
        ))
    });
    let accs = build_accs(&request.kinds, &group_sizes)
        .map_err(|err| EngineError::job(format!("invalid stats request: {err}")))?;
    let keys: Vec<StatKey> = request.kinds.iter().map(StatKind::key).collect();

    let abort = reed_engine::abort::abort_for(budget);
    let header = opened.facts.clone();
    if header.n_samples() != request.gt_groups.len() {
        return Err(EngineError::job(format!(
            "request assigns {} samples but the header declares {}",
            request.gt_groups.len(),
            header.n_samples()
        )));
    }
    let contig_len: Arc<HashMap<Vec<u8>, u64>> = Arc::new({
        let mut map: HashMap<Vec<u8>, u64> = request.contig_lengths;
        for (name, len) in header
            .contig_names
            .iter()
            .zip(header.contig_lengths.iter().copied())
        {
            if let Some(len) = len {
                map.insert(name.as_bytes().to_vec(), len);
            }
        }
        map
    });

    // The wrapper's mask decision follows the same static subset as the
    // group assignment.
    let keep_for_filter: Arc<[bool]> = Arc::from(
        request
            .keep_samples
            .clone()
            .unwrap_or_else(|| vec![true; request.gt_groups.len()]),
    );
    let job = StatsJob {
        width: binning.width(),
        n_groups,
        gt_groups: gt_groups.into(),
        n_samples: request.gt_groups.len(),
        seed: accs.clone().into(),
        contig_len: Arc::clone(&contig_len),
        targets: request.targets.clone().into(),
        targets_exclude: request.targets_exclude.clone().into(),
        // The native BCF path exists only for BCF input (02-bcf.md §2.4);
        // the table is built once per run.
        bcf: opened.bcf.as_ref().map(|dict| BcfNative {
            dict: Arc::clone(dict),
            gt_key: reed_io::bcf::format_key(dict, b"GT"),
            table: reed_io::bcf::GtTable::new(stats_policies()),
        }),
        sfs: sfs_plan,
        fstats: fstats_plan,
        ld_layout: ld_layout.clone(),
    };
    let out = Arc::new(std::sync::Mutex::new(None));
    let sink = StatsSink {
        sfs: None,
        fstats: None,
        ld_layout,
        ld_matrix: None,
        ld_jobs: Vec::new(),
        ld: ld_plan,
        size: request.window_size,
        step: request.step,
        width: binning.width(),
        kinds: Arc::from(request.kinds.into_boxed_slice()),
        seed: accs.clone(),
        levels: request.levels,
        skip_empty: request.skip_empty,
        layout: request.layout,
        windows: Vec::new(),
        n_invariant_sites: 0,
        n_no_alt_records: 0,
        contig_len: (*contig_len).clone(),
        current: None,
        closed: Vec::new(),
        genome_totals: accs,
        genome_span: 0,
        rows: String::new(),
        n_windows: 0,
        n_records: 0,
        n_multiallelic_skipped: 0,
        out: Arc::clone(&out),
    };
    // The site/threshold flags run in the same pass through the filter's own
    // record pipeline; the inner (statistics) job then sees exactly what a
    // physical `reed filter` + `reed stats` would have produced. Without a
    // plan, BCF input takes the native record path (02-bcf.md §2.3–§2.5):
    // the job decodes GT bytes instead of rendering every column.
    let report = match &request.filter {
        Some(plan) => {
            let job = crate::bcf_input::InputJob::new(
                crate::filtered::FilteredJob::with_records(
                    Some(Arc::clone(plan)),
                    Arc::clone(&keep_for_filter),
                    StatsJob::run_records,
                    job,
                ),
                &opened,
                &abort,
            );
            let mut source = opened.source;
            let report = run(source.as_mut(), &job, sink, budget, &abort)?;
            job.finish(report.n_chunks, true)?;
            report
        }
        None => {
            let job = crate::bcf_input::InputJob::with_bcf_records(
                job,
                &opened,
                &abort,
                StatsJob::run_bcf_records,
            );
            let mut source = opened.source;
            let report = run(source.as_mut(), &job, sink, budget, &abort)?;
            job.finish(report.n_chunks, true)?;
            report
        }
    };
    let output = out
        .lock()
        .expect("invariant: stats slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the output");

    let mut tsv = String::from("contig\tstart\tend");
    for key in &keys {
        tsv.push('\t');
        tsv.push_str(key.as_str());
    }
    tsv.push('\n');
    tsv.push_str(&output.rows);
    Ok(StatsReport {
        tsv,
        windows: output.windows,
        ld: output.ld.map(|ld| {
            let bin = request.ld.as_ref().map_or(1000, |l| l.decay_bin_bp.max(1));
            let bins = |acc: &reed_stats::LdDecayAcc| {
                let means = acc.mean_r2();
                let counts = acc.counts().to_vec();
                counts
                    .iter()
                    .zip(&means)
                    .enumerate()
                    .map(|(i, (&n, m))| {
                        let start = (i as u64) * bin;
                        (start, start + bin, n, *m)
                    })
                    .collect::<Vec<_>>()
            };
            LdReport {
                decay: bins(&ld.cohort),
                group_decay: ld.groups.iter().map(bins).collect(),
                pairs: ld
                    .pairs
                    .iter()
                    .map(|(c, a, b, r2)| (String::from_utf8_lossy(c).into_owned(), *a, *b, *r2))
                    .collect(),
            }
        }),
        sfs: output.sfs.as_ref().map(|partial| SfsReport {
            spectra: partial
                .per_group
                .iter()
                .map(ProjectedSfsAcc::spectrum)
                .collect(),
            projections: partial.per_group.iter().map(ProjectedSfsAcc::n).collect(),
            folded: partial
                .per_group
                .first()
                .is_some_and(ProjectedSfsAcc::folded),
            n_sites: partial.per_group.iter().map(|a| a.site_tally().0).collect(),
            n_skipped_m_lt_n: partial.per_group.iter().map(|a| a.site_tally().1).collect(),
        }),
        fstats: output.fstats.as_ref().map(|partial| FStatsReport {
            rows: partial.per_test.iter().filter_map(FStatAcc::row).collect(),
        }),
        n_invariant_sites: output.n_invariant_sites,
        n_no_alt_records: output.n_no_alt_records,
        n_windows: output.n_windows,
        n_records: output.n_records,
        n_multiallelic_skipped: output.n_multiallelic_skipped,
        n_chunks: report.n_chunks,
        metrics: report.metrics,
        elapsed: report.elapsed,
    })
}
