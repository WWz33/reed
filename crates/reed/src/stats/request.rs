//! Request, plan and report types of the statistics driver: what a caller asks
//! for (`PopgenRequest` and its sub-requests) and what comes back (`PopgenReport`).
//! The run driver lives in `run.rs`; the ordered sink in `sink.rs`.

use reed_core::count::SiteGroupCounts;
use reed_core::grid::WindowGrid;
pub use reed_stats::fst::FstKind;
use reed_stats::jackknife::BpBlockJackknife;
use reed_stats::patterson::{PattersonKind, f2_site, f3_site, f4_site};
pub use reed_stats::pbs::PbsKind;
use reed_stats::request::StatAcc;
pub use reed_stats::request::StatKind;
use reed_stats::sfs::{JointProjectedAcc, ProjectedSfsAcc};
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Duration;

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
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Levels {
    /// One row per window (the default).
    pub window: bool,
    /// One row per contig, emitted when the contig closes.
    pub contig: bool,
    /// One `genome` row, emitted at the end.
    pub genome: bool,
}

impl Default for Levels {
    /// The default is the window level (doc_v11 N1): `Levels`'s "zero value"
    /// is not the all-`false` bool tuple — a request built from
    /// `PopgenRequest::default()` must still emit its window rows.
    fn default() -> Self {
        Self {
            window: true,
            contig: false,
            genome: false,
        }
    }
}

impl Levels {
    /// Parses a comma-separated `window,contig,genome` list.
    ///
    /// # Errors
    ///
    /// A usage error on an unknown level token.
    pub fn parse(text: &str) -> Result<Self, crate::error::Error> {
        let mut levels = Self {
            window: false,
            contig: false,
            genome: false,
        };
        for token in text.split(',').map(str::trim).filter(|t| !t.is_empty()) {
            match token {
                "window" => levels.window = true,
                "contig" => levels.contig = true,
                "genome" => levels.genome = true,
                other => {
                    return Err(crate::error::Error::usage(format!(
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
    /// 1-based inclusive end (equivalently, the 0-based exclusive end).
    pub end: u64,
    /// Window span in bp (the per-bp denominator).
    pub bp: u64,
    /// Accumulators in `kinds` order.
    pub accs: Vec<StatAcc>,
}

/// The population assignment for a stats run (`--pop` vs the bare
/// `--sfs`/`--ld` one-group synthesis). The assignment against the header's
/// sample names happens inside `run_popgen_opened` (08 §6.2: request
/// construction in the front end, executable plans in the driver).
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub enum Populations {
    /// One group holding every sample, named "all".
    #[default]
    All,
    /// A `sample<TAB>population` file (`--pop`).
    File(PathBuf),
    /// An in-memory `(sample, population)` list (doc_v11 N11): the same rows
    /// a populations file carries, for callers that already hold the
    /// mapping. Group ids follow first-appearance order, as in the file.
    Assignment(Vec<(String, String)>),
}

impl Populations {
    /// The group names in assignment (first-appearance) order, without
    /// opening the input: kinds/fstats/sfs building resolves names to
    /// indices before the run. The run re-reads the file in
    /// [`assign_groups`]; both readers share `reed-core`'s line rules, so a
    /// successful run sees the same names in the same order.
    ///
    /// # Errors
    ///
    /// An unreadable or malformed populations file (a `--pop` argument
    /// error, cli.md §9 → usage).
    pub fn group_names(&self) -> Result<Vec<String>, crate::error::Error> {
        match self {
            Self::All => Ok(vec!["all".to_string()]),
            Self::File(path) => reed_core::group::GroupTable::names_from_tsv(path)
                .map_err(crate::error::Error::from),
            Self::Assignment(pairs) => Ok(first_appearance_names(
                pairs.iter().map(|(_, group)| group.as_str()),
            )),
        }
    }
}

/// The distinct names in first-appearance order (the assignment order).
fn first_appearance_names<'a>(names: impl Iterator<Item = &'a str>) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for name in names {
        if !out.iter().any(|n| n == name) {
            out.push(name.to_string());
        }
    }
    out
}

/// One `run_popgen` request: statistics, group assignment, window grid.
#[derive(Clone, Debug, Default)]
pub struct PopgenRequest {
    /// Statistics to compute, in output-column order.
    pub kinds: Vec<StatKind>,
    /// Site rules in typed form; compiled against the header inside
    /// `run_popgen_opened` (08 §6.2: the front end builds typed requests, it
    /// does not compile plans). `None` = no site-level filtering.
    pub site: Option<reed_filter::spec::FilterSpec>,
    /// Population assignment; resolved against the header's sample names
    /// inside `run_popgen_opened`.
    pub populations: Populations,
    /// `--force-samples`: population-file rows naming samples absent from
    /// the header are skipped with a warning instead of failing.
    pub force_samples: bool,
    /// The window grid: fixed bp windows (`--window/--step`) or explicit
    /// per-contig regions (`--bed-file`).
    pub grid: WindowGrid,
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
pub(super) struct LdPlan {
    pub(super) decay_max_bp: Option<u64>,
    pub(super) decay_bin_bp: u64,
    pub(super) pairs: bool,
}

impl LdPlan {
    /// Decay bin edges: `[0, bin), [bin, 2·bin), …` up to `decay_max_bp`.
    pub(super) fn decay_edges(&self) -> Vec<u64> {
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
pub(super) struct LdRun {
    pub(super) contig: Vec<u8>,
    pub(super) flat: Vec<u64>,
    pub(super) pos0: Vec<u64>,
}

/// Per-chunk LD partial: the runs the worker filled. The pairs are computed
/// position-block parallel from the accumulated matrices after the stream
/// ends (docs/rounds/v6/02-models.md §2.2), so no head/tail carry exists.
#[derive(Clone, Debug, Default)]
pub(super) struct LdPartial {
    pub(super) runs: Vec<LdRun>,
}

/// A projected site-frequency-spectrum request (docs/rounds/v5/05-algorithms.md
/// §3.1): one genome-wide spectrum per population; with two or more groups
/// the pairwise joint (2D) spectra are written alongside (cli.md §7.2,
/// easySFS `dadi_twoD_sfs_combinations`).
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
    /// Per pair: the joint (2D) spectra in `(i, j)` row-major order, one per
    /// population pair in input order (`None` when not requested).
    pub joint: Option<JointSfsReport>,
}

/// The pairwise joint spectra an [`SfsRequest`] with `joint` produced.
#[derive(Clone, Debug)]
pub struct JointSfsReport {
    /// Per pair: the two group indices `(g1, g2)`, input order.
    pub pairs: Vec<(usize, usize)>,
    /// Per pair: the spectrum, row-major `i * (n2 + 1) + j`.
    pub spectra: Vec<Vec<f64>>,
    /// Per pair: the projection sizes `(n1, n2)`.
    pub projections: Vec<(u32, u32)>,
    /// Whether the spectra are folded.
    pub folded: bool,
    /// Per pair: sites absorbed by both groups.
    pub n_sites: Vec<u64>,
    /// Per pair: sites skipped for `m < n` in either group.
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
pub(super) struct FStatAcc {
    pub(super) test: FStatTest,
    /// Σnum, Σden (the jackknife's inputs).
    pub(super) num: reed_core::math::sum::NeumaierSum,
    pub(super) den: reed_core::math::sum::NeumaierSum,
    pub(super) n_sites: u64,
    pub(super) jk: BpBlockJackknife,
}

impl FStatAcc {
    pub(super) fn new(test: FStatTest, block_bp: u64) -> Self {
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
    pub(super) fn site(&self, row: &[SiteGroupCounts]) -> Option<(f64, f64)> {
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

    pub(super) fn push(&mut self, contig: &[u8], pos0: u64, row: &[SiteGroupCounts]) {
        let Some((num, den)) = self.site(row) else {
            return;
        };
        self.num.add(num);
        self.den.add(den);
        self.n_sites += 1;
        self.jk.push(contig, pos0, num, den);
    }

    pub(super) fn merge(&mut self, other: &Self) {
        self.num.merge(&other.num);
        self.den.merge(&other.den);
        self.n_sites += other.n_sites;
        self.jk.merge(&other.jk);
    }

    pub(super) fn row(&self) -> Option<FStatRow> {
        // Additive statistics (f2/f3/f4) report Σnum with a sum-scale
        // delete-one-block SE — admixtools reports f2 as a (weighted) sum
        // over sites, so estimate and se share the sum scale (doc_v7 D7-9).
        // The ratios (D, f4-ratio) report Σnum/Σden over the pair.
        let (estimate, se, n_blocks) = match self.test {
            FStatTest::Dstat { .. } | FStatTest::F4Ratio { .. } => {
                let den = self.den.value();
                if den == 0.0 {
                    return None;
                }
                let (_, se, nb) =
                    self.jk
                        .jackknife()
                        .unwrap_or((self.num.value() / den, f64::NAN, 0));
                (self.num.value() / den, se, nb)
            }
            _ => {
                let (_, se, nb) =
                    self.jk
                        .jackknife_sum()
                        .unwrap_or((self.num.value(), f64::NAN, 0));
                (self.num.value(), se, nb)
            }
        };
        let z = if se.is_finite() && se > 0.0 {
            estimate / se
        } else {
            f64::NAN
        };
        Some(FStatRow {
            test: self.test.clone(),
            estimate,
            se,
            z,
            n_blocks,
            n_sites: self.n_sites,
        })
    }
}

/// The validated f-statistics plan shared by the job and the sink.
#[derive(Clone, Debug)]
pub(super) struct FStatsPlan {
    pub(super) tests: Arc<[FStatTest]>,
    pub(super) block_bp: u64,
}

/// Per-chunk f-statistics partials.
#[derive(Clone, Debug)]
pub(super) struct FStatsPartial {
    pub(super) per_test: Vec<FStatAcc>,
}

impl FStatsPartial {
    pub(super) fn new(plan: &FStatsPlan) -> Self {
        Self {
            per_test: plan
                .tests
                .iter()
                .map(|test| FStatAcc::new(test.clone(), plan.block_bp))
                .collect(),
        }
    }

    pub(super) fn merge(&mut self, other: &Self) {
        for (dst, src) in self.per_test.iter_mut().zip(&other.per_test) {
            dst.merge(src);
        }
    }
}

/// The validated SFS plan shared by the job and the sink (crate-visible:
/// the convert driver rebuilds it through the same `build_sfs_plan`, D13-6).
#[derive(Clone, Debug)]
pub(crate) struct SfsPlan {
    pub(crate) projections: Arc<[u32]>,
    pub(crate) folded: bool,
    /// Group pairs for the joint spectra, input order (`combinations`).
    pub(crate) joint_pairs: Arc<[(usize, usize)]>,
}

/// Per-chunk SFS partials: one projected spectrum per group, plus one joint
/// spectrum per pair when requested.
#[derive(Clone, Debug, Default)]
pub(crate) struct SfsPartial {
    pub(super) per_group: Vec<ProjectedSfsAcc>,
    pub(super) joint: Option<Vec<JointProjectedAcc>>,
}

impl SfsPartial {
    pub(crate) fn new(plan: &SfsPlan) -> Option<Self> {
        let mut per_group = Vec::with_capacity(plan.projections.len());
        for &n in plan.projections.iter() {
            per_group.push(ProjectedSfsAcc::new(n, plan.folded)?);
        }
        let joint = if plan.joint_pairs.is_empty() {
            None
        } else {
            let mut accs = Vec::with_capacity(plan.joint_pairs.len());
            for &(g1, g2) in plan.joint_pairs.iter() {
                accs.push(JointProjectedAcc::new(
                    plan.projections[g1],
                    plan.projections[g2],
                    plan.folded,
                )?);
            }
            Some(accs)
        };
        Some(Self { per_group, joint })
    }

    pub(crate) fn merge(&mut self, other: &Self) {
        for (dst, src) in self.per_group.iter_mut().zip(&other.per_group) {
            dst.merge(src);
        }
        if let (Some(dst), Some(src)) = (&mut self.joint, &other.joint) {
            for (d, s) in dst.iter_mut().zip(src) {
                d.merge(s);
            }
        }
    }
}

/// The run's output and counters.
#[derive(Clone, Debug)]
pub struct PopgenReport {
    /// Structured window rows, one per selected level; both the
    /// window-layout files and the flat table render from these.
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
    /// The projected spectra when [`PopgenRequest::sfs`] was set.
    pub sfs: Option<SfsReport>,
    /// The Patterson f-statistics when [`PopgenRequest::fstats`] was set.
    pub fstats: Option<FStatsReport>,
    /// The LD tables when [`PopgenRequest::ld`] was set.
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
/// A usage error for an unreadable file, a malformed row, or an unknown
/// sample (without `force_samples`).
pub fn assign_groups(
    groups_path: &Path,
    sample_names: &[String],
    force_samples: bool,
) -> Result<(Vec<Option<u32>>, Vec<String>), crate::error::Error> {
    let samples = sample_table(sample_names)?;
    let table = reed_core::group::GroupTable::from_tsv(
        groups_path,
        &samples,
        reed_core::group::GroupFileOptions {
            ignore_missing_samples: force_samples,
        },
    )
    .map_err(crate::error::Error::from)?;
    group_masks(&table, sample_names.len())
}

/// The in-memory twin of [`assign_groups`] (doc_v11 N11): the same
/// first-appearance group order and the same unknown-sample rule, over
/// `(sample, population)` pairs the caller already holds.
///
/// # Errors
///
/// A usage error for an unknown sample (without `force_samples`).
pub fn assign_groups_in_memory(
    pairs: &[(String, String)],
    sample_names: &[String],
    force_samples: bool,
) -> Result<(Vec<Option<u32>>, Vec<String>), crate::error::Error> {
    let samples = sample_table(sample_names)?;
    let rows: Vec<(&str, &str)> = if force_samples {
        pairs
            .iter()
            .filter(|(sample, _)| samples.get(sample.trim()).is_some())
            .map(|(sample, group)| (sample.trim(), group.trim()))
            .collect()
    } else {
        pairs
            .iter()
            .map(|(sample, group)| (sample.trim(), group.trim()))
            .collect()
    };
    let table = reed_core::group::GroupTable::from_assignment(&samples, &rows)
        .map_err(crate::error::Error::from)?;
    group_masks(&table, sample_names.len())
}

fn sample_table(
    sample_names: &[String],
) -> Result<reed_core::sample::SampleTable, crate::error::Error> {
    reed_core::sample::SampleTable::from_names(sample_names.iter().map(String::as_str))
        .map_err(crate::error::Error::from)
}

/// Flattens a [`GroupTable`] into the per-sample group slots the job reads,
/// plus the group names in table order.
fn group_masks(
    table: &reed_core::group::GroupTable,
    n_samples: usize,
) -> Result<(Vec<Option<u32>>, Vec<String>), crate::error::Error> {
    let mut gt_groups: Vec<Option<u32>> = vec![None; n_samples];
    let mut names: Vec<String> = Vec::with_capacity(table.n_groups());
    for (group_idx, group) in table.iter().enumerate() {
        names.push(group.name.clone());
        let group_idx = u32::try_from(group_idx)
            .map_err(|_| crate::error::Error::usage("more than u32::MAX groups"))?;
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
