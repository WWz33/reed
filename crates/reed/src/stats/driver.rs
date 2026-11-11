//! The statistics driver entry points: open-or-take the input, compile the
//! request into a plan, run the pipeline, and assemble the report.

use crate::stats::request::{
    FStatsPartial, FStatsPlan, FStatsReport, JointSfsReport, LdPlan, LdReport, Levels,
    PopgenReport, PopgenRequest, Populations, SfsPartial, SfsPlan, SfsReport, assign_groups,
    assign_groups_in_memory,
};
use crate::stats::run::{BcfNative, GridMode, StatsJob, stats_policies};
use crate::stats::sink::StatsSink;
use reed_core::grid::WindowGrid;
use reed_engine::budget::ThreadBudget;
use reed_engine::error::EngineError;
use reed_engine::pipeline::run;
use reed_stats::request::build_accs;
use reed_stats::sfs::{JointProjectedAcc, ProjectedSfsAcc};
use reed_stats::window::BpBinning;
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::sync::Arc;

/// Resolves the population request against the header's samples and applies
/// the static keep mask (stats/convert share this; doc_v13 D13-6). Returns
/// the per-sample group slots and the group names (first-appearance order).
///
/// # Errors
///
/// Group-file failures and a keep mask whose length differs from the header.
pub(crate) fn resolve_groups(
    populations: &Populations,
    force_samples: bool,
    sample_names: &[String],
    keep: Option<&[bool]>,
) -> Result<(Vec<Option<u32>>, Vec<String>), crate::error::Error> {
    let (gt_groups, group_names) = match populations {
        Populations::File(path) => assign_groups(path, sample_names, force_samples)?,
        Populations::Assignment(pairs) => {
            assign_groups_in_memory(pairs, sample_names, force_samples)?
        }
        Populations::All => (vec![Some(0); sample_names.len()], vec!["all".to_string()]),
    };
    let mut gt_groups = gt_groups;
    if let Some(keep) = keep {
        // One construction check for every driver (doc_v25 V25-07b).
        let keep = reed_core::sample::KeepMask::try_new(gt_groups.len(), keep.to_vec())
            .map_err(|e| crate::error::Error::usage(e.to_string()))?;
        for (slot, &kept) in gt_groups.iter_mut().zip(keep.iter()) {
            if !kept {
                *slot = None;
            }
        }
    }
    Ok((gt_groups, group_names))
}

/// Validates an SFS request against the group count and builds the plan
/// (the single implementation; stats and convert both call it, D13-6).
///
/// # Errors
///
/// A projection count that contradicts the groups, a projection below two
/// chromosomes, or one above a group's chromosome count (2 × samples).
pub(crate) fn build_sfs_plan(
    sfs: &crate::stats::request::SfsRequest,
    group_names: &[String],
    group_sizes: &[u32],
) -> Result<SfsPlan, crate::error::Error> {
    let n_groups = group_sizes.len();
    if sfs.projections.len() != n_groups {
        return Err(crate::error::Error::usage(format!(
            "the SFS requests {} group projection(s) but the sample assignment declares {n_groups} group(s)",
            sfs.projections.len()
        )));
    }
    for (group, &n) in sfs.projections.iter().enumerate() {
        if n < 2 {
            return Err(crate::error::Error::usage(format!(
                "projection size must be at least 2 chromosomes, got {n}"
            )));
        }
        // A projection above the group's chromosome count (2 x samples) can
        // never absorb a site: the run would quietly emit an all-zero
        // spectrum (doc_v15 D15-9; shared by popg and convert, D13-6).
        let max_chrom = 2 * u64::from(group_sizes[group]);
        if u64::from(n) > max_chrom {
            return Err(crate::error::Error::usage(format!(
                "--proj {n} exceeds population {}'s {max_chrom} chromosomes (2 x {} samples)",
                group_names[group], group_sizes[group]
            )));
        }
    }
    let joint_pairs: Vec<(usize, usize)> = (0..n_groups)
        .flat_map(|a| ((a + 1)..n_groups).map(move |b| (a, b)))
        .collect();
    Ok(SfsPlan {
        projections: sfs.projections.clone().into(),
        folded: sfs.folded,
        joint_pairs: joint_pairs.into(),
    })
}

/// Assembles an [`SfsReport`] from a merged SFS partial (stats and convert
/// share it, D13-6).
pub(crate) fn sfs_report(partial: &crate::stats::request::SfsPartial, plan: &SfsPlan) -> SfsReport {
    SfsReport {
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
        joint: partial.joint.as_ref().map(|accs| JointSfsReport {
            pairs: plan.joint_pairs.to_vec(),
            spectra: accs.iter().map(JointProjectedAcc::values).collect(),
            projections: accs.iter().map(|a| a.shape()).collect(),
            folded: accs.first().is_some_and(JointProjectedAcc::folded),
            n_sites: accs.iter().map(|a| a.site_tally().0).collect(),
            n_skipped_m_lt_n: accs.iter().map(|a| a.site_tally().1).collect(),
        }),
    }
}

/// Runs windowed statistics over a BGZF VCF (`docs/spec/architecture.md`
/// milestone M4: the full 3-population statistic set).
///
/// # Errors
///
/// Transport, job, and request failures — including a sample-count mismatch
/// against the header, an invalid grid, or a group index without samples.
pub fn run_popgen(
    path: Option<&Path>,
    request: &PopgenRequest,
    budget: &ThreadBudget,
) -> Result<PopgenReport, crate::error::Error> {
    let opened = crate::input::open(path).map_err(EngineError::from)?;
    run_popgen_opened(opened, request, budget)
}

/// Allocates the worker half of the budget to LD, leaving the other half
/// for input decoding. Reader + reducer account for the two non-worker
/// slots, so the whole run stays within `budget.threads` (V27-16).
fn ld_worker_share(budget: &ThreadBudget, requested: bool) -> usize {
    if requested && !budget.is_sequential() {
        budget.workers() / 2
    } else {
        0
    }
}

/// Runs windowed statistics over an already-opened input (the CLI opens
/// once so stdin streams; doc_v3 03-performance.md O5).
///
/// # Errors
///
/// As [`run_popgen`].
pub fn run_popgen_opened(
    opened: crate::input::OpenedInput,
    request: &PopgenRequest,
    budget: &ThreadBudget,
) -> Result<PopgenReport, crate::error::Error> {
    // ---- Request → executable plan (08 §6.2). Every failure up to the end
    // of the validation section below is request construction (Usage). ----
    let sample_names = opened.facts.sample_names.clone();
    let (gt_groups, group_names) = resolve_groups(
        &request.populations,
        request.force_samples,
        &sample_names,
        request.keep_samples.as_deref(),
    )?;
    // The biallelic gate is applied once in `StatsJob::site_prefix`, where
    // its filtered-record count reaches `PopgenReport` and the window log.
    // Do not also put it in this wrapper: a record dropped there would be
    // invisible to the statistics counters (doc_v27 V27-39). `-M 2` is the
    // same default gate; `-M > 2` is refused by the CLI.
    let filter: Option<Arc<reed_filter::plan::FilterPlan>> = match &request.site {
        Some(spec) => {
            let mut spec = spec.clone();
            spec.site.max_alt_alleles = None;
            Some(Arc::new(reed_filter::plan::compile(
                &spec,
                &reed_io::header::tag_table(&opened.facts.header),
                "",
                false,
            )?))
        }
        None => None,
    };
    let n_groups = group_names.len();

    request.grid.validate()?;
    // With no window-level statistic there is no window row at any level:
    // SFS, f-statistics and LD are genome-wide computations carried by their
    // own partials. Building a Bp grid here used to force a contig length at
    // `close_contig` even when every level was disabled (doc_v27 V27-07).
    // `GridMode::None` avoids both that false dependency and the unused bin
    // accumulators. A Regions grid still requires the window level, as
    // before; `--levels contig,genome` keeps the Bp span bookkeeping.
    let grid_mode = if request.kinds.is_empty() {
        if matches!(request.grid, WindowGrid::Regions(_)) {
            return Err(crate::error::Error::usage(
                "--bed-file needs a window-level statistic (for example --pi)",
            ));
        }
        GridMode::None
    } else {
        match &request.grid {
            WindowGrid::Bp { size, step } => {
                let width = BpBinning::from_grid(&request.grid)
                    .ok_or_else(|| crate::error::Error::usage("stats driver requires a Bp grid"))?
                    .width();
                GridMode::Bp {
                    size: *size,
                    step: *step,
                    width,
                }
            }
            WindowGrid::Regions(regions) => {
                if !matches!(request.levels, Levels { window: true, .. }) {
                    return Err(crate::error::Error::usage(
                        "--bed-file windows require the window level (--levels window)",
                    ));
                }
                // Region ends carry the clipping the Bp path takes from the
                // contig length; nothing else needs it.
                GridMode::Regions(Arc::clone(regions))
            }
            WindowGrid::Sites { .. } => {
                return Err(crate::error::Error::usage(
                    "site-count windows are not implemented (use --window/--step or --bed-file)",
                ));
            }
        }
    };

    // The SFS is genome-wide and per group: validate the projections once
    // (docs/rounds/v5/05-algorithms.md §3.1).
    let fstats_plan = match &request.fstats {
        None => None,
        Some(fstats) => {
            if fstats.tests.is_empty() {
                return Err(crate::error::Error::usage(
                    "the f-statistics request lists no tests",
                ));
            }
            for test in &fstats.tests {
                for g in test.pops() {
                    if g >= n_groups {
                        return Err(crate::error::Error::usage(format!(
                            "f-statistic {} names group {g} but the sample assignment declares {n_groups} group(s)",
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
        pairs_path: ld.pairs_path.clone(),
    });
    // `-s/-S/--keep/--remove`: the keep mask already emptied excluded
    // samples' group slots inside `resolve_groups`.
    let mut group_sizes = vec![0_u32; n_groups];
    for group in gt_groups.iter().flatten() {
        group_sizes[*group as usize] += 1;
    }
    let sfs_plan = request
        .sfs
        .as_ref()
        .map(|sfs| build_sfs_plan(sfs, &group_names, &group_sizes))
        .transpose()?;
    // The LD plane layout: samples packed group-major (word-aligned), so a
    // group's pair sums read one contiguous slice and the cohort sums are
    // the exact sum of the group sums (docs/rounds/v6/02-models.md §2.2).
    let ld_layout = ld_plan.is_some().then(|| {
        Arc::new(reed_stats::SampleLayout::from_assignments(
            &gt_groups,
            n_groups,
            request.keep_samples.as_deref(),
        ))
    });
    let accs = build_accs(&request.kinds, &group_sizes)
        .map_err(|err| crate::error::Error::usage(format!("invalid stats request: {err}")))?;
    let abort = reed_engine::abort::abort_for(budget);
    let header = opened.facts.clone();
    let contig_len: Arc<HashMap<Vec<u8>, u64>> = Arc::new({
        let mut map: HashMap<Vec<u8>, u64> = HashMap::new();
        for (name, len) in header
            .contig_names
            .iter()
            .zip(header.contig_lengths.iter().copied())
        {
            if let Some(len) = len {
                map.insert(name.as_bytes().to_vec(), len);
            }
        }
        // `--contig-lengths` is the operator's explicit value and wins over
        // the header (the CLI prelude already applies the same precedence).
        map.extend(request.contig_lengths.clone());
        map
    });

    // The wrapper's mask decision follows the same static subset as the
    // group assignment.
    let keep_for_filter: Arc<[bool]> = Arc::from(
        request
            .keep_samples
            .clone()
            .unwrap_or_else(|| vec![true; sample_names.len()]),
    );
    let active_samples = reed_core::sample::KeepMask::try_new(
        sample_names.len(),
        gt_groups.iter().map(Option::is_some).collect(),
    )
    .expect("invariant: one group slot per sample");
    let job = StatsJob {
        grid: grid_mode.clone(),
        stream_windows: request.compact_window_accumulators,
        n_groups,
        gt_groups: gt_groups.into(),
        active_samples,
        n_samples: sample_names.len(),
        seed: accs.clone().into(),
        contig_len: Arc::clone(&contig_len),
        select: Arc::new(request.select.clone()),
        // The native BCF path exists only for BCF input (02-bcf.md §2.4);
        // the table is built once per run.
        region_accs: Arc::new(std::sync::Mutex::new(HashMap::new())),
        bcf: opened.bcf.as_ref().map(|dict| BcfNative {
            dict: Arc::clone(dict),
            gt_key: reed_io::bcf::format_key(dict, b"GT"),
            table: reed_io::bcf::GtTable::new(stats_policies()),
        }),
        sfs: sfs_plan.clone(),
        fstats: fstats_plan.clone(),
        ld_layout: ld_layout.clone(),
    };
    let out = Arc::new(std::sync::Mutex::new(None));
    // This is an output-session promise: `--ld-pairs` always publishes a
    // table, including its header-only zero-pair form. The temporary exists
    // now; the final path changes only after the input checks below pass.
    let pairs_file = match ld_plan.as_ref().and_then(|plan| plan.pairs_path.as_deref()) {
        Some(path) => crate::stats::sink::open_pairs(path)?,
        None => None,
    };
    // `--threads` is one budget, not one budget per stage. Two slots belong
    // to the reader/reducer; the remaining workers are shared by the input
    // pipeline and the LD pool. This removes both the old quadratic growth
    // and the linear `2 * (threads - 2)` oversubscription (V27-16).
    let ld_workers = ld_worker_share(budget, ld_plan.is_some());
    let pipeline_workers = budget.workers() - ld_workers;
    let pipeline_threads = if ld_workers == 0 {
        budget.threads
    } else {
        pipeline_workers + 2
    };
    let pipeline_budget = ThreadBudget {
        threads: pipeline_threads,
        mem: budget.mem,
    };
    let ld_executor = match ld_plan.as_ref() {
        Some(_) => Some(crate::stats::sink::LdExecutor::new(
            ld_workers,
            abort.handle.clone(),
            abort.signal.clone(),
        )?),
        None => None,
    };
    let output_session = request
        .window_output
        .as_ref()
        .filter(|_| request.compact_window_accumulators && !request.kinds.is_empty())
        .map(|config| {
            crate::stats::writers::session::WindowOutputSession::new(
                config,
                &request.kinds,
                &group_names,
                request.layout.as_ref(),
            )
        })
        .transpose()?;
    let sink = StatsSink {
        // Request-scoped accumulators exist even when no record reaches
        // them. An empty but successful run must still produce a valid zero
        // SFS/f-statistics report, not `None` and a stale previous output
        // (doc_v27 V27-05).
        sfs: sfs_plan.as_ref().and_then(SfsPartial::new),
        fstats: fstats_plan.as_ref().map(FStatsPartial::new),
        pairs_file,
        ld_layout,
        ld_matrix: None,
        ld_executor,
        ld_merged: None,
        max_ld_jobs: budget.auxiliary(),
        ld: ld_plan,
        grid: grid_mode,
        kinds: request.kinds.clone().into(),
        seed: accs.clone(),
        levels: request.levels,
        compact_window_accumulators: request.compact_window_accumulators,
        skip_empty: request.skip_empty,
        layout: request.layout,
        windows: Vec::new(),
        window_scratch: Vec::new(),
        output_session,
        n_invariant_sites: 0,
        n_no_alt_records: 0,
        contig_len: (*contig_len).clone(),
        current: None,
        closed: HashSet::new(),
        genome_totals: accs,
        genome_span: 0,
        genome_denom_bp: 0,
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
    let report = match &filter {
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
                &pipeline_budget,
            );
            let mut source = opened.source;
            let report = run(source.as_mut(), &job, sink, &pipeline_budget, &abort)?;
            job.finish(report.n_chunks, true)?;
            report
        }
        None => {
            let job = crate::bcf_input::InputJob::with_bcf_records(
                job,
                &opened,
                &abort,
                &pipeline_budget,
                StatsJob::run_bcf_records,
            );
            let mut source = opened.source;
            let report = run(source.as_mut(), &job, sink, &pipeline_budget, &abort)?;
            job.finish(report.n_chunks, true)?;
            report
        }
    };
    let output = out
        .lock()
        .expect("invariant: stats slot not poisoned")
        .take()
        .expect("invariant: sink finish stored the output");

    let sfs_report = output.sfs.as_ref().map(|partial| {
        sfs_report(
            partial,
            sfs_plan
                .as_ref()
                .expect("invariant: an SFS partial implies a plan"),
        )
    });
    let window_files = match output.window_session {
        Some(session) => {
            let config = request
                .window_output
                .as_ref()
                .expect("invariant: a session implies a config");
            session.finish(
                output.n_windows,
                output.n_records,
                output.n_invariant_sites,
                output.n_multiallelic_skipped,
                request
                    .layout
                    .as_ref()
                    .map_or(crate::stats::request::InvariantMode::Require, |layout| {
                        layout.invariant
                    }),
                &config.prefix,
                sfs_report.as_ref(),
            )?
        }
        None => Vec::new(),
    };

    // Publish the streaming pair table only now: every input-integrity check
    // (`job.finish` above detects a record truncated at the end of a BCF)
    // has passed. Committing inside `sink.finish` published a truncated
    // `.ld` over the previous run's good file and then failed the run
    // (doc_v27 V27-02). Dropping the `AtomicFile` on any early return
    // removes the temp and leaves the old `.ld` untouched.
    if let Some(atomic) = output.pending_commit {
        atomic.commit().map_err(EngineError::Io)?;
    }

    let ld_report = request.ld.as_ref().map(|ld_request| {
        let bin = ld_request.decay_bin_bp.max(1);
        let n_bins = ld_request.decay_max_bp.map_or(0, |max| {
            usize::try_from(max.div_ceil(bin)).unwrap_or(usize::MAX)
        });
        // `n_bins` is derived from the request, not from how many contigs
        // yielded pairs. A zero-pair run therefore has the same valid table
        // shape (0 / NA in every requested bin) as a nonempty one, and the
        // writer can replace the old decay file unconditionally (V27-05).
        let bins = |acc: Option<&reed_stats::LdDecayAcc>| {
            let means = acc.map(reed_stats::LdDecayAcc::mean_r2);
            let counts = acc.map(reed_stats::LdDecayAcc::counts);
            (0..n_bins)
                .map(|i| {
                    let start = (i as u64) * bin;
                    (
                        start,
                        start + bin,
                        counts.as_ref().and_then(|v| v.get(i)).copied().unwrap_or(0),
                        means.as_ref().and_then(|v| v.get(i)).copied().flatten(),
                    )
                })
                .collect::<Vec<_>>()
        };
        LdReport {
            decay: bins(output.ld.as_ref().map(|ld| &ld.cohort)),
            group_decay: (0..n_groups)
                .map(|g| bins(output.ld.as_ref().and_then(|ld| ld.groups.get(g))))
                .collect(),
            pairs: output
                .ld
                .as_ref()
                .map(|ld| {
                    ld.pairs
                        .iter()
                        .map(|(c, a, b, r2)| (String::from_utf8_lossy(c).into_owned(), *a, *b, *r2))
                        .collect()
                })
                .unwrap_or_default(),
        }
    });

    Ok(PopgenReport {
        windows: output.windows,
        window_files,
        ld: ld_report,
        sfs: sfs_report,
        pairs_path: output.pairs_path,
        fstats: output.fstats.as_ref().map(|partial| FStatsReport {
            rows: partial
                .per_test
                .iter()
                .map(|acc| acc.row(&partial.grid))
                .collect(),
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

#[cfg(test)]
mod thread_allocation_tests {
    use super::*;

    #[test]
    fn ld_and_pipeline_workers_share_one_budget() {
        for threads in 1..=64 {
            let budget = ThreadBudget {
                threads,
                mem: reed_engine::budget::MemBudget::from_bytes(u64::MAX),
            };
            let ld = ld_worker_share(&budget, true);
            let pipeline_workers = budget.workers() - ld;
            let pipeline_threads = if ld == 0 {
                budget.threads
            } else {
                pipeline_workers + 2
            };
            assert!(pipeline_threads + ld <= threads, "threads={threads}");
            if threads <= 2 {
                assert_eq!(ld, 0);
                assert_eq!(pipeline_threads, threads);
            } else {
                assert_eq!(ld, budget.workers() / 2, "threads={threads}");
            }
        }
    }
}
