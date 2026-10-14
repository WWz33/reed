//! The statistics driver entry points: open-or-take the input, compile the
//! request into a plan, run the pipeline, and assemble the report.

use crate::stats::request::{
    FStatAcc, FStatsPlan, FStatsReport, JointSfsReport, LdPlan, LdReport, Levels, Populations,
    SfsPlan, SfsReport, StatsReport, StatsRequest, assign_groups, assign_groups_in_memory,
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
use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

/// Runs windowed statistics over a BGZF VCF (`docs/spec/architecture.md`
/// milestone M4: the full 3-population statistic set).
///
/// # Errors
///
/// Transport, job, and request failures — including a sample-count mismatch
/// against the header, an invalid grid, or a group index without samples.
pub fn run_stats(
    path: Option<&Path>,
    request: &StatsRequest,
    budget: &ThreadBudget,
) -> Result<StatsReport, crate::error::Error> {
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
    request: &StatsRequest,
    budget: &ThreadBudget,
) -> Result<StatsReport, crate::error::Error> {
    // ---- Request → executable plan (08 §6.2). Every failure up to the end
    // of the validation section below is request construction (Usage). ----
    let sample_names = opened.facts.sample_names.clone();
    let (gt_groups, group_names) = match &request.populations {
        Populations::File(path) => assign_groups(path, &sample_names, request.force_samples)?,
        Populations::Assignment(pairs) => {
            assign_groups_in_memory(pairs, &sample_names, request.force_samples)?
        }
        Populations::All => (vec![Some(0); sample_names.len()], vec!["all".to_string()]),
    };
    // No biallelic gate here: the stats reporters decide per statistic how
    // to treat multiallelic sites (doc_v3 05-roadmap R7); gating records
    // would silently change the summary counts.
    let filter: Option<Arc<reed_filter::plan::FilterPlan>> = match &request.site {
        Some(spec) => Some(Arc::new(reed_filter::plan::compile(
            spec,
            &reed_io::header::tag_table(&opened.facts.header),
            "",
            false,
        )?)),
        None => None,
    };
    let n_groups = group_names.len();

    request.grid.validate()?;
    let grid_mode = match &request.grid {
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
            if request.contig_lengths.is_empty() && !regions.is_empty() {
                // Region ends carry the clipping the Bp path takes from the
                // contig length; nothing else needs it.
            }
            GridMode::Regions(Arc::clone(regions))
        }
        WindowGrid::Sites { .. } => {
            return Err(crate::error::Error::usage(
                "site-count windows are not implemented (use --window/--step or --bed-file)",
            ));
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
    });
    let sfs_plan = match &request.sfs {
        None => None,
        Some(sfs) => {
            if sfs.projections.len() != n_groups {
                return Err(crate::error::Error::usage(format!(
                    "the SFS requests {} group projection(s) but the sample assignment declares {n_groups} group(s)",
                    sfs.projections.len()
                )));
            }
            for &n in &sfs.projections {
                if n < 2 {
                    return Err(crate::error::Error::usage(format!(
                        "projection size must be at least 2 chromosomes, got {n}"
                    )));
                }
            }
            let joint_pairs: Vec<(usize, usize)> = (0..n_groups)
                .flat_map(|a| ((a + 1)..n_groups).map(move |b| (a, b)))
                .collect();
            Some(SfsPlan {
                projections: sfs.projections.clone().into(),
                folded: sfs.folded,
                joint_pairs: joint_pairs.into(),
            })
        }
    };
    // `-s/-S/--keep/--remove`: excluded samples lose their group before the
    // group sizes are derived, so the denominators shrink with them.
    let mut gt_groups = gt_groups;
    if let Some(keep) = &request.keep_samples {
        if keep.len() != gt_groups.len() {
            return Err(crate::error::Error::usage(format!(
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
        let mut map: HashMap<Vec<u8>, u64> = request.contig_lengths.clone();
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
            .unwrap_or_else(|| vec![true; sample_names.len()]),
    );
    let job = StatsJob {
        grid: grid_mode.clone(),
        n_groups,
        gt_groups: gt_groups.into(),
        n_samples: sample_names.len(),
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
        sfs: sfs_plan.clone(),
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
        ld_merged: None,
        max_ld_jobs: budget.auxiliary(),
        ld_budget: *budget,
        ld: ld_plan,
        grid: grid_mode,
        seed: accs.clone(),
        levels: request.levels,
        skip_empty: request.skip_empty,
        layout: request.layout,        windows: Vec::new(),
        n_invariant_sites: 0,
        n_no_alt_records: 0,
        contig_len: (*contig_len).clone(),
        current: None,
        closed: Vec::new(),
        genome_totals: accs,
        genome_span: 0,
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

    Ok(StatsReport {
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
            joint: partial.joint.as_ref().map(|accs| {
                let plan = sfs_plan
                    .as_ref()
                    .expect("invariant: joint accs imply an SFS plan");
                JointSfsReport {
                    pairs: plan.joint_pairs.to_vec(),
                    spectra: accs.iter().map(JointProjectedAcc::values).collect(),
                    projections: accs.iter().map(|a| a.shape()).collect(),
                    folded: accs.first().is_some_and(JointProjectedAcc::folded),
                    n_sites: accs.iter().map(|a| a.site_tally().0).collect(),
                    n_skipped_m_lt_n: accs.iter().map(|a| a.site_tally().1).collect(),
                }
            }),
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
