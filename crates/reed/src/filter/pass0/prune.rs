//! The `--indep-pairwise` analysis over the unified pass0 matrices
//! (docs/rounds/v5/05-algorithms.md §3.4.1, doc_v6 03 §7): each contig's
//! [`SiteMatrix`] feeds one [`LdPruner`] in position order — plink2's
//! subcontig parallelism, strict plink2 order within a contig. The pruner
//! sees only the sites the tail rules kept: either the pass0 collected
//! them inline (no `--mind`) or the mid stage's survivor bitmap filters
//! here (`--mind` decides the keep set after the stream).

use rayon::prelude::*;

use reed_stats::ld::{LdPruner, LdWindow};

use super::matrix::SiteMatrix;
use reed_engine::budget::ThreadBudget;

/// The prune configuration.
#[derive(Clone, Copy, Debug)]
pub(crate) struct PruneConfig {
    pub threshold: f64,
    pub window: LdWindow,
    pub step: usize,
}

/// The kept sites, per contig, sorted by position (pass1 binary-searches
/// this: worker chunks are unordered, so an ordinal bitmap cannot be used).
///
/// Keyed by `(contig, pos0)`: multiple records at one position (rare;
/// e.g. a split multiallelic) share the prune decision — all kept or all
/// dropped (doc_v7 D7-25, the same position-keyed convention as
/// vcftools/plink2's `.prune.in`). A contig absent from `per_contig` counts
/// as fully pruned.
#[derive(Clone, Debug, Default)]
pub(crate) struct SiteKeepSet {
    /// `(contig, sorted 0-based positions)`.
    pub per_contig: Vec<(Vec<u8>, Vec<u64>)>,
}

impl SiteKeepSet {
    /// Whether the site at `(contig, pos0)` was kept.
    pub fn contains(&self, contig: &[u8], pos0: u64) -> bool {
        self.per_contig
            .iter()
            .filter(|(c, _)| c == contig)
            .any(|(_, positions)| positions.binary_search(&pos0).is_ok())
    }
}

/// The prune result: the kept sites and the counters.
#[derive(Clone, Debug, Default)]
pub(crate) struct PruneOutcome {
    /// The kept sites.
    pub keep: SiteKeepSet,
    /// Sites the pruner saw (the post-tail survivors).
    pub n_seen: u64,
    /// Sites the pruner removed.
    pub n_pruned: u64,
    /// Wall time of the parallel prune (the driver's timing report).
    pub prune_secs: f64,
}

/// One site's plink2 aggregates over all samples, from its four planes:
/// the major-allele frequency and the hom-ref / hom-alt / het counts
/// (`FillVaggs`).
/// One site's plink2 aggregates over the kept, unmasked samples, from its
/// four planes plus the per-sample exclusion bitmap (doc_v7 D7-4): the
/// major-allele frequency and the hom-ref / hom-alt / het counts
/// (`FillVaggs`).
pub(crate) fn site_prune_stats_ex(
    planes: [&[u64]; 4],
    n_samples: usize,
    exclude: &[u64],
) -> reed_stats::ld::SitePruneStats {
    let mut alt = 0_u64;
    let mut called = 0_u64;
    let mut n_hom_ref = 0_u64;
    let mut n_hom_alt = 0_u64;
    let mut n_het = 0_u64;
    for (w, _) in planes[0].iter().enumerate() {
        let real = if (w + 1) * 64 <= n_samples {
            !0_u64
        } else {
            (1_u64 << (n_samples % 64)) - 1
        };
        let excl = exclude.get(w).copied().unwrap_or(0);
        let miss = (planes[1][w] | planes[3][w] | excl) & real;
        let called_w = real & !miss;
        let alt0 = planes[0][w] & called_w;
        let alt1 = planes[2][w] & called_w;
        alt += alt0.count_ones() as u64 + alt1.count_ones() as u64;
        called += 2 * called_w.count_ones() as u64;
        let hom = called_w & !(planes[1][w] | planes[3][w]) & !(alt0 ^ alt1);
        n_hom_alt += (hom & alt0).count_ones() as u64;
        n_hom_ref += (hom & !alt0).count_ones() as u64;
        n_het += (called_w & !(planes[1][w] | planes[3][w]) & (alt0 ^ alt1)).count_ones() as u64;
    }
    let maj_freq = if called == 0 {
        1.0
    } else {
        let alt_freq = alt as f64 / called as f64;
        alt_freq.max(1.0 - alt_freq)
    };
    reed_stats::ld::SitePruneStats {
        maj_freq,
        n_hom_ref,
        n_hom_alt,
        n_het,
    }
}

/// Prunes one contig's matrix; returns `(contig, kept 0-based positions,
/// sites seen, sites pruned)`.
///
/// # Panics
///
/// Panics only if the matrix is corrupt (shorter rows than the stride
/// declares), which `push_block` maintains by construction.
pub(crate) fn prune_one(
    matrix: SiteMatrix,
    config: &PruneConfig,
    survivor: Option<&[u64]>,
    base: usize,
    n_samples: usize,
    keep_words: &[u64],
) -> PruneOutcome {
    let w = matrix.w;
    let mut pruner = LdPruner::new(config.threshold, config.window, config.step, w, n_samples);
    let mut n_seen = 0_u64;
    for (i, &pos0) in matrix.pos0().iter().enumerate() {
        if let Some(words) = survivor
            && !survivor_bit(words, base + i)
        {
            continue;
        }
        n_seen += 1;
        let planes = matrix.planes(i);
        // doc_v7 D7-4: the pruner sees the post-mask, kept-only genotypes.
        // The site's mask plane carries the fired mask; the keep bitmap
        // removes the samples the static subset (and `--mind`) dropped.
        let mask = matrix.mask_plane(i);
        let mut exclude = vec![0_u64; w];
        for (wi, e) in exclude.iter_mut().enumerate() {
            let keep = if keep_words.is_empty() {
                !0_u64
            } else {
                keep_words.get(wi).copied().unwrap_or(0)
            };
            *e = mask.get(wi).copied().unwrap_or(0) | !keep;
        }
        let stats = site_prune_stats_ex(planes, n_samples, &exclude);
        pruner.push_ex(i, pos0, planes, stats, exclude);
    }
    pruner.finish_contig();
    // A later pair can prune an earlier site, so `push`'s return is
    // provisional; `kept_sites()` after `finish_contig` is the answer.
    let pos = matrix.pos0();
    let kept: Vec<u64> = pruner.kept_sites().into_iter().map(|o| pos[o]).collect();
    let n_pruned = n_seen - kept.len() as u64;
    PruneOutcome {
        keep: SiteKeepSet {
            per_contig: vec![(matrix.contig.clone(), kept)],
        },
        n_seen,
        n_pruned,
        prune_secs: 0.0,
    }
}

/// Whether global site ordinal `o` survived the mid stage.
fn survivor_bit(words: &[u64], o: usize) -> bool {
    words
        .get(o / 64)
        .is_some_and(|word| word >> (o % 64) & 1 == 1)
}

/// Prunes every contig in parallel (plink2's subcontig parallelism) and
/// merges the outcomes in contig order — the deterministic result the
/// per-contig independence guarantees. The pool comes from the budget,
/// never the rayon global pool (doc_v8 D8-1).
pub(crate) fn prune_matrices(
    matrices: Vec<SiteMatrix>,
    survivor: Option<&[u64]>,
    bases: &[usize],
    config: &PruneConfig,
    n_samples: usize,
    keep_words: &[u64],
    budget: &ThreadBudget,
) -> PruneOutcome {
    let started = std::time::Instant::now();
    let per_contig: Vec<PruneOutcome> = reed_engine::budget::run_auxiliary(budget, |parallel| {
        if parallel {
            matrices
                .into_par_iter()
                .enumerate()
                .map(|(index, matrix)| {
                    let base = bases.get(index).copied().unwrap_or(0);
                    prune_one(matrix, config, survivor, base, n_samples, keep_words)
                })
                .collect()
        } else {
            matrices
                .into_iter()
                .enumerate()
                .map(|(index, matrix)| {
                    let base = bases.get(index).copied().unwrap_or(0);
                    prune_one(matrix, config, survivor, base, n_samples, keep_words)
                })
                .collect()
        }
    })
    .expect("invariant: the auxiliary pool builds");
    let mut outcome = PruneOutcome::default();
    for part in per_contig {
        let PruneOutcome {
            keep: mut site_keep,
            n_seen,
            n_pruned,
            ..
        } = part;
        outcome.n_seen += n_seen;
        outcome.n_pruned += n_pruned;
        outcome.keep.per_contig.append(&mut site_keep.per_contig);
    }
    outcome.prune_secs = started.elapsed().as_secs_f64();
    outcome
}
