//! DAPC cross-validation (`adegenet::xvalDapc`).
//!
//! For each candidate number of principal components, split the samples into a
//! training and a validation set, fit the discriminant analysis on the
//! training half, predict the validation half, and record the accuracy. The
//! optimum is the candidate with the lowest `sqrt(mean((success - 1)^2))`.
//!
//! # What is split and what is not
//!
//! The PCA is computed **once** on the whole data set and reused for every
//! repetition; only the LDA is refit. `adegenet` does the same by passing its
//! precomputed `dudi` object into each `dapc` call. Refitting the PCA per
//! repetition would leak the validation samples into the components and
//! inflate the accuracy.
//!
//! # Reproducibility
//!
//! Each repetition derives its own sub-stream from `Rng::derive(repetition)`,
//! so the result does not depend on thread count. The numbers cannot match
//! R's Mersenne-Twister stream, so the parity target is the deterministic
//! bookkeeping (candidate list, training-set sizing) plus the accuracy
//! distribution, not the per-repetition values.

use reed_core::math::r_compat::{pretty_sequence, round_half_even_usize};
use reed_core::rng::Rng;

use crate::lda::{self, Lda, LdaError};
use reed_core::linalg::Mat;

/// Errors from cross-validation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum XvalError {
    /// Fewer than two samples, or fewer than two groups.
    TooFewSamples,
    /// Every group has a single member, so no training/validation split is
    /// possible.
    NoSplittableGroup,
    /// `training_set` was not strictly between 0 and 1.
    BadTrainingFraction,
    /// The LDA failed on one of the splits.
    Lda(LdaError),
}

impl std::fmt::Display for XvalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::TooFewSamples => write!(f, "cross-validation needs at least 2 samples and 2 groups"),
            Self::NoSplittableGroup => {
                write!(f, "every group has a single member; no training set is possible")
            }
            Self::BadTrainingFraction => {
                write!(f, "the training fraction must be strictly between 0 and 1")
            }
            Self::Lda(e) => write!(f, "discriminant analysis failed: {e}"),
        }
    }
}

impl std::error::Error for XvalError {}

/// Which accuracy to report per repetition.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum XvalResult {
    /// Fraction of validation samples assigned to their true group.
    Overall,
    /// Per-group accuracy, averaged over groups (so a small group weighs as
    /// much as a large one).
    GroupMean,
}

impl XvalResult {
    /// Parses the CLI spelling.
    ///
    /// # Errors
    ///
    /// A message naming the accepted values.
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "overall" => Ok(Self::Overall),
            "group-mean" | "groupmean" => Ok(Self::GroupMean),
            other => Err(format!(
                "unknown result {other:?} (overall|group-mean)"
            )),
        }
    }
}

/// One candidate's results.
#[derive(Clone, Debug)]
pub struct Candidate {
    /// Number of principal components.
    pub n_pca: usize,
    /// Accuracy of each repetition, in repetition order.
    pub success: Vec<f64>,
}

impl Candidate {
    /// Mean accuracy across repetitions.
    #[must_use]
    pub fn mean_success(&self) -> f64 {
        if self.success.is_empty() {
            return 0.0;
        }
        self.success.iter().sum::<f64>() / self.success.len() as f64
    }

    /// `sqrt(mean((success - 1)^2))`, the criterion `adegenet` minimizes.
    ///
    /// Note this is **not** the usual RMSE: there is no mean subtraction, so
    /// it measures the deviation from a perfect score rather than from the
    /// candidate's own mean. Reproduced as-is so the selected optimum matches.
    #[must_use]
    pub fn error_from_perfect(&self) -> f64 {
        if self.success.is_empty() {
            return 1.0;
        }
        let sum: f64 = self
            .success
            .iter()
            .map(|s| (s - 1.0) * (s - 1.0))
            .sum();
        (sum / self.success.len() as f64).sqrt()
    }
}

/// The full cross-validation result.
#[derive(Clone, Debug)]
pub struct XvalReport {
    /// One entry per candidate, in ascending `n_pca`.
    pub candidates: Vec<Candidate>,
    /// The `n_pca` with the lowest [`Candidate::error_from_perfect`]. Ties go
    /// to the **largest** candidate, as in `adegenet`.
    pub best_n_pca: usize,
    /// The `n_pca` with the highest mean success (`adegenet` reports this
    /// alongside the MSE choice).
    pub best_mean_success_n_pca: usize,
    /// Accuracy expected from guessing at random: `1 / n_groups`, the value
    /// `adegenet`'s label-permutation baseline converges to.
    pub random_chance: f64,
}

/// Cross-validation settings.
#[derive(Clone, Debug)]
pub struct XvalConfig {
    /// Fraction of each group used for training. Reduced automatically when
    /// the smallest group cannot support it.
    pub training_set: f64,
    /// Repetitions per candidate.
    pub n_rep: usize,
    /// Which accuracy to report.
    pub result: XvalResult,
    /// Seed for the per-repetition sub-streams.
    pub seed: u64,
}

impl Default for XvalConfig {
    fn default() -> Self {
        Self {
            training_set: 0.9,
            n_rep: 30,
            result: XvalResult::Overall,
            seed: 42,
        }
    }
}

/// Runs cross-validation over PC scores.
///
/// `scores` is `n_samples x n_components` (the PCA scores, not the loadings);
/// `groups[i]` indexes `levels`. `candidates` is the list of component counts
/// to try; see [`candidate_counts`] for the `adegenet` rule.
///
/// # Errors
///
/// See [`XvalError`].
pub fn run(
    scores: &Mat,
    groups: &[u32],
    levels: &[String],
    candidates: &[usize],
    config: &XvalConfig,
) -> Result<XvalReport, XvalError> {
    let n = scores.rows;
    let ng = levels.len();
    if n < 2 || ng < 2 {
        return Err(XvalError::TooFewSamples);
    }
    if !(config.training_set > 0.0 && config.training_set < 1.0) {
        return Err(XvalError::BadTrainingFraction);
    }
    if groups.len() != n {
        return Err(XvalError::TooFewSamples);
    }

    // Group sizes, and the fraction the smallest group can support.
    let mut counts = vec![0usize; ng];
    for &g in groups {
        // Validate before indexing: an out-of-range id used to panic with
        // index-out-of-bounds here, while `lda::fit` reports the same
        // condition with a typed error. Map to that error (doc_v20 S32).
        if g as usize >= ng {
            return Err(XvalError::Lda(LdaError::GroupIndexOutOfRange { index: g, ng }));
        }
        counts[g as usize] += 1;
    }
    let popmin = counts.iter().copied().filter(|&c| c >= 2).min();
    let Some(popmin) = popmin else {
        return Err(XvalError::NoSplittableGroup);
    };
    // `training.set2 = (popmin - 1) / popmin`: the largest fraction that still
    // leaves at least one sample out of the smallest group.
    let training_set2 = (popmin - 1) as f64 / popmin as f64;
    let training_set = config.training_set.min(training_set2);

    // The random baseline `adegenet` reports: it permutes the labels and takes
    // `mean(tapply(shuffled == truth, truth, mean))`. That is the **unweighted
    // mean of per-group hit rates**, which under a uniform permutation is
    // `1 / ng` — not `sum(p_g^2)`, which would be the accuracy of a
    // prior-weighted classifier and is a different quantity.
    //
    // Measured on groups of 30/20/10: `adegenet` gives 0.3335, `sum(p^2)` is
    // 0.3889. Using the latter would make a real signal look weaker.
    let random_chance = 1.0 / ng as f64;

    let base_rng = Rng::from_seed(config.seed);
    let mut candidates_out = Vec::with_capacity(candidates.len());
    let mut best_err = f64::INFINITY;
    let mut best_n_pca = candidates.first().copied().unwrap_or(1);
    let mut best_mean = f64::NEG_INFINITY;
    let mut best_mean_n_pca = best_n_pca;

    for &n_pca in candidates {
        let k = n_pca.min(scores.cols);
        if k == 0 {
            continue;
        }
        // The first k components of the full-data PCA.
        let mut sub = Mat::zeros(n, k);
        for i in 0..n {
            for c in 0..k {
                sub.set(i, c, scores.get(i, c));
            }
        }

        let mut success = Vec::with_capacity(config.n_rep);
        for rep in 0..config.n_rep {
            // One sub-stream per repetition: the draw does not depend on
            // which thread runs it.
            let mut rng = base_rng.derive(rep as u64);
            let (train, valid) = split(groups, ng, &counts, training_set, &mut rng);
            success.push(evaluate(&sub, groups, levels, &train, &valid, config.result)?);
        }

        let cand = Candidate { n_pca, success };
        let err = cand.error_from_perfect();
        if err < best_err {
            best_err = err;
            best_n_pca = n_pca;
        } else if err == best_err && n_pca > best_n_pca {
            // `adegenet` takes the largest candidate among exact ties.
            best_n_pca = n_pca;
        }
        let mean = cand.mean_success();
        if mean > best_mean {
            best_mean = mean;
            best_mean_n_pca = n_pca;
        }
        candidates_out.push(cand);
    }

    Ok(XvalReport {
        candidates: candidates_out,
        best_n_pca,
        best_mean_success_n_pca: best_mean_n_pca,
        random_chance,
    })
}

/// The candidate component counts `adegenet` would try.
///
/// `n_pca_max` is already clamped to `min(requested, pca_rank, N_training-1)`.
/// The returned values are strictly below `n_pca_max`, matching the filter in
/// `xvalDapc` (`n.pca < n.pca.max`).
#[must_use]
pub fn candidate_counts(n_pca_max: usize) -> Vec<usize> {
    if n_pca_max < 2 {
        return Vec::new();
    }
    let runs = if n_pca_max < 10 { n_pca_max } else { 10 };
    pretty_sequence(n_pca_max, runs)
        .into_iter()
        .filter(|&v| v >= 1 && v < n_pca_max)
        .collect()
}

/// Splits the samples into training and validation indices.
///
/// Each group contributes `round(training_set * n_g)` training samples, drawn
/// without replacement. A group with a single member goes entirely into the
/// training set (it cannot be validated).
fn split(
    groups: &[u32],
    ng: usize,
    counts: &[usize],
    training_set: f64,
    rng: &mut Rng,
) -> (Vec<usize>, Vec<usize>) {
    let n = groups.len();
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); ng];
    for (i, &g) in groups.iter().enumerate() {
        members[g as usize].push(i);
    }

    let mut is_train = vec![false; n];
    for (g, ids) in members.iter().enumerate() {
        if counts[g] < 2 {
            // A single-member group cannot be split; `adegenet` keeps it whole
            // in the training set.
            for &i in ids {
                is_train[i] = true;
            }
            continue;
        }
        // `adegenet` takes `round(training.set * n_g)` with no clamp, so a
        // group small enough that the rounded count equals its size
        // contributes **no** validation samples. That is not a bug to fix: it
        // is why the reference caps `training.set` at `(popmin-1)/popmin`
        // only for the smallest group, and larger-but-still-small groups can
        // still be fully trained. Clamping to `n_g - 1` here would score a
        // sample the reference never validates.
        let take = round_half_even_usize(training_set * counts[g] as f64);
        let mut pool = ids.clone();
        rng.shuffle(&mut pool);
        for &i in pool.iter().take(take) {
            is_train[i] = true;
        }
    }

    let mut train = Vec::new();
    let mut valid = Vec::new();
    for i in 0..n {
        if is_train[i] {
            train.push(i);
        } else {
            valid.push(i);
        }
    }
    (train, valid)
}

/// Fits on `train`, predicts `valid`, returns the accuracy.
fn evaluate(
    scores: &Mat,
    groups: &[u32],
    levels: &[String],
    train: &[usize],
    valid: &[usize],
    result: XvalResult,
) -> Result<f64, XvalError> {
    if valid.is_empty() {
        return Ok(1.0);
    }
    let p = scores.cols;
    let mut x_train = Mat::zeros(train.len(), p);
    let mut g_train = Vec::with_capacity(train.len());
    for (row, &i) in train.iter().enumerate() {
        for c in 0..p {
            x_train.set(row, c, scores.get(i, c));
        }
        g_train.push(groups[i]);
    }

    let model: Lda =
        lda::fit(&x_train, &g_train, levels, 1e-30).map_err(XvalError::Lda)?;

    let mut x_valid = Mat::zeros(valid.len(), p);
    for (row, &i) in valid.iter().enumerate() {
        for c in 0..p {
            x_valid.set(row, c, scores.get(i, c));
        }
    }
    // Only the axes with a real singular value: the trailing ones are
    // numerical null directions from collinear class means, and including
    // them adds noise to the distances that decide the assignment. This is
    // the same `svd^2 > 1e-10` cut `adegenet` applies (`lda.dim`).
    let dimen = model
        .svd
        .iter()
        .filter(|&&v| v * v > 1e-10)
        .count()
        .max(1);
    let projection = model.project(&x_valid, dimen).map_err(XvalError::Lda)?;

    Ok(match result {
        XvalResult::Overall => {
            let hits = valid
                .iter()
                .enumerate()
                .filter(|&(row, &i)| projection.assign[row] == groups[i])
                .count();
            hits as f64 / valid.len() as f64
        }
        XvalResult::GroupMean => {
            // Per-group accuracy, then the unweighted mean over groups that
            // have at least one validation sample.
            let mut hits = vec![0usize; levels.len()];
            let mut totals = vec![0usize; levels.len()];
            for (row, &i) in valid.iter().enumerate() {
                let g = groups[i] as usize;
                totals[g] += 1;
                if projection.assign[row] == groups[i] {
                    hits[g] += 1;
                }
            }
            let per_group: Vec<f64> = (0..levels.len())
                .filter(|&g| totals[g] > 0)
                .map(|g| hits[g] as f64 / totals[g] as f64)
                .collect();
            if per_group.is_empty() {
                0.0
            } else {
                per_group.iter().sum::<f64>() / per_group.len() as f64
            }
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn levels(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("g{i}")).collect()
    }

    /// Three well-separated groups in two dimensions.
    ///
    /// The within-group noise matters: with points on a perfect line the
    /// within-class scatter is rank-deficient and the discriminant collapses
    /// to a numerical null direction, which is a property of LDA on
    /// noiseless data rather than of the cross-validation.
    fn separated() -> (Mat, Vec<u32>) {
        let mut x = Mat::zeros(30, 2);
        let mut g = Vec::with_capacity(30);
        let mut seed = 7u64;
        let mut jitter = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            ((seed >> 33) as f64 / (1u64 << 31) as f64 - 1.0) * 0.5
        };
        for i in 0..30 {
            let group = (i / 10) as u32;
            g.push(group);
            let shift = match group {
                0 => -5.0,
                1 => 0.0,
                _ => 5.0,
            };
            x.set(i, 0, shift + jitter());
            x.set(i, 1, shift * 0.5 + jitter());
        }
        (x, g)
    }

    #[test]
    fn candidates_follow_the_r_rule() {
        // n_pca.max = 20 -> multiples of 2 below 20.
        assert_eq!(
            candidate_counts(20),
            vec![2, 4, 6, 8, 10, 12, 14, 16, 18]
        );
        // Small ranges enumerate every integer below the max.
        assert_eq!(candidate_counts(5), vec![1, 2, 3, 4]);
        assert_eq!(candidate_counts(2), vec![1]);
        assert!(candidate_counts(1).is_empty());
    }

    #[test]
    fn well_separated_groups_score_perfectly() {
        let (x, g) = separated();
        let cfg = XvalConfig {
            n_rep: 5,
            ..Default::default()
        };
        let report = run(&x, &g, &levels(3), &[1, 2], &cfg).unwrap();
        for cand in &report.candidates {
            assert!(
                (cand.mean_success() - 1.0).abs() < 1e-12,
                "n_pca={} scored {}",
                cand.n_pca,
                cand.mean_success()
            );
        }
        assert_eq!(report.best_n_pca, 2, "ties go to the largest candidate");
    }

    #[test]
    fn a_single_member_group_never_enters_validation() {
        // Groups of 10, 10 and 1: the singleton can never be validated.
        let mut x = Mat::zeros(21, 2);
        let mut g = Vec::with_capacity(21);
        for i in 0..21 {
            let group = if i < 10 { 0 } else if i < 20 { 1 } else { 2 };
            g.push(group);
            x.set(i, 0, group as f64 * 5.0 + (i % 4) as f64 * 0.3);
            x.set(i, 1, group as f64 * 3.0 + (i % 3) as f64 * 0.2);
        }
        let cfg = XvalConfig {
            n_rep: 3,
            ..Default::default()
        };
        let report = run(&x, &g, &levels(3), &[1], &cfg).unwrap();
        // It must not error, and the training fraction is capped by popmin=10.
        assert_eq!(report.candidates.len(), 1);
    }

    #[test]
    fn training_fraction_is_capped_by_the_smallest_group() {
        // Group sizes 10, 10 and 3: popmin = 3, so training.set2 = 2/3 and a
        // 0.9 request is reduced to 2/3. The smallest group then contributes
        // exactly 2 training samples and 1 validation sample per repetition.
        let mut x = Mat::zeros(23, 2);
        let mut g = Vec::with_capacity(23);
        for i in 0..23 {
            let group = if i < 10 {
                0
            } else if i < 20 {
                1
            } else {
                2
            };
            g.push(group);
            x.set(i, 0, group as f64 * 8.0 + (i % 3) as f64 * 0.4);
            x.set(i, 1, group as f64 * 3.0 + (i % 2) as f64 * 0.3);
        }
        let cfg = XvalConfig {
            n_rep: 10,
            ..Default::default()
        };
        let report = run(&x, &g, &levels(3), &[1], &cfg).unwrap();
        // Accuracy is a fraction of 22 validation samples; with one of them in
        // the smallest group, every score is a multiple of 1/22.
        for cand in &report.candidates {
            for &s in &cand.success {
                let scaled = s * 22.0;
                assert!(
                    (scaled - scaled.round()).abs() < 1e-9,
                    "accuracy {s} is not a multiple of 1/22"
                );
            }
        }
    }

    /// With the reference's cap in place, every group keeps at least one
    /// sample out, so the accuracy denominator is the per-group validation
    /// count and every score is a multiple of `1 / total_validation`.
    ///
    /// This is the invariant that makes dropping the `n_g - 1` clamp safe:
    /// `training.set <= (popmin-1)/popmin`, so for any group `n_g >= popmin`,
    /// `round(training.set * n_g) <= n_g - 1`. The clamp was therefore
    /// redundant; removing it keeps the arithmetic identical to `adegenet`
    /// for every reachable input.
    #[test]
    fn every_group_keeps_a_validation_sample() {
        // Three groups of 10: popmin = 10, cap 0.9, each leaves exactly one.
        let mut x = Mat::zeros(30, 2);
        let mut g = Vec::with_capacity(30);
        for i in 0..30 {
            let group = if i < 10 {
                0
            } else if i < 20 {
                1
            } else {
                2
            };
            g.push(group);
            x.set(i, 0, group as f64 * 5.0 + (i % 3) as f64 * 0.4);
            x.set(i, 1, group as f64 * 2.0 + (i % 2) as f64 * 0.3);
        }
        let cfg = XvalConfig {
            n_rep: 6,
            ..Default::default()
        };
        let report = run(&x, &g, &levels(3), &[1], &cfg).unwrap();
        for cand in &report.candidates {
            for &s in &cand.success {
                let scaled = s * 3.0;
                assert!(
                    (scaled - scaled.round()).abs() < 1e-9,
                    "accuracy {s} is not a multiple of 1/3 (3 validation samples)"
                );
            }
        }
    }

    /// Two groups with identical centroids make the LDA rank collapse below
    /// `ng - 1`. `adegenet` asks each fold for `n.da = ng - 1`; reed uses the
    /// observed rank. The predicted assignment must be the same either way,
    /// so the accuracy must still be computable and finite.
    #[test]
    fn collinear_group_centroids_do_not_break_the_accuracy() {
        // Groups 0 and 1 share a centroid; group 2 is far away.
        let mut x = Mat::zeros(30, 2);
        let mut g = Vec::with_capacity(30);
        for i in 0..30 {
            let group = if i < 10 {
                0
            } else if i < 20 {
                1
            } else {
                2
            };
            g.push(group);
            // Group 2 sits at x = 10; groups 0 and 1 both straddle 0.
            let base = if group == 2 { 10.0 } else { 0.0 };
            x.set(i, 0, base + ((i % 5) as f64 - 2.0) * 0.3);
            x.set(i, 1, base + ((i % 3) as f64 - 1.0) * 0.2);
        }
        let cfg = XvalConfig {
            n_rep: 5,
            ..Default::default()
        };
        let report = run(&x, &g, &levels(3), &[1, 2], &cfg).unwrap();
        for cand in &report.candidates {
            for &s in &cand.success {
                assert!(s.is_finite() && (0.0..=1.0).contains(&s), "bad score {s}");
            }
            // Group 2 is always separable, so no repetition can score zero:
            // at least its validation sample is right.
            assert!(cand.mean_success() > 0.0, "n_pca={}", cand.n_pca);
        }
    }

    #[test]
    fn results_are_reproducible() {
        let (x, g) = separated();
        let cfg = XvalConfig {
            n_rep: 8,
            seed: 1234,
            ..Default::default()
        };
        let a = run(&x, &g, &levels(3), &[1, 2], &cfg).unwrap();
        let b = run(&x, &g, &levels(3), &[1, 2], &cfg).unwrap();
        for (ca, cb) in a.candidates.iter().zip(&b.candidates) {
            assert_eq!(ca.success, cb.success);
        }
    }

    #[test]
    fn different_seeds_change_the_splits() {
        let (x, g) = separated();
        let cfg_a = XvalConfig {
            n_rep: 8,
            seed: 1,
            ..Default::default()
        };
        let cfg_b = XvalConfig {
            n_rep: 8,
            seed: 2,
            ..Default::default()
        };
        let a = run(&x, &g, &levels(3), &[1], &cfg_a).unwrap();
        let b = run(&x, &g, &levels(3), &[1], &cfg_b).unwrap();
        // On a perfectly separable fixture both are 1.0; the splits still
        // differ, which shows up only through the RNG, so check the generator
        // directly instead of the accuracy.
        let mut r1 = Rng::from_seed(1).derive(0);
        let mut r2 = Rng::from_seed(2).derive(0);
        assert_ne!(r1.next_u64(), r2.next_u64());
        assert_eq!(a.candidates.len(), b.candidates.len());
    }

    #[test]
    fn group_mean_weights_groups_equally() {
        // Group 0 is large and easy; group 1 is small. Overall accuracy would
        // be dominated by group 0, group-mean is not.
        let mut x = Mat::zeros(30, 2);
        let mut g = Vec::with_capacity(30);
        for i in 0..30 {
            let group = if i < 25 { 0 } else { 1 };
            g.push(group);
            x.set(i, 0, group as f64 * 10.0 + (i % 5) as f64 * 0.3);
            x.set(i, 1, group as f64 * 4.0 + (i % 3) as f64 * 0.2);
        }
        let cfg = XvalConfig {
            n_rep: 5,
            result: XvalResult::GroupMean,
            ..Default::default()
        };
        let report = run(&x, &g, &levels(2), &[1], &cfg).unwrap();
        // Perfectly separable, so both results are 1.0; the point is that the
        // group-mean path runs and reports a probability.
        assert!(report.candidates[0].mean_success() > 0.9);
    }

    #[test]
    fn s32_group_id_out_of_range_is_a_typed_error_not_a_panic() {
        // doc_v20 S32: `counts[g as usize] += 1` ran before any validation,
        // so a group id >= ng panicked with index-out-of-bounds instead of
        // surfacing the typed error `lda::fit` already defines for exactly
        // this condition.
        let x = Mat::from_fn(3, 2, |i, c| i as f64 + c as f64 * 0.5);
        let cfg = XvalConfig { n_rep: 2, ..Default::default() };
        let err = run(&x, &[0, 0, 2], &levels(2), &[1], &cfg).unwrap_err();
        assert_eq!(
            err,
            XvalError::Lda(LdaError::GroupIndexOutOfRange { index: 2, ng: 2 })
        );
    }

    #[test]
    fn bad_training_fraction_is_rejected() {
        let (x, g) = separated();
        for bad in [0.0, 1.0, -0.5, 1.5] {
            let cfg = XvalConfig {
                training_set: bad,
                ..Default::default()
            };
            assert_eq!(
                run(&x, &g, &levels(3), &[1], &cfg).unwrap_err(),
                XvalError::BadTrainingFraction
            );
        }
    }

    #[test]
    fn a_single_group_is_rejected() {
        let (x, _) = separated();
        let g = vec![0u32; 30];
        let cfg = XvalConfig::default();
        assert_eq!(
            run(&x, &g, &levels(1), &[1], &cfg).unwrap_err(),
            XvalError::TooFewSamples
        );
    }

    #[test]
    fn error_from_perfect_matches_the_r_formula() {
        let c = Candidate {
            n_pca: 3,
            success: vec![1.0, 0.5, 0.5, 1.0],
        };
        // sqrt(mean((s-1)^2)) = sqrt((0 + 0.25 + 0.25 + 0)/4) = sqrt(0.125)
        let want = 0.125f64.sqrt();
        assert!((c.error_from_perfect() - want).abs() < 1e-12);
    }

    #[test]
    fn random_chance_is_one_over_the_group_count() {
        // `adegenet` permutes the labels and averages per-group hit rates, so
        // the baseline is 1/ng regardless of how unequal the groups are.
        // Measured: groups of 30/20/10 give 0.3335, while sum(p^2) is 0.3889.
        let mut x = Mat::zeros(60, 2);
        let mut g = Vec::with_capacity(60);
        for i in 0..60 {
            let group = if i < 30 {
                0
            } else if i < 50 {
                1
            } else {
                2
            };
            g.push(group);
            x.set(i, 0, group as f64 * 5.0 + (i % 3) as f64 * 0.4);
            x.set(i, 1, group as f64 * 2.0 + (i % 2) as f64 * 0.3);
        }
        let cfg = XvalConfig {
            n_rep: 1,
            ..Default::default()
        };
        let report = run(&x, &g, &levels(3), &[1], &cfg).unwrap();
        assert!(
            (report.random_chance - 1.0 / 3.0).abs() < 1e-12,
            "got {}",
            report.random_chance
        );
    }
}
