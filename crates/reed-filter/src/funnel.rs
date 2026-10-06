//! Funnel and report v2 (§11.5): stable string `RuleId`s (fixing A24),
//! additive per-chunk counters, mask attribution, and the JSON report.

use serde::Serialize;

/// A stable rule identifier (`site.info.QD`, `after_mask.maf`, …), JSON-
/// serialized as a string.
pub type RuleId = String;

/// One funnel stage: the rule and the count of records that passed it
/// (cumulative semantics: a record that failed rule i counts in stages
/// 1..i only).
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct FunnelStage {
    /// The stable rule id.
    pub rule: RuleId,
    /// Records that reached and passed this rule.
    pub n_passed: u64,
}

/// Mask attribution (§11.5 `masked_by_rule`): the first rule that fires
/// owns the masked transition.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct MaskedBy {
    /// The mask rule id (`genotype.mask.DP`, in declaration order).
    pub rule: RuleId,
    /// Called-to-missing transitions attributed to this rule.
    pub n_masked: u64,
}

/// The run report (§11.5).
#[derive(Clone, Debug, Serialize)]
pub struct FilterReport {
    /// Records read by the final pass (the funnel input).
    pub n_input_records: u64,
    /// Records written.
    pub n_output_records: u64,
    /// Sample columns before the cohort rule.
    pub n_samples_in: u32,
    /// Sample columns after the cohort rule.
    pub n_samples_out: u32,
    /// Funnel stages in evaluation order.
    pub funnel: Vec<FunnelStage>,
    /// Total GT masked to `./.` (called-to-missing transitions only).
    pub n_masked_gt: u64,
    /// Mask attribution per rule.
    pub masked_by_rule: Vec<MaskedBy>,
    /// Policy triggers (ploidy/multiallelic rejections).
    pub policy_counters: PolicyCounters,
    /// The chosen strategy (§11.4: choice affects speed only).
    pub strategy: String,
    /// Wall time (seconds).
    pub elapsed_secs: f64,
}

/// Policy trigger counters (§11.5).
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct PolicyCounters {
    /// Records rejected for haploid/polyploid calls under
    /// `require_diploid`.
    pub ploidy_rejections: u64,
    /// Records rejected as multiallelic under `reject`.
    pub multiallelic_rejections: u64,
}

/// Additive funnel partial (one per chunk; summed in chunk order).
#[derive(Clone, Debug, Default)]
pub struct FunnelPartial {
    /// Per-head-rule pass counts (aligned with the plan's head list).
    pub head_hits: Vec<u64>,
    /// Per-tail-rule pass counts (aligned with the plan's tail list).
    pub tail_hits: Vec<u64>,
    /// Per-mask-rule attribution counts.
    pub mask_hits: Vec<u64>,
    /// Called-to-missing transitions.
    pub masked_gt: u64,
    /// Records read.
    pub n_input: u64,
    /// Records written.
    pub n_output: u64,
    /// Policy triggers.
    pub policy: PolicyCounters,
}

impl FunnelPartial {
    /// Re-arms the counters for a fresh chunk (the worker's funnel is
    /// `mem::take`n at chunk end, leaving empty vecs behind).
    pub fn reset(&mut self, n_head: usize, n_tail: usize, n_mask: usize) {
        self.head_hits = vec![0; n_head];
        self.tail_hits = vec![0; n_tail];
        self.mask_hits = vec![0; n_mask];
        self.masked_gt = 0;
        self.n_input = 0;
        self.n_output = 0;
        self.policy = PolicyCounters::default();
    }

    /// Adds another partial (component-wise, order-independent for the
    /// integer components).
    pub fn merge(&mut self, other: &Self) {
        // The receiving side may start with unsized (default) vectors;
        // extend so every counter of `other` lands somewhere.
        self.head_hits.resize(self.head_hits.len().max(other.head_hits.len()), 0);
        self.tail_hits.resize(self.tail_hits.len().max(other.tail_hits.len()), 0);
        self.mask_hits.resize(self.mask_hits.len().max(other.mask_hits.len()), 0);
        for (dst, src) in self.head_hits.iter_mut().zip(&other.head_hits) {
            *dst += src;
        }
        for (dst, src) in self.tail_hits.iter_mut().zip(&other.tail_hits) {
            *dst += src;
        }
        for (dst, src) in self.mask_hits.iter_mut().zip(&other.mask_hits) {
            *dst += src;
        }
        self.masked_gt += other.masked_gt;
        self.n_input += other.n_input;
        self.n_output += other.n_output;
        self.policy.ploidy_rejections += other.policy.ploidy_rejections;
        self.policy.multiallelic_rejections += other.policy.multiallelic_rejections;
    }
}
