//! Funnel and report (§11.5): stable string `RuleId`s (fixing A24),
//! additive per-chunk counters, mask attribution, and the JSON report.

use serde::Serialize;

/// A stable rule identifier (`site.info.QD`, `after_mask.maf`, …), JSON-
/// serialized as a string.
pub type RuleId = String;

/// The site-level missing-rate rule's funnel id (`--geno`).
pub const AFTER_MASK_MISSING: &str = "after_mask.missing";
/// The sample-level missing-rate rule's id (`--mind`). The funnel stages
/// count records, so this rule is reported through
/// [`FilterReport::sample_rule`] instead of posing as a stage
/// (doc_v4 02-cli-semantics.md 回退版 §5.3).
pub const SAMPLE_MIND: &str = "sample.mind";
/// The LD-pruning rule's funnel id (`--indep-pairwise`,
/// docs/rounds/v5/05-algorithms.md §3.4.1).
pub const AFTER_MASK_LD_PRUNE: &str = "after_mask.ld_prune";

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

/// Phase wall times plus the final pass's engine metrics (`--timing`, doc_v3
/// O10; doc_v4 03 A4-5 extends the engine side). The CLI renders this for
/// `--timing FILE` and for `--verbosity 2`, so the driver never prints.
#[derive(Clone, Debug, Default, Serialize)]
pub struct TimingReport {
    /// Worker threads the run used.
    pub threads: usize,
    /// `(phase, seconds)` in signal order (open, pass0, mid, pass1).
    pub phases: Vec<(String, f64)>,
    /// Whole-run wall time.
    pub total_secs: f64,
    /// Chunks the final pass consumed.
    pub chunks: u64,
    /// Chunks read across all passes (doc_v4 03 A4-5).
    pub chunks_read: u64,
    /// Summed worker-busy time across all passes.
    pub worker_busy_secs: f64,
    /// Cumulative reader time (`next_chunk`) across all passes.
    pub reader_secs: f64,
    /// Cumulative reducer time (`accept`) across all passes.
    pub reduce_secs: f64,
    /// Cumulative time the reader spent blocked on tickets/permits/send (A6-4).
    pub reader_wait_secs: f64,
    /// Cumulative time the reducer spent blocked on the result channel (A6-4).
    pub reducer_wait_secs: f64,
    /// Cumulative worker inflate time across all passes (A6-4).
    pub inflate_secs: f64,
    /// Cumulative worker job time across all passes (deflate included; A6-4).
    pub job_secs: f64,
    /// Cumulative output-compression time reported by the sink (A6-4).
    pub deflate_secs: f64,
    /// Cumulative time workers waited on the BCF frame chain.
    pub chain_wait_secs: f64,
    /// Peak in-flight bytes, when the engine reported it.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub in_flight_peak_bytes: Option<u64>,
}

/// The dynamic sample rule's outcome (`sample.mind`, doc_v4
/// 02-cli-semantics.md §5.3). The funnel stages count records, so the
/// sample-level rule reports separately instead of posing as a stage.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct SampleRuleReport {
    /// The stable rule id (`sample.mind`).
    pub rule: RuleId,
    /// The threshold as rendered (exact decimal).
    pub threshold: String,
    /// Samples the rule dropped.
    pub n_dropped: u64,
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
    /// The `--max-missing-sample` outcome when the rule was requested.
    pub sample_rule: Option<SampleRuleReport>,
    /// `--indep-pairwise`: `(survivors, pruned)` when the rule was
    /// requested (docs/rounds/v5/05-algorithms.md §3.4.1).
    pub ld_prune: Option<(u64, u64)>,
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
    /// Phase breakdown for `--timing` / `--verbosity 2`.
    pub timing: TimingReport,
}

impl FilterReport {
    /// The `--timing` text (doc_v3 O10): one `key<TAB>value` line per
    /// counter, phases in signal order.
    #[must_use]
    pub fn timing_text(&self) -> String {
        use std::fmt::Write as _;
        let t = &self.timing;
        let mut text = String::from("# reed filter timing (03-performance.md O10)\n");
        let _ = writeln!(text, "strategy\t{}", self.strategy);
        let _ = writeln!(text, "threads\t{}", t.threads);
        let _ = writeln!(text, "n_records_in\t{}", self.n_input_records);
        let _ = writeln!(text, "n_records_out\t{}", self.n_output_records);
        for (phase, secs) in &t.phases {
            let _ = writeln!(text, "phase.{phase}_ms\t{:.3}", secs * 1e3);
        }
        let _ = writeln!(text, "total_ms\t{:.3}", t.total_secs * 1e3);
        let _ = writeln!(text, "chunks\t{}", t.chunks);
        let _ = writeln!(text, "chunks_read\t{}", t.chunks_read);
        let _ = writeln!(text, "worker_busy_ms\t{:.3}", t.worker_busy_secs * 1e3);
        let _ = writeln!(text, "reader_ms\t{:.3}", t.reader_secs * 1e3);
        let _ = writeln!(text, "reduce_ms\t{:.3}", t.reduce_secs * 1e3);
        let _ = writeln!(text, "reader_wait_ms\t{:.3}", t.reader_wait_secs * 1e3);
        let _ = writeln!(text, "reduce_wait_ms\t{:.3}", t.reducer_wait_secs * 1e3);
        let _ = writeln!(text, "inflate_ms\t{:.3}", t.inflate_secs * 1e3);
        let _ = writeln!(text, "job_ms\t{:.3}", t.job_secs * 1e3);
        let _ = writeln!(text, "deflate_ms\t{:.3}", t.deflate_secs * 1e3);
        let _ = writeln!(text, "chain_wait_ms\t{:.3}", t.chain_wait_secs * 1e3);
        if let Some(peak) = t.in_flight_peak_bytes {
            let _ = writeln!(text, "in_flight_peak_bytes\t{peak}");
        }
        text
    }
}

/// Policy trigger counters (§11.5): the decode-level transformations the
/// scanner applied, mirroring `reed_core::genotype::PolicyCounters`.
/// Rejections are errors, so they never appear in a written report.
#[derive(Clone, Copy, Debug, Default, Serialize)]
pub struct PolicyCounters {
    /// Haploid genotypes converted to half calls (`allow_haploid`).
    pub haploid_as_half_call: u64,
    /// Haploid genotypes duplicated into homozygous calls.
    pub haploid_as_homozygous: u64,
    /// Haploid genotypes zeroed to missing.
    pub haploid_as_missing: u64,
    /// Half calls (including haploid-as-half-call) zeroed to missing.
    pub half_call_as_missing: u64,
    /// ALT alleles with index > 1 collapsed onto the single ALT bit.
    pub multiallelic_collapsed: u64,
}

impl PolicyCounters {
    /// The counters the decoder produced for one record.
    #[must_use]
    pub fn from_decode(c: &reed_core::genotype::PolicyCounters) -> Self {
        Self {
            haploid_as_half_call: c.n_haploid_as_half_call,
            haploid_as_homozygous: c.n_haploid_as_homozygous,
            haploid_as_missing: c.n_haploid_as_missing,
            half_call_as_missing: c.n_half_call_as_missing,
            multiallelic_collapsed: c.n_multiallelic_collapsed,
        }
    }

    /// Adds the triggers of one decoded record.
    pub fn add_decode(&mut self, c: &reed_core::genotype::PolicyCounters) {
        self.haploid_as_half_call += c.n_haploid_as_half_call;
        self.haploid_as_homozygous += c.n_haploid_as_homozygous;
        self.haploid_as_missing += c.n_haploid_as_missing;
        self.half_call_as_missing += c.n_half_call_as_missing;
        self.multiallelic_collapsed += c.n_multiallelic_collapsed;
    }
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
    /// Records dropped by the LD-pruning rule (`after_mask.ld_prune`,
    /// docs/rounds/v5/05-algorithms.md §3.4.1).
    pub n_ld_pruned: u64,
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
        self.n_ld_pruned = 0;
        self.n_input = 0;
        self.n_output = 0;
        self.policy = PolicyCounters::default();
    }

    /// Adds another partial (component-wise, order-independent for the
    /// integer components).
    pub fn merge(&mut self, other: &Self) {
        // The receiving side may start with unsized (default) vectors;
        // extend so every counter of `other` lands somewhere.
        self.head_hits
            .resize(self.head_hits.len().max(other.head_hits.len()), 0);
        self.tail_hits
            .resize(self.tail_hits.len().max(other.tail_hits.len()), 0);
        self.mask_hits
            .resize(self.mask_hits.len().max(other.mask_hits.len()), 0);
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
        self.n_ld_pruned += other.n_ld_pruned;
        self.n_input += other.n_input;
        self.n_output += other.n_output;
        self.policy.haploid_as_half_call += other.policy.haploid_as_half_call;
        self.policy.haploid_as_homozygous += other.policy.haploid_as_homozygous;
        self.policy.haploid_as_missing += other.policy.haploid_as_missing;
        self.policy.half_call_as_missing += other.policy.half_call_as_missing;
        self.policy.multiallelic_collapsed += other.policy.multiallelic_collapsed;
    }
}

#[cfg(test)]
mod tests {
    use super::FunnelPartial;

    /// doc_v7 D7-28: `n_ld_pruned` must survive `merge` and be cleared by
    /// `reset` (it was omitted from both).
    #[test]
    fn merge_and_reset_carry_the_ld_prune_count() {
        let mut a = FunnelPartial {
            n_ld_pruned: 2,
            ..FunnelPartial::default()
        };
        let b = FunnelPartial {
            n_ld_pruned: 3,
            ..FunnelPartial::default()
        };
        a.merge(&b);
        assert_eq!(a.n_ld_pruned, 5);
        a.reset(1, 1, 1);
        assert_eq!(a.n_ld_pruned, 0, "reset must clear the LD-prune count");
    }
}
