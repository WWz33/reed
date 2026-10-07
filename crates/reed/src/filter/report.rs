//! The run report (`doc_v2/03-architecture.md` §8, §11.5): funnel stages,
//! mask attribution and the human stderr line.

use reed_filter::funnel::{FilterReport, FunnelPartial, FunnelStage, MaskedBy};

/// Funnel stages in rule order: head rules then tail rules.
#[must_use]
pub(super) fn funnel_stages(
    head_ids: &[String],
    tail_ids: &[String],
    funnel: &FunnelPartial,
) -> Vec<FunnelStage> {
    let mut stages: Vec<FunnelStage> = Vec::new();
    for (id, n) in head_ids.iter().zip(&funnel.head_hits) {
        stages.push(FunnelStage { rule: id.clone(), n_passed: *n });
    }
    for (id, n) in tail_ids.iter().zip(&funnel.tail_hits) {
        stages.push(FunnelStage { rule: id.clone(), n_passed: *n });
    }
    stages
}

/// Mask attribution, one entry per mask rule in declaration order.
#[must_use]
pub(super) fn masked_by_rule(mask_ids: &[String], funnel: &FunnelPartial) -> Vec<MaskedBy> {
    mask_ids
        .iter()
        .zip(&funnel.mask_hits)
        .map(|(id, n)| MaskedBy { rule: id.clone(), n_masked: *n })
        .collect()
}

/// The one-line human report on stderr.
pub(super) fn eprint_report(report: &FilterReport) {
    eprintln!(
        "reed filter: {} -> {} records ({} of {} samples), strategy {}, {:.1}s",
        report.n_input_records,
        report.n_output_records,
        report.n_samples_out,
        report.n_samples_in,
        report.strategy,
        report.elapsed_secs
    );
}
