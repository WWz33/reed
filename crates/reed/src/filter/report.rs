//! The run report (`doc_v2/03-architecture.md` §8, §11.5): funnel stages and
//! mask attribution. Printing lives in the CLI (doc_v4 03 §3.3).

use reed_filter::funnel::{FunnelPartial, FunnelStage, MaskedBy};

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
