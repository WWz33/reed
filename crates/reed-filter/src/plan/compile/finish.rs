//! The final stage: strategy choice, recompute/stale tags, decode policies,
//! provenance and the assembled [`FilterPlan`](crate::plan::FilterPlan).

use reed_core::tag_table::TagTable;

use super::{HeadOutcome, RuleSet};
use crate::FilterError;
use crate::plan::{ExecStrategy, FilterPlan, RecomputeFlags, TailRule, provenance};
use crate::spec::FilterSpec;

/// Assembles the plan from the three sections plus the finishing checks
/// (doc_v25 08 §3 ④).
///
/// # Errors
///
/// [`FilterError::Spec`] on strategy/recompute/provenance validation.
pub(super) fn compile(
    spec: &FilterSpec,
    tags: &TagTable,
    command_line: &str,
    head: HeadOutcome,
    rules: RuleSet,
    tail: Vec<TailRule>,
) -> Result<FilterPlan, FilterError> {
    let HeadOutcome {
        head,
        cohort_changed,
        format_tags: _,
    } = head;
    let RuleSet {
        mask,
        int_fields,
        float_fields,
        fmt_expr_fields,
    } = rules;
    // ---- strategy (§11.4; SketchThenEmit selection is M5b) ----
    // The sketch path evaluates tail rules in the mid stage, where the
    // record bytes are gone; a tail rule that reads them therefore forces
    // the single-pass path (02-cli.md §4: expressions run on the record).
    // The predicate mirrors `FilterPlan::tail_reads_record` (the plan does
    // not exist yet here), so the driver's later override asks the same
    // question (doc_v12 D12-1).
    let tail_reads_record = tail.iter().any(|r| {
        matches!(
            r,
            TailRule::Expr { .. } | TailRule::MinAf(..) | TailRule::MaxAf(..)
        )
    });
    let strategy = if tail_reads_record {
        ExecStrategy::SinglePass
    } else if spec.mind.is_some() {
        ExecStrategy::TwoPassReparse
    } else {
        ExecStrategy::SinglePass
    };

    // ---- recompute + stale ----
    // `--fill-tags` entries are validated here (not only in the CLI): an
    // unsupported tag is a spec error, like bcftools `+fill-tags`
    // ("the tag X is not supported", plugins/fill-tags.c:530-535). Silently
    // accepting one turned `--fill-tags AF` into a no-op that the default
    // `-x` then used to delete the input's AF (doc_v12 D12-11).
    let canonical_tags: Vec<&'static str> = spec
        .recompute
        .tags
        .iter()
        .map(|tag| {
            crate::spec::recompute::canonical_fill_tag(tag)
                .map_err(|e| FilterError::Spec(e.to_string()))
        })
        .collect::<Result<_, _>>()?;
    let no_update = spec.recompute.no_update;
    // Registry order is fill order (AC, AN, NS, MAF, HWE), independent of
    // the order the user listed them in.
    let recompute = RecomputeFlags {
        tags: reed_core::count::RECOMPUTE_TAGS
            .iter()
            .filter(|tag| canonical_tags.contains(&tag.name))
            .collect(),
    };
    // bcftools `view -I`: keep the input's AC/AN. The default `-x` list names
    // them because the *update* replaces them; with the update suppressed
    // there is nothing to replace, so they must stay. An explicit
    // `-x INFO/AC` is still the user's request and applies
    // (`drop_stale_explicit`, doc_v30 §4.4).
    let mut stale_tags = spec.recompute.stale_bytes();
    if no_update && !spec.recompute.drop_stale_explicit {
        stale_tags.retain(|t| t.as_slice() != b"AC" && t.as_slice() != b"AN");
    }

    // ---- ascertainment (§11.7) + provenance (§11.8; the type lives in
    // reed-core, doc_v25 V25-08c) ----
    let provenance = provenance::Provenance::new(command_line)
        .with_spec(spec, spec.to_args().join(" "))
        .map_err(|e| FilterError::Spec(format!("spec serialization failed: {e}")))?;

    let policies = reed_core::genotype::GenotypePolicies {
        ploidy: match spec.policy.ploidy {
            crate::spec::PloidyPolicy::RequireDiploid => {
                reed_core::genotype::PloidyPolicy::RequireDiploid
            }
            crate::spec::PloidyPolicy::AllowHaploid => {
                reed_core::genotype::PloidyPolicy::HaploidAsHalfCall
            }
        },
        multiallelic: reed_core::genotype::MultiallelicPolicy::CollapseAlt,
        ..reed_core::genotype::GenotypePolicies::default()
    };

    Ok(FilterPlan {
        head,
        mask,
        tail,
        strategy,
        int_fields,
        float_fields,
        fmt_expr_fields,
        policies,
        stale_tags,
        auto_stale_tags: Vec::new(),
        recompute,
        no_update,
        info_trusted: !cohort_changed,
        tags: std::sync::Arc::new(tags.clone()),
        provenance,
    })
}
