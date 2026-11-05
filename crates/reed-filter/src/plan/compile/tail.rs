//! Tail-rule compilation: `-q/-Q`, missing/MAF/MAC/HWE, genotype classes
//! and the `[site_after_mask]` expressions (including the `--set-gts .`
//! mask rules they produce).

use reed_core::tag_table::{TagNumber, TagType, TagTable};

use super::RuleSet;
use crate::FilterError;
use crate::plan::{FmtExprField, IntField, MaskRule, TailRule, af_threshold, validate_numeric, warn_stale_info};
use crate::spec::FilterSpec;

/// Compiles the tail section (doc_v25 08 §3 ④: moved from `compile()`).
pub(super) fn compile(
    spec: &FilterSpec,
    tags: &TagTable,
    format_tags: &[String],
    rules: &mut RuleSet,
) -> Result<Vec<TailRule>, FilterError> {
    // ---- tail rules (fixed order §11.2) ----
    let mut tail = Vec::new();
    // bcftools applies -q/-Q (site-level AF bounds) before the plink2-style
    // missing/MAF blocks.
    if let Some(f) = spec.site_after_mask.min_af {
        if !(0.0..=1.0).contains(&f) {
            return Err(FilterError::Spec(format!("min_af out of range [0,1]: {f}")));
        }
        tail.push(TailRule::MinAf(
            af_threshold(f),
            spec.site_after_mask
                .min_af_kind
                .unwrap_or(crate::spec::AfKind::Nref),
        ));
    }
    if let Some(f) = spec.site_after_mask.max_af {
        if !(0.0..=1.0).contains(&f) {
            return Err(FilterError::Spec(format!("max_af out of range [0,1]: {f}")));
        }
        tail.push(TailRule::MaxAf(
            af_threshold(f),
            spec.site_after_mask
                .max_af_kind
                .unwrap_or(crate::spec::AfKind::Nref),
        ));
    }
    if let Some(rate) = spec.site_after_mask.max_missing {
        // The range is validated where the rate is parsed
        // (spec/missing_rate.rs), so a compiled spec can only carry
        // `q ≤ 10^9`.
        tail.push(TailRule::MaxMissing(rate));
    }
    if let Some(m) = spec.site_after_mask.min_maf {
        if !(0.0..=0.5).contains(&m) {
            return Err(FilterError::Spec(format!(
                "min_maf out of range [0,0.5]: {m}"
            )));
        }
        tail.push(TailRule::MinMaf(m));
    }
    if let Some(m) = spec.site_after_mask.max_maf {
        if !(0.0..=0.5).contains(&m) {
            return Err(FilterError::Spec(format!(
                "max_maf out of range [0,0.5]: {m}"
            )));
        }
        tail.push(TailRule::MaxMaf(m));
    }
    // doc_v21 S13: `--maf 0.4 --max-maf 0.2` (and the MAC pair) is an empty
    // intersection — bcftools/plink2 reject it; reed silently emitted an
    // empty file.
    if let (Some(min), Some(max)) = (
        spec.site_after_mask.min_maf,
        spec.site_after_mask.max_maf,
    ) && min > max
    {
        return Err(FilterError::Spec(format!(
            "--maf {min} exceeds --max-maf {max} (empty intersection)"
        )));
    }
    if spec.site_after_mask.min_mac > 0 {
        tail.push(TailRule::MinMac(spec.site_after_mask.min_mac));
    }
    if spec.site_after_mask.max_mac > 0 {
        tail.push(TailRule::MaxMac(spec.site_after_mask.max_mac));
    }
    if spec.site_after_mask.max_mac > 0 && spec.site_after_mask.min_mac > spec.site_after_mask.max_mac
    {
        return Err(FilterError::Spec(format!(
            "--mac {} exceeds --max-mac {} (empty intersection)",
            spec.site_after_mask.min_mac, spec.site_after_mask.max_mac
        )));
    }
    if let Some(p) = spec.site_after_mask.min_hwe_p {
        if spec.site_after_mask.hwe_mode == crate::spec::HweMode::MidP {
            // doc_v21 S12: no internal milestone code; the exact test is the
            // only implemented mode and the message says so.
            return Err(FilterError::Spec(
                "--hwe midp is not implemented yet: only the exact test is available; \
                 drop the `midp` token"
                    .into(),
            ));
        }
        if !(0.0..=1.0).contains(&p) {
            return Err(FilterError::Spec(format!(
                "min_hwe_p out of range [0,1]: {p}"
            )));
        }
        tail.push(TailRule::MinHweP(p));
    }
    // `-g/--genotype` (repeated flags AND, unlike bcftools' last-wins):
    // exact duplicates collapse; requiring a class and its absence is
    // unsatisfiable — reject at compile time rather than silently drop
    // every record.
    let mut gt_rules: Vec<crate::spec::GtClassRule> =
        Vec::with_capacity(spec.site_after_mask.gt_class.len());
    for &rule in &spec.site_after_mask.gt_class {
        if gt_rules.contains(&rule) {
            continue; // an exact duplicate: AND is idempotent
        }
        if let Some(prior) = gt_rules.iter().find(|r| r.class() == rule.class()) {
            return Err(FilterError::Spec(format!(
                "-g/--genotype rules {:?} and {:?} are contradictory \
                 (a site cannot have and lack a {} call)",
                prior.id(),
                rule.id(),
                rule.class(),
            )));
        }
        gt_rules.push(rule);
    }
    tail.extend(gt_rules.into_iter().map(TailRule::GtClass));

    // `--set-gts .` turns the `FMT/` expressions into per-sample mask rules
    // (bcftools `filter -S .`); genotype-derived site fields stay tail rules.
    let set_gts = match spec.genotype.set_gts.as_deref() {
        None => false,
        Some(".") => true,
        Some(other) => {
            return Err(FilterError::Spec(format!(
                "--set-gts {other:?}: only `.` (the missing genotype) is supported"
            )));
        }
    };
    for (index, include, text) in spec
        .site_after_mask
        .include
        .iter()
        .enumerate()
        .map(|(i, t)| (i, true, t))
        .chain(
            spec.site_after_mask
                .exclude
                .iter()
                .enumerate()
                .map(|(i, t)| (i, false, t)),
        )
    {
        let flag = if include { "-i" } else { "-e" };
        let expr = crate::expr::parse_and_check_fmt(text, tags, format_tags)
            .map_err(|e| FilterError::Spec(format!("{flag} {text:?} (after mask): {e}")))?;
        warn_stale_info(&expr, text, tags);
        if !expr.references_genotypes() {
            return Err(FilterError::Spec(format!(
                "{flag} {text:?} does not read genotype-derived values and belongs in [site]"
            )));
        }
        let mut fmt_tags = Vec::new();
        expr.fmt_tags(&mut fmt_tags);
        if fmt_tags.is_empty() {
            tail.push(TailRule::Expr {
                id: if include {
                    format!("after_mask.include#{}", index + 1)
                } else {
                    format!("after_mask.exclude#{}", index + 1)
                },
                include,
                expr,
            });
            continue;
        }
        if !set_gts {
            return Err(FilterError::Spec(format!(
                "{flag} {text:?} reads FORMAT values: they are masks, not site rules — add --set-gts ."
            )));
        }
        if expr.reads_gt_fields() {
            return Err(FilterError::Spec(format!(
                "{flag} {text:?} mixes FORMAT values with genotype-derived fields; \
                 --set-gts . expressions may read FORMAT and site fields only"
            )));
        }
        for tag in &fmt_tags {
            let def = tags.format(tag).ok_or_else(|| {
                FilterError::Spec(format!("FORMAT tag {tag} is not declared in the header"))
            })?;
            if def.number != TagNumber::Count(1) {
                return Err(FilterError::Spec(format!(
                    "FORMAT tag {tag} is multi-valued (Number={:?}): subscripts are not supported yet",
                    def.number
                )));
            }
            validate_numeric(tag, def.ty)?;
            let is_float = def.ty == TagType::Float;
            let int_column_of = |fields: &mut Vec<IntField>| match fields
                .iter()
                .position(|f| f.name == tag.as_bytes() && f.fallback.is_none())
            {
                Some(i) => i,
                None => {
                    fields.push(IntField {
                        name: tag.as_bytes().to_vec(),
                        // An explicit `FMT/...` expression is literal: no
                        // depth fallback (doc_v22 D22-1).
                        fallback: None,
                    });
                    fields.len() - 1
                }
            };
            let float_column_of = |fields: &mut Vec<Vec<u8>>| match fields
                .iter()
                .position(|f| f.as_slice() == tag.as_bytes())
            {
                Some(i) => i,
                None => {
                    fields.push(tag.as_bytes().to_vec());
                    fields.len() - 1
                }
            };
            let (int_column, float_column) = if is_float {
                (0, Some(float_column_of(&mut rules.float_fields)))
            } else {
                (int_column_of(&mut rules.int_fields), None)
            };
            rules.fmt_expr_fields.push(FmtExprField {
                tag: tag.as_bytes().to_vec(),
                int_column,
                float_column,
            });
        }
        rules.mask.push(MaskRule {
            expr: Some(expr),
            include,
            field: Vec::new(),
            cmp: crate::spec::Cmp::parse(">= 0").expect("invariant: a constant comparison"),
            on_missing: crate::spec::OnMissing::Keep,
            int_column: 0,
            float_column: None,
        });
    }

    Ok(tail)
}
