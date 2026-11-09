//! Mask-rule compilation: `[[genotype.mask_if]]` field rules and the
//! FORMAT columns they scan.

use reed_core::tag_table::{TagTable, TagType};

use super::RuleSet;
use crate::FilterError;
use crate::plan::{FmtExprField, IntField, MaskRule, validate_numeric};
use crate::spec::{FilterSpec, OnMissing};

/// Compiles the mask section (doc_v25 08 §3 ④: moved from `compile()`).
pub(super) fn compile(spec: &FilterSpec, tags: &TagTable) -> Result<RuleSet, FilterError> {
    // ---- mask rules ----
    let fmt_expr_fields: Vec<FmtExprField> = Vec::new();
    let mut mask = Vec::new();
    let mut int_fields: Vec<IntField> = Vec::new();
    let mut float_fields: Vec<Vec<u8>> = Vec::new();
    for rule in &spec.genotype.mask_if {
        // `depth = true` is the explicit depth policy (doc_v22 D22-1): the
        // record's DP when it carries one, otherwise the sample's AD sum.
        // Only `DP` can be a depth field.
        if rule.depth && rule.field != "DP" {
            return Err(FilterError::Spec(format!(
                "mask rule {}: depth = true is only valid for FORMAT/DP",
                rule.field
            )));
        }
        let fallback = rule.depth.then(|| b"AD".to_vec());
        // The fallback reads AD as integers whenever a record's FORMAT has
        // no DP key. Validate that even when DP is declared: otherwise a
        // Float AD makes the mask a silent no-op for exactly those records
        // (doc_v22 V22-R7).
        if rule.depth
            && let Some(ad) = tags.format("AD")
            && ad.ty != TagType::Integer
        {
            return Err(FilterError::Spec(
                "FORMAT field AD must be Integer for the depth fallback (DP is read when the \
                 record carries it)"
                    .to_string(),
            ));
        }
        let def = match tags.format(&rule.field) {
            Some(def) => def,
            // A depth mask on a callset that only emits AD (`GT:PL:AD`,
            // common for BSA-seq/GATK-lite pipelines). The type check above
            // already rejected a non-Integer AD.
            None if rule.depth => tags.format("AD").ok_or_else(|| {
                FilterError::Spec(
                    "FORMAT field DP is not declared in the header, and neither is an \
                     Integer FORMAT/AD to fall back to"
                        .to_string(),
                )
            })?,
            None => {
                return Err(FilterError::Spec(format!(
                    "FORMAT field {} is not declared in the header",
                    rule.field
                )));
            }
        };
        validate_numeric(&rule.field, def.ty)?;
        let is_float = def.ty == TagType::Float;
        if is_float && rule.depth {
            return Err(FilterError::Spec(
                "FORMAT field DP is Float; the depth mask needs an Integer DP (or Integer AD \
                 to fall back to)"
                    .to_string(),
            ));
        }
        let cmp = crate::spec::predicate::Cmp::parse(&rule.when)?;
        // The mask default is keep (§11.1): an omitted on_missing
        // deserializes to the enum default (Pass) and normalizes to Keep.
        let on_missing = match rule.on_missing {
            OnMissing::Keep | OnMissing::Pass => OnMissing::Keep,
            OnMissing::Mask | OnMissing::Fail => OnMissing::Mask,
        };
        let int_column_of = |fields: &mut Vec<IntField>| match fields
            .iter()
            .position(|f| f.name == rule.field.as_bytes() && f.fallback == fallback)
        {
            Some(i) => i,
            None => {
                fields.push(IntField {
                    name: rule.field.as_bytes().to_vec(),
                    fallback: fallback.clone(),
                });
                fields.len() - 1
            }
        };
        let float_column_of = |fields: &mut Vec<Vec<u8>>| match fields
            .iter()
            .position(|f| f.as_slice() == rule.field.as_bytes())
        {
            Some(i) => i,
            None => {
                fields.push(rule.field.as_bytes().to_vec());
                fields.len() - 1
            }
        };
        let (int_column, float_column) = if is_float {
            (0, Some(float_column_of(&mut float_fields)))
        } else {
            (int_column_of(&mut int_fields), None)
        };
        mask.push(MaskRule {
            expr: None,
            include: true,
            field: rule.field.as_bytes().to_vec(),
            cmp,
            on_missing,
            int_column,
            float_column,
        });
    }

    Ok(RuleSet {
        mask,
        int_fields,
        float_fields,
        fmt_expr_fields,
    })
}
