//! Head-rule compilation: fixed columns, `-q/-Q` INFO bounds and the
//! `[site]` expressions (`-i/-e`, `[[site.info]]`).

use reed_core::tag_table::{TagNumber, TagTable};

use super::HeadOutcome;
use crate::FilterError;
use crate::plan::{HeadRule, InfoRule, af_info_sides, af_threshold, compile_regions, validate_numeric, warn_stale_af, warn_stale_info};
use crate::spec::{FilterSpec, OnMissing, Quantifier};

/// Compiles the head section (doc_v25 08 §3 ④: moved from `compile()`).
pub(super) fn compile(
    spec: &FilterSpec,
    tags: &TagTable,
    biallelic_gate: bool,
) -> Result<HeadOutcome, FilterError> {
    // ---- head rules (§11.2 order: targets → class/alt → FILTER → QUAL →
    // INFO in declaration order) ----
    let mut head = Vec::new();
    // The indexed family (`-r/-R`) and the streaming family (`-t/-T`) are two
    // separate predicates, so a site must satisfy both — bcftools'
    // intersection (doc_v23 V23-61).
    let indexed =
        crate::spec::region::parse_regions_strict(&spec.targets.indexed, tags.contigs())?;
    if !indexed.is_empty() {
        head.push(HeadRule::Regions {
            include: indexed,
            exclude: Vec::new(),
        });
    }
    let regions = compile_regions(&spec.targets, tags.contigs())?;
    if let Some(r) = regions {
        head.push(HeadRule::Regions {
            include: r.0,
            exclude: r.1,
        });
    }
    if !spec.site.variant_class.is_empty() {
        head.push(HeadRule::VariantClass(spec.site.variant_class.clone()));
    }
    if spec.site.just_acgt {
        head.push(HeadRule::JustAcgt);
    }
    if !spec.site.exclude_variant_class.is_empty() {
        head.push(HeadRule::ExcludeVariantClass(
            spec.site.exclude_variant_class.clone(),
        ));
    }
    if spec.site.min_alt_alleles > 0 {
        head.push(HeadRule::MinAltAlleles(spec.site.min_alt_alleles));
    }
    if let Some(max) = spec.site.max_alt_alleles {
        head.push(HeadRule::MaxAltAlleles(max));
    } else if biallelic_gate && !spec.multiallelics_allowed() {
        // The default gate (§11.6, plink2-style): without `multiallelic =
        // "allow"` or an explicit `max_alt_alleles > 1`, multiallelic
        // records do not reach the decoder — the collapsed per-record
        // counts (and the recomputed AC/MAF/HWE) are biallelic-valid only.
        head.push(HeadRule::BiallelicGate);
    }
    if !spec.site.require_filter.is_empty() {
        head.push(HeadRule::RequireFilter(
            spec.site
                .require_filter
                .iter()
                .map(String::as_bytes)
                .map(<[u8]>::to_vec)
                .collect(),
        ));
    }
    if let Some(q) = spec.site.min_qual {
        if q < 0.0 {
            return Err(FilterError::Spec(format!(
                "min_qual out of range [0,inf): {q}"
            )));
        }
        head.push(HeadRule::MinQual(q));
    }
    // `-q/-Q` (bcftools): the INFO AC/AN values decide at the site stage when
    // the header declares both AND the threshold is the `nref` form; the
    // other `:type`s are genotype-derived and live in the tail only.
    //
    // The INFO preference holds only while AC/AN still describe the cohort
    // being decided on. Once the sample set changes (a static subset, a
    // genotype mask, `--mind`, or `--set-gts .`), bcftools recomputes from
    // the subset GT before applying `-q/-Q` (vcfview.c:341-350) and gates
    // its own INFO preference on `!sample_mask_set` (filter.c:1471); reed
    // mirrors that by not pushing the INFO short-circuit at all, so the
    // decision runs through `af_of`, whose `info_trusted` flag is false in
    // every execution path once the cohort changed (doc_v12 D12-6).
    let info_af = tags.info("AC").is_some() && tags.info("AN").is_some();
    // An empty `SampleSelect` (JSON `"sample_select": {}`) selects nothing;
    // treating it as "cohort changed" silently diverted `-q/-Q` from the
    // INFO short-circuit (doc_v23 V23-25).
    let cohort_changed = spec
        .sample_select
        .as_ref()
        .is_some_and(|select| !select.is_empty())
        || spec.mind.is_some()
        || !spec.genotype.mask_if.is_empty()
        || spec.genotype.set_gts.is_some();
    let (min_af_info, max_af_info) = af_info_sides(&spec.site_after_mask, info_af, cohort_changed);
    // `-q/-Q` bounds are validated here but the RULES are pushed after the
    // `-i/-e` expressions (bcftools runs `filter_test` at vcfview.c:322,
    // before its `-q/-Q` block at :418): a record an expression excludes must
    // never reach the AF test, so a malformed `AC>AN` record cannot abort a
    // run whose surviving set is non-empty (`bcftools view -e 'QUAL<5' -q`
    // exits 0 on that input).
    let mut af_rules: Vec<HeadRule> = Vec::new();
    for bound in spec.site_after_mask.min_af.iter().copied() {
        if !(0.0..=1.0).contains(&bound) {
            return Err(FilterError::Spec(format!(
                "min_af out of range [0,1]: {bound}"
            )));
        }
        if min_af_info {
            af_rules.push(HeadRule::MinAfInfo(af_threshold(bound)));
        }
    }
    for bound in spec.site_after_mask.max_af.iter().copied() {
        if !(0.0..=1.0).contains(&bound) {
            return Err(FilterError::Spec(format!(
                "max_af out of range [0,1]: {bound}"
            )));
        }
        if max_af_info {
            af_rules.push(HeadRule::MaxAfInfo(af_threshold(bound)));
        }
    }
    // Warn only after both bounds validate, so an out-of-range threshold
    // errors without the stale-INFO noise first.
    if min_af_info || max_af_info {
        warn_stale_af();
    }
    // ---- expressions (02-cli.md §4): non-genotype-dependent ones run here,
    // genotype-dependent ones after the mask stage ----
    let format_tags: Vec<String> = tags.format_names();
    for (index, text) in spec.site.include.iter().enumerate() {
        let expr = crate::expr::parse_and_check(text, tags)
            .map_err(|e| FilterError::Spec(format!("-i {text:?}: {e}")))?;
        warn_stale_info(&expr, text, tags);
        if expr.references_genotypes() {
            return Err(FilterError::Spec(format!(
                "-i {text:?} reads genotype-derived values and belongs in [site_after_mask]"
            )));
        }
        head.push(HeadRule::Expr {
            id: format!("site.include#{}", index + 1),
            include: true,
            expr,
        });
    }
    for (index, text) in spec.site.exclude.iter().enumerate() {
        let expr = crate::expr::parse_and_check(text, tags)
            .map_err(|e| FilterError::Spec(format!("-e {text:?}: {e}")))?;
        warn_stale_info(&expr, text, tags);
        if expr.references_genotypes() {
            return Err(FilterError::Spec(format!(
                "-e {text:?} reads genotype-derived values and belongs in [site_after_mask]"
            )));
        }
        head.push(HeadRule::Expr {
            id: format!("site.exclude#{}", index + 1),
            include: false,
            expr,
        });
    }
    // Now the AF bounds, after every expression rule (bcftools' order).
    head.extend(af_rules);
    for rule in &spec.site.info {
        let def = tags.info(&rule.tag).ok_or_else(|| {
            FilterError::Spec(format!(
                "INFO tag {} is not declared in the header",
                rule.tag
            ))
        })?;
        validate_numeric(&rule.tag, def.ty)?;
        let number_a = def.number == TagNumber::A;
        if def.number != TagNumber::Count(1) && !number_a {
            return Err(FilterError::Spec(format!(
                "INFO tag {} has an unsupported Number (want 1 or A)",
                rule.tag
            )));
        }
        let quantifier = if number_a {
            rule.quantifier.ok_or_else(|| {
                FilterError::Spec(format!(
                    "INFO tag {} is Number=A: declare quantifier = \"any\" | \"all\"",
                    rule.tag
                ))
            })?
        } else {
            // A Number=1 tag with a declared quantifier would silently
            // change nothing; keep the declaration honest.
            // A declared quantifier on a Number=1 tag is harmless (the
            // value set has one element); accepted silently.
            Quantifier::Any
        };
        let cmp = crate::spec::predicate::Cmp::parse(&rule.keep_if)?;
        // Missing-tag behaviour: bcftools semantics default to pass.
        let on_missing = rule.on_missing;
        if on_missing == OnMissing::Keep || on_missing == OnMissing::Mask {
            return Err(FilterError::Spec(format!(
                "INFO rule {} uses a genotype on_missing value (pass|fail only)",
                rule.tag
            )));
        }
        head.push(HeadRule::Info(InfoRule {
            tag: rule.tag.as_bytes().to_vec(),
            cmp,
            on_missing,
            quantifier,
            number_a,
        }));
    }

    Ok(HeadOutcome {
        head,
        cohort_changed,
        format_tags,
    })
}
