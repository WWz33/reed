//! `FilterSpecV2::to_args`: the inverse of the CLI translation
//! (`doc_v2/03-architecture.md` §8, `spec/render.rs`).
//!
//! Round trip: `to_args()` renders a spec as the flag vector that
//! `reed-cli`'s `spec_from_flags` parses back into the same spec. Two spec
//! shapes have no flag spelling and are rendered as their equivalent
//! expression form instead (documented per tag): a typed `min_qual` becomes
//! `-i QUAL>=P`, and typed INFO predicates (`site.info`) become no flag at
//! all — the CLI carries INFO rules as expressions in `site.include`.

use super::{FilterSpecV2, HweMode, OnMissing, SampleList};

/// Appends a flag plus its value.
fn flag(args: &mut Vec<String>, name: &str, value: impl std::fmt::Display) {
    args.push(name.to_string());
    args.push(value.to_string());
}

/// Renders `value` the way the CLI parses floats back (`{}` round-trips).
fn f(value: f64) -> String {
    let mut s = format!("{value}");
    if !s.contains(['.', 'e', 'E']) {
        s.push_str(".0");
    }
    s
}

impl FilterSpecV2 {
    /// The flag vector equivalent to this spec (02-cli.md §3).
    ///
    /// The order follows the evaluation order of 02-cli.md §6 so the output
    /// reads like the documented command line.
    #[must_use]
    pub fn to_args(&self) -> Vec<String> {
        let mut args: Vec<String> = Vec::new();

        // ---- samples (static subset, then the dynamic --mind threshold) ----
        if let Some(select) = &self.sample_select {
            match &select.keep {
                Some(SampleList::Names(names)) => flag(&mut args, "-s", names.join(",")),
                Some(SampleList::File(path)) => flag(&mut args, "-S", path),
                None => {}
            }
            if let Some(remove) = &select.remove {
                match remove {
                    // The whole list carries one `^` prefix (bcftools syntax).
                    SampleList::Names(names) => {
                        flag(&mut args, "-s", format!("^{}", names.join(",")));
                    }
                    SampleList::File(path) => flag(&mut args, "--remove", path),
                }
            }
        }
        if let Some(threshold) = self.mind {
            flag(&mut args, "--mind", threshold);
        }

        // ---- regions (one `-t` list: `^` prefixes the excludes) ----
        if !self.targets.include.is_empty() || !self.targets.exclude.is_empty() {
            let mut list: Vec<String> = self.targets.include.clone();
            list.extend(self.targets.exclude.iter().map(|r| format!("^{r}")));
            flag(&mut args, "-t", list.join(","));
        }

        // ---- variant class / alleles / FILTER / QUAL ----
        if !self.site.variant_class.is_empty() {
            flag(&mut args, "-v", self.site.variant_class.join(","));
        }
        if !self.site.exclude_variant_class.is_empty() {
            flag(&mut args, "-V", self.site.exclude_variant_class.join(","));
        }
        if self.site.min_alt_alleles > 0 {
            // bcftools counts REF+ALT, so ALT count n renders as n+1.
            flag(&mut args, "-m", self.site.min_alt_alleles + 1);
        }
        if self.site.max_alt_alleles > 0 {
            flag(&mut args, "-M", self.site.max_alt_alleles + 1);
        }
        if !self.site.require_filter.is_empty() {
            flag(&mut args, "-f", self.site.require_filter.join(","));
        }
        if let Some(qual) = self.site.min_qual {
            flag(&mut args, "-i", format!("QUAL>={qual}"));
        }

        // ---- expressions (source order: includes then excludes) ----
        for text in &self.site.include {
            flag(&mut args, "-i", text);
        }
        for text in &self.site.exclude {
            flag(&mut args, "-e", text);
        }
        for text in &self.site_after_mask.include {
            flag(&mut args, "-i", text);
        }
        for text in &self.site_after_mask.exclude {
            flag(&mut args, "-e", text);
        }

        // ---- set-gts mode (the FMT expressions above become masks) ----
        if let Some(value) = &self.genotype.set_gts {
            flag(&mut args, "--set-gts", value);
        }

        // ---- genotype masking (the three CLI spellings) ----
        for rule in &self.genotype.mask_if {
            let when = rule.when.trim();
            let name = match (&*rule.field, when.starts_with('<'), when.starts_with('>')) {
                ("DP", true, _) if rule.on_missing == OnMissing::Keep => "--min-dp",
                ("DP", _, true) if rule.on_missing == OnMissing::Keep => "--max-dp",
                ("GQ", true, _) if rule.on_missing == OnMissing::Keep => "--min-gq",
                // Not one of the CLI spellings (the CLI never builds those).
                _ => continue,
            };
            let value = when.trim_start_matches(['<', '>']).trim();
            flag(&mut args, name, value);
        }

        // ---- after-mask thresholds (02-cli.md §6 step 6 order) ----
        if let Some(g) = self.site_after_mask.max_missing {
            flag(&mut args, "--geno", f(g));
        }
        if let Some(m) = self.site_after_mask.min_maf {
            flag(&mut args, "--maf", f(m));
        }
        if let Some(m) = self.site_after_mask.max_maf {
            flag(&mut args, "--max-maf", f(m));
        }
        if self.site_after_mask.min_mac > 0 {
            flag(&mut args, "--mac", self.site_after_mask.min_mac);
        }
        if self.site_after_mask.max_mac > 0 {
            flag(&mut args, "--max-mac", self.site_after_mask.max_mac);
        }
        if let Some(p) = self.site_after_mask.min_hwe_p {
            args.push("--hwe".to_string());
            args.push(p.to_string());
            if self.site_after_mask.hwe_mode == HweMode::MidP {
                args.push("midp".to_string());
            }
        }

        // ---- INFO ----
        if !self.recompute.drop_stale.is_empty() {
            let list: Vec<String> = self
                .recompute
                .drop_stale
                .iter()
                .map(|t| format!("INFO/{t}"))
                .collect();
            flag(&mut args, "-x", list.join(","));
        }
        if !self.recompute.tags.is_empty() {
            flag(&mut args, "--fill-tags", self.recompute.tags.join(","));
        }
        if self.recompute.no_update {
            args.push("-I".to_string());
        }

        args
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args_of(spec: &FilterSpecV2) -> Vec<String> {
        spec.to_args()
    }

    #[test]
    #[test]
    fn renders_a_combined_sample_selection() {
        let spec = FilterSpecV2 {
            sample_select: Some(super::super::SampleSelect {
                keep: Some(SampleList::Names(vec!["a".into(), "b".into()])),
                remove: Some(SampleList::Names(vec!["c".into()])),
            }),
            mind: Some(0.25),
            ..Default::default()
        };
        let args = args_of(&spec);
        assert_eq!(&args[..6], &["-s", "a,b", "-s", "^c", "--mind", "0.25"]);
    }

    #[test]
    fn renders_the_documented_flags() {
        let spec = FilterSpecV2 {
            site: super::super::SiteRulesV2 {
                variant_class: vec!["snv".into()],
                exclude_variant_class: vec!["indel".into()],
                max_alt_alleles: 1,
                require_filter: vec!["PASS".into()],
                ..Default::default()
            },
            site_after_mask: super::super::TailRules {
                max_missing: Some(0.1),
                min_maf: Some(0.05),
                max_maf: Some(0.4),
                min_mac: 3,
                max_mac: 30,
                min_hwe_p: Some(1e-6),
                ..Default::default()
            },
            recompute: super::super::RecomputeSpec {
                tags: vec!["AC".into(), "AN".into()],
                drop_stale: vec!["DP".into()],
                no_update: true,
            },
            ..Default::default()
        };
        assert_eq!(
            args_of(&spec),
            vec![
                "-v", "snv", "-V", "indel", "-M", "2", "-f", "PASS", "--geno", "0.1", "--maf",
                "0.05", "--max-maf", "0.4", "--mac", "3", "--max-mac", "30", "--hwe", "0.000001",
                "-x", "INFO/DP", "--fill-tags", "AC,AN", "-I",
            ]
        );
    }

    #[test]
    fn floats_always_render_as_parsable() {
        let spec = FilterSpecV2 {
            site_after_mask: super::super::TailRules {
                max_missing: Some(0.5),
                ..Default::default()
            },
            ..Default::default()
        };
        assert!(args_of(&spec).contains(&"0.5".to_string()));
    }
}
