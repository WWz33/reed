//! `reed filter` flag translation (02-cli.md §3, §6): the typed
//! [`reed_filter::spec::FilterSpecV2`] plus the driver call. The inverse
//! mapping (`FilterSpec::to_args`) lives in `reed-filter`.

use std::path::PathBuf;


use crate::cli::FilterCli;
use crate::units::input_path;

/// Translates the CLI flags into a typed [`FilterSpecV2`] (02-cli.md §3, §6).
///
/// # Errors
///
/// Unknown tokens (variant types, output type) and malformed expressions.
pub fn spec_from_flags(cli: &FilterCli) -> Result<reed_filter::spec::FilterSpecV2, String> {
    use reed_filter::spec::{
        FilterSpecV2, GenotypeMaskRule, OnMissing, SampleList, SampleSelect,
    };
    let mut spec = FilterSpecV2::default();

    // ---- samples: the static subset (04 §C0b). `-s` may appear twice —
    // once with names to keep, once with a `^`-prefixed removal list — but
    // two sources of the same polarity are a conflict, not a silent pick.
    let parse_list = |list: &str| -> Result<(bool, Vec<String>), String> {
        let exclude = list.trim_start().starts_with('^');
        let names: Vec<String> = list
            .trim_start_matches('^')
            .split(',')
            .map(str::trim)
            .filter(|n| !n.is_empty())
            .map(str::to_string)
            .collect();
        if names.is_empty() {
            return Err(format!("empty sample list {list:?}"));
        }
        Ok((exclude, names))
    };
    let mut inline_keep: Option<Vec<String>> = None;
    let mut inline_remove: Option<Vec<String>> = None;
    for list in &cli.samples {
        let (exclude, names) = parse_list(list)?;
        let slot = if exclude { &mut inline_remove } else { &mut inline_keep };
        if slot.is_some() {
            return Err(format!(
                "repeated -s list of the same polarity ({list:?}); combine the names into one list"
            ));
        }
        *slot = Some(names);
    }
    let file_keep = cli.samples_file.as_ref().or(cli.keep.as_ref());
    if inline_keep.is_some() && file_keep.is_some() {
        return Err("-s and -S/--keep are mutually exclusive: give one sample source".into());
    }
    if inline_remove.is_some() && cli.remove.is_some() {
        return Err("-s ^… and --remove are mutually exclusive: give one removal source".into());
    }
    let keep = inline_keep
        .map(SampleList::Names)
        .or_else(|| file_keep.map(|p| SampleList::File(p.clone())));
    let remove = inline_remove
        .map(SampleList::Names)
        .or_else(|| cli.remove.clone().map(SampleList::File));
    if keep.is_some() || remove.is_some() {
        spec.sample_select = Some(SampleSelect { keep, remove });
    }
    spec.mind = cli.mind;

    // ---- regions (02-cli.md §2): the spec predicate, validated here.
    // `-r` (indexed) and `-t` (streamed) share it; `-T`/`^` add excludes ----
    let mut add_regions = |list: &str, include: bool| -> Result<(), String> {
        for raw in list.split(',').map(str::trim).filter(|t| !t.is_empty()) {
            let (token, include) = match raw.strip_prefix('^') {
                Some(rest) => (rest, false),
                None => (raw, include),
            };
            reed_filter::spec::Region::parse(token)
                .map_err(|e| format!("invalid region {token:?}: {e}"))?;
            if include {
                spec.targets.include.push(token.to_string());
            } else {
                spec.targets.exclude.push(token.to_string());
            }
        }
        Ok(())
    };
    if let Some(list) = &cli.regions {
        add_regions(list, true)?;
    }
    if let Some(path) = &cli.regions_file {
        let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
        add_regions(&text.lines().collect::<Vec<_>>().join(","), true)?;
    }
    if let Some(list) = &cli.targets {
        add_regions(list, true)?;
    }
    if let Some(path) = &cli.targets_file {
        let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
        add_regions(&text.lines().collect::<Vec<_>>().join(","), true)?;
    }

    // ---- site: variant class / alleles / FILTER ----
    let variant_class = |token: &str| -> Result<&'static str, String> {
        Ok(match token.trim() {
            "snps" | "snp" | "snv" | "snvs" => "snv",
            "mnps" | "mnp" | "mnv" => "mnv",
            "indels" | "indel" => "indel",
            "other" | "symbolic" => "symbolic",
            "ref" => "ref",
            other => return Err(format!("unknown variant type {other:?}")),
        })
    };
    if let Some(list) = &cli.types {
        for token in list.split(',') {
            spec.site.variant_class.push(variant_class(token)?.to_string());
        }
    }
    if let Some(list) = &cli.exclude_types {
        for token in list.split(',') {
            spec.site
                .exclude_variant_class
                .push(variant_class(token)?.to_string());
        }
    }
    if let Some(modifier) = &cli.snps_only {
        if !matches!(modifier.as_str(), "just-acgt") {
            return Err(format!("unknown --snps-only modifier {modifier:?}"));
        }
        if !spec.site.variant_class.iter().any(|c| c == "snv") {
            spec.site.variant_class.push("snv".to_string());
        }
    }
    // bcftools counts REF+ALT, the spec counts ALTs.
    if let Some(m) = cli.min_alleles {
        spec.site.min_alt_alleles = m.saturating_sub(1);
    }
    if let Some(m) = cli.max_alleles {
        spec.site.max_alt_alleles = m.saturating_sub(1);
    }
    if let Some(list) = &cli.apply_filters {
        spec.site.require_filter = list
            .split(',')
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(str::to_string)
            .collect();
    }

    // ---- expressions (02-cli.md §4): route by genotype dependence ----
    for (text, include) in cli
        .include
        .iter()
        .map(|t| (t, true))
        .chain(cli.exclude.iter().map(|t| (t, false)))
    {
        let expr = reed_filter::expr::parse(text)
            .map_err(|e| format!("{} {text:?}: {e}", if include { "-i" } else { "-e" }))?;
        let genotype_derived = expr.references_genotypes();
        match (include, genotype_derived) {
            (true, false) => spec.site.include.push(text.clone()),
            (true, true) => spec.site_after_mask.include.push(text.clone()),
            (false, false) => spec.site.exclude.push(text.clone()),
            (false, true) => spec.site_after_mask.exclude.push(text.clone()),
        }
    }

    // ---- genotype masking (plink2 semantics: missing never masks) ----
    let mut mask = |field: &str, when: String| {
        spec.genotype.mask_if.push(GenotypeMaskRule {
            field: field.to_string(),
            when,
            on_missing: OnMissing::Keep,
        });
    };
    if let Some(v) = cli.min_dp {
        mask("DP", format!("< {v}"));
    }
    if let Some(v) = cli.max_dp {
        mask("DP", format!("> {v}"));
    }
    if let Some(v) = cli.min_gq {
        mask("GQ", format!("< {v}"));
    }
    if let Some(value) = &cli.set_gts {
        spec.genotype.set_gts = Some(value.clone());
    }

    // ---- after-mask site rules ----
    if let Some(g) = cli.geno {
        spec.site_after_mask.max_missing = Some(g);
    }
    if let Some(m) = cli.maf {
        spec.site_after_mask.min_maf = Some(m);
    }
    if let Some(m) = cli.max_maf {
        spec.site_after_mask.max_maf = Some(m);
    }
    if let Some(m) = cli.mac {
        spec.site_after_mask.min_mac = m;
    }
    if let Some(m) = cli.max_mac {
        spec.site_after_mask.max_mac = m;
    }
    spec.recompute.no_update = cli.no_update;
    if !cli.hwe.is_empty() {
        let p: f64 = cli.hwe[0]
            .parse()
            .map_err(|_| format!("invalid --hwe p-value {:?}", cli.hwe[0]))?;
        spec.site_after_mask.min_hwe_p = Some(p);
        if cli.hwe.get(1).is_some_and(|m| m == "midp") {
            spec.site_after_mask.hwe_mode = reed_filter::spec::HweMode::MidP;
        }
    }

    // ---- INFO ----
    if let Some(list) = &cli.remove_info {
        spec.recompute.drop_stale = list
            .split(',')
            .map(str::trim)
            .map(|t| t.strip_prefix("INFO/").unwrap_or(t).to_string())
            .filter(|t| !t.is_empty())
            .collect();
    }
    if let Some(list) = &cli.fill_tags {
        spec.recompute.tags = list
            .split(',')
            .map(str::trim)
            .filter(|t| !t.is_empty())
            .map(str::to_string)
            .collect();
    }
    Ok(spec)
}

pub fn run_filter_cli(cli: FilterCli) -> Result<(), Box<dyn std::error::Error + Send + Sync>> {
    let FilterCli {
        input,
        out,
        write_index,
        compress_level,
        report,
        strategy,
        threads,
        output_type,
        force_samples,
        ..
    } = &cli;
    let spec_v2 = spec_from_flags(&cli)?;
    let strategy = strategy
        .as_deref()
        .map(|s| match s {
            "single_pass" => Ok(reed_filter::plan::ExecStrategy::SinglePass),
            "two_pass_reparse" => Ok(reed_filter::plan::ExecStrategy::TwoPassReparse),
            "sketch_then_emit" => Ok(reed_filter::plan::ExecStrategy::SketchThenEmit),
            other => Err(format!("unknown strategy {other:?}")),
        })
        .transpose()?;
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");
    let threads = threads.unwrap_or_else(|| {
        std::thread::available_parallelism().expect("invariant: parallelism").get()
    });
    // -O v|z|b|u[level]: explicit output type, `-o` names the file (the
    // extension is only inferred when -O is absent).
    let mut output_bcf = None;
    let (output, compress_level, output_compressed) = match output_type.as_deref() {
        None => (out.clone().map(PathBuf::from), *compress_level, None),
        Some(spec) => {
            let (kind, level) = spec.split_at(1);
            let level = if level.is_empty() {
                *compress_level
            } else {
                level.parse::<u8>().map_err(|_| format!("invalid -O level in {spec:?}"))?
            };
            match kind {
                "v" => (out.clone().map(PathBuf::from), level, Some(false)),
                "z" => (out.clone().map(PathBuf::from), level, Some(true)),
                "b" => {
                    output_bcf = Some(level);
                    (out.clone().map(PathBuf::from), level, None)
                }
                // `-Ou`: BGZF level-0 blocks (htslib's uncompressed BCF).
                "u" => {
                    output_bcf = Some(0);
                    (out.clone().map(PathBuf::from), 0, None)
                }
                other => return Err(format!("unknown output type {other:?}").into()),
            }
        }
    };
    let mut regions = Vec::new();
    let indexed_list = match (&cli.regions, &cli.regions_file) {
        (Some(list), _) => Some(list.clone()),
        (None, Some(path)) => Some(
            std::fs::read_to_string(path)
                .map_err(|e| format!("{path}: {e}"))?
                .lines()
                .collect::<Vec<_>>()
                .join(","),
        ),
        (None, None) => None,
    };
    if let Some(list) = &indexed_list {
        for token in list.split(',').map(str::trim).filter(|t| !t.is_empty()) {
            let parsed = reed_filter::spec::Region::parse(token)
                .map_err(|e| format!("invalid region {token:?}: {e}"))?;
            regions.push(reed_io::chunk::region::Region {
                contig: parsed.contig.clone(),
                start1: parsed.start,
                end1: parsed.end,
            });
        }
    }
    let opts = reed::filter::FilterOptions {
        threads,
        output,
        write_index: reed::filter::WriteIndex::parse(write_index.as_deref())?,
        compress_level,
        strategy,
        command_line,
        quiet: false,
        report_json: report.as_ref().map(PathBuf::from),
        force_samples: *force_samples,
        regions,
        tmp_dir: cli.tmp_dir.as_ref().map(PathBuf::from),
        output_compressed,
        output_bcf,
    };
    let report = reed::filter::run_filter(input_path(input).as_deref(), &spec_v2, &opts)?;
    eprintln!(
        "reed filter: {} -> {} records in {:.1}s",
        report.n_input_records, report.n_output_records, report.elapsed_secs
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use reed_filter::spec::FilterSpecV2;

    /// Parses a flag vector through the real clap definitions.
    fn cli_from(args: &[String]) -> FilterCli {
        // `Cli` requires a subcommand; prepend the program name and `filter`.
        let mut argv = vec!["reed".to_string(), "filter".to_string(), "in.vcf.gz".to_string()];
        argv.extend_from_slice(args);
        match crate::cli::Cli::try_parse_from(argv).expect("invariant: flags parse").command {
            crate::cli::Command::Filter(cli) => cli,
            _ => unreachable!("invariant: the filter subcommand was given"),
        }
    }

    /// R5 gate (02-cli.md §3): `to_args` and the flag translation are
    /// inverses — every spec the CLI can build renders to flags that parse
    /// back into the same spec.
    #[test]
    fn spec_to_args_round_trips_through_the_flags() {
        let flags: Vec<&str> = vec![
            "-s", "s1,s2",
            "-t", "Chr01:1-1000,^Chr02",
            "-v", "snv", "-V", "indel",
            "-m", "2", "-M", "3",
            "-f", "PASS",
            "-i", "QD>=4",
            "-e", "FS>50",
            "-i", "F_MISSING<0.05",
            "-e", "N_MISSING>10",
            "--min-dp", "5", "--max-dp", "100", "--min-gq", "20",
            "--geno", "0.1", "--maf", "0.05", "--max-maf", "0.4",
            "--mac", "3", "--max-mac", "30",
            "--hwe", "0.000001", "midp",
            "-x", "INFO/DP,INFO/MLEAC",
            "--fill-tags", "AC,AN,NS",
            "-I",
        ];
        let first = spec_from_flags(&cli_from(&flags.iter().map(ToString::to_string).collect::<Vec<_>>()))
            .expect("invariant: the reference flags translate");
        let rendered = first.to_args();
        let second = spec_from_flags(&cli_from(&rendered))
            .expect("invariant: the rendered flags translate");
        assert_eq!(first, second, "round trip differs; rendered: {rendered:?}");
        assert_eq!(rendered, flags, "the rendering is the canonical spelling");
    }

    /// `--mind` is the dynamic sample rule the spec renders as a threshold
    /// (and the CLI parses back into `mind`).
    #[test]
    fn mind_round_trips() {
        let flags = vec!["--mind".to_string(), "0.25".to_string()];
        let first = spec_from_flags(&cli_from(&flags)).expect("invariant: --mind translates");
        assert_eq!(first.mind, Some(0.25));
        assert_eq!(first.sample_select, None);
        let rendered = first.to_args();
        let second = spec_from_flags(&cli_from(&rendered)).expect("invariant: renders back");
        assert_eq!(first, second);
        assert_eq!(&rendered[..flags.len()], &flags[..]);
        // The rest of the rendering is the default `-x` stale list.
        assert!(rendered[flags.len()] == "-x", "{rendered:?}");
    }

    /// C0b: a static subset (file) and an inline removal list plus `--mind`
    /// round-trip; the rendering is the canonical spelling.
    #[test]
    fn static_subset_and_mind_round_trip() {
        let flags: Vec<String> = ["-S", "keep.txt", "-s", "^drop1,drop2", "--mind", "0.2"]
            .iter()
            .map(ToString::to_string)
            .collect();
        let first = spec_from_flags(&cli_from(&flags)).expect("invariant: translate");
        let select = first.sample_select.as_ref().expect("invariant: subset");
        assert_eq!(select.keep, Some(reed_filter::spec::SampleList::File("keep.txt".into())));
        assert_eq!(
            select.remove,
            Some(reed_filter::spec::SampleList::Names(vec!["drop1".into(), "drop2".into()]))
        );
        assert_eq!(first.mind, Some(0.2));
        let rendered = first.to_args();
        assert_eq!(&rendered[..flags.len()], &flags[..]);
        let second = spec_from_flags(&cli_from(&rendered)).expect("invariant: renders back");
        assert_eq!(first, second);
    }

    /// C0b: same-polarity sample sources are a conflict, not a silent pick.
    #[test]
    fn conflicting_sample_sources_error() {
        for flags in [
            vec!["-s", "A", "-S", "f"],
            vec!["-s", "A", "--keep", "f"],
            vec!["-s", "^A", "--remove", "f"],
            vec!["-s", "A", "-s", "B"],
        ] {
            let args: Vec<String> = flags.iter().map(ToString::to_string).collect();
            let err = spec_from_flags(&cli_from(&args)).expect_err("must conflict");
            assert!(err.contains("mutually exclusive") || err.contains("repeated"), "{err}");
        }
    }

    /// The rendered flags are also what `##reed_filterCommand` reports, so
    /// the provenance line and the CLI agree on one spelling.
    #[test]
    fn no_update_renders_only_when_set() {
        let base = FilterSpecV2::default();
        assert!(!base.to_args().contains(&"-I".to_string()));
        let spec = FilterSpecV2 {
            recompute: reed_filter::spec::RecomputeSpec {
                no_update: true,
                ..reed_filter::spec::RecomputeSpec::default()
            },
            ..FilterSpecV2::default()
        };
        assert!(spec.to_args().contains(&"-I".to_string()));
    }
}
