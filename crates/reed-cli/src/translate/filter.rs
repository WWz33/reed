//! `reed filter` flag translation (02-cli.md §3, §6): the typed
//! [`reed::spec::FilterSpec`] plus the driver call. The inverse
//! mapping (`FilterSpec::to_args`) lives in `reed-filter`.

use std::path::PathBuf;

use crate::cli::FilterCli;
use crate::units::input_path;

/// Translates the shared site/threshold flags (`SiteFilterArgs`) into the
/// typed spec — the same flags serve `reed filter` and `reed stats`
/// (04-cleanup §2.1).
///
/// # Errors
///
/// Unknown variant types or malformed expressions.
/// Parses `--indep-pairwise W[kb] [S] R2` (docs/rounds/v5/05-algorithms.md
/// §3.4.1): `W` variants (or `W kb`), the step defaults to 1 and must be 1
/// when the window is in kb.
fn parse_indep_pairwise(args: &[String]) -> Result<Option<reed::spec::LdPruneSpec>, String> {
    if args.is_empty() {
        return Ok(None);
    }
    // `W[kb] [S] R2`: the r² threshold is always the LAST token, so parse by
    // position — an integer r² (`50 5 1`, `50 5 0`) must be accepted (doc_v7
    // D7-18; the old "does the token contain a dot" heuristic mistook the
    // integer for the step).
    let window_token = args.first().expect("invariant: args is non-empty").clone();
    let (step, r2_token) = match args.len() {
        2 => (1_usize, &args[1]),
        3 => (
            args[1]
                .parse()
                .map_err(|_| format!("invalid --indep-pairwise step {:?}", args[1]))?,
            &args[2],
        ),
        n => {
            return Err(format!(
                "--indep-pairwise wants `W[kb] [S] R2`, got {n} value(s)"
            ));
        }
    };
    let r2: f64 = r2_token
        .parse()
        .map_err(|_| "--indep-pairwise r² must be a number".to_string())?;
    // `parse_size` scales the full token itself (`50kb` → 50,000); stripping
    // the suffix first would silently read the window as bare bp. A suffixed
    // token is a bp window (plink2's `kb` form, extended to reed's other
    // length suffixes); a bare number is a variant count.
    let kb = !window_token.bytes().all(|b| b.is_ascii_digit());
    let window = crate::units::parse_size(&window_token)?;
    if kb && step != 1 {
        return Err("--indep-pairwise with a kb window requires step 1".into());
    }
    Ok(Some(reed::spec::LdPruneSpec {
        window,
        kb,
        step,
        r2,
    }))
}

pub(crate) fn site_spec(
    args: &crate::cli::filter_flags::SiteFilterArgs,
) -> Result<reed::spec::FilterSpec, String> {
    use reed::spec::{FilterSpec, GenotypeMaskRule, OnMissing};
    // `site_spec` serves both `reed filter` and `reed stats`; the removed
    // drop-direction flags never reach here (their value parsers reject them
    // with the conversion hint at parse time).
    let cli = args;
    let mut spec = FilterSpec::default();
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
            spec.site
                .variant_class
                .push(variant_class(token)?.to_string());
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
        spec.site.max_alt_alleles = Some(m.saturating_sub(1));
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
        let expr = reed::expr::parse(text)
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
    if let Some(text) = &cli.min_af {
        let (f, kind) = crate::cli::filter_flags::parse_af_with_type(text)?;
        spec.site_after_mask.min_af = Some(f);
        spec.site_after_mask.min_af_kind = Some(kind);
    }
    if let Some(text) = &cli.max_af {
        let (f, kind) = crate::cli::filter_flags::parse_af_with_type(text)?;
        spec.site_after_mask.max_af = Some(f);
        spec.site_after_mask.max_af_kind = Some(kind);
    }
    if let Some(rate) = cli.geno {
        spec.site_after_mask.max_missing = Some(rate);
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
    if !cli.hwe.is_empty() {
        let p: f64 = cli.hwe[0]
            .parse()
            .map_err(|_| format!("invalid --hwe p-value {:?}", cli.hwe[0]))?;
        spec.site_after_mask.min_hwe_p = Some(p);
        if cli.hwe.get(1).is_some_and(|m| m == "midp") {
            spec.site_after_mask.hwe_mode = reed::spec::HweMode::MidP;
        }
    }
    for text in &cli.genotype {
        spec.site_after_mask
            .gt_class
            .push(reed::spec::GtClassRule::parse(text).map_err(|e| e.to_string())?);
    }

    Ok(spec)
}

/// Translates the CLI flags into a typed [`reed::spec::FilterSpec`]
/// (02-cli.md §3, §6): the shared site/threshold flags (04-cleanup §2.1)
/// plus the sample, region and INFO-update parts that stay filter-only.
///
/// # Errors
///
/// Unknown tokens (variant types, output type) and malformed expressions.
pub(crate) fn spec_from_flags(cli: &FilterCli) -> Result<reed::spec::FilterSpec, String> {
    use reed::spec::{SampleList, SampleSelect};
    let mut spec = site_spec(&cli.site_args)?;

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
    for list in &cli.sample_args.samples {
        let (exclude, names) = parse_list(list)?;
        let slot = if exclude {
            &mut inline_remove
        } else {
            &mut inline_keep
        };
        if slot.is_some() {
            return Err(format!(
                "repeated -s list of the same polarity ({list:?}); combine the names into one list"
            ));
        }
        *slot = Some(names);
    }
    let file_keep = cli
        .sample_args
        .samples_file
        .as_ref()
        .or(cli.sample_args.keep.as_ref());
    if inline_keep.is_some() && file_keep.is_some() {
        return Err("-s and -S/--keep are mutually exclusive: give one sample source".into());
    }
    if inline_remove.is_some() && cli.sample_args.remove.is_some() {
        return Err("-s ^… and --remove are mutually exclusive: give one removal source".into());
    }
    let keep = inline_keep
        .map(SampleList::Names)
        .or_else(|| file_keep.map(|p| SampleList::File(p.clone())));
    let remove = inline_remove
        .map(SampleList::Names)
        .or_else(|| cli.sample_args.remove.clone().map(SampleList::File));
    if keep.is_some() || remove.is_some() {
        spec.sample_select = Some(SampleSelect { keep, remove });
    }
    spec.mind = cli.mind;
    spec.ld_prune = parse_indep_pairwise(&cli.site_args.indep_pairwise)?;

    // ---- regions (02-cli.md §2): the spec predicate, validated here.
    // `-r` (indexed) and `-t` (streamed) share it; `-T`/`^` add excludes ----
    let mut add_regions = |list: &str, include: bool| -> Result<(), String> {
        for raw in list.split(',').map(str::trim).filter(|t| !t.is_empty()) {
            let (token, include) = match raw.strip_prefix('^') {
                Some(rest) => (rest, false),
                None => (raw, include),
            };
            reed::spec::Region::parse(token)
                .map_err(|e| format!("invalid region {token:?}: {e}"))?;
            if include {
                spec.targets.include.push(token.to_string());
            } else {
                spec.targets.exclude.push(token.to_string());
            }
        }
        Ok(())
    };
    if let Some(list) = &cli.regions_args.regions {
        add_regions(list, true)?;
    }
    if let Some(path) = &cli.regions_args.regions_file {
        let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
        add_regions(&text.lines().collect::<Vec<_>>().join(","), true)?;
    }
    if let Some(list) = &cli.regions_args.targets {
        add_regions(list, true)?;
    }
    if let Some(path) = &cli.regions_args.targets_file {
        let text = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
        add_regions(&text.lines().collect::<Vec<_>>().join(","), true)?;
    }

    // ---- INFO (filter-only: `-x`, `--fill-tags`, `-I`) ----
    spec.recompute.no_update = cli.no_update;
    if let Some(list) = &cli.remove_info {
        spec.recompute.drop_stale_explicit = true;
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

pub(crate) fn run_filter_cli(cli: FilterCli) -> Result<(), reed::Error> {
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
        memory,
        verbosity,
        no_version,
        ..
    } = &cli;
    verbosity.init_tracing();
    let spec = spec_from_flags(&cli).map_err(reed::Error::usage)?;
    let strategy = strategy
        .as_deref()
        .map(|s| match s {
            "single_pass" => Ok(reed::plan::ExecStrategy::SinglePass),
            "two_pass_reparse" => Ok(reed::plan::ExecStrategy::TwoPassReparse),
            "sketch_then_emit" => Ok(reed::plan::ExecStrategy::SketchThenEmit),
            other => Err(reed::Error::usage(format!("unknown strategy {other:?}"))),
        })
        .transpose()?;
    let command_line = std::env::args().collect::<Vec<_>>().join(" ");
    let threads = threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .expect("invariant: parallelism")
            .get()
    });
    // -O v|z|b|u[level]: explicit output type, `-o` names the file (the
    // extension is only inferred when -O is absent).
    let mut output_bcf = None;
    let (output, compress_level, output_compressed) = match output_type.as_deref() {
        None => (out.as_deref().and_then(crate::units::output_path), *compress_level, None),
        Some(spec) => {
            let (kind, level) = spec.split_at(1);
            let level = if level.is_empty() {
                *compress_level
            } else {
                level
                    .parse::<u8>()
                    .map_err(|_| reed::Error::usage(format!("invalid -O level in {spec:?}")))?
            };
            match kind {
                "v" => (out.as_deref().and_then(crate::units::output_path), level, Some(false)),
                "z" => (out.as_deref().and_then(crate::units::output_path), level, Some(true)),
                "b" => {
                    output_bcf = Some(level);
                    (out.as_deref().and_then(crate::units::output_path), level, None)
                }
                // `-Ou`: BGZF level-0 blocks (htslib's uncompressed BCF).
                "u" => {
                    output_bcf = Some(0);
                    (out.as_deref().and_then(crate::units::output_path), 0, None)
                }
                other => return Err(reed::Error::usage(format!("unknown output type {other:?}"))),
            }
        }
    };
    let mut regions = Vec::new();
    let indexed_list = match (&cli.regions_args.regions, &cli.regions_args.regions_file) {
        (Some(list), _) => Some(list.clone()),
        (None, Some(path)) => Some(
            std::fs::read_to_string(path)
                .map_err(|e| reed::Error::usage(format!("{path}: {e}")))?
                .lines()
                .collect::<Vec<_>>()
                .join(","),
        ),
        (None, None) => None,
    };
    if let Some(list) = &indexed_list {
        for token in list.split(',').map(str::trim).filter(|t| !t.is_empty()) {
            let parsed = reed::spec::Region::parse(token)
                .map_err(|e| reed::Error::usage(format!("invalid region {token:?}: {e}")))?;
            regions.push(reed::chunk::region::Region {
                contig: parsed.contig.clone(),
                start1: parsed.start,
                end1: parsed.end,
            });
        }
    }
    let opts = reed::FilterOptions {
        threads,
        output,
        write_index: reed::WriteIndex::parse(write_index.as_deref())?,
        compress_level,
        strategy,
        command_line,
        timing: cli.timing.as_ref().map(PathBuf::from),
        report_json: report.as_ref().map(PathBuf::from),
        force_samples: *force_samples,
        regions,
        tmp_dir: cli.tmp_dir.as_ref().map(PathBuf::from),
        output_compressed,
        output_bcf,
        mem_bytes: memory.map(|m| m.saturating_mul(1 << 20)),
        no_version: *no_version,
        no_command: cli.no_command,
    };
    tracing::debug!(threads, strategy = ?strategy, "reed filter starting");
    let report = reed::run_filter(input_path(input).as_deref(), &spec, &opts)?;
    if verbosity.shows_summary() {
        eprint_report(&report);
    }
    if verbosity.shows_timing() {
        eprint!("{}", report.timing_text());
    }
    tracing::debug!(
        n_input = report.n_input_records,
        n_output = report.n_output_records,
        samples = report.n_samples_out,
        "reed filter finished"
    );
    Ok(())
}

/// The human summary (`--verbosity 1`): one line plus the funnel table
/// (doc_v4 03 §3.3). Moved out of the driver so the library never prints.
fn eprint_report(report: &reed::funnel::FilterReport) {
    eprintln!(
        "reed filter: {} -> {} records ({} of {} samples), strategy {}, {:.1}s",
        report.n_input_records,
        report.n_output_records,
        report.n_samples_out,
        report.n_samples_in,
        report.strategy,
        report.elapsed_secs
    );
    if !report.funnel.is_empty() {
        eprintln!("funnel:");
        for stage in &report.funnel {
            eprintln!("  {:<28} {}", stage.rule, stage.n_passed);
        }
    }
    if let Some(rule) = &report.sample_rule {
        eprintln!(
            "  {:<28} {} ({} dropped)",
            rule.rule, report.n_samples_out, rule.n_dropped
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use reed::spec::FilterSpec;

    /// Parses a flag vector through the real clap definitions.
    fn cli_from(args: &[String]) -> FilterCli {
        // `Cli` requires a subcommand; prepend the program name and `filter`.
        let mut argv = vec![
            "reed".to_string(),
            "filter".to_string(),
            "in.vcf.gz".to_string(),
        ];
        argv.extend_from_slice(args);
        match crate::cli::Cli::try_parse_from(argv)
            .expect("invariant: flags parse")
            .command
        {
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
            "-s",
            "s1,s2",
            "-t",
            "Chr01:1-1000,^Chr02",
            "-v",
            "snv",
            "-V",
            "indel",
            "-m",
            "2",
            "-M",
            "3",
            "-f",
            "PASS",
            "-i",
            "QD>=4",
            "-e",
            "FS>50",
            "-i",
            "F_MISSING<0.05",
            "-e",
            "N_MISSING>10",
            "--min-dp",
            "5",
            "--max-dp",
            "100",
            "--min-gq",
            "20",
            "--geno",
            "0.1",
            "--maf",
            "0.05",
            "--max-maf",
            "0.4",
            "--mac",
            "3",
            "--max-mac",
            "30",
            "--hwe",
            "0.000001",
            "midp",
            "-x",
            "INFO/DP,INFO/MLEAC",
            "--fill-tags",
            "AC,AN,NS",
            "-I",
        ];
        let first = spec_from_flags(&cli_from(
            &flags.iter().map(ToString::to_string).collect::<Vec<_>>(),
        ))
        .expect("invariant: the reference flags translate");
        let rendered = first.to_args();
        let second =
            spec_from_flags(&cli_from(&rendered)).expect("invariant: the rendered flags translate");
        assert_eq!(first, second, "round trip differs; rendered: {rendered:?}");
        assert_eq!(rendered, flags, "the rendering is the canonical spelling");
    }

    /// `-g/--genotype` (bcftools): repeated flags AND (bcftools keeps the
    /// last one); the spec renders them back in order, in the canonical
    /// lowercase spelling.
    #[test]
    fn genotype_rules_round_trip() {
        let flags: Vec<String> = ["-g", "het", "-g", "^miss", "-g", "Hom"]
            .iter()
            .map(ToString::to_string)
            .collect();
        let first = spec_from_flags(&cli_from(&flags)).expect("invariant: translates");
        use reed::spec::GtClassRule as R;
        assert_eq!(
            first.site_after_mask.gt_class,
            vec![R::Het, R::NoMiss, R::Hom]
        );
        let rendered = first.to_args();
        assert_eq!(&rendered[..4], &["-g", "het", "-g", "^miss",]);
        let second = spec_from_flags(&cli_from(&rendered)).expect("invariant: renders back");
        assert_eq!(first, second);
        // Bad tokens are usage errors at parse time.
        let err = spec_from_flags(&cli_from(&["-g".into(), "homo".into()]))
            .expect_err("unknown class must fail");
        assert!(err.contains("hom, het, miss"), "{err}");
    }

    /// `--mind` is the dynamic sample rule the spec renders as a threshold
    /// (and the CLI parses back into `mind`).
    #[test]
    fn mind_round_trips() {
        let flags = vec!["--mind".to_string(), "0.25".to_string()];
        let first = spec_from_flags(&cli_from(&flags)).expect("invariant: translates");
        assert_eq!(
            first.mind,
            Some(reed::spec::MissingRate::parse("0.25").expect("valid"))
        );
        assert_eq!(first.sample_select, None);
        let rendered = first.to_args();
        let second = spec_from_flags(&cli_from(&rendered)).expect("invariant: renders back");
        assert_eq!(first, second);
        assert_eq!(&rendered[..flags.len()], &flags[..]);
        // The default stale list is not rendered (doc_v12 D12-12): the CLI
        // applies it anyway, and rendering it would turn the round trip into
        // an explicit `-x` (`drop_stale_explicit`).
        assert_eq!(rendered, flags, "{rendered:?}");
    }

    /// doc_v4 02 回退版 §5.2: the withdrawn keep-direction flags error with
    /// the converted spelling (exact decimal), and no silent alias exists.
    #[test]
    fn withdrawn_flags_error_with_the_conversion() {
        use clap::Parser;
        let cases = [
            (vec!["--max-missing", "0.9"], "--geno 0.1"),
            (vec!["--max-missing", "0.8"], "--geno 0.2"),
            (vec!["--max-missing-sample", "0.9"], "--mind 0.1"),
            (
                vec!["--max-missing-count", "5"],
                "the missing-threshold family is --geno/--mind only",
            ),
        ];
        for (flags, expected) in cases {
            let mut argv = vec![
                "reed".to_string(),
                "filter".to_string(),
                "in.vcf.gz".to_string(),
            ];
            argv.extend(flags.iter().map(ToString::to_string));
            let err = crate::cli::Cli::try_parse_from(argv)
                .err()
                .expect("withdrawn flag must fail");
            assert_eq!(
                err.kind(),
                clap::error::ErrorKind::ValueValidation,
                "{flags:?}"
            );
            let text = err.to_string();
            assert!(text.contains("was removed"), "{flags:?}: {text}");
            assert!(text.contains(expected), "{flags:?}: {text}");
        }
    }

    /// C0b: a static subset (file) and an inline removal list plus the
    /// dynamic sample rule round-trip; the rendering is the canonical
    /// spelling.
    #[test]
    fn static_subset_and_mind_round_trip() {
        let flags: Vec<String> = ["-S", "keep.txt", "-s", "^drop1,drop2", "--mind", "0.2"]
            .iter()
            .map(ToString::to_string)
            .collect();
        let first = spec_from_flags(&cli_from(&flags)).expect("invariant: translate");
        let select = first.sample_select.as_ref().expect("invariant: subset");
        assert_eq!(
            select.keep,
            Some(reed::spec::SampleList::File("keep.txt".into()))
        );
        assert_eq!(
            select.remove,
            Some(reed::spec::SampleList::Names(vec![
                "drop1".into(),
                "drop2".into()
            ]))
        );
        assert_eq!(
            first.mind,
            Some(reed::spec::MissingRate::parse("0.2").expect("valid"))
        );
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
            assert!(
                err.contains("mutually exclusive") || err.contains("repeated"),
                "{err}"
            );
        }
    }

    /// The rendered flags are also what `##reed_filterCommand` reports, so
    /// the provenance line and the CLI agree on one spelling.
    #[test]
    fn no_update_renders_only_when_set() {
        let base = FilterSpec::default();
        assert!(!base.to_args().contains(&"-I".to_string()));
        let spec = FilterSpec {
            recompute: reed::spec::RecomputeSpec {
                no_update: true,
                ..reed::spec::RecomputeSpec::default()
            },
            ..FilterSpec::default()
        };
        assert!(spec.to_args().contains(&"-I".to_string()));
    }
}
