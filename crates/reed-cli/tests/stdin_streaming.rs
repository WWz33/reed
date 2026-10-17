//! O5 gates (03-performance.md): stdin streams for single-pass commands and
//! leaves no temporary file behind; multi-pass runs spill under `--tmp-dir`
//! and delete the spill on exit (04 §C6).

use std::fmt::Write as _;
use std::io::Write;
use std::process::{Command, Stdio};

const HEADER: &str = "\
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=100000>
##FILTER=<ID=PASS,Description=\"All filters passed\">
##INFO=<ID=QD,Number=1,Type=Float,Description=\"q\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"d\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2
";

/// 3 samples; s2 is missing in 9 of 10 records so `--mind` drops it.
fn fixture() -> Vec<u8> {
    let mut text = String::from(HEADER);
    for i in 0..10_u32 {
        let s2 = if i == 0 { "0/1" } else { "./." };
        let _ = writeln!(
            text,
            "Chr01\t{}\t.\tA\tG\t50\tPASS\tQD=10\tGT:DP\t0/1:30\t0/1:30\t{s2}:30",
            i + 1
        );
    }
    let mut deflater = reed::bgzf::deflate::Deflater::new(6).expect("invariant: level");
    let mut out = Vec::new();
    for block in text.as_bytes().chunks(4096) {
        deflater
            .deflate_block(block, &mut out)
            .expect("invariant: block");
    }
    out.extend_from_slice(&reed::bgzf::block::EOF_BLOCK);
    out
}

/// The non-comment lines of a BGZF VCF (the provenance lines differ between
/// the two command lines).
fn record_lines(path: &std::path::Path) -> Vec<String> {
    let bytes = std::fs::read(path).expect("invariant: read output");
    let mut inflater = reed::bgzf::inflate::Inflater::new();
    let mut text = Vec::new();
    let mut rest: &[u8] = &bytes;
    while let Some(total) = reed::bgzf::block::parse_block_size(rest) {
        inflater
            .inflate_block(&rest[..total], &mut text)
            .expect("invariant: inflate");
        rest = &rest[total..];
    }
    String::from_utf8_lossy(&text)
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(str::to_string)
        .collect()
}

fn run(args: &[&str], stdin_bytes: Option<&[u8]>) -> std::process::Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_reed"));
    command
        .args(args)
        .stdin(Stdio::piped())
        .stderr(Stdio::piped());
    let mut child = command.spawn().expect("invariant: spawn reed");
    if let Some(bytes) = stdin_bytes {
        child
            .stdin
            .take()
            .expect("invariant: piped stdin")
            .write_all(bytes)
            .expect("invariant: write stdin");
    }
    child.wait_with_output().expect("invariant: wait")
}

/// The fixture and its stdin form produce identical output; the run leaves
/// no reed file in the temp directory.
#[test]
fn single_pass_stdin_streams_and_leaves_nothing() {
    let dir = std::env::temp_dir().join(format!("reed-o5-single-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: write fixture");

    // The reports + summary pass reads the input exactly once: it must
    // stream from stdin and leave no spill.
    let p_file = dir.join("file");
    let p_stdin = dir.join("stdin");
    let file_out = run(
        &[
            "stats",
            input.to_str().unwrap(),
            "--summary",
            "-o",
            p_file.to_str().unwrap(),
            "--verbosity",
            "0",
        ],
        None,
    );
    assert!(
        file_out.status.success(),
        "{}",
        String::from_utf8_lossy(&file_out.stderr)
    );
    let stdin_out = run(
        &[
            "stats",
            "-",
            "--summary",
            "-o",
            p_stdin.to_str().unwrap(),
            "--verbosity",
            "0",
        ],
        Some(&fixture()),
    );
    assert!(
        stdin_out.status.success(),
        "{}",
        String::from_utf8_lossy(&stdin_out.stderr)
    );
    let file_stats =
        std::fs::read_to_string(format!("{}.stats", p_file.display())).expect("invariant: stats");
    let stdin_stats =
        std::fs::read_to_string(format!("{}.stats", p_stdin.display())).expect("invariant: stats");
    assert_eq!(
        file_stats, stdin_stats,
        "streamed stdin matches the file run"
    );

    let leftovers: Vec<_> = std::fs::read_dir(std::env::temp_dir())
        .expect("invariant: temp dir")
        .filter_map(Result::ok)
        .filter(|e| e.file_name().to_string_lossy().starts_with("reed-stdin-"))
        .collect();
    assert!(
        leftovers.is_empty(),
        "no spill for single-pass runs: {leftovers:?}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// `--mind` reads twice: stdin spills under `--tmp-dir` and
/// the spill is deleted on exit; the sample set matches the file run.
#[test]
fn multi_pass_stdin_spills_and_cleans_up() {
    let dir = std::env::temp_dir().join(format!("reed-o5-multi-{}", std::process::id()));
    let tmp = dir.join("tmp");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&tmp).expect("invariant: temp dir");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: write fixture");
    let out_file = dir.join("file.vcf.gz");
    let out_stdin = dir.join("stdin.vcf.gz");

    let from_file = run(
        &[
            "filter",
            input.to_str().unwrap(),
            "--mind",
            "0.5",
            "-Oz",
            "-o",
            out_file.to_str().unwrap(),
        ],
        None,
    );
    assert!(
        from_file.status.success(),
        "{}",
        String::from_utf8_lossy(&from_file.stderr)
    );
    let from_stdin = run(
        &[
            "filter",
            "-",
            "--mind",
            "0.5",
            "-Oz",
            "-o",
            out_stdin.to_str().unwrap(),
            "--tmp-dir",
            tmp.to_str().unwrap(),
        ],
        Some(&fixture()),
    );
    assert!(
        from_stdin.status.success(),
        "{}",
        String::from_utf8_lossy(&from_stdin.stderr)
    );

    // The provenance line records the command line, which differs; the
    // records (and the rest of the header) must be identical.
    assert_eq!(
        record_lines(&out_file),
        record_lines(&out_stdin),
        "streamed stdin matches the file run"
    );
    let leftovers: Vec<_> = std::fs::read_dir(&tmp)
        .expect("invariant: tmp dir")
        .filter_map(Result::ok)
        .map(|e| e.file_name())
        .collect();
    assert!(
        leftovers.is_empty(),
        "the spill is deleted on exit: {leftovers:?}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v7 D7-2: `--indep-pairwise` runs a prune pass0 that consumes the
/// stream, so stdin must spill even without `--mind`; otherwise pass1 reopens
/// stdin, reads EOF and the run fails with `input is empty`.
#[test]
fn indep_pairwise_stdin_spills_without_mind() {
    let dir = std::env::temp_dir().join(format!("reed-o5-indep-{}", std::process::id()));
    let tmp = dir.join("tmp");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&tmp).expect("invariant: temp dir");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: write fixture");
    let out_file = dir.join("file.vcf.gz");
    let out_stdin = dir.join("stdin.vcf.gz");

    let from_file = run(
        &[
            "filter",
            input.to_str().unwrap(),
            "--indep-pairwise",
            "50",
            "5",
            "0.2",
            "-Oz",
            "-o",
            out_file.to_str().unwrap(),
        ],
        None,
    );
    assert!(
        from_file.status.success(),
        "{}",
        String::from_utf8_lossy(&from_file.stderr)
    );
    let from_stdin = run(
        &[
            "filter",
            "-",
            "--indep-pairwise",
            "50",
            "5",
            "0.2",
            "-Oz",
            "-o",
            out_stdin.to_str().unwrap(),
            "--tmp-dir",
            tmp.to_str().unwrap(),
        ],
        Some(&fixture()),
    );
    assert!(
        from_stdin.status.success(),
        "{}",
        String::from_utf8_lossy(&from_stdin.stderr)
    );
    assert_eq!(
        record_lines(&out_file),
        record_lines(&out_stdin),
        "streamed stdin matches the file run"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// The spill guard deletes the file on early error returns too.
#[test]
fn spill_is_removed_on_error_exit() {
    let dir = std::env::temp_dir().join(format!("reed-o5-error-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    // An invalid spec fails after the input (and the spill) was opened.
    let out = run(
        &[
            "filter",
            "-",
            "--mind",
            "0.5",
            "-Oz",
            "-o",
            dir.join("o.vcf.gz").to_str().unwrap(),
            "--tmp-dir",
            dir.to_str().unwrap(),
            "--maf",
            "9",
        ],
        Some(&fixture()),
    );
    assert!(!out.status.success(), "the threshold is invalid");
    let leftovers: Vec<_> = std::fs::read_dir(&dir)
        .expect("invariant: dir")
        .filter_map(Result::ok)
        .filter(|e| e.file_name().to_string_lossy().starts_with("reed-stdin-"))
        .collect();
    assert!(
        leftovers.is_empty(),
        "spill deleted on the error path: {leftovers:?}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// 04 §2: `-O` picks the output format, `-o` names the file. `-O v` used to
/// force stdout and ignore `-o`, which sent a full 35 GB run into the caller's
/// pipe (and made the `-Ov` benchmarks measure the pipe consumer).
#[test]
fn output_type_and_path_compose() {
    let dir = std::env::temp_dir().join(format!("reed-o-args-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: write fixture");

    let plain = dir.join("out.vcf");
    let plain_run = run(
        &[
            "filter",
            input.to_str().unwrap(),
            "-Ov",
            "-o",
            plain.to_str().unwrap(),
        ],
        None,
    );
    assert!(
        plain_run.status.success(),
        "{}",
        String::from_utf8_lossy(&plain_run.stderr)
    );
    assert!(plain_run.stdout.is_empty(), "-o redirects the output");
    let text = std::fs::read_to_string(&plain).expect("invariant: plain output");
    assert!(text.starts_with("##fileformat="), "{text:.80}");

    let compressed = dir.join("out.dat");
    let gz_run = run(
        &[
            "filter",
            input.to_str().unwrap(),
            "-Oz1",
            "-o",
            compressed.to_str().unwrap(),
        ],
        None,
    );
    assert!(
        gz_run.status.success(),
        "{}",
        String::from_utf8_lossy(&gz_run.stderr)
    );
    assert_eq!(
        &std::fs::read(&compressed).expect("invariant: read")[..2],
        &[0x1f, 0x8b]
    );

    // Without -O the extension still infers the format.
    let inferred = dir.join("out.vcf.gz");
    let inf_run = run(
        &[
            "filter",
            input.to_str().unwrap(),
            "-o",
            inferred.to_str().unwrap(),
        ],
        None,
    );
    assert!(
        inf_run.status.success(),
        "{}",
        String::from_utf8_lossy(&inf_run.stderr)
    );
    assert_eq!(
        &std::fs::read(&inferred).expect("invariant: read")[..2],
        &[0x1f, 0x8b]
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// 04-cleanup §2.1: `stats` accepts the filter's region and sample flags.
/// `-r` (indexed) and `-t` (streamed) must agree, and `-s` must shrink the
/// group denominators.
#[test]
fn stats_accepts_filter_region_and_sample_flags() {
    let dir = std::env::temp_dir().join(format!("reed-stats-flags-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: write fixture");
    // A two-population groups file over the fixture's three samples.
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tpopA\ns1\tpopA\ns2\tpopB\n").expect("invariant: groups");

    // `-r` needs an index next to the input.
    let indexed = run(&["index", input.to_str().unwrap(), "-c"], None);
    assert!(
        indexed.status.success(),
        "{}",
        String::from_utf8_lossy(&indexed.stderr)
    );

    let prefix = |name: &str| dir.join(name).display().to_string();
    let run_stats = |extra: &[&str], out: &str| {
        let mut argv = vec![
            "popg",
            input.to_str().unwrap(),
            "--pop",
            pops.to_str().unwrap(),
            "--pi",
            "--variant-only",
            "-o",
            out,
        ];
        argv.extend_from_slice(extra);
        let output = run(&argv, None);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        std::fs::read_to_string(format!("{out}_pi.txt")).expect("invariant: pi file")
    };
    let all = run_stats(&[], &prefix("all"));
    let region_indexed = run_stats(&["-r", "Chr01:1-5"], &prefix("ri"));
    let region_streamed = run_stats(&["-t", "Chr01:1-5"], &prefix("rs"));
    assert_eq!(region_indexed, region_streamed, "-r and -t must agree");
    assert_ne!(region_indexed, all, "the region must restrict the windows");

    // popA holds s0 and s1; keeping only s0 must halve its comparison counts.
    let subset = run_stats(&["-s", "s0"], &prefix("sub"));
    let comparisons = |text: &str, pop: &str| -> u64 {
        text.lines()
            .skip(1)
            .filter(|l| l.starts_with(pop))
            .map(|l| {
                l.split('\t')
                    .nth(7)
                    .expect("comparisons")
                    .parse::<u64>()
                    .expect("number")
            })
            .sum()
    };
    assert!(
        comparisons(&subset, "popA") < comparisons(&all, "popA"),
        "-s must shrink the popA denominators"
    );
    // An unknown sample is an error, not a silent skip.
    let bad = run(
        &[
            "popg",
            input.to_str().unwrap(),
            "--pop",
            pops.to_str().unwrap(),
            "--pi",
            "-s",
            "nope",
            "-o",
            &prefix("bad"),
        ],
        None,
    );
    assert!(!bad.status.success());
    let _ = std::fs::remove_dir_all(&dir);
}

/// 04-cleanup §2.1: the site/threshold flags on `stats` produce exactly what
/// the composition (`filter` to a file, then `stats`) would, in one pass.
#[test]
fn stats_site_filters_match_the_filter_then_stats_composition() {
    let dir = std::env::temp_dir().join(format!("reed-stats-sites-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: write fixture");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tpopA\ns1\tpopA\ns2\tpopB\n").expect("invariant: groups");
    let filtered = dir.join("filtered.vcf.gz");
    let prefix = |name: &str| dir.join(name).display().to_string();
    for args in [
        vec!["--maf", "0.05"],
        vec!["--min-dp", "20"],
        vec!["--geno", "0.1"],
        vec!["-Q", "0.9"],
    ] {
        let mut filter_argv = vec!["filter", input.to_str().unwrap()];
        filter_argv.extend(args.iter().copied());
        filter_argv.extend(["-Oz", "-o", filtered.to_str().unwrap()]);
        let filtered_run = run(&filter_argv, None);
        assert!(
            filtered_run.status.success(),
            "{}",
            String::from_utf8_lossy(&filtered_run.stderr)
        );

        let stats = |source: &str, extra: &[&str], out: &str| -> String {
            let mut argv = vec![
                "popg",
                source,
                "--pop",
                pops.to_str().unwrap(),
                "--pi",
                "--variant-only",
                "-o",
                out,
            ];
            argv.extend_from_slice(extra);
            let output = run(&argv, None);
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            std::fs::read_to_string(format!("{out}_pi.txt")).expect("invariant: pi file")
        };
        let composed = stats(filtered.to_str().unwrap(), &[], &prefix("composed"));
        let single_pass = stats(input.to_str().unwrap(), &args, &prefix("single"));
        assert_eq!(composed, single_pass, "args {args:?} must agree");
    }
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v7 D7-6: report mode parsed `-s/-S/--keep/--remove` and the
/// site/threshold flags and silently ignored them. It must refuse them
/// instead; `--sites-file` now reaches report mode and narrows the rows.
#[test]
fn report_mode_rejects_unapplied_filters() {
    let dir = std::env::temp_dir().join(format!("reed-reports-flags-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: write fixture");
    let out = dir.join("rep");

    // A sample subset is not applied by report mode → usage error (exit 2).
    let kept = run(
        &[
            "stats",
            input.to_str().unwrap(),
            "--freq",
            "-s",
            "s0",
            "-o",
            out.to_str().unwrap(),
        ],
        None,
    );
    assert_eq!(
        kept.status.code(),
        Some(2),
        "{}",
        String::from_utf8_lossy(&kept.stderr)
    );
    assert!(
        String::from_utf8_lossy(&kept.stderr).contains("does not apply"),
        "{}",
        String::from_utf8_lossy(&kept.stderr)
    );

    // A site/threshold filter likewise.
    let maf = run(
        &[
            "stats",
            input.to_str().unwrap(),
            "--freq",
            "--maf",
            "0.01",
            "-o",
            out.to_str().unwrap(),
        ],
        None,
    );
    assert_eq!(maf.status.code(), Some(2));

    // `--sites-file` is applied: only the listed site's row appears.
    let sites = dir.join("sites.txt");
    std::fs::write(&sites, "Chr01\t3\n").expect("invariant: sites");
    let out2 = dir.join("rep2");
    let ok = run(
        &[
            "stats",
            input.to_str().unwrap(),
            "--freq",
            "--sites-file",
            sites.to_str().unwrap(),
            "-o",
            out2.to_str().unwrap(),
        ],
        None,
    );
    assert!(
        ok.status.success(),
        "{}",
        String::from_utf8_lossy(&ok.stderr)
    );
    let afreq = std::fs::read_to_string(format!("{}.afreq", out2.display())).expect("afreq");
    let rows: Vec<&str> = afreq.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(rows.len(), 1, "only the sites-file site: {afreq:?}");
    assert_eq!(rows[0].split('\t').nth(1), Some("Chr01_3"), "{:?}", rows[0]);
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v7 D7-8: `--f4-ratio A,B;C,D;X` names five populations; the parser
/// used to accept only three and hardcode b = a, c = o.
#[test]
fn f4_ratio_accepts_five_populations() {
    let dir = std::env::temp_dir().join(format!("reed-f4ratio-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let vcf = dir.join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=Chr01,length=100000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1\tP2\tP3\tP4\tP5\n\
         Chr01\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t0/1\t1/1\t0/0\n\
         Chr01\t2\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\t1/1\t0/0\n\
         Chr01\t3\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/1\t0/1\t1/1\n\
         Chr01\t4\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/0\t1/1\t0/1\t0/1\n",
    )
    .expect("invariant: fixture");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "P1\tA\nP2\tO\nP3\tX\nP4\tC\nP5\tB\n").expect("invariant: pops");
    let prefix = dir.join("out");
    let out = run(
        &[
            "popg",
            vcf.to_str().unwrap(),
            "--pop",
            pops.to_str().unwrap(),
            "--f4-ratio",
            "A,O;X,C;B",
            "-o",
            prefix.to_str().unwrap(),
        ],
        None,
    );
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let wrote = std::fs::read_dir(&dir)
        .expect("invariant: dir")
        .filter_map(Result::ok)
        .any(|e| e.file_name().to_string_lossy().starts_with("out"));
    assert!(wrote, "the f4-ratio output file must be written");
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v7 D7-12: a `-S`/`--keep`/`--remove` file that cannot be read is a
/// usage error, not a silent "no filter".
#[test]
fn missing_sample_file_errors() {
    let dir = std::env::temp_dir().join(format!("reed-keepfile-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: write fixture");
    let missing = dir.join("nope.txt");
    let out = run(
        &[
            "popg",
            input.to_str().unwrap(),
            "--pi",
            "-S",
            missing.to_str().unwrap(),
            "-o",
            dir.join("p").to_str().unwrap(),
        ],
        None,
    );
    assert_eq!(
        out.status.code(),
        Some(2),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(
        String::from_utf8_lossy(&out.stderr).contains("cannot read sample file"),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v7 D7-24: a missing input file must name the file in the error.
#[test]
fn missing_input_names_the_file() {
    let missing = "/definitely/not/here.vcf";
    let out = run(&["popg", missing, "--pi", "-o", "/tmp/reed-nope"], None);
    assert_eq!(
        out.status.code(),
        Some(1),
        "an unreadable input is a runtime error"
    );
    assert!(
        String::from_utf8_lossy(&out.stderr).contains(missing),
        "the error must name the file: {}",
        String::from_utf8_lossy(&out.stderr)
    );
}

/// doc_v7 D7-18: an integer r² threshold parses by position
/// (`--indep-pairwise 50 5 1`), instead of being mistaken for the step.
#[test]
fn indep_pairwise_accepts_integer_r2() {
    let dir = std::env::temp_dir().join(format!("reed-indep-int-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    let input = dir.join("in.vcf.gz");
    std::fs::write(&input, fixture()).expect("invariant: write fixture");
    let out = run(
        &[
            "filter",
            input.to_str().unwrap(),
            "--indep-pairwise",
            "50",
            "5",
            "1",
            "-Oz",
            "-o",
            dir.join("o.vcf.gz").to_str().unwrap(),
        ],
        None,
    );
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let _ = std::fs::remove_dir_all(&dir);
}

/// doc_v7 D7-13: `--skip-empty` reaches the windowed route (it was hardcoded
/// false) and `--levels contig` suppresses the window rows.
#[test]
fn skip_empty_and_levels_apply_to_windowed_stats() {
    let dir = std::env::temp_dir().join(format!("reed-lv-{}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("invariant: temp dir");
    // Sites at 1 and 10 leave eight empty 1 bp windows between them.
    let vcf = dir.join("in.vcf");
    std::fs::write(
        &vcf,
        "##fileformat=VCFv4.2\n\
         ##contig=<ID=Chr01,length=100000>\n\
         ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
         #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts0\ts1\ts2\n\
         Chr01\t1\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\t0/0\n\
         Chr01\t10\t.\tA\tG\t50\tPASS\t.\tGT\t0/0\t0/1\t0/1\n",
    )
    .expect("invariant: fixture");
    let pops = dir.join("pops.txt");
    std::fs::write(&pops, "s0\tpopA\ns1\tpopA\ns2\tpopA\n").expect("invariant: pops");
    let run_pi = |name: &str, extra: &[&str]| -> String {
        let prefix = dir.join(name);
        let mut argv = vec![
            "popg",
            vcf.to_str().unwrap(),
            "--pop",
            pops.to_str().unwrap(),
            "--pi",
            "--variant-only",
            "--size",
            "1",
            "--step",
            "1",
            "-o",
        ];
        let prefix_s = prefix.display().to_string();
        argv.push(&prefix_s);
        argv.extend_from_slice(extra);
        let out = run(&argv, None);
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        std::fs::read_to_string(format!("{prefix_s}_pi.txt")).expect("invariant: pi file")
    };
    let data_rows = |text: &str| text.lines().filter(|l| !l.starts_with("pop\t")).count();

    let all = run_pi("all", &["--levels", "window,contig"]);
    let empty = run_pi("empty", &["--levels", "window,contig", "--skip-empty"]);
    let contig = run_pi("contig", &["--levels", "contig"]);
    assert!(
        data_rows(&empty) < data_rows(&all),
        "--skip-empty must drop the empty windows: {} vs {}",
        data_rows(&empty),
        data_rows(&all)
    );
    assert_eq!(
        data_rows(&contig),
        1,
        "--levels contig writes only the contig row: {contig:?}"
    );
    let _ = std::fs::remove_dir_all(&dir);
}
