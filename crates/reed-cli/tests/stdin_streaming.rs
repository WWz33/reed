//! O5 gates (03-performance.md): stdin streams for single-pass commands and
//! leaves no temporary file behind; multi-pass runs spill under `--tmp-dir`
//! and delete the spill on exit (04 §C6).

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

/// 3 samples; s2 is missing in 9 of 10 records so the sample-level call-rate
/// rule drops it.
fn fixture() -> Vec<u8> {
    let mut text = String::from(HEADER);
    for i in 0..10_u32 {
        let s2 = if i == 0 { "0/1" } else { "./." };
        text.push_str(&format!(
            "Chr01\t{}\t.\tA\tG\t50\tPASS\tQD=10\tGT:DP\t0/1:30\t0/1:30\t{s2}:30\n",
            i + 1
        ));
    }
    let mut deflater = reed_io::bgzf::deflate::Deflater::new(6).expect("invariant: level");
    let mut out = Vec::new();
    for block in text.as_bytes().chunks(4096) {
        deflater.deflate_block(block, &mut out).expect("invariant: block");
    }
    out.extend_from_slice(&reed_io::bgzf::block::EOF_BLOCK);
    out
}

/// The non-comment lines of a BGZF VCF (the provenance lines differ between
/// the two command lines).
fn record_lines(path: &std::path::Path) -> Vec<String> {
    let bytes = std::fs::read(path).expect("invariant: read output");
    let mut inflater = reed_io::bgzf::inflate::Inflater::new();
    let mut text = Vec::new();
    let mut rest: &[u8] = &bytes;
    while let Some(total) = reed_io::bgzf::block::parse_block_size(rest) {
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
    command.args(args).stdin(Stdio::piped()).stderr(Stdio::piped());
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

    // qc reads the input exactly once: it must stream from stdin.
    let file_out = run(&["qc", input.to_str().unwrap()], None);
    assert!(file_out.status.success(), "{}", String::from_utf8_lossy(&file_out.stderr));
    let stdin_out = run(&["qc", "-"], Some(&fixture()));
    assert!(stdin_out.status.success(), "{}", String::from_utf8_lossy(&stdin_out.stderr));
    assert_eq!(file_out.stdout, stdin_out.stdout, "streamed stdin matches the file run");

    let leftovers: Vec<_> = std::fs::read_dir(std::env::temp_dir())
        .expect("invariant: temp dir")
        .filter_map(Result::ok)
        .filter(|e| e.file_name().to_string_lossy().starts_with("reed-stdin-"))
        .collect();
    assert!(leftovers.is_empty(), "no spill for single-pass runs: {leftovers:?}");
    let _ = std::fs::remove_dir_all(&dir);
}

/// `--max-missing-sample` reads twice: stdin spills under `--tmp-dir` and
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
            "--max-missing-sample",
            "0.5",
            "-Oz",
            "-o",
            out_file.to_str().unwrap(),
        ],
        None,
    );
    assert!(from_file.status.success(), "{}", String::from_utf8_lossy(&from_file.stderr));
    let from_stdin = run(
        &[
            "filter",
            "-",
            "--max-missing-sample",
            "0.5",
            "-Oz",
            "-o",
            out_stdin.to_str().unwrap(),
            "--tmp-dir",
            tmp.to_str().unwrap(),
        ],
        Some(&fixture()),
    );
    assert!(from_stdin.status.success(), "{}", String::from_utf8_lossy(&from_stdin.stderr));

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
    assert!(leftovers.is_empty(), "the spill is deleted on exit: {leftovers:?}");
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
            "--max-missing-sample",
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
    assert!(leftovers.is_empty(), "spill deleted on the error path: {leftovers:?}");
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
        &["filter", input.to_str().unwrap(), "-Ov", "-o", plain.to_str().unwrap()],
        None,
    );
    assert!(plain_run.status.success(), "{}", String::from_utf8_lossy(&plain_run.stderr));
    assert!(plain_run.stdout.is_empty(), "-o redirects the output");
    let text = std::fs::read_to_string(&plain).expect("invariant: plain output");
    assert!(text.starts_with("##fileformat="), "{text:.80}");

    let compressed = dir.join("out.dat");
    let gz_run = run(
        &["filter", input.to_str().unwrap(), "-Oz1", "-o", compressed.to_str().unwrap()],
        None,
    );
    assert!(gz_run.status.success(), "{}", String::from_utf8_lossy(&gz_run.stderr));
    assert_eq!(&std::fs::read(&compressed).expect("invariant: read")[..2], &[0x1f, 0x8b]);

    // Without -O the extension still infers the format.
    let inferred = dir.join("out.vcf.gz");
    let inf_run = run(
        &["filter", input.to_str().unwrap(), "-o", inferred.to_str().unwrap()],
        None,
    );
    assert!(inf_run.status.success(), "{}", String::from_utf8_lossy(&inf_run.stderr));
    assert_eq!(&std::fs::read(&inferred).expect("invariant: read")[..2], &[0x1f, 0x8b]);
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
    assert!(indexed.status.success(), "{}", String::from_utf8_lossy(&indexed.stderr));

    let prefix = |name: &str| dir.join(name).display().to_string();
    let run_stats = |extra: &[&str], out: &str| {
        let mut argv = vec![
            "stats",
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
        assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
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
            .map(|l| l.split('\t').nth(7).expect("comparisons").parse::<u64>().expect("number"))
            .sum()
    };
    assert!(
        comparisons(&subset, "popA") < comparisons(&all, "popA"),
        "-s must shrink the popA denominators"
    );
    // An unknown sample is an error, not a silent skip.
    let bad = run(
        &["stats", input.to_str().unwrap(), "--pop", pops.to_str().unwrap(), "--pi", "-s", "nope",
          "-o", &prefix("bad")],
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
        vec!["--max-missing", "0.9"],
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
                "stats",
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
            assert!(output.status.success(), "{}", String::from_utf8_lossy(&output.stderr));
            std::fs::read_to_string(format!("{out}_pi.txt")).expect("invariant: pi file")
        };
        let composed = stats(filtered.to_str().unwrap(), &[], &prefix("composed"));
        let single_pass = stats(input.to_str().unwrap(), &args, &prefix("single"));
        assert_eq!(composed, single_pass, "args {args:?} must agree");
    }
    let _ = std::fs::remove_dir_all(&dir);
}
